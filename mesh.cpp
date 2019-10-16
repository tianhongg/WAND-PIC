//----------------------------------------------------------------------------------||
//-------------------                mesh.cpp                    -------------------||
//----------------------------------------------------------------------------------||
//                                                                                  ||
//               __        ___    _   _ ____        ____ ___ ____                   ||
//               \ \      / / \  | \ | |  _ \      |  _ \_ _/ ___|                  ||
//                \ \ /\ / / _ \ |  \| | | | |_____| |_) | | |                      ||
//                 \ V  V / ___ \| |\  | |_| |_____|  __/| | |___                   ||
//                  \_/\_/_/   \_\_| \_|____/      |_|  |___\____|                  ||
//                                                                                  ||
//----------------------------------------------------------------------------------||
//--  (W)akefield (A)cceleration a(n)d (D)LA - (P)article (i)n (C)ell Simulation  --||
//----------------------------------------------------------------------------------||
//---Author-----------           : Tianhong Wang                --------------------||
//---Starting---------           : Jan-24-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||


#include "wand_PIC.h"
#include <exception>



//---------------------------- Mesh::Mesh --------------------
Mesh::Mesh(int XGridN, int YGridN, int ZGridN, FILE *f): NList ("Plasma")
{
	p_Trajectory = NULL;
	p_Particle	 = NULL;
	p_CellArray  = NULL;
	dAx1 = dAx2 = dAy1 = dAy2 = NULL;

	Rank = p_domain()->p_Partition()->Get_Rank();

	dx = p_domain()->Get_dx();
	dy = p_domain()->Get_dy();
	dz = p_domain()->Get_dz();
	dt = p_domain()->Get_dt();
	dzz = dz;

	GridX = XGridN;
	GridY = YGridN;
	GridZ = ZGridN;
	CellNum = (GridX+2)*(GridY+2)*ZGridN;  //Extra Cells Overlapping With Neighboring Partitions;

	//Big Coordinates of the Rank
	RankIdx_X = p_domain()->p_Partition()->RankIdx_X(); 
	RankIdx_Y = p_domain()->p_Partition()->RankIdx_Y(); 

	// Coordinates of the Mesh Bottom_Left Corner. 
	//(0,0) origin is defined at the domain Bottom_Left Corner. 
	Offset_X = (RankIdx_X-1)*(GridX*dx)-p_domain()->Get_Xmax()*0.5;
	Offset_Y = (RankIdx_Y-1)*(GridY*dy)-p_domain()->Get_Ymax()*0.5;

	//Seed Cells
	try
  	{
	p_CellArray = new Cell[CellNum];
	}
  	catch (std::exception& e)
  	{
  	std::cout << "==== Mesh: Creating Cell Failed.         ====\n";
  	std::cout << "==== Maybe Too Many Cells per Processor, ====\n";
  	std::cout << "==== Or Wrong Cell Number                ====\n";
    std::cout << e.what() << '\n';
    exit(1);
  	}
	// Pointer to the Cell Object.


  	// set shifted position for ions
  	for (int k=0; k<GridZ; k++)
	{
		for (int j=1; j<=GridY; j++)
		{
			for (int i=1; i<=GridX; i++)
			{
				Cell &c = GetCell(i,  j,  k);
				c.Z_shifted = CellZ(k);
			}
		}
	}


  	AddEntry((char*)"PlasProfile_T", 	&PProfileT, 1);
  	AddEntry((char*)"PlasProfile_L", 	&PProfileL, 1);
  	AddEntry((char*)"Traj_per_Cellx", 	&TpCellx, 1);
  	AddEntry((char*)"Traj_per_Celly", 	&TpCelly, 1);
  	AddEntry((char*)"Push_Traj_Order",	&PushOrder, 1);
  	AddEntry((char*)"AdaptiveStep",		&AdaptiveStep,0);
  	AddEntry((char*)"Threshold_V",		&V_thresh,1.0);
  	AddEntry((char*)"IfIonization",	&if_ioniz,0);
  	AddEntry((char*)"Ion_Radius",	&Ion_R, 0.0);
  	
  	AddEntry((char*)"Dop_Begin",	&Dop_TB,0.0);
  	AddEntry((char*)"Dop_End",		&Dop_TE,0.0);
  	
  	AddEntry((char*)"PlasmaDen",	&Pla_ne,1e18);
  	AddEntry((char*)"DopingRate",	&Dop_ne,0.01);


  	AddEntry((char*)"PlateauBegin",	&PlateauBegin,0);
  	AddEntry((char*)"PlasmaBegin",	&PlasmaBegin, 0);
  	AddEntry((char*)"PlateauEnd",	&PlateauEnd,  1e10);
  	AddEntry((char*)"PlasmaEnd",	&PlasmaEnd,	  1e10);



	if (f)
    {
      rewind(f);
      read(f);
    }


    if(TpCellx == 0 || TpCelly == 0)
  	{
   		if (Rank==0)  std::cout << "==== Mesh: Wrong Trajectory Per Cell     ====\n";
  		exit(-10);
  	}


  	if(PushOrder<0 ||PushOrder>2)
  	{
   		if (Rank==0)  std::cout << "==== Mesh: Wrong Trajectory Pushing Order====\n";
  		exit(-10);
  	}



}


void Mesh::SeedTrajectory()
{
	int Xpa= p_domain()->p_Partition()->GetXpart();
	int Ypa= p_domain()->p_Partition()->GetYpart();

	double ztime = p_domain()->Get_RunTime();

	double SeedOffX;
	double SeedOffY;

	int TgridX;
	int TgridY;

	double xt;
	double yt;
	SeedOffX = Offset_X;
	SeedOffY = Offset_Y;
	TgridX= GridX;
	TgridY= GridY;
	
	TrajNum = TgridX*TgridY*TpCellx*TpCelly;

	Trajectory *p =NULL;

	for (int j=0; j<TgridY*TpCelly; j++)
	{
		for (int i=0; i<TgridX*TpCellx; i++)
		{
			xt = SeedOffX + double(i + 0.5)*dx/TpCellx;
			yt = SeedOffY + double(j + 0.5)*dy/TpCelly;
			p = new Trajectory(xt, yt, ztime, TpCellx, TpCelly);
		}
	}

	Vlim=V_thresh*dx/dz;
	Vmax=0.0;
	SetIonDensity();
	return;
}

void Mesh::SetIonDensity()
{

	double ztime=p_domain()->Get_RunTime();

	for (int k=0; k<GridZ; k++) 
	{
		for (int j=0; j<GridY+2; j++) 
		{
			double y = CellY(j);

			for (int i=0; i<GridX+2; i++)
			{
				
				double x = CellX(i);
				Cell &c = GetCell(i, j, k);
				c.W_Deni = ProfileLongi(x,y,ztime)*ProfileTrans(x,y,ztime);
			}
		}
	}

	return;

}


//=============================================================
//=================Trajectory Chain============================
//=============================================================
//==== p0 -> p1 -> p2 -> p3 -> p4 -> p5 -> p6 -> p7 ...............
void Mesh::AddTrajectory(Trajectory* p_Traj)
{
	p_Traj->p_PrevTraj = p_Trajectory;
	p_Trajectory = p_Traj;  

	if(p_Traj->p_PrevTraj)
	{
		p_Traj->p_PrevTraj->p_NextTraj = p_Traj;
	}
	return;
}


void Mesh::AddParticle(Particle* p_Part)
{
	p_Part->p_PrevPart = p_Particle; 
	p_Particle = p_Part;  

	if(p_Part->p_PrevPart)
	{
		p_Part->p_PrevPart->p_NextPart = p_Part;
	}
	return;
}

void Mesh::SetdAs()
{
  	//========allocate matrix for dA=======
	int NF = p_domain()->NFreqs;
	int nn = (GridX+2)*(GridY+2)*NF;

  	dAx1 = new dcomplex[nn];
	dAx2 = new dcomplex[nn];
	for(int i=0; i<nn; i++)
	{	
		dAx1[i] = 0.0;
		dAx2[i] = 0.0;	
	}

  	dAy1 = new dcomplex[nn];
  	dAy2 = new dcomplex[nn];
  	for(int i=0; i<nn; i++)
	{	
		dAy1[i] = 0.0;
		dAy2[i] = 0.0;	
	}
  	dA0 = 1.5/dz/dt;
  	//======================================
  	return;
}

void Mesh::ResetPlasma()
{
	Trajectory *p = NULL;
	Trajectory *pm = NULL;
	p = p_Trajectory;

	//====delete all trajectories;
	while (p)
	{
		pm = p->p_PrevTraj;
		delete p;
		p = pm;
	}

	//====Seed New Trajectories=====
	p_Trajectory = NULL;
	SeedTrajectory();


	int NF = p_domain()->NFreqs;
	int nn = (GridX+2)*(GridY+2)*NF;
	//========= Set d_A=0  =========

	for(int i=0; i<nn; i++)
	{	
		dAx1[i] = 0.0; 
		dAx2[i] = 0.0;	
	}
	

  	for(int i=0; i<nn; i++)
	{	
		dAy1[i] = 0.0;
		dAy2[i] = 0.0;	
	}
  	
	return;

}

 // initial stage 
// simple version of Ionization block
// tianhong sep-13-2019
// test.


void Mesh::Ionization()
{

	double ExR, EyR, EzR, EL_Inver, EL;
	int i, j, k, k2, NF;
	double x,y,z;
	double wm, wp;
	double OmegaL;
	double proba;
	double xtemp,ytemp,ztemp, weight;
	int NFreqs=p_domain()->NFreqs;


	//====temporary ion data=================//
	double UH  =13.6;		
	double UIon=871.4;  //ionzaition potential eV
	int Z_Ion=8; 		  //ion charge after ionization
	//========================================


	double E_natlog=2.718281828459;
	double omega_alpha= 733.3*sqrt(1e18/Pla_ne);// atomic unit frequency in kp
	double E_alpha=omega_alpha/137.0;

	double neff =Z_Ion*sqrt(UH/UIon);  //n_star;

	srand(time(NULL));
	Particle *p =NULL;


	for (k=0; k<GridZ; k++)
	{
		for (j=1; j<=GridY; j++)
		{
			for (i=1; i<=GridX; i++)
			{
				
				x = CellX(i);
				y = CellY(j);
				Cell &c = GetCell(i,  j,  k);
				z = c.Z_shifted;// shifted_Position;

				k2 = floor(z/dz);
				if( k2<0 || k2>(GridZ-2)) continue;

				Cell &cm = GetCell(i,  j,  k2);		wm = (k2+1-z/dz);
				Cell &cp = GetCell(i,  j,  k2+1);	wp = (z/dz-k2);

         		EL  = 0.0;
         		ExR = 0.0;
         		EyR = 0.0;
         		EzR = 0.0;
         	
         		if((x*x+y*y)<Ion_R*Ion_R&&c.InoState)
				{
				///////
					for(NF=0; NF<NFreqs; NF++)
					{
						OmegaL = p_domain()->OmegaL[NF];
						EL += sqrt(abs((wm*cm.L_Ex[NF]+wp*cp.L_Ex[NF]))*abs((wm*cm.L_Ex[NF]+wp*cp.L_Ex[NF]))
						          +abs((wm*cm.L_Ey[NF]+wp*cp.L_Ey[NF]))*abs((wm*cm.L_Ey[NF]+wp*cp.L_Ey[NF]))); 
						ExR += ( (wm*cm.L_Ex[NF]+wp*cp.L_Ex[NF]) *exp(-ci*OmegaL*z)).real();
						EyR += ( (wm*cm.L_Ey[NF]+wp*cp.L_Ey[NF]) *exp(-ci*OmegaL*z)).real();
						EzR += ( (wm*cm.L_Ez[NF]+wp*cp.L_Ez[NF]) *exp(-ci*OmegaL*z)).real();
					}
					//===caculate ionization probability===
					if(EL>1e-2)
					{
						EL_Inver=E_alpha/EL;
						proba=dt*omega_alpha*0.311*pow(E_natlog,2*neff)*Z_Ion*Z_Ion/pow(neff,4.5)
		 				*pow(4*EL_Inver*Z_Ion*Z_Ion*Z_Ion/neff/neff/neff/neff,2*neff-1.5)*exp(-2/3.*EL_Inver
		 				*(Z_Ion*Z_Ion*Z_Ion/neff/neff/neff));

		 				double seedr=(double) rand() / (RAND_MAX);
		 				if(1-exp(-proba)>seedr)
		 				{
							xtemp = Offset_X + (i - 0.5)*dx;
							ytemp = Offset_Y + (j - 0.5)*dy;
							ztemp = z;
							weight = Dop_ne;
							p = new Electron(xtemp, ytemp, ztemp, 0.0, 0.0, 0.0, ExR, EyR, EzR, 1.0, weight);
							c.InoState=0;
 						}
				////////
					}

				}

			}
		}
	}


	//shift all ion position;
	for (k=0; k<GridZ; k++)
	{
		for (j=1; j<=GridY; j++)
		{
			for (i=1; i<=GridX; i++)
			{
				
				Cell &c = GetCell(i,  j,  k);
				c.Z_shifted += dt;

				if(c.Z_shifted>(GridZ-1)*dz)
				{
					double dzz=c.Z_shifted-GridZ*dz;
					c.Z_shifted=dzz;
					c.InoState =1;

				}

			}
		}
	}

	
	


	return;
}



Mesh::~Mesh() 
{
	delete[] p_CellArray;
	delete 	 p_Trajectory;
	delete	 p_Particle;
	delete[] dAx1;
	delete[] dAx2;
	delete[] dAy1;
	delete[] dAy2;
};
