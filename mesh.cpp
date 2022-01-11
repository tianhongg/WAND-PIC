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
	XRayDetector = NULL;
	dAx1 = dAx2 = dAy1 = dAy2 = NULL;

	Rank = p_domain()->p_Partition()->Get_Rank();

	dx = p_domain()->Get_dx();
	dy = p_domain()->Get_dy();

	dz = p_domain()->Get_dz();
	dt = p_domain()->Get_dt();
	dt0 = dt;
	dzz = dz;

	minGamma=1.0;
	GridX = XGridN;
	GridY = YGridN;
	GridZ = ZGridN;
	CellNum = (GridX+2)*(GridY+2)*ZGridN;  //Extra Cells Overlapping With Neighboring Partitions;


	int N_Xpart=p_domain()->p_Partition()->GetXpart();
	int N_Ypart=p_domain()->p_Partition()->GetYpart();

	//
	GridsTmp=std::vector<WDOUBLE>(XGridN*N_Xpart+2,0);
	std::vector<WDOUBLE> GridsAcc(XGridN*N_Xpart+2,0);

	WDOUBLE ddxx=0;
	for(int i=0;i<XGridN*N_Xpart/2+1;i++)
	{	
		WDOUBLE dd=p_domain()->CustomGrid(ddxx);
		dd= floor(dd*1e4)/1e4;
		GridsTmp[XGridN*N_Xpart/2-i]   =dd;
		GridsTmp[XGridN*N_Xpart/2+i+1] =dd;
		ddxx+=dd;
	}

	for(int i=1;i<XGridN*N_Xpart+2;i++)
	{
		GridsAcc[i]=GridsAcc[i-1]+GridsTmp[i];
	}

	//set up the new size of the domain;
	p_domain()->Set_Xmax(GridsAcc[XGridN*N_Xpart]);
	p_domain()->Set_Ymax(GridsAcc[XGridN*N_Xpart]);

	//Big Coordinates of the Rank
	RankIdx_X = p_domain()->p_Partition()->RankIdx_X(); 
	RankIdx_Y = p_domain()->p_Partition()->RankIdx_Y(); 

	// Coordinates of the Mesh Bottom_Left Corner. 
	Offset_X = GridsAcc[(RankIdx_X-1)*(GridX)] - p_domain()->Get_Xmax()*0.5;
	Offset_Y = GridsAcc[(RankIdx_Y-1)*(GridY)] - p_domain()->Get_Ymax()*0.5;


	if(Rank==0)
	{
		char name[128];
  		sprintf(name,"Grids.ini");
		FILE * dFile;
		dFile = fopen (name,"w");
		for (int j=1; j<=XGridN*N_Xpart; j++) fprintf(dFile, "%f ", GridsTmp[j]);
		fclose (dFile);
	}


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
		for (int j=0; j<=GridY+1; j++)
		{
			for (int i=0; i<=GridX+1; i++)
			{
				Cell &c = GetCell(i,  j,  k);

				//assign grid size....
				c.dx=GridsTmp[(RankIdx_X-1)*(GridX)+i];
				c.dy=GridsTmp[(RankIdx_Y-1)*(GridY)+j];
				// coordinates of the center of cell...
				c.Xcord = GridsAcc[(RankIdx_X-1)*(GridX)+i]- p_domain()->Get_Xmax()*0.5-c.dx*0.5;
				c.Ycord = GridsAcc[(RankIdx_Y-1)*(GridY)+j]- p_domain()->Get_Ymax()*0.5-c.dy*0.5;
				c.Z_shifted = CellZ(k);
			}
		}
	}

  	AddEntry((char*)"PlasProfile_T", 	&PProfileT, 1);
  	AddEntry((char*)"PlasProfile_L", 	&PProfileL, 1);
  	AddEntry((char*)"Traj_per_Cellx", 	&TpCellx, 1);
  	AddEntry((char*)"Traj_per_Celly", 	&TpCelly, 1);
  	AddEntry((char*)"Delta_P", 			&Delta_P, 0.0);
  	AddEntry((char*)"Push_Traj_Order",	&PushOrder, 1);
  	AddEntry((char*)"AdaptiveStep",		&AdaptiveStep,0);
  	AddEntry((char*)"Threshold_V",		&V_thresh,1.0);
  	AddEntry((char*)"AdjustPsi", 		&IfAdjustPsi,0);
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
  	AddEntry((char*)"PlasRadius",	&PlasRadius,  10000);

/*
  	char name[128];
  	sprintf(name,"Grids_%d_.dg",Rank);
	FILE * dFile;
	dFile = fopen (name,"w");

	for (int j=0; j<=GridY+1; j++)
	{
		for (int i=0; i<=GridX+1; i++)
		{
			Cell &c = GetCell(i,  j,  1);
			fprintf(dFile, "%f ", c.Xcord);
		}
		fprintf(dFile, "\n");
	}
	fclose (dFile);
*/
  	

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


  	// Radiation Detector
  	XRayDetector = new Detector(f);

}


void Mesh::SeedTrajectory()
{
	int Xpa= p_domain()->p_Partition()->GetXpart();
	int Ypa= p_domain()->p_Partition()->GetYpart();

	WDOUBLE ztime = p_domain()->Get_RunTime();

	WDOUBLE xt;
	WDOUBLE yt;

	WDOUBLE dxp;
	WDOUBLE dyp;

	srand(time(NULL));
	
	TrajNum = GridX*GridY*TpCellx*TpCelly;

	Trajectory *p =NULL;


	// char name[128];
 	// sprintf(name,"Trajs_%d_.dg",Rank);
	// FILE * dFile;
	// dFile = fopen (name,"w");

	//loop cell
	for (int j=1; j<=GridY; j++)
	{
		for (int i=1; i<=GridX; i++)
		{
			//seed
			for(int sj=0;sj<TpCelly;sj++)
			{
				for(int si=0;si<TpCellx;si++)
				{
					Cell &c = GetCell(i, j, 0); 

					dxp = c.dx/TpCellx;
					dyp = c.dy/TpCelly;
					
					xt = c.Xcord-c.dx*0.5 + WDOUBLE(si + 0.5)*dxp;
					yt = c.Ycord-c.dy*0.5 + WDOUBLE(sj + 0.5)*dyp;

					p = new Trajectory(xt, yt, ztime, TpCellx, TpCelly, dxp, dyp);
					p->idx_i=i;
					p->idx_j=j;

					// // finite temperature section: test
					WDOUBLE r1 = rand_gaussian(Delta_P);
					WDOUBLE r2 = rand_gaussian(Delta_P);

 					p->Vx=p->old_vx=r1;  
 					p->Vy=p->old_vy=r2;    
                    p->Vxx=r1*r1; 
                    p->Vyy=r2*r2; 
                    p->Vxy=r1*r2;

				}

			}

		}
	}

	// fclose (dFile);
	Vlim=V_thresh*dx/dz;
	Vmax=0.0;
	SetIonDensity();
	return;
}

void Mesh::SetIonDensity()
{
	WDOUBLE ztime=p_domain()->Get_RunTime();

	for (int k=0; k<GridZ; k++) 
	{
		for (int j=0; j<GridY+2; j++) 
		{
			for (int i=0; i<GridX+2; i++)
			{
				Cell &c = GetCell(i, j, k);
				c.W_Deni = ProfileLongi(c.Xcord,c.Ycord,ztime)*ProfileTrans(c.Xcord,c.Ycord,ztime);
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

	WDOUBLE ExR, EyR, EzR, EL_Inver, EL;
	int i, j, k, k2, NF;
	WDOUBLE x,y,z;
	WDOUBLE wm, wp;
	WDOUBLE OmegaL;
	WDOUBLE proba;
	WDOUBLE xtemp,ytemp,ztemp, weight;
	int NFreqs=p_domain()->NFreqs;


	//====temporary ion data=================//
	WDOUBLE UH  =13.6;		
	WDOUBLE UIon=871.4;  //ionzaition potential eV
	int Z_Ion=8; 		  //ion charge after ionization
	//========================================


	WDOUBLE E_natlog=2.718281828459;
	WDOUBLE omega_alpha= 733.3*sqrt(1e18/Pla_ne);// atomic unit frequency in kp
	WDOUBLE E_alpha=omega_alpha/137.0;

	WDOUBLE neff =Z_Ion*sqrt(UH/UIon);  //n_star;

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

		 				WDOUBLE seedr=(WDOUBLE) rand() / (RAND_MAX);
		 				if(1-exp(-proba)>seedr)
		 				{
							xtemp = Offset_X + (i - 0.5)*dx;
							ytemp = Offset_Y + (j - 0.5)*dy;
							ztemp = z;
							weight = c.W_Deni*Dop_ne;
							p = new Electron(xtemp, ytemp, ztemp, 0.0, 0.0, 0.0, ExR, EyR, EzR, 1.0, weight);
							c.InoState=0;
 						}

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
					WDOUBLE dzz=c.Z_shifted-GridZ*dz;
					c.Z_shifted=dzz;
					c.InoState =1;

				}

			}
		}
	}

	
	


	return;
}




Detector::Detector(FILE *f): NList ("Detector")
{
	p_DetectArray=NULL;


	int rank = p_domain()->p_Partition()->Get_Rank();
	
  	AddEntry((char*)"Radiation", 	&IfRadiation, 0);

  	AddEntry((char*)"ThetaMax", 	&ThetaMax, 0.2);
  	AddEntry((char*)"ThetaGrid", 	&NTheta, 150);

  	AddEntry((char*)"PhiMax", 		&PhiMax, 0.2);
  	AddEntry((char*)"PhiGrid",		&NPhi, 150);

  	AddEntry((char*)"OmegaMin",		&OmegaMin, 0);
  	AddEntry((char*)"OmegaMax",		&OmegaMax, 50);
  	AddEntry((char*)"OmegaGrid",	&NOmega,   200);


	if (f)
    {
      rewind(f);
      read(f);
    }


    if(IfRadiation>0 && NTheta>0 && NPhi>0 &&NOmega>0)
    {

    	try
    	{
  
			p_DetectArray = new WDOUBLE[NTheta*NPhi*NOmega];
			if (rank==0)  std::cout << "==== Mesh: Radiation Detector Created.   ====\n";
		}
  		catch (std::exception& e)
  		{
  			std::cout << "==== Mesh: Creating Cell Failed.         ====\n";
  			std::cout << "==== Maybe Too Many Points per Processor ====\n";
    		std::cout << e.what() << '\n';
    		exit(1);
  		}
    		
    }
    else
    {
    		if (rank==0)  std::cout << "==== Mesh: No Radiation Detector.        ====\n";
    }

 	if(IfRadiation>0 && NTheta>0 && NPhi>0 &&NOmega>0)
    {

 		for (int k=0; k<NPhi; k++)
		{
			for (int j=0; j<NTheta; j++)
			{
				for (int i=0; i<NOmega; i++)
				{
					p_DetectArray[GetDectN(i,j,k)]=0;
				}
			}
		}

	}

}


WDOUBLE Mesh::rand_gaussian (WDOUBLE sigma) //sigma=standard deviation 
{
	WDOUBLE x, y, r2;
	do
	{
		x = (2.*rand()-RAND_MAX)/RAND_MAX;
		y = (2.*rand()-RAND_MAX)/RAND_MAX;
		r2 = x * x + y * y;
	}
	while (r2 > 1.0 || r2 == 0);
	return sigma * y * sqrt (-2.0 * log (r2) / r2);
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
	delete XRayDetector;
};

Detector::~Detector() 
{
	if(p_DetectArray) delete[] p_DetectArray;
};
