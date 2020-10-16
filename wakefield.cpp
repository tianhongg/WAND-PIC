//----------------------------------------------------------------------------------||
//-------------------                wakefield.cpp               -------------------||
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
//---Starting---------           : Feb-01-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#include "wand_PIC.h"



void Mesh::MacroSource(bool exbeam, int k)
{

	Trajectory *p = NULL;
	double ddx;
	double ddy;

	double wmm, wmp, wpm, wpp;

	double dxdy = dx*dy;
	int i,j;


	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;


	SetSourceZero(k);
	
	p = p_Trajectory;
	while (p)
	{
		double xt 		=  p->x;
		double yt 		=  p->y;
		double massweig = (p->mass)*(p->Weight);

		//------------------------------
		//      _______
		//     |mp | pp|
		//     |___|___|
		//     |mm | pm|
		//     |___|___|
		//
		//------------------------------
		ddx = xt-(Offset_X-dx*0.5);
		ddy = yt-(Offset_Y-dy*0.5);
		// idex of the corner cell
		i = floor(ddx/dx);
		j = floor(ddy/dy);


		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(i < 0 || i > GridX || j < 0 || j > GridY)
		{
			p = p->p_PrevTraj;
			continue;
		}
		//==================================================

		
		wmm = (i+1-ddx/dx)*(j+1-ddy/dy);
		wmp = (i+1-ddx/dx)*(ddy/dy-j);
		wpm = (ddx/dx-i)*(j+1-ddy/dy);
		wpp = (ddx/dx-i)*(ddy/dy-j);

		Cell &cmm = GetCell(i,j,k);
		Cell &cmp = GetCell(i,j+1,k);
		Cell &cpm = GetCell(i+1,j,k);
		Cell &cpp = GetCell(i+1,j+1,k);

		for (int n=0; n<SOU_DIM; n++)
		{
		cmm.W_Source[n] += massweig * wmm * p->T_Source[n];
		cmp.W_Source[n] += massweig * wmp * p->T_Source[n];
		cpm.W_Source[n] += massweig * wpm * p->T_Source[n];
		cpp.W_Source[n] += massweig * wpp * p->T_Source[n];
		}
		p = p->p_PrevTraj;

	}
	
	AdjustSource(exbeam,k);

	return;

}

void Mesh::SetSourceZero(int k)
{
	for(int j=0; j<=GridY+1; j++)
	{
		for(int i=0; i<=GridX+1; i++)
		{
			Cell &ccc = GetCell(i,j,k);
			for (int n=0; n<SOU_DIM; n++)
			{
				ccc.W_Source[n] =0.0;
			}

		}
	}	
	return;
}



void Mesh::AdjustSource(bool exbeam, int k)
{
	int Xpa = p_domain()->p_Partition()->GetXpart();
	int Ypa = p_domain()->p_Partition()->GetYpart();

	//=====mm-corner======================
	Cell &c1  = GetCell( 1,1,k);
	Cell &co1 = GetCell( 0,0,k);

	int b=0;
	if(exbeam) b=1;

	for (int i=0; i<SOU_DIM+BEA_DIM*b; i++)
	{ c1.W_Source[i] += co1.W_Source[i]; }

	//=====mp-corner======================
	Cell &c2  = GetCell( 1,GridY,k);
	Cell &co2 = GetCell( 0,GridY+1,k);

	for (int i=0; i<SOU_DIM+BEA_DIM*b; i++)
	{ c2.W_Source[i] += co2.W_Source[i]; }

	//=====pm-corner======================
	Cell &c3  = GetCell( GridX,1,k);
	Cell &co3 = GetCell( GridX+1,0,k);

	for (int i=0; i<SOU_DIM+BEA_DIM*b; i++)
	{ c3.W_Source[i] += co3.W_Source[i]; }

	//=====pp-corner======================
	Cell &c4  = GetCell( GridX,GridY,k);
	Cell &co4 = GetCell( GridX+1,GridY+1,k);

	for (int i=0; i<SOU_DIM+BEA_DIM*b; i++)
	{	c4.W_Source[i] += co4.W_Source[i]; }


	return;

}



//===============================
//======== djx/dx+djy/dy ========
//===============================
double Mesh::Dive_J(int i, int j, int k)
{


	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);
	return (cxp.W_Jx-cxm.W_Jx + cxp.B_Jx-cxm.B_Jx)*0.5/dx
		  +(cyp.W_Jy-cym.W_Jy + cyp.B_Jy-cym.B_Jy)*0.5/dy;
	
}



//===============================
//======== djy/dx-djx/dy ========
//===============================
double Mesh::Curl_J(int i, int j, int k)
{

	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	return (cxp.W_Jy-cxm.W_Jy + cxp.B_Jy-cxm.B_Jy)*0.5/dx
		  -(cyp.W_Jx-cym.W_Jx + cyp.B_Jx-cym.B_Jx)*0.5/dy;
}



void Mesh::Put_Jz(int k)
{

	double Asq;
	int i,j;

	for (j=0; j<=GridY+1; j++)
	{
		for (i=0; i<=GridX+1; i++)
		{
			Cell &ccc = GetCell(i, j, k);
			Asq  = ccc.W_Asq;
			ccc.W_Jz = 0.5*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn*( (1.0+0.5*Asq)/(1+ccc.W_Psi)/(1+ccc.W_Psi) -1 ) );

		}
	}

	return;
}





//=============================================
//======== Source For Magnetic Field: Sx ======
//=============================================
double Mesh::SourceX(int i, int j, int k)
{
	double gamma, ax, sx, Asq;
	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);


	Asq  = ccc.W_Asq;
	gamma = 0.5*(1+ccc.W_Psi)*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*Asq)/(1+ccc.W_Psi);

	ax = ((gamma*ccc.W_Ex-0.25*ccc.W_Ponx*ccc.W_Denn)/(1+ccc.W_Psi)-ccc.W_Jx*ccc.W_Ez
		-(ccc.W_Jxx*ccc.W_Ex+ccc.W_Jxy*ccc.W_Ey))/(1+ccc.W_Psi);

	sx = -ccc.W_Jy*ccc.W_Bz/(1+ccc.W_Psi) + ax - (cxp.W_Jxx-cxm.W_Jxx)*0.5/dx
		 -(cyp.W_Jxy-cym.W_Jxy)*0.5/dy + (cxp.W_Jz-cxm.W_Jz + cxp.B_Jz-cxm.B_Jz )*0.5/dx;

	if(k>0 && k< GridZ-1)
	{
		Cell &czp = GetCell(i,  j,k-1);
		Cell &czm = GetCell(i,  j,k+1);
		sx +=  (czp.B_Jx-czm.B_Jx)*0.5/dz;
	}
	return sx;
}




//=============================================
//======== Source For Magnetic Field: Sy ======
//=============================================
double Mesh::SourceY(int i, int j, int k)
{

	double gamma, ay, sy, Asq;

	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);


	Asq  = ccc.W_Asq;
	gamma = 0.5*(1+ccc.W_Psi)*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*Asq)/(1+ccc.W_Psi);

	ay = ((gamma*ccc.W_Ey-0.25*ccc.W_Pony*ccc.W_Denn)/(1+ccc.W_Psi)-ccc.W_Jy*ccc.W_Ez
		-(ccc.W_Jyy*ccc.W_Ey+ccc.W_Jxy*ccc.W_Ex))/(1+ccc.W_Psi);

	sy = ccc.W_Jx*ccc.W_Bz/(1+ccc.W_Psi) + ay - (cyp.W_Jyy-cym.W_Jyy)*0.5/dy
		 -(cxp.W_Jxy-cxm.W_Jxy)*0.5/dx + ( cyp.W_Jz-cym.W_Jz + cyp.B_Jz-cym.B_Jz )*0.5/dy;

	if(k>0 && k< GridZ-1)
	{
		Cell &czp = GetCell(i,  j,k-1);
		Cell &czm = GetCell(i,  j,k+1);
		sy +=  (czp.B_Jy-czm.B_Jy)*0.5/dz;
	}
	return sy;
}




//======================================
//======== Put Chi: n*/(1+Psi) =========
//======================================
void Mesh::Put_Chi(int k)
{

	int i,j;

	for (j=1; j<=GridY; j++)
	{
		for (i=1; i<=GridX; i++)
		{
			Cell &cc = GetCell(i,j,k);
			cc.W_Chi = cc.W_Denn/(1.0+cc.W_Psi) + cc.B_Chi;

		}

	}


	return;
}


//======================================
//======== Put dPsi/dx into W_Ex =======
//======== Put dPsi/dy into W_Ey =======
//======================================

void Mesh::Partial_Psi(int k)
{

	int i,j;

	for (j=1; j<=GridY; j++)
	{
		for (i=1; i<=GridX; i++)
		{
			Cell &ccc = GetCell(i,  j,k);
			Cell &cxm = GetCell(i-1,j,k);
			Cell &cxp = GetCell(i+1,j,k);
			Cell &cym = GetCell(i,j-1,k);
			Cell &cyp = GetCell(i,j+1,k);

			ccc.W_Ex = (cxp.W_Psi - cxm.W_Psi)*0.5/dx;
			ccc.W_Ey = (cyp.W_Psi - cym.W_Psi)*0.5/dy;

		}

	}

	return;
}



//=============================================
//======== Put d|A|^2/dx into W_Ponx ==========
//======== Put d|A|^2/dy into W_Pony ==========
//=============================================

void Mesh::Pondermotive(int k)
{

	int i,j, NF;
	double Asq_xm, Asq_ym, Asq_xp, Asq_yp;
	int NFreqs = p_domain()->NFreqs;


	for (j=1; j<=GridY; j++)
	{
		for (i=1; i<=GridX; i++)
		{
			Cell &ccc = GetCell(i,  j,k);
			Cell &cxm = GetCell(i-1,j,k);
			Cell &cxp = GetCell(i+1,j,k);
			Cell &cym = GetCell(i,j-1,k);
			Cell &cyp = GetCell(i,j+1,k);
			Asq_xm = Asq_ym = Asq_xp = Asq_yp = 0.0;
			for (NF=0; NF<NFreqs; NF++)
			{
				Asq_xm += abs(cxm.Acomx[NF])*abs(cxm.Acomx[NF]) + abs(cxm.Acomy[NF])*abs(cxm.Acomy[NF]);
				Asq_xp += abs(cxp.Acomx[NF])*abs(cxp.Acomx[NF]) + abs(cxp.Acomy[NF])*abs(cxp.Acomy[NF]);

				Asq_ym += abs(cym.Acomx[NF])*abs(cym.Acomx[NF]) + abs(cym.Acomy[NF])*abs(cym.Acomy[NF]);
				Asq_yp += abs(cyp.Acomx[NF])*abs(cyp.Acomx[NF]) + abs(cyp.Acomy[NF])*abs(cyp.Acomy[NF]);
			}

			ccc.W_Ponx = (Asq_xp - Asq_xm)*0.5/dx;
			ccc.W_Pony = (Asq_yp - Asq_ym)*0.5/dy;

		
		}

	}

	for (j=0; j<=GridY+1; j++)
	{
		for (i=0; i<=GridX+1; i++)
		{
			Cell &c = GetCell(i, j, k);
			c.W_Asq = 0.0;
			for (NF=0; NF<NFreqs; NF++)
			{
 				c.W_Asq += abs(c.Acomx[NF])*abs(c.Acomx[NF]) + abs(c.Acomy[NF])*abs(c.Acomy[NF]);
			}		
		}

	}


			

	return;
}

void Mesh::AdjustFields(int k)
{
	int i;
	int Xpa = p_domain()->p_Partition()->GetXpart();
	int Ypa = p_domain()->p_Partition()->GetYpart();

	if(RankIdx_X != 1 & RankIdx_Y != 1)
	{
	//=====mm-corner======================
		Cell &c  = GetCell( 0,0,k);
		Cell &c1 = GetCell( 1,0,k);
		Cell &c2 = GetCell( 0,1,k);
		Cell &c3 = GetCell( 1,1,k);
		for (i=0; i<WAK_DIM; i++)
		{ c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;

	}


	if(RankIdx_X != 1 & RankIdx_Y != Ypa)
	{ 	
	//=====mp-corner======================
		Cell c  = GetCell( 0,GridY+1,k);
		Cell c1 = GetCell( 1,GridY+1,k);
		Cell c2 = GetCell( 0,GridY,k);
		Cell c3 = GetCell( 1,GridY,k);
		for (i=0; i<WAK_DIM; i++)
		{ c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}


	//=====pm-corner======================
	if(RankIdx_X !=Xpa & RankIdx_Y != 1)
	{ 
		Cell c  = GetCell( GridX+1,0,k);
		Cell c1 = GetCell( GridX,  0,k);
		Cell c2 = GetCell( GridX+1,1,k);
		Cell c3 = GetCell( GridX,   1,k);
		for (i=0; i<WAK_DIM; i++)
		{ c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}

	if(RankIdx_X !=Xpa & RankIdx_Y != Ypa)
	{
	//=====pp-corner======================
		Cell c  = GetCell( GridX+1,GridY+1,k);
		Cell c1 = GetCell( GridX,  GridY+1,k);
		Cell c2 = GetCell( GridX+1,GridY,k);
		Cell c3 = GetCell( GridX,  GridY,k);
		for (i=0; i<WAK_DIM; i++)
		{ c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}
	return;

}



dcomplex Mesh::SourceAx(int i, int j, int k, int NF)
{
	dcomplex laplace;
	double k0 = p_domain()->OmegaL[NF];
	int n = i+(GridX+2)*j;
	int na= (GridX+2)*(GridY+2)*NF;
			
	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	laplace = (cxp.Acomx[NF]-2*ccc.Acomx[NF]+cxm.Acomx[NF])/dx/dx
			 +(cyp.Acomx[NF]-2*ccc.Acomx[NF]+cym.Acomx[NF])/dy/dy;
	
	return (ccc.W_Chi+4*ci*k0/dt-4*dA0)*ccc.Acomx[NF]-laplace+4*dAx1[n+na]+4*dAx2[n+na];

}


dcomplex Mesh::SourceAy(int i, int j, int k, int NF)
{
	dcomplex laplace;
	double k0 = p_domain()->OmegaL[NF];
	int n = i+(GridX+2)*j;
	int na= (GridX+2)*(GridY+2)*NF;
			
	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	laplace = (cxp.Acomy[NF]-2*ccc.Acomy[NF]+cxm.Acomy[NF])/dx/dx
			 +(cyp.Acomy[NF]-2*ccc.Acomy[NF]+cym.Acomy[NF])/dy/dy;
	
	return (ccc.W_Chi+4*ci*k0/dt-4*dA0)*ccc.Acomy[NF]-laplace+4*dAy1[n+na]+4*dAy2[n+na];
}

void Mesh::Put_dA12(int what, int k, int NF)
{
	int i, j, n;
	int na= (GridX+2)*(GridY+2)*NF;
	switch(what)
	{

		case 5:
		for (j=1; j<=GridY; j++)
		{
			for (i=1; i<=GridX; i++)
			{
				Cell &ccc = GetCell(i,j,k);
				n = i+(GridX+2)*j;
				dAx2[n+na] = dAx1[n+na]*(-0.25);
				dAx1[n+na] = -(ccc.Acomx[NF]-ccc.Acomxm[NF])*2/dz/dt;
			}
		}
		break;

		case 6:
		for (j=1; j<=GridY; j++)
		{
			for (i=1; i<=GridX; i++)
			{
				Cell &ccc = GetCell(i,j,k);
				n = i+(GridX+2)*j;
				dAy2[n+na] = dAy1[n+na]*(-0.25);
				dAy1[n+na] = -(ccc.Acomy[NF]-ccc.Acomym[NF])*2/dz/dt;
			}
		}
		break;

	}
	return;
}


void Mesh::SetFieldZeroAfter(int k0)
{

	int i,j,k,n;
	for (k=k0; k<GridZ; k++)
	{
		for (j=0; j<=GridY+1; j++)
		{
			for (i=0; i<=GridX+1; i++)
			{
				Cell &ccc = GetCell(i, j, k);
				for (n=0; n<=7; n++)
				{
					ccc.W_Fields[n]=0.0;

				}
					ccc.W_Denn=0.0;
					ccc.W_Jxx=0.0;
					ccc.W_Jyy=0.0;
			}
		}
	}

return;
}

// only for very rare situations you need set Psi-limit
// most of the divergence can be solved by change adaptive step.
void Mesh::AdjustPsi(int k)
{
	int i,j;
	double psimax=0.8/dzz; 
	double psimin=-psimax/(1.+psimax);
	for (j=0; j<=GridY+1; j++)
	{
		for (i=0; i<=GridX+1; i++)
		{
			Cell &ccc = GetCell(i, j, k);
		//	if(ccc.W_Psi<psimin) ccc.W_Psi=psimin;
		}
	}


return;

}




