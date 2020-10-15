//----------------------------------------------------------------------------------||
//-------------------                 pushparts.cpp              -------------------||
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
//---Starting---------           : Apr-29-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||


#include "wand_PIC.h"


void Mesh::PushParticle()
{


	Particle *p = NULL;

	double ddx;
	double ddy;
	double ddz;
	double ddz2;

	double dt  = p_domain()->Get_dt();
	int Nt  = p_domain()->Get_SubCycle();
	double dtt = dt/Nt;

	double Ex, Ey, Ez, Psi, Bx, By, Bz, Pondx, Pondy, Asq;

	dcomplex Exl, Eyl, Ezl, Bxl, Byl, Bzl;

	double ExlR, EylR, EzlR, BxlR, BylR, BzlR;

	double x, y, z, x0, y0, z0, px, py, pz, px0, py0, pz0, gamma;

	double pxp, pyp, pzp;

	double q2m;
	double wmmm, wmpm, wpmm, wppm;
	double wmmp, wmpp, wpmp, wppp;
	double cc0,  cc1,  cc2;

	int i, j, k, k2, n;

	int type;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;

	int NFreqs =  p_domain()->NFreqs;
	int NF;
	double OmegaL;


	//===radiation part=========
	double reinkp=sqrt(3)*(1e-9)*sqrt(Pla_ne/1.40251e20);//  sqrt(3)*re/lambda_L; classical electron radius in lambda_p
	double theta_x;
	double phi_y;
	double Rcurve;   // radius of curvature.
	double omegac;   //critical frequency.

 	double ThetaMax = XRayDetector->ThetaMax;
 	double PhiMax   = XRayDetector->PhiMax;
	double OmegaMax = XRayDetector->OmegaMax;
	double OmegaMin = XRayDetector->OmegaMin;

	int NOmega = XRayDetector->NOmega;
	int NTheta = XRayDetector->NTheta;
	int NPhi   = XRayDetector->NPhi;

	double dtheta =2*ThetaMax/NTheta;
	double dphi   =2*PhiMax/NPhi;
	double domega =(OmegaMax-OmegaMin)/NOmega;  //in KeV


	p = p_Particle;

	while (p)
	{

		q2m  = p-> q2m;
		type = p-> type;

		x  = p->  x;
		y  = p->  y;
		z  = p->  z;

		px = p-> px;
		py = p-> py;
		pz = p-> pz;

		gamma = sqrt(1 + px*px + py*py + pz*pz);

		for (n=0; n<Nt; n++)
		{
			//=======================================
			//= Step-(1) half-position push in time =
			x0=x;
			y0=y;
			z0=z;
			px0=px;
			py0=py;
			pz0=pz;

			x += 0.5*dtt*   px/gamma;
			y += 0.5*dtt*   py/gamma;
			z += 0.5*dtt*(1-pz/gamma);
			//=======================================


			//=============================//
			if(RankIdx_X == 1	&& x<=Offset_X) break;
			if(RankIdx_X == Xpa && x>=Xmax)	    break; 
			if(RankIdx_Y == 1	&& y<=Offset_Y) break;
			if(RankIdx_Y == Ypa && y>=Ymax)     break;
			//=============================//


			ddx = x-(Offset_X-dx*0.5);
			ddy = y-(Offset_Y-dy*0.5);
			ddz = z;
			ddz2= z-0.5*dz;

			i = floor(ddx/dx);
			j = floor(ddy/dy);
			k = floor(ddz/dz);
			k2= floor(ddz2/dz);

			//=============================//
			if(i<0) i = 0;
			if(i>GridX) i = GridX;
			if(j<0) j = 0;
			if(j>GridY) j = GridY;
			if(k<0 || k > GridZ-2 || k2<0 || k2>GridZ-2) break;
			//=============================//

			//==============================
			//=====wakefields===============
			wmmm = (i+1-ddx/dx) *(j+1-ddy/dy) *(k+1-ddz/dz);
			wmpm = (i+1-ddx/dx) *(ddy/dy-j)   *(k+1-ddz/dz);
			wpmm = (ddx/dx-i)   *(j+1-ddy/dy) *(k+1-ddz/dz);
			wppm = (ddx/dx-i)   *(ddy/dy-j)   *(k+1-ddz/dz);

			wmmp = (i+1-ddx/dx) *(j+1-ddy/dy) *(ddz/dz-k);
			wmpp = (i+1-ddx/dx) *(ddy/dy-j)   *(ddz/dz-k);
			wpmp = (ddx/dx-i)   *(j+1-ddy/dy) *(ddz/dz-k);
			wppp = (ddx/dx-i)   *(ddy/dy-j)   *(ddz/dz-k);

			Cell &cmmm = GetCell(i,   j,   k);
			Cell &cmpm = GetCell(i,   j+1, k);
			Cell &cpmm = GetCell(i+1, j,   k);
			Cell &cppm = GetCell(i+1, j+1, k);

			Cell &cmmp = GetCell(i,   j,   k+1);
			Cell &cmpp = GetCell(i,   j+1, k+1);
			Cell &cpmp = GetCell(i+1, j,   k+1);
			Cell &cppp = GetCell(i+1, j+1, k+1);


			Bx = wmmm*cmmm.W_Bx  + wmpm*cmpm.W_Bx  + wpmm*cpmm.W_Bx  + wppm*cppm.W_Bx
			   + wmmp*cmmp.W_Bx  + wmpp*cmpp.W_Bx  + wpmp*cpmp.W_Bx  + wppp*cppp.W_Bx;
			
			By = wmmm*cmmm.W_By  + wmpm*cmpm.W_By  + wpmm*cpmm.W_By  + wppm*cppm.W_By
			   + wmmp*cmmp.W_By  + wmpp*cmpp.W_By  + wpmp*cpmp.W_By  + wppp*cppp.W_By;
			
			Bz = wmmm*cmmm.W_Bz  + wmpm*cmpm.W_Bz  + wpmm*cpmm.W_Bz  + wppm*cppm.W_Bz
			   + wmmp*cmmp.W_Bz  + wmpp*cmpp.W_Bz  + wpmp*cpmp.W_Bz  + wppp*cppp.W_Bz;
			

			// Ex= - dpsi/dy + By
			Ex =-(wmmm*cmmm.W_Ex  + wmpm*cmpm.W_Ex  + wpmm*cpmm.W_Ex  + wppm*cppm.W_Ex
			    + wmmp*cmmp.W_Ex  + wmpp*cmpp.W_Ex  + wpmp*cpmp.W_Ex  + wppp*cppp.W_Ex) + By;
			
			// Ey= - dpsi/dy - Bx
			Ey =-(wmmm*cmmm.W_Ey  + wmpm*cmpm.W_Ey  + wpmm*cpmm.W_Ey  + wppm*cppm.W_Ey
			    + wmmp*cmmp.W_Ey  + wmpp*cmpp.W_Ey  + wpmp*cpmp.W_Ey  + wppp*cppp.W_Ey) - Bx;
			

			Ez = wmmm*cmmm.W_Ez  + wmpm*cmpm.W_Ez  + wpmm*cpmm.W_Ez  + wppm*cppm.W_Ez
			   + wmmp*cmmp.W_Ez  + wmpp*cmpp.W_Ez  + wpmp*cpmp.W_Ez  + wppp*cppp.W_Ez;
			//=====wakefields===============
			//==============================



			//===============================
			//=====laserfields===============
			ExlR = EylR = EzlR = 0.0;
			BxlR = BylR = BzlR = 0.0;
			if(NFreqs>0)
			{

			wmmm = (i+1-ddx/dx) *(j+1-ddy/dy) *(k2+1-ddz2/dz);
			wmpm = (i+1-ddx/dx) *(ddy/dy-j)   *(k2+1-ddz2/dz);
			wpmm = (ddx/dx-i)   *(j+1-ddy/dy) *(k2+1-ddz2/dz);
			wppm = (ddx/dx-i)   *(ddy/dy-j)   *(k2+1-ddz2/dz);

			wmmp = (i+1-ddx/dx) *(j+1-ddy/dy) *(ddz2/dz-k2);
			wmpp = (i+1-ddx/dx) *(ddy/dy-j)   *(ddz2/dz-k2);
			wpmp = (ddx/dx-i)   *(j+1-ddy/dy) *(ddz2/dz-k2);
			wppp = (ddx/dx-i)   *(ddy/dy-j)   *(ddz2/dz-k2);

			Cell &cmmm2 = GetCell(i,   j,   k2);
			Cell &cmpm2 = GetCell(i,   j+1, k2);
			Cell &cpmm2 = GetCell(i+1, j,   k2);
			Cell &cppm2 = GetCell(i+1, j+1, k2);

			Cell &cmmp2 = GetCell(i,   j,   k2+1);
			Cell &cmpp2 = GetCell(i,   j+1, k2+1);
			Cell &cpmp2 = GetCell(i+1, j,   k2+1);
			Cell &cppp2 = GetCell(i+1, j+1, k2+1);

			for(NF=0; NF<NFreqs; NF++)
			{
				
				Exl = wmmm*cmmm2.L_Ex[NF] + wmpm*cmpm2.L_Ex[NF] + wpmm*cpmm2.L_Ex[NF] + wppm*cppm2.L_Ex[NF]
					+ wmmp*cmmp2.L_Ex[NF] + wmpp*cmpp2.L_Ex[NF] + wpmp*cpmp2.L_Ex[NF] + wppp*cppp2.L_Ex[NF];

				Eyl = wmmm*cmmm2.L_Ey[NF] + wmpm*cmpm2.L_Ey[NF] + wpmm*cpmm2.L_Ey[NF] + wppm*cppm2.L_Ey[NF]
					+ wmmp*cmmp2.L_Ey[NF] + wmpp*cmpp2.L_Ey[NF] + wpmp*cpmp2.L_Ey[NF] + wppp*cppp2.L_Ey[NF];

				Ezl = wmmm*cmmm2.L_Ez[NF] + wmpm*cmpm2.L_Ez[NF] + wpmm*cpmm2.L_Ez[NF] + wppm*cppm2.L_Ez[NF]
					+ wmmp*cmmp2.L_Ez[NF] + wmpp*cmpp2.L_Ez[NF] + wpmp*cpmp2.L_Ez[NF] + wppp*cppp2.L_Ez[NF];

				Bxl = wmmm*cmmm2.L_Bx[NF] + wmpm*cmpm2.L_Bx[NF] + wpmm*cpmm2.L_Bx[NF] + wppm*cppm2.L_Bx[NF]
					+ wmmp*cmmp2.L_Bx[NF] + wmpp*cmpp2.L_Bx[NF] + wpmp*cpmp2.L_Bx[NF] + wppp*cppp2.L_Bx[NF];

				Byl = wmmm*cmmm2.L_By[NF] + wmpm*cmpm2.L_By[NF] + wpmm*cpmm2.L_By[NF] + wppm*cppm2.L_By[NF]
					+ wmmp*cmmp2.L_By[NF] + wmpp*cmpp2.L_By[NF] + wpmp*cpmp2.L_By[NF] + wppp*cppp2.L_By[NF];

				OmegaL = p_domain()->OmegaL[NF];

				ExlR += (Exl*exp(-ci*OmegaL*z)).real();
				EylR += (Eyl*exp(-ci*OmegaL*z)).real();
				EzlR += (Ezl*exp(-ci*OmegaL*z)).real();
				BxlR += (Bxl*exp(-ci*OmegaL*z)).real();
				BylR += (Byl*exp(-ci*OmegaL*z)).real();	
				BzlR  =0.0;
			}

			}
			//=====laserfields===============
			//===============================


			//===========All fields==============
			ExlR *=(-q2m); EylR *=(-q2m); EzlR *=(-q2m);
			BxlR *=(-q2m); BylR *=(-q2m); BzlR *=(-q2m);

			Ex   *=(-q2m); Ey   *=(-q2m); Ez   *=(-q2m);
			Bx   *=(-q2m); By   *=(-q2m); Bz   *=(-q2m);

			Bx +=BxlR; By +=BylR; Bz +=BzlR;


			//=======================================
			//= Step-(2) half-Electric push in time =
			pxp = px + 0.5*dtt*(Ex+ExlR);
			pyp = py + 0.5*dtt*(Ey+EylR);
			pzp = pz + 0.5*dtt*(Ez+EzlR);
			//=======================================


			//=======================================
			//= Step-(3) Full-Magnetic push in time =
			cc0 = dtt*0.5/sqrt(1+pxp*pxp+pyp*pyp+pzp*pzp); 
			cc1 = 2*cc0/(1+cc0*cc0*(Bx*Bx+By*By+Bz*Bz));
			cc2 = cc0*cc1;

			px = (1-(By*By+Bz*Bz)*cc2)*pxp + ( Bz*cc1 + Bx*By*cc2)*pyp + (-By*cc1 + Bx*Bz*cc2)*pzp;
			py = (-Bz*cc1 + Bx*By*cc2)*pxp + (1-(Bx*Bx+Bz*Bz)*cc2)*pyp + ( Bx*cc1 + By*Bz*cc2)*pzp;
			pz = ( By*cc1 + Bx*Bz*cc2)*pxp + (-Bx*cc1 + By*Bz*cc2)*pyp + (1-(Bx*Bx+By*By)*cc2)*pzp;
			//=======================================

			//=======================================
			//= Step-(4) half-Electric push in time =
			px = px + 0.5*dtt*(Ex+ExlR);
			py = py + 0.5*dtt*(Ey+EylR);
			pz = pz + 0.5*dtt*(Ez+EzlR);
			//=======================================


			//=======================================
			//= Step-(5) half-position push in time =
			gamma = sqrt(1 + px*px + py*py + pz*pz);

			x += 0.5*dtt*   px/gamma;
			y += 0.5*dtt*   py/gamma;
			z += 0.5*dtt*(1-pz/gamma);

			p->Wxw += 	   (x-x0) *Ex;
			p->Wyw +=      (y-y0) *Ey;
			p->Wzw += (dtt-(z-z0))*Ez;

			p->Wxl +=     ( x-x0) *ExlR;
			p->Wyl += 	   (y-y0) *EylR;
			p->Wzl += (dtt-(z-z0))*EzlR;





			//======================================	
			// Raidation Calculation
			//
			if(XRayDetector->IfRadiation==1 && gamma>2 && pz>0)
			{

				theta_x = px/pz;
				phi_y   = py/pz;
				double vs =(px*px + py*py + pz*pz)/gamma/gamma;
				double vvs=( (px-px0)*(px-px0) + (py-py0)*(py-py0) + (pz-pz0)*(pz-pz0) )/gamma/gamma/dtt/dtt;
				double vdvv= ( (px-px0)*px+ (py-py0)*py + (pz-pz0)*pz )/gamma/gamma/dtt;
				if(vs*vvs-vdvv*vdvv==0)
				{
					Rcurve  = 1e20;
					omegac  = 0.0;
					
				}
				else
				{
					Rcurve  = pow(vs, 1.5)/sqrt(vs*vvs-vdvv*vdvv);
					omegac  = 1.5*gamma*gamma*gamma/Rcurve;
				}
				
				// calculate wavepacket energy
				double k1=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(dtt-(z-z0))*(dtt-(z-z0)))/Rcurve*gamma*reinkp;
				k1 *=p->weight;

				omegac *= 1.2398*sqrt(Pla_ne/1.114855e21)/1000;

				// deposit the packet on detector
				int idx_i=floor(omegac/domega);
				int idx_j=floor((theta_x+ThetaMax)/dtheta-0.5);
				int idx_k=floor((phi_y+PhiMax)/dphi-0.5);
				
				if(idx_i>=0&&idx_i<NOmega
				 &&idx_j>=0&&idx_j<NTheta
				 &&idx_k>=0&&idx_k<NPhi) XRayDetector->PutPacket(idx_i, idx_j, idx_k, k1);

			}




		}

		p-> x  = x;
		p-> y  = y;
		p-> z  = z;

		p-> px = px;
		p-> py = py;
		p-> pz = pz;

		p-> gamma = gamma;

		p = p->p_PrevPart;

	}

	ExchangeP();

	return;
}



Particle* Mesh::Reconnect(Particle* p_Part)
{

	Particle* p_temp;
	if(p_Part->p_PrevPart)
	{
		p_Part->p_PrevPart->p_NextPart = p_Part->p_NextPart;
		p_temp = p_Part->p_PrevPart;
	}
	else
	{
		p_temp = NULL;
	}

	if(p_Part->p_NextPart)
	{
		p_Part->p_NextPart->p_PrevPart = p_Part->p_PrevPart;
	}
	else
	{
		p_Particle = p_Part->p_PrevPart;

	}

	delete p_Part;
	return p_temp;

}




void Mesh::ExchangeP()
{
	Particle *p = NULL;
	
	//=========Send and Receive Buf Size===============
	double bufsize = p_domain()->p_Com()->Get_bufsize();
	bufsize *= (GridX*SOU_DIM*2.0/SDP_DIM);
	//=================================================

	double xp, yp, zp;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;
	double Zmax = (GridZ-2)*dz;

	double Xsize = GridX*dx*Xpa;
	double Ysize = GridY*dy*Ypa;

	int Sendxm, Sendxp, Sendym, Sendyp;
	int S_SUM, A_SUM;

	while(1)
	{

		Sendxm = Sendxp = Sendym = Sendyp = 0; 
		p = p_Particle;
		while (p)
		{
			xp = p-> x;
			yp = p-> y;
			zp = p-> z;

		//=================================================
		//============Partectory Outside Boundary =========
		//=================================================
		if(p_domain()->Get_BC()==1)
		{
			
		if(zp<0 || zp>Zmax)				     { p = Reconnect(p); continue;}
		if(RankIdx_X ==1	&& xp<=Offset_X) { p = Reconnect(p); continue;}
		if(RankIdx_X == Xpa && xp>=Xmax)     { p = Reconnect(p); continue;}
		if(RankIdx_Y == 1	&& yp<=Offset_Y) { p = Reconnect(p); continue;}
		if(RankIdx_Y == Ypa && yp>=Ymax)     { p = Reconnect(p); continue;}

		}

		//==================================================

		//====================================
		//====== Send to left Neighbor =======
		//====================================
		if(xp < Offset_X)
		{
			if(RankIdx_X ==1) 	p-> x = xp + Xsize;
			Sendxm +=1;
			PackP(p, Sendxm, 0);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to right Neighbor ======
		//====================================
		else if(xp > Xmax)
		{
			if(RankIdx_X ==Xpa) p-> x = xp - Xsize;
			Sendxp +=1;
			PackP(p, Sendxp, 1);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to Neighbor Below ======
		//====================================
		else if(yp < Offset_Y)
		{
			if(RankIdx_Y ==1) 	p-> y = yp + Ysize;
			Sendym +=1;
			PackP(p, Sendym, 2);
			p = Reconnect(p);

		}
		//====================================
		//====== Send to Neighbor Above ======
		//====================================
		else if(yp > Ymax)
		{
			if(RankIdx_Y ==Ypa) p-> y = yp - Ysize;
			Sendyp +=1;
			PackP(p, Sendyp, 3);
			p = Reconnect(p);
		}
		else
		{
			p = p->p_PrevPart;
		}

		}


		//=====================================================================
		//======== Exchange Partectoryies with Neighboring Processors =========
		//=====================================================================
		int maxsend=Sendxm;
		if(bufsize<Sendxm || bufsize<Sendxp || bufsize<Sendym || bufsize<Sendyp)
		{	
			if(Sendxp>maxsend) maxsend=Sendxp;
			if(Sendyp>maxsend) maxsend=Sendyp;
			if(Sendym>maxsend) maxsend=Sendym;

			printf("==== Mesh: At Rank: %5d. ==================\n",Rank);
			std::cout << "==== Mesh: Send Too Many Particles.      ====\n";
			std::cout << "==== Mesh: May Cause Memory Problems.    ====\n";
			printf("==== Mesh: Try Buf Size: %3d           ====\n",maxsend);
			
		}

		S_SUM = Sendxm+Sendxp+Sendym+Sendyp;
		MPI_Allreduce(&S_SUM, &A_SUM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( A_SUM == 0 ) {break;};
		p_domain()->p_Com()->DoCommuteT(COMMU_P, Sendxm, Sendxp, Sendym, Sendyp);

	}

	return;
}



void Mesh::PackP(Particle* p_Part, int Sendn, int where)
{
	Commute *p_COMM = p_domain()->p_Com();

	switch(where)
	{

		case 0:
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 0] = p_Part-> x;   	
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 1] = p_Part-> y;   	
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 2] = p_Part-> z;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 3] = p_Part-> x0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 4] = p_Part-> y0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 5] = p_Part-> z0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 6] = p_Part-> px;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 7] = p_Part-> py;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 8] = p_Part-> pz;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM + 9] = p_Part-> Ex0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +10] = p_Part-> Ey0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +11] = p_Part-> Ez0;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +12] = p_Part-> type*1.0;	
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +13] = p_Part-> q2m;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +14] = p_Part-> weight;

		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +15] = p_Part-> Wxw;	
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +16] = p_Part-> Wyw;	
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +17] = p_Part-> Wzw;
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +18] = p_Part-> Wxl;
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +19] = p_Part-> Wyl;		
		p_COMM->SendSourceXm[(Sendn-1)*SDP_DIM +20] = p_Part-> Wzl;

		break;
	
		case 1:
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 0] = p_Part-> x;   	
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 1] = p_Part-> y;   	
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 2] = p_Part-> z;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 3] = p_Part-> x0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 4] = p_Part-> y0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 5] = p_Part-> z0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 6] = p_Part-> px;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 7] = p_Part-> py;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 8] = p_Part-> pz;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM + 9] = p_Part-> Ex0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +10] = p_Part-> Ey0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +11] = p_Part-> Ez0;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +12] = p_Part-> type*1.0;	
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +13] = p_Part-> q2m;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +14] = p_Part-> weight;

		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +15] = p_Part-> Wxw;	
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +16] = p_Part-> Wyw;	
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +17] = p_Part-> Wzw;
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +18] = p_Part-> Wxl;
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +19] = p_Part-> Wyl;		
		p_COMM->SendSourceXp[(Sendn-1)*SDP_DIM +20] = p_Part-> Wzl;


		break;

		case 2:
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 0] = p_Part-> x;   	
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 1] = p_Part-> y;   	
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 2] = p_Part-> z;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 3] = p_Part-> x0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 4] = p_Part-> y0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 5] = p_Part-> z0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 6] = p_Part-> px;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 7] = p_Part-> py;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 8] = p_Part-> pz;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM + 9] = p_Part-> Ex0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +10] = p_Part-> Ey0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +11] = p_Part-> Ez0;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +12] = p_Part-> type*1.0;	
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +13] = p_Part-> q2m;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +14] = p_Part-> weight;

		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +15] = p_Part-> Wxw;	
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +16] = p_Part-> Wyw;	
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +17] = p_Part-> Wzw;
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +18] = p_Part-> Wxl;
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +19] = p_Part-> Wyl;		
		p_COMM->SendSourceYm[(Sendn-1)*SDP_DIM +20] = p_Part-> Wzl;


		break;

		case 3:
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 0] = p_Part-> x;   	
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 1] = p_Part-> y;   	
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 2] = p_Part-> z;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 3] = p_Part-> x0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 4] = p_Part-> y0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 5] = p_Part-> z0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 6] = p_Part-> px;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 7] = p_Part-> py;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 8] = p_Part-> pz;
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM + 9] = p_Part-> Ex0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +10] = p_Part-> Ey0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +11] = p_Part-> Ez0;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +12] = p_Part-> type*1.0;	
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +13] = p_Part-> q2m;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +14] = p_Part-> weight;

		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +15] = p_Part-> Wxw;	
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +16] = p_Part-> Wyw;	
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +17] = p_Part-> Wzw;
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +18] = p_Part-> Wxl;
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +19] = p_Part-> Wyl;		
		p_COMM->SendSourceYp[(Sendn-1)*SDP_DIM +20] = p_Part-> Wzl;

		break;
	}

	return;
}