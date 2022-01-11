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


void Mesh::PushParticle()//v
{


	Particle *p = NULL;

	WDOUBLE ddx;
	WDOUBLE ddy;
	WDOUBLE ddz;
	WDOUBLE ddz2;
	WDOUBLE sx,sy,sxy;

	WDOUBLE dt  = p_domain()->Get_dt();
	int Nt  = p_domain()->Get_SubCycle();
	WDOUBLE dtt = dt/Nt;

	WDOUBLE Ex, Ey, Ez, Psi, Bx, By, Bz, Pondx, Pondy, Asq;

	dcomplex Exl, Eyl, Ezl, Bxl, Byl, Bzl;

	WDOUBLE ExlR, EylR, EzlR, BxlR, BylR, BzlR;

	WDOUBLE x, y, z, x0, y0, z0, px, py, pz, px0, py0, pz0, gamma;

	WDOUBLE pxp, pyp, pzp;

	WDOUBLE q2m;

	WDOUBLE cc0,  cc1,  cc2;

	int i, j, k, k2, n;

	int type;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	WDOUBLE Xmax = Offset_X+GridX*dx;
	WDOUBLE Ymax = Offset_Y+GridY*dy;

	Cell  **c = new Cell*[18];
	WDOUBLE weight[18];

	int NFreqs =  p_domain()->NFreqs;
	int NF;
	WDOUBLE OmegaL;


	//===radiation part=========
	WDOUBLE reinkp=sqrt(3)*(1e-9)*sqrt(Pla_ne/1.40251e20);//  sqrt(3)*re/lambda_L; classical electron radius in lambda_p
	WDOUBLE theta_x;
	WDOUBLE phi_y;
	WDOUBLE Rcurve;   // radius of curvature.
	WDOUBLE omegac;   //critical frequency.

 	WDOUBLE ThetaMax = XRayDetector->ThetaMax;
 	WDOUBLE PhiMax   = XRayDetector->PhiMax;
	WDOUBLE OmegaMax = XRayDetector->OmegaMax;
	WDOUBLE OmegaMin = XRayDetector->OmegaMin;

	int NOmega = XRayDetector->NOmega;
	int NTheta = XRayDetector->NTheta;
	int NPhi   = XRayDetector->NPhi;

	WDOUBLE dtheta =2*ThetaMax/NTheta;
	WDOUBLE dphi   =2*PhiMax/NPhi;
	WDOUBLE domega =(OmegaMax-OmegaMin)/NOmega;  //in KeV

	WDOUBLE minga=1e20;


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

			i=p->idx_i;
			j=p->idx_j;

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

			Cell *ctmp = &GetCell(i,j,0);
			while(x>(ctmp->Xcord+ctmp->dx*0.5) && i<GridX) { p->idx_i++;  i++; ctmp=&GetCell(i,j,0); }
			while(x<(ctmp->Xcord-ctmp->dx*0.5) && i>1) 	   { p->idx_i--;  i--; ctmp=&GetCell(i,j,0); }
			while(y>(ctmp->Ycord+ctmp->dy*0.5) && j<GridY) { p->idx_j++;  j++; ctmp=&GetCell(i,j,0); }
			while(y<(ctmp->Ycord-ctmp->dy*0.5) && j>1) 	   { p->idx_j--;  j--; ctmp=&GetCell(i,j,0); }

			//=======================================

			//=============================//
			if(RankIdx_X == 1	&& i==0) 		break;
			if(RankIdx_X == Xpa && i==GridX+1)	break; 
			if(RankIdx_Y == 1	&& j==0) 		break;
			if(RankIdx_Y == Ypa && j==GridY+1)  break;
			//=============================//

			ddz = z;
			ddz2= z-0.5*dz;

			k = floor(ddz/dz);
			k2= floor(ddz2/dz);
			if(k<0 || k > GridZ-2 || k2<0 || k2>GridZ-2) break;
	

			//==============================
			//=====wakefields===============
			c[0] = &GetCell(i-1,j-1,k); c[9] =  &GetCell(i-1,j-1,k+1);
			c[1] = &GetCell(i-1,j  ,k); c[10] = &GetCell(i-1,j  ,k+1);
			c[2] = &GetCell(i-1,j+1,k); c[11] = &GetCell(i-1,j+1,k+1);

			c[3] = &GetCell(i,  j-1,k); c[12] = &GetCell(i,  j-1,k+1);
			c[4] = &GetCell(i,  j  ,k); c[13] = &GetCell(i,  j  ,k+1);
			c[5] = &GetCell(i,  j+1,k); c[14] = &GetCell(i,  j+1,k+1);
		
			c[6] = &GetCell(i+1,j-1,k); c[15] = &GetCell(i+1,j-1,k+1);
			c[7] = &GetCell(i+1,j  ,k); c[16] = &GetCell(i+1,j  ,k+1);
			c[8] = &GetCell(i+1,j+1,k); c[17] = &GetCell(i+1,j+1,k+1);
			Cell &ccc= *c[4];

			ddx=ccc.dx;
			ddy=ccc.dy;
			sx=ddx;  //- re-size
			sy=ddy;  //- re-size
			sxy=sx*sy;

			WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+x-ccc.Xcord),0.0);
			WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+y-ccc.Ycord),0.0);

			WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-x+ccc.Xcord),0.0);
			WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-y+ccc.Ycord),0.0);

			WDOUBLE deltaxc=sx-deltaxm-deltaxp;
			WDOUBLE deltayc=sy-deltaym-deltayp;

			WDOUBLE deltaz=(ddz/dz-k);
		
			weight[0] = deltaxm*deltaym/sxy;
			weight[1] = deltaxm*deltayc/sxy;
			weight[2] = deltaxm*deltayp/sxy;

			weight[3] = deltaxc*deltaym/sxy;
			weight[4] = deltaxc*deltayc/sxy;
			weight[5] = deltaxc*deltayp/sxy;

			weight[6] = deltaxp*deltaym/sxy;
			weight[7] = deltaxp*deltayc/sxy;
			weight[8] = deltaxp*deltayp/sxy;

			for(int n=0;n<9;n++)
			{
				weight[n+9]=weight[n]*deltaz;
				weight[n] *=(1-deltaz);
			}

			Bx=By=Bz=Ex=Ey=Ez=0.0;
			for(int n=0;n<18;n++)
			{
				Bx+=c[n]->W_Bx*weight[n];
				By+=c[n]->W_By*weight[n];
				Bz+=c[n]->W_Bz*weight[n];

				Ex-=c[n]->W_Ex*weight[n];
				Ey-=c[n]->W_Ey*weight[n];
				Ez+=c[n]->W_Ez*weight[n];
			}
			// Ex= - dpsi/dy + By
			Ex += By;
			// Ey= - dpsi/dy - Bx
			Ey -= Bx;
			
			//===============================
			//=====laserfields===============
			ExlR = EylR = EzlR = 0.0;
			BxlR = BylR = BzlR = 0.0;
			
			if(NFreqs>0)
			{

			c[0] = &GetCell(i-1,j-1,k2); c[9]  = &GetCell(i-1,j-1,k2+1);
			c[1] = &GetCell(i-1,j  ,k2); c[10] = &GetCell(i-1,j  ,k2+1);
			c[2] = &GetCell(i-1,j+1,k2); c[11] = &GetCell(i-1,j+1,k2+1);

			c[3] = &GetCell(i,  j-1,k2); c[12] = &GetCell(i,  j-1,k2+1);
			c[4] = &GetCell(i,  j  ,k2); c[13] = &GetCell(i,  j  ,k2+1);
			c[5] = &GetCell(i,  j+1,k2); c[14] = &GetCell(i,  j+1,k2+1);
		
			c[6] = &GetCell(i+1,j-1,k2); c[15] = &GetCell(i+1,j-1,k2+1);
			c[7] = &GetCell(i+1,j  ,k2); c[16] = &GetCell(i+1,j  ,k2+1);
			c[8] = &GetCell(i+1,j+1,k2); c[17] = &GetCell(i+1,j+1,k2+1);

			WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+x-ccc.Xcord),0.0);
			WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+y-ccc.Ycord),0.0);

			WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-x+ccc.Xcord),0.0);
			WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-y+ccc.Ycord),0.0);

			WDOUBLE deltaxc=sx-deltaxm-deltaxp;
			WDOUBLE deltayc=sy-deltaym-deltayp;

			WDOUBLE deltaz=(ddz2/dz-k2);
		
			weight[0] = deltaxm*deltaym/sxy;
			weight[1] = deltaxm*deltayc/sxy;
			weight[2] = deltaxm*deltayp/sxy;

			weight[3] = deltaxc*deltaym/sxy;
			weight[4] = deltaxc*deltayc/sxy;
			weight[5] = deltaxc*deltayp/sxy;

			weight[6] = deltaxp*deltaym/sxy;
			weight[7] = deltaxp*deltayc/sxy;
			weight[8] = deltaxp*deltayp/sxy;

			for(int n=0;n<9;n++)
			{
				weight[n+9]=weight[n]*deltaz;
				weight[n] *=weight[n]*(1-deltaz);
			}

			for(NF=0; NF<NFreqs; NF++)
			{
				Exl=Eyl=Ezl=Bxl=Byl=Bzl=0.0;
				for(int n=0;n<18;n++)
				{
					Exl+=c[n]->L_Ex[NF] *weight[n];
					Eyl+=c[n]->L_Ey[NF] *weight[n];
					Ezl+=c[n]->L_Ez[NF] *weight[n];
					Bxl+=c[n]->L_Bx[NF] *weight[n];
					Byl+=c[n]->L_By[NF] *weight[n];
					Bzl+=c[n]->L_Bz[NF] *weight[n];
				}

				OmegaL = p_domain()->OmegaL[NF];
				ExlR += (Exl*exp(-ci*OmegaL*z)).real();
				EylR += (Eyl*exp(-ci*OmegaL*z)).real();
				EzlR += (Ezl*exp(-ci*OmegaL*z)).real();
				BxlR += (Bxl*exp(-ci*OmegaL*z)).real();
				BylR += (Byl*exp(-ci*OmegaL*z)).real();	
				BzlR += (Bzl*exp(-ci*OmegaL*z)).real();	
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
			if(XRayDetector->IfRadiation==1 &&type==ELECTRON && gamma>2 && pz>0)
			{

				theta_x = px/pz;
				phi_y   = py/pz;
				WDOUBLE vs =(px*px + py*py + pz*pz)/gamma/gamma;
				WDOUBLE vvs=( (px-px0)*(px-px0) + (py-py0)*(py-py0) + (pz-pz0)*(pz-pz0) )/gamma/gamma/dtt/dtt;
				WDOUBLE vdvv= ( (px-px0)*px+ (py-py0)*py + (pz-pz0)*pz )/gamma/gamma/dtt;
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
				WDOUBLE k1=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(dtt-(z-z0))*(dtt-(z-z0)))/Rcurve*gamma*reinkp;
				k1 *=p->weight*p->sx*p->sy;

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

		if(pz>2&&gamma<minga) minga=gamma;

		Cell *ctmp = &GetCell(i,j,0);
		while(x>(ctmp->Xcord+ctmp->dx*0.5) && i<GridX+1) { p->idx_i++;  i++; ctmp=&GetCell(i,j,0); }
		while(x<(ctmp->Xcord-ctmp->dx*0.5) && i>0) 	     { p->idx_i--;  i--; ctmp=&GetCell(i,j,0); }
		while(y>(ctmp->Ycord+ctmp->dy*0.5) && j<GridY+1) { p->idx_j++;  j++; ctmp=&GetCell(i,j,0); }
		while(y<(ctmp->Ycord-ctmp->dy*0.5) && j>0) 	     { p->idx_j--;  j--; ctmp=&GetCell(i,j,0); }
		p = p->p_PrevPart;

	}

	minGamma=minga;

	ExchangeP();

	return;
}


void Mesh::SetNewTimeStep()
{
	WDOUBLE gmin;
	MPI_Allreduce(&minGamma, &gmin, 1, MPI_WDOUBLE, MPI_MIN, MPI_COMM_WORLD);
	dt= dt0*sqrt(gmin);
	p_domain()->set_new_dt(dt);

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
	Commute *p_COMM = p_domain()->p_Com();
	//=========Send and Receive Buf Size===============
	WDOUBLE bufsize = p_domain()->p_Com()->Get_bufsize();
	bufsize *= (GridX*SOU_DIM*2.0/SDP_DIM);
	//=================================================

	WDOUBLE xp, yp, zp;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();


	WDOUBLE Zmax = (GridZ-2)*dz;

	std::vector<int> SendN(8,0);//Sendmm, Sendmp, Sendpm, Sendpp; Sendxm, Sendxp, Sendym, Sendyp;

	int S_SUM, A_SUM;

	int i;
	int j;

	while(1)
	{

		p = p_Particle;

		SendN[0]=SendN[1]=SendN[2]=SendN[3]=SendN[4]=SendN[5]=SendN[6]=SendN[7]=0;
		
		WDOUBLE* SeXm=p_COMM->SendSourceXm;
		WDOUBLE* SeXp=p_COMM->SendSourceXp;
		WDOUBLE* SeYm=p_COMM->SendSourceYm;
		WDOUBLE* SeYp=p_COMM->SendSourceYp;

		WDOUBLE* Semm=p_COMM->SendSourcemm;
		WDOUBLE* Semp=p_COMM->SendSourcemp;
		WDOUBLE* Sepm=p_COMM->SendSourcepm;
		WDOUBLE* Sepp=p_COMM->SendSourcepp;

		while (p)
		{
			i=p->idx_i;
			j=p->idx_j;
			zp = p-> z;

			//=================================================
			//============Partectory Outside Boundary =========
			//=================================================
			if(zp<=dz || zp>=Zmax)	{ p = Reconnect(p); continue;}

			if(p_domain()->Get_BC()==1)
			{
				if(RankIdx_X ==1	&& i==0) 	   { p = Reconnect(p); continue;}
				if(RankIdx_X == Xpa && i==GridX+1) { p = Reconnect(p); continue;}
				if(RankIdx_Y == 1	&& j==0) 	   { p = Reconnect(p); continue;}
				if(RankIdx_Y == Ypa && j==GridY+1) { p = Reconnect(p); continue;}
			}

			//==================================================
			if(i==0&&j==0)
			{
				SendN[0] +=1;
				PackP(p, Semm);
				p = Reconnect(p);
				continue;
			}

			if(i==0&&j==GridY+1)
			{
				SendN[1] +=1;
				PackP(p, Semp);
				p = Reconnect(p);
				continue;
			}
			if(i==GridX+1&&j==0)
			{
				SendN[2] +=1;
				PackP(p, Sepm);
				p = Reconnect(p);
				continue;
			}
			if(i==GridX+1&&j==GridY+1)
			{
				SendN[3] +=1;
				PackP(p, Sepp);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to left Neighbor =======
			//====================================
			if(i==0)
			{
				SendN[4] +=1;
				PackP(p, SeXm);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to right Neighbor ======
			//====================================
			if(i==GridX+1)
			{
				SendN[5] +=1;
				PackP(p, SeXp);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to Neighbor Below ======
			//====================================
			if(j==0)
			{
				SendN[6] +=1;
				PackP(p, SeYm);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to Neighbor Above ======
			//====================================
			if(j==GridY+1)
			{
				SendN[7] +=1;
				PackP(p, SeYp);
				p = Reconnect(p);
				continue;
			}
			p = p->p_PrevPart;
		

		}


		//=====================================================================
		//======== Exchange Partectoryies with Neighboring Processors =========
		//=====================================================================
		if(bufsize/GridX<*std::max_element(SendN.begin(), SendN.begin()+4)||bufsize<*std::max_element(SendN.begin()+5, SendN.end()))
		{	
			printf("==== Mesh: At Rank: %5d. ==================\n",Rank);
			std::cout << "==== Mesh: Send Too Many Particles.      ====\n";
			std::cout << "==== Mesh: May Cause Memory Problems.    ====\n";
			// printf("==== Mesh: Try Buf Size: %3d           ====\n",maxsend/GridX);
			
		}

		S_SUM = SendN[0]+SendN[1]+SendN[2]+SendN[3]+SendN[4]+SendN[5]+SendN[6]+SendN[7];
		MPI_Allreduce(&S_SUM, &A_SUM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( A_SUM == 0 ) {break;};
		p_domain()->p_Com()->DoCommuteT(COMMU_P, SendN);
	}

	return;
}



void Mesh::PackP(Particle* p_Part, WDOUBLE* &Se)
{
	
	*Se = p_Part-> x;  Se++;
	*Se = p_Part-> y;  Se++;
	*Se = p_Part-> z;  Se++;

	*Se = p_Part-> x0; Se++;
	*Se = p_Part-> y0; Se++;
	*Se = p_Part-> z0; Se++;

	*Se = p_Part-> px; Se++;
	*Se = p_Part-> py; Se++;
	*Se = p_Part-> pz; Se++;

	*Se = p_Part-> Ex0; Se++;
	*Se = p_Part-> Ey0; Se++;
	*Se = p_Part-> Ez0; Se++;

	*Se = p_Part-> type*1.0; Se++;
	*Se = p_Part-> q2m; 	 Se++;
	*Se = p_Part-> weight;   Se++;

	*Se = p_Part-> Wxw; Se++;
	*Se = p_Part-> Wyw; Se++;
	*Se = p_Part-> Wzw; Se++;

	*Se = p_Part-> Wxl; Se++;
	*Se = p_Part-> Wyl; Se++;
	*Se = p_Part-> Wzl; Se++;

	*Se = p_Part-> sx; Se++;
	*Se = p_Part-> sy; Se++;

	return;
}