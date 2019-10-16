//----------------------------------------------------------------------------------||
//-------------------                pushtrajs.cpp               -------------------||
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
//---Starting---------           : Feb-22-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#include "wand_PIC.h"


void Mesh::PushTrajectory(double k0, int k, int step)
{

	Trajectory *p = NULL;
	double ddx;
	double ddy;
	double Ex, Ey, Ez, Psi, Bx, By, Bz, Pondx, Pondy, Asq;
	double Fx, Fy, gamma;

	double xt, yt, Vx, Vy, Vxp, Vyp, xtp, ytp;
	double wmm, wmp, wpm, wpp;
	double dxdy = dx*dy;
	double dztmp;
	int i,j;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;

	p = p_Trajectory;
	//====== for adpative z step========
	Vmax = 0.0;
	while (p)
	{

		xt = p-> x;
		yt = p-> y;
		Vx = p-> Vx;
		Vy = p-> Vy;

		ddx = xt-(Offset_X-dx*0.5);
		ddy = yt-(Offset_Y-dy*0.5);

		i = floor(ddx/dx);
		j = floor(ddy/dy);

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	& xt<=Offset_X)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa & xt>=Xmax)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	& yt<=Offset_Y)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa & yt>=Ymax)
		{ p = p->p_PrevTraj; continue;}
		//==================================================


		if(i < 0 || i > GridX || j < 0 || j > GridY)
		{
			std::cout <<"Rank:"<<Rank<<"==== Mesh: Wrong Trajectory Position. ====\n";
			std::cout <<"Rank:"<<Rank<<"==== Mesh: Check Trajs' Send & Rece.  ====\n";
			std::cout <<"Rank:"<<Rank<<"==== Mesh: Check Fields' Resolutions. ====\n";
			std::cout <<"Rank:"<<Rank<<"==== Mesh: Check Domain Resolutions.  ====\n";
			exit(21);
		}

		wmm = (i+1-ddx/dx)*(j+1-ddy/dy);
		wmp = (i+1-ddx/dx)*(ddy/dy-j);
		wpm = (ddx/dx-i)*(j+1-ddy/dy);
		wpp = (ddx/dx-i)*(ddy/dy-j);

		Cell &cmm = GetCell(i,j,  	k);
		Cell &cmp = GetCell(i,j+1,  k);
		Cell &cpm = GetCell(i+1,j,  k);
		Cell &cpp = GetCell(i+1,j+1,k);

		Ex  = wmm*cmm.W_Ex  + wmp*cmp.W_Ex  + wpm*cpm.W_Ex  + wpp*cpp.W_Ex;
		Ey  = wmm*cmm.W_Ey  + wmp*cmp.W_Ey  + wpm*cpm.W_Ey  + wpp*cpp.W_Ey;
		Ez  = wmm*cmm.W_Ez  + wmp*cmp.W_Ez  + wpm*cpm.W_Ez  + wpp*cpp.W_Ez;

		Bx  = wmm*cmm.W_Bx  + wmp*cmp.W_Bx  + wpm*cpm.W_Bx  + wpp*cpp.W_Bx;
		By  = wmm*cmm.W_By  + wmp*cmp.W_By  + wpm*cpm.W_By  + wpp*cpp.W_By;
		Bz  = wmm*cmm.W_Bz  + wmp*cmp.W_Bz  + wpm*cpm.W_Bz  + wpp*cpp.W_Bz;

		Psi   = wmm*cmm.W_Psi  + wmp*cmp.W_Psi  + wpm*cpm.W_Psi  + wpp*cpp.W_Psi;
		Pondx = wmm*cmm.W_Ponx + wmp*cmp.W_Ponx + wpm*cpm.W_Ponx + wpp*cpp.W_Ponx;
		Pondy = wmm*cmm.W_Pony + wmp*cmp.W_Pony + wpm*cpm.W_Pony + wpp*cpp.W_Pony;
		
		Asq   = wmm*cmm.W_Asq  + wmp*cmp.W_Asq  + wpm*cpm.W_Asq  + wpp*cpp.W_Asq;

		gamma = 0.5*(1+Psi)*(Vx*Vx+Vy*Vy+1)+0.5*(1.0+0.5*Asq)/(1+Psi);
		
		Fx = ((gamma*Ex-Pondx*0.25)/(1+Psi) - Vy*Bz - By - Vx*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);
		Fy = ((gamma*Ey-Pondy*0.25)/(1+Psi) + Vx*Bz + Bx - Vy*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);

		switch(step)
		{
		case 0:

			dztmp=dzz;
			Vxp = p-> old_vx + Fx*dztmp;
			Vyp = p-> old_vy + Fy*dztmp;

			xtp = p-> old_x  + Vx*dztmp;
			ytp = p-> old_y  + Vy*dztmp;

		break;
		case 1:

			dztmp=dzz*0.5;
			Vxp = Vx + Fx*dztmp;
			Vyp = Vy + Fy*dztmp;

			xtp = xt + Vx*dztmp;
			ytp = yt + Vy*dztmp;

		break;
		}

		//==========================================
		//=========== Adapteive Z Step =============
		double Vr = sqrt(Vxp*Vxp+Vyp*Vyp);
		if(AdaptiveStep>0 & Vr >=Vlim*AdaptiveStep)
		{
			Vxp = Vlim*AdaptiveStep*Vxp/Vr;
			Vyp = Vlim*AdaptiveStep*Vyp/Vr;
		}


		p-> x = xtp;
		p-> y = ytp;
		p-> Vx = Vxp;
		p-> Vy = Vyp;

		if(step==0)
		{
		p-> old_x = xtp;
		p-> old_y = ytp;
		p-> old_vx = Vxp;
		p-> old_vy = Vyp;
		}
		
		p-> Vxx = Vxp*Vxp;
		p-> Vxy = Vxp*Vyp;
		p-> Vyy = Vyp*Vyp;

		p = p->p_PrevTraj;

		double Vrr=sqrt(Vxp*Vxp+Vyp*Vyp);
		if(Vrr>Vmax) Vmax = Vrr;

	}
	//============================================
	//=========== Exchange Particles =============
	ExchangeT();
	//============================================
	return;

}

Trajectory* Mesh::Reconnect(Trajectory* p_Traj)
{

	Trajectory* p_temp;

	if(p_Traj->p_PrevTraj)
	{
		p_Traj->p_PrevTraj->p_NextTraj = p_Traj->p_NextTraj;
		p_temp = p_Traj->p_PrevTraj;
	}
	else
	{
		p_temp = NULL;
	}

	if(p_Traj->p_NextTraj)
	{
		p_Traj->p_NextTraj->p_PrevTraj = p_Traj->p_PrevTraj;
	}
	else
	{
		p_Trajectory = p_Traj->p_PrevTraj;

	}

	delete p_Traj;
	return p_temp;

}




void Mesh::PackT(Trajectory* p_Traj, int Sendn, int where)
{
	Commute *p_COMM = p_domain()->p_Com();

	switch(where)
	{

		case 0:
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;

		break;
	
		case 1:
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;
		break;

		case 2:
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;
		break;

		case 3:
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;	
		break;
	}

	return;
}



void Mesh::ExchangeT()
{
	Trajectory *p = NULL;
	
	//=========Send and Receive Buf Size===============
	double bufsize = p_domain()->p_Com()->Get_bufsize();
	bufsize *= (GridX*SOU_DIM*2.0/SDT_DIM);
	//=================================================

	double xtp, ytp;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;
	double Xsize = GridX*dx*Xpa;
	double Ysize = GridY*dy*Ypa;

	int Sendxm, Sendxp, Sendym, Sendyp;
	int S_SUM, A_SUM;

	while(1)
	{

		Sendxm = Sendxp = Sendym = Sendyp = 0; 
		p = p_Trajectory;
		while (p)
		{
			xtp = p-> x;
			ytp = p-> y;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(p_domain()->Get_BC()==1)
		{
		if(RankIdx_X ==1	&& xtp<=Offset_X)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && xtp>=Xmax)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& ytp<=Offset_Y)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && ytp>=Ymax)
		{ p = p->p_PrevTraj; continue;}
		}
		//==================================================

		//====================================
		//====== Send to left Neighbor =======
		//====================================
		if(xtp < Offset_X)
		{
			if(RankIdx_X ==1) 	p-> x = xtp + Xsize;
			Sendxm +=1;
			PackT(p, Sendxm, 0);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to right Neighbor ======
		//====================================
		else if(xtp > Xmax)
		{
			if(RankIdx_X ==Xpa) p-> x = xtp - Xsize;
			Sendxp +=1;
			PackT(p, Sendxp, 1);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to Neighbor Below ======
		//====================================
		else if(ytp < Offset_Y)
		{
			if(RankIdx_Y ==1) 	p-> y = ytp + Ysize;
			Sendym +=1;
			PackT(p, Sendym, 2);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to Neighbor Above ======
		//====================================
		else if(ytp > Ymax)
		{
			if(RankIdx_Y ==Ypa) p-> y = ytp - Ysize;
			Sendyp +=1;
			PackT(p, Sendyp, 3);
			p = Reconnect(p);
		}
		else
		{
			p = p->p_PrevTraj;
		}

		}

		//=====================================================================
		//======== Exchange Trajectoryies with Neighboring Processors =========
		//=====================================================================
		if(bufsize<Sendxm || bufsize<Sendxp || bufsize<Sendym || bufsize<Sendyp)
		{	
			printf("==== Mesh: At Rank: %5d. ==================\n",Rank);
			std::cout << "==== Mesh: Send Too Many Trajectoryies.  ====\n";
			std::cout << "==== Mesh: May Cause Memory Problems.    ====\n";
			std::cout << "==== Mesh: Try to Increase the Buf Size. ====\n";
			
		}

		S_SUM = Sendxm+Sendxp+Sendym+Sendyp;

		MPI_Allreduce(&S_SUM, &A_SUM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( A_SUM == 0 ) {break;};

		p_domain()->p_Com()->DoCommuteT(COMMU_T, Sendxm, Sendxp, Sendym, Sendyp);

	}

	return;
}

void Mesh::AdjustZstep(double k0, int k, double &dz2dz)
{
	//==========================================
	//=========== Adapteive Z Step =============
	if(AdaptiveStep>0)
	{	
		double A_Vmax;
		MPI_Allreduce(&Vmax, &A_Vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if(A_Vmax>Vlim) 
 		{
			dzz = dz/(A_Vmax/Vlim);
		}
		else
		{
			dzz = dz;
		}
		if(k0 + dzz/dz > k &dzz<dz &(k!=k0)) dzz=(k-k0)*dz;
	}

	dz2dz= dzz/dz;
	//==========================================
	return;
}
