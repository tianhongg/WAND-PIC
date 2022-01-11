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


void Mesh::PushTrajectory(WDOUBLE k0, int k, int step)
{

	// temporary remove
	return;

}

void Mesh::PushTrajectory_Half()
{

	Trajectory *p = NULL;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();
	
	WDOUBLE xt, yt, Vx, Vy, xtp, ytp;
	WDOUBLE dztmp;

	p = p_Trajectory;

	int i;
	int j;

	while (p)
	{

		xt = p-> x;
		yt = p-> y;

		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}

		if(RankIdx_X == Xpa && i==GridX+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}

		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}

		if(RankIdx_Y == Ypa && j==GridY+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		// Push  position only......
		dztmp=dzz;
		xtp = p-> old_x  + Vx*dztmp;
		ytp = p-> old_y  + Vy*dztmp;

		p-> x = xtp;
		p-> y = ytp;

		p-> old_x = xtp;
		p-> old_y = ytp;

		// update the cell index;
		Cell *c = &GetCell(i,j,0);
		while(xtp>(c->Xcord+c->dx*0.5) && i<GridX+1) { p->idx_i++;  i++; c=&GetCell(i,j,0); }
		while(xtp<(c->Xcord-c->dx*0.5) && i>0) 		 { p->idx_i--; 	i--; c=&GetCell(i,j,0); }
		while(ytp>(c->Ycord+c->dy*0.5) && j<GridY+1) { p->idx_j++;  j++; c=&GetCell(i,j,0); }
		while(ytp<(c->Ycord-c->dy*0.5) && j>0) 		 { p->idx_j--;  j--; c=&GetCell(i,j,0); }
		//
		p = p->p_PrevTraj;

	}
	//============================================
	//=========== Exchange Particles =============
	ExchangeT();
	//============================================
	return;

}

void Mesh::PushTrajectory_HalfE(int k) 
{

	Trajectory *p = NULL;

	WDOUBLE ddx, ddy;

	WDOUBLE Ex, Ey, Ez, Psi, Pondx, Pondy, Asq;
	WDOUBLE Fx, Fy, gamma;

	WDOUBLE xt, yt, Vx, Vy, Vxp, Vyp, xtp, ytp;
	
	WDOUBLE wmm,wmc,wmp;
	WDOUBLE wcm,wcc,wcp;
	WDOUBLE wpm,wpc,wpp;

	WDOUBLE dztmp;
	WDOUBLE sx, sy, sxy;

	int i,j;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	p = p_Trajectory;
	Vmax = 0.0;

	while (p)
	{

		xt = p-> x;
		yt = p-> y;

		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		Cell &cmm = GetCell(i-1,j-1,k);
		Cell &cmc = GetCell(i-1,j  ,k);
		Cell &cmp = GetCell(i-1,j+1,k);

		Cell &ccm = GetCell(i,  j-1,k);
		Cell &ccc = GetCell(i,  j  ,k);
		Cell &ccp = GetCell(i,  j+1,k);
		
		Cell &cpm = GetCell(i+1,j-1,k);
		Cell &cpc = GetCell(i+1,j  ,k);
		Cell &cpp = GetCell(i+1,j+1,k);

		ddx=ccc.dx;
		ddy=ccc.dy;

		// sx=ddx/TpCellx;  //- re-size
		// sy=ddy/TpCelly;  //- re-size


		// this one is better
		sx=ddx;  //- re-size
		sy=ddy;  //- re-size

		sxy = sx*sy;

		WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+xt-ccc.Xcord),0.0);
		WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+yt-ccc.Ycord),0.0);

		WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-xt+ccc.Xcord),0.0);
		WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-yt+ccc.Ycord),0.0);

		WDOUBLE deltaxc=sx-deltaxm-deltaxp;
		WDOUBLE deltayc=sy-deltaym-deltayp;

		wmm = deltaxm*deltaym/sxy;
		wmc = deltaxm*deltayc/sxy;
		wmp = deltaxm*deltayp/sxy;

		wcm = deltaxc*deltaym/sxy;
		wcc = deltaxc*deltayc/sxy;
		wcp = deltaxc*deltayp/sxy;

		wpm = deltaxp*deltaym/sxy;
		wpc = deltaxp*deltayc/sxy;
		wpp = deltaxp*deltayp/sxy;

		Ex  = wmm*cmm.W_Ex  + wmc*cmc.W_Ex  + wmp*cmp.W_Ex  
			+ wcm*ccm.W_Ex  + wcc*ccc.W_Ex  + wcp*ccp.W_Ex
			+ wpm*cpm.W_Ex  + wpc*cpc.W_Ex  + wpp*cpp.W_Ex;

		Ey  = wmm*cmm.W_Ey  + wmc*cmc.W_Ey  + wmp*cmp.W_Ey  
			+ wcm*ccm.W_Ey  + wcc*ccc.W_Ey  + wcp*ccp.W_Ey
			+ wpm*cpm.W_Ey  + wpc*cpc.W_Ey  + wpp*cpp.W_Ey;

		Ez  = wmm*cmm.W_Ez  + wmc*cmc.W_Ez  + wmp*cmp.W_Ez  
			+ wcm*ccm.W_Ez  + wcc*ccc.W_Ez  + wcp*ccp.W_Ez
			+ wpm*cpm.W_Ez  + wpc*cpc.W_Ez  + wpp*cpp.W_Ez;

		Psi = wmm*cmm.W_Psi  + wmc*cmc.W_Psi  + wmp*cmp.W_Psi  
			+ wcm*ccm.W_Psi  + wcc*ccc.W_Psi  + wcp*ccp.W_Psi
			+ wpm*cpm.W_Psi  + wpc*cpc.W_Psi  + wpp*cpp.W_Psi;

		Asq = wmm*cmm.W_Asq  + wmc*cmc.W_Asq  + wmp*cmp.W_Asq  
			+ wcm*ccm.W_Asq  + wcc*ccc.W_Asq  + wcp*ccp.W_Asq
			+ wpm*cpm.W_Asq  + wpc*cpc.W_Asq  + wpp*cpp.W_Asq;

		Pondx= wmm*cmm.W_Ponx  + wmc*cmc.W_Ponx  + wmp*cmp.W_Ponx  
			 + wcm*ccm.W_Ponx  + wcc*ccc.W_Ponx  + wcp*ccp.W_Ponx
			 + wpm*cpm.W_Ponx  + wpc*cpc.W_Ponx  + wpp*cpp.W_Ponx;

		Pondy= wmm*cmm.W_Pony  + wmc*cmc.W_Pony  + wmp*cmp.W_Pony  
			 + wcm*ccm.W_Pony  + wcc*ccc.W_Pony  + wcp*ccp.W_Pony
			 + wpm*cpm.W_Pony  + wpc*cpc.W_Pony  + wpp*cpp.W_Pony;

		gamma = 0.5*(1+Psi)*(Vx*Vx+Vy*Vy+1)+0.5*(1.0+0.5*Asq)/(1+Psi);
		Fx = ((gamma*Ex-Pondx*0.25)/(1+Psi) - Vx*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);
		Fy = ((gamma*Ey-Pondy*0.25)/(1+Psi) - Vy*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);

		dztmp=dzz;
		Vxp = p-> old_vx + Fx*dztmp*0.5;
		Vyp = p-> old_vy + Fy*dztmp*0.5;

		//==========================================
		//=========== Adapteive Z Step =============
		WDOUBLE Vr = sqrt(Vxp*Vxp+Vyp*Vyp);
		if(AdaptiveStep>0 && Vr >=Vlim*AdaptiveStep)
		{
			Vxp = Vlim*AdaptiveStep*Vxp/Vr;
			Vyp = Vlim*AdaptiveStep*Vyp/Vr;
		}

		p-> Vx = Vxp;
		p-> Vy = Vyp;

		p-> old_vx = Vxp;
		p-> old_vy = Vyp;
		
		p-> Vxx = Vxp*Vxp;
		p-> Vxy = Vxp*Vyp;
		p-> Vyy = Vyp*Vyp;

		p = p->p_PrevTraj;

		WDOUBLE Vrr=std::max(abs(Vxp)*dz/ccc.dx,abs(Vyp)*dz/ccc.dy);
		Vmax = std::max(Vmax, Vrr); // how many grids it can cross

	}
	return;

}


void Mesh::PushTrajectory_HalfB(int k)
{
	Trajectory *p = NULL;

	WDOUBLE ddx,ddy;
	WDOUBLE Psi, Bx, By, Bz;
	WDOUBLE Fx, Fy, gamma;

	WDOUBLE xt, yt, Vx, Vy, Vxp, Vyp, xtp, ytp;
	
	WDOUBLE wmm,wmc,wmp;
	WDOUBLE wcm,wcc,wcp;
	WDOUBLE wpm,wpc,wpp;

	WDOUBLE dztmp;
	int i,j;
	WDOUBLE sx,sy, sxy;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	p = p_Trajectory;
	//====== for adpative z step========
	Vmax = 0.0;
	while (p)
	{

		xt = p-> x;
		yt = p-> y;

		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;
		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY+1)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		Cell &cmm = GetCell(i-1,j-1,k);
		Cell &cmc = GetCell(i-1,j  ,k);
		Cell &cmp = GetCell(i-1,j+1,k);

		Cell &ccm = GetCell(i,  j-1,k);
		Cell &ccc = GetCell(i,  j  ,k);
		Cell &ccp = GetCell(i,  j+1,k);
		
		Cell &cpm = GetCell(i+1,j-1,k);
		Cell &cpc = GetCell(i+1,j  ,k);
		Cell &cpp = GetCell(i+1,j+1,k);

		ddx=ccc.dx;
		ddy=ccc.dy;

		// sx=ddx/TpCellx;  //- re-size
		// sy=ddy/TpCelly;  //- re-size

		// this one is better
		sx=ddx;  //- re-size
		sy=ddy;  //- re-size

		sxy = sx*sy;

		WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+xt-ccc.Xcord),0.0);
		WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+yt-ccc.Ycord),0.0);

		WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-xt+ccc.Xcord),0.0);
		WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-yt+ccc.Ycord),0.0);

		WDOUBLE deltaxc=sx-deltaxm-deltaxp;
		WDOUBLE deltayc=sy-deltaym-deltayp;

		wmm = deltaxm*deltaym/sxy;
		wmc = deltaxm*deltayc/sxy;
		wmp = deltaxm*deltayp/sxy;

		wcm = deltaxc*deltaym/sxy;
		wcc = deltaxc*deltayc/sxy;
		wcp = deltaxc*deltayp/sxy;

		wpm = deltaxp*deltaym/sxy;
		wpc = deltaxp*deltayc/sxy;
		wpp = deltaxp*deltayp/sxy;

		Bx  = wmm*cmm.W_Bx  + wmc*cmc.W_Bx  + wmp*cmp.W_Bx  
			+ wcm*ccm.W_Bx  + wcc*ccc.W_Bx  + wcp*ccp.W_Bx
			+ wpm*cpm.W_Bx  + wpc*cpc.W_Bx  + wpp*cpp.W_Bx;

		By  = wmm*cmm.W_By  + wmc*cmc.W_By  + wmp*cmp.W_By  
			+ wcm*ccm.W_By  + wcc*ccc.W_By  + wcp*ccp.W_By
			+ wpm*cpm.W_By  + wpc*cpc.W_By  + wpp*cpp.W_By;

		Bz  = wmm*cmm.W_Bz  + wmc*cmc.W_Bz  + wmp*cmp.W_Bz  
			+ wcm*ccm.W_Bz  + wcc*ccc.W_Bz  + wcp*ccp.W_Bz
			+ wpm*cpm.W_Bz  + wpc*cpc.W_Bz  + wpp*cpp.W_Bz;

		Psi = wmm*cmm.W_Psi  + wmc*cmc.W_Psi  + wmp*cmp.W_Psi  
			+ wcm*ccm.W_Psi  + wcc*ccc.W_Psi  + wcp*ccp.W_Psi
			+ wpm*cpm.W_Psi  + wpc*cpc.W_Psi  + wpp*cpp.W_Psi;
		
		Fx = (- Vy*Bz - By )/(1+Psi);
		Fy = (  Vx*Bz + Bx )/(1+Psi);

		dztmp=dzz;
		Vxp = p-> old_vx + Fx*dztmp;
		Vyp = p-> old_vy + Fy*dztmp;

		p-> old_vx = Vxp;
		p-> old_vy = Vyp;
		
		p-> Vxx = Vxp*Vxp;
		p-> Vxy = Vxp*Vyp;
		p-> Vyy = Vyp*Vyp;

		p = p->p_PrevTraj;


	}

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


void Mesh::PackT(Trajectory* p_Traj, WDOUBLE* &Se)
{

	*Se = p_Traj-> x;  Se++;
	*Se = p_Traj-> y;  Se++;
	*Se = p_Traj-> x0; Se++;
	*Se = p_Traj-> y0; Se++;
	*Se = p_Traj-> z0; Se++;
	*Se = p_Traj-> Vx; Se++;
	*Se = p_Traj-> Vy; Se++;
	*Se = p_Traj-> old_x; Se++;
	*Se = p_Traj-> old_y; Se++;
	*Se = p_Traj-> old_vx; Se++;
	*Se = p_Traj-> old_vy; Se++;
	*Se = p_Traj-> sx; Se++;
	*Se = p_Traj-> sy; Se++;

	return;
}



void Mesh::ExchangeT()
{
	Trajectory *p = NULL;

	Commute *p_COMM = p_domain()->p_Com();
	
	//=========Send and Receive Buf Size===============
	WDOUBLE bufsize = p_domain()->p_Com()->Get_bufsize();
	bufsize *= (GridX*SOU_DIM*2.0/SDT_DIM);
	//=================================================

	WDOUBLE xtp, ytp;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	std::vector<int> SendN(8,0);//Sendmm, Sendmp, Sendpm, Sendpp; Sendxm, Sendxp, Sendym, Sendyp;
	

	int S_SUM, A_SUM;
	int i;
	int j;

		//      ___________
		//     |mp | yp| pp|
		//     |___|___|___|
		//     |xm |   | xp|
		//     |___|___|___|
		//     |mm | ym| pm|
		//     |___|___|___|
		//
	while(1)
	{

		p = p_Trajectory;

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

			//=================================================
			//============Trajectory Outside Boundary =========
			//=================================================
			if(p_domain()->Get_BC()==1)
			{
				if(RankIdx_X ==1	&& i==0)
				{ p = p->p_PrevTraj; continue;}

				if(RankIdx_X == Xpa && i==GridX+1)
				{ p = p->p_PrevTraj; continue;}

				if(RankIdx_Y == 1	&& j==0)
				{ p = p->p_PrevTraj; continue;}

				if(RankIdx_Y == Ypa && j==GridY+1)
				{ p = p->p_PrevTraj; continue;}
			}
			//==================================================
			if(i==0&&j==0)
			{
				SendN[0] +=1;
				PackT(p, Semm);
				p = Reconnect(p);
				continue;
			}

			if(i==0&&j==GridY+1)
			{
				SendN[1] +=1;
				PackT(p, Semp);
				p = Reconnect(p);
				continue;
			}
			if(i==GridX+1&&j==0)
			{
				SendN[2] +=1;
				PackT(p, Sepm);
				p = Reconnect(p);
				continue;
			}
			if(i==GridX+1&&j==GridY+1)
			{
				SendN[3] +=1;
				PackT(p, Sepp);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to left Neighbor =======
			//====================================
			if(i==0)
			{
				SendN[4] +=1;
				PackT(p, SeXm);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to right Neighbor ======
			//====================================
			if(i==GridX+1)
			{
				SendN[5] +=1;
				PackT(p, SeXp);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to Neighbor Below ======
			//====================================
			if(j==0)
			{
				SendN[6] +=1;
				PackT(p, SeYm);
				p = Reconnect(p);
				continue;
			}
			//====================================
			//====== Send to Neighbor Above ======
			//====================================
			if(j==GridY+1)
			{
				SendN[7] +=1;
				PackT(p, SeYp);
				p = Reconnect(p);
				continue;
			}
			p = p->p_PrevTraj;
		
		}

		//=====================================================================
		//======== Exchange Trajectoryies with Neighboring Processors =========
		//=====================================================================
		if(bufsize/GridX<*std::max_element(SendN.begin(), SendN.begin()+4)||bufsize<*std::max_element(SendN.begin()+5, SendN.end()))
		{	
			printf("==== Mesh: At Rank: %5d. ==================\n",Rank);
			std::cout << "==== Mesh: Send Too Many Trajectoryies.  ====\n";
			std::cout << "==== Mesh: May Cause Memory Problems.    ====\n";
			std::cout << "==== Mesh: Try to Increase the Buf Size. ====\n";
		}

		S_SUM = SendN[0]+SendN[1]+SendN[2]+SendN[3]+SendN[4]+SendN[5]+SendN[6]+SendN[7];
		MPI_Allreduce(&S_SUM, &A_SUM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( A_SUM == 0 ) {break;};
		p_domain()->p_Com()->DoCommuteT(COMMU_T, SendN);
	}
	return;
}

void Mesh::AdjustZstep(WDOUBLE k0, int k, WDOUBLE &dz2dz)
{
	//==========================================
	//=========== Adapteive Z Step =============
	if(AdaptiveStep>0)
	{	
		WDOUBLE A_Vmax;
		MPI_Allreduce(&Vmax, &A_Vmax, 1, MPI_WDOUBLE, MPI_MAX, MPI_COMM_WORLD);
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
