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



void Mesh::MacroSource(int k) //v
{

	Trajectory *p = NULL;

	WDOUBLE wmm,wmc,wmp;
	WDOUBLE wcm,wcc,wcp;
	WDOUBLE wpm,wpc,wpp;

	int i,j;
	WDOUBLE ddx; //cell size
	WDOUBLE ddy;

	WDOUBLE sx; // particle size;
	WDOUBLE sy;


	SetSourceZero(k);
	
	p = p_Trajectory;
	while (p)
	{
		WDOUBLE xt 		=  p->x;
		WDOUBLE yt 		=  p->y;
		WDOUBLE massweig = (p->Weight);

		//------------------------------
		//      ___________
		//     |mp | cp| pp|
		//     |___|___|___|
		//     |mc | cc| pc|
		//     |___|___|___|
		//     |mm | cm| pm|
		//     |___|___|___|
		//
		//------------------------------

		// idex of the cell
		i = p->idx_i;
		j = p->idx_j;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		// 0 - Grid+1
		if(i < 1 || i > GridX || j < 1 || j > GridY)
		{
			p = p->p_PrevTraj;
			continue;
		}
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

		sx=ddx;  //- re-size to make the results more smooth
		sy=ddy;  //- re-size

		massweig *= (p->sx)*(p->sy)/sx/sy; // re-weight

		WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+xt-ccc.Xcord),0.0);
		WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+yt-ccc.Ycord),0.0);

		WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-xt+ccc.Xcord),0.0);
		WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-yt+ccc.Ycord),0.0);

		WDOUBLE deltaxc=sx-deltaxm-deltaxp;
		WDOUBLE deltayc=sy-deltaym-deltayp;
		
		wmm = deltaxm*deltaym/cmm.dx/cmm.dy;
		wmc = deltaxm*deltayc/cmc.dx/cmc.dy;
		wmp = deltaxm*deltayp/cmp.dx/cmp.dy;

		wcm = deltaxc*deltaym/ccm.dx/ccm.dy;
		wcc = deltaxc*deltayc/ccc.dx/ccc.dy;
		wcp = deltaxc*deltayp/ccp.dx/ccp.dy;

		wpm = deltaxp*deltaym/cpm.dx/cpm.dy;
		wpc = deltaxp*deltayc/cpc.dx/cpc.dy;
		wpp = deltaxp*deltayp/cpp.dx/cpp.dy;

		for (int n=0; n<SOU_DIM; n++)
		{
			cmm.W_Source[n] += massweig * wmm * p->T_Source[n];
			cmc.W_Source[n] += massweig * wmc * p->T_Source[n];
			cmp.W_Source[n] += massweig * wmp * p->T_Source[n];

			ccm.W_Source[n] += massweig * wcm * p->T_Source[n];
			ccc.W_Source[n] += massweig * wcc * p->T_Source[n];
			ccp.W_Source[n] += massweig * wcp * p->T_Source[n];

			cpm.W_Source[n] += massweig * wpm * p->T_Source[n];
			cpc.W_Source[n] += massweig * wpc * p->T_Source[n];
			cpp.W_Source[n] += massweig * wpp * p->T_Source[n];

		}
		p = p->p_PrevTraj;

	}
	
	return;

}

void Mesh::SetSourceZero(int k) //v
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


//deprecated//May-24-tianhong
void Mesh::AdjustSource(int k) //x
{
	
	return;

}



//===============================
//======== djx/dx+djy/dy ========
//===============================
WDOUBLE Mesh::Dive_J(int i, int j, WDOUBLE k0, int k) //v
{

	int km=k-1;
	if(k==0) km=k;
	WDOUBLE kzz=k-k0;

	WDOUBLE hp,hm,ha;

	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	Cell &cccm = GetCell(i, j,km);
	Cell &cxmm = GetCell(i-1,j,km);
	Cell &cxpm = GetCell(i+1,j,km);
	Cell &cymm = GetCell(i,j-1,km);
	Cell &cypm = GetCell(i,j+1,km);

	hp=(cxp.dx+ccc.dx)*0.5;
	hm=(ccc.dx+cxm.dx)*0.5;
	ha=hp+hm;

	WDOUBLE cccB_Jx=ccc.B_Jx*(1-kzz)+cccm.B_Jx*kzz;
	WDOUBLE cxmB_Jx=cxm.B_Jx*(1-kzz)+cxmm.B_Jx*kzz;
	WDOUBLE cxpB_Jx=cxp.B_Jx*(1-kzz)+cxpm.B_Jx*kzz;
	WDOUBLE cymB_Jx=cym.B_Jx*(1-kzz)+cymm.B_Jx*kzz;
	WDOUBLE cypB_Jx=cyp.B_Jx*(1-kzz)+cypm.B_Jx*kzz;

	WDOUBLE cccB_Jy=ccc.B_Jy*(1-kzz)+cccm.B_Jy*kzz;
	WDOUBLE cxmB_Jy=cxm.B_Jy*(1-kzz)+cxmm.B_Jy*kzz;
	WDOUBLE cxpB_Jy=cxp.B_Jy*(1-kzz)+cxpm.B_Jy*kzz;
	WDOUBLE cymB_Jy=cym.B_Jy*(1-kzz)+cymm.B_Jy*kzz;
	WDOUBLE cypB_Jy=cyp.B_Jy*(1-kzz)+cypm.B_Jy*kzz;

	WDOUBLE DJ = (cxp.W_Jx - ccc.W_Jx + cxpB_Jx - cccB_Jx)*hm/hp/ha + (ccc.W_Jx - cxm.W_Jx + cccB_Jx - cxmB_Jx)*hp/hm/ha;
	
	hp=(cyp.dy+ccc.dy)*0.5;
	hm=(ccc.dy+cym.dy)*0.5;
	ha=hp+hm;

	return DJ + (cyp.W_Jy - ccc.W_Jy + cypB_Jy - cccB_Jy)*hm/hp/ha + (ccc.W_Jy - cym.W_Jy + cccB_Jy - cymB_Jy)*hp/hm/ha;
	
}



//===============================
//======== djy/dx-djx/dy ========
//===============================
WDOUBLE Mesh::Curl_J(int i, int j, WDOUBLE k0, int k) //v
{

	int km=k-1;
	if(k==0) km=k;

	WDOUBLE kzz=k-k0;
	WDOUBLE hp,hm,ha;

	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	Cell &cccm = GetCell(i,  j,km);
	Cell &cxmm = GetCell(i-1,j,km);
	Cell &cxpm = GetCell(i+1,j,km);
	Cell &cymm = GetCell(i,j-1,km);
	Cell &cypm = GetCell(i,j+1,km);

	hp=(cxp.dx+ccc.dx)*0.5;
	hm=(ccc.dx+cxm.dx)*0.5;
	ha=hp+hm;

	WDOUBLE cccB_Jx=ccc.B_Jx*(1-kzz)+cccm.B_Jx*kzz;
	WDOUBLE cxmB_Jx=cxm.B_Jx*(1-kzz)+cxmm.B_Jx*kzz;
	WDOUBLE cxpB_Jx=cxp.B_Jx*(1-kzz)+cxpm.B_Jx*kzz;
	WDOUBLE cymB_Jx=cym.B_Jx*(1-kzz)+cymm.B_Jx*kzz;
	WDOUBLE cypB_Jx=cyp.B_Jx*(1-kzz)+cypm.B_Jx*kzz;

	WDOUBLE cccB_Jy=ccc.B_Jy*(1-kzz)+cccm.B_Jy*kzz;
	WDOUBLE cxmB_Jy=cxm.B_Jy*(1-kzz)+cxmm.B_Jy*kzz;
	WDOUBLE cxpB_Jy=cxp.B_Jy*(1-kzz)+cxpm.B_Jy*kzz;
	WDOUBLE cymB_Jy=cym.B_Jy*(1-kzz)+cymm.B_Jy*kzz;
	WDOUBLE cypB_Jy=cyp.B_Jy*(1-kzz)+cypm.B_Jy*kzz;

	WDOUBLE CJ = (cxp.W_Jy - ccc.W_Jy + cxpB_Jy - cccB_Jy)*hm/hp/ha + (ccc.W_Jy - cxm.W_Jy + cccB_Jy - cxmB_Jy)*hp/hm/ha;

	hp=(cyp.dy+ccc.dy)*0.5;
	hm=(ccc.dy+cym.dy)*0.5;
	ha=hp+hm;

	return CJ - (cyp.W_Jx - ccc.W_Jx + cypB_Jx - cccB_Jx)*hm/hp/ha - (ccc.W_Jx - cym.W_Jx + cccB_Jx - cymB_Jx)*hp/hm/ha;
}



void Mesh::Put_Jz(int k) //v
{

	WDOUBLE Asq;
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
WDOUBLE Mesh::SourceX(int i, int j, WDOUBLE k0, int k) //v
{
	int km=k-1;
	if(k==0) km=k;
	WDOUBLE kzz=k-k0;
	WDOUBLE hxp,hxm,hxa;
	WDOUBLE hyp,hym,hya;

	WDOUBLE gamma, ax, sx, Asq;

	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	Cell &cccm = GetCell(i,  j,km);
	Cell &cxmm = GetCell(i-1,j,km);
	Cell &cxpm = GetCell(i+1,j,km);

	WDOUBLE cccW_Jz=ccc.W_Jz + ccc.B_Jz*(1-kzz)+cccm.B_Jz*kzz;
	WDOUBLE cxmW_Jz=cxm.W_Jz + cxm.B_Jz*(1-kzz)+cxmm.B_Jz*kzz;
	WDOUBLE cxpW_Jz=cxp.W_Jz + cxp.B_Jz*(1-kzz)+cxpm.B_Jz*kzz;

	hxp=(cxp.dx+ccc.dx)*0.5;
	hxm=(ccc.dx+cxm.dx)*0.5;
	hxa=hxp+hxm;

	hyp=(cyp.dy+ccc.dy)*0.5;
	hym=(ccc.dy+cym.dy)*0.5;
	hya=hyp+hym;

	Asq  = ccc.W_Asq;
	gamma = 0.5*(1+ccc.W_Psi)*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*Asq)/(1+ccc.W_Psi);

	ax = ((gamma*ccc.W_Ex-0.25*ccc.W_Ponx*ccc.W_Denn)/(1+ccc.W_Psi)-ccc.W_Jx*ccc.W_Ez
		-(ccc.W_Jxx*ccc.W_Ex+ccc.W_Jxy*ccc.W_Ey))/(1+ccc.W_Psi);


	sx = -ccc.W_Jy*ccc.W_Bz/(1+ccc.W_Psi) + ax - ( (cxp.W_Jxx - ccc.W_Jxx - cxpW_Jz + cccW_Jz)*hxm/hxp/hxa + (ccc.W_Jxx - cxm.W_Jxx - cccW_Jz + cxmW_Jz)*hxp/hxm/hxa )
		 -((cyp.W_Jxy - ccc.W_Jxy)*hym/hyp/hya + (ccc.W_Jxy - cym.W_Jxy)*hyp/hym/hya);

	if(k>1 && k< GridZ-1)
	{
		Cell &czp  = GetCell(i, j,k+1);
		Cell &czm  = GetCell(i, j,k-1);
		Cell &cz   = GetCell(i, j,km+1);
		Cell &czmm = GetCell(i, j,km-1);

		sx += ( (czp.B_Jx-czm.B_Jx)*(1-kzz)+(cz.B_Jx-czmm.B_Jx)*kzz )*0.5/dz;
	}
	return sx;
}




//=============================================
//======== Source For Magnetic Field: Sy ======
//=============================================
WDOUBLE Mesh::SourceY(int i, int j, WDOUBLE k0, int k) //v
{

	int km=k-1;
	if(k==0) km=k;
	WDOUBLE kzz=k-k0;
	WDOUBLE hxp,hxm,hxa;
	WDOUBLE hyp,hym,hya;

	WDOUBLE gamma, ay, sy, Asq;

	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	Cell &cccm = GetCell(i,  j,km);
	Cell &cymm = GetCell(i,j-1,km);
	Cell &cypm = GetCell(i,j+1,km);

	WDOUBLE cccW_Jz=ccc.W_Jz + ccc.B_Jz*(1-kzz)+cccm.B_Jz*kzz;
	WDOUBLE cymW_Jz=cym.W_Jz + cym.B_Jz*(1-kzz)+cymm.B_Jz*kzz;
	WDOUBLE cypW_Jz=cyp.W_Jz + cyp.B_Jz*(1-kzz)+cypm.B_Jz*kzz;

	hxp=(cxp.dx+ccc.dx)*0.5;
	hxm=(ccc.dx+cxm.dx)*0.5;
	hxa=hxp+hxm;

	hyp=(cyp.dy+ccc.dy)*0.5;
	hym=(ccc.dy+cym.dy)*0.5;
	hya=hyp+hym;


	Asq  = ccc.W_Asq;
	gamma = 0.5*(1+ccc.W_Psi)*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*Asq)/(1+ccc.W_Psi);

	ay = ((gamma*ccc.W_Ey-0.25*ccc.W_Pony*ccc.W_Denn)/(1+ccc.W_Psi)-ccc.W_Jy*ccc.W_Ez
		-(ccc.W_Jyy*ccc.W_Ey+ccc.W_Jxy*ccc.W_Ex))/(1+ccc.W_Psi);

	sy = ccc.W_Jx*ccc.W_Bz/(1+ccc.W_Psi) + ay - ( (cyp.W_Jyy - ccc.W_Jyy - cypW_Jz + cccW_Jz)*hym/hyp/hya + (ccc.W_Jyy - cym.W_Jyy - cccW_Jz + cymW_Jz)*hyp/hym/hya )
		 -((cxp.W_Jxy - ccc.W_Jxy)*hxm/hxp/hxa + (ccc.W_Jxy - cxm.W_Jxy)*hxp/hxm/hxa);//-(cxp.W_Jxy-cxm.W_Jxy)*0.5/dx;

	if(k>1 && k< GridZ-1)
	{
		Cell &czp = GetCell(i,  j,k+1);
		Cell &czm = GetCell(i,  j,k-1);
		Cell &cz  = GetCell(i,  j,km+1);
		Cell &czmm = GetCell(i, j,km-1);
		sy +=  ( (czp.B_Jy-czm.B_Jy)*(1-kzz)+(cz.B_Jy-czmm.B_Jy)*kzz )*0.5/dz;
	}
	return sy;
}



//======================================
//======== Put Chi: n*/(1+Psi) =========
//======================================
void Mesh::Put_Chi(WDOUBLE k0, int k)//v
{
	int km=k-1;
	if(k==0) km=k;
	WDOUBLE kzz=k-k0;

	int i,j;

	for (j=1; j<=GridY; j++)
	{
		for (i=1; i<=GridX; i++)
		{
			Cell &cc = GetCell(i,j,k);
			Cell &cm = GetCell(i,j,km);
			cc.W_Chi = cc.W_Denn/(1.0+cc.W_Psi) + cc.B_Chi*(1-kzz)+cm.B_Chi*(kzz);
		}
	}


	return;
}


//======================================
//======== Put dPsi/dx into W_Ex =======
//======== Put dPsi/dy into W_Ey =======
//======================================

void Mesh::Partial_Psi(int k)//v
{

	int i,j;
	WDOUBLE hp,hm,ha;

	for (j=1; j<=GridY; j++)
	{
		for (i=1; i<=GridX; i++)
		{
			Cell &ccc = GetCell(i,  j,k);
			Cell &cxm = GetCell(i-1,j,k);
			Cell &cxp = GetCell(i+1,j,k);
			Cell &cym = GetCell(i,j-1,k);
			Cell &cyp = GetCell(i,j+1,k);

			hp=(cxp.dx+ccc.dx)*0.5;
			hm=(ccc.dx+cxm.dx)*0.5;
			ha=hp+hm;

			ccc.W_Ex = (cxp.W_Psi - ccc.W_Psi)*hm/hp/ha + (ccc.W_Psi - cxm.W_Psi)*hp/hm/ha;

			hp=(cyp.dy+ccc.dy)*0.5;
			hm=(ccc.dy+cym.dy)*0.5;
			ha=hp+hm;

			ccc.W_Ey = (cyp.W_Psi - ccc.W_Psi)*hm/hp/ha + (ccc.W_Psi - cym.W_Psi)*hp/hm/ha;

		}

	}

	return;
}



//=============================================
//======== Put d|A|^2/dx into W_Ponx ==========
//======== Put d|A|^2/dy into W_Pony ==========
//=============================================

void Mesh::Pondermotive(int k)//v
{

	int i,j, NF;
	WDOUBLE Asq_xm, Asq_ym, Asq_xp, Asq_yp;
	WDOUBLE Asq_cc;
	WDOUBLE hp,hm,ha;
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
			Asq_xm = Asq_ym = Asq_xp = Asq_yp = Asq_cc = 0.0;
			for (NF=0; NF<NFreqs; NF++)
			{
				Asq_cc += abs(ccc.Acomx[NF])*abs(ccc.Acomx[NF]) + abs(ccc.Acomy[NF])*abs(ccc.Acomy[NF]);

				Asq_xm += abs(cxm.Acomx[NF])*abs(cxm.Acomx[NF]) + abs(cxm.Acomy[NF])*abs(cxm.Acomy[NF]);
				Asq_xp += abs(cxp.Acomx[NF])*abs(cxp.Acomx[NF]) + abs(cxp.Acomy[NF])*abs(cxp.Acomy[NF]);

				Asq_ym += abs(cym.Acomx[NF])*abs(cym.Acomx[NF]) + abs(cym.Acomy[NF])*abs(cym.Acomy[NF]);
				Asq_yp += abs(cyp.Acomx[NF])*abs(cyp.Acomx[NF]) + abs(cyp.Acomy[NF])*abs(cyp.Acomy[NF]);
			}

			ccc.W_Asq = Asq_cc;

			hp=(cxp.dx+ccc.dx)*0.5;
			hm=(ccc.dx+cxm.dx)*0.5;
			ha=hp+hm;
			ccc.W_Ponx = (Asq_xp-Asq_cc)*hm/hp/ha + (Asq_cc-Asq_xm)*hp/hm/ha;
			
			hp=(cyp.dy+ccc.dy)*0.5;
			hm=(ccc.dy+cym.dy)*0.5;
			ha=hp+hm;
			ccc.W_Pony = (Asq_yp-Asq_cc)*hm/hp/ha + (Asq_cc-Asq_ym)*hp/hm/ha;
		}

	}

	for (j=0; j<=GridY+1; j+=GridY+1)
	{
		for (i=0; i<=GridX+1; i+=1)
		{
			Cell &c = GetCell(i, j, k);
			c.W_Asq = 0.0;
			for (NF=0; NF<NFreqs; NF++)
			{
 				c.W_Asq += abs(c.Acomx[NF])*abs(c.Acomx[NF]) + abs(c.Acomy[NF])*abs(c.Acomy[NF]);
			}		
		}

	}


	for (j=0; j<=GridY+1; j+=1)
	{
		for (i=0; i<=GridX+1; i+=GridX+1)
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

void Mesh::AdjustFields(int k) //?
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
		// for (i=0; i<WAK_DIM-WAK_DIM2; i++)
		// { c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}


	if(RankIdx_X != 1 & RankIdx_Y != Ypa)
	{ 	
	//=====mp-corner======================
		Cell c  = GetCell( 0,GridY+1,k);
		Cell c1 = GetCell( 1,GridY+1,k);
		Cell c2 = GetCell( 0,GridY,k);
		Cell c3 = GetCell( 1,GridY,k);
		// for (i=0; i<WAK_DIM-WAK_DIM2; i++)
		// { c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}


	//=====pm-corner======================
	if(RankIdx_X !=Xpa & RankIdx_Y != 1)
	{ 
		Cell c  = GetCell( GridX+1,0,k);
		Cell c1 = GetCell( GridX,  0,k);
		Cell c2 = GetCell( GridX+1,1,k);
		Cell c3 = GetCell( GridX,  1,k);
		// for (i=0; i<WAK_DIM-WAK_DIM2; i++)
		// { c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}

	if(RankIdx_X !=Xpa & RankIdx_Y != Ypa)
	{
	//=====pp-corner======================
		Cell c  = GetCell( GridX+1,GridY+1,k);
		Cell c1 = GetCell( GridX,  GridY+1,k);
		Cell c2 = GetCell( GridX+1,GridY,k);
		Cell c3 = GetCell( GridX,  GridY,k);
		// for (i=0; i<WAK_DIM-WAK_DIM2; i++)
		// { c.W_Fields[i] = c1.W_Fields[i]*0.4+c2.W_Fields[i]*0.4+c3.W_Fields[i]*0.2 ; }
		  c.W_Asq=c1.W_Asq*0.4+c2.W_Asq*0.4+c3.W_Asq*0.2;
	}
	return;

}


dcomplex Mesh::SourceAx(int i, int j, int k, int NF)//v
{
	dcomplex laplace;
	WDOUBLE k0 = p_domain()->OmegaL[NF];
	int n = i+(GridX+2)*j;
	int na= (GridX+2)*(GridY+2)*NF;

	WDOUBLE hxp,hxm,hxa;
	WDOUBLE hyp,hym,hya;
			
	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);


	hxp=(cxp.dx+ccc.dx)*0.5;
	hxm=(ccc.dx+cxm.dx)*0.5;
	hxa=hxp+hxm;

	hyp=(cyp.dy+ccc.dy)*0.5;
	hym=(ccc.dy+cym.dy)*0.5;
	hya=hyp+hym;

	laplace = (hxm*cxp.Acomx[NF]-hxa*ccc.Acomx[NF]+hxp*cxm.Acomx[NF])*2/hxa/hxp/hxm
			 +(hym*cyp.Acomx[NF]-hya*ccc.Acomx[NF]+hyp*cym.Acomx[NF])*2/hya/hyp/hym;
	
	return (ccc.W_Chi+4*ci*k0/dt-4*dA0)*ccc.Acomx[NF]-laplace+4*dAx1[n+na]+4*dAx2[n+na];

}


dcomplex Mesh::SourceAy(int i, int j, int k, int NF)//v
{
	dcomplex laplace;
	WDOUBLE k0 = p_domain()->OmegaL[NF];
	int n = i+(GridX+2)*j;
	int na= (GridX+2)*(GridY+2)*NF;

	WDOUBLE hxp,hxm,hxa;
	WDOUBLE hyp,hym,hya;
			
	Cell &ccc = GetCell(i,  j,k);
	Cell &cxm = GetCell(i-1,j,k);
	Cell &cxp = GetCell(i+1,j,k);
	Cell &cym = GetCell(i,j-1,k);
	Cell &cyp = GetCell(i,j+1,k);

	hxp=(cxp.dx+ccc.dx)*0.5;
	hxm=(ccc.dx+cxm.dx)*0.5;
	hxa=hxp+hxm;

	hyp=(cyp.dy+ccc.dy)*0.5;
	hym=(ccc.dy+cym.dy)*0.5;
	hya=hyp+hym;

	laplace = (hxm*cxp.Acomy[NF]-hxa*ccc.Acomy[NF]+hxp*cxm.Acomy[NF])*2/hxa/hxp/hxm
			 +(hym*cyp.Acomy[NF]-hya*ccc.Acomy[NF]+hyp*cym.Acomy[NF])*2/hya/hyp/hym;
	
	return (ccc.W_Chi+4*ci*k0/dt-4*dA0)*ccc.Acomy[NF]-laplace+4*dAy1[n+na]+4*dAy2[n+na];
}

void Mesh::Put_dA12(int what, int k, int NF)//v
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


void Mesh::SetFieldZeroAfter(int k0)//v
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
void Mesh::AdjustPsi(int k)//v
{

	if(IfAdjustPsi==0) return;
	int i,j;
	WDOUBLE psimax=0.8/dzz; 
	WDOUBLE psimin=-psimax/(1.+psimax);
	for (j=0; j<=GridY+1; j++)
	{
		for (i=0; i<=GridX+1; i++)
		{
			Cell &ccc = GetCell(i, j, k);
			if(ccc.W_Psi<psimin) ccc.W_Psi=psimin;
		}
	}


return;

}




