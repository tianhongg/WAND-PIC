//----------------------------------------------------------------------------------||
//-------------------                multigridC.cpp              -------------------||
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
//---Starting---------           : Mar-14-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||


#include "wand_PIC.h"


// Full Weight Restriction
void MultiGrid::RestrictionC(int send, int rece, int tolayer, int where)
{

	int i,j;
	WDOUBLE wmm, wxm, wmp;
	WDOUBLE wym, wcc, wyp;
	WDOUBLE wpm, wxp, wpp;
	WDOUBLE wa;

switch(where)
{

	case 0:
	for(j=1; j<=LayerGridY[tolayer]; j++)
	{
		for(i=1; i<=LayerGridX[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,tolayer);

			wmm=mgc.p_Res_mm->dx*mgc.p_Res_mm->dy/4;
			wxm=mgc.p_Res_xm->dx*mgc.p_Res_xm->dy/2;
			wmp=mgc.p_Res_mp->dx*mgc.p_Res_mp->dy/4;

			wym=mgc.p_Res_ym->dx*mgc.p_Res_ym->dy/2;
			wcc=mgc.p_Res_cc->dx*mgc.p_Res_cc->dy;
			wyp=mgc.p_Res_yp->dx*mgc.p_Res_yp->dy/2;

			wpm=mgc.p_Res_pm->dx*mgc.p_Res_pm->dy/4;
			wxp=mgc.p_Res_xp->dx*mgc.p_Res_xp->dy/2;
			wpp=mgc.p_Res_pp->dx*mgc.p_Res_pp->dy/4;

			wa=wmm+wxm+wmp + wym+wcc+wyp + wpm+wxp+wpp;

			mgc.C_value[rece] = (
				mgc.p_Res_mm->C_value[send]*wmm + mgc.p_Res_xm->C_value[send]*wxm + mgc.p_Res_mp->C_value[send]*wmp
			  + mgc.p_Res_ym->C_value[send]*wym + mgc.p_Res_cc->C_value[send]*wcc + mgc.p_Res_yp->C_value[send]*wyp	
			  + mgc.p_Res_pm->C_value[send]*wpm + mgc.p_Res_xp->C_value[send]*wxp + mgc.p_Res_pp->C_value[send]*wpp
			)/wa;

		}

	}
	break;

	case 1:
	// for(j=1; j<=BLayerGrid[tolayer]; j++)
	// {
	// 	for(i=1; i<=BLayerGrid[tolayer]; i++)
	// 	{

	// 		MG_Cell &mgc = GetMGBCell(i,j,tolayer);

	// 		mgc.C_value[rece] = 			(mgc.p_Res_cc)->C_value[send]*0.25
	// 		+((mgc.p_Res_xm)->C_value[send]+(mgc.p_Res_xp)->C_value[send]
	// 		 +(mgc.p_Res_ym)->C_value[send]+(mgc.p_Res_yp)->C_value[send])*0.125
	// 		+((mgc.p_Res_mm)->C_value[send]+(mgc.p_Res_mp)->C_value[send]
	// 		 +(mgc.p_Res_pm)->C_value[send]+(mgc.p_Res_pp)->C_value[send])*0.0625;

	// 	}
	// }
	break;


}

	return;
}



void MultiGrid::RestrictionBC(int send, int rece, int tolayer, int where)
{

	int i,j;

switch(where)
{

	case 0:
	for(j=1; j<=LayerGridY[tolayer]; j++)
	{
		for(i=1; i<=LayerGridX[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,tolayer);
			mgc.C_value[rece] = (mgc.p_Res_cc)->C_value[send];
		}

	}
	break;

	case 1:
	for(j=1; j<=BLayerGrid[tolayer]; j++)
	{
		for(i=1; i<=BLayerGrid[tolayer]; i++)
		{
			MG_Cell &mgc = GetMGBCell(i,j,tolayer);
			mgc.C_value[rece] = (mgc.p_Res_cc)->C_value[send];
		}

	}
	break;
}

	return;
}


// Four type prolongation
void MultiGrid::ProlongationC(int send, int rece, int tolayer, int where)
{

	int i,j;
	WDOUBLE dxm, dxp, dym, dyp;

switch(where)
{

	case 0:
	for(j=1; j<=LayerGridY[tolayer]; j++)
	{
		for(i=1; i<=LayerGridX[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,tolayer);

			switch (mgc.Protype)
			{
    			case 0:
    				mgc.C_value[rece] = (mgc.p_Pro_xm)->C_value[send];
					break;

				case 1: 
					dxm=mgc.p_Pro_xm->dx;
					dxp=mgc.p_Pro_xp->dx;
					mgc.C_value[rece] = ((mgc.p_Pro_xm)->C_value[send]*dxp
										+(mgc.p_Pro_xp)->C_value[send]*dxm)/(dxm+dxp);
					break;
					
				case 2: 
					dym=mgc.p_Pro_ym->dy;
					dyp=mgc.p_Pro_yp->dy;
					mgc.C_value[rece] = ((mgc.p_Pro_ym)->C_value[send]*dyp
										+(mgc.p_Pro_yp)->C_value[send]*dym)/(dym+dyp);
					break;

				case 3: 
					dxm=mgc.p_Pro_xm->dx;
					dxp=mgc.p_Pro_xp->dx;
					dym=mgc.p_Pro_ym->dy;
					dyp=mgc.p_Pro_yp->dy;

					mgc.C_value[rece] = ((mgc.p_Pro_xm)->C_value[send]*dxp*dyp
					+(mgc.p_Pro_xp)->C_value[send]*dxm*dyp+(mgc.p_Pro_ym)->C_value[send]*dxp*dym
					+(mgc.p_Pro_yp)->C_value[send]*dxm*dym)/(dxm+dxp)/(dym+dyp);
					break;

			}

		}

	}
	break;

	case 1:
	

	break;



}



	return;
}


//deprecated//May-25-tianhong
void MultiGrid::SendtoBottomC(int what)
{
	
	return;

}


//deprecated//May-25-tianhong
void MultiGrid::BottomSendBackC(int what)
{
	


	return;

}



//deprecated//May-25-tianhong
void MultiGrid::MG_BottomLayerC(int field)
{

	

	return;
}


void MultiGrid::RelaxationC(int field, int layer, int where)
{

	int i,j;
	int nx,ny, amp;

	WDOUBLE hxp,hxm,hxa,hxd,h2x;
	WDOUBLE hyp,hym,hya,hyd,h2y;
	WDOUBLE wmm,wmp,wpm,wpp;

	dcomplex wcc,wxm,wxp,wym,wyp;
	dcomplex d2xS, d2yS, dxS, dyS;
	dcomplex kcc,kxm,kxp,kym,kyp;

	switch(RelaxType)
	{

	// case 0: Hybird Gauss-Siedel Relaxation
	// case 1: Red-Black Relaxation
	// case 2: 
	case 0:
	default:
	nx=LayerGridX[layer];
	ny=LayerGridY[layer];
		// amp=MeshAmplif[layer];

			//Exchange boundary conditions;
			// case 0: Hybird Gauss-Shield Relaxatio
	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			
			MG_Cell &ccc = GetMGCell(i,   j, layer);

			MG_Cell &cxm = GetMGCell(i-1, j, layer);
			MG_Cell &cxp = GetMGCell(i+1, j, layer);
			MG_Cell &cym = GetMGCell(i, j-1, layer);
			MG_Cell &cyp = GetMGCell(i, j+1, layer);

			MG_Cell &cmm = GetMGCell(i-1, j-1, layer);
			MG_Cell &cmp = GetMGCell(i-1, j+1, layer);
			MG_Cell &cpm = GetMGCell(i+1, j-1, layer);
			MG_Cell &cpp = GetMGCell(i+1, j+1, layer);

			hxp=(cxp.dx+ccc.dx)*0.5;
			hxm=(ccc.dx+cxm.dx)*0.5;
			hxa=hxp+hxm; hxd=hxp-hxm;
			h2x=hxa*hxa-hxp*hxm*3;

			hyp=(cyp.dy+ccc.dy)*0.5;
			hym=(ccc.dy+cym.dy)*0.5;
			hya=hyp+hym; hyd=hyp-hym;
			h2y=hya*hya-hyp*hym*3;

			kcc=ccc.C_value[4];
			kxm=cxm.C_value[4];
			kxp=cxp.C_value[4];
			kym=cym.C_value[4];
			kyp=cyp.C_value[4];

			wmm = (h2x + h2y - 2*hxd*hxp - 2*hyd*hyp)/(3*hxa*hxm*hya*hym);
			wmp = (h2x + h2y - 2*hxd*hxp + 2*hyd*hym)/(3*hxa*hxm*hya*hyp);
			wpm = (h2x + h2y + 2*hxd*hxm - 2*hyd*hyp)/(3*hxa*hxp*hya*hym);
			wpp = (h2x + h2y + 2*hxd*hxm + 2*hyd*hym)/(3*hxa*hxp*hya*hyp);

			wxm = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x-2*hxd*hxp)*(2+hym*hyp*kxm) )/(6*hxa*hxm*hym*hyp);
			wxp = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x+2*hxd*hxm)*(2+hym*hyp*kxp) )/(6*hxa*hxp*hym*hyp);
			wym = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y-2*hyd*hyp)*(2+hxm*hxp*kym) )/(6*hxm*hxp*hya*hym);
			wyp = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y+2*hyd*hym)*(2+hxm*hxp*kyp) )/(6*hxm*hxp*hya*hyp);

			wcc = ( h2y*(-2+hxm*hxp*kcc) + h2x*(-2+hym*hyp*kcc) - 2*(hym*hyp*(4+hxd*hxd*kcc) + hxm*hxp*(4+hyd*hyd*kcc)) )/(6*hxm*hxp*hym*hyp);

			dxS = (cxp.C_value[1]-ccc.C_value[1])*hxm/hxp/hxa + (ccc.C_value[1]-cxm.C_value[1])*hxp/hxm/hxa;
			dyS = (cyp.C_value[1]-ccc.C_value[1])*hym/hyp/hya + (ccc.C_value[1]-cym.C_value[1])*hyp/hym/hya;

			d2xS = (hxm*cxp.C_value[1]-hxa*ccc.C_value[1]+hxp*cxm.C_value[1])*2/hxa/hxp/hxm;
			d2yS = (hym*cyp.C_value[1]-hya*ccc.C_value[1]+hyp*cym.C_value[1])*2/hya/hyp/hym;

			ccc.C_value[0]=(1-omega)*ccc.C_value[0]+omega*
			( ccc.C_value[1] + h2x/12*d2xS + h2y/12*d2yS + hxd/3*dxS + hyd/3*dyS
			 - wxm*cxm.C_value[0]- wxp*cxp.C_value[0]- wym*cym.C_value[0]- wyp*cyp.C_value[0]
			 - wmm*cmm.C_value[0]- wpm*cpm.C_value[0]- wmp*cmp.C_value[0]- wpp*cpp.C_value[0])/(wcc-kcc) ;
		}
	}


		break;

	}

	return;
}




void MultiGrid::ResidualC(int field, int layer, int where)
{

	int i,j;
	int nx,ny, amp;


	WDOUBLE hxp,hxm,hxa,hxd,h2x;
	WDOUBLE hyp,hym,hya,hyd,h2y;
	WDOUBLE wmm,wmp,wpm,wpp;

	dcomplex wcc,wxm,wxp,wym,wyp;
	dcomplex d2xS, d2yS, dxS, dyS;
	dcomplex kcc,kxm,kxp,kym,kyp;

	nx=LayerGridX[layer];
	ny=LayerGridY[layer];


	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{	

			MG_Cell &ccc = GetMGCell(i,   j, layer);

			MG_Cell &cxm = GetMGCell(i-1, j, layer);
			MG_Cell &cxp = GetMGCell(i+1, j, layer);
			MG_Cell &cym = GetMGCell(i, j-1, layer);
			MG_Cell &cyp = GetMGCell(i, j+1, layer);

			MG_Cell &cmm = GetMGCell(i-1, j-1, layer);
			MG_Cell &cmp = GetMGCell(i-1, j+1, layer);
			MG_Cell &cpm = GetMGCell(i+1, j-1, layer);
			MG_Cell &cpp = GetMGCell(i+1, j+1, layer);

			hxp=(cxp.dx+ccc.dx)*0.5;
			hxm=(ccc.dx+cxm.dx)*0.5;
			hxa=hxp+hxm; hxd=hxp-hxm;
			h2x=hxa*hxa-hxp*hxm*3;

			hyp=(cyp.dy+ccc.dy)*0.5;
			hym=(ccc.dy+cym.dy)*0.5;
			hya=hyp+hym; hyd=hyp-hym;
			h2y=hya*hya-hyp*hym*3;

			kcc=ccc.C_value[4];
			kxm=cxm.C_value[4];
			kxp=cxp.C_value[4];
			kym=cym.C_value[4];
			kyp=cyp.C_value[4];

			wmm = (h2x + h2y - 2*hxd*hxp - 2*hyd*hyp)/(3*hxa*hxm*hya*hym);
			wmp = (h2x + h2y - 2*hxd*hxp + 2*hyd*hym)/(3*hxa*hxm*hya*hyp);
			wpm = (h2x + h2y + 2*hxd*hxm - 2*hyd*hyp)/(3*hxa*hxp*hya*hym);
			wpp = (h2x + h2y + 2*hxd*hxm + 2*hyd*hym)/(3*hxa*hxp*hya*hyp);

			wxm = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x-2*hxd*hxp)*(2+hym*hyp*kxm) )/(6*hxa*hxm*hym*hyp);
			wxp = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x+2*hxd*hxm)*(2+hym*hyp*kxp) )/(6*hxa*hxp*hym*hyp);
			wym = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y-2*hyd*hyp)*(2+hxm*hxp*kym) )/(6*hxm*hxp*hya*hym);
			wyp = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y+2*hyd*hym)*(2+hxm*hxp*kyp) )/(6*hxm*hxp*hya*hyp);

			wcc = ( h2y*(-2+hxm*hxp*kcc) + h2x*(-2+hym*hyp*kcc) - 2*(hym*hyp*(4+hxd*hxd*kcc) + hxm*hxp*(4+hyd*hyd*kcc)) )/(6*hxm*hxp*hym*hyp);

			dxS = (cxp.C_value[1]-ccc.C_value[1])*hxm/hxp/hxa + (ccc.C_value[1]-cxm.C_value[1])*hxp/hxm/hxa;
			dyS = (cyp.C_value[1]-ccc.C_value[1])*hym/hyp/hya + (ccc.C_value[1]-cym.C_value[1])*hyp/hym/hya;

			d2xS = (hxm*cxp.C_value[1]-hxa*ccc.C_value[1]+hxp*cxm.C_value[1])*2/hxa/hxp/hxm;
			d2yS = (hym*cyp.C_value[1]-hya*ccc.C_value[1]+hyp*cym.C_value[1])*2/hya/hyp/hym;

			ccc.C_value[2]=ccc.C_value[1]-
			(  wxm*cxm.C_value[0]+ wxp*cxp.C_value[0]+ wym*cym.C_value[0]+ wyp*cyp.C_value[0]
			 + wmm*cmm.C_value[0]+ wpm*cpm.C_value[0]+ wmp*cmp.C_value[0]+ wpp*cpp.C_value[0]
			 + (wcc-kcc)*ccc.C_value[0] -h2x/12*d2xS -h2y/12*d2yS - hxd/3*dxS - hyd/3*dyS
			);
		}

	}


	return;
}





void MultiGrid::SetZeroC(int what, int layer, int where)
{
	int i,j;
	int nx,ny;



	nx=LayerGridX[layer];
	ny=LayerGridY[layer];


	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{

			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.C_value[what] = 0.0;

		}

	}	



	return;
}



void MultiGrid::ExchangeC(int what, int layer)
{

	switch(what)
	{
		case MG_Phi:
			p_domain()->p_Com()->DoCommute(COMMU_MG_P_C, layer);
			break;
		case MG_Sou:
			p_domain()->p_Com()->DoCommute(COMMU_MG_S_C, layer);
			break;
		case MG_Res:
			p_domain()->p_Com()->DoCommute(COMMU_MG_R_C, layer);
		case MG_Chi:
			p_domain()->p_Com()->DoCommute(COMMU_MG_C_C, layer);
			break;
	}
	return;
}



void MultiGrid::AddCorrectionC(int layer, int where)
{	

	int i,j;


	int nx,ny;
	nx=LayerGridX[layer];
	ny=LayerGridY[layer];

	// switch(where)
	// {

	// case 0:
	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.C_value[0] +=ccc.C_value[2];
		
		}

	}
	

	return;
}


WDOUBLE MultiGrid::FindErrorC(WDOUBLE &maxall)
{

	int i,j;

	dcomplex epsn;
	WDOUBLE   epsp;
	WDOUBLE   eps;

	WDOUBLE maxp;

	epsp = 0.0;
	maxp=0;


	for (j=1; j<=LayerGridY[1]; j++)
	{
		for (i=1; i<=LayerGridX[1]; i++)
		{
			MG_Cell &ccc = GetMGCell(i, j, 1);
			epsn = (ccc.C_value[0] -ccc.C_value[3]);
			epsp += abs(epsn)*abs(epsn);
			if(abs(ccc.C_value[0])>maxp) maxp=abs(ccc.C_value[0]);
		}
	}

	MPI_Allreduce(&epsp, &eps,    1, MPI_WDOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&maxp, &maxall, 1, MPI_WDOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return eps;
}

void MultiGrid::Put_SourceC(int field, int k, int NF)
{

	int i,j;
	dcomplex dA0=p_Meshs->GetdA0();
	dcomplex kk =4*ci*(p_domain()->OmegaL[NF])/(p_domain()->Get_dt());
	WDOUBLE time =p_domain()->Get_RunTime();

	for (j=0; j<=LayerGridY[1]+1; j++)
	{
		for (i=0; i<=LayerGridX[1]+1; i++)
		{

			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);
			switch(field)
			{	
				case 5:
				mgc.C_value[0] = 2*c.Acomx[NF]-c.Acomxm[NF];	//initial guess
				c.Acomxm[NF] = c.Acomx[NF];
				break;

				case 6:
				mgc.C_value[0] = 2*c.Acomy[NF]-c.Acomym[NF];
				c.Acomym[NF] = c.Acomy[NF];
				break;
			}

		}
	}


	for (j=1; j<=LayerGridY[1]; j++)
	{
		for (i=1; i<=LayerGridX[1]; i++)
		{
			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);
			//================================
			switch(field)
			{	
				case 5:
				mgc.C_value[1] = p_Meshs->SourceAx(i,j,k,NF);
				mgc.C_value[4] = (4*dA0+c.W_Chi-kk);
				break;

				case 6:
				mgc.C_value[1] = p_Meshs->SourceAy(i,j,k,NF);
				mgc.C_value[4] = (4*dA0+c.W_Chi-kk);
				break;
			}

		}

	}

	return;
}



void MultiGrid::Put_FieldsC(int field, int k, int NF)
{

	int i,j;

	int nx=LayerGridX[1];
	int ny=LayerGridY[1];

	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{
			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);

			switch(field)
			{	
				case 5:
				c.Acomx[NF] = mgc.C_value[0];  
				break;
				case 6:
				c.Acomy[NF] = mgc.C_value[0];  
				break;
			} 
				
		}

	}
	return;
}



int MultiGrid::MG_V_cycleC(int field, int k, int NF)
{


	int i,j,n;
	WDOUBLE eps;
	WDOUBLE maxall=1.0;

//============================================================
//==============   Put Source For Different Equation =========
//============================================================
	Put_SourceC(field, k, NF);
	ExchangeC(MG_Sou, 1);
	ExchangeC(MG_Chi, 1);

//============================================================
//==============    Restrict Additional Coefficient  =========
//==============   	 for Helmholtz Equation 		 =========
//============================================================

		for (n=1; n<MPI_Layer; n++)
		{
			RestrictionBC(MG_Chi,MG_Chi,n+1, 0);
			ExchangeC(MG_Chi, n+1);
		}

//============================================================
//==============     V-Cycle Start          ==================
//============================================================


	eps = 100.0;
	int iter = 0;


while(eps > EpsLim*maxall)
{

	iter++;
	//record old value
	for (j=1; j<=LayerGridY[1]; j++)
		for (i=1; i<=LayerGridX[1]; i++)
			GetMGCell(i, j, 1).C_value[3] = GetMGCell(i, j, 1).C_value[0];






//============================================================
//==============     V-Cycle Moving Down    ==================
//============================================================


	for (n=1; n<MPI_Layer; n++)
	{
		//relaxation u1 times
		for(int m=0; m< u1; m++)
		{
			RelaxationC(field, n, 0);
			ExchangeC(MG_Phi, n);
		}

		//find the residual;
   		ResidualC(field, n, 0);
		ExchangeC(MG_Res, n);


		//restrict to next layer;
		RestrictionC(MG_Res,MG_Sou,n+1, 0);
		ExchangeC(MG_Sou, n+1);

		SetZeroC(MG_Phi, n+1, 0);

	}
	
//============================================================
//==============  V-Cycle Layers Below MPI Layer  ============
//============================================================


	switch(BottomType)
	{

		case 0:

			for(int m=0; m< u1; m++)
			{
				RelaxationC(field, n, 0);
				ExchangeC(MG_Phi, MPI_Layer);
			}

		break;

		case 1:

			// SendtoBottomC(MG_Sou);
			// if(Rank == Worker) 
			// {
			// 	MG_BottomLayerC(field);
			// }
			// BottomSendBackC(MG_Phi);

		break;

	}

//============================================================
//==============     V-Cycle Moving Up      ==================
//============================================================


	for (n=MPI_Layer; n>1; n--)
	{

		ProlongationC(MG_Phi, MG_Res, n-1, 0);
		AddCorrectionC(n-1, 0);

		for(int m=0; m<u2; m++)
		{
			ExchangeC(MG_Phi, n-1);
			RelaxationC(field,n-1, 0);
		}

	}
//============================================================
//==============     Error Estimation      ==================
//============================================================
	eps = FindErrorC(maxall)/GridXY;

	//======== test==========
	//if(Rank == 0) std::cout<<"iter: "<<iter<<";   eps:"<<eps<<'\n';
	//======== test==========

//============================================================
//==============     Fail to Converge.      ==================
//============================================================
	if(iter>1e2 || eps>1e5)
	{
		if (Rank==0)  std::cout <<"==== Multigrid: Equation Type: "<<field<< " Failed To Converge at k(z) ="<<k<<'\n';
		return 1;
	}

}


//============================================================
//==============   Put Solution Back to Cell         =========
//============================================================
	Put_FieldsC(field, k, NF);

	return 0;



}



//=========================================================







