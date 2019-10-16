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
// Don't forget take care the corner values at exchange();
void MultiGrid::RestrictionC(int send, int rece, int tolayer, int where)
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

			mgc.C_value[rece] = 			(mgc.p_Res_cc)->C_value[send]*0.25
			+((mgc.p_Res_xm)->C_value[send]+(mgc.p_Res_xp)->C_value[send]
			 +(mgc.p_Res_ym)->C_value[send]+(mgc.p_Res_yp)->C_value[send])*0.125
			+((mgc.p_Res_mm)->C_value[send]+(mgc.p_Res_mp)->C_value[send]
			 +(mgc.p_Res_pm)->C_value[send]+(mgc.p_Res_pp)->C_value[send])*0.0625;

		}

	}
	break;

	case 1:
	for(j=1; j<=BLayerGrid[tolayer]; j++)
	{
		for(i=1; i<=BLayerGrid[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGBCell(i,j,tolayer);

			mgc.C_value[rece] = 			(mgc.p_Res_cc)->C_value[send]*0.25
			+((mgc.p_Res_xm)->C_value[send]+(mgc.p_Res_xp)->C_value[send]
			 +(mgc.p_Res_ym)->C_value[send]+(mgc.p_Res_yp)->C_value[send])*0.125
			+((mgc.p_Res_mm)->C_value[send]+(mgc.p_Res_mp)->C_value[send]
			 +(mgc.p_Res_pm)->C_value[send]+(mgc.p_Res_pp)->C_value[send])*0.0625;

		}
	}
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
					mgc.C_value[rece] = ((mgc.p_Pro_xm)->C_value[send]
										+(mgc.p_Pro_xp)->C_value[send])*0.5;
					break;
					
				case 2: 
					mgc.C_value[rece] = ((mgc.p_Pro_ym)->C_value[send]
										+(mgc.p_Pro_yp)->C_value[send])*0.5;

					break;

				case 3: 
					mgc.C_value[rece] =			  ((mgc.p_Pro_xm)->C_value[send]
					+(mgc.p_Pro_xp)->C_value[send]+(mgc.p_Pro_ym)->C_value[send]
					+(mgc.p_Pro_yp)->C_value[send])*0.25;
					break;

			}

		}

	}


	break;

	case 1:
	for(j=1; j<=BLayerGrid[tolayer]; j++)
	{
		for(i=1; i<=BLayerGrid[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGBCell(i,j,tolayer);

			switch (mgc.Protype)
			{
    			case 0:
    				mgc.C_value[rece] = (mgc.p_Pro_xm)->C_value[send];
					break;

				case 1: 
					mgc.C_value[rece] = ((mgc.p_Pro_xm)->C_value[send]
										+(mgc.p_Pro_xp)->C_value[send])*0.5;
					break;
					
				case 2: 
					mgc.C_value[rece] = ((mgc.p_Pro_ym)->C_value[send]
										+(mgc.p_Pro_yp)->C_value[send])*0.5;

					break;

				case 3: 
					mgc.C_value[rece] =			  ((mgc.p_Pro_xm)->C_value[send]
					+(mgc.p_Pro_xp)->C_value[send]+(mgc.p_Pro_ym)->C_value[send]
					+(mgc.p_Pro_yp)->C_value[send])*0.25;
					break;

			}

		}

	}

	break;



}



	return;
}



void MultiGrid::SendtoBottomC(int what)
{
	int i, j, n;

	int nsend;
	nsend = LayerGridX[MPI_Layer]*LayerGridY[MPI_Layer];
	double mysend[nsend];
	double myreceive[BottomCells];
	double mysend2[nsend];
	double myreceive2[BottomCells];

	n = 0;
	for(j=1; j<=LayerGridY[MPI_Layer]; j++)
	{
		for(i=1; i<=LayerGridX[MPI_Layer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,MPI_Layer);
			mysend[n]  = (mgc.C_value[what]).real();
			mysend2[n] = (mgc.C_value[what]).imag();
			n++;

		}
	}

	MPI_Gatherv(&mysend,  nsend, MPI_DOUBLE, &myreceive,  BottomSend, Bottomdisp,
				MPI_DOUBLE, Worker, MPI_COMM_WORLD);
	MPI_Gatherv(&mysend2, nsend, MPI_DOUBLE, &myreceive2, BottomSend, Bottomdisp,
				MPI_DOUBLE, Worker, MPI_COMM_WORLD);

	//put cells into Worker processor;
	int ii = 0;
	int jj = 0;

	int offx[Xpa];
	int offy[Ypa];

	offx[0] = 0;
	offy[0] = 0;

	for(n=0; n<Xpa*Ypa; n++)
	{
	   if(n%Xpa+1<Xpa)		{offx[n%Xpa+1] 		= recebottomX[n];};
	   if(int(n/Ypa)+1<Ypa) {offy[int(n/Ypa)+1] = recebottomY[n];};
	}


	for(n = 1; n<Xpa; n++)
	{
		offx[n] += offx[n-1];
		offy[n] += offy[n-1];
	}


	if(Rank == Worker)
	{

		for (n=0; n<Xpa*Ypa; n++)
		{

			for(j=1; j<=recebottomY[n]; j++)
			{

				for(i=1; i<=recebottomX[n]; i++)
				{

					ii = offx[ n%Xpa ]+i;
					jj = offy[ int(n/Xpa) ]+j;

					MG_Cell &mgc = GetMGBCell(ii,jj,1);
					mgc.C_value[what] =  myreceive[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)]
								   + ci*myreceive2[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)];
				}

			}

		}

	}


	return;

}



void MultiGrid::BottomSendBackC(int what)
{
	int i, j, n;

	int nrece;
	nrece = LayerGridX[MPI_Layer]*LayerGridY[MPI_Layer];

	double myrece[nrece];
	double mysend[BottomCells];
	double myrece2[nrece];
	double mysend2[BottomCells];


	//put cells into Worker processor;
	int ii = 0;
	int jj = 0;



	int offx[Xpa];
	int offy[Xpa];

	offx[0] = 0;
	offy[0] = 0;

	for(n=0; n<Xpa*Ypa; n++)
	{
	   if(n%Xpa+1<Xpa)		{offx[n%Xpa+1] 		= recebottomX[n];};
	   if(int(n/Ypa)+1<Ypa) {offy[int(n/Ypa)+1] = recebottomY[n];};
	}


	for(n = 1; n<Xpa; n++)
	{
		offx[n] += offx[n-1];
		offy[n] += offy[n-1];
	}

	if(Rank == Worker)
	{

		for (n=0; n<Xpa*Ypa; n++)
		{

			for(j=1; j<=recebottomY[n]; j++)
			{

				for(i=1; i<=recebottomX[n]; i++)
				{

					ii = offx[ n%Xpa ]+i;
					jj = offy[ int(n/Xpa) ]+j;

					MG_Cell &mgc = GetMGBCell(ii,jj,1);
					 mysend[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)] = (mgc.C_value[what]).real();
					mysend2[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)] = (mgc.C_value[what]).imag();
				
				}

			}

		}

	}

	MPI_Scatterv(&mysend,  BottomSend, Bottomdisp, MPI_DOUBLE,  &myrece, nrece,
				MPI_DOUBLE, Worker, MPI_COMM_WORLD);
	MPI_Scatterv(&mysend2, BottomSend, Bottomdisp, MPI_DOUBLE, &myrece2, nrece,
				MPI_DOUBLE, Worker, MPI_COMM_WORLD);


	n = 0;
	for(j=1; j<=LayerGridY[MPI_Layer]; j++)
	{
		for(i=1; i<=LayerGridX[MPI_Layer]; i++)
		{
			MG_Cell &mgc = GetMGCell(i,j,MPI_Layer);
			mgc.C_value[what] = myrece[n]+ci*myrece2[n];
			n++;
		}
	}



	return;

}




void MultiGrid::MG_BottomLayerC(int field)
{

	int n;
	switch(BottomType)
	{
	// case 0: Do not do the Bottom Layer;
	// case 1: Multigrid 
	// case 2: Direct Method


	case 1:

//============================================================
//==============     V-Cycle Moving Down    ==================
//============================================================
	for (n=1; n<SER_Layer; n++)
	{
		//relaxation u1 times
		for(int m=0; m< u1; m++)
		{
			RelaxationC(field,n,1);
		}

		//find the residual;
		ResidualC(field,n,1);

		//restrict to next layer;
		RestrictionC(MG_Res,MG_Sou,n+1,1);

		SetZeroC(MG_Phi, n+1,1);

	}

//============================================================
//==============  V-Cycle Layers Below MPI Layer  ============
//============================================================

	MG_Cell &c1 = GetMGBCell(1, 1, SER_Layer);
	MG_Cell &c2 = GetMGBCell(1, 2, SER_Layer);
	MG_Cell &c3 = GetMGBCell(2, 1, SER_Layer);
	MG_Cell &c4 = GetMGBCell(2, 2, SER_Layer);

	c1.C_value[0] = -(7*c1.C_value[1]+2*c2.C_value[1]+2*c3.C_value[1]+c4.C_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c2.C_value[0] = -(7*c2.C_value[1]+2*c1.C_value[1]+2*c4.C_value[1]+c3.C_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c3.C_value[0] = -(7*c3.C_value[1]+2*c4.C_value[1]+2*c1.C_value[1]+c2.C_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c4.C_value[0] = -(7*c4.C_value[1]+2*c3.C_value[1]+2*c2.C_value[1]+c1.C_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;

//============================================================
//==============     V-Cycle Moving Up      ==================
//============================================================


	for (n=SER_Layer; n>1; n--)
	{

		ProlongationC(MG_Phi, MG_Res, n-1,1);

		AddCorrectionC(n-1,1);

		for(int m=0; m< u2; m++)
		{
			RelaxationC(field,n-1,1);
		}

	}


	break;



	case 2:
	break;


	}


	return;
}


void MultiGrid::RelaxationC(int field, int layer, int where)
{

	int i,j;

	switch(RelaxType)
	{

	// case 0: Hybird Gauss-Siedel Relaxation
	// case 1: Red-Black Relaxation
	// case 2: 

		case 0:
		default:
			
			// case 0: Hybird Gauss-Shield Relaxation
	switch(where)
	{

		case 0:
		//MPI layer
			for (j=1; j<=LayerGridY[layer]; j++)
			{
				for (i=1; i<=LayerGridX[layer]; i++)
				{
					MG_Cell &ccc = GetMGCell(i,   j, layer);
					MG_Cell &cxm = GetMGCell(i-1, j, layer);
					MG_Cell &cxp = GetMGCell(i+1, j, layer);
					MG_Cell &cym = GetMGCell(i, j-1, layer);
					MG_Cell &cyp = GetMGCell(i, j+1, layer);

					ccc.C_value[0]=(1-omega)*ccc.C_value[0]+omega*(cxm.C_value[0]+cxp.C_value[0]
					+cym.C_value[0]+cyp.C_value[0]-ccc.C_value[1]*MeshAmplif[layer])/(4.0+ccc.C_value[4]*MeshAmplif[layer]);

				}
			}
	
		break;


		case 1:
		//Bottom layer
			for (j=1; j<=BLayerGrid[layer]; j++)
			{
				for (i=1; i<=BLayerGrid[layer]; i++)
				{
					MG_Cell &ccc = GetMGBCell(i,   j, layer);
					MG_Cell &cxm = GetMGBCell(i-1, j, layer);
					MG_Cell &cxp = GetMGBCell(i+1, j, layer);
					MG_Cell &cym = GetMGBCell(i, j-1, layer);
					MG_Cell &cyp = GetMGBCell(i, j+1, layer);

					ccc.C_value[0]=(1-omega)*ccc.C_value[0]+omega*(cxm.C_value[0]+cxp.C_value[0]
					+cym.C_value[0]+cyp.C_value[0]-ccc.C_value[1]*MeshAmplif[MPI_Layer-1+layer])
					/(4.0+ccc.C_value[4]*MeshAmplif[MPI_Layer-1+layer]);

				}

			}

		break;

		}


		break;


	}

	return;
}




void MultiGrid::ResidualC(int field, int layer, int where)
{

	int i,j;

	switch(where)
	{

	case 0:
	for (j=1; j<=LayerGridY[layer]; j++)
	{
		for (i=1; i<=LayerGridX[layer]; i++)
		{
			MG_Cell &ccc = GetMGCell(i,   j, layer);
			MG_Cell &cxm = GetMGCell(i-1, j, layer);
			MG_Cell &cxp = GetMGCell(i+1, j, layer);
			MG_Cell &cym = GetMGCell(i, j-1, layer);
			MG_Cell &cyp = GetMGCell(i, j+1, layer);


			ccc.C_value[2]=ccc.C_value[1]-(cxm.C_value[0]+cxp.C_value[0]
			+cym.C_value[0]+cyp.C_value[0]-(4+ccc.C_value[4]*MeshAmplif[layer])*ccc.C_value[0])/MeshAmplif[layer];
		}

	}

	break;

	case 1:

	for (j=1; j<=BLayerGrid[layer]; j++)
	{
		for (i=1; i<=BLayerGrid[layer]; i++)
		{
			MG_Cell &ccc = GetMGBCell(i,   j, layer);
			MG_Cell &cxm = GetMGBCell(i-1, j, layer);
			MG_Cell &cxp = GetMGBCell(i+1, j, layer);
			MG_Cell &cym = GetMGBCell(i, j-1, layer);
			MG_Cell &cyp = GetMGBCell(i, j+1, layer);

			 ccc.C_value[2]=ccc.C_value[1]-(cxm.C_value[0]+cxp.C_value[0]
			+cym.C_value[0]+cyp.C_value[0]-(4+ccc.C_value[4]*MeshAmplif[MPI_Layer-1+layer])*ccc.C_value[0])/MeshAmplif[MPI_Layer-1+layer];

		}

	}

	break;

}

	return;
}





void MultiGrid::SetZeroC(int what, int layer, int where)
{
	int i,j;


	switch(where)
	{

	case 0:

	for (j=0; j<=LayerGridY[layer]+1; j++)
	{
		for (i=0; i<=LayerGridX[layer]+1; i++)
		{

			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.C_value[what] = 0;

		}

	}	

	break;

	case 1:

	for (j=0; j<=BLayerGrid[layer]+1; j++)
	{
		for (i=0; i<=BLayerGrid[layer]+1; i++)
		{

			MG_Cell &ccc = GetMGBCell(i, j, layer);
			ccc.C_value[what] = 0;

		}

	}	

	break;


}

	return;
}



void MultiGrid::ExchangeC(int what, int layer)
{

	switch(what)
	{
		case 0:
			p_domain()->p_Com()->DoCommute(COMMU_MG_P_C, layer);
			break;
		case 1:
			p_domain()->p_Com()->DoCommute(COMMU_MG_S_C, layer);
			break;
		case 2:
			p_domain()->p_Com()->DoCommute(COMMU_MG_R_C, layer);
			break;

	}


	//adjust corner cell;
	if(what == 2)
	{

		MG_Cell &cmm = GetMGCell(0, 0, layer);
		cmm.C_value[what] = (GetMGCell(0, 1, layer).C_value[what] 
							+GetMGCell(1, 0, layer).C_value[what])*0.5;

		MG_Cell &cpm = GetMGCell(LayerGridX[layer]+1, 0, layer);
		cpm.C_value[what] = (GetMGCell(LayerGridX[layer]+1, 1, layer).C_value[what] 
							+GetMGCell(LayerGridX[layer],   0, layer).C_value[what])*0.5;

		MG_Cell &cmp = GetMGCell(0, LayerGridY[layer]+1, layer);
		cmp.C_value[what] = (GetMGCell(0, LayerGridY[layer],   layer).C_value[what] 
							+GetMGCell(1, LayerGridY[layer]+1, layer).C_value[what])*0.5;

		MG_Cell &cpp = GetMGCell(LayerGridX[layer]+1, LayerGridY[layer]+1, layer);
		cpp.C_value[what] = (GetMGCell(LayerGridX[layer],   LayerGridY[layer]+1, layer).C_value[what] 
							+GetMGCell(LayerGridX[layer]+1, LayerGridY[layer],   layer).C_value[what])*0.5;
	}

	return;
}



void MultiGrid::AddCorrectionC(int layer, int where)
{	

	int i,j;


	switch(where)
	{

	case 0:

	for (j=1; j<=LayerGridY[layer]; j++)
	{
		for (i=1; i<=LayerGridX[layer]; i++)
		{
			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.C_value[0] +=ccc.C_value[2];
		
		}

	}

	break;

	case 1:
	for (j=1; j<=BLayerGrid[layer]; j++)
	{
		for (i=1; i<=BLayerGrid[layer]; i++)
		{
			MG_Cell &ccc = GetMGBCell(i, j, layer);
			ccc.C_value[0] +=ccc.C_value[2];
		
		}

	}

	break;


}



	return;
}


double MultiGrid::FindErrorC(double &maxall)
{

	int i,j;

	dcomplex epsn;
	double   epsp;
	double   eps;

	double maxp;

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

	MPI_Allreduce(&epsp, &eps,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&maxp, &maxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return eps;
}

void MultiGrid::Put_SourceC(int field, int k, int NF)
{

	int i,j;
	dcomplex dA0=p_Meshs->GetdA0();
	dcomplex kk =4*ci*(p_domain()->OmegaL[NF])/(p_domain()->Get_dt());
	
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
				mgc.C_value[1] = p_Meshs->SourceAx(i,j,k,NF)*dxdy;
				mgc.C_value[4] = (4*dA0+c.W_Chi-kk)*dxdy;
				break;

				case 6:
				mgc.C_value[1] = p_Meshs->SourceAy(i,j,k,NF)*dxdy;
				mgc.C_value[4] = (4*dA0+c.W_Chi-kk)*dxdy;
				break;
			}

		}

	}

	return;
}



void MultiGrid::Put_FieldsC(int field, int k, int NF)
{

	int i,j;

	for (j=0; j<=LayerGridY[1]+1; j++)
	{
		for (i=0; i<=LayerGridX[1]+1; i++)
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
	double eps;
	double maxall=1.0;

//============================================================
//==============   Put Source For Different Equation =========
//============================================================
	Put_SourceC(field, k, NF);

//============================================================
//==============    Restrict Additional Coefficient  =========
//==============   	 for Helmholtz Equation 		 =========
//============================================================

		for (n=1; n<MPI_Layer; n++)
		{
			RestrictionBC(MG_Chi,MG_Chi,n+1, 0);
		}
		switch(BottomType)
		{
			case 1:
				SendtoBottomC(MG_Chi);
				for (n=1; n<SER_Layer; n++)
				{
					RestrictionBC(MG_Chi,MG_Chi,n+1,1);
				}
				break;
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
	{
		for (i=1; i<=LayerGridX[1]; i++)
		{
			GetMGCell(i, j, 1).C_value[3] = GetMGCell(i, j, 1).C_value[0];

		}

	}


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

			SendtoBottomC(MG_Sou);
			if(Rank == Worker) 
			{
				MG_BottomLayerC(field);
			}
			BottomSendBackC(MG_Phi);

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







