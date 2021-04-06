//----------------------------------------------------------------------------------||
//-------------------                multigrid.cpp               -------------------||
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
//---Starting---------           : Feb-07-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||


#include "wand_PIC.h"


// Initiate at the beginning of the simulation 
// Define the rule of Restriction and Prolongation
// Define the way of Restriction and Prolongation
// Define the 2D Cell Used for Multigrid Relaxation


//=======Multigrid Cell Class

MG_Cell::MG_Cell()
{
	  field  =   source =   Residu =   savefield =   Chi = 0.0;
	C_field  = C_source = C_Residu = C_savefield = C_Chi = 0.0;

	p_Res_cc = p_Res_xm = p_Res_xp = p_Res_ym = p_Res_yp = NULL;
	p_Res_mm = p_Res_mp = p_Res_pm = p_Res_pp = NULL;
	p_Pro_xm = p_Pro_xp = p_Pro_ym = p_Pro_yp = NULL;

	Protype = 0;
};




MultiGrid::MultiGrid(int rank, int XGridN, int YGridN, FILE *f):NList("MultiGrid")
{


	p_MGBCell = p_MGCell = NULL;
	BottomSend = Bottomdisp = NULL;
	recebottomX = recebottomY = NULL;

	p_Meshs = p_domain()->p_Mesh();

	AddEntry((char*)"OrderDown", 	&u1, 3);
	AddEntry((char*)"OrderUp", 	 	&u2, 3);
	AddEntry((char*)"RelaxMethod",	&RelaxType, 0);
	AddEntry((char*)"BottomMethod", &BottomType, 0);
	AddEntry((char*)"SOR_Omega", 	&omega, 1.5);
	AddEntry((char*)"ErrorLimit", 	&EpsLim, 1e-6);
	

	  if (f)
    {
      rewind(f);
      read(f);
    }


   
    Rank = rank;

	int i,j,n;

	//  Partition;
	int RankIdx_X = p_domain()->p_Partition()->RankIdx_X(); 
	int RankIdx_Y = p_domain()->p_Partition()->RankIdx_Y(); 

	Xpa = p_domain()->p_Partition()->GetXpart();
	Ypa = p_domain()->p_Partition()->GetYpart();

	dxdy = (p_domain()-> Get_dx())*(p_domain()-> Get_dy());

	GridXY = XGridN*YGridN*Xpa*Ypa;

	// initialize
	for (i=0; i<Layer_buf; i++)
	{
    LayerGridX[i] = 0;
    LayerGridY[i] = 0;
	}
	// Layer index starts from 1;

	LayerGridX[1] = XGridN;
	LayerGridY[1] = YGridN;


// ---------------------------------------------------------
	// Define the Restriction rule and layers in X direction
	// by a very simple and stupid way: walking and jumping.

// ---------------------------------------------------------

// ---------------------------------------------------------
	// define layer and grid number at each layer//
// ---------------------------------------------------------

//==========================================================//
//================= X Direction ============================//
	int gridmin = 2;
	MPI_Layer = 1;

	while(gridmin>1)
	{
    	i = 1;
    	MPI_Layer += 1;
    	while(i<=XGridN*Xpa)
    		{
    			if(i >= ((RankIdx_X-1)*XGridN+1) & i<= (RankIdx_X)*XGridN )
    			{
					LayerGridX[MPI_Layer] += 1;
        		}
				i += pow(2,(MPI_Layer-1));
			}

			MPI_Allreduce(&LayerGridX[MPI_Layer], &gridmin, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
	}
//================= X Direction ============================//
//==========================================================//




//==========================================================//
//================= Y Direction ============================//

	gridmin = 2;
	MPI_Layer = 1;

	while(gridmin>1)
	{
    	i = 1;
    	MPI_Layer += 1;
    	while(i<=YGridN*Ypa)
    		{
    			if(i >= ((RankIdx_Y-1)*YGridN+1) & i<= (RankIdx_Y)*YGridN )
    			{
					LayerGridY[MPI_Layer] += 1;
        		}
				i += pow(2,(MPI_Layer-1));
			}

			MPI_Allreduce(&LayerGridY[MPI_Layer], &gridmin, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
	}
//================= Y Direction ============================//
//==========================================================//


	int Tot_GridX = 0;
	int Tot_GridY = 0;

	for (i=1; i<MPI_Layer+1; i++)
	{
    	Tot_GridX += LayerGridX[i];
    	Tot_GridY += LayerGridY[i];
	}
	Tot_GridX += MPI_Layer*2;
	Tot_GridY += MPI_Layer*2;


// ---------------------------------------------------------
	// Find the 2D Cells Sum of all Layers above this Layer
// ---------------------------------------------------------
	Layer_Sum[1] =0;
	for (i=2; i<MPI_Layer+1; i++)
	{
		Layer_Sum[i] =Layer_Sum[i-1]+(LayerGridX[i-1]+2)*(LayerGridY[i-1]+2);
	}

// ---------------------------------------------------------
	// Find the 1D Cell Sum of all Layers above this Layer
// ---------------------------------------------------------
	int LayerSum1dX[MPI_Layer+1];
	int LayerSum1dY[MPI_Layer+1];

	LayerSum1dX[1] =0;
	LayerSum1dY[1] =0;

	for (i=2; i<=MPI_Layer; i++)
	{
		LayerSum1dX[i] =LayerSum1dX[i-1]+(LayerGridX[i-1]+2);
		LayerSum1dY[i] =LayerSum1dY[i-1]+(LayerGridY[i-1]+2);
	}


// ----------------------------------------------------------------
	// Find the Cell Grand_Index in X(or Y) direction at each layer
// ----------------------------------------------------------------
	int Cellidx_i[Tot_GridX];
	int Cellidx_j[Tot_GridY];

	int Celln = XGridN+2;

	for (i=0; i<Celln; i++) 
	{
		Cellidx_i[i] = i+((RankIdx_X-1)*XGridN);
		Cellidx_j[i] = i+((RankIdx_Y-1)*YGridN);
	}


	for(n=2;n<=MPI_Layer;n++)
	{
		i = 1;
		while(i<=XGridN*Xpa)
		{
			if(i >= ((RankIdx_X-1)*XGridN+1) & i<= (RankIdx_X)*XGridN )
			{
				Celln += 1;
				Cellidx_i[Celln] = i;
			}
		i += pow(2,(n-1));
		}

		Cellidx_i[Celln+1] = Cellidx_i[Celln]+pow(2,(n-1));
		Cellidx_i[Celln-LayerGridX[n]] = Cellidx_i[Celln-LayerGridX[n]+1]-pow(2,(n-1));
		Celln += 2;

	}

	Celln = YGridN+2;

	for(n=2;n<=MPI_Layer;n++)
	{
		i = 1;
		while(i<=YGridN*Ypa)
		{
			if(i >= ((RankIdx_Y-1)*YGridN+1) & i<= (RankIdx_Y)*YGridN )
			{
				Celln += 1;
				Cellidx_j[Celln] = i;
			}
		i += pow(2,(n-1));
		}


		Cellidx_j[Celln+1] = Cellidx_j[Celln]+pow(2,(n-1));
		Cellidx_j[Celln-LayerGridY[n]] = Cellidx_j[Celln-LayerGridY[n]+1]-pow(2,(n-1));
		Celln += 2;
	}




// ----------------------------------------------------------------
	// Creat All 2D Cells at All Layers
// ----------------------------------------------------------------

	int Tot_MGCell = 0;

	for (n=1; n<=MPI_Layer; n++)
	{
    	Tot_MGCell += (LayerGridX[n]+2)*(LayerGridY[n]+2);
	}

	p_MGCell = new MG_Cell[Tot_MGCell];

	for (n=1; n<=MPI_Layer; n++)
	{
		for (j=0; j<=LayerGridY[n]+1; j++)
		{
			for (i=0; i<=LayerGridX[n]+1; i++)
			{
				MG_Cell &mgc = GetMGCell(i,j,n);
				mgc.Grandidx_X = Cellidx_i[LayerSum1dX[n]+i];
				mgc.Grandidx_Y = Cellidx_j[LayerSum1dY[n]+j];
			}
		}
	}


// ----------------------------------------------------------------
	// Relate Cells at different layers
	// Find	restriction rules and prolongation rules
// ----------------------------------------------------------------

	int ii, jj;

	// Restriction starts at the second layer//
	for (n=2; n<=MPI_Layer;n++)
	{
		for (j=1; j<=LayerGridY[n]; j++)
		{
			for (i=1; i<=LayerGridX[n]; i++)
			{

				MG_Cell &mgc = GetMGCell(i,j,n);

				ii = 1;
				jj = 1;

				while (Cellidx_i[LayerSum1dX[n-1]+ii]<Cellidx_i[LayerSum1dX[n]+i])
				{ ii++; }
				while (Cellidx_j[LayerSum1dY[n-1]+jj]<Cellidx_j[LayerSum1dY[n]+j])
				{ jj++; }

				mgc.p_Res_cc = &GetMGCell(ii,  jj,n-1);
				mgc.p_Res_xm = &GetMGCell(ii-1,jj,n-1);
				mgc.p_Res_xp = &GetMGCell(ii+1,jj,n-1);
				mgc.p_Res_ym = &GetMGCell(ii,jj-1,n-1);
				mgc.p_Res_yp = &GetMGCell(ii,jj+1,n-1);

				mgc.p_Res_mm = &GetMGCell(ii-1,jj-1,n-1);
				mgc.p_Res_mp = &GetMGCell(ii-1,jj+1,n-1);
				mgc.p_Res_pm = &GetMGCell(ii+1,jj-1,n-1);
				mgc.p_Res_pp = &GetMGCell(ii+1,jj+1,n-1);
			}
		}
	}




	// Prolongation starts at the layer above last layer//
	int findx; 
	int findy;
	for (n=1; n<MPI_Layer;n++)
	{
		for (j=1; j<=LayerGridY[n]; j++)
		{
			for (i=1; i<=LayerGridX[n]; i++)
			{

				MG_Cell &mgc = GetMGCell(i,j,n);
				findx = 0; findy = 0;


				for (ii = 0; ii<=LayerGridX[n+1]+1; ii++)
				{
					if(Cellidx_i[LayerSum1dX[n+1]+ii]==Cellidx_i[LayerSum1dX[n]+i])
					{
						findx = 1;
						break;
					}
					

					if(Cellidx_i[LayerSum1dX[n+1]+ii]<Cellidx_i[LayerSum1dX[n]+i]&
					 Cellidx_i[LayerSum1dX[n+1]+ii+1]>Cellidx_i[LayerSum1dX[n]+i])
					{
						break;
					}


				} 


				for (jj = 0; jj<=LayerGridY[n+1]+1; jj++)
				{
					if(Cellidx_j[LayerSum1dY[n+1]+jj]==Cellidx_j[LayerSum1dY[n]+j])
					{
						findy = 1;
						break;
					}

					if(Cellidx_j[LayerSum1dY[n+1]+jj]<Cellidx_j[LayerSum1dY[n]+j]&
					 Cellidx_j[LayerSum1dY[n+1]+jj+1]>Cellidx_j[LayerSum1dY[n]+j])
					{
						break;
					}

				}

				//type 0 prolongation
				if(findx==1 &findy == 1)
				{
					mgc.Protype = 0;
					mgc.p_Pro_xm = &GetMGCell(ii,jj,n+1);

				}

				//type 1 prolongation
				if(findx==0 &findy == 1)
				{
					mgc.Protype = 1;
					mgc.p_Pro_xm = &GetMGCell(ii,  jj,n+1);
					mgc.p_Pro_xp = &GetMGCell(ii+1,jj,n+1);
				}
				//type 2 prolongation
				if(findx==1 &findy == 0)
				{
					mgc.Protype = 2;
					mgc.p_Pro_ym = &GetMGCell(ii,  jj,n+1);
					mgc.p_Pro_yp = &GetMGCell(ii,jj+1,n+1);
				}
				//type 3 prolongation
				if(findx==0 &findy == 0)
				{
					mgc.Protype = 3;
					mgc.p_Pro_xm = &GetMGCell(ii,    jj,n+1);
					mgc.p_Pro_xp = &GetMGCell(ii+1,  jj,n+1);
					mgc.p_Pro_ym = &GetMGCell(ii,  jj+1,n+1);
					mgc.p_Pro_yp = &GetMGCell(ii+1,jj+1,n+1);
				}

			
			}
		}
	}





// ----------------------------------------------------------------
//           Layers Below the bottom Layer
// ----------------------------------------------------------------
// ----------------------------------------------------------------
//           Layers Below the bottom Layer
// ----------------------------------------------------------------
// ----------------------------------------------------------------
//           Layers Below the bottom Layer
// ----------------------------------------------------------------
// ----------------------------------------------------------------
//           Layers Below the bottom Layer
// ----------------------------------------------------------------

	// collect botton grid number at each processor;

	SER_Layer = 0;
	Worker = floor(Ypa*Xpa/2-1);

	int sendx = LayerGridX[MPI_Layer];  //1d grid at bottom layer at different rank
	int sendy = LayerGridY[MPI_Layer]; 

	recebottomX = new int[Xpa*Ypa];
	recebottomY = new int[Xpa*Ypa];

	MPI_Gather(&sendx, 1, MPI_INT, recebottomX, 1, MPI_INT, Worker, MPI_COMM_WORLD);
 	MPI_Gather(&sendy, 1, MPI_INT, recebottomY, 1, MPI_INT, Worker, MPI_COMM_WORLD);




    // MPI_GatherV receive displacement for each sending processor;
    if(Rank == Worker)
	{	
    	Bottomdisp =  new int[Xpa*Ypa];
    	BottomSend =  new int[Xpa*Ypa];

    	Bottomdisp[0] = 0;
    	for(i=1; i<Xpa*Ypa; i++)
    	{
    		Bottomdisp[i] = Bottomdisp[i-1] + recebottomX[i-1]*recebottomY[i-1];
    	}



    	//total cells at the bottom layer without the bound cells.
    	BottomCells = Bottomdisp[Xpa*Ypa-1]+recebottomX[Xpa*Ypa-1]*recebottomY[Xpa*Ypa-1];
    
    	for(i=0; i<Xpa*Ypa; i++)
    	{
    		BottomSend[i] = recebottomX[i]*recebottomY[i];
    	}


		int BottomGrid = 0;
		for(i=0; i<Xpa; i++)
		{
    		BottomGrid+= recebottomX[i];
		}

		switch(BottomType)
		{
//====================multigrid===========
			case 1: 

			n = BottomGrid;

			SER_Layer = 1;
			BLayerGrid[SER_Layer] = BottomGrid;

			while (n>2)
			{

				SER_Layer++;
				n = ceil(n*0.5);
				BLayerGrid[SER_Layer] = n;

			}


			// tot grid(1d) at bottom layers alltogether
			int Tot_Grid = 0;

			for (n=1; n<=SER_Layer; n++)
			{
    			Tot_Grid += BLayerGrid[n];
			}

			Tot_Grid += SER_Layer*2;



			int  BCellidx_i[Tot_Grid];

			n = 0;
			for (j = 1; j<=SER_Layer; j++)
			{
				for (i = 0; i<=BLayerGrid[j]+1; i++)
				{
					BCellidx_i[n] = (i-1)*pow(2,j-1)+1;
					n++;
				}
				
			}


			BLayer_Sum[1] =0;
			for (i=2; i<SER_Layer+1; i++)
			{
				BLayer_Sum[i] =BLayer_Sum[i-1]+(BLayerGrid[i-1]+2)*(BLayerGrid[i-1]+2);
			}


			int BLayerSum1d[SER_Layer+1];
			BLayerSum1d[1] =0;
			for (i=2; i<SER_Layer+1; i++)
			{
				BLayerSum1d[i] =BLayerSum1d[i-1]+(BLayerGrid[i-1]+2);
			}


			Tot_MGCell = 0;

			for (i=1; i<SER_Layer+1; i++)
			{
    			Tot_MGCell += (BLayerGrid[i]+2)*(BLayerGrid[i]+2);
			}


			p_MGBCell = new MG_Cell[Tot_MGCell];




			for (n=1; n<=SER_Layer;n++)
			{
				for (j=0; j<=BLayerGrid[n]+1; j++)
				{
					for (i=0; i<=BLayerGrid[n]+1; i++)
					{
						MG_Cell &mgc = GetMGBCell(i,j,n);
						mgc.Grandidx_X = BCellidx_i[BLayerSum1d[n]+i];
						mgc.Grandidx_Y = BCellidx_i[BLayerSum1d[n]+j];
					}
				}
			}


// ----------------------------------------------------------------
	// Relate Cells at different layers
	// Find	restriction rules and prolongation rules
// ----------------------------------------------------------------

			// Restriction starts at the second layer//
			for (n=2; n<=SER_Layer;n++)
			{
				for (j=1; j<=BLayerGrid[n]; j++)
				{
					for (i=1; i<=BLayerGrid[n]; i++)
					{

						MG_Cell &mgc = GetMGBCell(i,j,n);

						ii = 1;
						jj = 1;

						while (BCellidx_i[BLayerSum1d[n-1]+ii]<BCellidx_i[BLayerSum1d[n]+i])
						{ ii++; }
						while (BCellidx_i[BLayerSum1d[n-1]+jj]<BCellidx_i[BLayerSum1d[n]+j])
						{ jj++; }

						mgc.p_Res_cc = &GetMGBCell(ii,  jj,n-1);
						mgc.p_Res_xm = &GetMGBCell(ii-1,jj,n-1);
						mgc.p_Res_xp = &GetMGBCell(ii+1,jj,n-1);
						mgc.p_Res_ym = &GetMGBCell(ii,jj-1,n-1);
						mgc.p_Res_yp = &GetMGBCell(ii,jj+1,n-1);

						mgc.p_Res_mm = &GetMGBCell(ii-1,jj-1,n-1);
						mgc.p_Res_mp = &GetMGBCell(ii-1,jj+1,n-1);
						mgc.p_Res_pm = &GetMGBCell(ii+1,jj-1,n-1);
						mgc.p_Res_pp = &GetMGBCell(ii+1,jj+1,n-1);
					}
				}
			}





	// Prolongation starts at the layer above last layer//

			for (n=1; n<SER_Layer;n++)
			{
				for (j=1; j<=BLayerGrid[n]; j++)
				{
					for (i=1; i<=BLayerGrid[n]; i++)
					{

						MG_Cell &mgc = GetMGBCell(i,j,n);
						findx = 0; findy = 0;


						for (ii = 0; ii<=BLayerGrid[n+1]+1; ii++)
						{
							if(BCellidx_i[BLayerSum1d[n+1]+ii]==BCellidx_i[BLayerSum1d[n]+i])
							{
								findx = 1;
								break;
							}
					

							if(BCellidx_i[BLayerSum1d[n+1]+ii]<BCellidx_i[BLayerSum1d[n]+i]&
							 BCellidx_i[BLayerSum1d[n+1]+ii+1]>BCellidx_i[BLayerSum1d[n]+i])
							{
								break;
							}


						} 


						for (jj = 0; jj<=BLayerGrid[n+1]+1; jj++)
						{
							if(BCellidx_i[BLayerSum1d[n+1]+jj]==BCellidx_i[BLayerSum1d[n]+j])
							{
								findy = 1;
								break;
							}

							if(BCellidx_i[BLayerSum1d[n+1]+jj]<BCellidx_i[BLayerSum1d[n]+j]&
							 BCellidx_i[BLayerSum1d[n+1]+jj+1]>BCellidx_i[BLayerSum1d[n]+j])
							{
								break;
							}

						}

						//type 0 prolongation
						if(findx==1 &findy == 1)
						{
							mgc.Protype = 0;
							mgc.p_Pro_xm = &GetMGBCell(ii,jj,n+1);

						}

						//type 1 prolongation
						if(findx==0 &findy == 1)
						{
							mgc.Protype = 1;
							mgc.p_Pro_xm = &GetMGBCell(ii,  jj,n+1);
							mgc.p_Pro_xp = &GetMGBCell(ii+1,jj,n+1);
						}
						//type 2 prolongation
						if(findx==1 &findy == 0)
						{
							mgc.Protype = 2;
							mgc.p_Pro_ym = &GetMGBCell(ii,  jj,n+1);
							mgc.p_Pro_yp = &GetMGBCell(ii,jj+1,n+1);
						}
						//type 3 prolongation
						if(findx==0 &findy == 0)
						{
							mgc.Protype = 3;
							mgc.p_Pro_xm = &GetMGBCell(ii,    jj,n+1);
							mgc.p_Pro_xp = &GetMGBCell(ii+1,  jj,n+1);
							mgc.p_Pro_ym = &GetMGBCell(ii,  jj+1,n+1);
							mgc.p_Pro_yp = &GetMGBCell(ii+1,jj+1,n+1);
						}

			
					}
				}
			}
			//====================multigrid===========
		
			break;

		}

	}

	//==== end worker rank=====//

	for (i=1; i<=MPI_Layer+SER_Layer; i++)
	{
    	MeshAmplif[i] = pow(2,i-1)*pow(2,i-1);
	}

	if (Rank==0)  std::cout << "==== Multigrid: Multigrid Meshes Done.   ====\n";


};
	




// Full Weight Restriction
// Don't forget take care the corner values at exchange();
void MultiGrid::Restriction(int send, int rece, int tolayer, int where)
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

			mgc.M_value[rece] = 			(mgc.p_Res_cc)->M_value[send]*0.25
			+((mgc.p_Res_xm)->M_value[send]+(mgc.p_Res_xp)->M_value[send]
			 +(mgc.p_Res_ym)->M_value[send]+(mgc.p_Res_yp)->M_value[send])*0.125
			+((mgc.p_Res_mm)->M_value[send]+(mgc.p_Res_mp)->M_value[send]
			 +(mgc.p_Res_pm)->M_value[send]+(mgc.p_Res_pp)->M_value[send])*0.0625;

		}

	}
	break;

	case 1:
	for(j=1; j<=BLayerGrid[tolayer]; j++)
	{
		for(i=1; i<=BLayerGrid[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGBCell(i,j,tolayer);

			mgc.M_value[rece] = 			(mgc.p_Res_cc)->M_value[send]*0.25
			+((mgc.p_Res_xm)->M_value[send]+(mgc.p_Res_xp)->M_value[send]
			 +(mgc.p_Res_ym)->M_value[send]+(mgc.p_Res_yp)->M_value[send])*0.125
			+((mgc.p_Res_mm)->M_value[send]+(mgc.p_Res_mp)->M_value[send]
			 +(mgc.p_Res_pm)->M_value[send]+(mgc.p_Res_pp)->M_value[send])*0.0625;

		}

	}
	break;

}

	return;
}



void MultiGrid::RestrictionB(int send, int rece, int tolayer, int where)
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
			mgc.M_value[rece] = (mgc.p_Res_cc)->M_value[send];
			
		}

	}

	break;

	case 1:
	for(j=1; j<=BLayerGrid[tolayer]; j++)
	{
		for(i=1; i<=BLayerGrid[tolayer]; i++)
		{

			MG_Cell &mgc = GetMGBCell(i,j,tolayer);

			mgc.M_value[rece] = (mgc.p_Res_cc)->M_value[send];


		}

	}
	break;


}

	return;
}


// Four type prolongation
void MultiGrid::Prolongation(int send, int rece, int tolayer, int where)
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
    				mgc.M_value[rece] = (mgc.p_Pro_xm)->M_value[send];
					break;

				case 1: 
					mgc.M_value[rece] = ((mgc.p_Pro_xm)->M_value[send]
										+(mgc.p_Pro_xp)->M_value[send])*0.5;
					break;
					
				case 2: 
					mgc.M_value[rece] = ((mgc.p_Pro_ym)->M_value[send]
										+(mgc.p_Pro_yp)->M_value[send])*0.5;

					break;

				case 3: 
					mgc.M_value[rece] =			  ((mgc.p_Pro_xm)->M_value[send]
					+(mgc.p_Pro_xp)->M_value[send]+(mgc.p_Pro_ym)->M_value[send]
					+(mgc.p_Pro_yp)->M_value[send])*0.25;
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
    				mgc.M_value[rece] = (mgc.p_Pro_xm)->M_value[send];
					break;

				case 1: 
					mgc.M_value[rece] = ((mgc.p_Pro_xm)->M_value[send]
										+(mgc.p_Pro_xp)->M_value[send])*0.5;
					break;
					
				case 2: 
					mgc.M_value[rece] = ((mgc.p_Pro_ym)->M_value[send]
										+(mgc.p_Pro_yp)->M_value[send])*0.5;

					break;

				case 3: 
					mgc.M_value[rece] =			  ((mgc.p_Pro_xm)->M_value[send]
					+(mgc.p_Pro_xp)->M_value[send]+(mgc.p_Pro_ym)->M_value[send]
					+(mgc.p_Pro_yp)->M_value[send])*0.25;
					break;

			}

		}

	}


	break;


}



	return;
}



void MultiGrid::SendtoBottom(int what)
{
	int i, j, n;

	int nsend;
	nsend = LayerGridX[MPI_Layer]*LayerGridY[MPI_Layer];
	double mysend[nsend];
	double myreceive[BottomCells];

	n = 0;
	for(j=1; j<=LayerGridY[MPI_Layer]; j++)
	{
		for(i=1; i<=LayerGridX[MPI_Layer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,MPI_Layer);
			mysend[n] = mgc.M_value[what];
			n++;

		}
	}

	MPI_Gatherv(&mysend, nsend, MPI_DOUBLE, &myreceive, BottomSend, Bottomdisp,
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
					mgc.M_value[what] = myreceive[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)];
				}

			}

		}

	}


	return;

}



void MultiGrid::BottomSendBack(int what)
{
	int i, j, n;

	int nrece;
	nrece = LayerGridX[MPI_Layer]*LayerGridY[MPI_Layer];

	double myrece[nrece];
	double mysend[BottomCells];


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
					mysend[ Bottomdisp[n]+(j-1)*recebottomX[n]+(i-1)] = mgc.M_value[what];
				
				}

			}

		}

	}

	MPI_Scatterv(&mysend, BottomSend, Bottomdisp, MPI_DOUBLE, &myrece, nrece,
				MPI_DOUBLE, Worker, MPI_COMM_WORLD);


	n = 0;
	for(j=1; j<=LayerGridY[MPI_Layer]; j++)
	{
		for(i=1; i<=LayerGridX[MPI_Layer]; i++)
		{

			MG_Cell &mgc = GetMGCell(i,j,MPI_Layer);
			mgc.M_value[what] = myrece[n];
			n++;

		}
	}



	return;

}




void MultiGrid::MG_BottomLayer(int field)
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
			Relaxation(field,n,1);
		}

		//find the residual;
		Residual(field,n,1);

		//restrict to next layer;
		Restriction(MG_Res,MG_Sou,n+1,1);

		SetZero(MG_Phi, n+1,1);

	}

//============================================================
//==============  V-Cycle Layers Below MPI Layer  ============
//============================================================

	MG_Cell &c1 = GetMGBCell(1, 1, SER_Layer);
	MG_Cell &c2 = GetMGBCell(1, 2, SER_Layer);
	MG_Cell &c3 = GetMGBCell(2, 1, SER_Layer);
	MG_Cell &c4 = GetMGBCell(2, 2, SER_Layer);

	c1.M_value[0] = -(7*c1.M_value[1]+2*c2.M_value[1]+2*c3.M_value[1]+c4.M_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c2.M_value[0] = -(7*c2.M_value[1]+2*c1.M_value[1]+2*c4.M_value[1]+c3.M_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c3.M_value[0] = -(7*c3.M_value[1]+2*c4.M_value[1]+2*c1.M_value[1]+c2.M_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;
	c4.M_value[0] = -(7*c4.M_value[1]+2*c3.M_value[1]+2*c2.M_value[1]+c1.M_value[1])
					*MeshAmplif[MPI_Layer+SER_Layer-1]/24;

//============================================================
//==============     V-Cycle Moving Up      ==================
//============================================================


	for (n=SER_Layer; n>1; n--)
	{

		Prolongation(MG_Phi, MG_Res, n-1,1);

		AddCorrection(n-1,1);

		for(int m=0; m< u2; m++)
		{
			Relaxation(field,n-1,1);
		}

	}


	break;



	case 2:
	
	break;


	}


	return;
}


void MultiGrid::Relaxation(int field, int layer, int where)
{

	int i,j;
	int nx,ny, amp;


	switch(RelaxType)
	{

	// case 0: Hybird Gauss-Siedel Relaxation
	// case 1: Red-Black Relaxation
	// case 2: 
		case 0:
		default:
		nx=LayerGridX[layer];
		ny=LayerGridY[layer];
		amp=MeshAmplif[layer];

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

					ccc.M_value[0]=(1-omega)*ccc.M_value[0]+omega*(cxm.M_value[0]+cxp.M_value[0]
					+cym.M_value[0]+cyp.M_value[0]-ccc.M_value[1]*amp)/(4.0+ccc.M_value[4]*amp);
				}

			}
			

		// //Bottom layer
		// 	for (j=1; j<=BLayerGrid[layer]; j++)
		// 	{
		// 		for (i=1; i<=BLayerGrid[layer]; i++)
		// 		{
		// 			MG_Cell &ccc = GetMGBCell(i,   j, layer);
		// 			MG_Cell &cxm = GetMGBCell(i-1, j, layer);
		// 			MG_Cell &cxp = GetMGBCell(i+1, j, layer);
		// 			MG_Cell &cym = GetMGBCell(i, j-1, layer);
		// 			MG_Cell &cyp = GetMGBCell(i, j+1, layer);

		// 			ccc.M_value[0]=(1-omega)*ccc.M_value[0]+omega*(cxm.M_value[0]+cxp.M_value[0]
		// 			+cym.M_value[0]+cyp.M_value[0]-ccc.M_value[1]*MeshAmplif[MPI_Layer-1+layer])
		// 			/(4.0+ccc.M_value[4]*MeshAmplif[MPI_Layer-1+layer]);
		// 		}
		// 	}

		// break;
		// }

		break;

	}

	return;
}




void MultiGrid::Residual(int field, int layer, int where)
{

	int i,j;
	int nx,ny, amp;


// switch(where)
// {

// 	case 0:
	nx=LayerGridX[layer];
	ny=LayerGridY[layer];
	amp=MeshAmplif[layer];

	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{	
			MG_Cell &ccc = GetMGCell(i,   j, layer);
			MG_Cell &cxm = GetMGCell(i-1, j, layer);
			MG_Cell &cxp = GetMGCell(i+1, j, layer);
			MG_Cell &cym = GetMGCell(i, j-1, layer);
			MG_Cell &cyp = GetMGCell(i, j+1, layer);

			ccc.M_value[2]=ccc.M_value[1]-(cxm.field+cxp.field
			+cym.field+cyp.field-(4+ccc.M_value[4]*amp)*ccc.field)/amp;
		}

	}
	// break;


// 	case 1:

// 	for (j=1; j<=BLayerGrid[layer]; j++)
// 	{
// 		for (i=1; i<=BLayerGrid[layer]; i++)
// 		{
// 			MG_Cell &ccc = GetMGBCell(i,   j, layer);
// 			MG_Cell &cxm = GetMGBCell(i-1, j, layer);
// 			MG_Cell &cxp = GetMGBCell(i+1, j, layer);
// 			MG_Cell &cym = GetMGBCell(i, j-1, layer);
// 			MG_Cell &cyp = GetMGBCell(i, j+1, layer);

// 			ccc.M_value[2]=ccc.M_value[1]-(cxm.M_value[0]+cxp.M_value[0]
// 			+cym.M_value[0]+cyp.M_value[0]-(4+ccc.M_value[4]*MeshAmplif[MPI_Layer-1+layer])*ccc.M_value[0])/MeshAmplif[MPI_Layer-1+layer];
// 		}

// 	}
// 	break;

// }

	return;
}





void MultiGrid::SetZero(int what, int layer, int where)
{
	int i,j;
	int nx,ny;


// switch(where)
// {

// 	case 0:

	nx=LayerGridX[layer];
	ny=LayerGridY[layer];

	// switch(where)
	// {

	// case 0:
	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{

			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.M_value[what] = 0;

		}

	}	

// 	break;

// 	case 1:

// 	for (j=0; j<=BLayerGrid[layer]+1; j++)
// 	{
// 		for (i=0; i<=BLayerGrid[layer]+1; i++)
// 		{

// 			MG_Cell &ccc = GetMGBCell(i, j, layer);
// 			ccc.M_value[what] = 0;

// 		}

// 	}

// 	break;

// }




	return;
}



void MultiGrid::Exchange(int what, int layer)
{

	switch(what)
	{
		case 0:
			p_domain()->p_Com()->DoCommute(COMMU_MG_P, layer);
			break;
		case 1:
			p_domain()->p_Com()->DoCommute(COMMU_MG_S, layer);
			break;
		case 2:
			p_domain()->p_Com()->DoCommute(COMMU_MG_R, layer);
			break;

	}


	//adjust corner cell;
	if(what == 2)
	{

		MG_Cell &cmm = GetMGCell(0, 0, layer);
		cmm.M_value[what] = (GetMGCell(0, 1, layer).M_value[what] 
							+GetMGCell(1, 0, layer).M_value[what])*0.5;

		MG_Cell &cpm = GetMGCell(LayerGridX[layer]+1, 0, layer);
		cpm.M_value[what] = (GetMGCell(LayerGridX[layer]+1, 1, layer).M_value[what] 
							+GetMGCell(LayerGridX[layer],   0, layer).M_value[what])*0.5;

		MG_Cell &cmp = GetMGCell(0, LayerGridY[layer]+1, layer);
		cmp.M_value[what] = (GetMGCell(0, LayerGridY[layer],   layer).M_value[what] 
							+GetMGCell(1, LayerGridY[layer]+1, layer).M_value[what])*0.5;

		MG_Cell &cpp = GetMGCell(LayerGridX[layer]+1, LayerGridY[layer]+1, layer);
		cpp.M_value[what] = (GetMGCell(LayerGridX[layer],   LayerGridY[layer]+1, layer).M_value[what] 
							+GetMGCell(LayerGridX[layer]+1, LayerGridY[layer],   layer).M_value[what])*0.5;
	}

	

	return;
}



void MultiGrid::AddCorrection(int layer, int where)
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
			ccc.M_value[0] +=ccc.M_value[2];
		
		}

	}

	// break;

	// case 1:

	// for (j=1; j<=BLayerGrid[layer]; j++)
	// {
	// 	for (i=1; i<=BLayerGrid[layer]; i++)
	// 	{
	// 		MG_Cell &ccc = GetMGBCell(i, j, layer);
	// 		ccc.M_value[0] +=ccc.M_value[2];
		
	// 	}

	// }

	// break;

	// }



	return;
}


double MultiGrid::FindError(double &maxall)
{

	int i,j;
	int nx,ny;


// switch(where)
// {

// 	case 0:

	nx=LayerGridX[1];
	ny=LayerGridY[1];

	double epsn;
	double epsp;
	double eps;

	double maxp;

	epsp = 0.0;
	maxp=0;


	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			MG_Cell &ccc = GetMGCell(i, j, 1);
			epsn = (ccc.M_value[0] -ccc.M_value[3]);
			epsp += epsn*epsn;

			if(abs(ccc.M_value[0])>maxp) maxp=abs(ccc.M_value[0]);
		
		}

	}

	MPI_Allreduce(&epsp, &eps,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&maxp, &maxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return eps;
}

void MultiGrid::Put_Source(int field, double k0, int k)
{

	int i,j;
	int nx,ny;

	nx=LayerGridX[1];
	ny=LayerGridY[1];

	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{

			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);

			if(k>1)
			{
				Cell &cm  = p_Meshs->GetCell(i, j, k-1);
				Cell &cmm = p_Meshs->GetCell(i, j, k-2);
				mgc.M_value[0] = 2*cm.W_Fields[field]-cmm.W_Fields[field];	//initial guess
			}
			else
			{
				mgc.M_value[0] = 0.0;				//initial guess
			}
		}
	}

	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{

			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);
			switch(field)
			{
				//===Psi==================================
				case 0:
				mgc.M_value[1] = (c.W_Denn-c.W_Deni + c.B_Den - c.B_Jz)*dxdy;
				mgc.M_value[4] = 0.0;
				break;
				//===Ez==================================
				case 1:
				mgc.M_value[1] = -p_Meshs->Dive_J(i, j, k0, k)*dxdy;
				mgc.M_value[4] = 0.0;
				break;
				//===Bz==================================
				case 2:
				mgc.M_value[1] =  p_Meshs->Curl_J(i, j, k0, k)*dxdy;
				mgc.M_value[4] = 0.0;
				break;
				//===Bx==================================
				case 3:
				mgc.M_value[1] =  p_Meshs->SourceY(i, j, k0, k)*dxdy;
				mgc.M_value[4] =  c.W_Chi*dxdy;
				break;
				//===By==================================
				case 4:
				mgc.M_value[1] = -p_Meshs->SourceX(i, j, k0, k)*dxdy;
				mgc.M_value[4] =  c.W_Chi*dxdy;
				break;

			}
		}

	}

	return;
}



void MultiGrid::Put_Fields(int field, int k)
{

	int i,j;
	int nx,ny;

	nx=LayerGridX[1];
	ny=LayerGridY[1];


	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{
			Cell 	  &c = p_Meshs->GetCell(i, j, k);
			MG_Cell &mgc = GetMGCell(i, j, 1);
			c.W_Fields[field] = mgc.M_value[0];   
				
		}

	}
	return;
}



int MultiGrid::MG_V_cycle(int field, double k0, int k)
{


	int i,j,n;
	double eps;
	double maxall=1.0;

//============================================================
//==============   Put Source For Different Equation =========
//============================================================
	Put_Source(field, k0, k);


//============================================================
//==============    Restrict Additional Coefficient  =========
//==============   	 for Helmholtz Equation 		 =========
//============================================================

	for (n=1; n<MPI_Layer; n++)
	{
		RestrictionB(MG_Chi,MG_Chi,n+1, 0);
	}

	switch(BottomType)
	{
		case 1:
		SendtoBottom(MG_Chi);
		for (n=1; n<SER_Layer; n++)
		{
			RestrictionB(MG_Chi,MG_Chi,n+1,1);
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
	//if(Rank==0) std::cout<<k<<'\n';

	iter++;
	//record old value
	for (j=1; j<=LayerGridY[1]; j++)
	{
		for (i=1; i<=LayerGridX[1]; i++)
		{
			GetMGCell(i, j, 1).M_value[3] = GetMGCell(i, j, 1).M_value[0];

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
			Relaxation(field, n, 0);
			Exchange(MG_Phi, n);
		}
		
		//find the residual;
   		Residual(field, n, 0);
		Exchange(MG_Res, n);
		
		Restriction(MG_Res,MG_Sou,n+1, 0);
		
		SetZero(MG_Phi, n+1, 0);



	}
	
//============================================================
//==============  V-Cycle Layers Below MPI Layer  ============
//============================================================


	switch(BottomType)
	{

		case 0:

			for(int m=0; m< u1; m++)
			{
				Relaxation(field, n, 0);
				Exchange(MG_Phi, MPI_Layer);
			}

		break;

		case 1:

			SendtoBottom(MG_Sou);
			if(Rank == Worker) 
			{
				MG_BottomLayer(field);
			}
			BottomSendBack(MG_Phi);

		break;

	}

//============================================================
//==============     V-Cycle Moving Up      ==================
//============================================================


	for (n=MPI_Layer; n>1; n--)
	{


		Prolongation(MG_Phi, MG_Res, n-1, 0);
		AddCorrection(n-1, 0);


		for(int m=0; m<u2; m++)
		{
			Exchange(MG_Phi, n-1);
			Relaxation(field,n-1, 0);
		}
		// sprintf(name,"After_%d_Relax2_%d.dg",n,Rank);
		// DebugWrite(n-1, 0, name, 0);
	}
//============================================================
//==============     Error Estimation      ==================
//============================================================
	eps = FindError(maxall)/GridXY;

	//======== test==========
	//if(Rank == 0) std::cout<<"iter: "<<iter<<";   eps:"<<eps<<'\n';
	//======== test==========

//============================================================
//==============     Fail to Converge.      ==================
//============================================================
	if(iter>1e3 || eps>1e10)
	{
		if (Rank==0)  std::cout <<"==== Multigrid: Equation Type: "<<field<< " Failed To Converge at k(z) ="<<k<<'\n';
		return 1;
	}

}


//============================================================
//==============   Put Solution Back to Cell         =========
//============================================================
	Put_Fields(field, k);

	return 0;



}

void MultiGrid::DebugWrite(int layer, int what, char *name, int where)
{

	  FILE * dFile;
      int i,j;
	if(where==0)
	{
		dFile = fopen (name,"w");
		for (j=1; j<=LayerGridY[layer]; j++)
		{
			for (i=1; i<=LayerGridX[layer]; i++)
			{
				MG_Cell &c = GetMGCell(i,j,layer);
				fprintf(dFile, "%12.10f ", c.M_value[what]);
			}
			fprintf(dFile, "\n");
      	}
	fclose (dFile);
	}
	else
	{
		dFile = fopen (name,"w");
		for (j=1; j<=BLayerGrid[layer]; j++)
		{
			for (i=1; i<=BLayerGrid[layer]; i++)
			{
				MG_Cell &c = GetMGBCell(i,j,layer);
				fprintf(dFile, "%12.10f ", c.M_value[what]);
			}
			fprintf(dFile, "\n");
      	}

	fclose (dFile);
	}
	return;
}

//=========================================================

MultiGrid::~MultiGrid()
{
	delete[] p_MGCell;
	delete[] p_MGBCell;
	delete[] Bottomdisp;
	delete[] BottomSend;
	delete[] recebottomX;
	delete[] recebottomY;
};







