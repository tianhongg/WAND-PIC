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
	
	BottomType=0; //May-25-2021; Tianhong.

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


	//May-25-2021- Tianhong updates; 
	std::vector< std::vector< std::pair<WDOUBLE,int> > > Dx;
	std::vector< std::pair<WDOUBLE,int> > tmp;

	for(i=0;i<p_Meshs->GridsTmp.size();i++)
	{
		tmp.push_back({p_Meshs->GridsTmp[i],i}); // cell size and grand index
	}
	Dx.push_back(tmp);

	int gridmin=2;

	MPI_Layer = 1;
	n=1;
	while(gridmin>1)
	{
		tmp.clear();
		tmp.push_back({0.0,0}); //left ghost

		MPI_Layer++;
		n++;
		for(int i=1;i<Dx[n-2].size()-1;i+=2)
		{
			WDOUBLE nextsize=Dx[n-2][i].first+(Dx[n-2][i-1].first+Dx[n-2][i+1].first)*0.5;

			tmp.push_back( {nextsize, Dx[n-2][i].second} );


			if(Dx[n-2][i].second >= ((RankIdx_X-1)*XGridN+1) && Dx[n-2][i].second<= (RankIdx_X)*XGridN )
    		{
				LayerGridX[MPI_Layer] += 1;
        	}

        	if(Dx[n-2][i].second  >= ((RankIdx_Y-1)*YGridN+1) && Dx[n-2][i].second <= (RankIdx_Y)*YGridN )
    		{
				LayerGridY[MPI_Layer] += 1;
        	}

		}

		tmp[0].first=tmp[1].first;
		tmp.push_back( { tmp[tmp.size()-1].first, p_Meshs->GridsTmp.size()-1} ); //right ghost
		Dx.push_back(tmp);

		int minGrid=std::min(LayerGridY[MPI_Layer],LayerGridX[MPI_Layer]);
		MPI_Allreduce(&minGrid, &gridmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	}


	// //=====
	// if(Rank==0)
	
	// {FILE * dFile;
	// char name[128];
 //  	sprintf(name,"AllLayers_%d_.dg",Rank);
	// dFile = fopen (name,"w");

	// for (n=1; n<=MPI_Layer; n++)
	// {
	// 	for (j=0; j<Dx[n-1].size(); j++)
	// 	{
		
	// 		fprintf(dFile, "[%d] ", Dx[n-1][j].second);
	// 	}
	// 	fprintf(dFile, "\n");
	// }
	// //====
	// fclose (dFile);
	// }


	int Tot_GridX = 0;
	int Tot_GridY = 0;


	Layer_Sum[1] =0;
	for (i=1; i<=MPI_Layer; i++) //layer starts with 1
	{
    	Tot_GridX += LayerGridX[i]+2;
    	Tot_GridY += LayerGridX[i]+2;
    	
    	if(i>1)
    	{
    		// Find the 2D Cells Sum of all Layers above this Layer
    		Layer_Sum[i]   = Layer_Sum[i-1]+(LayerGridX[i-1]+2)*(LayerGridY[i-1]+2);
    	}
	}


// ----------------------------------------------------------------
	// Find the Cell Grand_Index in X(or Y) direction at each layer
// ----------------------------------------------------------------
	
	std::vector< std::vector<int> >    Cellidx_i(MPI_Layer,std::vector<int>() );
	std::vector< std::vector<int> >    Cellidx_j(MPI_Layer,std::vector<int>() );

	std::vector< std::vector<WDOUBLE> > Celldx(MPI_Layer,std::vector<WDOUBLE>() );
	std::vector< std::vector<WDOUBLE> > Celldy(MPI_Layer,std::vector<WDOUBLE>() );


	for (n=1; n<=MPI_Layer; n++)
	{
    	int iflag=0;
    	int jflag=0;
		for(i=0;i<Dx[n-1].size();i++)
		{
			WDOUBLE ddx = Dx[n-1][i].first;
			int    idx = Dx[n-1][i].second;

			//x-dir
			if(idx >= ((RankIdx_X-1)*XGridN+1) && idx<= (RankIdx_X)*XGridN )
    		{
				if(iflag==0)
				{	
					Celldx[n-1].push_back(Dx[n-1][i-1].first);
					Cellidx_i[n-1].push_back(Dx[n-1][i-1].second);
					iflag=1;
				}

				Celldx[n-1].push_back(ddx);
				Cellidx_i[n-1].push_back(idx);
        	}
        	else if(iflag)
        	{	iflag=0;
        		Celldx[n-1].push_back(ddx);
				Cellidx_i[n-1].push_back(idx);
        	}


    		//y-dir
        	if(idx >= ((RankIdx_Y-1)*YGridN+1) && idx<= (RankIdx_Y)*YGridN )
    		{
				if(jflag==0)
				{	
					Celldy[n-1].push_back(Dx[n-1][i-1].first);
					Cellidx_j[n-1].push_back(Dx[n-1][i-1].second);
					jflag=1;
				}

				Celldy[n-1].push_back(ddx);
				Cellidx_j[n-1].push_back(idx);
				
				
        	}
        	else if(jflag)
			{
				jflag=0;
				Celldy[n-1].push_back(ddx);
				Cellidx_j[n-1].push_back(idx);
			}


		}


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
				mgc.Grandidx_X = Cellidx_i[n-1][i];
				mgc.Grandidx_Y = Cellidx_j[n-1][j];

				mgc.dx = Celldx[n-1][i];
				mgc.dy = Celldy[n-1][j];

			}
		}
	}

	// //=====
	// FILE * dFile;
	// char name[128];
 //  	sprintf(name,"MultiGrids_%d_.dg",Rank);
	// dFile = fopen (name,"w");

	// for (n=1; n<=MPI_Layer; n++)
	// {
	// 	for (j=0; j<=LayerGridY[n]+1; j++)
	// 	{

	// 		MG_Cell &mgc = GetMGCell(0,j,n);
	// 		fprintf(dFile, "[%f] ", mgc.dy);
	// 	}
	// 	fprintf(dFile, "\n");
	// }
	// //====
	// fclose (dFile);


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

				while (Cellidx_i[n-2][ii]<Cellidx_i[n-1][i])
				{ ii++; }
				while (Cellidx_j[n-2][jj]<Cellidx_j[n-1][j])
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
					if(Cellidx_i[n][ii]==Cellidx_i[n-1][i])
					{
						findx = 1;
						break;
					}

					if(Cellidx_i[n][ii]<Cellidx_i[n-1][i]&&Cellidx_i[n][ii+1]>Cellidx_i[n-1][i])
					{
						break;
					}
				} 


				for (jj = 0; jj<=LayerGridY[n+1]+1; jj++)
				{
					if(Cellidx_j[n][jj]==Cellidx_j[n-1][j])
					{
						findy = 1;
						break;
					}

					if(Cellidx_j[n][jj]<Cellidx_j[n-1][j]&&Cellidx_j[n][jj+1]>Cellidx_j[n-1][j])
					{
						break;
					}

				}

				//type 0 prolongation
				if(findx==1 &&findy == 1)
				{
					mgc.Protype = 0;
					mgc.p_Pro_xm = &GetMGCell(ii,jj,n+1);

				}

				//type 1 prolongation
				if(findx==0 &&findy == 1)
				{
					mgc.Protype = 1;
					mgc.p_Pro_xm = &GetMGCell(ii,  jj,n+1);
					mgc.p_Pro_xp = &GetMGCell(ii+1,jj,n+1);
				}
				//type 2 prolongation
				if(findx==1 &&findy == 0)
				{
					mgc.Protype = 2;
					mgc.p_Pro_ym = &GetMGCell(ii,  jj,n+1);
					mgc.p_Pro_yp = &GetMGCell(ii,jj+1,n+1);
				}
				//type 3 prolongation
				if(findx==0 &&findy == 0)
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



/*  removed by Tianhong; May-25
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

*/
	//==== end worker rank=====//

	for (i=1; i<=MPI_Layer+SER_Layer; i++)
	{
    	MeshAmplif[i] = pow(2,i-1)*pow(2,i-1);
	}

	if (Rank==0)  std::cout << "==== Multigrid: Multigrid Meshes Done.   ====\n";


};
	




// Full Weight Restriction
void MultiGrid::Restriction(int send, int rece, int tolayer, int where) //v
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

			mgc.M_value[rece] = (
				mgc.p_Res_mm->M_value[send]*wmm + mgc.p_Res_xm->M_value[send]*wxm + mgc.p_Res_mp->M_value[send]*wmp
			  + mgc.p_Res_ym->M_value[send]*wym + mgc.p_Res_cc->M_value[send]*wcc + mgc.p_Res_yp->M_value[send]*wyp	
			  + mgc.p_Res_pm->M_value[send]*wpm + mgc.p_Res_xp->M_value[send]*wxp + mgc.p_Res_pp->M_value[send]*wpp
			)/wa;

		}

	}
	break;

	case 1:
	/*
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
	*/
	break;

}

	return;
}



void MultiGrid::RestrictionB(int send, int rece, int tolayer, int where)//?
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


void MultiGrid::Prolongation(int send, int rece, int tolayer, int where) //v
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
    				mgc.M_value[rece] = (mgc.p_Pro_xm)->M_value[send];
					break;

				case 1: 
					dxm=mgc.p_Pro_xm->dx;
					dxp=mgc.p_Pro_xp->dx;
					mgc.M_value[rece] = ((mgc.p_Pro_xm)->M_value[send]*dxp
										+(mgc.p_Pro_xp)->M_value[send]*dxm)/(dxm+dxp);
					break;
					
				case 2: 
					dym=mgc.p_Pro_ym->dy;
					dyp=mgc.p_Pro_yp->dy;
					mgc.M_value[rece] = ((mgc.p_Pro_ym)->M_value[send]*dyp
										+(mgc.p_Pro_yp)->M_value[send]*dym)/(dym+dyp);
					break;

				case 3: 
					dxm=mgc.p_Pro_xm->dx;
					dxp=mgc.p_Pro_xp->dx;
					dym=mgc.p_Pro_ym->dy;
					dyp=mgc.p_Pro_yp->dy;

					mgc.M_value[rece] = ((mgc.p_Pro_xm)->M_value[send]*dxp*dyp
					+(mgc.p_Pro_xp)->M_value[send]*dxm*dyp+(mgc.p_Pro_ym)->M_value[send]*dxp*dym
					+(mgc.p_Pro_yp)->M_value[send]*dxm*dym)/(dxm+dxp)/(dym+dyp);
					break;

			}

		}

	}

	break;


	case 1:
	/*
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
	*/


	break;


}



	return;
}



//deprecated//May-25-tianhong
void MultiGrid::SendtoBottom(int what)
{
	


	return;

}


//deprecated//May-25-tianhong
void MultiGrid::BottomSendBack(int what)
{
	


	return;

}



//deprecated//May-25-tianhong
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


void MultiGrid::Relaxation(int field, int layer, int where) //v 
{

	int i,j;
	int nx,ny, amp;

	WDOUBLE hxp,hxm,hxa,hxd,h2x;
	WDOUBLE hyp,hym,hya,hyd,h2y;

	WDOUBLE wcc,wxm,wxp,wym,wyp,wmm,wmp,wpm,wpp;
	WDOUBLE d2xS, d2yS, dxS, dyS;
	WDOUBLE kcc,kxm,kxp,kym,kyp;

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

			kcc=ccc.M_value[4];
			kxm=cxm.M_value[4];
			kxp=cxp.M_value[4];
			kym=cym.M_value[4];
			kyp=cyp.M_value[4];

			wmm = (h2x + h2y - 2*hxd*hxp - 2*hyd*hyp)/(3*hxa*hxm*hya*hym);
			wmp = (h2x + h2y - 2*hxd*hxp + 2*hyd*hym)/(3*hxa*hxm*hya*hyp);
			wpm = (h2x + h2y + 2*hxd*hxm - 2*hyd*hyp)/(3*hxa*hxp*hya*hym);
			wpp = (h2x + h2y + 2*hxd*hxm + 2*hyd*hym)/(3*hxa*hxp*hya*hyp);

			wxm = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x-2*hxd*hxp)*(2+hym*hyp*kxm) )/(6*hxa*hxm*hym*hyp);
			wxp = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x+2*hxd*hxm)*(2+hym*hyp*kxp) )/(6*hxa*hxp*hym*hyp);
			wym = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y-2*hyd*hyp)*(2+hxm*hxp*kym) )/(6*hxm*hxp*hya*hym);
			wyp = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y+2*hyd*hym)*(2+hxm*hxp*kyp) )/(6*hxm*hxp*hya*hyp);

			wcc = ( h2y*(-2+hxm*hxp*kcc) + h2x*(-2+hym*hyp*kcc) - 2*(hym*hyp*(4+hxd*hxd*kcc) + hxm*hxp*(4+hyd*hyd*kcc)) )/(6*hxm*hxp*hym*hyp);

			dxS = (cxp.M_value[1]-ccc.M_value[1])*hxm/hxp/hxa + (ccc.M_value[1]-cxm.M_value[1])*hxp/hxm/hxa;
			dyS = (cyp.M_value[1]-ccc.M_value[1])*hym/hyp/hya + (ccc.M_value[1]-cym.M_value[1])*hyp/hym/hya;

			d2xS = (hxm*cxp.M_value[1]-hxa*ccc.M_value[1]+hxp*cxm.M_value[1])*2/hxa/hxp/hxm;
			d2yS = (hym*cyp.M_value[1]-hya*ccc.M_value[1]+hyp*cym.M_value[1])*2/hya/hyp/hym;

			ccc.M_value[0]=(1-omega)*ccc.M_value[0]+omega*
			( ccc.M_value[1] + h2x/12*d2xS + h2y/12*d2yS + hxd/3*dxS + hyd/3*dyS
			 - wxm*cxm.M_value[0]- wxp*cxp.M_value[0]- wym*cym.M_value[0]- wyp*cyp.M_value[0]
			 - wmm*cmm.M_value[0]- wpm*cpm.M_value[0]- wmp*cmp.M_value[0]- wpp*cpp.M_value[0])/(wcc-kcc) ;
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




void MultiGrid::Residual(int field, int layer, int where) //v
{

	int i,j;
	int nx,ny, amp;


	WDOUBLE hxp,hxm,hxa,hxd,h2x;
	WDOUBLE hyp,hym,hya,hyd,h2y;

	WDOUBLE wcc,wxm,wxp,wym,wyp,wmm,wmp,wpm,wpp;
	WDOUBLE d2xS, d2yS, dxS, dyS;
	WDOUBLE kcc,kxm,kxp,kym,kyp;

// switch(where)
// {

// 	case 0:
	nx=LayerGridX[layer];
	ny=LayerGridY[layer];
	// amp=MeshAmplif[layer];

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

			kcc=ccc.M_value[4];
			kxm=cxm.M_value[4];
			kxp=cxp.M_value[4];
			kym=cym.M_value[4];
			kyp=cyp.M_value[4];

			wmm = (h2x + h2y - 2*hxd*hxp - 2*hyd*hyp)/(3*hxa*hxm*hya*hym);
			wmp = (h2x + h2y - 2*hxd*hxp + 2*hyd*hym)/(3*hxa*hxm*hya*hyp);
			wpm = (h2x + h2y + 2*hxd*hxm - 2*hyd*hyp)/(3*hxa*hxp*hya*hym);
			wpp = (h2x + h2y + 2*hxd*hxm + 2*hyd*hym)/(3*hxa*hxp*hya*hyp);

			wxm = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x-2*hxd*hxp)*(2+hym*hyp*kxm) )/(6*hxa*hxm*hym*hyp);
			wxp = -( 2*h2y-4*(hyd*hyd+3*hym*hyp) + (h2x+2*hxd*hxm)*(2+hym*hyp*kxp) )/(6*hxa*hxp*hym*hyp);
			wym = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y-2*hyd*hyp)*(2+hxm*hxp*kym) )/(6*hxm*hxp*hya*hym);
			wyp = -( 2*h2x-4*(hxd*hxd+3*hxm*hxp) + (h2y+2*hyd*hym)*(2+hxm*hxp*kyp) )/(6*hxm*hxp*hya*hyp);

			wcc = ( h2y*(-2+hxm*hxp*kcc) + h2x*(-2+hym*hyp*kcc) - 2*(hym*hyp*(4+hxd*hxd*kcc) + hxm*hxp*(4+hyd*hyd*kcc)) )/(6*hxm*hxp*hym*hyp);

			dxS = (cxp.M_value[1]-ccc.M_value[1])*hxm/hxp/hxa + (ccc.M_value[1]-cxm.M_value[1])*hxp/hxm/hxa;
			dyS = (cyp.M_value[1]-ccc.M_value[1])*hym/hyp/hya + (ccc.M_value[1]-cym.M_value[1])*hyp/hym/hya;

			d2xS = (hxm*cxp.M_value[1]-hxa*ccc.M_value[1]+hxp*cxm.M_value[1])*2/hxa/hxp/hxm;
			d2yS = (hym*cyp.M_value[1]-hya*ccc.M_value[1]+hyp*cym.M_value[1])*2/hya/hyp/hym;

			ccc.M_value[2]=ccc.M_value[1]-
			(  wxm*cxm.M_value[0]+ wxp*cxp.M_value[0]+ wym*cym.M_value[0]+ wyp*cyp.M_value[0]
			 + wmm*cmm.M_value[0]+ wpm*cpm.M_value[0]+ wmp*cmp.M_value[0]+ wpp*cpp.M_value[0]
			 + (wcc-kcc)*ccc.M_value[0] -h2x/12*d2xS -h2y/12*d2yS - hxd/3*dxS - hyd/3*dyS
			);
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



	nx=LayerGridX[layer];
	ny=LayerGridY[layer];


	for (j=0; j<=ny+1; j++)
	{
		for (i=0; i<=nx+1; i++)
		{

			MG_Cell &ccc = GetMGCell(i, j, layer);
			ccc.M_value[what] = 0.0;

		}

	}	



	return;
}

void MultiGrid::Exchange(int what, int layer)
{

	switch(what)
	{
		case MG_Phi:
			p_domain()->p_Com()->DoCommute(COMMU_MG_P, layer);
			break;
		case MG_Sou:
			p_domain()->p_Com()->DoCommute(COMMU_MG_S, layer);
			break;
		case MG_Res:
			p_domain()->p_Com()->DoCommute(COMMU_MG_R, layer);
		case MG_Chi:
			p_domain()->p_Com()->DoCommute(COMMU_MG_C, layer);
			break;
	}
	return;
}



void MultiGrid::AddCorrection(int layer, int where)//v
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
	

	return;
}


WDOUBLE MultiGrid::FindError(WDOUBLE &maxall) //v
{

	int i,j;
	int nx,ny;


// switch(where)
// {

// 	case 0:
	nx=LayerGridX[1];
	ny=LayerGridY[1];

	WDOUBLE epsn;
	WDOUBLE epsp;
	WDOUBLE eps;

	WDOUBLE maxp;

	epsp = 0.0;
	maxp=0;

	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			MG_Cell &ccc = GetMGCell(i, j, 1);
			epsn = (ccc.M_value[0] -ccc.M_value[3]);
			epsp += epsn*epsn;
			maxp=std::max(maxp,abs(ccc.M_value[0]));
		
		}

	}

	MPI_Allreduce(&epsp, &eps,    1, MPI_WDOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&maxp, &maxall, 1, MPI_WDOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return eps;
}

void MultiGrid::Put_Source(int field, WDOUBLE k0, int k) //v
{

	int i,j;
	int nx,ny;

	// int time =p_domain()->Get_Step();


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
				// else mgc.M_value[0] = c.W_Fields[field];
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
				mgc.M_value[1] = (c.W_Denn-c.W_Deni + c.B_Den - c.B_Jz);
				mgc.M_value[4] = 0.0;
				break;
				//===Ez==================================
				case 1:
				mgc.M_value[1] = -p_Meshs->Dive_J(i, j, k0, k);
				mgc.M_value[4] = 0.0;
				break;
				//===Bz==================================
				case 2:
				mgc.M_value[1] =  p_Meshs->Curl_J(i, j, k0, k);
				mgc.M_value[4] = 0.0;
				break;
				//===Bx==================================
				case 3:
				mgc.M_value[1] =  p_Meshs->SourceY(i, j, k0, k);
				mgc.M_value[4] =  c.W_Chi;
				break;
				//===By==================================
				case 4:
				mgc.M_value[1] = -p_Meshs->SourceX(i, j, k0, k);
				mgc.M_value[4] =  c.W_Chi;
				break;

			}
		}

	}

	return;
}



void MultiGrid::Put_Fields(int field, int k)//v
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



int MultiGrid::MG_V_cycle(int field, WDOUBLE k0, int k)
{


	int i,j,n;
	WDOUBLE eps;
	WDOUBLE maxall=1.0;

//============================================================
//==============   Put Source For Different Equation =========
//============================================================
	Put_Source(field, k0, k);
	Exchange(MG_Sou, 1);
	if(field>2) Exchange(MG_Chi,1);



//============================================================
//==============    Restrict Additional Coefficient  =========
//==============   	 for Helmholtz Equation 		 =========
//============================================================

	for (n=1; n<MPI_Layer; n++)
	{
		RestrictionB(MG_Chi,MG_Chi,n+1, 0); 
		Exchange(MG_Chi, n+1);
	}
	
	// switch(BottomType)
	// {
	// 	case 1:
	// 	// SendtoBottom(MG_Chi);
	// 	// for (n=1; n<SER_Layer; n++)
	// 	// {
	// 	// 	RestrictionB(MG_Chi,MG_Chi,n+1,1);
	// 	// }
	// 	break;
	// }

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
		for (i=1; i<=LayerGridX[1]; i++) 
			GetMGCell(i, j, 1).M_value[3] = GetMGCell(i, j, 1).M_value[0];



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
		Exchange(MG_Sou, n+1);

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
			//deprecated//May-25-tianhong
			// SendtoBottom(MG_Sou);
			// if(Rank == Worker) 
			// {
			// 	MG_BottomLayer(field);
			// }
			// BottomSendBack(MG_Phi);

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







