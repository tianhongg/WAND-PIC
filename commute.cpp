//----------------------------------------------------------------------------------||
//-------------------                commute.cpp                 -------------------||
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
//---Starting---------           : Feb-05-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||




#include "wand_PIC.h"



Commute::Commute(int XGridN, int YGridN)
{

	bufsize = ceil((p_domain()->Get_Buffersize())
		*p_domain()->p_Partition()->GetXpart());

	SendSourceXm = NULL;
	SendSourceXp = NULL;
	SendSourceYm = NULL;
	SendSourceYp = NULL;

	ReceSourceXm = NULL;
	ReceSourceXp = NULL;
	ReceSourceYm = NULL;
	ReceSourceYp = NULL;


	Rank = p_domain()->p_Partition()->Get_Rank();
	//Big Coordinates of the Rank
	RankIdx_X = p_domain()->p_Partition()->RankIdx_X(); 
	RankIdx_Y = p_domain()->p_Partition()->RankIdx_Y(); 

	//Defind the size of SendArray
	Xpa  = p_domain()->p_Partition()->GetXpart();
	Ypa  = p_domain()->p_Partition()->GetYpart();

	XmPE = p_domain()->p_Partition()->GetXmPE(); 
	XpPE = p_domain()->p_Partition()->GetXpPE(); 
	YmPE = p_domain()->p_Partition()->GetYmPE(); 
	YpPE = p_domain()->p_Partition()->GetYpPE(); 

	mmPE = p_domain()->p_Partition()->GetmmPE(); 
	mpPE = p_domain()->p_Partition()->GetmpPE(); 
	pmPE = p_domain()->p_Partition()->GetpmPE(); 
	ppPE = p_domain()->p_Partition()->GetppPE(); 

	GridX = XGridN;
	GridY = YGridN;

	kold=-10;

	if( (p_domain()->Get_Nbeam()) > 0 )
	{ 
		SendSouSizeX = YGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in X directions: left and right
		SendSouSizeY = XGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in Y directions: up   and down
	}
	else
	{ 
		//SOU_DIM = denn jx jy jxx jyy jxy
		SendSouSizeX = YGridN * SOU_DIM * 2;   //send in X directions: left and right
		SendSouSizeY = XGridN * SOU_DIM * 2;   //send in Y directions: up   and down	
	}

	SendSourceXm = new WDOUBLE[SendSouSizeX*bufsize];
	SendSourceXp = new WDOUBLE[SendSouSizeX*bufsize];

	ReceSourceXm = new WDOUBLE[SendSouSizeX*bufsize];
	ReceSourceXp = new WDOUBLE[SendSouSizeX*bufsize];

	SendSourceYm = new WDOUBLE[SendSouSizeY*bufsize];
	SendSourceYp = new WDOUBLE[SendSouSizeY*bufsize];

	ReceSourceYm = new WDOUBLE[SendSouSizeY*bufsize];
	ReceSourceYp = new WDOUBLE[SendSouSizeY*bufsize];

	// diagonal direction
	SendSourcemm = new WDOUBLE[SendSouSizeX/YGridN*bufsize];
	SendSourcemp = new WDOUBLE[SendSouSizeX/YGridN*bufsize];

	ReceSourcemm = new WDOUBLE[SendSouSizeX/YGridN*bufsize];
	ReceSourcemp = new WDOUBLE[SendSouSizeX/YGridN*bufsize];

	SendSourcepm = new WDOUBLE[SendSouSizeY/XGridN*bufsize];
	SendSourcepp = new WDOUBLE[SendSouSizeY/XGridN*bufsize];

	ReceSourcepm = new WDOUBLE[SendSouSizeY/XGridN*bufsize];
	ReceSourcepp = new WDOUBLE[SendSouSizeY/XGridN*bufsize];


	//Cell Position Accumulative
	CellAccX = std::vector<WDOUBLE> (XGridN+3,0.0);
	CellAccY = std::vector<WDOUBLE> (YGridN+3,0.0);

	for(int i=0;i<XGridN+2;i++)
	{	
		Cell &ccc = p_domain()->p_Mesh()->GetCell(i,0,0);
		CellAccX[i]= ccc.Xcord-ccc.dx*0.5;
	}
	Cell &c1 = p_domain()->p_Mesh()->GetCell(XGridN+1,0,0);
	CellAccX[XGridN+2]=c1.Xcord+c1.dx*0.5;

	for(int i=0;i<YGridN+2;i++)
	{	
		Cell &ccc = p_domain()->p_Mesh()->GetCell(0,i,0);
		CellAccY[i]= ccc.Ycord-ccc.dy*0.5;
	}
	Cell &c2 = p_domain()->p_Mesh()->GetCell(0,YGridN+1,0);
	CellAccY[YGridN+2]=c2.Ycord+c2.dy*0.5;


};

void Commute::DoCommute(exchange what, int k)
{


	//===============================================================
	//=====================           Pack         ==================
	//===============================================================
	//pack the fields or source all together in order to send.
	DoPack(what, k);
	int ssx;
	int ssy;

	int ssxd;
	int ssyd;

	MultiGrid *p_Multi = NULL;

	switch(what)
	{

		case COMMU_S:
		case COMMU_SO:
			ssx = SendSouSizeX;
			ssy = SendSouSizeY;
			ssxd= ssx/GridY;
			ssyd= ssy/GridX;
		break;

		case COMMU_F:
			ssx = GridY * WAK_DIM2;
			ssy = GridX * WAK_DIM2;
			ssxd= WAK_DIM2;
			ssyd= WAK_DIM2;
		break;

		case COMMU_A:
			ssx = GridY * 4;
			ssy = GridX * 4;
		break;

		case COMMU_MG_P:
		case COMMU_MG_S:
		case COMMU_MG_R:
		case COMMU_MG_C:
			p_Multi = p_domain()->p_MG();
			ssx = p_Multi->GetLayerGridY(k);
			ssy = p_Multi->GetLayerGridX(k);
			ssxd= 1;
			ssyd= 1;
		break;

		case COMMU_MG_P_C:
		case COMMU_MG_S_C:
		case COMMU_MG_R_C:
		case COMMU_MG_C_C:
			p_Multi = p_domain()->p_MG();
			ssx = p_Multi->GetLayerGridY(k)*2;
			ssy = p_Multi->GetLayerGridX(k)*2;
			ssxd= 2;
			ssyd= 2;
		break;
	}
		//===============================================================
		//=====================   ISend and IRecev     ==================
		//===============================================================
    	MPI_Request Request[16];
    	MPI_Status 	 Status[16];

 		MPI_Irecv(ReceSourceXm, ssx,  MPI_WDOUBLE, XmPE, 0, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(ReceSourceXp, ssx,  MPI_WDOUBLE, XpPE, 1, MPI_COMM_WORLD, &Request[1]);
		MPI_Irecv(ReceSourceYm, ssy,  MPI_WDOUBLE, YmPE, 2, MPI_COMM_WORLD, &Request[2]);
		MPI_Irecv(ReceSourceYp, ssy,  MPI_WDOUBLE, YpPE, 3, MPI_COMM_WORLD, &Request[3]);

		MPI_Irecv(ReceSourcemm, ssxd, MPI_WDOUBLE, mmPE, 4, MPI_COMM_WORLD, &Request[4]);
		MPI_Irecv(ReceSourcemp, ssxd, MPI_WDOUBLE, mpPE, 5, MPI_COMM_WORLD, &Request[5]);
		MPI_Irecv(ReceSourcepm, ssyd, MPI_WDOUBLE, pmPE, 6, MPI_COMM_WORLD, &Request[6]);
		MPI_Irecv(ReceSourcepp, ssyd, MPI_WDOUBLE, ppPE, 7, MPI_COMM_WORLD, &Request[7]);

		MPI_Isend(SendSourceXp, ssx,  MPI_WDOUBLE, XpPE, 0, MPI_COMM_WORLD, &Request[8]);
		MPI_Isend(SendSourceXm, ssx,  MPI_WDOUBLE, XmPE, 1, MPI_COMM_WORLD, &Request[9]);
 		MPI_Isend(SendSourceYp, ssy,  MPI_WDOUBLE, YpPE, 2, MPI_COMM_WORLD, &Request[10]);
		MPI_Isend(SendSourceYm, ssy,  MPI_WDOUBLE, YmPE, 3, MPI_COMM_WORLD, &Request[11]);

		MPI_Isend(SendSourcepp, ssxd, MPI_WDOUBLE, ppPE, 4, MPI_COMM_WORLD, &Request[12]);
		MPI_Isend(SendSourcepm, ssxd, MPI_WDOUBLE, pmPE, 5, MPI_COMM_WORLD, &Request[13]);
 		MPI_Isend(SendSourcemp, ssyd, MPI_WDOUBLE, mpPE, 6, MPI_COMM_WORLD, &Request[14]);
		MPI_Isend(SendSourcemm, ssyd, MPI_WDOUBLE, mmPE, 7, MPI_COMM_WORLD, &Request[15]);



 		int ierr;
		ierr = MPI_Waitall(16, Request,Status);

		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
		if (ierr!=0) 
		{
   			char errtxt[200];
			for (int i=0; i<16; i++) 
			{
				int err = Status[i].MPI_ERROR; 
				int len=200;
				MPI_Error_string(err,errtxt,&len);
				printf("%s; \n",errtxt);
			}
			MPI_Abort(MPI_COMM_WORLD,0);
		}

	//===============================================================
	//=============== For Conduction Boundary Condition =============
	//=============== This Part May Change for Periodic BC ==========
	//===============================================================
	switch(what)
	{
		case COMMU_F:

		int i = 0;
		if (RankIdx_X == 1) 
		{ 
			for (i = 0; i < ssx;  i++) {ReceSourceXm[i] = 0.0;}
			for (i = 0; i < ssxd; i++) {ReceSourcemm[i] = 0.0;ReceSourcemp[i] = 0.0;}
		};

		if (RankIdx_X == Xpa) 
		{ 
			for (i = 0; i < ssx;  i++) {ReceSourceXp[i] = 0.0;}
			for (i = 0; i < ssxd; i++) {ReceSourcepm[i] = 0.0;ReceSourcepp[i] = 0.0;}
		};

		if (RankIdx_Y == 1) 
		{ 
			for (i = 0; i < ssy;  i++) {ReceSourceYm[i] = 0.0;}
			for (i = 0; i < ssyd; i++) {ReceSourcemm[i] = 0.0;ReceSourcepm[i] = 0.0;}
		};

		if (RankIdx_Y == Ypa) 
		{
			for (i = 0; i < ssy;  i++) {ReceSourceYp[i] = 0.0;}
			for (i = 0; i < ssyd; i++) {ReceSourcemp[i] = 0.0;ReceSourcepp[i] = 0.0;}
		};
		break;

	}

	//===============================================================
	//=====================           UnPack       ==================
	//===============================================================
	//unpack the fields or source.
	UnPack(what, k);



	return;

}


void Commute::DoCommuteT(exchange what, std::vector<int> SendN)
{

	//===============================================================
	//=====================   ISend and IRecev     ==================
	//===============================================================
	int SendDIM;

	std::vector<int> ReceN(8,0);

    MPI_Request Request[16];
    MPI_Status 	 Status[16];

 	MPI_Irecv(&(*(ReceN.begin()+0)), 1, MPI_INT, mmPE, 0, MPI_COMM_WORLD, &Request[0]);
	MPI_Irecv(&(*(ReceN.begin()+1)), 1, MPI_INT, mpPE, 1, MPI_COMM_WORLD, &Request[1]);
	MPI_Irecv(&(*(ReceN.begin()+2)), 1, MPI_INT, pmPE, 2, MPI_COMM_WORLD, &Request[2]);
	MPI_Irecv(&(*(ReceN.begin()+3)), 1, MPI_INT, ppPE, 3, MPI_COMM_WORLD, &Request[3]);

	MPI_Irecv(&(*(ReceN.begin()+4)), 1, MPI_INT, XmPE, 4, MPI_COMM_WORLD, &Request[4]);
	MPI_Irecv(&(*(ReceN.begin()+5)), 1, MPI_INT, XpPE, 5, MPI_COMM_WORLD, &Request[5]);
	MPI_Irecv(&(*(ReceN.begin()+6)), 1, MPI_INT, YmPE, 6, MPI_COMM_WORLD, &Request[6]);
	MPI_Irecv(&(*(ReceN.begin()+7)), 1, MPI_INT, YpPE, 7, MPI_COMM_WORLD, &Request[7]);
	//-------
	MPI_Isend(&(*(SendN.begin()+3)), 1, MPI_INT, ppPE, 0, MPI_COMM_WORLD, &Request[8]);
	MPI_Isend(&(*(SendN.begin()+2)), 1, MPI_INT, pmPE, 1, MPI_COMM_WORLD, &Request[9]);
 	MPI_Isend(&(*(SendN.begin()+1)), 1, MPI_INT, mpPE, 2, MPI_COMM_WORLD, &Request[10]);
	MPI_Isend(&(*(SendN.begin()+0)), 1, MPI_INT, mmPE, 3, MPI_COMM_WORLD, &Request[11]);

	MPI_Isend(&(*(SendN.begin()+5)), 1, MPI_INT, XpPE, 4, MPI_COMM_WORLD, &Request[12]);
	MPI_Isend(&(*(SendN.begin()+4)), 1, MPI_INT, XmPE, 5, MPI_COMM_WORLD, &Request[13]);
 	MPI_Isend(&(*(SendN.begin()+7)), 1, MPI_INT, YpPE, 6, MPI_COMM_WORLD, &Request[14]);
	MPI_Isend(&(*(SendN.begin()+6)), 1, MPI_INT, YmPE, 7, MPI_COMM_WORLD, &Request[15]);


 	int ierr;
	ierr = MPI_Waitall(16, Request,Status);

	//==========test=========================
	//if(Rank == 9) printf("%d=%d\n",Rank,Sendym);
	//if(Rank == 5) printf("%d=%d\n",Rank,Receyp);
	//==========test=========================

    MPI_Request Request2[16];
    MPI_Status 	 Status2[16];

    switch(what)
    {
    	case COMMU_T:
    	SendDIM = SDT_DIM;
    	break;
    	case COMMU_P:
    	SendDIM = SDP_DIM;
    	break;
    }

    MPI_Irecv(ReceSourcemm, ReceN[0]*SendDIM, MPI_WDOUBLE, mmPE, 0, MPI_COMM_WORLD, &Request2[0]);
	MPI_Irecv(ReceSourcemp, ReceN[1]*SendDIM, MPI_WDOUBLE, mpPE, 1, MPI_COMM_WORLD, &Request2[1]);
	MPI_Irecv(ReceSourcepm, ReceN[2]*SendDIM, MPI_WDOUBLE, pmPE, 2, MPI_COMM_WORLD, &Request2[2]);
	MPI_Irecv(ReceSourcepp, ReceN[3]*SendDIM, MPI_WDOUBLE, ppPE, 3, MPI_COMM_WORLD, &Request2[3]);
 	
 	MPI_Irecv(ReceSourceXm, ReceN[4]*SendDIM, MPI_WDOUBLE, XmPE, 4, MPI_COMM_WORLD, &Request2[4]);
	MPI_Irecv(ReceSourceXp, ReceN[5]*SendDIM, MPI_WDOUBLE, XpPE, 5, MPI_COMM_WORLD, &Request2[5]);
	MPI_Irecv(ReceSourceYm, ReceN[6]*SendDIM, MPI_WDOUBLE, YmPE, 6, MPI_COMM_WORLD, &Request2[6]);
	MPI_Irecv(ReceSourceYp, ReceN[7]*SendDIM, MPI_WDOUBLE, YpPE, 7, MPI_COMM_WORLD, &Request2[7]);

	MPI_Isend(SendSourcepp, SendN[3]*SendDIM, MPI_WDOUBLE, ppPE, 0, MPI_COMM_WORLD, &Request2[8]);
	MPI_Isend(SendSourcepm, SendN[2]*SendDIM, MPI_WDOUBLE, pmPE, 1, MPI_COMM_WORLD, &Request2[9]);
 	MPI_Isend(SendSourcemp, SendN[1]*SendDIM, MPI_WDOUBLE, mpPE, 2, MPI_COMM_WORLD, &Request2[10]);
	MPI_Isend(SendSourcemm, SendN[0]*SendDIM, MPI_WDOUBLE, mmPE, 3, MPI_COMM_WORLD, &Request2[11]);
	
	MPI_Isend(SendSourceXp, SendN[5]*SendDIM, MPI_WDOUBLE, XpPE, 4, MPI_COMM_WORLD, &Request2[12]);
	MPI_Isend(SendSourceXm, SendN[4]*SendDIM, MPI_WDOUBLE, XmPE, 5, MPI_COMM_WORLD, &Request2[13]);
 	MPI_Isend(SendSourceYp, SendN[7]*SendDIM, MPI_WDOUBLE, YpPE, 6, MPI_COMM_WORLD, &Request2[14]);
	MPI_Isend(SendSourceYm, SendN[6]*SendDIM, MPI_WDOUBLE, YmPE, 7, MPI_COMM_WORLD, &Request2[15]);

	ierr = MPI_Waitall(16, Request2,Status2);

	if (ierr!=0) 
	{
		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
   		char errtxt[200];
		for (int i=0; i<16; i++) 
		{
			int err = Status2[i].MPI_ERROR; 
			int len=200;
			MPI_Error_string(err,errtxt,&len);
			printf("%s; \n",errtxt);
		}
		MPI_Abort(MPI_COMM_WORLD,0);
	}


	UnPackT(what, ReceN);
	return;
}

		//      ___________
		//     |mp | yp| pp|
		//     |___|___|___|
		//     |xm |   | xp|
		//     |___|___|___|
		//     |mm | ym| pm|
		//     |___|___|___|
		//

void Commute::DoPack(exchange what, int k)

{

	Mesh *p_Meshs  = p_domain()->p_Mesh();
	int i,j,n,m;
	int LayerGridX;
	int LayerGridY;
	MultiGrid *p_Multi = NULL;


	WDOUBLE* SeXm=SendSourceXm;
	WDOUBLE* SeXp=SendSourceXp;
	WDOUBLE* SeYm=SendSourceYm;
	WDOUBLE* SeYp=SendSourceYp;

	WDOUBLE* Semm=SendSourcemm;
	WDOUBLE* Semp=SendSourcemp;
	WDOUBLE* Sepm=SendSourcepm;
	WDOUBLE* Sepp=SendSourcepp;

	switch (what)
	{

		//===============================================================
		//===================== Pack Plasma Source ======================
		//===============================================================
		case COMMU_S:
		case COMMU_SO:

		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0&&kold!=k) NSource=SOU_DIM+BEA_DIM;
		
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: Y: up and down
		for (m = 0 ; m<=1; m++)
		{
			Cell &mp = p_Meshs->GetCell(m,GridY+1-m,k); 
			Cell &pp = p_Meshs->GetCell(GridX+1-m,GridY+1-m,k);

			for (i=0; i < GridX; i++)
			{

				Cell &cm = p_Meshs->GetCell(i+1,m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+1-m,k);

				for (n = 0; n < NSource; n++)
				{
					*SeYm = cm.W_Source[n]; SeYm++;
					*SeYp = cp.W_Source[n]; SeYp++;
					//diagonal
					if(i==0)
					{
						*Semp = mp.W_Source[n]; Semp++;
						*Sepp = pp.W_Source[n]; Sepp++;
					}
				}

			}
		}
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: X: left and right
		for (m = 0 ; m<=1; m++)
		{
			Cell &mm = p_Meshs->GetCell(m,m,k);
			Cell &pm = p_Meshs->GetCell(GridX+1-m,m,k);

			for (j=0; j < GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+1-m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
					*SeXm = cm.W_Source[n]; SeXm++;
					*SeXp = cp.W_Source[n]; SeXp++;
					//diagonal
					if(j==0)
					{
						*Semm = mm.W_Source[n]; Semm++;
						*Sepm = pm.W_Source[n]; Sepm++;
					}
				}
			}
		}

		break;



		//===============================================================
		//===================== Pack Wakefields    ======================
		//===============================================================
		case COMMU_F:

			Cell &mp = p_Meshs->GetCell(1,GridY,k);
			Cell &pp = p_Meshs->GetCell(GridX,GridY,k);
			Cell &mm = p_Meshs->GetCell(1,1,k);
			Cell &pm = p_Meshs->GetCell(GridX,1,k);

			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	  1,k);
				Cell &cp = p_Meshs->GetCell(i,GridY,k);

				for (n = 0; n < WAK_DIM2; n++)
				{
					*SeYm = cm.W_Fields[n+5]; SeYm++;
					*SeYp = cp.W_Fields[n+5]; SeYp++;

					if(i==1)
					{
						*Semp = mp.W_Fields[n+5]; Semp++;
						*Sepp = pp.W_Fields[n+5]; Sepp++;
					}
				}
			}

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(1,	  j,k);
				Cell &cp = p_Meshs->GetCell(GridX,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					*SeXm = cm.W_Fields[n+5]; SeXm++;
					*SeXp = cp.W_Fields[n+5]; SeXp++;

					if(j==1)
					{
						*Semm = mm.W_Fields[n+5]; Semm++;
						*Sepm = pm.W_Fields[n+5]; Sepm++;
					}
				}

			}

		break;


		//===============================================================
		//===================== Pack Vector Potential  ==================
		//===============================================================




		//===============================================================
		//===================== Pack Multigrid Field  ===================
		//===============================================================
		//Exchange multigrid potential
		case COMMU_MG_P:
		case COMMU_MG_S:
		case COMMU_MG_R:
		case COMMU_MG_C:

		int field=what-COMMU_MG_P;

		p_Multi = p_domain()->p_MG();
		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			*SeXm = cxm.M_value[field]; SeXm++;
			*SeXp = cxp.M_value[field]; SeXp++;
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			*SeYm = cym.M_value[field]; SeYm++;
			*SeYp = cyp.M_value[field]; SeYp++;
		}

		MG_Cell &cmm = p_Multi->GetMGCell(1,		  1, 			 k);
		MG_Cell &cmp = p_Multi->GetMGCell(1, 		  LayerGridY,    k);
		MG_Cell &cpm = p_Multi->GetMGCell(LayerGridX, 1, 			 k);
		MG_Cell &cpp = p_Multi->GetMGCell(LayerGridX, LayerGridY,	 k);

		*Semm = cmm.M_value[field];
		*Semp = cmp.M_value[field];
		*Sepm = cpm.M_value[field];
		*Sepp = cpp.M_value[field];

		break;

		case COMMU_MG_P_C:
		case COMMU_MG_S_C:
		case COMMU_MG_R_C:
		case COMMU_MG_C_C:

		field=what-COMMU_MG_P_C;

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			*SeXm = (cxm.C_value[field]).real(); SeXm++;
			*SeXm = (cxm.C_value[field]).imag(); SeXm++;
			*SeXp = (cxp.C_value[field]).real(); SeXp++;
			*SeXp = (cxp.C_value[field]).imag(); SeXp++;
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			*SeYm = (cym.C_value[field]).real(); SeYm++;
			*SeYm = (cym.C_value[field]).imag(); SeYm++;
			*SeYp = (cyp.C_value[field]).real(); SeYp++;
			*SeYp = (cyp.C_value[field]).imag(); SeYp++;
		}

		MG_Cell &Cmm = p_Multi->GetMGCell(1,		 	1, 			k);
		MG_Cell &Cmp = p_Multi->GetMGCell(1, 			LayerGridY, k);
		MG_Cell &Cpm = p_Multi->GetMGCell(LayerGridX,	1, 			k);
		MG_Cell &Cpp = p_Multi->GetMGCell(LayerGridX, 	LayerGridY, k);

		*Semm = Cmm.C_value[field].real(); Semm++; *Semm = Cmm.C_value[field].imag();
		*Semp = Cmp.C_value[field].real(); Semp++; *Semp = Cmp.C_value[field].imag();
		*Sepm = Cpm.C_value[field].real(); Sepm++; *Sepm = Cpm.C_value[field].imag();
		*Sepp = Cpp.C_value[field].real(); Sepp++; *Sepp = Cpp.C_value[field].imag();

		break;

	}



	return;
}


void Commute::UnPack(exchange what, int k)
{
	Mesh *p_Meshs  = p_domain()->p_Mesh();

	int i,j,n,m;
	MultiGrid *p_Multi = NULL;
	int LayerGridX;
	int LayerGridY;

	WDOUBLE* ReXm=ReceSourceXm;
	WDOUBLE* ReXp=ReceSourceXp;
	WDOUBLE* ReYm=ReceSourceYm;
	WDOUBLE* ReYp=ReceSourceYp;

	WDOUBLE* Remm=ReceSourcemm;
	WDOUBLE* Remp=ReceSourcemp;
	WDOUBLE* Repm=ReceSourcepm;
	WDOUBLE* Repp=ReceSourcepp;

	switch (what)
	{
		//===============================================================
		//===================== Unpack Plasma Source  ===================
		//===============================================================
		case COMMU_S:
		case COMMU_SO:
		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0&&kold!=k) 
		{
			NSource=SOU_DIM+BEA_DIM;
			kold=k;
		}
		// Pull sources from the Rece Array, and add on the edging cells ;
		// Receive Direction: Y
		for (m = 0; m<=1; m++)
		{
			Cell &mm = p_Meshs->GetCell(1-m,1-m,k);
			Cell &pm = p_Meshs->GetCell(GridX+m,1-m,k);

			for (i = 0; i < GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i+1,1-m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+m,k);
				for (n = 0; n < NSource; n++)
				{
					cm.W_Source[n] += *ReYm; ReYm++; 
					cp.W_Source[n] += *ReYp; ReYp++; 
				
					if(i==0)
					{
						mm.W_Source[n] += *Remm; Remm++;
						pm.W_Source[n] += *Repm; Repm++;
					}
				}

			}
		}

		for (m = 0; m<=1; m++)
		{

			Cell &mp = p_Meshs->GetCell(1-m,    GridY+m, k);
			Cell &pp = p_Meshs->GetCell(GridX+m,GridY+m, k);

			for (j = 0; j < GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(1-m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
					cm.W_Source[n] += *ReXm; ReXm++;
					cp.W_Source[n] += *ReXp; ReXp++;

					if(j==0)
					{
						mp.W_Source[n] += *Remp; Remp++;
						pp.W_Source[n] += *Repp; Repp++;
					}
				}



			}
		}

		//======== if Conduction Boundary Condition ======
		if(p_domain()->Get_BC()==1)
		{	
			if (RankIdx_X == 1) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	
					for(int b=0;b<2;b++)
					{
						Cell &c = p_Meshs->GetCell(b, i,k); c.W_Denn  = c.W_Deni;
						for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
					}
					
				}
			}

			if (RankIdx_X == Xpa) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	
					for(int b=0;b<2;b++)
					{
						Cell &c = p_Meshs->GetCell(GridX+b,i,k); c.W_Denn  = c.W_Deni;
						for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
					}
				}
			
			}

			if (RankIdx_Y == 1) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	
					for(int b=0;b<2;b++)
					{
						Cell &c = p_Meshs->GetCell(i,b,k); c.W_Denn  = c.W_Deni;
						for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
					}
				}
			}

			if (RankIdx_Y == Ypa) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	
					for(int b=0;b<2;b++)
					{
						Cell &c = p_Meshs->GetCell(i,GridY+b,k); c.W_Denn  = c.W_Deni;
						for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
					}
				}
			}

		}
		break;


		//===============================================================
		//===================== Unpack Wakefields     ===================
		//===============================================================
		case COMMU_F:

			Cell &mm = p_Meshs->GetCell(0,0,k);
			Cell &pm = p_Meshs->GetCell(GridX+1,0,k);

			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	    0,k);
				Cell &cp = p_Meshs->GetCell(i,GridY+1,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = *ReYm; ReYm++;
					cp.W_Fields[n+5] = *ReYp; ReYp++;

					if(i==1)
					{
						mm.W_Fields[n+5] = *Remm; Remm++;
						pm.W_Fields[n+5] = *Repm; Repm++;
					}
				}
			}

			Cell &mp = p_Meshs->GetCell(0,GridY+1,k);
			Cell &pp = p_Meshs->GetCell(GridX+1,GridY+1,k);

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(0,	    j,k);
				Cell &cp = p_Meshs->GetCell(GridX+1,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = *ReXm; ReXm++;
					cp.W_Fields[n+5] = *ReXp; ReXp++;

					if(j==1)
					{
						mp.W_Fields[n+5] = *Remp; Remp++;
						pp.W_Fields[n+5] = *Repp; Repp++;
					}
				}

			}
		break;


		//===============================================================
		//===================== Unpack Multigird Field ==================
		//===============================================================
		case COMMU_MG_P:
		case COMMU_MG_S:
		case COMMU_MG_R:
		case COMMU_MG_C:

		int field=what-COMMU_MG_P;

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.M_value[field] = *ReXm; ReXm++;
			cxp.M_value[field] = *ReXp; ReXp++;
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.M_value[field] = *ReYm; ReYm++;
			cyp.M_value[field] = *ReYp; ReYp++;
		}

		MG_Cell &cmm = p_Multi->GetMGCell(0,		    0, k);
		MG_Cell &cmp = p_Multi->GetMGCell(0, LayerGridY+1, k);
		MG_Cell &cpm = p_Multi->GetMGCell(LayerGridX+1, 0, k);
		MG_Cell &cpp = p_Multi->GetMGCell(LayerGridX+1, LayerGridY+1, k);

		cmm.M_value[field]=*Remm;
		cmp.M_value[field]=*Remp;
		cpm.M_value[field]=*Repm;
		cpp.M_value[field]=*Repp;

		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).M_value[field]			  = 0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).M_value[field] = 0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).M_value[field] 			  = 0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).M_value[field] = 0; } };

		break;

		case COMMU_MG_P_C:
		case COMMU_MG_S_C:
		case COMMU_MG_R_C:
		case COMMU_MG_C_C:

		field=what-COMMU_MG_P_C;

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.C_value[field] = ReceSourceXm[(j-1)*2+0]+ci*ReceSourceXm[(j-1)*2+1];
			cxp.C_value[field] = ReceSourceXp[(j-1)*2+0]+ci*ReceSourceXp[(j-1)*2+1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.C_value[field] = ReceSourceYm[(i-1)*2+0]+ci*ReceSourceYm[(i-1)*2+1];
			cyp.C_value[field] = ReceSourceYp[(i-1)*2+0]+ci*ReceSourceYp[(i-1)*2+1];
		}

		MG_Cell &Cmm = p_Multi->GetMGCell(0,		  0, k);
		MG_Cell &Cmp = p_Multi->GetMGCell(0, LayerGridY+1, k);
		MG_Cell &Cpm = p_Multi->GetMGCell(LayerGridX+1, 0, k);
		MG_Cell &Cpp = p_Multi->GetMGCell(LayerGridX+1, LayerGridY+1, k);

		Cmm.C_value[field]=ReceSourcemm[0]+ci*ReceSourcemm[1];
		Cmp.C_value[field]=ReceSourcemp[0]+ci*ReceSourcemp[1];
		Cpm.C_value[field]=ReceSourcepm[0]+ci*ReceSourcepm[1];
		Cpp.C_value[field]=ReceSourcepp[0]+ci*ReceSourcepp[1];

		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).C_value[field]			  = 0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).C_value[field] = 0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).C_value[field] 			  = 0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).C_value[field] = 0; } };
		break;	
	}



	return;
}



void Commute::UnPackT(exchange what, std::vector<int> ReceN)
{

	int n;
	WDOUBLE x0,y0,z0,px,py,pz;
	WDOUBLE xt,yt,zt, sx, sy;
	WDOUBLE vx,vy,old_x,old_y,old_vx,old_vy; 
	WDOUBLE q2m, weight;
	WDOUBLE Ex0,Ey0,Ez0;
	WDOUBLE Wxw,Wyw,Wzw;
	WDOUBLE Wxl,Wyl,Wzl;
	int type;

	Trajectory *p =NULL;
	Particle  *pp =NULL;


	int TpCellx =  p_domain()->p_Mesh()->Get_TpCellx();
	int TpCelly =  p_domain()->p_Mesh()->Get_TpCelly();


	switch(what)
	{
		
		case COMMU_T: // 
		

		for(int dir=0; dir<8; dir++)
		{
			WDOUBLE* Re=NULL;

			if(dir==0) Re=ReceSourcemm;
			if(dir==1) Re=ReceSourcemp;
			if(dir==2) Re=ReceSourcepm;
			if(dir==3) Re=ReceSourcepp;
			if(dir==4) Re=ReceSourceXm;
			if(dir==5) Re=ReceSourceXp;
			if(dir==6) Re=ReceSourceYm;
			if(dir==7) Re=ReceSourceYp;

			for(n=0; n<ReceN[dir]; n++)
			{
				xt = *Re; Re++; 
				yt = *Re; Re++;

				x0 = *Re; Re++; 
				y0 = *Re; Re++; 
				z0 = *Re; Re++; 

				vx = *Re; Re++; 
				vy = *Re; Re++; 

				old_x  = *Re; Re++; 
				old_y  = *Re; Re++; 
				old_vx = *Re; Re++; 
				old_vy = *Re; Re++; 

				sx = *Re; Re++; 
				sy = *Re; Re++; 

				p = new Trajectory(x0, y0, z0, TpCellx, TpCelly,sx,sy);

				p->x = xt;
				p->y = yt;
				p->Vx= vx;
				p->Vy= vy;

				p->old_x  = old_x;
				p->old_y  = old_y;
				p->old_vx = old_vx;
				p->old_vy = old_vy;

				p->Vxx = (p->Vx)*(p->Vx);
				p->Vyy = (p->Vy)*(p->Vy);
				p->Vxy = (p->Vx)*(p->Vy);

				auto upper=std::upper_bound(CellAccX.begin(),CellAccX.end(),xt);
				p->idx_i= (upper-CellAccX.begin()-1);
				upper=std::upper_bound(CellAccY.begin(),CellAccY.end(),yt);
				p->idx_j=(upper-CellAccY.begin()-1);
			}
		}
		break; //COMMU_T: // 



		case COMMU_P:

		for(int dir=0; dir<8; dir++)
		{
			WDOUBLE* Re=NULL;

			if(dir==0) Re=ReceSourcemm;
			if(dir==1) Re=ReceSourcemp;
			if(dir==2) Re=ReceSourcepm;
			if(dir==3) Re=ReceSourcepp;
			if(dir==4) Re=ReceSourceXm;
			if(dir==5) Re=ReceSourceXp;
			if(dir==6) Re=ReceSourceYm;
			if(dir==7) Re=ReceSourceYp;

			for(n=0; n<ReceN[dir]; n++)
			{
				xt = *Re; Re++; 
				yt = *Re; Re++;
				zt = *Re; Re++;

				x0 = *Re; Re++; 
				y0 = *Re; Re++; 
				z0 = *Re; Re++; 

				px = *Re; Re++; 
				py = *Re; Re++; 
				pz = *Re; Re++; 

				Ex0 = *Re; Re++; 
				Ey0 = *Re; Re++; 
				Ez0 = *Re; Re++; 

				type=(int)*Re; Re++; 
				q2m = 	  *Re; Re++; 
				weight =  *Re; Re++; 

				Wxw = *Re; Re++; 
				Wyw = *Re; Re++; 
				Wzw = *Re; Re++; 

				Wxl = *Re; Re++; 
				Wyl = *Re; Re++; 
				Wzl = *Re; Re++; 

				sx = *Re; Re++; 
				sy = *Re; Re++; 

				switch(type)
				{
					case ELECTRON:
						pp = new Electron(x0, y0, z0, px, py, pz, Ex0, Ey0, Ez0, q2m, weight);
					break;

					case ION:
						pp = new Ion(x0, y0, z0, px, py, pz, Ex0, Ey0, Ez0, q2m, weight);
					break;
				}

				pp->x =xt;
				pp->y =yt;
				pp->z =zt;

				pp->Wxw = Wxw;
				pp->Wyw = Wyw;
				pp->Wzw = Wzw;
				pp->Wxl = Wxl;
				pp->Wyl = Wyl;
				pp->Wzl = Wzl;

				pp->sx=sx;
				pp->sy=sy;

				auto upper=std::upper_bound(CellAccX.begin(),CellAccX.end(),xt);
				pp->idx_i= (upper-CellAccX.begin()-1);
				upper=std::upper_bound(CellAccY.begin(),CellAccY.end(),yt);
				pp->idx_j=(upper-CellAccY.begin()-1);

			}
		}
		break;
		
	

	}

	return;
}


Commute::~Commute()
{
	delete[] SendSourceXm;  //Array to send for the sources;
	delete[] SendSourceYm;  //Array to send for the sources;
   	delete[] SendSourceXp;  //Array to send for the sources;
  	delete[] SendSourceYp;  //Array to send for the sources;

   	delete[] ReceSourceXm;  //Array for Rece 
  	delete[] ReceSourceYm;  //Array for Rece 
	delete[] ReceSourceXp;  //Array for Rece 
	delete[] ReceSourceYp;  //Array for Rece 


}




