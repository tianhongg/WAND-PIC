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

	GridX = XGridN;
	GridY = YGridN;

	if( (p_domain()->Get_Nbeam()) > 0 )
	{ 
	SendSouSizeX = YGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in X directions: left and right
	SendSouSizeY = XGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in Y directions: up   and down
	}
	else
	{
	SendSouSizeX = YGridN * SOU_DIM * 2;   //send in X directions: left and right
	SendSouSizeY = XGridN * SOU_DIM * 2;   //send in Y directions: up   and down	
	}

	SendSourceXm = new double[SendSouSizeX*bufsize];
	SendSourceXp = new double[SendSouSizeX*bufsize];
	ReceSourceXm = new double[SendSouSizeX*bufsize];
	ReceSourceXp = new double[SendSouSizeX*bufsize];

	SendSourceYm = new double[SendSouSizeY*bufsize];
	SendSourceYp = new double[SendSouSizeY*bufsize];
	ReceSourceYm = new double[SendSouSizeY*bufsize];
	ReceSourceYp = new double[SendSouSizeY*bufsize];


};

void Commute::DoCommute(int what, int k)
{


	//===============================================================
	//=====================           Pack         ==================
	//===============================================================
	//pack the fields or source all together in order to send.
	DoPack(what, k);
	int ssx;
	int ssy;

	MultiGrid *p_Multi = NULL;

	switch(what)
	{

		case COMMU_S:
		ssx = SendSouSizeX;
		ssy = SendSouSizeY;
		break;

		case COMMU_F:
		ssx = GridY * WAK_DIM2;
		ssy = GridX * WAK_DIM2;
		break;

		case COMMU_A:
		ssx = GridY * 4;
		ssy = GridX * 4;
		break;

		case COMMU_MG_P:
		case COMMU_MG_R:
		p_Multi = p_domain()->p_MG();
		ssx = p_Multi->GetLayerGridY(k);
		ssy = p_Multi->GetLayerGridX(k);
		break;

		case COMMU_MG_P_C:
		case COMMU_MG_R_C:
		p_Multi = p_domain()->p_MG();
		ssx = p_Multi->GetLayerGridY(k)*2;
		ssy = p_Multi->GetLayerGridX(k)*2;
		break;
	


	}


		//===============================================================
		//=====================   ISend and IRecev     ==================
		//===============================================================
    	MPI_Request Request[8];
    	MPI_Status 	 Status[8];
 		MPI_Irecv(ReceSourceXm, ssx, MPI_DOUBLE, XmPE, 0, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(ReceSourceXp, ssx, MPI_DOUBLE, XpPE, 1, MPI_COMM_WORLD, &Request[1]);
		MPI_Irecv(ReceSourceYm, ssy, MPI_DOUBLE, YmPE, 2, MPI_COMM_WORLD, &Request[2]);
		MPI_Irecv(ReceSourceYp, ssy, MPI_DOUBLE, YpPE, 3, MPI_COMM_WORLD, &Request[3]);

		MPI_Isend(SendSourceXp, ssx, MPI_DOUBLE, XpPE, 0, MPI_COMM_WORLD, &Request[4]);
		MPI_Isend(SendSourceXm, ssx, MPI_DOUBLE, XmPE, 1, MPI_COMM_WORLD, &Request[5]);
 		MPI_Isend(SendSourceYp, ssy, MPI_DOUBLE, YpPE, 2, MPI_COMM_WORLD, &Request[6]);
		MPI_Isend(SendSourceYm, ssy, MPI_DOUBLE, YmPE, 3, MPI_COMM_WORLD, &Request[7]);

 		int ierr;
		ierr = MPI_Waitall(8, Request,Status);

		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
		if (ierr!=0) 
		{
   			char errtxt[200];
			for (int i=0; i<8; i++) 
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
		{ for (i = 0; i < ssx; i++) {ReceSourceXm[i] = 0.0;}};

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i < ssx; i++) {ReceSourceXp[i] = 0.0;}};

		if (RankIdx_Y == 1) 
		{ for (i = 0; i < ssy; i++) {ReceSourceYm[i] = 0.0;}};

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i < ssy; i++) {ReceSourceYp[i] = 0.0;}};
		break;

	}

	//===============================================================
	//=====================           UnPack       ==================
	//===============================================================
	//unpack the fields or source.
	UnPack(what, k);

	return;

}


void Commute::DoCommuteT(int what,
int &Sendxm, int &Sendxp, int &Sendym, int &Sendyp)
{

	//===============================================================
	//=====================   ISend and IRecev     ==================
	//===============================================================
	int SendDIM;
	int Recexm, Receym, Recexp, Receyp;
    MPI_Request Request[8];
    MPI_Status 	 Status[8];
 	MPI_Irecv(&Recexm, 1, MPI_INT, XmPE, 0, MPI_COMM_WORLD, &Request[0]);
	MPI_Irecv(&Recexp, 1, MPI_INT, XpPE, 1, MPI_COMM_WORLD, &Request[1]);
	MPI_Irecv(&Receym, 1, MPI_INT, YmPE, 2, MPI_COMM_WORLD, &Request[2]);
	MPI_Irecv(&Receyp, 1, MPI_INT, YpPE, 3, MPI_COMM_WORLD, &Request[3]);

	MPI_Isend(&Sendxp, 1, MPI_INT, XpPE, 0, MPI_COMM_WORLD, &Request[4]);
	MPI_Isend(&Sendxm, 1, MPI_INT, XmPE, 1, MPI_COMM_WORLD, &Request[5]);
 	MPI_Isend(&Sendyp, 1, MPI_INT, YpPE, 2, MPI_COMM_WORLD, &Request[6]);
	MPI_Isend(&Sendym, 1, MPI_INT, YmPE, 3, MPI_COMM_WORLD, &Request[7]);
 	int ierr;
	ierr = MPI_Waitall(8, Request,Status);

	//==========test=========================
	//if(Rank == 9) printf("%d=%d\n",Rank,Sendym);
	//if(Rank == 5) printf("%d=%d\n",Rank,Receyp);
	//==========test=========================

    MPI_Request Request2[8];
    MPI_Status 	 Status2[8];

    switch(what)
    {
    	case COMMU_T:
    	SendDIM = SDT_DIM;
    	break;
    	case COMMU_P:
    	SendDIM = SDP_DIM;
    	break;
    }

 	MPI_Irecv(ReceSourceXm, Recexm*SendDIM, MPI_DOUBLE, XmPE, 0, MPI_COMM_WORLD, &Request2[0]);
	MPI_Irecv(ReceSourceXp, Recexp*SendDIM, MPI_DOUBLE, XpPE, 1, MPI_COMM_WORLD, &Request2[1]);
	MPI_Irecv(ReceSourceYm, Receym*SendDIM, MPI_DOUBLE, YmPE, 2, MPI_COMM_WORLD, &Request2[2]);
	MPI_Irecv(ReceSourceYp, Receyp*SendDIM, MPI_DOUBLE, YpPE, 3, MPI_COMM_WORLD, &Request2[3]);

	MPI_Isend(SendSourceXp, Sendxp*SendDIM, MPI_DOUBLE, XpPE, 0, MPI_COMM_WORLD, &Request2[4]);
	MPI_Isend(SendSourceXm, Sendxm*SendDIM, MPI_DOUBLE, XmPE, 1, MPI_COMM_WORLD, &Request2[5]);
 	MPI_Isend(SendSourceYp, Sendyp*SendDIM, MPI_DOUBLE, YpPE, 2, MPI_COMM_WORLD, &Request2[6]);
	MPI_Isend(SendSourceYm, Sendym*SendDIM, MPI_DOUBLE, YmPE, 3, MPI_COMM_WORLD, &Request2[7]);
	ierr = MPI_Waitall(8, Request2,Status2);

	if (ierr!=0) 
	{
		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
   		char errtxt[200];
		for (int i=0; i<8; i++) 
		{
			int err = Status2[i].MPI_ERROR; 
			int len=200;
			MPI_Error_string(err,errtxt,&len);
			printf("%s; \n",errtxt);
		}
		MPI_Abort(MPI_COMM_WORLD,0);
	}


	UnPackT(what, Recexm, Recexp, Receym, Receyp);
	return;
}



void Commute::DoPack(int what, int k)

{

	Mesh *p_Meshs  = p_domain()->p_Mesh();
	int i,j,n,m;
	int LayerGridX;
	int LayerGridY;
	MultiGrid *p_Multi = NULL;

	switch (what)
	{

		//===============================================================
		//===================== Pack Plasma Source ======================
		//===============================================================
		case COMMU_S:

		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0) NSource=SOU_DIM+BEA_DIM;
		
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: Y: up and down
		for (m = 0 ; m<=1; m++)
		{
			for (i=0; i < GridX; i++)
			{

				Cell &cm = p_Meshs->GetCell(i+1,m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+1-m,k);

				for (n = 0; n < NSource; n++)
				{
					SendSourceYm[ (GridX*m+i)*NSource + n ] = cm.W_Source[n];
					SendSourceYp[ (GridX*m+i)*NSource + n ] = cp.W_Source[n];
				}

			}
		}
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: X: left and right
		for (m = 0 ; m<=1; m++)
		{
			for (j=0; j < GridY; j++)
			{

				Cell &cm = p_Meshs->GetCell(m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+1-m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
					SendSourceXm[ (GridY*m+j)*NSource + n] = cm.W_Source[n];
					SendSourceXp[ (GridY*m+j)*NSource + n] = cp.W_Source[n];
				}

			}
		}

		break;



		//===============================================================
		//===================== Pack Wakefields    ======================
		//===============================================================
		case COMMU_F:
			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	  1,k);
				Cell &cp = p_Meshs->GetCell(i,GridY,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					SendSourceYm[ (i-1)*WAK_DIM2 + n ] = cm.W_Fields[n+5];
					SendSourceYp[ (i-1)*WAK_DIM2 + n ] = cp.W_Fields[n+5];
				}

			}

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(1,	  j,k);
				Cell &cp = p_Meshs->GetCell(GridX,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					SendSourceXm[ (j-1)*WAK_DIM2 + n ] = cm.W_Fields[n+5];
					SendSourceXp[ (j-1)*WAK_DIM2 + n ] = cp.W_Fields[n+5];
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

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);
			SendSourceXm[j-1] = cxm.M_value[0];
			SendSourceXp[j-1] = cxp.M_value[0];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);


			SendSourceYm[i-1] = cym.M_value[0];
			SendSourceYp[i-1] = cyp.M_value[0];
		}
		break;


		case COMMU_MG_P_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			SendSourceXm[(j-1)*2+0] = (cxm.C_value[0]).real();
			SendSourceXm[(j-1)*2+1] = (cxm.C_value[0]).imag();

			SendSourceXp[(j-1)*2+0] = (cxp.C_value[0]).real();
			SendSourceXp[(j-1)*2+1] = (cxp.C_value[0]).imag();
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			SendSourceYm[(i-1)*2+0] = (cym.C_value[0]).real();
			SendSourceYm[(i-1)*2+1] = (cym.C_value[0]).imag();
			SendSourceYp[(i-1)*2+0] = (cyp.C_value[0]).real();
			SendSourceYp[(i-1)*2+1] = (cyp.C_value[0]).imag();
		}
		break;


		//===============================================================
		//===================== Pack Multigrid Source  ==================
		//===============================================================
		//Exchange multigrid residual
		case COMMU_MG_R:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);
			SendSourceXm[j-1] = cxm.M_value[2];
			SendSourceXp[j-1] = cxp.M_value[2];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);


			SendSourceYm[i-1] = cym.M_value[2];
			SendSourceYp[i-1] = cyp.M_value[2];
		}
		break;



		case COMMU_MG_R_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			SendSourceXm[(j-1)*2+0] = (cxm.C_value[2]).real();
			SendSourceXm[(j-1)*2+1] = (cxm.C_value[2]).imag();
			SendSourceXp[(j-1)*2+0] = (cxp.C_value[2]).real();
			SendSourceXp[(j-1)*2+1] = (cxp.C_value[2]).imag();
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			SendSourceYm[(i-1)*2+0] = (cym.C_value[2]).real();
			SendSourceYm[(i-1)*2+1] = (cym.C_value[2]).imag();
			SendSourceYp[(i-1)*2+0] = (cyp.C_value[2]).real();
			SendSourceYp[(i-1)*2+1] = (cyp.C_value[2]).imag();
		}
		break;


	}



	return;
}




void Commute::UnPack(int what, int k)
{
	Mesh *p_Meshs  = p_domain()->p_Mesh();

	int i,j,n,m;
	MultiGrid *p_Multi = NULL;
	int LayerGridX;
	int LayerGridY;

	switch (what)
	{
		//===============================================================
		//===================== Unpack Plasma Source  ===================
		//===============================================================
		case COMMU_S:
		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0) NSource=SOU_DIM+BEA_DIM;
		// Pull sources from the Rece Array, and add on the edging cells ;
		// Receive Direction: Y
		for (m = 0; m<=1; m++)
		{
			for (i = 0; i < GridX; i++)
			{

				Cell &cm = p_Meshs->GetCell(i+1,1-m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+m,k);
				for (n = 0; n < NSource; n++)
				{
					cm.W_Source[n] += ReceSourceYm[ (GridX*m+i)*NSource+ n ];
					cp.W_Source[n] += ReceSourceYp[ (GridX*m+i)*NSource+ n ];
				}

			}
		}

		for (m = 0; m<=1; m++)
		{
			for (j = 0; j < GridY; j++)
			{

				Cell &cm = p_Meshs->GetCell(1-m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
				cm.W_Source[n] += ReceSourceXm[ (GridY*m+j)*NSource + n];
				cp.W_Source[n] += ReceSourceXp[ (GridY*m+j)*NSource + n];
				}

			}
		}

		//======== if Conduction Boundary Condition ======
		if(p_domain()->Get_BC()==1)
		{	
			if (RankIdx_X == 1) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	Cell &c = p_Meshs->GetCell(1, i,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};

				}
			}

			if (RankIdx_X == Xpa) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	Cell &c = p_Meshs->GetCell(GridX,i,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			
			}

			if (RankIdx_Y == 1) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	Cell &c = p_Meshs->GetCell(i,1,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			}

			if (RankIdx_Y == Ypa) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	Cell &c = p_Meshs->GetCell(i,GridY,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			}

		}


		break;


		//===============================================================
		//===================== Unpack Wakefields     ===================
		//===============================================================
		case COMMU_F:
			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	    0,k);
				Cell &cp = p_Meshs->GetCell(i,GridY+1,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = ReceSourceYm[ (i-1)*WAK_DIM2 + n ];
					cp.W_Fields[n+5] = ReceSourceYp[ (i-1)*WAK_DIM2 + n ];
				}

			}

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(0,	    j,k);
				Cell &cp = p_Meshs->GetCell(GridX+1,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = ReceSourceXm[ (j-1)*WAK_DIM2 + n ];
					cp.W_Fields[n+5] = ReceSourceXp[ (j-1)*WAK_DIM2 + n ];
				}

			}
		break;

		//===============================================================
		//===================== Unpack Vector Potential   ===============
		//===============================================================







		//===============================================================
		//===================== Unpack Multigird Field ==================
		//===============================================================
		case COMMU_MG_P:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.M_value[0] = ReceSourceXm[j-1];
			cxp.M_value[0] = ReceSourceXp[j-1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.M_value[0] = ReceSourceYm[i-1];
			cyp.M_value[0] = ReceSourceYp[i-1];
		}


		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(0, 			  i, k)).M_value[0] = 0.0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(LayerGridX+1, i, k)).M_value[0] = 0.0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, 			  0, k)).M_value[0] = 0.0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, LayerGridY+1, k)).M_value[0] = 0.0; } };

		break;


		case COMMU_MG_P_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.C_value[0] = ReceSourceXm[(j-1)*2+0]+ci*ReceSourceXm[(j-1)*2+1];
			cxp.C_value[0] = ReceSourceXp[(j-1)*2+0]+ci*ReceSourceXp[(j-1)*2+1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.C_value[0] = ReceSourceYm[(i-1)*2+0]+ci*ReceSourceYm[(i-1)*2+1];
			cyp.C_value[0] = ReceSourceYp[(i-1)*2+0]+ci*ReceSourceYp[(i-1)*2+1];
		}

		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(0, 			  i, k)).C_value[0] = 0.0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(LayerGridX+1, i, k)).C_value[0] = 0.0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, 			  0, k)).C_value[0] = 0.0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, LayerGridY+1, k)).C_value[0] = 0.0; } };

		break;

		//===============================================================
		//===================== Unpack Multigird Source =================
		//===============================================================
		case COMMU_MG_R:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.M_value[2] = ReceSourceXm[j-1];
			cxp.M_value[2] = ReceSourceXp[j-1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.M_value[2] = ReceSourceYm[i-1];
			cyp.M_value[2] = ReceSourceYp[i-1];
		}



		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).M_value[2]			  = (p_Multi->GetMGCell(1, i, k)).M_value[2]; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).M_value[2] = (p_Multi->GetMGCell(LayerGridX, i, k)).M_value[2]; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).M_value[2] 			  = (p_Multi->GetMGCell(i, 1, k)).M_value[2]; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).M_value[2] = (p_Multi->GetMGCell(i, LayerGridY, k)).M_value[2]; } };

		break;


		case COMMU_MG_R_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.C_value[2] = ReceSourceXm[(j-1)*2+0]+ci*ReceSourceXm[(j-1)*2+1];
			cxp.C_value[2] = ReceSourceXp[(j-1)*2+0]+ci*ReceSourceXp[(j-1)*2+1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.C_value[2] = ReceSourceYm[(i-1)*2+0]+ci*ReceSourceYm[(i-1)*2+1];
			cyp.C_value[2] = ReceSourceYp[(i-1)*2+0]+ci*ReceSourceYp[(i-1)*2+1];
		}



		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).C_value[2]			  = (p_Multi->GetMGCell(1, i, k)).C_value[2]; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).C_value[2] = (p_Multi->GetMGCell(LayerGridX, i, k)).C_value[2]; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).C_value[2] 			  = (p_Multi->GetMGCell(i, 1, k)).C_value[2]; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).C_value[2] = (p_Multi->GetMGCell(i, LayerGridY, k)).C_value[2]; } };

		break;


	}




	return;
}



void Commute::UnPackT(int what, int Recexm, int Recexp, int Receym, int Receyp)
{

	int n;
	double x0,y0,z0,px,py,pz;
	double q2m, weight;
	double Ex0,Ey0,Ez0;
	int type;

	Trajectory *p =NULL;
	Particle  *pp =NULL;


	int TpCellx =  p_domain()->p_Mesh()->Get_TpCellx();
	int TpCelly =  p_domain()->p_Mesh()->Get_TpCelly();


	switch(what)
	{
		case COMMU_T:
		
		for(n=0; n<Recexm; n++)
		{
			x0 = ReceSourceXm[n*SDT_DIM + 2];
			y0 = ReceSourceXm[n*SDT_DIM + 3];
			z0 = ReceSourceXm[n*SDT_DIM + 4];

			p = new Trajectory(x0, y0, z0, TpCellx, TpCelly);
			p->x = ReceSourceXm[n*SDT_DIM + 0];
			p->y = ReceSourceXm[n*SDT_DIM + 1];
			p->Vx= ReceSourceXm[n*SDT_DIM + 5];
			p->Vy= ReceSourceXm[n*SDT_DIM + 6];

			p->old_x  = ReceSourceXm[n*SDT_DIM + 7];
			p->old_y  = ReceSourceXm[n*SDT_DIM + 8];
			p->old_vx = ReceSourceXm[n*SDT_DIM + 9];
			p->old_vy = ReceSourceXm[n*SDT_DIM + 10];

			p->Vxx = (p->Vx)*(p->Vx);
			p->Vyy = (p->Vy)*(p->Vy);
			p->Vxy = (p->Vx)*(p->Vy);
		}


		for(n=0; n<Recexp; n++)
		{
			x0 = ReceSourceXp[n*SDT_DIM + 2];
			y0 = ReceSourceXp[n*SDT_DIM + 3];
			z0 = ReceSourceXp[n*SDT_DIM + 4];

			p = new Trajectory(x0, y0, z0, TpCellx, TpCelly);
			p->x = ReceSourceXp[n*SDT_DIM + 0];
			p->y = ReceSourceXp[n*SDT_DIM + 1];
			p->Vx= ReceSourceXp[n*SDT_DIM + 5];
			p->Vy= ReceSourceXp[n*SDT_DIM + 6];

			p->old_x  = ReceSourceXp[n*SDT_DIM + 7];
			p->old_y  = ReceSourceXp[n*SDT_DIM + 8];
			p->old_vx = ReceSourceXp[n*SDT_DIM + 9];
			p->old_vy = ReceSourceXp[n*SDT_DIM + 10];

			p->Vxx = (p->Vx)*(p->Vx);
			p->Vyy = (p->Vy)*(p->Vy);
			p->Vxy = (p->Vx)*(p->Vy);
		}

		for(n=0; n<Receym; n++)
		{
			x0 = ReceSourceYm[n*SDT_DIM + 2];
			y0 = ReceSourceYm[n*SDT_DIM + 3];
			z0 = ReceSourceYm[n*SDT_DIM + 4];

			p = new Trajectory(x0, y0, z0, TpCellx, TpCelly);
			p->x = ReceSourceYm[n*SDT_DIM + 0];
			p->y = ReceSourceYm[n*SDT_DIM + 1];
			p->Vx= ReceSourceYm[n*SDT_DIM + 5];
			p->Vy= ReceSourceYm[n*SDT_DIM + 6];

			p->old_x  = ReceSourceYm[n*SDT_DIM + 7];
			p->old_y  = ReceSourceYm[n*SDT_DIM + 8];
			p->old_vx = ReceSourceYm[n*SDT_DIM + 9];
			p->old_vy = ReceSourceYm[n*SDT_DIM + 10];

			p->Vxx = (p->Vx)*(p->Vx);
			p->Vyy = (p->Vy)*(p->Vy);
			p->Vxy = (p->Vx)*(p->Vy);
		}

		for(n=0; n<Receyp; n++)
		{
			x0 = ReceSourceYp[n*SDT_DIM + 2];
			y0 = ReceSourceYp[n*SDT_DIM + 3];
			z0 = ReceSourceYp[n*SDT_DIM + 4];

			p = new Trajectory(x0, y0, z0, TpCellx, TpCelly);
			p->x = ReceSourceYp[n*SDT_DIM + 0];
			p->y = ReceSourceYp[n*SDT_DIM + 1];
			p->Vx= ReceSourceYp[n*SDT_DIM + 5];
			p->Vy= ReceSourceYp[n*SDT_DIM + 6];

			p->old_x  = ReceSourceYp[n*SDT_DIM + 7];
			p->old_y  = ReceSourceYp[n*SDT_DIM + 8];
			p->old_vx = ReceSourceYp[n*SDT_DIM + 9];
			p->old_vy = ReceSourceYp[n*SDT_DIM + 10];

			p->Vxx = (p->Vx)*(p->Vx);
			p->Vyy = (p->Vy)*(p->Vy);
			p->Vxy = (p->Vx)*(p->Vy);
		}



		break;



		case COMMU_P:
		

		for(n=0; n<Recexm; n++)
		{
			x0 = ReceSourceXm[n*SDP_DIM + 3];
			y0 = ReceSourceXm[n*SDP_DIM + 4];
			z0 = ReceSourceXm[n*SDP_DIM + 5];
			px = ReceSourceXm[n*SDP_DIM + 6];
			py = ReceSourceXm[n*SDP_DIM + 7];
			pz = ReceSourceXm[n*SDP_DIM + 8];
			Ex0= ReceSourceXm[n*SDP_DIM + 9];
			Ey0= ReceSourceXm[n*SDP_DIM +10];
			Ez0= ReceSourceXm[n*SDP_DIM +11];

			type=(int)ReceSourceXm[n*SDP_DIM +12];
			q2m    =  ReceSourceXm[n*SDP_DIM +13];
			weight =  ReceSourceXm[n*SDP_DIM +14];

			switch(type)
			{
			case ELECTRON:
			pp = new Electron(x0, y0, z0, px, py, pz,
							Ex0, Ey0, Ez0, q2m, weight);
			break;

			case ION:
			break;
			}

			pp->x = ReceSourceXm[n*SDP_DIM + 0];
			pp->y = ReceSourceXm[n*SDP_DIM + 1];
			pp->z = ReceSourceXm[n*SDP_DIM + 2];

			pp->Wxw = ReceSourceXm[n*SDP_DIM + 15];
			pp->Wyw = ReceSourceXm[n*SDP_DIM + 16];
			pp->Wzw = ReceSourceXm[n*SDP_DIM + 17];
			pp->Wxl = ReceSourceXm[n*SDP_DIM + 18];
			pp->Wyl = ReceSourceXm[n*SDP_DIM + 19];
			pp->Wzl = ReceSourceXm[n*SDP_DIM + 20];


		}


		for(n=0; n<Recexp; n++)
		{
			x0 = ReceSourceXp[n*SDP_DIM + 3];
			y0 = ReceSourceXp[n*SDP_DIM + 4];
			z0 = ReceSourceXp[n*SDP_DIM + 5];
			px = ReceSourceXp[n*SDP_DIM + 6];
			py = ReceSourceXp[n*SDP_DIM + 7];
			pz = ReceSourceXp[n*SDP_DIM + 8];
			Ex0= ReceSourceXp[n*SDP_DIM + 9];
			Ey0= ReceSourceXp[n*SDP_DIM +10];
			Ez0= ReceSourceXp[n*SDP_DIM +11];

			type=(int)ReceSourceXp[n*SDP_DIM +12];
			q2m=      ReceSourceXp[n*SDP_DIM +13];
			weight =  ReceSourceXp[n*SDP_DIM +14];
			switch(type)
			{
			case ELECTRON:
			pp = new Electron(x0, y0, z0, px, py, pz,
							Ex0, Ey0, Ez0, q2m, weight);
			break;

			case ION:
			break;
			}

			pp->x = ReceSourceXp[n*SDP_DIM + 0];
			pp->y = ReceSourceXp[n*SDP_DIM + 1];
			pp->z = ReceSourceXp[n*SDP_DIM + 2];

			pp->Wxw = ReceSourceXp[n*SDP_DIM + 15];
			pp->Wyw = ReceSourceXp[n*SDP_DIM + 16];
			pp->Wzw = ReceSourceXp[n*SDP_DIM + 17];
			pp->Wxl = ReceSourceXp[n*SDP_DIM + 18];
			pp->Wyl = ReceSourceXp[n*SDP_DIM + 19];
			pp->Wzl = ReceSourceXp[n*SDP_DIM + 20];
		}

		for(n=0; n<Receym; n++)
		{
			x0 = ReceSourceYm[n*SDP_DIM + 3];
			y0 = ReceSourceYm[n*SDP_DIM + 4];
			z0 = ReceSourceYm[n*SDP_DIM + 5];
			px = ReceSourceYm[n*SDP_DIM + 6];
			py = ReceSourceYm[n*SDP_DIM + 7];
			pz = ReceSourceYm[n*SDP_DIM + 8];
			Ex0= ReceSourceYm[n*SDP_DIM + 9];
			Ey0= ReceSourceYm[n*SDP_DIM +10];
			Ez0= ReceSourceYm[n*SDP_DIM +11];

			type=(int)ReceSourceYm[n*SDP_DIM +12];
			q2m =     ReceSourceYm[n*SDP_DIM +13];
			weight =  ReceSourceYm[n*SDP_DIM +14];
			switch(type)
			{
			case ELECTRON:
			pp = new Electron(x0, y0, z0, px, py, pz,
							Ex0, Ey0, Ez0, q2m, weight);
			break;

			case ION:
			break;
			}

			pp->x = ReceSourceYm[n*SDP_DIM + 0];
			pp->y = ReceSourceYm[n*SDP_DIM + 1];
			pp->z = ReceSourceYm[n*SDP_DIM + 2];

			pp->Wxw = ReceSourceYm[n*SDP_DIM + 15];
			pp->Wyw = ReceSourceYm[n*SDP_DIM + 16];
			pp->Wzw = ReceSourceYm[n*SDP_DIM + 17];
			pp->Wxl = ReceSourceYm[n*SDP_DIM + 18];
			pp->Wyl = ReceSourceYm[n*SDP_DIM + 19];
			pp->Wzl = ReceSourceYm[n*SDP_DIM + 20];
		}

		for(n=0; n<Receyp; n++)
		{
			x0 = ReceSourceYp[n*SDP_DIM + 3];
			y0 = ReceSourceYp[n*SDP_DIM + 4];
			z0 = ReceSourceYp[n*SDP_DIM + 5];
			px = ReceSourceYp[n*SDP_DIM + 6];
			py = ReceSourceYp[n*SDP_DIM + 7];
			pz = ReceSourceYp[n*SDP_DIM + 8];
			Ex0= ReceSourceYp[n*SDP_DIM + 9];
			Ey0= ReceSourceYp[n*SDP_DIM +10];
			Ez0= ReceSourceYp[n*SDP_DIM +11];

			type=(int)ReceSourceYp[n*SDP_DIM +12];
			q2m=      ReceSourceYp[n*SDP_DIM +13];
			weight =  ReceSourceYp[n*SDP_DIM +14];
			switch(type)
			{
			case ELECTRON:
			pp = new Electron(x0, y0, z0, px, py, pz,
							Ex0, Ey0, Ez0, q2m, weight);
			break;

			case ION:
			break;
			}

			pp->x = ReceSourceYp[n*SDP_DIM + 0];
			pp->y = ReceSourceYp[n*SDP_DIM + 1];
			pp->z = ReceSourceYp[n*SDP_DIM + 2];

			pp->Wxw = ReceSourceYp[n*SDP_DIM + 15];
			pp->Wyw = ReceSourceYp[n*SDP_DIM + 16];
			pp->Wzw = ReceSourceYp[n*SDP_DIM + 17];
			pp->Wxl = ReceSourceYp[n*SDP_DIM + 18];
			pp->Wyl = ReceSourceYp[n*SDP_DIM + 19];
			pp->Wzl = ReceSourceYp[n*SDP_DIM + 20];
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




