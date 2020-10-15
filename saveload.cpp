//----------------------------------------------------------------------------------||
//-------------------                saveload.cpp                -------------------||
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
//---Starting---------           : Feb-27-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||


#include <pnetcdf.h>
#include "wand_PIC.h"


int Domain::Save(int nt)
{

	char sFile[128];
	char vname1[24];
	char vname2[24];
	static const int NC_ERR = 1;

	int Xpa = p_Partition()->GetXpart();
	int Ypa = p_Partition()->GetYpart();

	int Idx_x = p_Partition()->RankIdx_X();
	int Idx_y = p_Partition()->RankIdx_Y();

	int i,j,k, NF;

	//========================================
	//========= Pnetcdf Variables ============
	int retval, ncid; 

	int nx_id,  ny_id,  nz_id, F_DIM[3];
	int n_e_id, n_n_id, Psi_id, WEx_id, WEy_id, WEz_id, WBx_id, WBy_id, WBz_id;
	int AxR_id[NFreqs],	AxI_id[NFreqs],	AyR_id[NFreqs],	AyI_id[NFreqs];
	int n_b_id;

	int LExR_id[NFreqs],	LEyR_id[NFreqs],	LEzR_id[NFreqs];
	int LBxR_id[NFreqs],	LByR_id[NFreqs],	LBzR_id[NFreqs];

	int LExI_id[NFreqs],	LEyI_id[NFreqs],	LEzI_id[NFreqs];
	int LBxI_id[NFreqs],	LByI_id[NFreqs],	LBzI_id[NFreqs];

	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Offset fstart[3], fcount[3];

	double  ne[XGridN][YGridN][ZGridN],  nn[XGridN][YGridN][ZGridN], Psi[XGridN][YGridN][ZGridN];
	double WEx[XGridN][YGridN][ZGridN], WEy[XGridN][YGridN][ZGridN], WEz[XGridN][YGridN][ZGridN];
	double WBx[XGridN][YGridN][ZGridN], WBy[XGridN][YGridN][ZGridN], WBz[XGridN][YGridN][ZGridN];
	double AxR[XGridN][YGridN][ZGridN], AxI[XGridN][YGridN][ZGridN];
	double AyR[XGridN][YGridN][ZGridN], AyI[XGridN][YGridN][ZGridN];

	double LExR[XGridN][YGridN][ZGridN], LExI[XGridN][YGridN][ZGridN];
	double LEyR[XGridN][YGridN][ZGridN], LEyI[XGridN][YGridN][ZGridN];
	double LEzR[XGridN][YGridN][ZGridN], LEzI[XGridN][YGridN][ZGridN];
	double LBxR[XGridN][YGridN][ZGridN], LBxI[XGridN][YGridN][ZGridN];
	double LByR[XGridN][YGridN][ZGridN], LByI[XGridN][YGridN][ZGridN];
	double LBzR[XGridN][YGridN][ZGridN], LBzI[XGridN][YGridN][ZGridN];

	double  nb[XGridN][YGridN][ZGridN];


	fstart[0] = (Idx_x-1)*XGridN;
	fstart[1] = (Idx_y-1)*YGridN;
	fstart[2] = 0;

	fcount[0] = XGridN;
	fcount[1] = YGridN;
	fcount[2] = ZGridN;


	//========================================
	//========= Put Fields for Output ========
	for(k=0; k<ZGridN; k++)
	{
		for(j=1; j<=YGridN; j++)
		{
			for(i=1; i<=XGridN; i++)
			{
				Cell &ccc = p_Mesh()->GetCell(i,j,k);

				 ne[i-1][j-1][k] = 0.5*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*ccc.W_Asq)/(1+ccc.W_Psi)/(1+ccc.W_Psi);
				 nn[i-1][j-1][k] = ccc.W_Denn;
				 nb[i-1][j-1][k] = ccc.B_Den;
				Psi[i-1][j-1][k] = ccc.W_Psi;

				WEx[i-1][j-1][k] = -ccc.W_Ex+ccc.W_By;
				WEy[i-1][j-1][k] = -ccc.W_Ey-ccc.W_Bx;
				WEz[i-1][j-1][k] = ccc.W_Ez;

				WBx[i-1][j-1][k] = ccc.W_Bx;
				WBy[i-1][j-1][k] = ccc.W_By;
				WBz[i-1][j-1][k] = ccc.W_Bz;

			}
		}
	}

	//========================================
	//========= Create nc File ===============
	sprintf(sFile,"Fields_%d.nc",nt);
	if((retval = ncmpi_create(MPI_COMM_WORLD, sFile, NC_CLOBBER|NC_64BIT_OFFSET, info, &ncid)) ) 
	{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while creating " << sFile << "." <<'\n';
    	return NC_ERR;
	}


	//========================================
	//========= Define Dimension =============
	if ( (retval = ncmpi_def_dim(ncid, "nx", XGridN*Xpa, &nx_id)) )	return NC_ERR;
	if ( (retval = ncmpi_def_dim(ncid, "ny", YGridN*Ypa, &ny_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_dim(ncid, "nz", ZGridN, 	 &nz_id)) ) return NC_ERR;

	F_DIM[0] = nx_id;
  	F_DIM[1] = ny_id;
  	F_DIM[2] = nz_id;

	//========================================
	//========= Define Variables =============
	if ( (retval = ncmpi_def_var(ncid, "n_e", NC_DOUBLE, 3, F_DIM, &n_e_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "n_n", NC_DOUBLE, 3, F_DIM, &n_n_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "Psi", NC_DOUBLE, 3, F_DIM, &Psi_id)) ) return NC_ERR;

	if ( (retval = ncmpi_def_var(ncid, "WEx", NC_DOUBLE, 3, F_DIM, &WEx_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WEy", NC_DOUBLE, 3, F_DIM, &WEy_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WEz", NC_DOUBLE, 3, F_DIM, &WEz_id)) ) return NC_ERR;

	if ( (retval = ncmpi_def_var(ncid, "WBx", NC_DOUBLE, 3, F_DIM, &WBx_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WBy", NC_DOUBLE, 3, F_DIM, &WBy_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WBz", NC_DOUBLE, 3, F_DIM, &WBz_id)) ) return NC_ERR;

	// Beam Density;
	if(Nbeam){
	if ( (retval = ncmpi_def_var(ncid, "n_b", NC_DOUBLE, 3, F_DIM, &n_b_id)) ) return NC_ERR;}

	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{
		if(ifAx[NF])
		{
		sprintf(vname1,"AxR_%d",NF);sprintf(vname2,"AxI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &AxR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &AxI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LExR_%d",NF);sprintf(vname2,"LExI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &LExR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &LExI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LByR_%d",NF);sprintf(vname2,"LByI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &LByR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &LByI_id[NF])) ) return NC_ERR;


		}

		if(ifAy[NF])
		{
		sprintf(vname1,"AyR_%d",NF);sprintf(vname2,"AyI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &AyR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &AyI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LEyR_%d",NF);sprintf(vname2,"LEyI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &LEyR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &LEyI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LBxR_%d",NF);sprintf(vname2,"LBxI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &LBxR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &LBxI_id[NF])) ) return NC_ERR;

		}

		sprintf(vname1,"LEzR_%d",NF);sprintf(vname2,"LEzI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 3, F_DIM, &LEzR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 3, F_DIM, &LEzI_id[NF])) ) return NC_ERR;
	}


	if ( (retval = ncmpi_enddef(ncid)) ) return NC_ERR;

	//========================================
	//============ Put Variables =============
	if ( (retval = ncmpi_put_vara_double_all(ncid, n_e_id, fstart, fcount,  &ne[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, n_n_id, fstart, fcount,  &nn[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, Psi_id, fstart, fcount, &Psi[0][0][0])) ) return NC_ERR;

  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEx_id, fstart, fcount, &WEx[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEy_id, fstart, fcount, &WEy[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEz_id, fstart, fcount, &WEz[0][0][0])) ) return NC_ERR;

  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBx_id, fstart, fcount, &WBx[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBy_id, fstart, fcount, &WBy[0][0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBz_id, fstart, fcount, &WBz[0][0][0])) ) return NC_ERR;
	
	if(Nbeam){
	if ( (retval = ncmpi_put_vara_double_all(ncid, n_b_id, fstart, fcount,  &nb[0][0][0])) ) return NC_ERR;}

	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{

		for(k=0; k<ZGridN; k++)
		{
			for(j=1; j<=YGridN; j++)
			{
				for(i=1; i<=XGridN; i++)
				{
					Cell &ccc = p_Mesh()->GetCell(i,j,k);
					if(ifAx[NF])
					{
						
						 AxR[i-1][j-1][k] = ccc.Acomx[NF].real();
						 AxI[i-1][j-1][k] = ccc.Acomx[NF].imag();

						LExR[i-1][j-1][k] =  ccc.L_Ex[NF].real();
						LExI[i-1][j-1][k] =  ccc.L_Ex[NF].imag();

						LByR[i-1][j-1][k] =  ccc.L_By[NF].real();
						LByI[i-1][j-1][k] =  ccc.L_By[NF].imag();

					}

					if(ifAy[NF])
					{
		
						 AyR[i-1][j-1][k] = ccc.Acomy[NF].real();
						 AyI[i-1][j-1][k] = ccc.Acomy[NF].imag();

						LEyR[i-1][j-1][k] =  ccc.L_Ey[NF].real();
						LEyI[i-1][j-1][k] =  ccc.L_Ey[NF].imag();

						LBxR[i-1][j-1][k] =  ccc.L_Bx[NF].real();
						LBxI[i-1][j-1][k] =  ccc.L_Bx[NF].imag();
					}

						LEzR[i-1][j-1][k] =  ccc.L_Ez[NF].real();
						LEzI[i-1][j-1][k] =  ccc.L_Ez[NF].imag();
				}
			}
		}


		if(ifAx[NF])
		{ 
			if ( (retval = ncmpi_put_vara_double_all(ncid, AxR_id[NF], fstart, fcount,   &AxR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, AxI_id[NF], fstart, fcount,   &AxI[0][0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LExR_id[NF], fstart, fcount, &LExR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LExI_id[NF], fstart, fcount, &LExI[0][0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LByR_id[NF], fstart, fcount, &LByR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LByI_id[NF], fstart, fcount, &LByI[0][0][0])) ) return NC_ERR;
	
		}
	
				
		if(ifAy[NF])
		{
			if ( (retval = ncmpi_put_vara_double_all(ncid, AyR_id[NF], fstart, fcount,   &AyR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, AyI_id[NF], fstart, fcount,   &AyI[0][0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LEyR_id[NF], fstart, fcount, &LEyR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LEyI_id[NF], fstart, fcount, &LEyI[0][0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LBxR_id[NF], fstart, fcount, &LBxR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LBxI_id[NF], fstart, fcount, &LBxI[0][0][0])) ) return NC_ERR;
		}

			if ( (retval = ncmpi_put_vara_double_all(ncid, LEzR_id[NF], fstart, fcount, &LEzR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LEzI_id[NF], fstart, fcount, &LEzI[0][0][0])) ) return NC_ERR;
	}



	if ( (retval = ncmpi_sync(ncid)) )
  	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval <<  " while syncing " << sFile << "." <<'\n';
		return NC_ERR; 
	}
  
	if ( (retval = ncmpi_close(ncid)) ) 
	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
		return NC_ERR;
	}

	if(Nbeam) SaveP(nt);
	SaveXray(nt);
	return 0;

}



int Domain::Save2D(int nt, int savedim)
{

	char sFile[128];
	char vname1[24];
	char vname2[24];

	int Xpa   = p_Partition()->GetXpart();
	int Ypa   = p_Partition()->GetYpart();
	int Idx_x = p_Partition()->RankIdx_X();
	int Idx_y = p_Partition()->RankIdx_Y();

	int i,j,k, NF;
	int IndexY, IndexX,nx, nnn;
	//========================================
	//========= Pnetcdf Variables ============
	static const int NC_ERR = 1;
	int retval, ncid; 
	int nx_id,  ny_id,  nz_id, F_DIM[2];
	int n_b_id;
	int n_e_id, n_n_id, Psi_id, WEx_id, WEy_id, WEz_id, WBx_id, WBy_id, WBz_id;

	int AxR_id[NFreqs],	AxI_id[NFreqs],	AyR_id[NFreqs],	AyI_id[NFreqs];

	int LExR_id[NFreqs],	LEyR_id[NFreqs],	LEzR_id[NFreqs];
	int LBxR_id[NFreqs],	LByR_id[NFreqs],	LBzR_id[NFreqs];

	int LExI_id[NFreqs],	LEyI_id[NFreqs],	LEzI_id[NFreqs];
	int LBxI_id[NFreqs],	LByI_id[NFreqs],	LBzI_id[NFreqs];

	


	switch(savedim)
	{

	case 1:
		nx=XGridN;
		IndexY = floor(Ypa/2.0)+1;
	break;

	case 2:
		nx=YGridN;
		IndexX = floor(Xpa/2.0)+1;
	break;

	}

	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Offset fstart[2], fcount[2];

	fcount[0]=0;
	fstart[0]=0;

	fcount[1] = 0;
	fstart[1] = 0;
	

	double  ne[nx][ZGridN],  nn[nx][ZGridN], Psi[nx][ZGridN];
	double WEx[nx][ZGridN], WEy[nx][ZGridN], WEz[nx][ZGridN];
	double WBx[nx][ZGridN], WBy[nx][ZGridN], WBz[nx][ZGridN];

	double AxR[nx][ZGridN], AxI[nx][ZGridN];
	double AyR[nx][ZGridN], AyI[nx][ZGridN];

	double LExR[nx][ZGridN], LEyR[nx][ZGridN], LEzR[nx][ZGridN];
	double LBxR[nx][ZGridN], LByR[nx][ZGridN], LBzR[nx][ZGridN];

	double LExI[nx][ZGridN], LEyI[nx][ZGridN], LEzI[nx][ZGridN];
	double LBxI[nx][ZGridN], LByI[nx][ZGridN], LBzI[nx][ZGridN];

	double  nb[nx][ZGridN];

	switch(savedim)
	{
	case 1:

		if(IndexY==Idx_y)
		{
			fstart[0] = (Idx_x-1)*XGridN;
			fcount[0] = XGridN;
			fcount[1] = ZGridN;
		}
		j = int(YGridN*(Ypa/2.0+1-IndexY))+1;
		nnn=j;
		//========================================
		//========= Put Fields for Output ========
		for(k=0; k<ZGridN; k++)
		{
			for(i=1; i<=XGridN; i++)
			{
				Cell &ccc = p_Mesh()->GetCell(i,j,k);
				 ne[i-1][k] = 0.5*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*ccc.W_Asq)/(1+ccc.W_Psi)/(1+ccc.W_Psi);
				 nn[i-1][k] = ccc.W_Denn;
				 nb[i-1][k] = ccc.B_Den;
				Psi[i-1][k] = ccc.W_Psi;

				WEx[i-1][k] = -ccc.W_Ex+ccc.W_By;
				WEy[i-1][k] = -ccc.W_Ey-ccc.W_Bx;
				WEz[i-1][k] = ccc.W_Ez;

				WBx[i-1][k] = ccc.W_Bx;
				WBy[i-1][k] = ccc.W_By;
				WBz[i-1][k] = ccc.W_Bz;
			}
		}


	break;

	case 2:

		if(IndexX==Idx_x)
		{
			fstart[0] = (Idx_y-1)*YGridN;;
			fcount[0] = YGridN;
			fcount[1] = ZGridN;
		}
		i = int(XGridN*(Xpa/2.0+1-IndexX))+1;
		nnn=i;

		//========================================
		//========= Put Fields for Output ========
		for(k=0; k<ZGridN; k++)
		{
			for(j=1; j<=YGridN; j++)
			{
				Cell &ccc = p_Mesh()->GetCell(i,j,k);
				 ne[j-1][k] = 0.5*(ccc.W_Jxx + ccc.W_Jyy + ccc.W_Denn)+0.5*ccc.W_Denn*(1.0+0.5*ccc.W_Asq)/(1+ccc.W_Psi)/(1+ccc.W_Psi);
				 nn[j-1][k] = ccc.W_Denn;
				 nb[j-1][k] = ccc.B_Den;
				Psi[j-1][k] = ccc.W_Psi;

				WEx[j-1][k] = -ccc.W_Ex+ccc.W_By;
				WEy[j-1][k] = -ccc.W_Ey-ccc.W_Bx;
				WEz[j-1][k] = ccc.W_Ez;

				WBx[j-1][k] = ccc.W_Bx;
				WBy[j-1][k] = ccc.W_By;
				WBz[j-1][k] = ccc.W_Bz;

			}
		}

	break;

	}


	//========================================
	//========= Create nc File ===============
	sprintf(sFile,"Fields2d_%d.nc",nt);
	if((retval = ncmpi_create(MPI_COMM_WORLD, sFile, NC_CLOBBER|NC_64BIT_OFFSET, info, &ncid)) ) 
	{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while creating " << sFile << "." <<'\n';
    	return NC_ERR;
	}

	//========================================
	//========= Define Dimension =============
	if ( (retval = ncmpi_def_dim(ncid, "nz", ZGridN, 	 &nz_id)) ) return NC_ERR;

	switch(savedim)
	{
		case 1:
		if ( (retval = ncmpi_def_dim(ncid, "nx", XGridN*Xpa, &nx_id)) )	return NC_ERR;
		F_DIM[0] = nx_id;
  		F_DIM[1] = nz_id;
		break;

		case 2:
		if ( (retval = ncmpi_def_dim(ncid, "ny", YGridN*Ypa, &ny_id)) ) return NC_ERR;
		F_DIM[0] = ny_id;
  		F_DIM[1] = nz_id;
		break;
	}

	
	//========================================
	//========= Define Variables =============
	if ( (retval = ncmpi_def_var(ncid, "n_e", NC_DOUBLE, 2, F_DIM, &n_e_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "n_n", NC_DOUBLE, 2, F_DIM, &n_n_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "Psi", NC_DOUBLE, 2, F_DIM, &Psi_id)) ) return NC_ERR;

	if ( (retval = ncmpi_def_var(ncid, "WEx", NC_DOUBLE, 2, F_DIM, &WEx_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WEy", NC_DOUBLE, 2, F_DIM, &WEy_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WEz", NC_DOUBLE, 2, F_DIM, &WEz_id)) ) return NC_ERR;

	if ( (retval = ncmpi_def_var(ncid, "WBx", NC_DOUBLE, 2, F_DIM, &WBx_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WBy", NC_DOUBLE, 2, F_DIM, &WBy_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_var(ncid, "WBz", NC_DOUBLE, 2, F_DIM, &WBz_id)) ) return NC_ERR;

	if(Nbeam){
	if ( (retval = ncmpi_def_var(ncid, "n_b", NC_DOUBLE, 2, F_DIM, &n_b_id)) ) return NC_ERR;}


	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{
		if(ifAx[NF])
		{
		sprintf(vname1,"AxR_%d",NF);sprintf(vname2,"AxI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &AxR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &AxI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LExR_%d",NF);sprintf(vname2,"LExI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &LExR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &LExI_id[NF])) ) return NC_ERR;

		sprintf(vname1,"LByR_%d",NF);sprintf(vname2,"LByI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &LByR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &LByI_id[NF])) ) return NC_ERR;

		}

		if(ifAy[NF])
		{
		sprintf(vname1,"AyR_%d",NF);sprintf(vname2,"AyI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &AyR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &AyI_id[NF])) ) return NC_ERR;
		
		sprintf(vname1,"LEyR_%d",NF);sprintf(vname2,"LEyI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &LEyR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &LEyI_id[NF])) ) return NC_ERR;
		
		sprintf(vname1,"LBxR_%d",NF);sprintf(vname2,"LBxI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &LBxR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &LBxI_id[NF])) ) return NC_ERR;
		}

		sprintf(vname1,"LEzR_%d",NF);sprintf(vname2,"LEzI_%d",NF);
		if ( (retval = ncmpi_def_var(ncid, vname1, NC_DOUBLE, 2, F_DIM, &LEzR_id[NF])) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, vname2, NC_DOUBLE, 2, F_DIM, &LEzI_id[NF])) ) return NC_ERR;
	}



	if ( (retval = ncmpi_enddef(ncid)) ) return NC_ERR;

	//========================================
	//============ Put Variables =============


	if ( (retval = ncmpi_put_vara_double_all(ncid, n_e_id, fstart, fcount,  &ne[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, n_n_id, fstart, fcount,  &nn[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, Psi_id, fstart, fcount, &Psi[0][0])) ) return NC_ERR;

  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEx_id, fstart, fcount, &WEx[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEy_id, fstart, fcount, &WEy[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WEz_id, fstart, fcount, &WEz[0][0])) ) return NC_ERR;

  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBx_id, fstart, fcount, &WBx[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBy_id, fstart, fcount, &WBy[0][0])) ) return NC_ERR;
  	if ( (retval = ncmpi_put_vara_double_all(ncid, WBz_id, fstart, fcount, &WBz[0][0])) ) return NC_ERR;
	
	if(Nbeam){
	if ( (retval = ncmpi_put_vara_double_all(ncid, n_b_id, fstart, fcount,  &nb[0][0])) ) return NC_ERR;}


	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{

		switch(savedim)
		{
			case 1:
			j = int(YGridN*(Ypa/2.0+1-IndexY))+1;
			for(k=0; k<ZGridN; k++)
			{
				for(i=1; i<=XGridN; i++)
				{
					Cell &ccc = p_Mesh()->GetCell(i,j,k);
					if(ifAx[NF])
					{
						AxR[i-1][k] = ccc.Acomx[NF].real();
						AxI[i-1][k] = ccc.Acomx[NF].imag();

						LExR[i-1][k] = ccc.L_Ex[NF].real();
						LExI[i-1][k] = ccc.L_Ex[NF].imag();

						LByR[i-1][k] = ccc.L_By[NF].real();
						LByI[i-1][k] = ccc.L_By[NF].imag();

					}
					if(ifAy[NF])
					{
						AyR[i-1][k] = ccc.Acomy[NF].real();
						AyI[i-1][k] = ccc.Acomy[NF].imag();

						LEyR[i-1][k] = ccc.L_Ey[NF].real();
						LEyI[i-1][k] = ccc.L_Ey[NF].imag();

						LBxR[i-1][k] = ccc.L_Bx[NF].real();
						LBxI[i-1][k] = ccc.L_Bx[NF].imag();
					}

						LEzR[i-1][k] = ccc.L_Ez[NF].real();
						LEzI[i-1][k] = ccc.L_Ez[NF].imag();
				}
			}
			break;
			case 2:
			i = int(XGridN*(Xpa/2.0+1-IndexX))+1;
			for(k=0; k<ZGridN; k++)
			{
				for(j=1; j<=YGridN; j++)
				{
					Cell &ccc = p_Mesh()->GetCell(i,j,k);
					if(ifAx[NF])
					{
						AxR[j-1][k] = ccc.Acomx[NF].real();
						AxI[j-1][k] = ccc.Acomx[NF].imag();

						LExR[j-1][k] = ccc.L_Ex[NF].real();
						LExI[j-1][k] = ccc.L_Ex[NF].imag();

						LByR[j-1][k] = ccc.L_By[NF].real();
						LByI[j-1][k] = ccc.L_By[NF].imag();
					}
					if(ifAy[NF])
					{
						AyR[j-1][k] = ccc.Acomy[NF].real();
						AyI[j-1][k] = ccc.Acomy[NF].imag();

						LEyR[j-1][k] = ccc.L_Ey[NF].real();
						LEyI[j-1][k] = ccc.L_Ey[NF].imag();

						LBxR[j-1][k] = ccc.L_Bx[NF].real();
						LBxI[j-1][k] = ccc.L_Bx[NF].imag();
					}
						LEzR[j-1][k] = ccc.L_Ez[NF].real();
						LEzI[j-1][k] = ccc.L_Ez[NF].imag();
				}
			}
			break;
		}

		if(ifAx[NF]){
			if ( (retval = ncmpi_put_vara_double_all(ncid, AxR_id[NF],  fstart, fcount,  &AxR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, AxI_id[NF],  fstart, fcount,  &AxI[0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LExR_id[NF], fstart, fcount, &LExR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LExI_id[NF], fstart, fcount, &LExI[0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LByR_id[NF], fstart, fcount, &LByR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LByI_id[NF], fstart, fcount, &LByI[0][0])) ) return NC_ERR;}


		if(ifAy[NF]){
			if ( (retval = ncmpi_put_vara_double_all(ncid, AyR_id[NF],  fstart, fcount,  &AyR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, AyI_id[NF],  fstart, fcount,  &AyI[0][0])) ) return NC_ERR;


			if ( (retval = ncmpi_put_vara_double_all(ncid, LEyR_id[NF], fstart, fcount, &LEyR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LEyI_id[NF], fstart, fcount, &LEyI[0][0])) ) return NC_ERR;

			if ( (retval = ncmpi_put_vara_double_all(ncid, LBxR_id[NF], fstart, fcount, &LBxR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LBxI_id[NF], fstart, fcount, &LBxI[0][0])) ) return NC_ERR;}
			if ( (retval = ncmpi_put_vara_double_all(ncid, LEzR_id[NF], fstart, fcount, &LEzR[0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_put_vara_double_all(ncid, LEzI_id[NF], fstart, fcount, &LEzI[0][0])) ) return NC_ERR;
	
	}
	
	if ( (retval = ncmpi_sync(ncid)) )
  	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval <<  " while syncing " << sFile << "." <<'\n';
		return NC_ERR; 
	}
  
	if ( (retval = ncmpi_close(ncid)) ) 
	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
		return NC_ERR;
	}


	
	if(Nbeam) SaveP(nt);

	SaveXray(nt);

	return 0;

}







int Domain::SaveP(int nt)
{
	int N_processor;
	MPI_Comm_size(MPI_COMM_WORLD, &N_processor);

	char sFile[128];
	char pFile[128];
	static const int NC_ERR = 1;

	int Xpa = p_Partition()->GetXpart();
	int Ypa = p_Partition()->GetYpart();

	int Idx_x = p_Partition()->RankIdx_X();
	int Idx_y = p_Partition()->RankIdx_Y();

	int n, type, Npart, NAll, Nparts[N_processor];

	//========================================
	//========= Pnetcdf Variables ============
	int retval, ncid; 
	int np_id, F_DIM[1];
	int x_id, y_id, z_id, x0_id, y0_id, z0_id;
	int px_id, py_id, pz_id, Ex0_id, Ey0_id, Ez0_id;
	int q2m_id, wei_id;

	int Wxw_id, Wyw_id, Wzw_id;
	int Wxl_id, Wyl_id, Wzl_id;

	int PP_id; //particle partition information;
	int PDim_id, PSC_id, P_DIM[2]; //start and count//dim
	MPI_Offset pstart[2], pcount[2];


	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Offset fstart[1], fcount[1];
	int SC[1][2];




	for (n=0; n<NSpecie; n++)
	{


		type = SpecieType[n];

		Npart= Get_NSpecie(type);


		double   *x_P, *y_P, *z_P, *x0_P, *y0_P,  *z0_P;	
		double   *px_P,  *py_P,  *pz_P;	
		double   *Ex0_P, *Ey0_P, *Ez0_P;
		double   *Wxw_P, *Wyw_P, *Wzw_P;
		double   *Wxl_P, *Wyl_P, *Wzl_P;
		double   *q2m_P, *Weight_P;

		x_P  = new double[Npart+1];
		y_P  = new double[Npart+1];  
		z_P  = new double[Npart+1];
		x0_P = new double[Npart+1];  
		y0_P = new double[Npart+1];  
		z0_P = new double[Npart+1];
		px_P = new double[Npart+1];
		py_P = new double[Npart+1];
		pz_P = new double[Npart+1];
		Ex0_P= new double[Npart+1];
		Ey0_P= new double[Npart+1];
		Ez0_P= new double[Npart+1];

		Wxw_P = new double[Npart+1];
		Wyw_P = new double[Npart+1];
		Wzw_P = new double[Npart+1];
		Wxl_P = new double[Npart+1];
		Wyl_P = new double[Npart+1];
		Wzl_P = new double[Npart+1];

		q2m_P    = new double[Npart+1];
		Weight_P = new double[Npart+1];




		Particle *p = p_Meshes -> p_Particle;
		int npp=0;
		while(p)
		{
			if(p->type==type)
			{
				x_P[npp]  = p->x;
				y_P[npp]  = p->y;
				z_P[npp]  = p->z;

				x0_P[npp] = p->x0;
				y0_P[npp] = p->y0;
				z0_P[npp] = p->z0;

				px_P[npp] = p->px;
				py_P[npp] = p->py;
				pz_P[npp] = p->pz;

				Ex0_P[npp]= p->Ex0;
				Ey0_P[npp]= p->Ey0;
				Ez0_P[npp]= p->Ez0;

				Wxw_P[npp]= p->Wxw;
				Wyw_P[npp]= p->Wyw;
				Wzw_P[npp]= p->Wzw;

				Wxl_P[npp]= p->Wxl;
				Wyl_P[npp]= p->Wyl;
				Wzl_P[npp]= p->Wzl;

				q2m_P[npp]= p->q2m;
				Weight_P[npp]=p->weight;

				npp++;
			}

			p = p-> p_PrevPart;

		}


		MPI_Allgather(&Npart, 1, MPI_INT, &Nparts, 1, MPI_INT, MPI_COMM_WORLD);

		fstart[0] = 0;
		for(int i=0; i<Rank; i++) {fstart[0]+= Nparts[i];};
		fcount[0] = Npart;

		NAll = 0;
		for(int i=0; i<N_processor; i++) {NAll += Nparts[i];};

		// //save partition for the particles;
		// sprintf(pFile,"PPartition_Rank_%d_sp_%d.dt",rank,n);
		// FILE * dFile;
		// dFile = fopen (pFile,"w");
		// fprintf(dFile, "%d %d\n", fstart[0],fcount[0]);
		// fclose (dFile);


		sprintf(sFile,"Parts_%d_sp_%d.nc",nt,n);

		if((retval = ncmpi_create(MPI_COMM_WORLD, sFile, NC_CLOBBER|NC_64BIT_OFFSET, info, &ncid)) ) 
		{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while creating " << sFile << "." <<'\n';
    	return NC_ERR;
		}



		//========================================
		//========= Define Dimension =============
		if ( (retval = ncmpi_def_dim(ncid, "np", NAll, &np_id)) )	return NC_ERR;
		F_DIM[0] = np_id;

		//particle partition information;
		if ( (retval = ncmpi_def_dim(ncid, "npro", N_processor, &PDim_id)) )	return NC_ERR;
		if ( (retval = ncmpi_def_dim(ncid, "SC", 	2, 			&PSC_id)) )		return NC_ERR;

		P_DIM[0] = PDim_id; //Number of processors
		P_DIM[1] = PSC_id;	//start and count

		pstart[0]=Rank;
		pstart[1]=0;
		pcount[0]=1;
		pcount[1]=2;

		SC[0][0]=fstart[0];
  		SC[0][1]=fcount[0];

		

		//========================================
		//========= Define Variables =============
		if ( (retval = ncmpi_def_var(ncid, "partition", NC_INT, 2, P_DIM, &PP_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "x0", NC_DOUBLE, 1, F_DIM, &x0_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "y0", NC_DOUBLE, 1, F_DIM, &y0_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "z0", NC_DOUBLE, 1, F_DIM, &z0_id)) ) return NC_ERR;

		if ( (retval = ncmpi_def_var(ncid, "xx", NC_DOUBLE, 1, F_DIM, &x_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "yy", NC_DOUBLE, 1, F_DIM, &y_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "zz", NC_DOUBLE, 1, F_DIM, &z_id)) ) return NC_ERR;

		if ( (retval = ncmpi_def_var(ncid, "px", NC_DOUBLE, 1, F_DIM, &px_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "py", NC_DOUBLE, 1, F_DIM, &py_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "pz", NC_DOUBLE, 1, F_DIM, &pz_id)) ) return NC_ERR;

		if ( (retval = ncmpi_def_var(ncid, "Ex", NC_DOUBLE, 1, F_DIM, &Ex0_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Ey", NC_DOUBLE, 1, F_DIM, &Ey0_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Ez", NC_DOUBLE, 1, F_DIM, &Ez0_id)) ) return NC_ERR;
		
		if ( (retval = ncmpi_def_var(ncid, "Wxw", NC_DOUBLE, 1, F_DIM, &Wxw_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Wyw", NC_DOUBLE, 1, F_DIM, &Wyw_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Wzw", NC_DOUBLE, 1, F_DIM, &Wzw_id)) ) return NC_ERR;

		if ( (retval = ncmpi_def_var(ncid, "Wxl", NC_DOUBLE, 1, F_DIM, &Wxl_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Wyl", NC_DOUBLE, 1, F_DIM, &Wyl_id)) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "Wzl", NC_DOUBLE, 1, F_DIM, &Wzl_id)) ) return NC_ERR;

		if ( (retval = ncmpi_def_var(ncid, "q2m",    NC_DOUBLE, 1, F_DIM, &q2m_id )) ) return NC_ERR;
		if ( (retval = ncmpi_def_var(ncid, "weight", NC_DOUBLE, 1, F_DIM, &wei_id )) ) return NC_ERR;

		if ( (retval = ncmpi_enddef(ncid)) ) return NC_ERR;



		//========================================
		//============ Put Variables =============
		if ( (retval = ncmpi_put_vara_double_all(ncid, x0_id, fstart, fcount,  &x0_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, y0_id, fstart, fcount,  &y0_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, z0_id, fstart, fcount,  &z0_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, x_id,  fstart, fcount,   &x_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, y_id,  fstart, fcount,   &y_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, z_id,  fstart, fcount,   &z_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, px_id, fstart, fcount,  &px_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, py_id, fstart, fcount,  &py_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, pz_id, fstart, fcount,  &pz_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, Ex0_id,fstart, fcount, &Ex0_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Ey0_id,fstart, fcount, &Ey0_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Ez0_id,fstart, fcount, &Ez0_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wxw_id,fstart, fcount, &Wxw_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wyw_id,fstart, fcount, &Wyw_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wzw_id,fstart, fcount, &Wzw_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wxl_id,fstart, fcount, &Wxl_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wyl_id,fstart, fcount, &Wyl_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, Wzl_id,fstart, fcount, &Wzl_P[0]))) return NC_ERR;

  		if ( (retval = ncmpi_put_vara_double_all(ncid, q2m_id,fstart, fcount, &q2m_P[0]))) return NC_ERR;
  		if ( (retval = ncmpi_put_vara_double_all(ncid, wei_id,fstart, fcount, &Weight_P[0]))) return NC_ERR;
  		
		if ( (retval = ncmpi_put_vara_int_all(ncid, PP_id, pstart, pcount,  &SC[0][0]))) return NC_ERR;




		if ( (retval = ncmpi_sync(ncid)) )
  		{
			if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval <<  " while syncing " << sFile << "." <<'\n';
			return NC_ERR; 
		}
  
		if ( (retval = ncmpi_close(ncid)) ) 
		{
			if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
			return NC_ERR;
		}



		delete [] x_P; 	 delete [] y_P;   delete [] z_P; 
		delete [] x0_P;  delete [] y0_P;  delete [] z0_P;	
		delete [] px_P;  delete [] py_P;  delete [] pz_P;	
		delete [] Ex0_P; delete [] Ey0_P; delete [] Ez0_P;
		delete [] Wxw_P; delete [] Wyw_P; delete [] Wzw_P;
		delete [] Wxl_P; delete [] Wyl_P; delete [] Wzl_P;
		delete [] q2m_P; delete [] Weight_P;




	}


	return 0;

}



int Domain::SaveXray(int nt)
{

	int type;
	int Npart=0;
	int Nparts;

	if(p_Meshes->XRayDetector->IfRadiation==0)   return 0;

	for (int n=0; n<NSpecie; n++)
	{
		type   = SpecieType[n];
		Npart += Get_NSpecie(type);
	}
	MPI_Allreduce(&Npart, &Nparts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	
 	if( Nparts==0) return 0;
 	
	int NOmega = p_Meshes->XRayDetector->NOmega;
	int NTheta = p_Meshes->XRayDetector->NTheta;
	int NPhi   = p_Meshes->XRayDetector->NPhi;

	double Xray[NOmega][NTheta][NPhi];
	double XraySpatial[NTheta][NPhi];
	double XraySpatialAll[NTheta][NPhi];

	for(int i=0;i<NOmega;i++)
	{

		for(int j=0;j<NTheta;j++)
		{
			for(int k=0;k<NPhi;k++)
			{
				XraySpatial[j][k]   =p_Meshes->XRayDetector->GetDetector(i,j,k);
				XraySpatialAll[j][k]=0.0;
			}
		}
		MPI_Reduce(&XraySpatial[0][0], &XraySpatialAll[0][0],NTheta*NPhi, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(Rank==0) 
		{


			for(int j=0;j<NTheta;j++)
			{
				for(int k=0;k<NPhi;k++)
				{
					Xray[i][j][k]=XraySpatialAll[j][k];

				}
			}


		}

	}



	
	char sFile[128];
	static const int NC_ERR = 1;
	int retval, ncid; 
	int no_id,  nt_id,  np_id, F_DIM[3];
	int Xray_id;
	MPI_Info info;
	MPI_Info_create(&info);
	MPI_Offset fstart[3], fcount[3];

	//========================================
	//========= Create nc File ===============
	sprintf(sFile,"XRay_%d.nc",nt);
	if((retval = ncmpi_create(MPI_COMM_WORLD, sFile, NC_CLOBBER|NC_64BIT_OFFSET, info, &ncid)) ) 
	{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while creating " << sFile << "." <<'\n';
    	return NC_ERR;
	}

	if ( (retval = ncmpi_def_dim(ncid, "n_omega",    NOmega, &no_id)) )	return NC_ERR;
	if ( (retval = ncmpi_def_dim(ncid, "n_theta(x)", NTheta, &nt_id)) ) return NC_ERR;
	if ( (retval = ncmpi_def_dim(ncid, "n_theta(y)", NPhi, 	 &np_id)) ) return NC_ERR;


	F_DIM[0] = no_id;
  	F_DIM[1] = nt_id;
  	F_DIM[2] = np_id;

  	fstart[0]=0;fstart[1]=0;fstart[2]=0;
  	fcount[0]=0;fcount[1]=0;fcount[2]=0;

  	if(Rank==0) {fcount[0]=NOmega;fcount[1]=NTheta;fcount[2]=NPhi;}
  	//========================================
	//========= Define Variables =============
	if ( (retval = ncmpi_def_var(ncid, "X_Ray", NC_DOUBLE, 3, F_DIM, &Xray_id)) ) return NC_ERR;
	if ( (retval = ncmpi_enddef(ncid)) ) return NC_ERR;

	if ( (retval = ncmpi_put_vara_double_all(ncid, Xray_id, fstart, fcount,  &Xray[0][0][0])) ) return NC_ERR;

	if ( (retval = ncmpi_sync(ncid)) )
  	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval <<  " while syncing " << sFile << "." <<'\n';
		return NC_ERR; 
	}
  
	if ( (retval = ncmpi_close(ncid)) ) 
	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
		return NC_ERR;
	}

	

	



	return 0;


}




int Domain::LoadPulse(int nt)
{	
	//---set pulse amplitude zeros-------
	for (int k=0; k<ZGridN; k++) 
	{
		for (int j=0; j<YGridN+2; j++) 
		{
			for (int i=0; i<XGridN+2; i++)
			{
				Cell &c = p_Mesh()->GetCell(i, j, k);
				for (int f=0;f<NFreqs;f++)
				{
					c.Acomx[f]  = c.Acomy[f]  = c.Acomxm[f] = c.Acomym[f] = 0.0;
				}
               
			}
		}
	}


	//-------pnetcdf part--------------

	char sFile[128];
	char vname1[24];
	char vname2[24];
	static const int NC_ERR = 1;

	int Xpa = p_Partition()->GetXpart();
	int Ypa = p_Partition()->GetYpart();

	int Idx_x = p_Partition()->RankIdx_X();
	int Idx_y = p_Partition()->RankIdx_Y();

	int i,j,k, NF;

	//========================================
	//========= Pnetcdf Variables ============
	int retval, ncid; 
	int AxR_id[NFreqs],		AxI_id[NFreqs],		AyR_id[NFreqs],		AyI_id[NFreqs];
	int AxR_natts[NFreqs],	AxI_natts[NFreqs],	AyR_natts[NFreqs],	AyI_natts[NFreqs];


	MPI_Offset fstart[3], fcount[3];

	// pulse info
	int nxx = XGridN+ (Idx_x>1) + (Idx_x<Xpa);
	int nyy = YGridN+ (Idx_y>1) + (Idx_y<Ypa);

	int offx = -(Idx_x>1);
	int offy = -(Idx_y>1);


	double AxR[nxx][nyy][ZGridN], AxI[nxx][nyy][ZGridN];
	double AyR[nxx][nyy][ZGridN], AyI[nxx][nyy][ZGridN];



	// reading start and count
	fstart[0] = (Idx_x-1)*XGridN+offx;
	fstart[1] = (Idx_y-1)*YGridN+offy;
	fstart[2] = 0;

	fcount[0] = nxx;
	fcount[1] = nyy;
	fcount[2] = ZGridN;


	//========================================
	//========= open existing nc File ===============
	sprintf(sFile,"Fields_%d.nc",nt);
	if((retval = ncmpi_open(MPI_COMM_WORLD, sFile, NC_NOWRITE, MPI_INFO_NULL, &ncid)) ) 
	{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while opening " << sFile << "." <<'\n';
    	return NC_ERR;
	}

	//========================================


	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{
		if(ifAx[NF])
		{
		sprintf(vname1,"AxR_%d",NF);sprintf(vname2,"AxI_%d",NF);
			if ( (retval = ncmpi_inq_varid(ncid, vname1 , &AxR_id[NF]))) return NC_ERR;
			if ( (retval = ncmpi_inq_varid(ncid, vname2 , &AxI_id[NF]))) return NC_ERR;
		}

		if(ifAy[NF])
		{
		sprintf(vname1,"AyR_%d",NF);sprintf(vname2,"AyI_%d",NF);
			if ( (retval = ncmpi_inq_varid(ncid, vname1 , &AyR_id[NF]))) return NC_ERR;
			if ( (retval = ncmpi_inq_varid(ncid, vname2 , &AyI_id[NF]))) return NC_ERR;
		}

	}

	//===== Laser Wake potential =======
	for (NF=0; NF<NFreqs; NF++)
	{

		for(k=0; k<ZGridN; k++)
		{
			for(j=0; j<nyy; j++)
			{
				for(i=0; i<nxx; i++)
				{
					AxR[i][j][k]=AxI[i][j][k]=AyR[i][j][k]=AyI[i][j][k]=0.0;
				}
			}
		}
		
		if(ifAx[NF])
		{ 
			if ( (retval = ncmpi_get_vara_double_all(ncid, AxR_id[NF], fstart, fcount,   &AxR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, AxI_id[NF], fstart, fcount,   &AxI[0][0][0])) ) return NC_ERR;
	
		}
	
				
		if(ifAy[NF])
		{
			if ( (retval = ncmpi_get_vara_double_all(ncid, AyR_id[NF], fstart, fcount,   &AyR[0][0][0])) ) return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, AyI_id[NF], fstart, fcount,   &AyI[0][0][0])) ) return NC_ERR;
		}


		// put vector potential on the cell.
		for(k=0; k<ZGridN; k++)
		{
			for(j=1; j<=nyy; j++)
			{
				for(i=1; i<=nxx; i++)
				{
					Cell &ccc = p_Mesh()->GetCell(i+offx,j+offy,k);
					ccc.AddAComs(AxR[i-1][j-1][k]+ci*AxI[i-1][j-1][k], AyR[i-1][j-1][k]+ci*AyI[i-1][j-1][k], NF);
				}
			}
		}

	}

	if ( (retval = ncmpi_close(ncid)) ) 
	{
		if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
		return NC_ERR;
	}

	return 0;

}

int  Domain::LoadParti(int nt)
{

	int N_processor;
	MPI_Comm_size(MPI_COMM_WORLD, &N_processor);

	char sFile[128];
	static const int NC_ERR = 1;

	double Offset_X = p_Mesh()->GetOffset_X();
	double Offset_Y = p_Mesh()->GetOffset_Y();


	double XXmax = Offset_X+XGridN*dx;
	double YYmax = Offset_Y+YGridN*dy;
	double ZZmax = (ZGridN-2)*dz;

	double XXmin = Offset_X;
	double YYmin = Offset_Y;
	double ZZmin = 0.0;


	int n, type;

	//========================================
	//========= Pnetcdf Variables ============
	int retval, ncid; 
	int np_id, F_DIM[1];

	int x_id, y_id, z_id, x0_id, y0_id, z0_id;
	int px_id, py_id, pz_id, Ex0_id, Ey0_id, Ez0_id;
	int q2m_id, wei_id;
	int Wxw_id, Wyw_id, Wzw_id;
	int Wxl_id, Wyl_id, Wzl_id;

	//=======
	int PP_id;
	int SCAll[N_processor][2];
	MPI_Offset pstart[2], pcount[2];
	pstart[0]=pstart[1]=0;
	pcount[0]=N_processor; pcount[1]=2;


	MPI_Offset fstart[1], fcount[1];
	int Npart;
	Particle *p =NULL;

	for (n=0; n<NSpecie; n++)
	{

		type = SpecieType[n];

		sprintf(sFile,"Parts_%d_sp_%d.nc",nt,n);
		if((retval = ncmpi_open(MPI_COMM_WORLD, sFile, NC_NOWRITE, MPI_INFO_NULL, &ncid)) ) 
		{
    	if(Rank==0) std::cout<< "Domain: pnetcdf error:" << retval << " while opening " << sFile << "." <<'\n';
    	return NC_ERR;
		}

		if ( (retval = ncmpi_inq_varid(ncid,  "x0" , 	 &x0_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "y0" , 	 &y0_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "z0" , 	 &z0_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid,  "xx" , 	  &x_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "yy" , 	  &y_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "zz" , 	  &z_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid,  "px" , 	 &px_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "py" , 	 &py_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "pz" , 	 &pz_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid,  "Ex" , 	&Ex0_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "Ey" , 	&Ey0_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid,  "Ez" , 	&Ez0_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid, "Wxw" , 	&Wxw_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid, "Wyw" , 	&Wyw_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid, "Wzw" , 	&Wzw_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid, "Wxl" , 	&Wxl_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid, "Wyl" , 	&Wyl_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid, "Wzl" , 	&Wzl_id))) return NC_ERR;

		if ( (retval = ncmpi_inq_varid(ncid, "q2m" , 	&q2m_id))) return NC_ERR;
		if ( (retval = ncmpi_inq_varid(ncid, "weight" , &wei_id))) return NC_ERR;

		// read particle partitioning info

		if ( (retval = ncmpi_inq_varid(ncid, "partition" , &PP_id))) return NC_ERR;
		if ( (retval = ncmpi_get_vara_int_all(ncid, PP_id, pstart, pcount,  &SCAll[0][0])) ) return NC_ERR;
		// number of particles in each processor.
		Npart=SCAll[Rank][1];

			double   *x_P, *y_P, *z_P, *x0_P, *y0_P,  *z0_P;	
			double   *px_P,  *py_P,  *pz_P;	
			double   *Ex0_P, *Ey0_P, *Ez0_P;
			double   *Wxw_P, *Wyw_P, *Wzw_P;
			double   *Wxl_P, *Wyl_P, *Wzl_P;
			double   *Weight_P;

			x_P  = new double[Npart];
			y_P  = new double[Npart];  
			z_P  = new double[Npart];
			x0_P = new double[Npart];  
			y0_P = new double[Npart];  
			z0_P = new double[Npart];
			px_P = new double[Npart];
			py_P = new double[Npart];
			pz_P = new double[Npart];
			Ex0_P= new double[Npart];
			Ey0_P= new double[Npart];
			Ez0_P= new double[Npart];

			Wxw_P = new double[Npart];
			Wyw_P = new double[Npart];
			Wzw_P = new double[Npart];
			Wxl_P = new double[Npart];
			Wyl_P = new double[Npart];
			Wzl_P = new double[Npart];
			Weight_P = new double[Npart];

			fstart[0]=SCAll[Rank][0];
			fcount[0]=SCAll[Rank][1];

			if ( (retval = ncmpi_get_vara_double_all(ncid, x_id, fstart, fcount,  &x_P[0])))   return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, y_id, fstart, fcount,  &y_P[0])))   return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, z_id, fstart, fcount,  &z_P[0])))   return NC_ERR;

			if ( (retval = ncmpi_get_vara_double_all(ncid, x0_id, fstart, fcount, &x0_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, y0_id, fstart, fcount, &y0_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, z0_id, fstart, fcount, &z0_P[0])))  return NC_ERR;

			if ( (retval = ncmpi_get_vara_double_all(ncid, px_id, fstart, fcount, &px_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, py_id, fstart, fcount, &py_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, pz_id, fstart, fcount, &pz_P[0])))  return NC_ERR;

			if ( (retval = ncmpi_get_vara_double_all(ncid, Ex0_id, fstart, fcount, &Ex0_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Ey0_id, fstart, fcount, &Ey0_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Ez0_id, fstart, fcount, &Ez0_P[0])))  return NC_ERR;


			if ( (retval = ncmpi_get_vara_double_all(ncid, Wxw_id, fstart, fcount, &Wxw_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Wyw_id, fstart, fcount, &Wyw_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Wzw_id, fstart, fcount, &Wzw_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Wxl_id, fstart, fcount, &Wxl_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Wyl_id, fstart, fcount, &Wyl_P[0])))  return NC_ERR;
			if ( (retval = ncmpi_get_vara_double_all(ncid, Wzl_id, fstart, fcount, &Wzl_P[0])))  return NC_ERR;

			if ( (retval = ncmpi_get_vara_double_all(ncid, wei_id, fstart, fcount, &Weight_P[0])))  return NC_ERR;



			for(int i=0;i<Npart;i++)
			{	

				switch(type)
				{

					case ELECTRON:
					p = new Electron(x0_P[i], y0_P[i], z0_P[i], px_P[i], py_P[i], pz_P[i], Ex0_P[i], Ey0_P[i], Ez0_P[i], 1.0, Weight_P[i]);
					p->x=x_P[i];
					p->y=y_P[i];
					p->z=z_P[i];

					p->Wxw = Wxw_P[i];
					p->Wyw = Wyw_P[i];
					p->Wzw = Wzw_P[i];

					p->Wxl = Wxl_P[i];
					p->Wyl = Wyl_P[i];
					p->Wzl = Wzl_P[i];
					break;
				}
			}


		if ( (retval = ncmpi_close(ncid)) ) 
		{
			if(Rank==0) std::cout<< "Domain: pnetcdf error " << retval << " while closing " << sFile << "." <<'\n';
			return NC_ERR;
		}


	}
		





	return 0;

}







