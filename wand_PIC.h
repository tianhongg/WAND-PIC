//----------------------------------------------------------------------------------||
//-------------------                wand_PIC.h                  -------------------||
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
//---Starting---------           : Jan-11-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||



#ifndef H_WAND
#define H_WAND


// switch the data type
#ifdef _WTYPE
	#define WDOUBLE double
	#define NC_WDOUBLE NC_DOUBLE
	#define MPI_WDOUBLE MPI_DOUBLE 
	#define ncmpi_put_att_WDOUBLE ncmpi_put_att_double
	#define ncmpi_get_att_WDOUBLE ncmpi_get_att_double
	#define ncmpi_put_vara_WDOUBLE_all ncmpi_put_vara_double_all
	#define ncmpi_get_vara_WDOUBLE_all ncmpi_get_vara_double_all
#else
	#define WDOUBLE float
	#define NC_WDOUBLE NC_FLOAT
	#define MPI_WDOUBLE MPI_FLOAT
	#define ncmpi_put_att_WDOUBLE ncmpi_put_att_float
	#define ncmpi_get_att_WDOUBLE ncmpi_get_att_float
	#define ncmpi_put_vara_WDOUBLE_all ncmpi_put_vara_float_all
	#define ncmpi_get_vara_WDOUBLE_all ncmpi_get_vara_float_all
#endif


#define dcomplex std::complex<WDOUBLE>
#define ci dcomplex(0.0,1.0)



#define L_MPI

#define FLD_DIM 6
#define WAK_DIM 9
#define BEA_DIM 5
#define SOU_DIM 6
#define SDT_DIM 13  // x-y-z-x0-y0-z0-vx-vy-old_x-old_y-old_vx-old_vy-sx-sy
#define SDP_DIM 23
#define WAK_DIM2 4

#define MG_Phi 0
#define MG_Sou 1
#define MG_Res 2
#define MG_Chi 4


// exchange
enum exchange
{
	COMMU_S,   //source
	COMMU_F,   // field
	COMMU_A,   // a0
	COMMU_T,   //trajectory
	COMMU_P,   //particle
	COMMU_SO,  //source

	COMMU_MG_P,  //potential
	COMMU_MG_S,  //source
	COMMU_MG_R,  //residual
	COMMU_MG_F, 
	COMMU_MG_C,  //\chi

	COMMU_MG_P_C,//potential
	COMMU_MG_S_C,//source
	COMMU_MG_R_C, //residual
	COMMU_MG_F_C, 
	COMMU_MG_C_C
};


#include <mpi.h>
#include <complex>
#include <cmath>
#include <time.h> 
#include <vector>
#include "namelist.h"
#include "domain.h"
#include "cell.h"
#include "pulse.h"
#include "partition.h"
#include "trajectory.h"
#include "particles.h"
#include "multigrid.h"
#include "mesh.h"
#include "commute.h"
#include <algorithm> 


#define Banner "\n||----------------------------------------------------------------------------------||\n||----------------------------------------------------------------------------------||\n||               __        ___    _   _ ____        ____ ___ ____                   ||\n||               \\ \\      / / \\  | \\ | |  _ \\      |  _ \\_ _/ ___|                  ||\n||                \\ \\ /\\ / / _ \\ |  \\| | | | |_____| |_) | | |                      ||\n||                 \\ V  V / ___ \\| |\\  | |_| |_____|  __/| | |___                   ||\n||                  \\_/\\_/_/   \\_\\_| \\_|____/      |_|  |___\\____|                  ||\n||                                                                                  ||\n||----------------------------------------------------------------------------------||\n||--  (W)akefield (A)cceleration a(n)d (D)LA - (P)article (i)n (C)ell Simulation  --||\n||----------------------------------------------------------------------------------||\n\n"


#endif


