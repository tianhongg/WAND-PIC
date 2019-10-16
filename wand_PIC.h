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


#define dcomplex std::complex<double>
#define ci dcomplex(0.0,1.0)


#define FLD_DIM 6
#define WAK_DIM 9
#define BEA_DIM 5
#define SOU_DIM 6
#define SDT_DIM 11
#define SDP_DIM 21
#define WAK_DIM2 4


#define L_MPI

#define COMMU_S 11
#define COMMU_F 22
#define COMMU_A 33
#define COMMU_T 44
#define COMMU_P 55

#define COMMU_MG_P 66
#define COMMU_MG_S 77
#define COMMU_MG_R 88

#define COMMU_MG_P_C 99
#define COMMU_MG_S_C 100
#define COMMU_MG_R_C 110

#define MG_Phi 0
#define MG_Sou 1
#define MG_Res 2
#define MG_Chi 4




#include <mpi.h>
#include <complex>
#include <cmath>
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


#define Banner "\n||----------------------------------------------------------------------------------||\n||----------------------------------------------------------------------------------||\n||               __        ___    _   _ ____        ____ ___ ____                   ||\n||               \\ \\      / / \\  | \\ | |  _ \\      |  _ \\_ _/ ___|                  ||\n||                \\ \\ /\\ / / _ \\ |  \\| | | | |_____| |_) | | |                      ||\n||                 \\ V  V / ___ \\| |\\  | |_| |_____|  __/| | |___                   ||\n||                  \\_/\\_/_/   \\_\\_| \\_|____/      |_|  |___\\____|                  ||\n||                                                                                  ||\n||----------------------------------------------------------------------------------||\n||--  (W)akefield (A)cceleration a(n)d (D)LA - (P)article (i)n (C)ell Simulation  --||\n||----------------------------------------------------------------------------------||\n\n"


#endif


