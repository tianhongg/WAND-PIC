//----------------------------------------------------------------------------------||
//-------------------                wand_PIC.cpp                -------------------||
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


#include <stdio.h>
#include <stdlib.h>
#include "wand_PIC.h"


int main(int argc, char** argv)
{

#ifdef L_MPI
//====MPI=====initialization====
MPI_Init(&argc,&argv);
// Get the number of processes
int N_processor;
MPI_Comm_size(MPI_COMM_WORLD, &N_processor);
// Get the rank of the process
int Rank;
MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#endif

	if(Rank==0) 
	{	
		std::cout << Banner;
		std::cout << "\n=============================================\n";
		std::cout <<"==== Starting Program: Wand-PIC.         ====\n";
		printf("==== %8d Processors Initiated.      ====\n",N_processor);
	};

	Domain *domain = new Domain((char*)"WAND.ini",Rank);
	
	domain-> Run();
	

	delete domain;



#ifdef L_MPI
// Finalize the MPI environment.
MPI_Finalize();
#endif
return 0;
}