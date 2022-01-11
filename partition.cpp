//----------------------------------------------------------------------------------||
//-------------------                partition.cpp               -------------------||
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


//---------------------------- Partition::Partition -----------------------
Partition::Partition (FILE *f) : NList ("Partition")
{
  
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &N_Processor);

  AddEntry((char*)"Xpartition", &N_Xpart, 1);
  AddEntry((char*)"Ypartition", &N_Ypart, 1);

  if (f)
    {
      rewind(f);
      read(f);
    }

  if (N_Processor != N_Xpart*N_Ypart)
    {
      if(Rank==0) std::cout<< "Error Partition Number. \n";
      exit(-10);
    }

// coordinates of the processor
  idx_Y= (Rank/N_Xpart)+1;
  idx_X= (Rank%N_Xpart)+1;  //tianhong May-18-2021

//Neighbor Processors: m==minus, p=plus.

  idx_XmPE=Rank-1;
  idx_XpPE=Rank+1;


 //tianhong May-18-2021
  if(idx_X==1)idx_XmPE+=N_Xpart; 
  if(idx_X==N_Xpart) idx_XpPE-=N_Xpart;

  idx_YmPE=Rank-N_Xpart;
  idx_YpPE=Rank+N_Xpart;

 //tianhong May-18-2021
  if(idx_Y==1) idx_YmPE+=N_Xpart*(N_Ypart);
  if(idx_Y==N_Ypart) idx_YpPE-=N_Xpart*(N_Ypart);


  idx_mmPE=idx_YmPE-1;
  idx_pmPE=idx_YmPE+1;

  idx_mpPE=idx_YpPE-1;
  idx_ppPE=idx_YpPE+1;

 //tianhong May-18-2021
  if(idx_X==1)
  {
    idx_mmPE+=N_Xpart;
    idx_mpPE+=N_Xpart; 
  }
  if(idx_X==N_Xpart)
  {
    idx_pmPE-=N_Xpart;
    idx_ppPE-=N_Xpart; 
  }
  
  // printf("%d: %d-%d-%d-%d-%d-%d-%d-%d-\n", Rank,idx_XmPE,idx_XpPE,idx_YmPE,idx_YpPE,idx_mmPE,idx_pmPE,idx_mpPE,idx_ppPE);
  

}