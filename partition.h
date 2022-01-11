//----------------------------------------------------------------------------------||
//-------------------                partition.h                 -------------------||
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


#ifndef H_PARTITION
#define H_PARTITION


//---------------------------- Partition class -----------------------

class Partition : public NList {

  friend class Domain;

private:

  Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.

  int N_Xpart, N_Ypart;   // partition number;
  int N_Processor;        // Total Number of Processor;
  int Rank;

  int idx_X,idx_Y;   // Index of processor in X and Y direction, start with 1;

  //int *p_myPE;      // pointer to my processor.
  int idx_XmPE;     // idx of processor left to me.
  int idx_XpPE;     // idx of processor right to me.
  int idx_YmPE;     // idx of processor below.
  int idx_YpPE;     // idx of processor above.

  int idx_mmPE;     // idx of processor left-blew
  int idx_pmPE;     // idx of processor right-blew
  int idx_mpPE;     // idx of processor left-above
  int idx_ppPE;     // idx of processor right-above




 public:

  int GetXpart() {return N_Xpart;};
  int GetYpart() {return N_Ypart;};

  int GetXmPE()  {return idx_XmPE;};
  int GetXpPE()  {return idx_XpPE;};
  int GetYmPE()  {return idx_YmPE;};
  int GetYpPE()  {return idx_YpPE;};

  int GetmmPE()  {return idx_mmPE;};
  int GetmpPE()  {return idx_mpPE;};
  int GetpmPE()  {return idx_pmPE;};
  int GetppPE()  {return idx_ppPE;};

  int RankIdx_X() {return idx_X;};
  int RankIdx_Y() {return idx_Y;};
  int  Get_Rank() {return Rank;};


  Partition(FILE *f);

  ~Partition() { 
  };
};

#endif
