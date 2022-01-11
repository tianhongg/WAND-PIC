//-------------------                pulse.h                     -------------------||
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
//---Starting---------           : Jan-22-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#ifndef H_PULSE
#define H_PULSE

#include <stdio.h>

//#define EX 1

//---------------------------- Pulse class -----------------------

class Pulse : public NList {
  friend class Domain;

private:

  Domain *p_domain() {return Domain::p_Domain;}; 

  //std::complex<WDOUBLE> Acom(1.0, 0.0);  
  //complex type of normalized vector potential

  WDOUBLE a0;

  WDOUBLE XSpotSize;
  WDOUBLE YSpotSize;
  WDOUBLE ZLength;

  WDOUBLE Xcenter;
  WDOUBLE Ycenter;
  WDOUBLE Zcenter;

  WDOUBLE ZFocal;

  WDOUBLE ZRise;
  WDOUBLE ZDrop;

  WDOUBLE Xphase;
  WDOUBLE Yphase;

  WDOUBLE Chirp;

  

  int Tprofile;
  int Lprofile;

  int i_Pulse;
  


//=====Laser Fields=======

  dcomplex Acomx;
  dcomplex Acomy;


public:

  WDOUBLE Xpol;
  WDOUBLE Ypol;
  WDOUBLE Omega;


  void Put_AComplex(WDOUBLE x, WDOUBLE y, WDOUBLE z, Cell &c, int NF);
  WDOUBLE    ProfileLongi(WDOUBLE x, WDOUBLE y, WDOUBLE z);
  WDOUBLE    ProfileTrans(WDOUBLE x, WDOUBLE y, WDOUBLE z);



//===


  Pulse(char *name, FILE *f);

  ~Pulse() { 
  };

};
#endif
