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

  //std::complex<double> Acom(1.0, 0.0);  
  //complex type of normalized vector potential

  double a0;

  double XSpotSize;
  double YSpotSize;
  double ZLength;

  double Xcenter;
  double Ycenter;
  double Zcenter;

  double ZFocal;

  double ZRise;
  double ZDrop;

  double Xphase;
  double Yphase;

  double Chirp;

  

  int Tprofile;
  int Lprofile;

  int i_Pulse;
  


//=====Laser Fields=======

  dcomplex Acomx;
  dcomplex Acomy;


public:

  double Xpol;
  double Ypol;
  double Omega;


  void Put_AComplex(double x, double y, double z, Cell &c, int NF);
  double    ProfileLongi(double x, double y, double z);
  double    ProfileTrans(double x, double y, double z);



//===


  Pulse(char *name, FILE *f);

  ~Pulse() { 
  };

};
#endif
