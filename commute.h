//----------------------------------------------------------------------------------||
//-------------------                commute.h                   -------------------||
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



#ifndef H_COMM
#define H_COMM

#include<vector>

//---------------------------- Mesh class -----------------------
class Commute
   {
   friend class Domain;
   friend class Mesh;

private:
   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.


   int bufsize;

   WDOUBLE *SendSourceXm;  //Array to send for the sources;
   WDOUBLE *SendSourceYm;  //Array to send for the sources;
   WDOUBLE *SendSourceXp;  //Array to send for the sources;
   WDOUBLE *SendSourceYp;  //Array to send for the sources;

   WDOUBLE *ReceSourceXm;  //Array for Rece 
   WDOUBLE *ReceSourceYm;  //Array for Rece 
   WDOUBLE *ReceSourceXp;  //Array for Rece 
   WDOUBLE *ReceSourceYp;  //Array for Rece 



   // diagonal
   WDOUBLE *SendSourcemm;  //Array to send for the sources;
   WDOUBLE *SendSourcemp;  //Array to send for the sources;
   WDOUBLE *SendSourcepm;  //Array to send for the sources;
   WDOUBLE *SendSourcepp;  //Array to send for the sources;

   WDOUBLE *ReceSourcemm;  //Array for Rece 
   WDOUBLE *ReceSourcemp;  //Array for Rece 
   WDOUBLE *ReceSourcepm;  //Array for Rece 
   WDOUBLE *ReceSourcepp;  //Array for Rece 

   //WDOUBLE *SendFields;  //Array to send for the fields;

   int SendSouSizeX;
   int SendSouSizeY;
   int SendFieSizeX;

   int Rank;

   int RankIdx_X;
   int RankIdx_Y;

   int Xpa;
   int Ypa;

   int XmPE; 
   int XpPE;
   int YmPE;
   int YpPE; 

   int mmPE; 
   int mpPE;
   int pmPE;
   int ppPE; 

   int GridX;
   int GridY;

   int kold;

   std::vector<WDOUBLE> CellAccX; 
   std::vector<WDOUBLE> CellAccY; 

public:

   
   void DoCommuteT(exchange what, std::vector<int> SendN);
   void    UnPackT(exchange what, std::vector<int> ReceN);

   void DoCommute(exchange what, int k);
   void    DoPack(exchange what, int k);
   void    UnPack(exchange what, int k);

   int Get_bufsize() {return bufsize;};



   Commute(int XGridN, int YGridN);
   ~Commute();


};

#endif
