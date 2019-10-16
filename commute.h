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



//---------------------------- Mesh class -----------------------
class Commute
   {
   friend class Domain;
   friend class Mesh;

private:
   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.


   int bufsize;

   double *SendSourceXm;  //Array to send for the sources;
   double *SendSourceYm;  //Array to send for the sources;
   double *SendSourceXp;  //Array to send for the sources;
   double *SendSourceYp;  //Array to send for the sources;

   double *ReceSourceXm;  //Array for Rece 
   double *ReceSourceYm;  //Array for Rece 
   double *ReceSourceXp;  //Array for Rece 
   double *ReceSourceYp;  //Array for Rece 

   //double *SendFields;  //Array to send for the fields;


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

   int GridX;
   int GridY;





public:


   
   void DoCommuteT(int what, int &Sendxm, int &Sendxp, int &Sendym, int &Sendyp);
   void    UnPackT(int what, int  Recexm, int  Recexp, int  Receym, int  Receyp);

   void DoCommute(int what, int k);
   void    DoPack(int what, int k);
   void    UnPack(int what, int k);

   int Get_bufsize() {return bufsize;};



   Commute(int XGridN, int YGridN);
   ~Commute();


};

#endif
