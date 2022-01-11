//----------------------------------------------------------------------------------||
//-------------------                cell.h                      -------------------||
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
//---Starting---------           : Jan-27-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||



#ifndef H_CELLS
#define H_CELLS


#include "wand_PIC.h"


//class Domain; // forward declaration of Domain
class Mesh;


class Cell
{

   friend class Mesh;
   friend class Commute;
   friend class Domain;
   friend class MultiGrid;


private:

   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.
//=====Wake Fields=======
   union
   {
      WDOUBLE W_Fields[1];
      WDOUBLE W_Psi;
   };
   // Ex and Ey here is actually d\psi/dx and d\psi/dy
   // Ponx is the pondermotive force.
   // Psi Ez, Bz, Bx, By is solved by MG.
   WDOUBLE W_Ez, W_Bz, W_Bx, W_By;
   WDOUBLE W_Ex, W_Ey, W_Ponx, W_Pony, W_Asq;

//=====Wake Currents and Densities=======
   union
   {
      WDOUBLE W_Source[1];
      WDOUBLE W_Denn;
   };
   WDOUBLE W_Jx, W_Jy, W_Jxx, W_Jyy, W_Jxy;
   WDOUBLE B_Den, B_Jx, B_Jy, B_Jz, B_Chi;   // beam density and current
   WDOUBLE W_Deni; // ion macro density
   WDOUBLE W_Chi;
   WDOUBLE W_Jz;

   // Variable Grid Size.
   WDOUBLE dx;
   WDOUBLE dy;

   WDOUBLE Xcord;
   WDOUBLE Ycord;


   dcomplex *Acomx;
   dcomplex *Acomy;
   dcomplex *Acomxm;  //old value;
   dcomplex *Acomym;  //old value;

   //=====Laser Fields=======
   dcomplex *L_Ex;
   dcomplex *L_Ey, *L_Ez, *L_Bx, *L_By, *L_Bz;



public:

   int InoState;
   WDOUBLE Z_shifted;
   
   void AddAComs(dcomplex acomx, dcomplex acomy, int NF)
   {
      Acomx[NF] += acomx;
      Acomy[NF] += acomy;
      Acomxm[NF] = Acomx[NF];
      Acomym[NF] = Acomy[NF];
   }



   Cell();
   ~Cell();
};

#endif
