//----------------------------------------------------------------------------------||
//-------------------                cell.cpp                    -------------------||
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



#include <stdio.h>
#include <stdlib.h>

#include "wand_PIC.h"

//=======Cell Class
Cell::Cell()
{



      W_Ex = W_Ey = W_Ez = 0.0;
      W_Bx = W_By = W_Bz = 0.0;

      W_Denn = W_Jx  = W_Jy  = 0.0;
      W_Jxx  = W_Jyy = W_Jxy = 0.0;

      B_Den = B_Jx = B_Jy = B_Jz = B_Chi = 0.0; 

      W_Deni = 1.0;
      W_Psi  = 0.0;

      W_Ponx = W_Pony =0.0;


//===========================================
      int NF =  p_domain()->NFreqs;
      L_Ex =  L_Ey =  L_Ez = NULL;
      L_Bx =  L_By =  L_Bz = NULL;
      Acomx  = Acomy  = NULL;
      Acomxm = Acomym = NULL;

      if(NF>0)
      {
         L_Ex   = new dcomplex[ NF];
         L_Ey   = new dcomplex[ NF];
         L_Ez   = new dcomplex[ NF];

         L_Bx   = new dcomplex[ NF];
         L_By   = new dcomplex[ NF];
         L_Bz   = new dcomplex[ NF];

         Acomx  = new dcomplex[ NF];
         Acomy  = new dcomplex[ NF];

         Acomxm = new dcomplex[ NF];
         Acomym = new dcomplex[ NF];

         for (int i=0;i<NF;i++)
         {
            L_Ex[i] =  L_Ey[i] =  L_Ez[i] = 0.0;
            L_Bx[i] =  L_By[i] =  L_Bz[i] = 0.0;
            Acomx[i]  = Acomy[i]  = Acomxm[i] = Acomym[i] = 0.0;

         }
      }

      InoState=1;


};



//---Cell::~Cell--------------------------------------------->
Cell::~Cell()
{}
