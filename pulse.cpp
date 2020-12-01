//----------------------------------------------------------------------------------||
//-------------------                pulse.cpp                   -------------------||
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

#include <stdio.h>
#include <stdlib.h>
#include "wand_PIC.h"


//---------------------------- Pulse::Pulse -----------------------
Pulse::Pulse (char *name, FILE *f) : NList (name)
{

   AddEntry("a0", &a0, 1.0);
   AddEntry("Xpol", &Xpol, 0.);
   AddEntry("Ypol", &Ypol, 1.);

   AddEntry("XSpotSize", &XSpotSize, 1.0);
   AddEntry("YSpotSize", &YSpotSize, 1.0);
   AddEntry("Length",    &ZLength,   1.0);

   AddEntry("RiseTime", &ZRise, 0.);
   AddEntry("DropTime", &ZDrop, 0.);

   AddEntry("Xcenter", &Xcenter, 0.);
   AddEntry("Ycenter", &Ycenter, 0.);
   AddEntry("Zcenter", &Zcenter, 0.);

   AddEntry("FocalPlane", &ZFocal, 0.);

   AddEntry("Xphase", &Xphase, 0.);
   AddEntry("Yphase", &Yphase, 0.);

   AddEntry("Tprofile", &Tprofile, 0);
   AddEntry("Lprofile", &Lprofile, 0);

   AddEntry("Omega", &Omega, 10.);
   AddEntry("Chirp", &Chirp, 0.);


   if (f)
   {
      rewind(f);
      read(f);
   };

   int i, NF, ifp;
   ifp=0;

   //====== store frequency of the pulse.
   if(p_domain()->NFreqs==0)
   {
      p_domain()->NFreqs =1;
      NF=0;

   }
   else
   {

      for(i=0; i<p_domain()->NFreqs; i++)
      {
         if(p_domain()->OmegaL[i] == Omega)
         {
            ifp=1; NF=i;
            break;
         }

      }
      if(ifp==0) {NF = p_domain()->NFreqs; p_domain()->NFreqs++;}

   }

   p_domain()->OmegaL[NF]=Omega;
   p_domain()->ifAx[NF] += Xpol;
   p_domain()->ifAy[NF] += Ypol;


}

 
//---------------------------- Mesh::Mesh --------------------
void Mesh::InitPulse(Pulse *pulse)
{

   double Omega =pulse->Omega;
   int i, NF;
   //====== store frequency of the pulse.
   for(i=0; i<p_domain()->NFreqs; i++)
   {
        if(p_domain()->OmegaL[i] == Omega)
       {
          NF=i;
          break;
      }

   }

   for (int k=0; k<GridZ; k++) 
   {
      double z = CellZ(k);
      for (int j=0; j<GridY+2; j++) 
      {
         double y = CellY(j);
         for (int i=0; i<GridX+2; i++)
         {
               double x = CellX(i);
               Cell &c = GetCell(i, j, k);
               pulse->Put_AComplex(x,y,z,c, NF);
         }
      }
   }


   return;
}



void Pulse::Put_AComplex(double x, double y, double z, Cell &c, int NF)
{

   x -= Xcenter;
   y -= Ycenter;
   z -= Zcenter;
   
   double phasex;
   double phasey;
   double phase0;


   double ZRx = 0.5*XSpotSize*XSpotSize*Omega;
   double ZRy = 0.5*YSpotSize*YSpotSize*Omega;


   double RCurva =ZFocal-Zcenter;


   if (RCurva==0)
   {
      phase0 = (Chirp*z)*z;

   }
   else
   {
      phase0 = (Chirp*z)*z - Omega*(x*x+y*y)/RCurva/2;
 
   }

   phasex = phase0 + Xphase;
   phasey = phase0 + Yphase;


   Acomx = a0*Xpol*ProfileLongi(x,y,z)*ProfileTrans(x,y,z)*exp(ci*phasex);
   Acomy = a0*Ypol*ProfileLongi(x,y,z)*ProfileTrans(x,y,z)*exp(ci*phasey);


   c.AddAComs(Acomx, Acomy, NF);

   return;
}



double Pulse::ProfileLongi(double x, double y, double z)
{

   if (ZLength<=0) return 0.;

   double arg = 0.;
   double amp = 0.;

   arg = (z/ZLength)*(z/ZLength);

   switch(Lprofile) 
   {

      case 0:
      default:

      if ( arg>4.) 
      {
        return 0.;
      }
     
      amp = exp(-arg);       // default Gaussian pulse
      return amp;

   }
   //  }
   return amp;
}


double Pulse::ProfileTrans(double x, double y, double z)
{

   if (XSpotSize<=0 || YSpotSize<=0) return 0.;
   
   double arg = 0.;
   double amp = 0.;

   arg = (x/XSpotSize)*(x/XSpotSize)+(y/YSpotSize)*(y/YSpotSize);

   switch(Tprofile) 
   {
      case 0:
      default:

      if (arg>8.)
      {
         return 0.;
      }

      amp = exp(-arg);       // default Gaussian pulse
      return amp;

   }

   //  }
   return amp;
}



void Mesh::LaserFields()
{

   int i,j,k, NF;
   double k0;
   int ifax, ifay;
   int NFreqs=p_domain()->NFreqs;


   for (k=0; k<GridZ-1; k++)
   {
      for (j=1; j<=GridY; j++)
      {
         for (i=1; i<=GridX; i++)
         {

         //  Fields are shifted in z by +dz/2;
         Cell &ccc   = GetCell(i,   j, k);
         Cell &czp   = GetCell(i,   j, k+1);

         Cell &cxm   = GetCell(i-1, j, k);
         Cell &cxp   = GetCell(i+1, j, k);
         Cell &cym   = GetCell(i, j-1, k);
         Cell &cyp   = GetCell(i, j+1, k);

         Cell &cxmzp = GetCell(i-1, j, k+1);
         Cell &cxpzp = GetCell(i+1, j, k+1);
         Cell &cymzp = GetCell(i, j-1, k+1);
         Cell &cypzp = GetCell(i, j+1, k+1);

         for(NF=0; NF<NFreqs; NF++)
         {
            ccc.L_Ex[NF] = ccc.L_Ey[NF] = ccc.L_Ez[NF] = 0.0;
            ccc.L_Bx[NF] = ccc.L_By[NF] = ccc.L_Bz[NF] = 0.0;

            ifax = p_domain()->ifAx[NF];
            ifay = p_domain()->ifAy[NF];
            k0   = p_domain()->OmegaL[NF];
         //===============================
            if(ifax)
            {
            ccc.L_By[NF] =ci*k0*0.25*( ccc.Acomx[NF]+ccc.Acomxm[NF] + czp.Acomx[NF]+czp.Acomxm[NF] )
                                -0.5*( czp.Acomx[NF]+czp.Acomxm[NF] - ccc.Acomx[NF]-ccc.Acomxm[NF] )/dz;

            ccc.L_Ex[NF] = ccc.L_By[NF]-0.5*( ccc.Acomx[NF]-ccc.Acomxm[NF] + czp.Acomx[NF]-czp.Acomxm[NF] )/dt;
   
            ccc.L_Ez[NF] = -0.125*( cxp.Acomx[NF]+  cxp.Acomxm[NF] -   cxm.Acomx[NF]-  cxm.Acomxm[NF]
                                 +cxpzp.Acomx[NF]+cxpzp.Acomxm[NF] - cxmzp.Acomx[NF]-cxmzp.Acomxm[NF])/dx;

            ccc.L_Bz[NF] = -0.125*( cyp.Acomx[NF]+  cyp.Acomxm[NF] -   cym.Acomx[NF]-  cym.Acomxm[NF]
                                 +cypzp.Acomx[NF]+cypzp.Acomxm[NF] - cymzp.Acomx[NF]-cymzp.Acomxm[NF])/dy;
            }
         //===============================

         //===============================
            if(ifay)
            {
            ccc.L_Bx[NF] = -ci*k0*0.25*( ccc.Acomy[NF]+ccc.Acomym[NF] + czp.Acomy[NF]+czp.Acomym[NF] )
                                  +0.5*( czp.Acomy[NF]+czp.Acomym[NF] - ccc.Acomy[NF]-ccc.Acomym[NF] )/dz;
   
            ccc.L_Ey[NF] =-ccc.L_Bx[NF]-0.5*( ccc.Acomy[NF]-ccc.Acomym[NF] + czp.Acomy[NF]-czp.Acomym[NF] )/dt;

            ccc.L_Ez[NF]+= -0.125*( cyp.Acomy[NF]+  cyp.Acomym[NF] -   cym.Acomy[NF]-  cym.Acomym[NF]
                                 +cypzp.Acomy[NF]+cypzp.Acomym[NF] - cymzp.Acomy[NF]-cymzp.Acomym[NF])/dy;

            ccc.L_Bz[NF]+=  0.125*( cxp.Acomy[NF]+  cxp.Acomym[NF] -   cxm.Acomy[NF]-  cxm.Acomym[NF]
                                 +cxpzp.Acomy[NF]+cxpzp.Acomym[NF] - cxmzp.Acomy[NF]-cxmzp.Acomym[NF])/dx;
            }
          }
         //===============================

         }
      }
   }


//adjust//

   for (k=0; k<GridZ-1; k++)
   {

      //=============================================================
         for (i=1; i<=GridX; i++)
         {
            Cell &c1 = GetCell(i, 0, k);  Cell &c2 = GetCell(i, 1, k);
            for(NF=0; NF<NFreqs; NF++)
            {
               c1.L_Ex[NF]=c2.L_Ex[NF]; c1.L_Ey[NF]=c2.L_Ey[NF]; c1.L_Ez[NF]=c2.L_Ez[NF];
               c1.L_Bx[NF]=c2.L_Bx[NF]; c1.L_By[NF]=c2.L_By[NF]; c1.L_Bz[NF]=c2.L_Bz[NF];
            }

         }
      //=============================================================
         for (i=1; i<=GridX; i++)
         {
            Cell &c1 = GetCell(i, GridY+1, k);  Cell &c2 = GetCell(i, GridY, k);
            for(NF=0; NF<NFreqs; NF++)
            {
               c1.L_Ex[NF]=c2.L_Ex[NF]; c1.L_Ey[NF]=c2.L_Ey[NF]; c1.L_Ez[NF]=c2.L_Ez[NF];
               c1.L_Bx[NF]=c2.L_Bx[NF]; c1.L_By[NF]=c2.L_By[NF]; c1.L_Bz[NF]=c2.L_Bz[NF];
            }

         }

      //=============================================================
         for (j=1; j<=GridY; j++)
         {
            Cell &c1 = GetCell(GridX+1, j, k);  Cell &c2 = GetCell(GridX, j, k);
            for(NF=0; NF<NFreqs; NF++)
            {
               c1.L_Ex[NF]=c2.L_Ex[NF]; c1.L_Ey[NF]=c2.L_Ey[NF]; c1.L_Ez[NF]=c2.L_Ez[NF];
               c1.L_Bx[NF]=c2.L_Bx[NF]; c1.L_By[NF]=c2.L_By[NF]; c1.L_Bz[NF]=c2.L_Bz[NF];
            }

         }
      //=============================================================
         for (j=1; j<=GridY; j++)
         {
            Cell &c1 = GetCell(0, j, k);  Cell &c2 = GetCell(1, j, k);
            for(NF=0; NF<NFreqs; NF++)
            {
               c1.L_Ex[NF]=c2.L_Ex[NF]; c1.L_Ey[NF]=c2.L_Ey[NF]; c1.L_Ez[NF]=c2.L_Ez[NF];
               c1.L_Bx[NF]=c2.L_Bx[NF]; c1.L_By[NF]=c2.L_By[NF]; c1.L_Bz[NF]=c2.L_Bz[NF];
            }

         }

      //=============================================================
      Cell &a  = GetCell(0, 0, k);  Cell &a1 = GetCell(1, 0, k);
      Cell &a2 = GetCell(0, 1, k);  Cell &a3 = GetCell(1, 1, k);
      for(NF=0; NF<NFreqs; NF++)
      {
         a.L_Ex[NF] = a1.L_Ex[NF]*0.4+a2.L_Ex[NF]*0.4+a3.L_Ex[NF]*0.2;
         a.L_Ey[NF] = a1.L_Ey[NF]*0.4+a2.L_Ey[NF]*0.4+a3.L_Ey[NF]*0.2;
         a.L_Ez[NF] = a1.L_Ez[NF]*0.4+a2.L_Ez[NF]*0.4+a3.L_Ez[NF]*0.2;
         a.L_Bx[NF] = a1.L_Bx[NF]*0.4+a2.L_Bx[NF]*0.4+a3.L_Bx[NF]*0.2;
         a.L_By[NF] = a1.L_By[NF]*0.4+a2.L_By[NF]*0.4+a3.L_By[NF]*0.2;
         a.L_Bz[NF] = a1.L_Bz[NF]*0.4+a2.L_Bz[NF]*0.4+a3.L_Bz[NF]*0.2;
      }

      //=============================================================
      Cell &b  = GetCell(0, GridY+1, k);  Cell &b1 = GetCell(1, GridY+1, k);
      Cell &b2 = GetCell(0, GridY,   k);  Cell &b3 = GetCell(1, GridY,   k);
      for(NF=0; NF<NFreqs; NF++)
      {
         b.L_Ex[NF] = b1.L_Ex[NF]*0.4+b2.L_Ex[NF]*0.4+b3.L_Ex[NF]*0.2;
         b.L_Ey[NF] = b1.L_Ey[NF]*0.4+b2.L_Ey[NF]*0.4+b3.L_Ey[NF]*0.2;
         b.L_Ez[NF] = b1.L_Ez[NF]*0.4+b2.L_Ez[NF]*0.4+b3.L_Ez[NF]*0.2;
         b.L_Bx[NF] = b1.L_Bx[NF]*0.4+b2.L_Bx[NF]*0.4+b3.L_Bx[NF]*0.2;
         b.L_By[NF] = b1.L_By[NF]*0.4+b2.L_By[NF]*0.4+b3.L_By[NF]*0.2;
         b.L_Bz[NF] = b1.L_Bz[NF]*0.4+b2.L_Bz[NF]*0.4+b3.L_Bz[NF]*0.2;
      }

//=============================================================
      Cell &c  = GetCell(GridX+1, 0, k);  Cell &c1 = GetCell(GridX, 0, k);
      Cell &c2 = GetCell(GridX+1, 1, k);  Cell &c3 = GetCell(GridX, 1, k);
      for(NF=0; NF<NFreqs; NF++)
      {
         c.L_Ex[NF] = c1.L_Ex[NF]*0.4+c2.L_Ex[NF]*0.4+c3.L_Ex[NF]*0.2;
         c.L_Ey[NF] = c1.L_Ey[NF]*0.4+c2.L_Ey[NF]*0.4+c3.L_Ey[NF]*0.2;
         c.L_Ez[NF] = c1.L_Ez[NF]*0.4+c2.L_Ez[NF]*0.4+c3.L_Ez[NF]*0.2;
         c.L_Bx[NF] = c1.L_Bx[NF]*0.4+c2.L_Bx[NF]*0.4+c3.L_Bx[NF]*0.2;
         c.L_By[NF] = c1.L_By[NF]*0.4+c2.L_By[NF]*0.4+c3.L_By[NF]*0.2;
         c.L_Bz[NF] = c1.L_Bz[NF]*0.4+c2.L_Bz[NF]*0.4+c3.L_Bz[NF]*0.2;
      }

//=============================================================
      Cell &d  = GetCell(GridX+1, GridY+1, k);  Cell &d1 = GetCell(GridX+1, GridY, k);
      Cell &d2 = GetCell(GridX,   GridY+1, k);  Cell &d3 = GetCell(GridX,   GridY, k);
      for(NF=0; NF<NFreqs; NF++)
      {
         d.L_Ex[NF] = d1.L_Ex[NF]*0.4+d2.L_Ex[NF]*0.4+d3.L_Ex[NF]*0.2;
         d.L_Ey[NF] = d1.L_Ey[NF]*0.4+d2.L_Ey[NF]*0.4+d3.L_Ey[NF]*0.2;
         d.L_Ez[NF] = d1.L_Ez[NF]*0.4+d2.L_Ez[NF]*0.4+d3.L_Ez[NF]*0.2;
         d.L_Bx[NF] = d1.L_Bx[NF]*0.4+d2.L_Bx[NF]*0.4+d3.L_Bx[NF]*0.2;
         d.L_By[NF] = d1.L_By[NF]*0.4+d2.L_By[NF]*0.4+d3.L_By[NF]*0.2;
         d.L_Bz[NF] = d1.L_Bz[NF]*0.4+d2.L_Bz[NF]*0.4+d3.L_Bz[NF]*0.2;
      }




   }

   return;
}


