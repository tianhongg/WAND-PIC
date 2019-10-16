//----------------------------------------------------------------------------------||
//-------------------                stepscpp                    -------------------||
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
//---Starting---------           : Apr-08-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#include "wand_PIC.h"


int Domain::RK2(double &k0, int &k)
{
   int ierr=0;
   double dz2dz;
   int NF;


   p_Meshes ->AdjustZstep(k0, k, dz2dz);


   //================================
   //========  RK2  ======
   if((ierr = PushWakeFields(k))) return 0;
   p_Meshes ->PushTrajectory(k0, k, 1);
   if((ierr = PushWakeFields(k))) return 0;
   p_Meshes ->PushTrajectory(k0, k, 0);

   // p_Meshes ->AdjustZstep(k0, k);

   //================================
   k  = ceil(k0 + dz2dz);
	k0 = 	 (k0 + dz2dz); 



   //================================
   //================================
   //======== Push Laser Ax, Ay ======
   if(NFreqs>0 &&  k0==k)
   {
      for(NF=0; NF<NFreqs; NF++)
      {
         if((ierr = PushPulses(k-1, NF))) return 0;
      }
   }
   //================================

   return 1;

}


int Domain::RK1(double &k0, int &k)
{
   int ierr=0;
   int NF;
   double dz2dz;


   p_Meshes ->AdjustZstep(k0, k, dz2dz);

   //================================
   //======== Push Wakefields  ======
   if((ierr = PushWakeFields(k))) return 0;
   //================================
   //======== Push Trajectory  ======
   //======== and Exchange Traj =====
   
   p_Meshes ->PushTrajectory(k0, k,0);
   // p_Meshes ->AdjustZstep(k0, k);

   //================================
  	k  = ceil(k0 + dz2dz);
  	k0 = 	   (k0 + dz2dz); 

   //================================
   //================================
   //======== Push Laser Ax, Ay ======
   if(NFreqs>0 &&  k0==k)
   {
      for(NF=0; NF<NFreqs; NF++)
      {
         if((ierr = PushPulses(k-1, NF))) return 0;
      }
   }
   //================================

   return 1;

}


int Domain:: PushPulses(int k, int NF)
{
   int ierr = 0;

   if(ifAx[NF])
   {
      ierr = p_Multi-> MG_V_cycleC(5, k, NF);
      p_Meshes ->         Put_dA12(5, k, NF);
   }
   if(ifAy[NF])
   {
      ierr = p_Multi-> MG_V_cycleC(6, k, NF);
      p_Meshes ->         Put_dA12(6, k, NF);
   }

   return ierr;

}


int Domain:: PushWakeFields(int k)
{
   int ierr=0;
   //================================
   //========  Plasma Source  =======
   p_Meshes -> MacroSource(k);
   p_Comm   -> DoCommute(COMMU_S, k);
   //================================
   //========  Wake: Psi   ==========
   if( (ierr = p_Multi  -> MG_V_cycle(0, k)) )  return 1;
   p_Meshes -> AdjustPsi(k);
   //==== Psi Exchanged After MG=====
   //================================
   //========  Wake: Chi   ==========
   p_Meshes -> Put_Chi(k);
   //================================
   //========  Wake: dPsi/dx ========
   p_Meshes -> Partial_Psi(k);
   //================================
   //========  Wake: Ez    ==========
   if( (ierr = p_Multi  -> MG_V_cycle(1, k)) )  return 1;
   //==== Ez Exchanged After MG======
   //================================
   //========  Wake: Bz    ==========
   if( (ierr = p_Multi  -> MG_V_cycle(2, k)) )  return 1;
   //==== Bz Exchanged After MG======
   //================================
   //======= Wake: d|A|^2/dx ========

   if(NFreqs>0) p_Meshes -> Pondermotive(k);
   //================================
   //========  Wake: Jz  ============
   p_Meshes -> Put_Jz(k);
   //================================
   //========  Wake: Bx    ==========
   if( (ierr = p_Multi  -> MG_V_cycle(3, k)) )  return 1;
   //==== Bx Exchanged After MG======
   //================================
   //========  Wake: By    ==========
   if( (ierr = p_Multi  -> MG_V_cycle(4, k)) )  return 1;
   //==== By Exchanged After MG======
   //================================
   //======== Exchange Wakefield  ===
   //======== Ex Ey Ponx Pony =======
   p_Comm   -> DoCommute(COMMU_F, k);
   p_Meshes ->AdjustFields(k);


   return 0;


}