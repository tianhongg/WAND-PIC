//----------------------------------------------------------------------------------||
//-------------------                domain.cpp                  -------------------||
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

#include "wand_PIC.h"
#include <ctime>

//Domain::p_D = NULL;
Domain* Domain::p_Domain= NULL;
//---------------------------- Domain::Domain -----------------------
Domain::Domain (char * infile, int rank) : NList("Domain") 
{

   Time = 0.0;

   ifAx = ifAy = 0;

   Domain::p_Domain = this; 

   AddEntry((char*)"XStep",		&dx,	1.0);
   AddEntry((char*)"YStep",		&dy,	1.0);
   AddEntry((char*)"ZStep",		&dz,	1.0);

   AddEntry((char*)"TStep",		&dt,	1.0);
   AddEntry((char*)"SubCycle",   &Ndt, 1  );

   AddEntry((char*)"XLength", 	&Xmax,	10.0);
   AddEntry((char*)"YLength", 	&Ymax,	10.0);
   AddEntry((char*)"ZLength", 	&Zmax,	10.0);

   AddEntry((char*)"TimeMax",      &Tmax,    10.0);
   AddEntry((char*)"OutputPeriod", &OutDt,   1.0);
   AddEntry((char*)"SaveDim",      &savedim,   0);

   AddEntry((char*)"Npulse",       &Npulse,   1);
   AddEntry((char*)"Nbeam",        &Nbeam,    0);
   AddEntry((char*)"Boundary",     &BC,       0);
   AddEntry((char*)"BufferSize",   &Buffersize,  2.0);
   AddEntry((char*)"Restart",      &IfRestart,0);

   Rank = rank;

   p_File = fopen(infile,"rt");

   if (p_File)
   {
      rewind(p_File);
      read(p_File);
   }

   //===============================================================
   //=======================Creating Partition======================
   //===============================================================
   if (rank==0)  std::cout << "==== Domain: Creating Partition.         ====\n";
   p_MPP = new Partition(p_File); //pointer point to partition class
   if (rank==0)  std::cout << "==== Domain: Partition Created.          ====\n";


   //===============================================================
   //=======================Creating Mesh     ======================
   //===============================================================
   // set number of grids
   XGridN = round(Xmax/(p_MPP->GetXpart()*1.0)/dx);
   YGridN = round(Ymax/(p_MPP->GetYpart()*1.0)/dy);
   ZGridN = round(Zmax/dz);


   //Round XY-Mesh to even number
   XGridN = ceil(XGridN/2.0)*2; 
   YGridN = ceil(YGridN/2.0)*2; 


   if( XGridN != YGridN || ( p_MPP->GetXpart() ) != ( p_MPP->GetYpart() ) )
   {
   if (rank==0)  std::cout << "==== Please Make X, Y Direction Identical====\n";
   exit(0);
   }


   //Resize the domain;
   Xmax=XGridN*dx*p_MPP->GetXpart();
   Ymax=YGridN*dy*p_MPP->GetYpart();
   Zmax=ZGridN*dz;


   //===============================================================
   //=======================Initiate Laser Pulse in the Domain======
   //===============================================================
   if (rank==0)  std::cout << "==== Domain: Initiating Pulses.          ====\n";
   int i;
   if (Npulse > 0) 
   {
      
      pp_Pulses = new Pulse*[Npulse];  
      //====================================
      //How Many frequencies possible in the domain.
      OmegaL = new double[Npulse];
      ifAx = new int[Npulse];
      ifAy = new int[Npulse];
      NFreqs = 0;

      for (i=0;i<Npulse;i++) 
      {ifAx[i]=ifAy[i]=OmegaL[i]=0.0;}
      //====================================

      char name[128];
      for (i=0; i<Npulse; i++)
      {
         sprintf(name,"Pulse%d",i);
         pp_Pulses[i] = new Pulse(name, p_File);
      }
      if (rank==0)  std::cout << "==== Domain: Pulses Parameters Are Read. ====\n";      

   }
   else 
   {
      pp_Pulses = NULL;
      OmegaL = NULL; ifAx=ifAy=NULL;
      if (rank==0)  std::cout << "==== Domain: No Pulse In Domain.         ====\n";

   }

   

   if (rank==0)  std::cout << "==== Domain: Creating Mesh.              ====\n";
   p_Meshes = new Mesh(XGridN,YGridN,ZGridN,p_File);
   if (rank==0)  std::cout << "==== Domain: Mesh Created.               ====\n";
 


   if (Npulse > 0) 
   {
      //Put Vector Potential A in the Domain
      for (i=0; i<Npulse; i++) 
      {
         p_Meshes->InitPulse(pp_Pulses[i]);
      }

      p_Meshes->SetdAs();
      if (rank==0)  std::cout << "==== Domain: Pulses are Initiated.       ====\n";



      //------Restart--------------
      if(IfRestart)
      {
         if (rank==0)  std::cout << "==== Domain: Restart... Loading files.   ====\n";
         LoadPulse(IfRestart);
         if (rank==0)        printf("==== Domain: Pusle File No.%4d Loaded.  ====\n",IfRestart);
         // reset the timer.
         Time=OutDt*IfRestart;

      }

   }




   //===============================================================
   //=======================Initiate Beam in the Domain  ===========
   //===============================================================
   NSpecie = 0;
   if(Nbeam>0)
   {
      SpecieType = new int[Nbeam];
      pp_Species = new Specie*[Nbeam]; 
      char name[128];
      for (int i=0; i<Nbeam; i++)
      {
         sprintf(name,"Beam%d",i);
         pp_Species[i] = new Specie(name, p_File);
         AddSpecie(pp_Species[i]->P_type);
      }

      
      if (rank==0)  std::cout << "==== Domain: Beam Parameters Are Read.   ====\n";

      for (int i=0; i<Nbeam; i++)
      {

         if(IfRestart==0)
         {
            p_Meshes->SeedParticles(pp_Species[i]);
            if (rank==0)  std::cout << "==== Domain: Particles Are Seeded.       ====\n";
         }
         
      }
      
   }
   else
   {
      if(p_Meshes->Ifioniz())
      {
         Nbeam=1;
         char name[128];
         SpecieType = new int[1];
         pp_Species = new Specie*[1];
         sprintf(name,"Beam%d",0);
         pp_Species[0] = new Specie(name, p_File);
         pp_Species[0]->density=0;
         AddSpecie(0);


      }
      pp_Species = NULL;
      if (rank==0)  std::cout << "==== Domain: No Beam In Domain.          ====\n";

   }
   if(IfRestart)
   {
      LoadParti(IfRestart);
      if (rank==0)  printf("==== Domain: ParticleFile No.%4d Loaded.====\n",IfRestart);
      Time=OutDt*IfRestart;
   }

   //===============================================================
   //=======================Seed Trajectory in the Domain===========
   //===============================================================
   if (rank==0)  std::cout << "==== Domain: Seed Trajectories in Mesh.  ====\n";
   p_Meshes->SeedTrajectory();
   if (rank==0)  std::cout << "==== Domain: Trajectories Are Seeded.    ====\n";

   trajorder=p_Meshes->GetPushOrder();

   if(trajorder==0 || trajorder==1 || trajorder==2)
   {   
   }
   else
   {
   if (rank==0)  std::cout << "==== Domain: Wrong Traj Pushing Order.   ====\n";
   exit(12);
   }

   //===============================================================
   //=======================  Create Communication.      ===========
   //===============================================================
   p_Comm = new Commute(XGridN,YGridN);
   if (rank==0)  std::cout << "==== Domain: Communications Established. ====\n";

   //===============================================================
   //=======================  Create MultiGrid Meshes    ===========
   //===============================================================

   p_Multi = new MultiGrid(rank, XGridN, YGridN, p_File);
   if (rank==0)  std::cout << "==== Domain: MultiGrid Meshes Created.   ====\n";
   if (rank==0)  std::cout << "=============================================\n\n\n\n";


}



void Domain::Run()
{

   //=======Timer=======
   std::clock_t start;
   start = std::clock();
   double duration;
   //====================

   int n=0;
   int nc=0;
   int k=0;
   double k0=0;
   int ierr;
   
   // restart;
   if(IfRestart) nc=IfRestart+1;

//======================================
//============= Main Routine ===========
while(Time<Tmax)
{

   if(Nbeam) p_Meshes -> BeamSource();
   
   //===================================
   //======== Push Wakefield in z=======
   k=0; k0=0; ierr=1;
   switch(trajorder)
   {
      case 0:
      while (k<ZGridN && ierr) { ierr = Boris(k0, k);};
      break;

      case 1:
      while (k<ZGridN && ierr) { ierr = RK1(k0, k);};
      break;

      case 2:
      while (k<ZGridN && ierr) { ierr = RK2(k0, k);};
      break;
   }
   if(ierr == 0)  p_Meshes -> SetFieldZeroAfter(k);
   //===================================


   //===================================
   //======get laser field =============
   if(NFreqs && Nbeam) p_Meshes -> LaserFields();
   //===================================



   
   //===================================
   //=== here do ionization block ======
   if(p_Meshes->Ifioniz() && Time>p_Meshes->DopeBegin() 
      && Time<p_Meshes->DopeEnd())
   {
      p_Meshes ->Ionization(); //tianhong_sep_15 test_version
   }
   //===================================


   //==========push particles ==========
   if(Nbeam || p_Meshes->Ifioniz()) p_Meshes -> PushParticle();
   //===================================

   //==================save==============
   duration=(std::clock()-start)/(double)CLOCKS_PER_SEC/60;
   if(Rank==0) printf("==== Step: %8d ----  %7.3f Mins.  ====\n",n,duration);
   Time += dt;
   n++;
   //===========================


   //===========================
   //======== Output ===========
   if(Time>=OutDt*nc)
   {
      if(savedim==0)
      {Save(nc);}
      else
      {Save2D(nc,savedim);}

      nc++;
      duration=(std::clock()-start)/(double)CLOCKS_PER_SEC/60;

      if(Rank==0)
      { 
         printf("\n================ Output: %4d.  =============\n",nc-1);
         printf("==== %4.1f%% == Time-Elapsed: %5.1f Mins.  ====\n",(Time/Tmax)*100,duration);
         printf("==== Estimated Time Remain: %4.2f Hours.  ====\n\n",duration*(Tmax-Time)/(Time-OutDt*IfRestart)/60);
      }

   }
   //==================save================

   //===========================
   //======== Reset Plasma =====
   p_Meshes ->ResetPlasma();
   //===========================

   
}

//============= Main Routine ===========
//======================================


   return;
}



void Domain::AddSpecie(int P_type)
{

   int i;
   int temp=1;

   if(NSpecie == 0)
   {
      SpecieType[0] = P_type;
      NSpecie++;
   }

   else
   {
      for(i=0; i<NSpecie; i++)
      {
         if(SpecieType[i]==P_type)
         {
            temp=0;
         }

      }

      if(temp==1)
      {
         NSpecie++;
         SpecieType[NSpecie-1]=P_type;
      }
   }

   return;
}


int Domain::Get_NSpecie(int SpecieType)
{
   int num=0;
   Particle *p = NULL;
   p = p_Mesh() -> p_Particle;

   while(p)
   {

      if(p->type==SpecieType) num++;
      p = p->p_PrevPart;

   }
   return num;

}

//---------------------------- Domain::~Domain() -----------------------
Domain::~Domain()
{
      if(p_File) fclose(p_File); 
      if(Rank_cout) Rank_cout.close();
      delete[] pp_Pulses;
      delete p_Comm;
      delete p_Multi;
      delete p_Meshes;
      delete p_MPP;
      delete[] pp_Species;
      delete[] SpecieType;

};
