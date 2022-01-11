//----------------------------------------------------------------------------------||
//-------------------                mesh.h                      -------------------||
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
//---Starting---------           : Jan-24-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||



#ifndef H_MESH
#define H_MESH

#include <stdlib.h>
#include<vector>
class Detector;
//---------------------------- Mesh class -----------------------
class Mesh : public NList {

   friend class Domain;

private:
   
   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.

   WDOUBLE dx, dy, dz, dt, dt0;
   WDOUBLE dzz;


   int GridX,GridY,GridZ; // Grid number in X, Y, Z;
   int CellNum;           // Cell number in X-Y-Z domain;
   int TrajNum;           // Traj number in X-Y transverse plane;

   int RankIdx_X;
   int RankIdx_Y;

   int Rank;

   WDOUBLE Offset_X;     //Defined as the Mesh Bottom_Left Corner.
   WDOUBLE Offset_Y;     //Defined as the Mesh Bottom_Left Corner.

//====plasma=======
   int PProfileT;       //Plasma transverse profile
   int PProfileL;       //Plasma longitudinal profile
   int TpCellx;         //Trajectory per cell in x direction
   int TpCelly;         //Trajectory per cell in y direction
   int PushOrder;       //Push Trajectory Order in z direction
   WDOUBLE Delta_P;
   WDOUBLE AdaptiveStep;    //Adaptive Z-step for pushing trajectories
   WDOUBLE V_thresh;
   WDOUBLE Vmax;
   WDOUBLE Vlim;

   int IfAdjustPsi;

//===Possible Ionization
   int if_ioniz;      // if there is a need for ionization.
   WDOUBLE Ion_R;
   WDOUBLE Pla_ne;
   WDOUBLE Dop_ne;
   WDOUBLE Dop_TB;   //time begine
   WDOUBLE Dop_TE;   //time end

   dcomplex *dAx1;     //store the difference A(t+1)-A(t)
   dcomplex *dAx2;     //store the difference A(t+1)-A(t)
   dcomplex *dAy1;     //store the difference A(t+1)-A(t)
   dcomplex *dAy2;     //store the difference A(t+1)-A(t)
   WDOUBLE dA0;

   Cell*        p_CellArray;
   Trajectory*  p_Trajectory;
   Particle*    p_Particle;
   Detector*    XRayDetector;

   WDOUBLE minGamma; // Smallest Gamma in domain
   
public:


   WDOUBLE PlateauBegin;
   WDOUBLE PlasmaBegin; 
   WDOUBLE PlateauEnd; 
   WDOUBLE PlasmaEnd; 
   WDOUBLE PlasRadius;  
   std::vector<WDOUBLE> GridsTmp;

   int GetRankIdx_X() {return RankIdx_X;};
   int GetRankIdx_Y() {return RankIdx_Y;};

   WDOUBLE GetOffset_X() {return Offset_X;};
   WDOUBLE GetOffset_Y() {return Offset_Y;};

   int GetPushOrder() {return PushOrder;};
   
   inline int GetCellN(int i, int j, int k)  //inline avoid function calling overhead.
   {
      return i+(GridX+2)*(j+(GridY+2)*k);
   };

   inline Cell& GetCell(int i, int j, int k)
   {
      return p_CellArray[GetCellN(i,j,k)];
   };

   inline Cell& GetCell(int n)
   {
       return p_CellArray[n];
   };


   //Get Physical Coordinates of the Cells in the Mesh;
   inline WDOUBLE CellX(int i) //not useful in this version
   {
      return dx*(i - 0.5)+Offset_X;
   };

   inline WDOUBLE CellY(int j) //not useful in this version
   {
      return dy*(j - 0.5)+Offset_Y;
   };



   inline WDOUBLE CellZ(int k)
   {
      return dz*k;
   };

   void InitPulse(Pulse*);
   void SetdAs();
   void SeedTrajectory();
   void SetIonDensity();
   void SeedParticles(Specie *specie);
   void ResetPlasma();
   void AddTrajectory(Trajectory* p_Traj);
   void   AddParticle(Particle*   p_Part);


   void  MacroSource(int k);
   void AdjustSource(int k);
   void AdjustFields(int k);
   void AdjustPsi(int k);
   void SetSourceZero(int k);
   void BeamSource();
   void AdjustBSource();
   void SetBSourceZero();
   void SetFieldZeroAfter(int k);



   void     Put_Chi(WDOUBLE k0, int k); 
   void Partial_Psi(int k); 
   void Pondermotive(int k);
   void       Put_Jz(int k); 
   void  LaserFields();
   void  Ionization();  // initial stage simple version
   int   Ifioniz(void) {return if_ioniz;};

   WDOUBLE DopeBegin(void) {return Dop_TB;};
   WDOUBLE DopeEnd(void)   {return Dop_TE;};

   WDOUBLE  Dive_J(int i, int j, WDOUBLE k0, int k); // divergence of transverse current
   WDOUBLE  Curl_J(int i, int j, WDOUBLE k0, int k); // Curl of transverse current
   
   WDOUBLE SourceX(int i, int j, WDOUBLE k0, int k);
   WDOUBLE SourceY(int i, int j, WDOUBLE k0, int k);

   int Get_TpCellx()
   {
      return TpCellx;
   }

   int Get_TpCelly()
   {
      return TpCelly;
   }

   WDOUBLE ProfileLongi(WDOUBLE xt, WDOUBLE yt, WDOUBLE zt);
   WDOUBLE ProfileTrans(WDOUBLE xt, WDOUBLE yt, WDOUBLE zt);
   void PushTrajectory(WDOUBLE k0, int k, int step);
   void PushTrajectory_Half();
   void PushTrajectory_HalfE(int k);
   void PushTrajectory_HalfB(int k);
   void PushParticle();
   void SetNewTimeStep();
   void AdjustZstep(WDOUBLE k0, int k, WDOUBLE &dz2dz);
   void ExchangeT();
   void ExchangeP();
   void PackT(Trajectory* p_Traj, WDOUBLE* &Se);
   void PackP(Particle*   p_Part, WDOUBLE* &Se);
   Trajectory* Reconnect(Trajectory* p_Traj);
   Particle*   Reconnect(Particle*   p_Part);
   dcomplex SourceAx(int i, int j, int k, int NF);
   dcomplex SourceAy(int i, int j, int k, int NF);
   void Put_dA12(int what,  int k, int NF);
   dcomplex GetdA0() {return dA0;}
   WDOUBLE rand_gaussian (WDOUBLE sigma);

   Mesh(int XGridN, int YGridN, int ZGridN, FILE *f);
   ~Mesh();
};




class Detector : public NList {

   friend class Domain;
   friend class Mesh;

private:
   
   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.

public:


   WDOUBLE* p_DetectArray;
   int IfRadiation;

   
   WDOUBLE ThetaMax;
   int NTheta;

   WDOUBLE PhiMax;
   int NPhi;

   WDOUBLE OmegaMin;
   WDOUBLE OmegaMax;
   int NOmega;  


   // i=omega j=theta k=phi
   inline int GetDectN(int i, int j, int k) 
   {
      return i+(NOmega)*(j+(NTheta)*k);
   };

   inline WDOUBLE GetDetector(int i, int j, int k)
   {
      return p_DetectArray[GetDectN(i,j,k)];
   };

   inline void PutPacket(int i, int j, int k, WDOUBLE EE)
   {
      p_DetectArray[GetDectN(i,j,k)]+=EE;
      return;
   };

   Detector(FILE *f);
   ~Detector();
};
#endif
