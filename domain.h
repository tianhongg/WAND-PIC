//----------------------------------------------------------------------------------||
//-------------------                domain.h                    -------------------||
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



#ifndef H_DOMAIN
#define H_DOMAIN



//=====include======//
#include <stdio.h>
#include <iostream>
#include <fstream>




#define PI 3.1415926536
//=====Domain class=====//

class Partition;
class Pulse;
class Mesh;
class Commute;
class MultiGrid;
class Specie;

class Domain : public NList{

  public:
    static Domain *p_Domain;
    std::ofstream Rank_cout;
 
  private: 
  //variables
    WDOUBLE dx;
    WDOUBLE dy;
    WDOUBLE dz;
    WDOUBLE dt;

    // for variable mesh size
    int MeshType;
    WDOUBLE delta;
    WDOUBLE dxRefine;
    WDOUBLE radius0;
    WDOUBLE order;


    int Adap_dt;

    int Ndt;

    WDOUBLE Xmax;
    WDOUBLE Ymax;
    WDOUBLE Zmax;

    int XGridN;
    int YGridN;
    int ZGridN;

    WDOUBLE Tmax;
    WDOUBLE OutDt;
    WDOUBLE Time;
    // int Tstep;
    int savedim;

    int Npulse;
    int Nbeam;
    int BC;
    WDOUBLE Buffersize;
    int IfRestart; 
    int trajorder;

    int Rank;

    
    int NSpecie;    
    int *SpecieType;

    FILE *p_File;
    Partition *p_MPP;
    Pulse **pp_Pulses;
    Specie **pp_Species;
    Mesh *p_Meshes;
    Commute *p_Comm;
    MultiGrid *p_Multi;

  public:    // Public functions

    int *ifAx;
    int *ifAy;
    WDOUBLE *OmegaL;
    int NFreqs;
    WDOUBLE Q1;
    WDOUBLE Get_dx() {return dx;};
    WDOUBLE Get_dy() {return dy;};
    WDOUBLE Get_dz() {return dz;};
    WDOUBLE Get_dt() {return dt;};
    void set_new_dt(WDOUBLE newdt) {dt=newdt;};
    
    int Get_SubCycle() {return Ndt;};


    WDOUBLE Get_Xmax() {return Xmax;};
    WDOUBLE Get_Ymax() {return Ymax;};

    WDOUBLE Set_Xmax(WDOUBLE x) {return Xmax=x;};
    WDOUBLE Set_Ymax(WDOUBLE y) {return Ymax=y;};

    WDOUBLE Get_Zmax() {return Zmax;};
    WDOUBLE Get_RunTime() {return Time;};
    // int Get_Step() {return Tstep;};
    
    int Get_BC()         {return BC;}
    WDOUBLE Get_Buffersize() {return Buffersize;}

    int Get_Nbeam()   {return Nbeam;};

    // ========run===========
    int RK1(WDOUBLE &k0, int &k);
    int RK2(WDOUBLE &k0, int &k);
    int Boris(WDOUBLE &k0, int &k);
    int PushWakeFields(WDOUBLE k0, int k);
    int PushWakeFieldsE(WDOUBLE k0, int k);
    int PushWakeFieldsEz(WDOUBLE k0, int k);
    int PushWakeFieldsB(WDOUBLE k0, int k);
    int PushPulses(int k, int NF);
    void Run();



    Partition *p_Partition()    {return p_MPP;}; 
    Mesh      *p_Mesh()         {return p_Meshes;}; 
    MultiGrid *p_MG()           {return p_Multi;}; 
    Pulse    **GetPulses(void)  {return pp_Pulses;};
    Commute    *p_Com()         {return p_Comm;}


    int  Save(int nt);
    int  Save2D(int nt, int savedim, bool part);
    int  SaveP(int nt);
    int  SaveT(int nt, int k);
    int  SaveXray(int nt);
    int  LoadPulse(int nt);
    int  LoadParti(int nt);
    void AddSpecie(int P_type);

    int Get_NSpecie(int SpecieType);
    WDOUBLE CustomGrid(WDOUBLE r);




  //=====constructor-destructor=====//
    Domain(char *infile, int rank);  
    ~Domain();        

};

#endif
