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
    double dx;
    double dy;
    double dz;
    double dt;

    int Ndt;

    double Xmax;
    double Ymax;
    double Zmax;

    int XGridN;
    int YGridN;
    int ZGridN;

    double Tmax;
    double OutDt;
    double Time;
    int savedim;

    int Npulse;
    int Nbeam;
    int BC;
    double Buffersize;
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
    double *OmegaL;
    int NFreqs;
    
    double Get_dx() {return dx;};
    double Get_dy() {return dy;};
    double Get_dz() {return dz;};
    double Get_dt() {return dt;};
    
    int Get_SubCycle() {return Ndt;};


    double Get_Xmax() {return Xmax;};
    double Get_Ymax() {return Ymax;};
    double Get_Zmax() {return Zmax;};

    double Get_RunTime() {return Time;};
    int Get_BC()         {return BC;}
    double Get_Buffersize() {return Buffersize;}

    int Get_Nbeam()   {return Nbeam;};

    // ========run===========
    int RK1(double &k0, int &k);
    int RK2(double &k0, int &k);
    int Boris(double &k0, int &k);
    int PushWakeFields( int k);
    int PushWakeFieldsE( int k);
    int PushWakeFieldsEz( int k);
    int PushWakeFieldsB( int k);
    int PushPulses(int k, int NF);
    void Run();



    Partition *p_Partition()    {return p_MPP;}; 
    Mesh      *p_Mesh()         {return p_Meshes;}; 
    MultiGrid *p_MG()           {return p_Multi;}; 
    Pulse    **GetPulses(void)  {return pp_Pulses;};
    Commute    *p_Com()         {return p_Comm;}


    int  Save(int nt);
    int  Save2D(int nt, int savedim);
    int  SaveP(int nt);
    int  SaveXray(int nt);
    int  LoadPulse(int nt);
    int  LoadParti(int nt);
    void AddSpecie(int P_type);

    int Get_NSpecie(int SpecieType);




  //=====constructor-destructor=====//
    Domain(char *infile, int rank);  
    ~Domain();        

};

#endif
