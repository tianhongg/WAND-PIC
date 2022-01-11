//-------------------                trajectory.cpp              -------------------||
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
//---Starting---------           : Jan-29-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||
#include "wand_PIC.h"


Trajectory::Trajectory(WDOUBLE xt, WDOUBLE yt, WDOUBLE ztime, int TpCellx, int TpCelly, WDOUBLE sxx, WDOUBLE syy)
  {
  p_PrevTraj = p_NextTraj = NULL;
  x = old_x = x0 = xt;
  y = old_y = y0 = yt;
  z0 = ztime;

  Vx = Vy = Vxx = Vyy = Vxy = 0.0;
  old_vx = old_vy = 0.0;
  
  unitm = 1.0;
  Gamma = 1.0;

  Weight = p_domain()->p_Mesh()->ProfileLongi(x0,y0,z0)
          *p_domain()->p_Mesh()->ProfileTrans(x0,y0,z0);
          
  mass = 1.0/(TpCellx*TpCelly);

  sx=sxx;
  sy=syy;

  p_domain()->p_Mesh()->AddTrajectory(this);

  }; 



WDOUBLE Mesh::ProfileLongi(WDOUBLE xt, WDOUBLE yt, WDOUBLE zt)
{
    WDOUBLE tmp;
 

    switch(PProfileL)
    {

      case 1:

        if (zt<PlasmaBegin || zt>PlasmaEnd) return 0.0;

        if (zt<PlateauBegin)
        {
          tmp = (zt-PlasmaBegin)/(PlateauBegin-PlasmaBegin);
          return tmp;
        }

        if (zt>PlateauEnd) 
        {
          tmp = (PlasmaEnd - zt)/(PlasmaEnd-PlateauEnd);
          return tmp;
        
        }

      break;

      //
      case 2:

      break;




    }



    return 1.0;
}; 


WDOUBLE Mesh::ProfileTrans(WDOUBLE xt, WDOUBLE yt, WDOUBLE zt)
{


    switch(PProfileT)
    {

        case 1:
        if (xt*xt+yt*yt>PlasRadius*PlasRadius) return 0.0;
        break;

    }

    return 1.0;
}; 




//---Cell::~Cell--------------------------------------------->
Trajectory::~Trajectory()
{;};
