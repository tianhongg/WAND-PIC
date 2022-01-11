//-------------------                trajectory.h                -------------------||
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


#ifndef H_TRAJS
#define H_TRAJS

class Trajectory{

  friend class Domain; 
  friend class Mesh;
  friend class Commute;

private:

  WDOUBLE x;
  WDOUBLE y;

  WDOUBLE x0;
  WDOUBLE y0;
  WDOUBLE z0;

  WDOUBLE mass;
  WDOUBLE Gamma;
  WDOUBLE Weight;

  // cell index 
  int idx_i;
  int idx_j;
  //shape
  WDOUBLE sx; // shape sx;
  WDOUBLE sy; // shape sy;

   union
   {
      WDOUBLE T_Source[1];
      WDOUBLE unitm;
   };
   WDOUBLE Vx, Vy, Vxx, Vyy, Vxy;

   WDOUBLE old_x,old_y, old_vx, old_vy;
   
//=================================
//=======  Need to Send    ========
//=======  x0,y0, x, y, z  ========
//======= Vx, Vy, mass, weight ====
//=================================

  Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.
  Trajectory *p_PrevTraj;
  Trajectory *p_NextTraj;

public:

  Trajectory(WDOUBLE xt, WDOUBLE yt, WDOUBLE ztime, int TpCellx, int TpCelly, WDOUBLE sxx, WDOUBLE syy);
  ~Trajectory();

};


#endif