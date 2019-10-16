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

  double x;
  double y;

  double x0;
  double y0;
  double z0;

  double mass;
  double Gamma;
  double Weight;

   union
   {
      double T_Source[1];
      double unitm;
   };
   double Vx, Vy, Vxx, Vyy, Vxy;

   double old_x,old_y, old_vx, old_vy;
   
//=================================
//=======  Need to Send    ========
//=======  x0,y0, x, y, z  ========
//======= Vx, Vy, mass, weight ====
//=================================

  Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.
  Trajectory *p_PrevTraj;
  Trajectory *p_NextTraj;

public:
  
  Trajectory(double xt, double yt, double ztime, int TpCellx, int TpCelly);
  ~Trajectory();

};


#endif