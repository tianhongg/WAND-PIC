//----------------------------------------------------------------------------------||
//-------------------                 particles.h                -------------------||
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
//---Starting---------           : Apr-23-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||



#ifndef H_PARTS
#define H_PARTS

#define ELECTRON 0
#define ION 1




class Particle{

	friend class Domain; 
	friend class Mesh;
	friend class Commute;

public:

	WDOUBLE x;
	WDOUBLE y;
	WDOUBLE z;

	WDOUBLE x0;
	WDOUBLE y0;
	WDOUBLE z0;

	WDOUBLE Ex0;
	WDOUBLE Ey0;
	WDOUBLE Ez0;

	int type;

	int idx_i;
	int idx_j;

	WDOUBLE sx;
	WDOUBLE sy;

	union
	{
		WDOUBLE P_Source[1];
		WDOUBLE unitm;
	};
	WDOUBLE px, py, pz, gamma, q2m, weight, Q;

	WDOUBLE Wxw, Wyw, Wzw, Wxl, Wyl, Wzl;

  

	Domain *p_domain() {return Domain::p_Domain;}; 

	Particle *p_PrevPart;
	Particle *p_NextPart;


public:

	Particle(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
			 WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
  			 WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
  			 WDOUBLE q2mp, WDOUBLE weightp);

  virtual ~Particle() {;};

};



class Electron : public Particle {
	public:
	Electron(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
			 WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
  			 WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
  			 WDOUBLE q2mp, WDOUBLE weightp);
	~Electron();

};



class Ion : public Particle {
	public:
	Ion(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
		WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
  		WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
  		WDOUBLE q2mp, WDOUBLE weightp);
	~Ion();

};




class Specie : public NList {

	friend class Domain;
	friend class Mesh;
	
public:

	Domain *p_domain() {return Domain::p_Domain;}; 


	// 0 - Electron; 1 - ion
	int P_type;

	int P_profile;
	WDOUBLE density;

	int PpCellx;        
	int PpCelly;
	int PpCellz;   

	int Seed_type;

	WDOUBLE P_Centerx;
	WDOUBLE P_Centery;
	WDOUBLE P_Centerz;

	WDOUBLE P_Sizex;
	WDOUBLE P_Sizey;
	WDOUBLE P_Sizez;


	int P_order;
	WDOUBLE P_deltaZ;


	WDOUBLE P_px0;
	WDOUBLE P_py0;
	WDOUBLE P_pz0;

	WDOUBLE pxspread;
	WDOUBLE pyspread;
	WDOUBLE pzspread;

	WDOUBLE p_q2m;

public:

	WDOUBLE Density(WDOUBLE x0, WDOUBLE y0, WDOUBLE z0);
	int Get_NParts();

	Specie(char *name, FILE *f);
	~Specie();



};






#endif