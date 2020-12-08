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

	double x;
	double y;
	double z;

	double x0;
	double y0;
	double z0;

	double Ex0;
	double Ey0;
	double Ez0;

	int type;

	union
	{
		double P_Source[1];
		double unitm;
	};
	double px, py, pz, gamma, q2m, weight, Q;

	double Wxw, Wyw, Wzw, Wxl, Wyl, Wzl;

  

	Domain *p_domain() {return Domain::p_Domain;}; 

	Particle *p_PrevPart;
	Particle *p_NextPart;


public:

	Particle(double x0p,  double y0p,  double z0p,
			 double pxp,  double pyp,  double pzp,
  			 double Ex0p, double Ey0p, double Ez0p,
  			 double q2mp, double weightp);

  virtual ~Particle() {;};

};



class Electron : public Particle {
	public:
	Electron(double x0p,  double y0p,  double z0p,
			 double pxp,  double pyp,  double pzp,
  			 double Ex0p, double Ey0p, double Ez0p,
  			 double q2mp, double weightp);
	~Electron();

};



class Ion : public Particle {
	public:
	Ion(double x0p,  double y0p,  double z0p,
		double pxp,  double pyp,  double pzp,
  		double Ex0p, double Ey0p, double Ez0p,
  		double q2mp, double weightp);
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
	double density;

	int PpCellx;        
	int PpCelly;
	int PpCellz;   

	int Seed_type;

	double P_Centerx;
	double P_Centery;
	double P_Centerz;

	double P_Sizex;
	double P_Sizey;
	double P_Sizez;

	double P_deltaZ;


	double P_px0;
	double P_py0;
	double P_pz0;

	double pxspread;
	double pyspread;
	double pzspread;

	double p_q2m;

public:

	double Density(double x0, double y0, double z0);
	int Get_NParts();

	Specie(char *name, FILE *f);
	~Specie();



};






#endif