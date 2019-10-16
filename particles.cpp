//----------------------------------------------------------------------------------||
//-------------------                 particles.cpp              -------------------||
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

#include "wand_PIC.h"
#include <time.h> 
double rand_gaussian (double sigma);




Particle::Particle(double x0p,  double y0p,  double z0p,
			 double pxp,  double pyp,  double pzp,
  			 double Ex0p, double Ey0p, double Ez0p,
  			 double q2mp, double weightp)
{
	x = x0 = x0p;	y = y0 = y0p;	z = z0 = z0p;
	px = pxp;		py = pyp;		pz = pzp;
	Ex0 = Ex0p;		Ey0 = Ey0p;		Ez0 = Ez0p;
	q2m = q2mp;		weight = weightp;
	unitm = 1.0;	gamma=sqrt(1+px*px+py*py+pz*pz);
	
	Wxw = Wyw = Wzw = Wxl = Wyl = Wzl=0.0;

	p_PrevPart = p_NextPart = NULL;

	p_domain()->p_Mesh()->AddParticle(this);
}



Electron::Electron(double x0p,  double y0p,  double z0p,
				   double pxp,  double pyp,  double pzp,
  			 	   double Ex0p, double Ey0p, double Ez0p,
  			 	   double q2mp, double weightp)
		:Particle(x0p, y0p, z0p, pxp, pyp, pzp, Ex0p,Ey0p,Ez0p, q2mp,weightp)
{
	type = ELECTRON;
	Q    = 1.0;
}
Electron::~Electron()
{;};

Ion::Ion(double x0p,  double y0p,  double z0p,
		 double pxp,  double pyp,  double pzp,
		 double Ex0p, double Ey0p, double Ez0p,
		 double q2mp, double weightp)
		:Particle(x0p, y0p, z0p, pxp, pyp, pzp, Ex0p,Ey0p,Ez0p, q2mp,weightp)
{
	type = ION;
	Q    = -1.0;
}
Ion::~Ion()
{;};




//============================================================================
//=============================  Species   ===================================
//============================================================================
Specie::Specie(char *name, FILE *f) : NList(name)
{

	AddEntry("Specie", &P_type, 0);

	AddEntry("Profile", &P_profile, 0);
	AddEntry("Density", &density, 1.0);

	AddEntry("Part_per_Cellx", &PpCellx, 1);
	AddEntry("Part_per_Celly", &PpCelly, 1);
	AddEntry("Part_per_Cellz", &PpCellz, 1);
	AddEntry("SeedMethod", 	   &Seed_type,0);

	AddEntry("X0", &P_Centerx, 0.);
	AddEntry("Y0", &P_Centery, 0.);
	AddEntry("Z0", &P_Centerz, 0.);

	AddEntry("Sizex", &P_Sizex, 1.);
	AddEntry("Sizey", &P_Sizey, 1.);
	AddEntry("Sizez", &P_Sizez, 1.);

	AddEntry("Px0", &P_px0, 0.0);
	AddEntry("Py0", &P_py0, 0.0);
	AddEntry("Pz0", &P_pz0, 0.0);

	AddEntry("PxSpread", &pxspread, 0.0);
	AddEntry("PySpread", &pyspread, 0.0);
	AddEntry("PzSpread", &pzspread, 0.0);


	if(f) 
	{
		rewind(f);
		read(f);
	}
}

Specie::~Specie()
{;};

double Specie::Density(double x0, double y0, double z0)
{

	double dentemp;
	double arg;
	arg=0.0;

	switch(P_profile)
	{

		// Gaussian ellipsoid
		case 0:
			arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
			arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;
			arg += (z0-P_Centerz)*(z0-P_Centerz)/P_Sizez/P_Sizez;

			if(arg<4.0)
			{
				return exp(-arg);
			}

			return 0;
		break;

		// Transvesely gaussian and longitudinally triangle;
		case 1:
			arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
			arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;

			if(arg<4.0 && abs(z0-P_Centerz)*2<P_Sizez)
			{
				return exp(-arg)*((z0-P_Centerz)/P_Sizez+0.5);
			}

			return 0;
		break;

		// Rectangle;
		case 2:

			if( abs(x0-P_Centerx)*2<P_Sizex && abs(y0-P_Centery)*2<P_Sizey && abs(z0-P_Centerz)*2<P_Sizez )
			{
				return 1;
			}

			return 0;
		break;


	}

	return 0.0;

}


void Mesh::SeedParticles(Specie *specie)
{

	int i,j,k;

	int S_type = specie->Seed_type;
	int P_type = specie->P_type;

	double dxp = dx/(specie->PpCellx);
	double dyp = dy/(specie->PpCelly);
	double dzp = dz/(specie->PpCellz);

	double wtemp=(specie->PpCellx)*(specie->PpCelly)*(specie->PpCellz);;

	int GridXp = GridX*(specie->PpCellx);
	int GridYp = GridY*(specie->PpCelly);
	int GridZp = GridZ*(specie->PpCellz);

	double x0,y0,z0;
	double px,  py,  pz;
  	double Ex0, Ey0, Ez0;
  	double q2m, weight, Den;

	srand(time(NULL));

	Particle *p =NULL;

	for(k=0; k<GridZp; k++)
	{
		for(j=0; j<GridYp; j++)
		{
			for(i=0; i<GridXp; i++)
			{

				switch(S_type)
				{

					case 0:
						x0 = Offset_X + double(i + 0.5)*dxp;
						y0 = Offset_Y + double(j + 0.5)*dyp;
						z0 = 0		  + double(k + 0.5)*dzp;
					break;

					case 1:
						double r1 = ((double) rand() / (RAND_MAX));
						double r2 = ((double) rand() / (RAND_MAX));
						double r3 = ((double) rand() / (RAND_MAX));

						x0 = Offset_X + (i + r1)*dxp;
						y0 = Offset_Y + (j + r2)*dyp;
						z0 = 0		  + (k + r3)*dzp;
					break;

				}

				Den = (specie->density)*specie->Density(x0, y0, z0);
				if(Den>0) 
				{
					switch(P_type)
					{
						case ELECTRON:
							q2m = 1.0;
							weight = Den/wtemp;
							px = specie->P_px0 + rand_gaussian (specie->pxspread);
							py = specie->P_py0 + rand_gaussian (specie->pyspread);
							pz = specie->P_pz0 + rand_gaussian (specie->pzspread);
							Ex0 = Ey0 = Ez0 = 0.0;

							p = new Electron(x0, y0, z0, px, py, pz,
											Ex0, Ey0, Ez0, q2m, weight);
						break;

						case ION:
							//p = new Ion();
						break;

					}

				}

			}
		}

	}


	return;

}



void Mesh::BeamSource()
{

	Particle *p = NULL;

	double ddx;
	double ddy;
	double ddz;

	double xt;
	double yt;
	double zt;

	double q2m, massweig, Q;

	double weight[8];
	Cell  *ccc = NULL;

	double dxdy = dx*dy;
	int i,j,k;
	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;

	SetBSourceZero();
	
	p = p_Particle;

	while (p)
	{
		xt =  p->x;
		yt =  p->y;
		zt =  p->z;
		
		q2m = p->q2m;
		Q   = p->Q;
		massweig = p->weight;

		ddx = xt-(Offset_X-dx*0.5);
		ddy = yt-(Offset_Y-dy*0.5);
		ddz = zt;

		// idex of the corner cell
		i = floor(ddx/dx);
		j = floor(ddy/dy);
		k = floor(ddz/dz);

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(i < 0 || i > GridX || j < 0 || j > GridY || k < 0 || k > GridZ-2 )
		{
			p = p->p_PrevPart;
			continue;
		}
		//==================================================

		weight[0] = (i+1-ddx/dx) *(j+1-ddy/dy) *(k+1-ddz/dz);
		weight[1] = (i+1-ddx/dx) *(ddy/dy-j)   *(k+1-ddz/dz);
		weight[2] = (ddx/dx-i)   *(j+1-ddy/dy) *(k+1-ddz/dz);
		weight[3] = (ddx/dx-i)   *(ddy/dy-j)   *(k+1-ddz/dz);

		weight[4] = (i+1-ddx/dx) *(j+1-ddy/dy) *(ddz/dz-k);
		weight[5] = (i+1-ddx/dx) *(ddy/dy-j)   *(ddz/dz-k);
		weight[6] = (ddx/dx-i)   *(j+1-ddy/dy) *(ddz/dz-k);
		weight[7] = (ddx/dx-i)   *(ddy/dy-j)   *(ddz/dz-k);
		
		for (int n=0; n<8; n++)
		{
			switch(n)
			{
				case 0:  ccc = &GetCell(i,   j,   k);   break;//cmmm
				case 1:  ccc = &GetCell(i,   j+1, k);   break;//cmpm
				case 2:  ccc = &GetCell(i+1, j,   k);   break;//cpmm
				case 3:  ccc = &GetCell(i+1, j+1, k);   break;//cppm
				case 4:  ccc = &GetCell(i,   j,   k+1); break;//cmmp
				case 5:  ccc = &GetCell(i,   j+1, k+1); break;//cmpp
				case 6:  ccc = &GetCell(i+1, j,   k+1); break;//cpmp
				case 7:  ccc = &GetCell(i+1, j+1, k+1); break;//cppp
			}

			ccc -> W_Source[6]  += massweig * weight[n] * Q; //density
			ccc -> W_Source[7]  += massweig * weight[n] * Q * (p->px)/(p->gamma); //jx
			ccc -> W_Source[8]  += massweig * weight[n] * Q * (p->py)/(p->gamma); //jy
			ccc -> W_Source[9]  += massweig * weight[n] * Q * (p->pz)/(p->gamma); //jz
			ccc -> W_Source[10] += massweig * weight[n] * abs(Q*q2m) /(p->gamma); //chi

		}
		p = p->p_PrevPart;

	}

	return;

}

void Mesh::SetBSourceZero()
{

	for(int k=0; k<GridZ; k++)
	{
		for(int j=0; j<=GridY+1; j++)
		{
			for(int i=0; i<=GridX+1; i++)
			{
				Cell &ccc = GetCell(i,j,k);
				for (int n=6; n<11; n++)
				{
					ccc.W_Source[n] =0.0;
				}
	
			}
		}	
	
	}

	return;
}



double rand_gaussian (double sigma)
{
	double x, y, r2;
	do
	{
		x = (2.*rand()-RAND_MAX)/RAND_MAX;
		y = (2.*rand()-RAND_MAX)/RAND_MAX;
		r2 = x * x + y * y;
	}
	while (r2 > 1.0 || r2 == 0);
	return sigma * y * sqrt (-2.0 * log (r2) / r2);
}



