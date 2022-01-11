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




Particle::Particle(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
			 WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
  			 WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
  			 WDOUBLE q2mp, WDOUBLE weightp)
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



Electron::Electron(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
				   WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
  			 	   WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
  			 	   WDOUBLE q2mp, WDOUBLE weightp)
		:Particle(x0p, y0p, z0p, pxp, pyp, pzp, Ex0p,Ey0p,Ez0p,q2mp,weightp)
{
	type = ELECTRON;
	Q    = 1.0;
}
Electron::~Electron()
{;};

Ion::Ion(WDOUBLE x0p,  WDOUBLE y0p,  WDOUBLE z0p,
		 WDOUBLE pxp,  WDOUBLE pyp,  WDOUBLE pzp,
		 WDOUBLE Ex0p, WDOUBLE Ey0p, WDOUBLE Ez0p,
		 WDOUBLE q2mp, WDOUBLE weightp)
		:Particle(x0p, y0p, z0p, pxp, pyp, pzp, Ex0p,Ey0p,Ez0p,q2mp,weightp)
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

	AddEntry("Q2M", &p_q2m, 1.0);

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

	AddEntry("Order",  &P_order, 2);
	AddEntry("deltaZ",  &P_deltaZ, 0.5);

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

WDOUBLE Specie::Density(WDOUBLE x0, WDOUBLE y0, WDOUBLE z0)
{

	WDOUBLE dentemp;
	WDOUBLE arg;
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

		// transversely rectangle and longitudinally trapezoid
		case 3:

			if( abs(x0-P_Centerx)*2<P_Sizex && abs(y0-P_Centery)*2<P_Sizey && abs(z0-P_Centerz)*2<P_Sizez )
			{
				return ((P_deltaZ-1)/P_Sizez*(z0-P_Centerz)*2+(P_deltaZ+1))/2.0;
			}

			return 0;
		break;


		// Transvesely gaussian and longitudinally high-order poly;
		case 4:
            arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
            arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;
            WDOUBLE arg2=0.0;

            if(z0-P_Centerz<0) return 0;

            arg2=1.0/(pow((z0-P_Centerz)/P_Sizez,P_order)+1);

            if(arg<4.0&&arg2>0.01)
            {
                    return exp(-arg)*arg2;
            }

            return 0;

        break;


        // Transves Donut-like
        case 5:
        	WDOUBLE r0=sqrt(x0*x0+y0*y0);
            arg += (r0-P_Centerx)*(r0-P_Centerx)/P_Sizex/P_Sizex;
            arg += (z0-P_Centerz)*(z0-P_Centerz)/P_Sizez/P_Sizez;
      		
            if(arg<4.0)
			{
				return exp(-arg);
			}

			return 0;

        break;


         // Half Gaussian ellipsoid
        case 7:
        	arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
        	arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;
        	arg += (z0-P_Centerz)*(z0-P_Centerz)/P_Sizez/P_Sizez;
        	if(z0-P_Centerz<0) return 0;
        	if(arg<16)
        	{
            	return exp(-arg);
        	}

        	return 0;
		break;

		// transversely Gaussian and longitudinally trapezoid
		case 8:

			arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
			arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;


			if( arg<4.0&&abs(z0-P_Centerz)*2<P_Sizez )
			{
				return exp(-arg)*( (P_deltaZ-1)/P_Sizez*(z0-P_Centerz)*2 + (P_deltaZ+1) )/2.0;
			}

			return 0;
		break;

		// transversely Gaussian and x
		case 9:

			arg += (x0-P_Centerx)*(x0-P_Centerx)/P_Sizex/P_Sizex;
			arg += (y0-P_Centery)*(y0-P_Centery)/P_Sizey/P_Sizey;

			if(z0-P_Centerz<0)
				arg += (z0-P_Centerz)*(z0-P_Centerz)/P_Sizez/P_Sizez;
			else
				arg += (z0-P_Centerz)*(z0-P_Centerz)/P_Sizez/P_Sizez/P_deltaZ/P_deltaZ;


			if( arg<4.0)
			{
				return exp(-arg);
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

	WDOUBLE dxp = dx/(specie->PpCellx);
	WDOUBLE dyp = dy/(specie->PpCelly);
	WDOUBLE dzp = dz/(specie->PpCellz);

	WDOUBLE wtemp;

	int GridXp = GridX*(specie->PpCellx);
	int GridYp = GridY*(specie->PpCelly);
	int GridZp = GridZ*(specie->PpCellz);

	WDOUBLE x0,y0,z0;
	WDOUBLE px,  py,  pz;
  	WDOUBLE Ex0, Ey0, Ez0;
  	WDOUBLE q2m, weight, Den;

	srand(time(NULL));

	Particle *p =NULL;

	// char name[128];
 // 	sprintf(name,"Parts_%d_.dg",Rank);
	// FILE * dFile;
	// dFile = fopen (name,"w");


	// loop cells
	for(k=0; k<GridZp; k++)
	{
		for(j=1; j<=GridY; j++)
		{
			for(i=1; i<=GridX; i++)
			{

				// seed particles in cell
				for(int sj=0;sj<specie->PpCelly;sj++)
				{
					for(int si=0;si<specie->PpCellx;si++)
					{

						Cell &c = GetCell(i, j, 0); 

						dxp = c.dx/(specie->PpCellx);
						dyp = c.dy/(specie->PpCelly);

						switch(S_type)
						{

							case 0:
								x0 = c.Xcord-c.dx*0.5 + WDOUBLE(si + 0.5)*dxp;
								y0 = c.Ycord-c.dy*0.5 + WDOUBLE(sj + 0.5)*dyp;
								z0 = WDOUBLE(k + 0.5)*dzp;

								
							break;

							case 1:
								WDOUBLE r1 = ((WDOUBLE) rand() / (RAND_MAX));
								WDOUBLE r2 = ((WDOUBLE) rand() / (RAND_MAX));
								WDOUBLE r3 = ((WDOUBLE) rand() / (RAND_MAX));

								x0 = c.Xcord-c.dx*0.5 + (si + r1)*dxp;
								y0 = c.Ycord-c.dy*0.5 + (sj + r2)*dyp;
								z0 = 0		  + (k + r3)*dzp;
							break;

						}

						Den = (specie->density)*specie->Density(x0, y0, z0);
						
						// wtemp=c.dx*c.dy/(specie->PpCellx)/(specie->PpCelly)/(specie->PpCellz);

						if(Den>0) 
						{

							// fprintf(dFile, "%f,%f,%f\n", x0,y0,z0);

							switch(P_type)
							{
								case ELECTRON:
									q2m = specie->p_q2m;
									weight = Den/(specie->PpCellz);
									px = specie->P_px0 + rand_gaussian (specie->pxspread);
									py = specie->P_py0 + rand_gaussian (specie->pyspread);
									pz = specie->P_pz0 + rand_gaussian (specie->pzspread);
									Ex0 = Ey0 = Ez0 = 0.0;
									p = new Electron(x0, y0, z0, px, py, pz, Ex0, Ey0, Ez0, q2m, weight);
									p->idx_i=i;
									p->idx_j=j;
									p->sx=dxp;
									p->sy=dyp;
								break;

								case ION:
									q2m = specie->p_q2m;
									weight = Den/(specie->PpCellz);
									px = specie->P_px0 + rand_gaussian (specie->pxspread);
									py = specie->P_py0 + rand_gaussian (specie->pyspread);
									pz = specie->P_pz0 + rand_gaussian (specie->pzspread);
									Ex0 = Ey0 = Ez0 = 0.0;
									p = new Ion(x0, y0, z0, px, py, pz, Ex0, Ey0, Ez0, q2m, weight);
									p->idx_i=i;
									p->idx_j=j;
									p->sx=dxp;
									p->sy=dyp;
								break;

							}

						}

					}//si

				}//sj

				

			}
		}

	}


	// fclose(dFile);

	return;

}



void Mesh::BeamSource()//v
{

	Particle *p = NULL;

	WDOUBLE ddx;
	WDOUBLE ddy;
	WDOUBLE ddz;

	WDOUBLE sx,sy;

	WDOUBLE xt;
	WDOUBLE yt;
	WDOUBLE zt;

	WDOUBLE q2m, massweig, Q;

	Cell  **c = new Cell*[18];
	WDOUBLE weight[18];

	int i,j,k;

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

		ddz = zt;
		i = p->idx_i;
		j = p->idx_j;
		k = floor(ddz/dz);

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(i < 1 || i > GridX || j < 1 || j > GridY || k < 0 || k > GridZ-2 )
		{
			p = p->p_PrevPart;
			continue;
		}

		c[0] = &GetCell(i-1,j-1,k); c[9] =  &GetCell(i-1,j-1,k+1);
		c[1] = &GetCell(i-1,j  ,k); c[10] = &GetCell(i-1,j  ,k+1);
		c[2] = &GetCell(i-1,j+1,k); c[11] = &GetCell(i-1,j+1,k+1);

		c[3] = &GetCell(i,  j-1,k); c[12] = &GetCell(i,  j-1,k+1);
		c[4] = &GetCell(i,  j  ,k); c[13] = &GetCell(i,  j  ,k+1);
		c[5] = &GetCell(i,  j+1,k); c[14] = &GetCell(i,  j+1,k+1);
		
		c[6] = &GetCell(i+1,j-1,k); c[15] = &GetCell(i+1,j-1,k+1);
		c[7] = &GetCell(i+1,j  ,k); c[16] = &GetCell(i+1,j  ,k+1);
		c[8] = &GetCell(i+1,j+1,k); c[17] = &GetCell(i+1,j+1,k+1);

		Cell &ccc= *c[4];

		ddx=ccc.dx;
		ddy=ccc.dy;

		sx=ddx;  //- re-size
		sy=ddy;  //- re-size

		massweig *= (p->sx)*(p->sy)/sx/sy; // re-weight

		WDOUBLE deltaxm=std::max(sx*0.5-(ddx*0.5+xt-ccc.Xcord),0.0);
		WDOUBLE deltaym=std::max(sy*0.5-(ddy*0.5+yt-ccc.Ycord),0.0);

		WDOUBLE deltaxp=std::max(sx*0.5-(ddx*0.5-xt+ccc.Xcord),0.0);
		WDOUBLE deltayp=std::max(sy*0.5-(ddy*0.5-yt+ccc.Ycord),0.0);

		WDOUBLE deltaxc=sx-deltaxm-deltaxp;
		WDOUBLE deltayc=sy-deltaym-deltayp;

		WDOUBLE deltaz=(ddz/dz-k);
		
		weight[0] = deltaxm*deltaym/c[0]->dx/c[0]->dy;
		weight[1] = deltaxm*deltayc/c[1]->dx/c[1]->dy;
		weight[2] = deltaxm*deltayp/c[2]->dx/c[2]->dy;

		weight[3] = deltaxc*deltaym/c[3]->dx/c[3]->dy;
		weight[4] = deltaxc*deltayc/c[4]->dx/c[4]->dy;
		weight[5] = deltaxc*deltayp/c[5]->dx/c[5]->dy;

		weight[6] = deltaxp*deltaym/c[6]->dx/c[6]->dy;
		weight[7] = deltaxp*deltayc/c[7]->dx/c[7]->dy;
		weight[8] = deltaxp*deltayp/c[8]->dx/c[8]->dy;

		for(int n=0;n<9;n++)
		{
			weight[n+9]=weight[n]*deltaz;
			weight[n] *=(1-deltaz);
		}
		for (int n=0; n<18; n++)
		{
			c[n] -> W_Source[6]  += massweig * weight[n] * Q; //density
			c[n] -> W_Source[7]  += massweig * weight[n] * Q * (p->px)/(p->gamma); //jx
			c[n] -> W_Source[8]  += massweig * weight[n] * Q * (p->py)/(p->gamma); //jy
			c[n] -> W_Source[9]  += massweig * weight[n] * Q * (p->pz)/(p->gamma); //jz
			c[n] -> W_Source[10] += massweig * weight[n] * abs(Q*q2m) /(p->gamma); //chi

		}
		p = p->p_PrevPart;

	}

	return;

}

//deprecated//May-27-tianhong
void Mesh::AdjustBSource()
{
	
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







