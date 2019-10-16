//----------------------------------------------------------------------------------||
//-------------------                multigrid.h                 -------------------||
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
//---Starting---------           : Feb-07-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||








#ifndef H_MULTI
#define H_MULTI




class MG_Cell
{

   friend class MultiGrid;
   friend class Commute;

private:

	// value at this cell
	union
	{	
		double M_value[1];
		double field;
	};
		double source;   //1
		double Residu;   //2
		double savefield; //3
		double Chi;  //4


	union
	{	
		dcomplex C_value[1];
		dcomplex C_field;
	};
		dcomplex C_source;   //1
		dcomplex C_Residu;   //2
		dcomplex C_savefield; //3
		dcomplex C_Chi;  //4



	//point to the four cells we need to restrict from
	MG_Cell *p_Res_cc;
	MG_Cell *p_Res_xm;
	MG_Cell *p_Res_xp;
	MG_Cell *p_Res_ym;
	MG_Cell *p_Res_yp;
	MG_Cell *p_Res_mm;
	MG_Cell *p_Res_mp;
	MG_Cell *p_Res_pm;
	MG_Cell *p_Res_pp;


	//point to the four cells we need to prolongate from
	MG_Cell *p_Pro_xm;
	MG_Cell *p_Pro_xp;
	MG_Cell *p_Pro_ym;
	MG_Cell *p_Pro_yp;

	int Protype;  // four types of prolongation from 0 to 3;
	int Grandidx_X;
	int Grandidx_Y;

public:


	MG_Cell(); 
	~MG_Cell()
	{;};// do not release any pointer, this class persists to the end of simulation anyway


};




class MultiGrid :public NList {

	friend class Mesh;
private:


	Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.


	static const int Layer_buf = 20;  //Layer Number Buffer


	int LayerGridX[Layer_buf];
	int LayerGridY[Layer_buf];

	int Layer_Sum[Layer_buf];
	int MeshAmplif[Layer_buf+10];

	int BLayerGrid[Layer_buf];
	int BLayer_Sum[Layer_buf];

	int MPI_Layer;  //MPI Level layers idx starts from one !!!!!!!
	int SER_Layer; 	//Single rank layers
	int u1;
	int u2;

	int RelaxType;
	int BottomType;

	double dxdy;
	int GridXY;
	double EpsLim;
	double omega;


	int BottomCells;  //total cells at the bottom layer without the bound cells.
	int Worker;
	int *Bottomdisp;
	int *BottomSend;
	int *recebottomX;
	int *recebottomY;


	int Rank;
	int Xpa, Ypa;

	MG_Cell*  p_MGCell;
	MG_Cell*  p_MGBCell;
	Mesh *p_Meshs;


public:


	inline int GetMG_N(int i, int j, int layer)
	{
		return Layer_Sum[layer]+(LayerGridX[layer]+2)*j+i;
	};



	inline MG_Cell& GetMGCell(int i, int j, int layer)
	{
		return p_MGCell[GetMG_N(i,j,layer)];
	};


	inline int GetMGB_N(int i, int j, int layer)
	{
		return BLayer_Sum[layer]+(BLayerGrid[layer]+2)*j+i;
	};

	inline MG_Cell& GetMGBCell(int i, int j, int layer)
	{
		return p_MGBCell[GetMGB_N(i,j,layer)];
	};

	inline int GetLayerGridX(int layer)
	{
		return LayerGridX[layer];
	};

	inline int GetLayerGridY(int layer)
	{
		return LayerGridY[layer];
	};


	//======== Double Type ==========================================
	void		Exchange(int what, int layer);
	void 	 Restriction(int send, int rece, int tolayer, int where);
	void 	RestrictionB(int send, int rece, int tolayer, int where);
	void 	Prolongation(int send, int rece, int tolayer, int where);
	void		 SetZero(int what, int layer, int where);
	void	  Relaxation(int field, int layer, int where);
	void		Residual(int field, int layer, int where);
	void	SendtoBottom(int what);
	void  BottomSendBack(int what);
	void   AddCorrection(int layer, int where);
	void  MG_BottomLayer(int field);
	int	  MG_V_cycle(int field, int k);
	void	  Put_Source(int field, int k);
	void	  Put_Fields(int field, int k);
	double	   FindError(double &maxall);

	//======== Complex Type ==========================================
	void 	 RestrictionC(int send, int rece, int tolayer, int where);
	void 	RestrictionBC(int send, int rece, int tolayer, int where);
	void 	ProlongationC(int send, int rece, int tolayer, int where);
	void	SendtoBottomC(int what);
	void  BottomSendBackC(int what);
	void  MG_BottomLayerC(int field);
	void	  RelaxationC(int field, int layer, int where);
	void		ResidualC(int field, int layer, int where);
	void		 SetZeroC(int what, int layer, int where);
	void		ExchangeC(int what, int layer);
	void   AddCorrectionC(int layer, int where);
	double	   FindErrorC(double &maxall);
	void	  Put_SourceC(int field, int k, int NF);
	void	  Put_FieldsC(int field, int k, int NF);
	int	  MG_V_cycleC(int field, int k, int NF);

	//====== debug ======================
	void	DebugWrite(int layer, int what, char *name, int where);


	MultiGrid(int rank, int  XGridN, int  YGridN, FILE *f);
	~MultiGrid();




};


#endif