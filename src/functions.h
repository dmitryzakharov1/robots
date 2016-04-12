using namespace std;

//константы
double deltt=2.5;
//?
double infinity=1e8;
double eps=1e-8;
double eps1=1e-2;
double pokmax=16;
double epstermc=1e-2;
double eps2c=0.5e-2;

//почему не integer
double x0c[18] = {0,2,0, 0,4,0, 0,6,0, 0,0,0, 0,0,0, 0,0,0};
//почему не integer
double xfc[9] = {-5,0,0, 0,0,0, 5,0,0};
double uminc[6] = {-5,-1, -5,-1, -5,-1};
double umaxc[6] = {5,1, 5,1, 5,1};
int O1sc[20] = {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16, 17,18,19,23};
int O2sc[2] = {1,2};
double qc[9] = {1,1,1, 1,1,1, 1,1,1};
double qminc[9] = {-10,-10,-10, -10,-10,-10, -10,-10,-10};
double qmaxc[9] = {10,10,10, 10,10,10, 10,10,10};
double qyminc[3] = {-10,-10,-10};
double qymaxc[3] = {10,10,10};
double stepsqyc[3] = {20,20,20};
int Pnumc[9] = {0,1,2, 3,4,5, 6,7,8};
int Rnumc[9] = {9,10,11, 12,13,14, 15,16,17};
double prepc[2][5][2] = {
	{{-8,1},{-20,1},{-20,1},{-8,-1},{-8,1}},
		{{20,1},{8,1},{8,-1},{20,-1},{20,1}}
};

//время окончания
double tfc = 8;
//шаг интегрирования
double dtc = 0.01;

int PsiMultBasc[4][16][16] = {
	//1
	{{0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},
	{0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},
	{0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},
	{0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},
				
	{0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0},

    {0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1}},
		
	//2
	{{0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0}, 
    {0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0},

    {0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1}},
	
	//3
	{{0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},  
    {0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0},

    {0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1}},
	
	//4
	
    {{0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0}, 
    {0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,0},
    {0,0,0,0, 0,0,0,1, 0,0,1,0, 0,0,0,0},

    {0,0,0,0, 0,0,0,0, 1,0,0,1, 0,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0},

    {0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0},
    {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1}}
};

int M_Entrc[4][12] = {{0,0,0,9, 0,1,0,10, 0,2,0,11},
{0,3,0,12, 0,4,0,13, 0,5,0,14},
{0,6,0,15, 0,7,0,16, 0,8,0,17},
{1,14,1,15, 2,14,2,15, 3,14,3,15}};

int M_Outc[9][2] = {{1,15},
{1,13},
{4,15},
{2,15},
{2,13},
{4,14},
{3,15},
{3,13},
{4,13}};

//переменные
// Parameters of NOP
int kL=4; // number of network operators (layers) in MNOP
int kInp=6; // number of inputs in each layer
int L=16; // //dimension of network operator matrices
int kp=9; //size of the set of variables
int kr=9; //size of the set of parameters
int kw=20; //size of the set of unary operations
int kv=2; //size of the set of binary operations
int Mout=9; // number of outputs
vector<double> Vs; //set of variables
vector<double> Cs; //set of parameters
vector<int> O1s; //set of unary operations
vector<int> O2s; //set of binary operations
int kOut;// quantity of exit in each matrix
vector<vector <int> > M_Entr; // matrix of connects
vector<double> V_Entr; //vector of inputs
vector<int> Pnum; //vector of number of positions in input vector for variables
vector<int> Rnum; //vector of number of positions in input vector for parameters
vector<vector <int> > M_Out; //matrix of outputs
vector<vector <double> > z; //vector of nodes
vector<vector <string> > zs; //string for mathematical expression
vector<vector <vector <int> > > Psi,Psi0; //Network  operator matrices

// Parameters of GA
int HH=1024; //number of chromosomes in the initial population
int PP=64; // number of generations (loops of GA)
int RR=64; // number of couples in one generation
int lchr=8; //length of structural part of chromosome
int Epo=24;  // number of generations between epochs
int kel=8;  // number of elitaring chromosomes
double alfa=0.4;  // parameter to calculate probability of crossover
double pmut=0.7;  // probability of mutation
vector<vector <vector <int> > >PopChrStr[5]; //array for structural parts of chromosomes
vector<vector <vector <int> > >PopChrPar; //array for parametrical parts of chromosomes

int nfu=2;  // number of functionals
vector<vector <double> > Fuh; // values of functionals for each chromosome
vector<vector <double> > FuhNorm;
vector<double> Shtraf;
vector<int> Lh; // values distance to Pareto set
vector<int> Pareto; // Pareto set

vector<vector <int> > Son1s[5],Son2s[5],Son3s[5],Son4s[5]; //structural part of sons
vector<int> Son1p,Son2p,Son3p,Son4p; //parametrical part of sons
int L1,L2,L3,L4; //values distance to Pareto set for sons
vector<double> Fu1,Fu2,Fu3,Fu4; // values of functionals for sons

int p=9;    //number of serching parameters
int c=4;    //number of bits for integer part
int d=12;    //number of bits for fractional part
vector<double> q; //vector of parameters
vector<double> qmax,qmin; // vectors of constraints for parameters
vector<int> zb; //additional vector

// Undefined parameters
int ny=3;// dimention of vector of undefine parameters
vector<double> qy; //vextor of undefined parameters
vector<double> qymax,qymin;
vector<double> stepsqy;//vector of steps of undefine parameters

// Parameters of model
int kRob=3;// number of robots;
int n=18;//dimention of system;
int m=6;//dimention of control;
int ll=42;//dimention of viewing;

vector<double> x; // vector of state
vector<double> xs;
vector<double> x0; // vector of initial condition
vector<double> xf;
vector<double> fb;
vector<double> fa;
vector<double> u; // vector of control
vector<double> umin; // vector of control
vector<double> umax; // vector of control
int lv;
double t;
double dt=0.01; //step of integration;
double tf=8;//terminal time;
double dtp=0.1;//step of print;

int kChoose;// number of choose chromosome
double f0,f1;
vector<double> xf1;
vector<int> Prior;
vector<bool> FlagStop;
double epsterm=0.01;
int ix,iy;
vector<int> EAix,EAixmax;
double sumixmax;

int i,j,k, pt,rt,k1,k2,lmax,imax,ks1,ks2;
double ksi,su, su1, Fumax,Fumin;

 //************************  1.   *******************************
//ОПИСАНИЕ ВЛОЖЕННЫХ ПРОЦЕДУР И ФУНКЦИЙ
//**************************************************************
void LexPM(vector<int> EAix, bool flag) {
int i,j;
i=ny-1;
while((i>=0) && (EAix[i]==EAixmax[i])) { i--; }
if (i>=0) {
	EAix[i]=EAix[i]+1;
	for (j=i+1; j<(ny-1); j++) {
		EAix[j]=0;
	}
	flag = true;
} else {
	flag = false;
}
	
}