using namespace std;

//���������
double deltt=2.5;
//?
double infinity=1e8;
double eps=1e-8;
double eps1=1e-2;
double pokmax=16;
double epstermc=1e-2;
double eps2c=0.5e-2;

//������ �� integer
double x0c[18] = {0,2,0, 0,4,0, 0,6,0, 0,0,0, 0,0,0, 0,0,0};
//������ �� integer
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
int arr[] = {0,1,2, 3,4,5, 6,7,8};
vector<int> Pnumc(arr, arr+9);
int arr1[] = {9,10,11, 12,13,14, 15,16,17};
vector<int> Rnumc(arr1, arr1+9);
double prepc[2][5][2] = {
	{{-8,1},{-20,1},{-20,1},{-8,-1},{-8,1}},
		{{20,1},{8,1},{8,-1},{20,-1},{20,1}}
};

//����� ���������
double tfc = 8;
//��� ��������������
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

//����������
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
//�������� ��������� �������� � �������
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

//**************************************************************
void SetEAixmax(vector<int> EAix1) {
int i;
bool flag;
sumixmax=0;
for (i=0; i<(ny-1); i++) {
	EAixmax[i]=EAix1[i];
}
for (i=0; i<(ny-1); i++) {
	EAix[i]=0;
}

while(not flag) {
	sumixmax++;
	for(i=0; i<(ny-1); i++) {
		sumixmax=sumixmax+EAix[i];
	}
	LexPM(EAix,flag);
}

}

//*************************************************************
void SetPsi(vector<vector <vector <int> > > Psi1) {
int i,j,k;
for(k=0; k<(kL-1); k++)
{
	for (i=0; i<(L-1); i++) {
		for (j=0; j<(L-1); j++) {
			Psi[k,i,j]=Psi1[k,i,j];
		}
	}
}	
}

//*************************************************************
//��������� ������ ���� � ����
void SetPsiBas(vector<vector <vector <int> > > Psi1) {
int i,j,k;
for(k=0; k<(kL-1); k++)
{
	for (i=0; i<(L-1); i++) {
		for (j=0; j<(L-1); j++) {
			Psi[k,i,j]=Psi1[k,i,j];
		}
	}
}		  
}
//*************************************************************
void VectortoGrey(vector<int> y) {
int x,i,j,k;
double r,g1;
g1=1;
for(i=0; i<(c-1); i++) { g1=g1*2; }
for(i=0; i<(p-1); i++) { q[i]=(q[i]-qmin[i])*g1/(qmax[i]-qmin[i]); }
for(i=0; i<(p*(c+d)-1); i++) { zb[i]=0; }
for(j=0; j<(p-1); j++) {
	x=trunc(q[j]);
    r=q[j]-x;
    k=c+j*(c+d)-1;
	while(k>=j*(c+d)) {
		zb[k]=x%2;
        x=x/2;
        k=k-1;
	}
	k=c+j*(c+d);
	while(k<(c+d)*(j+1)) {
		r=2*r;
        x=floor(r);
        zb[k]=x;
        r=r-x;
        k=k+1;
	}
	y[j*(c+d)]=zb[j*(c+d)];
	for(i=(j*(c+d)+1); i<((j+1)*(c+d)-1); i++) { y[i]=zb[i] xor zb[i-1]; }
}
}
//*************************************************************
// ��������� ������������ ��������
bool TestSource(int j) {
// ���� j-����� ���� ���������, �� ���������� false
int i,flag;
flag=true;
i=0;

while(i<=(*max_element(Pnum.begin(), Pnum.end())) && (j!=Pnum[i])) { i++; }

if (i<=*max_element(Rnum.begin(), Rnum.end())) { flag=false; } else {
	i=0;
	while((i<=*max_element(Rnum.begin(), Rnum.end())) and (j!=Rnum[i])) { i++; }
	if (i<=*max_element(Rnum.begin(), Rnum.end())) {flag=false;}
}
return flag;	
}


void GenVar(vector<int> w[5]) {
w[0].push_back(rand()%kL+1); //����� ����
w[1].push_back(rand()%4+1);
int w2 = (*w)[2];
switch((*w)[1]) {
	case 0:
	w[2].push_back(rand()%(L-1)+1);
	w[3].push_back(rand()%(L-w2-1)+1+w2+1);
	w[4].push_back(O1s[rand()%kw+1]);
	case '2':
	w[2].push_back(rand()%(L-1)+1);
	w[3].push_back(rand()%(L-w2-1)+1+w2+1);
	w[4].push_back(O1s[rand()%kw+1]);
	case '3':
	w[2].push_back(rand()%(L-1)+1);
	w[3].push_back(rand()%(L-w2-1)+1+w2+1);
	w[4].push_back(O1s[rand()%kw+1]);
	case '1':
	w[2].push_back(rand()%L+1);
	w[3].push_back((*w)[2]);
	w[4].push_back(O2s[rand()%kv+1]);
 }
}

//*************************************************************
void Variations(vector<vector <vector <vector <vector <int> > > > > w) {
// ������������ ��������
// 0 - ������ ��������������� ��������
// 1 - ������ ������������� ��������
// 2 - ���������� ����
// 3 - �������� ����
int i,j,s1,s2;

if (w[1]!=0) && (w[2]!=0) && (w[3]!=0) {
	switch ((*w)[1]) {
		case 0: if (Psi[w[0],w[2],w[3]]!=0) { Psi[w[0],w[2],w[3]]=w[4]; }
		case 1: if (Psi[w[0],w[2],w[2]]!=0) { Psi[w[0],w[2],w[3]]=w[4]; }
		case 2: if (Psi[w[0],w[2],w[3]]==0) { if (Psi[w[0],w[3],w[3]]!=0) {Psi[w[0],w[2],w[3]]=w[4];} }
		case 3: 
		s1=0;
		for (i=0; i<w[3]-1; i++) { if (Psi[w[0],i,w[3]]!=0) { s1++; } }
		s2=0;
		for (j=w[2]+1; j<(L-1); j++) { if (Psi[w[0],w[2],j]!=0) { s2++; } }
		if (s1>1) { if (s2>1) { Psi[w[0],w[1],w[2]]=0; } }
	}
}
}
//*************************************************************
void SetV_Entr() {
int i;
for (i=0; i< kP-1; i++) { V_Entr[Pnum[i]]=Vs[i]; }	
for (i=0; i< kR-1; i++) { V_Entr[Rnum[i]]=Cs[i]; }
}
//*************************************************************
double Ro_1(double z) {	return z; }
//*************************************************************
double Ro_2(double z) {	if (abs(z)>sqrt(infinity)) { return infinity; } }
//*************************************************************
double Ro_3(double z) {	return -z; }
//*************************************************************
double Ro_10(double z) { if (z>=0) { return 1; } else { return -1; } }
//*************************************************************
double Ro_4(double z) { return Ro_10(z)*sqrt(abs(z)); }
//*************************************************************
double Ro_5(double z) { if(abs(z)>eps) { return 1/z; } else { return Ro_10(z)/eps; } }
//*************************************************************
double Ro_6(double z) { if (z>-ln(eps)) { return -ln(eps); } else { return exp(z); } }
//*************************************************************
double Ro_7(double z) { if (abs(z)<exp(-pokmax)) { return ln(eps); } else { return ln(abs(z)); } }
//*************************************************************
double Ro_8(double z) { if (abs(z)>-ln(eps)) { return Ro_10(z); } else { return (1-exp(-z))/(1+exp(-z)); } }
//*************************************************************
double Ro_9(double z) { if (z>=0) {return 1; } else { return 0; } }
//*************************************************************
double Ro_11(double z) { return cos(z); }
//*************************************************************
double Ro_12(double z) { return sin(z); }
//*************************************************************
double Ro_13(double z) { return arctan(z); }
//*************************************************************
double Ro_15(double z) { if (abs(z)<eps) { return Ro_10(z)*eps; } else { return Ro_10(z)*exp(ln(abs(z))/3); } }
//*************************************************************
double Ro_14(double z) { if (abs(z)>Ro_15(infinity)) { return Ro_10(z)*infinity; } else { return sqr(z)*z; } }
//*************************************************************
double Ro_16(double z) { if (abs(z)<1) {return z;} else { return Ro_10(z); } }
//*************************************************************
double Ro_17(double z) { return Ro_10(z)*ln(abs(z)+1); }
//*************************************************************
double Ro_18(double z) { if (abs(z)>-ln(eps)) { return Ro_10(z)*infinity; } else { return Ro_10(z)*(exp(abs(z))-1); } }
//*************************************************************
double Ro_19(double z) { if (abs(z)>1/eps) { return Ro_10(z)*eps; } else { return Ro_10(z)*exp(-abs(z)); } }
//*************************************************************
double Ro_20(double z) { return z/2; }
//*************************************************************
double Ro_21(double z) { return 2*z; }
//*************************************************************
double Ro_22(double z) { if (abs(z)>-ln(eps)) { return 0; } else { return exp(abs(z)); } }
//*************************************************************
double Ro_23(double z) { if (abs(z)>1/eps) { return -Ro_10(z)/eps; } else { return z-z*sqr(z); } }
//*************************************************************
double Ro_24(double z) { if (z>-ln(eps)) { return eps/(1+eps); } else { return 1/(1+exp(-z)); } }
//*************************************************************
double Ro_25(double z) { if (z>0) { return 1; } else { return 0; } }
//*************************************************************
double Ro_26(double z) { if (abs(z)<eps1) { return 0; } else { return Ro_10(z); } }
//*************************************************************
double Ro_27(double z) { if (abs(z)>1) { return Ro_10(z); } else { return Ro_10(z)*(1-sqrt(1-sqr(z))); } }
//*************************************************************
double Ro_28(double z) { if (z*z>ln(infinity)) { return z*(1-eps); } else { return z*(1-exp(-sqr(z))); } }
//*************************************************************
double Xi_0(double z1, double z2) { return z1+z2; }
//*************************************************************
double Xi_1(double z1, double z2) { return z1*z2; }
//*************************************************************
double Xi_2 (double z1, double z2) { if (z1>=z2) { return z1; } else { return z2; } }
//*************************************************************
double Xi_3(double z1, double z2) { if (z1<z2) { return z1; } else { return z2; } }
//*************************************************************
double Xi_4(double z1, double z2) { return z1+z2-z1*z2; }
//*************************************************************
double Xi_5(double z1, double z2) { return Ro_10(z1+z2)*sqrt(sqr(z1)+sqr(z2)); }
//*************************************************************
double Xi_6(double z1, double z2) { return Ro_10(z1+z2)*(abs(z1)+abs(z2)); }
//*************************************************************
double Xi_7(double z1, double z2) { return Ro_10(z1+z2)*Xi_2(abs(z1),abs(z2)); }
//*************************************************************
void RPControl() {
int k,i,j;
double zz;
SetV_Entr();
for (k=0; k<(kL-1); k++) {
	for (i=0; i<(L-1); i++) {
		switch (psi[k,i,i]) {
			case 2:  z[k,i]=1; 
			case 3:  z[k,i]=-infinity; 
			case 4:  z[k,i]:=infinity; 
			default: z[k,i]:=0;
		}
	}
	for (i=0; i<kInP-1; i++) {
		if (M_Entr[k,2*i]=0) { z[k,i]=V_Entr[M_Entr[k,2*i+1]]; } else { z[k,i]:=z[M_Entr[k,2*i]-1,M_Entr[k,2*i+1]]; }
	}
	for (i=0; i<(L-2); i++) {
		for (j=i+1; j<(L-1); j++) {
			if (Psi[k,i,j]!=0) { 
				if (Psi[k,j,j]!=0) {
					switch (Psi[k,i,j]) {
						case 1: zz=Ro_1(z[k,i]);
						case 2: zz=Ro_2(z[k,i]);
						case 3: zz=Ro_3(z[k,i]);
						case 4: zz=Ro_4(z[k,i]);
						case 5: zz=Ro_5(z[k,i]);
						case 6: zz=Ro_6(z[k,i]);
						case 7: zz=Ro_7(z[k,i]);
						case 8: zz=Ro_8(z[k,i]);
						case 9: zz=Ro_9(z[k,i]);
						case 10: zz=Ro_10(z[k,i]);
						case 11: zz=Ro_11(z[k,i]);
						case 12: zz=Ro_12(z[k,i]);
						case 13: zz=Ro_13(z[k,i]);
						case 14: zz=Ro_14(z[k,i]);
						case 15: zz=Ro_15(z[k,i]);
						case 16: zz=Ro_16(z[k,i]);
						case 17: zz=Ro_17(z[k,i]);
						case 18: zz=Ro_18(z[k,i]);
						case 19: zz=Ro_19(z[k,i]);
						case 20: zz=Ro_20(z[k,i]);
						case 21: zz=Ro_21(z[k,i]);
						case 22: zz=Ro_22(z[k,i]);
						case 23: zz=Ro_23(z[k,i]);
						case 24: zz=Ro_24(z[k,i]);
						case 25: zz=Ro_25(z[k,i]);
						case 26: zz=Ro_26(z[k,i]);
						case 27: zz=Ro_27(z[k,i]); 
						case 28: zz=Ro_28(z[k,i]);
					}
					switch (Psi[k,j,j]) {
						case 1: z[k,j]=Xi_0(z[k,j],zz);
						case 2: z[k,j]=Xi_1(z[k,j],zz);
						case 3: z[k,j]=Xi_2(z[k,j],zz);
						case 4: z[k,j]=Xi_3(z[k,j],zz);
						case 5: z[k,j]=Xi_4(z[k,j],zz);
						case 6: z[k,j]=Xi_5(z[k,j],zz);
						case 7: z[k,j]=Xi_6(z[k,j],zz);
						case 8: z[k,j]=Xi_7(z[k,j],zz);
					}
				}
			}
		}
	} 
}
}
