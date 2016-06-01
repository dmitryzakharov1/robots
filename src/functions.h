using namespace std;

//êîíñòàíòû
double deltt=2.5;
//?
double infinity=1e8;
double eps=1e-8;
double eps1=1e-2;
double pokmax=16;
double epstermc=1e-2;
double eps2c=0.5e-2;
int Lbasc=2;
int kLc=4;
int kInpc=6;
int L1c=16;
int Mout1c=9;

//ïî÷åìó íå integer
double x0c[18] = {0,2,0, 0,4,0, 0,6,0, 0,0,0, 0,0,0, 0,0,0};
//ïî÷åìó íå integer
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

//âðåìÿ îêîí÷àíèÿ
double tfc = 8;
//øàã èíòåãðèðîâàíèÿ
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

//ïåðåìåííûå
// Parameters of NOP
int kL=4; // number of network operators (layers) in MNOP
int kInP=6; // number of inputs in each layer
int L=16; // //dimension of network operator matrices
int kP=9; //size of the set of variables
int kR=9; //size of the set of parameters
int kW=20; //size of the set of unary operations
int kV=2; //size of the set of binary operations
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
vector<vector<vector <int> > > Psi,Psi0; //Network  operator matrices

// Parameters of GA
int HH=102; //number of chromosomes in the initial population 1024
int PP=6; // number of generations (loops of GA) 64
int RR=6; // number of couples in one generation 64
int lchr=8; //length of structural part of chromosome
int Epo=24;  // number of generations between epochs
int kel=8;  // number of elitaring chromosomes
double alfa=0.4;  // parameter to calculate probability of crossover
double pmut=0.7;  // probability of mutation
vector<vector <vector <int> > >PopChrStr; //array for structural parts of chromosomes
vector<vector <int> > PopChrPar; //array for parametrical parts of chromosomes

int nfu=2;  // number of functionals
vector<vector <double> > Fuh; // values of functionals for each chromosome
vector<vector <double> > FuhNorm;
vector<double> Shtraf;
vector<int> Lh; // values distance to Pareto set
vector<int> Pareto; // Pareto set

vector<vector <int> > Son1s,Son2s,Son3s,Son4s; //structural part of sons
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
int lv=42;//dimention of viewing;

vector <double> y;// vector of vewing

vector<double> x; // vector of state
vector<double> xs;
vector<double> x0; // vector of initial condition
vector<double> xf;
vector<double> fb;
vector<double> fa;
vector<double> u; // vector of control
vector<double> umin; // vector of control
vector<double> umax; // vector of control

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
double ksi, Fumax,Fumin;
double suGA, su1GA;
vector <double> su, su1;
int kol;

 //************************  1.   *******************************
//    
//**************************************************************
void LexPM(vector<int> EAix, bool flag) {
int i,j;
i=ny-1;
while((i>=0) && (EAix[i]==EAixmax[i])) { i--; }
if (i>=0) {
	EAix[i]=EAix[i]+1;
	for (j=i+1; j<=(ny-1); j++) {
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
//ny=3
for (i=0; i<=(ny-1); i++) {
	EAixmax[i]=EAix1[i];
}
for (i=0; i<=(ny-1); i++) {
	EAix[i]=0;
}

do {
	sumixmax++;
	for(i=0; i<=(ny-1); i++) {
		sumixmax=sumixmax+EAix[i];
	}
	//cout << "Start LexPM(EAix,flag)" << endl;
	LexPM(EAix,flag);
	//cout << "End of LexPM(EAix,flag)" << endl;
}
while (flag);

}

//*************************************************************
void SetPsi(vector<vector <vector <int> > > Psi1) {
//int i,j,k;
//for(k=0; k<(kL-1); k++)
//{
//	for (i=0; i<(L-1); i++) {
//		for (j=0; j<(L-1); j++) {
//			Psi[k,i,j]=Psi1[k,i,j];
//		}
//	}
//}
Psi=Psi1;	
}

//*************************************************************
//ïðîöåäóðû äåëàþò îäíî è òîæå
void SetPsiBas(vector<vector <vector <int> > > Psi1) {
//int i,j,k;
//for(k=0; k<(kL-1)-1; k++)
//{
//	for (i=0; i<(L-1)-1; i++) {
//		for (j=0; j<(L-1)-1; j++) {
//			Psi[k,i,j]=Psi1[k,i,j];
//		}
//	}
//}
Psi=Psi1;		  
}
//*************************************************************
void VectortoGrey(vector<int> y) {
int x,i,j,k;
double r,g1;
g1=1;
for(i=0; i<=(c-1); i++) { g1=g1*2;  }
for(i=0; i<=(p-1); i++) { q[i]=(q[i]-qmin[i])*g1/(qmax[i]-qmin[i]); }
for(i=0; i<=(p*(c+d)-1); i++) { zb[i]=0; }



for(j=0; j<=(p-1); j++) {
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
	
	//for(int jj=0; jj<zb.size(); jj++) {
	//	cout << "zb[" << jj << "]=" << fixed << zb[jj] << endl;
	//	}
	
	y[j*(c+d)]=zb[j*(c+d)];
	
	for(i=(j*(c+d)+1); i<=((j+1)*(c+d)-1); i++) { y[i]=zb[i] xor zb[i-1]; }
}

		//for(int jj=0; jj<y.size(); jj++) {
		//cout << "y[" << jj << "]=" << fixed << y[jj] << endl;
		//}
		//getchar();

}
//*************************************************************
// Ãåíåðàöèÿ ýëåìåíòàðíîé îïåðàöèè
bool TestSource(int j) {
// åñëè j-íîìåð óçëà èñòî÷íèêà, òî âîçâðàùàåò false
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


void GenVar(vector<int> w) {
w.push_back(rand()%kL+1); //íîìåð ñëîÿ
w.push_back(rand()%4+1);
int w2 = w[2];
switch(w[1]) {
	case 0:
	w.push_back(rand()%(L-1)+1);
	w.push_back(rand()%(L-w2-1)+1+w2+1);
	w.push_back(O1s[rand()%kW+1]);
	case '2':
	w.push_back(rand()%(L-1)+1);
	w.push_back(rand()%(L-w2-1)+1+w2+1);
	w.push_back(O1s[rand()%kW+1]);
	case '3':
	w.push_back(rand()%(L-1)+1);
	w.push_back(rand()%(L-w2-1)+1+w2+1);
	w.push_back(O1s[rand()%kW+1]);
	case '1':
	w.push_back(rand()%L+1);
	w.push_back(w[2]);
	w.push_back(O2s[rand()%kV+1]);
 }
}

//*************************************************************
void Variations(vector<int> w) {
// Ýëåìåíòàðíûå îïåðàöèè
// 0 - çàìåíà íåäèàãîíàëüíîãî ýëåìåíòà
// 1 - çàìåíà äèàãîíàëüíîãî ýëåìåíòà
// 2 - äîáàâëåíèå äóãè
// 3 - óäàëåíèå äóãè
int i,j,s1,s2;

if ((w[1] != 0) || (w[2] != 0) || (w[3] != 0)) {
	switch (w[1]) {
		case 0: if (Psi[w[0]][w[2]][w[3]] != 0) { Psi[w[0]][w[2]][w[3]]=w[4]; }
		case 1: if (Psi[w[0]][w[2]][w[2]]!=0) { Psi[w[0]][w[2]][w[2]]=w[4]; }
		case 2: if (Psi[w[0]][w[2]][w[3]]==0) { if (Psi[w[0]][w[3]][w[3]]!=0) {Psi[w[0]][w[2]][w[3]]=w[4];} }
		case 3: 
		s1=0;
		for (i=0; i<=w[3]-1; i++) { if (Psi[w[0]][i][w[3]]!=0) { s1++; } }
		s2=0;
		for (j=w[2]+1; j<=(L-1); j++) { if (Psi[w[0]][w[2]][j]!=0) { s2++; } }
		if (s1>1) { if (s2>1) { Psi[w[0]][w[1]][w[2]]=0; } } 
	}
}
}
//*************************************************************
void SetV_Entr() {
int i;
for (i=0; i<=kP-1; i++) { V_Entr[Pnum[i]]=Vs[i]; }	
for (i=0; i<=kR-1; i++) { V_Entr[Rnum[i]]=Cs[i]; }
}
//*************************************************************
double Ro_1(double z) {	return z; }
//*************************************************************
double Ro_2(double z) {	if (fabs(z)>sqrt(numeric_limits<double>::infinity( ))) { return numeric_limits<double>::infinity( ); } }
//*************************************************************
double Ro_3(double z) {	return -z; }
//*************************************************************
double Ro_10(double z) { if (z>=0) { return 1; } else { return -1; } }
//*************************************************************
double Ro_4(double z) { return Ro_10(z)*sqrt(fabs(z)); }
//*************************************************************
double Ro_5(double z) { if(fabs(z)>eps) { return 1/z; } else { return Ro_10(z)/eps; } }
//*************************************************************
double Ro_6(double z) { if (z>-log(eps)) { return -log(eps); } else { return exp(z); } }
//*************************************************************
double Ro_7(double z) { if (fabs(z)<exp(-pokmax)) { return log(eps); } else { return log(fabs(z)); } }
//*************************************************************
double Ro_8(double z) { if (fabs(z)>-log(eps)) { return Ro_10(z); } else { return (1-exp(-z))/(1+exp(-z)); } }
//*************************************************************
double Ro_9(double z) { if (z>=0) {return 1; } else { return 0; } }
//*************************************************************
double Ro_11(double z) { return cos(z); }
//*************************************************************
double Ro_12(double z) { return sin(z); }
//*************************************************************
double Ro_13(double z) { return atan(z); }
//*************************************************************
double Ro_15(double z) { if (fabs(z)<eps) { return Ro_10(z)*eps; } else { return Ro_10(z)*exp(log(fabs(z))/3); } }
//*************************************************************
double Ro_14(double z) { if (fabs(z)>Ro_15(numeric_limits<double>::infinity( ))) { return Ro_10(z)*numeric_limits<double>::infinity( ); } else { return sqrt(z)*z; } }
//*************************************************************
double Ro_16(double z) { if (fabs(z)<1) {return z;} else { return Ro_10(z); } }
//*************************************************************
double Ro_17(double z) { return Ro_10(z)*log(fabs(z)+1); }
//*************************************************************
double Ro_18(double z) { if (fabs(z)>-log(eps)) { return Ro_10(z)*numeric_limits<double>::infinity( ); } else { return Ro_10(z)*(exp(fabs(z))-1); } }
//*************************************************************
double Ro_19(double z) { if (fabs(z)>1/eps) { return Ro_10(z)*eps; } else { return Ro_10(z)*exp(-fabs(z)); } }
//*************************************************************
double Ro_20(double z) { return z/2; }
//*************************************************************
double Ro_21(double z) { return 2*z; }
//*************************************************************
double Ro_22(double z) { if (fabs(z)>-log(eps)) { return 0; } else { return exp(fabs(z)); } }
//*************************************************************
double Ro_23(double z) { if (fabs(z)>1/eps) { return -Ro_10(z)/eps; } else { return z-z*sqrt(z); } }
//*************************************************************
double Ro_24(double z) { if (z>-log(eps)) { return eps/(1+eps); } else { return 1/(1+exp(-z)); } }
//*************************************************************
double Ro_25(double z) { if (z>0) { return 1; } else { return 0; } }
//*************************************************************
double Ro_26(double z) { if (fabs(z)<eps1) { return 0; } else { return Ro_10(z); } }
//*************************************************************
double Ro_27(double z) { if (fabs(z)>1) { return Ro_10(z); } else { return Ro_10(z)*(1-sqrt(1-sqrt(z))); } }
//*************************************************************
double Ro_28(double z) { if (z*z>log(numeric_limits<double>::infinity( ))) { return z*(1-eps); } else { return z*(1-exp(-sqrt(z))); } }
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
double Xi_5(double z1, double z2) { return Ro_10(z1+z2)*sqrt(sqrt(z1)+sqrt(z2)); }
//*************************************************************
double Xi_6(double z1, double z2) { return Ro_10(z1+z2)*(fabs(z1)+fabs(z2)); }
//*************************************************************
double Xi_7(double z1, double z2) { return Ro_10(z1+z2)*Xi_2(fabs(z1),fabs(z2)); }
//*************************************************************
void RPControl() {
int k,i,j;
double zz;
SetV_Entr();
//cout << "At RPControl z.size()=" << z.size() << " z[0].size()=" << z[0].size() << endl;
for (k=0; k<=kL-1; k++) {
	for (i=0; i<=L-1; i++) {
	 
		int temp;
		temp=Psi[k][i][i];
		//cout << "point" << endl;
		switch (temp) {
			case 2:   z[k][i]=1;  //cout << "Choosed 2a" << endl; //cout << "Choosed 2a z[" << k << "][" << i << "]=" << z[k][i] << endl; 
			case 3:   z[k][i]=-numeric_limits<double>::infinity( ); // cout << "Choosed 3a" << endl;
			case 4:   z[k][i]=numeric_limits<double>::infinity( );  //cout << "Choosed 4a" << endl;
			default:  z[k][i]=0; //cout << "Choosed default" << endl;
		}
		//cout <<  "k=" << k << " i=" << i << endl;
		//cout << "RPControl. Debug. Psi[" << k << "][" << i << "][" << i <<"] " << " temp=" << Psi[k][i][i] << " z[" << k << "][" << i << "]=" << z[k][i] << endl;
	}
	
		/*for (int ii=0; ii<=kL-1; ii++) {
			for (int jj=0; jj <= L-1; jj++) {
				cout << "z[" << ii << "][" << jj << "]=" << z[ii][jj] << endl;
			}
		}*/	 
	
	for (i=0; i<=kInP-1; i++) {
	//cout <<  "k=" << k << " i=" << i << endl;
		//if (M_Entr[k][2*i]==0) { z[k][i]=V_Entr[M_Entr[k][2*i+1]]; } else { z[k][i]=z[M_Entr[k][2*i]-1][M_Entr[k][2*i+1]]; }
		if (M_Entr[k][2*i]==0) { 
		z[k][i]=V_Entr[M_Entr[k][2*i+1]]; 
		} else {
		z[k][i]=z[M_Entr[k][2*i]-1][M_Entr[k][2*i+1]]; 
		}
	}

	//cout << "point" << endl;
	
	for (i=0; i<=(L-2); i++) {

		for (j=i+1; j<=(L-1); j++) {
			if (Psi[k][i][j]!=0) { 
				if (Psi[k][j][j]!=0) {
					switch (Psi[k][i][j]) {
						case 1: zz=Ro_1(z[k][i]);
						case 2: zz=Ro_2(z[k][i]);
						case 3: zz=Ro_3(z[k][i]);
						case 4: zz=Ro_4(z[k][i]);
						case 5: zz=Ro_5(z[k][i]);
						case 6: zz=Ro_6(z[k][i]);
						case 7: zz=Ro_7(z[k][i]);
						case 8: zz=Ro_8(z[k][i]);
						case 9: zz=Ro_9(z[k][i]);
						case 10: zz=Ro_10(z[k][i]);
						case 11: zz=Ro_11(z[k][i]);
						case 12: zz=Ro_12(z[k][i]);
						case 13: zz=Ro_13(z[k][i]);
						case 14: zz=Ro_14(z[k][i]);
						case 15: zz=Ro_15(z[k][i]);
						case 16: zz=Ro_16(z[k][i]);
						case 17: zz=Ro_17(z[k][i]);
						case 18: zz=Ro_18(z[k][i]);
						case 19: zz=Ro_19(z[k][i]);
						case 20: zz=Ro_20(z[k][i]);
						case 21: zz=Ro_21(z[k][i]);
						case 22: zz=Ro_22(z[k][i]);
						case 23: zz=Ro_23(z[k][i]);
						case 24: zz=Ro_24(z[k][i]);
						case 25: zz=Ro_25(z[k][i]);
						case 26: zz=Ro_26(z[k][i]);
						case 27: zz=Ro_27(z[k][i]); 
						case 28: zz=Ro_28(z[k][i]);
					}
					switch (Psi[k][j][j]) {
						case 1: z[k][j]=Xi_0(z[k][j],zz);
						case 2: z[k][j]=Xi_1(z[k][j],zz); //cout << "Choosed 2b z[" << k << "][" << j << "]=" << z[k][i] << endl;
						case 3: z[k][j]=Xi_2(z[k][j],zz);
						case 4: z[k][j]=Xi_3(z[k][j],zz);
						case 5: z[k][j]=Xi_4(z[k][j],zz);
						case 6: z[k][j]=Xi_5(z[k][j],zz);
						case 7: z[k][j]=Xi_6(z[k][j],zz);
						case 8: z[k][j]=Xi_7(z[k][j],zz);
					}
				}
			}
		}
	}
}	

}
//*************************************************************
int Rast(vector < double > Fu) {
    int i, j, k, count;
    count = 0;
	
	for(int jj=0; jj<Fu.size(); jj++) {
	cout << "Fu[" << jj << "]=" << Fu[jj] << endl;
	}
	
	for(int jj=0; jj<Fuh[jj].size(); jj++) {
	cout << "Fuh[0][" << jj << "]=" << Fuh[0][jj] << endl;
	}
	
    for (i = 0; i <=(HH - 1); i++) {
        j = 0;
        while ((j < nfu) && (Fu[j] >= Fuh[i][j])) { j++; }
        if (j >= nfu) {
            k = 0;
            while ((k < nfu) && (Fu[k] = Fuh[i][k])) { k++; }
                if (k < nfu) { count++; }
            
        }
    }
	cout << "count=" << count << endl;
	getchar();
    return count;
}
//*************************************************************
void ChoosePareto() {
	int i,j;
	j=0;
	
	for(int jj=0; jj<Lh.size(); jj++) {
	cout << "Lh[" << jj << "]=" << Lh[jj] << endl;
	}
	
	for (i=0; i<=(HH-1); i++) {
		if (Lh[i]==0) {
			j++;
			//setlength(Pareto,j);
			Pareto.resize(j);
			Pareto[j-1]=i;
		}
	}
	for(int jj=0; jj<Pareto.size(); jj++) {
	cout << "Pareto[" << jj << "]=" << Pareto[jj] << endl;
	}
	getchar();
}
//*************************************************************
void GreytoVector(vector < int > y) {
	int i,j,l1,l;
	double g,g1;
	l=c+d;
	//l1=high(y)+1;
	
	for (int jj=0; jj <y.size(); jj++) {
		cout << "y[" << jj << "]=" << y[jj] << endl;
		}
	
	l1=y.size();
	for (i=0; i<y.size(); i++) {
		if (i%l==0) { zb[i]=y[i]; } else { zb[i]=zb[i-1] xor y[i]; }
	}
	j=-1;
	g1=1;
	g=1;
	for (i=0; i<=(c-2); i++) { g1=g1*2; }
	for (i=0; i<=(l1-1); i++) {
		if ((i % l)==0) {
			j++;
			q[j]=0;
			g=g1;
		}
		q[j]=q[j]+g*zb[i];
		
		
		//for (int jj=0; jj <zb.size(); jj++) {
		//cout << "zb[" << jj << "]=" << zb[jj] << endl;
		//}
		//for (int jj=0; jj <q.size(); jj++) {
		//cout << "q[" << jj << "]=" << q[jj] << endl;
		//}
		getchar();
		
		g=g/2;
	}
	g1=g1*2;
	
	for (j=0; j<=(p-1); j++) {
	cout << "q[" << j << "]=" << q[j];	
	q[j]=(qmax[j]-qmin[j])*q[j]/g1+qmin[j]; 
	cout << " qmax[" << j << "]=" << qmax[j] << " qmin[" << j << "]=" << qmin[j] << " q[" << j << "]=" << q[j] << " g1=" << g1 << endl; 
	}
	
	//for (int jj=0; jj <qmin.size(); jj++) {
	//	cout << "qmin[" << jj << "]=" << qmin[jj] << endl;
	//}
	
	//getchar();
}
//*********************************************************
void Initial() {
x[0]=x0[0]+qy[0];
x[1]=x0[1];
x[2]=x0[2];

x[3]=x0[3]+qy[1];
x[4]=x0[4];
x[5]=x0[5];

x[6]=x0[6]+qy[2];
x[7]=x0[7];
x[8]=x0[8];

x[9]=xf[0]-x0[0];
x[10]=xf[1]-x0[1];
x[11]=xf[2]-x0[2];

x[12]=xf[3]-x0[3];
x[13]=xf[4]-x0[4];
x[14]=xf[5]-x0[5];

x[15]=xf[6]-x0[6];
x[16]=xf[7]-x0[7];
x[17]=xf[8]-x0[8];

u[0]=0;
u[1]=0;
u[2]=0;
u[3]=0;
u[4]=0;
u[5]=0;

t=0;

//for (int jj=0; jj <x.size(); jj++) {
//		cout << "x[" << jj << "]=" << x[jj] << endl;
//	}
	
//	getchar();

}
//*********************************************************
void Viewer() {
	double xlf,xlb,xrf,xrb,ylf,ylb,yrf,yrb;
	double steta,cteta;
	int i;
	
	// ïåðâûé ðîáîò
  y[0]=x[0];
  y[1]=x[1];
  y[2]=x[2];
  steta=sin(x[2]);
  cteta=cos(x[2]);
  xlf=x[0]*cteta+x[1]*steta+Lbasc;//A1ac;
  ylf=-x[0]*steta+x[1]*cteta+Lbasc/2;//H1hc;

  xlb=x[0]*cteta+x[1]*steta-Lbasc;//1bc;
  ylb=ylf;

  xrf=xlf;
  yrf=-x[0]*steta+x[1]*cteta-Lbasc/2;

  xrb=xlb;
  yrb=yrf;

  y[3]=xlf*cteta-ylf*steta;
  y[4]=xlf*steta+ylf*cteta;

  y[5]=xlb*cteta-ylb*steta;
  y[6]=xlb*steta+ylb*cteta;

  y[7]=xrb*cteta-yrb*steta;
  y[8]=xrb*steta+yrb*cteta;

  y[9]=xrf*cteta-yrf*steta;
  y[10]=xrf*steta+yrf*cteta;

// âòîðîé ðîáîò
  y[11]=x[3];
  y[12]=x[4];
  y[13]=x[5];
  steta=sin(x[5]);
  cteta=cos(x[5]);
  xlf=x[3]*cteta+x[4]*steta+Lbasc;//A1ac;
  ylf=-x[3]*steta+x[4]*cteta+Lbasc/2;//H1hc;

  xlb=x[3]*cteta+x[4]*steta-Lbasc;//1bc;
  ylb=ylf;

  xrf=xlf;
  yrf=-x[3]*steta+x[4]*cteta-Lbasc/2;

  xrb=xlb;
  yrb=yrf;

  y[14]=xlf*cteta-ylf*steta;
  y[15]=xlf*steta+ylf*cteta;

  y[16]=xlb*cteta-ylb*steta;
  y[17]=xlb*steta+ylb*cteta;

  y[18]=xrb*cteta-yrb*steta;
  y[19]=xrb*steta+yrb*cteta;

  y[20]=xrf*cteta-yrf*steta;
  y[21]=xrf*steta+yrf*cteta;

// òðåòèé ðîáîò
  y[22]=x[6];
  y[23]=x[7];
  y[24]=x[8];
  steta=sin(x[8]);
  cteta=cos(x[8]);
  xlf=x[6]*cteta+x[7]*steta+Lbasc;//A1ac;
  ylf=-x[6]*steta+x[7]*cteta+Lbasc/2;//H1hc;

  xlb=x[6]*cteta+x[7]*steta-Lbasc;//1bc;
  ylb=ylf;

  xrf=xlf;
  yrf=-x[6]*steta+x[7]*cteta-Lbasc/2;

  xrb=xlb;
  yrb=yrf;

  y[25]=xlf*cteta-ylf*steta;
  y[26]=xlf*steta+ylf*cteta;

  y[27]=xlb*cteta-ylb*steta;
  y[28]=xlb*steta+ylb*cteta;

  y[29]=xrb*cteta-yrb*steta;
  y[30]=xrb*steta+yrb*cteta;

  y[31]=xrf*cteta-yrf*steta;
  y[32]=xrf*steta+yrf*cteta;
  
  for (i=0; i<=8; i++) { y[33+i]=x[9+i]; }
}
//*********************************************************
double Check1_4(double xt,double yt, double x1,double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
double d1,d2,d3,d4;
d1=(x1-xt)*(x2-x1)+(y1-yt)*(y2-y1);
d2=(x2-xt)*(x3-x2)+(y2-yt)*(y3-y2);
d3=(x3-xt)*(x4-x3)+(y3-yt)*(y4-y3);
d4=(x4-xt)*(x1-x4)+(y4-yt)*(y1-y4);
if ((d1*d2>0) && (d2*d3>0) && (d3*d4>0))
{
	d1=sqrt(sqrt(x1-xt)+sqrt(y1-yt));
    d2=sqrt(sqrt(x2-xt)+sqrt(y2-yt));
    d3=sqrt(sqrt(x3-xt)+sqrt(y3-yt));
    d4=sqrt(sqrt(x4-xt)+sqrt(y4-yt));
    if (d2<d1) { d1=d2; }
    if (d3<d1) { d1=d3; }
    if (d4<d1) { d1=d4; }
    return d1;
} else {
	return 0;
}	
}
//*********************************************************
double Check2(vector <double> y) {
	int i,j, l1;
	suGA=0;
//узнать что в эту переменную попадает
l1=1;
	for (i=0; i<l1; i++) {
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(y[3+2*j],y[4+2*j],prepc[i][0][0],prepc[i][0][1],prepc[i][1][0],
          prepc[i][1][1],prepc[i][2][0],prepc[i][2][1],prepc[i][3][0],prepc[i][3][1]);
		}
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(prepc[i][j][0],prepc[i][j][1],y[3],y[4],y[5],y[6],y[7],y[8],
           y[9],y[10]);
		}
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(y[14+2*j],y[15+2*j],prepc[i][0][0],prepc[i][0][1],prepc[i][1][0],
          prepc[i][1][1],prepc[i][2][0],prepc[i][2][1],prepc[i][3][0],prepc[i][3][1]);
		}
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(prepc[i][j][0],prepc[i][j][1],y[14],y[15],y[16],y[17],y[18],y[19],
           y[20],y[21]);
		}
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(y[25+2*j],y[26+2*j],prepc[i][0][0],prepc[i][0][1],prepc[i][1][0],
          prepc[i][1][1],prepc[i][2][0],prepc[i][2][1],prepc[i][3][0],prepc[i][3][1]);
		}
		for (j=0; j<=3; j++) {
			suGA=suGA+Check1_4(prepc[i][j][0],prepc[i][j][1],y[25],y[26],y[27],y[28],y[29],y[30],
           y[31],y[32]);
		}
	}
	return suGA;
}
//*********************************************************
double Check3(vector <double> y) {
	int i;
	suGA=0;
	//first with second
	su1GA=0;
	for (i=0; i<=(kRob - 1); i++) { FlagStop[i]=true; }
	for (i=0; i<=3; i++) {
		 su1GA=su1GA+Check1_4(y[3+2*i],y[4+2*i],y[14],y[15], y[16],y[17],
                    y[18],y[19], y[20],y[21]);
	}
	for (i=0; i<=3; i++) {
		su1GA=su1GA+Check1_4(y[14+2*i],y[15+2*i],y[3],y[4], y[5],y[6],
                    y[7],y[8], y[9],y[10]);
	}
	if (su1GA>0) {
		if (Prior[0]>Prior[1]) { FlagStop[1]=false; } else { FlagStop[0]=false; }
	}
	suGA=suGA+su1GA;
	//first with third
	su1GA=0;
	for (i=0; i<=3; i++) {
		su1GA=su1GA+Check1_4(y[3+2*i],y[4+2*i],y[25],y[26], y[27],y[28],
                    y[29],y[30], y[31],y[32]);
	}
	for (i=0; i<=3; i++) {
		su1GA=su1GA+Check1_4(y[25+2*i],y[26+2*i],y[3],y[4], y[5],y[6],
                    y[7],y[8], y[9],y[10]);
	}
	if (su1GA>0) {
		if (Prior[0]>Prior[2]) { FlagStop[2]=false; } else { FlagStop[0]=false; }
	}
	suGA=suGA+su1GA;
	//second with third
	su1GA=0;
	for (i=0; i<=3; i++) {
		su1GA=su1GA+Check1_4(y[14+2*i],y[15+2*i],y[25],y[26], y[27],y[28],
                    y[29],y[30], y[31],y[32]);
	}
	for (i=0; i<=3; i++) {
		 su1GA=su1GA+Check1_4(y[25+2*i],y[26+2*i],y[14],y[15], y[16],y[17],
                    y[18],y[19], y[20],y[21]);
	}
	if (su1GA>0) {
		if (Prior[1]>Prior[2]) { FlagStop[2]=false; } else { FlagStop[1]=false; }
	}
	suGA=suGA+su1GA;
	return suGA;
}
//*************************************************************
void OgrUpr() {
	int i;
	for (i=0; i<=(m-1); i++) {
		if (u[i]>umax[i]) { u[i]=umax[i]; } else { if (u[i]<umin[i]) { u[i]=umin[i]; } }
	}
}
//*********************************************************
void RP(double t1, vector <double> x1, vector <double> f1) {
const double q_0=10.218;
const double q_1=0.44775390625;
const double q_2=1.4932;
const double q_3=0.42098;
const double q_4=14.37744140625;
const double q_5=8.479736328125;
const double q_6=0.28297;

double dx1,dy1,dteta1,dx2,dy2,dteta2,dx3,dy3,dteta3;
double z_0_9,z_0_15,z_1_10,z_1_13,z_1_15;
int i;

Vs[0]=x1[9];
Vs[1]=x1[10];
Vs[2]=x1[11];

Vs[3]=x1[12];
Vs[4]=x1[13];
Vs[5]=x1[14];

Vs[6]=x1[15];
Vs[7]=x1[16];
Vs[8]=x1[17];

//cout << "start RPControl()" << endl;
RPControl();
//cout << "end of RPControl()" << endl;

dx1=x1[9]-x1[0];
dy1=x1[10]-x1[1];
dteta1=x1[11]-x1[2];

dx2=x1[12]-x1[3];
dy2=x1[13]-x1[4];
dteta2=x1[14]-x1[5];

dx3=x1[15]-x1[6];
dy3=x1[16]-x1[7];
dteta3=x1[17]-x1[8];

z_0_9=-q_2*dteta1*dy1+q_1*dy1+q_0*dx1;
z_0_15=Ro_20(z_0_9*q_3+Ro_19(q_0*dx1))+Ro_8(z_0_9*q_3)+Ro_14(dx1);
z_1_10=Ro_4(q_6*dx1)*q_5*dteta1+q_4*dy1+q_6*dx1+
        Ro_5(Ro_4(q_6*dx1)*q_5*dteta1)+Ro_10(q_4*dy1)+Ro_3(q_6*dx1);
  z_1_13=z_1_10+Ro_4(z_1_10)+Ro_18(Ro_4(q_6*dx1)*q_5*dteta1);
  z_1_15=z_1_13*Ro_11(dx1)+Ro_18(z_1_13)+Ro_5(z_1_10);
  u[0]=3*z_0_15/2;
  u[1]=Ro_18(3*z_1_15+Ro_19(3*z_0_15));
  
  z_0_9=-q_2*dteta2*dy2+q_1*dy2+q_0*dx2;
  z_0_15=Ro_20(z_0_9*q_3+Ro_19(q_0*dx2))+Ro_8(z_0_9*q_3)+Ro_14(dx2);
  z_1_10=Ro_4(q_6*dx2)*q_5*dteta2+q_4*dy2+q_6*dx2+
        Ro_5(Ro_4(q_6*dx2)*q_5*dteta2)+Ro_10(q_4*dy2)+Ro_3(q_6*dx2);
  z_1_13=z_1_10+Ro_4(z_1_10)+Ro_18(Ro_4(q_6*dx2)*q_5*dteta2);
  z_1_15=z_1_13*Ro_11(dx2)+Ro_18(z_1_13)+Ro_5(z_1_10);
  u[2]=3*z_0_15/2;
  u[3]=Ro_18(3*z_1_15+Ro_19(3*z_0_15));
  z_0_9=-q_2*dteta3*dy3+q_1*dy3+q_0*dx3;
  z_0_15=Ro_20(z_0_9*q_3+Ro_19(q_0*dx3))+Ro_8(z_0_9*q_3)+Ro_14(dx3);
  z_1_10=Ro_4(q_6*dx3)*q_5*dteta3+q_4*dy3+q_6*dx3+
        Ro_5(Ro_4(q_6*dx3)*q_5*dteta3)+Ro_10(q_4*dy3)+Ro_3(q_6*dx3);
  z_1_13=z_1_10+Ro_4(z_1_10)+Ro_18(Ro_4(q_6*dx3)*q_5*dteta3);
  z_1_15=z_1_13*Ro_11(dx3)+Ro_18(z_1_13)+Ro_5(z_1_10);
  u[4]=3*z_0_15/2;
  u[5]=Ro_18(3*z_1_15+Ro_19(3*z_0_15));

//cout << "start OgrUpr()" << endl;
  OgrUpr();
//cout << "end of OgrUpr()" << endl;

  
  if (!(FlagStop[0])) {
	  u[0]=0;
    u[1]=0;
  }
  
  if (!(FlagStop[1])) {
	   u[2]=0;
    u[3]=0;
  }
  
  if (!(FlagStop[2])) {
	  u[4]=0;
    u[5]=0;
  }
  
  
  f1[0]=u[0]*cos(x1[2]);
  f1[1]=u[0]*sin(x1[2]);
  f1[2]=(u[0]/Lbasc)*sin(u[1])/cos(u[1]);

  f1[3]=u[2]*cos(x1[5]);
  f1[4]=u[2]*sin(x1[5]);
  f1[5]=(u[2]/Lbasc)*sin(u[3])/cos(u[3]);

  f1[6]=u[4]*cos(x1[8]);
  f1[7]=u[4]*sin(x1[8]);
  f1[8]=(u[4]/Lbasc)*sin(u[5])/cos(u[5]);

  f1[9]=z[M_Out[0][0]-1][M_Out[0][1]];
  f1[10]=z[M_Out[1][0]-1][M_Out[1][1]];
  f1[11]=z[M_Out[2][0]-1][M_Out[2][1]];

  f1[12]=z[M_Out[3][0]-1][M_Out[3][1]]; 
  f1[13]=z[M_Out[4][0]-1][M_Out[4][1]];
  f1[14]=z[M_Out[5][0]-1][M_Out[5][1]];


  f1[15]=z[M_Out[6][0]-1][M_Out[6][1]];
  f1[16]=z[M_Out[7][0]-1][M_Out[7][1]];  
  f1[17]=z[M_Out[8][0]-1][M_Out[8][1]];

  for (i=0; i<=(n-1); i++) { if (fabs(f1[i])>numeric_limits<double>::infinity( )) { f1[i]=Ro_10(f1[i])*numeric_limits<double>::infinity( ); } }
}
//*************************************************************
void Euler2() {
//t не меняется
	int i;
	//cout << "start rp1" << endl;
	RP(t,x,fa);
	//cout << "end of rp1" << endl;
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*fa[i]; }
	//cout << "start rp2" << endl;
	RP(t+dt,xs,fb);
	//cout << "end of rp2" << endl;
	for (i=0; i<=(n-1); i++) { x[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	t=t+dt;
	//cout << "t=" << t << endl;
	//getchar();
}
//*************************************************************
void Euler3() {
	int i;
	RP(t,x,fa);
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*fa[i]; }
	RP(t+dt,xs,fb);
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	RP(t+dt,xs,fb);
	for (i=0; i<=(n-1); i++) { x[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	t=t+dt;
}
//*************************************************************
void Euler4() {
	int i;
	RP(t,x,fa);
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*fa[i]; }
	RP(t+dt,xs,fb);
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	RP(t+dt,xs,fb);
	for (i=0; i<=(n-1); i++) { xs[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	RP(t+dt,xs,fb);
	for (i=0; i<=(n-1); i++) { x[i]=x[i]+dt*(fa[i]+fb[i])/2; }
	t=t+dt;
}
//*********************************************************
double Normdist(vector <double> x1, vector <double> xf1) {
	double sum,aa;
	int i;
	sum=0;
	//cout << "x1.size()=" << x1.size() << " xf1.size()=" << xf1.size() << endl;
	
		//for(int jj=0; jj<xf1.size(); jj++) {
		//cout << "xf1[" << jj << "]=" << fixed << xf1[jj] << endl;
		//}
		//for(int jj=0; jj<x1.size(); jj++) {
		//cout << "x1[" << jj << "]=" << fixed << x1[jj] << endl;
		//}
	
 	for (i=0;i<xf1.size(); i++) {
	//cout << "i=" << i << endl;
		aa=fabs(xf1[i]-x1[i]);
		

		
		//cout << "aa=" << aa << " sum=" << fixed <<  sum << endl;
		//getchar();
		if (aa>sum) { sum=aa; }
	}
	return sum;
}
//*********************************************************
void Func(vector <double> Fu) {
	const double shtraf=0.5;
	double sumpen,promah,pr1;
	int i;
	//cout << "Start Initial()" << endl;
	Initial();
	//cout << "End of Initial()" << endl;
	sumpen=0;
	//for (i=0; i<kR-1; i++) { Cs[i]=q[i]; }
	Cs=q;
	
	
	//for(int jj=0; jj<y.size(); jj++) {
	//	cout << "y[" << jj << "]=" << fixed << y[jj] << endl;
	//	}
	
	//getchar();
	
	do {
		Viewer();
		sumpen=sumpen+Check2(y);
    	sumpen=sumpen+Check3(y);
		
		cout << "sumpen=" << fixed << sumpen << endl;
		
		//cout << "Start Euler2()" << endl;
    	Euler2();
		//cout << "End of Euler2()" << endl; 
		//cout << "t=" << t << " tf=" << tf << " Normdist(x,xf)=" << fixed << setprecision(3) << Normdist(x,xf) << " epsterm=" << fixed << setprecision(3) << epsterm << endl;
}
	while (!((t>tf) || (Normdist(x,xf)<epsterm)));
	promah=0;
	//cout << "point" << endl;
	for (i=0; i<xf1.size(); i++) { 
		pr1=fabs(x[i]-xfc[i]);
		if (pr1>promah) { promah=pr1; }
	}
	Fu[0]=t+shtraf*sumpen*dt;
	cout << "t=" << t << " shtraf=" << shtraf << " sumpen=" << sumpen << " dt=" << dt << endl;
	Fu[1]=promah+shtraf*sumpen*dt;
	cout << "promah=" << promah << " shtraf=" << shtraf << " sumpen=" << sumpen << " dt=" << dt << endl;
    //cout << "Fu[0]=" << Fu[0] << " Fu[1]=" << Fu[1] << endl;
	getchar();
}
//*************************************************************
void Integr() {
	int i,j;
	bool flag;
	for (i=0; i<=ny-1; i++) { EAix[i]=0; }
	for (i=0; i<=nfu-1; i++) { su[i]=0; }
	do {
	for (i=0; i<=ny-1; i++) { qy[i]=qymin[i]+stepsqy[i]*EAix[i]; }
//cout << "Start Func(su1)" << endl;
	Func(su1);
//cout << "End of Func(su1)" << endl;
	for (i=0; i<=nfu-1; i++) { su[i]=su[i]+su1[i]; }
	LexPM(EAix,flag); }
	while (flag);
}
//*************************************************************
void Func0(vector <double> Fu) {
	int i;
//cout << "Start Integr()" << endl;
	Integr();
//cout << "End of Integr()" << endl;
	for (i=0; i<=(nfu-1); i++) { Fu[i]=su[i]; }
}
//*************************************************************
void ImproveChrom(vector <double> q, vector <vector <int> > StrChrom) {
int i,j,k;
bool flag;
SetPsi(Psi0);
Func0(Fu1);
k=-1;
for (i=0; i<=lchr-1; i++) { 
 Variations(StrChrom[i]);
 Func0(Fu2);
 flag=true;
 for (j=0; j<=nfu-1; j++) { if (Fu2[j]>Fu1[j]) { flag=false; } }
 if (flag) {
	 for (j=0; j<=nfu-1; j++) { Fu1[j]=Fu2[j]; }
	 k=i;
 }
 for (i=k+1; i<=lchr-1;i++) {
	 for (j=0; j<=3; j++) { StrChrom[i][j]=0; }
 }
}
}

