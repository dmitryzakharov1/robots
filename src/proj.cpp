#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include "functions.h"



int main(int argc, const char* argv[])
{
//************************  2.  *******************************
//Определяем размерности и задаем значения параметров алгоритма
//*************************************************************

//номера унарных операций
//SetLength(O1s,kW);
O1s.resize(kW);
for (i=0; i<=kW-1; i++) { O1s[i]=O1sc[i]; }
//номера бинарных операций
//SetLength(O2s,kV);
O2s.resize(kV);
for (i=0; i<=kV-1;i++) { O2s[i]=O2sc[i]; }	
//базисная матрица
//SetLength(Psi,kL,L,L);
//Psi.resize(kL,L,L);
Psi.resize(kL);
for (int ii=0; ii< 3; ii++) {
	Psi[ii].resize(L);
	for (int jj=0; jj<L; jj++) {
	Psi[ii][jj].resize(L);
	}
}



//SetLength(Psi0,kL,L,L);
//Psi0.resize(kL,L,L);
Psi0.resize(kL);
for (int ii=0; ii< 3; ii++) {
	Psi0[ii].resize(L);
	for (int jj=0; jj<L; jj++) {
	Psi0[ii][jj].resize(L);
	}
}

for (k=0; k<kL-1; k++) {
	for (i=0; i<L-1; i++) {
		for (j=0; j<(L-1); j++) {
			Psi[k][i][j]=PsiMultBasc[k][i][j];
           Psi0[k][i][j]=PsiMultBasc[k][i][j];
		   //Psi[k][0][0]=1;
		}
	}
}



//матрица входов
//SetLength(M_Entr,kL,2*kInp);
//M_Entr.resize(kL,2*kInP);
M_Entr.resize(kL);
for (int ii=0; ii<2; ii++) {
	M_Entr[ii].resize(2*kInP);
}

//cout << M_Entr.size() << " " << M_Entr[1].size() << " | " << kL-1 << " " << 2*kInP - 1 << endl;
for (k=0; k<kL-1; k++) {
	for (i=0; i<(2*kInP - 1); i++) {  
	//M_Entr[k][i]=M_Entrc[k][i];
	//M_Entr[k][i]=1;
	cout << "M_Entrc[" << k << "][" << i << "] = " << M_Entrc[k][i] << endl;
	}
}
/*
//матрица выходов
//SetLength(M_Out,Mout,2);
//M_Out.resize(Mout,2);
M_Out.resize(Mout);
for (int ii=0; ii< 2; ii++) {
	M_Out[ii].resize(2);
}


for (i=0; i<=(Mout - 1); i++) {
	for (j=0; j<=1; j++) { M_Out[i][j]=M_Outc[i][j]; }
}
//параметры модели
dt=dtc;
tf=tfc;
//SetLength(qmax,p);
qmax.resize(p);
//SetLength(qmin,p);
qmin.resize(p);
//SetLength(q,p);
q.resize(p);
for (i=0; i<=(p-1); i++) {
	qmax[i]=qmaxc[i];
    qmin[i]=qminc[i];
    q[i]=qc[i];
}
//неопределенные параметры
//Setlength(qymax,ny);
qymax.resize(ny);
//Setlength(qymin,ny);
qymin.resize(ny);
//Setlength(stepsqy,ny);
stepsqy.resize(ny);
//SetLength(qy,ny);
qy.resize(ny);
//SetLength(EAixmax,ny);
EAixmax.resize(ny);
//SetLength(EAix,ny);
EAix.resize(ny);
  
for (i=0; i<=ny-1; i++) {
	qymin[i]=qyminc[i];
    qymax[i]=qymaxc[i];
    stepsqy[i]=stepsqyc[i];
    qy[i]=qymin[i];
}
for (i=0; i<=ny-1; i++) { EAix[i]=trunc((qymax[i]-qymin[i])/stepsqy[i]); }
SetEAixmax(EAix);

//вектор состояния объекта
//Setlength(x0,n);
x0.resize(n);
//SetLength(x,n);
x.resize(n);
//SetLength(xs,n);
xs.resize(n);
//SetLength(fa,n);
fa.resize(n);
//SetLength(fb,n);
fb.resize(n);
//Setlength(xf,high(xfc)+1); 
xf.resize(10); 
for (i=0; i<=n-1; i++) { x0[i]=x0c[i]; }
for (i=0; i<=(xf.size()-1); i++) { xf[i]=xfc[i]; }
//ограничения на управление
//SetLength(u,m);
u.resize(m);
//Setlength(umin,m);
umin.resize(m);
//Setlength(umax,m);
umax.resize(m);  
for (i=0; i<=m-1; i++) { umin[i]=uminc[i]; }
for (i=0; i<=m-1; i++) { umax[i]=umaxc[i]; }
//вектор наблюдения
//SetLength(y,lv);
y.resize(lv);
  //хранение промежуточных значений функционалов
//SetLength(su,nfu);
su.resize(nfu);
//SetLength(su1,nfu);
su1.resize(nfu);
  
  //0-ой слой (переменных и параметров)
//SetLength(V_Entr,kP+kR);
V_Entr.resize(kP+kR);
//SetLength(Pnum,kP);
Pnum.resize(kP);
//SetLength(Rnum,kR);
Rnum.resize(kR);
//SetLength(Vs,kP);
Vs.resize(kP);
//SetLength(Cs,kR);
Cs.resize(kR);
	  
for (i=0; i<=kP-1; i++) { Pnum[i]=Pnumc[i]; }	
for (i=0; i<=kR-1; i++) { Rnum[i]=Rnumc[i]; } 
for (i=0; i<=kR-1; i++) { Cs[i]=q[i]; } 

  SetV_Entr;
 // SetLength(Prior,kRob);
Prior.resize(kRob);
  //SetLength(FlagStop,kRob);
FlagStop.resize(kRob);
  
 for (i=0; i<=kRob-1; k++) {
	FlagStop[i]=true;
    Prior[i]=kRob-i;
 }
  
//  Setlength(PopChrStr,HH,lchr);
PopChrStr.resize(HH);
for (int ii=0; ii< 2; ii++) {
	PopChrStr[ii].resize(lchr);
	for (int jj=0; jj<lchr; jj++) {
	PopChrStr[ii][jj].resize(5);
	}
}


  //Setlength(PopChrPar,HH,p*(c+d));
PopChrPar.resize(HH);
for (int ii=0; ii< 2; ii++) {
	PopChrPar[ii].resize(p*(c+d));
}

  //Setlength(Fuh,HH,nfu);
Fuh.resize(HH);
for (int ii=0; ii< 2; ii++) {
	Fuh[ii].resize(nfu);
}

//  SetLength(FuhNorm,HH,nfu);
FuhNorm.resize(HH);
for (int ii=0; ii< 2; ii++) {
	FuhNorm[ii].resize(nfu);
}

//Setlength(Shtraf,nfu);
Shtraf.resize(nfu);
	
for (i=0; i<=nfu-1; i++) { Shtraf[i]=1; }
  //Setlength(Lh,HH);
Lh.resize(HH);
  //Setlength(Fu1,nfu);
Fu1.resize(nfu);
  //Setlength(Fu2,nfu);
Fu2.resize(nfu);
  //Setlength(Fu3,nfu);
Fu3.resize(nfu);
  //Setlength(Fu4,nfu);
Fu4.resize(nfu);
  //Setlength(Son1s,lchr);
Son1s.resize(lchr);
for (int ii=0; ii< lchr; ii++) {
	Son1s[ii].resize(5);
}
  //Setlength(Son2s,lchr);
Son2s.resize(lchr);
for (int ii=0; ii< lchr; ii++) {
	Son2s[ii].resize(5);
}

  //Setlength(Son3s,lchr);
Son3s.resize(lchr);
for (int ii=0; ii< lchr; ii++) {
	Son3s[ii].resize(5);
}
//  Setlength(Son4s,lchr);
Son4s.resize(lchr);
for (int ii=0; ii< lchr; ii++) {
	Son4s[ii].resize(5);
}
//Setlength(Son1p,p*(c+d));
Son1p.resize(p*(c+d));
//  Setlength(Son2p,p*(c+d));
Son2p.resize(p*(c+d));
//  Setlength(Son3p,p*(c+d));
Son3p.resize(p*(c+d));
  //Setlength(Son4p,p*(c+d));
Son4p.resize(p*(c+d));
//SetLength(zb,p*(c+d));
zb.resize(p*(c+d)); 
//SetLength(z,kL,L);
z.resize(kL);
for (int ii=0; ii< 2; ii++) {
	z[ii].resize(L);
}


//SetLength(zs,kL,L);
zs.resize(kL);
for (int ii=0; ii< 2; ii++) {
	zs[ii].resize(L);
}
  
  //************************  2.  *******************************
  //                 ГЕНЕТИЧЕСКИЙ АЛГОРИТМ
  //*************************************************************

  //задаем базисную матрицу сетевого оператора
  SetPsiBas(Psi);
  //генерируем начальную (нулевую) популяцию векторов вариаций
  VectortoGrey(PopChrPar[0]); 
for (i=0; i<= lchr-1; i++) {
	for (j=0; j<=4; j++) { PopChrStr[0][i][j]=0; }
}
//наполняем вектора случайными возможными значениями
for (i=1; i<=HH-1; i++) {
	for (j=0; j<=lchr-1; j++) { GenVar(PopChrStr[i][j]); }
	for (j=0; j<=(p*(c+d)-1); j++) { PopChrPar[i][j]=rand()%2+1; }
}
// считаем значения функционалов для каждой хромосомы
for (i=0; i<=(HH-1); i++) {
//для структурной части: берем текущую базисную матрицу
SetPsi(Psi0);
//вносим в нее вариации;
for (j=0; j<=(lchr-1); j++) { Variations(PopChrStr[i][j]); }
//для параметрической:
//текущие значения параметров в коде Грея переводим в вектор параметров
GreytoVector(PopChrPar[i]);
// для каждой полученной матрицы сетевого оператора и вектора параметров,
//определяющих искомое математическое выражение, оцениваем
//функции приспособленности
Func0(Fuh[i]);
}
//caculating distances to Pareto set
for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
//Start of cycle for generations
pt=1;  // first current generation

do {
    //start of cycle for crossovering
    rt=1;//first couple for crossoving
    do {
      //select of two parents
      k1=rand()% HH + 1;
      k2=rand()% HH + 1;
      ksi=rand();
      if ((ksi<(1+alfa*Lh[k1])/(1+Lh[k1])) ||
         (ksi<(1+alfa*Lh[k2])/(1+Lh[k2]))) {
      
        //if true
        ks1=random() % lchr+1;
        ks2=random() % (p*(c+d))+1;
        //crossover creating four sons
	
		for (i=0; i<=lchr-1; i++) {
		Son1s[i]=PopChrStr[k1][i];
        Son2s[i]=PopChrStr[k2][i];
		}
		
		for (i=0; i<=ks2-1; i++) {
		Son1p[i]=PopChrPar[k1][i];
        Son2p[i]=PopChrPar[k2][i];
        Son3p[i]=PopChrPar[k1][i];
        Son4p[i]=PopChrPar[k2][i];
		}
		
		for (i=ks2; i<=(p*(c+d)-1); i++) {
			Son1p[i]=PopChrPar[k2][i];
          Son2p[i]=PopChrPar[k1][i];
          Son3p[i]=PopChrPar[k2][i];
          Son4p[i]=PopChrPar[k1][i];
		}
				
		for (i=0; i<=(ks1-1); i++) {
			Son3s[i]=PopChrStr[k1][i];
          Son4s[i]=PopChrStr[k2][i];
		}
		
		for (i=ks1; i<=lchr-1; i++) {
		Son3s[i]=PopChrStr[k2][i];
        Son4s[i]=PopChrStr[k1][i];
		}
		
        //mutation for 1st son	
		if (rand()<pmut) {
		Son1p[random()%(p*(c+d))+1]=rand()%2+1;
          GenVar(Son1s[rand()%lchr+1]);
		}
		
        //functional for 1st son
        SetPsi(Psi0);		  
		for (j=0; j<=lchr-1; j++) { Variations(Son1s[j]); }
        GreytoVector(Son1p);
        Func0(Fu1);
        //Distance for 1st son
        L1=Rast(Fu1);
        //Chromosome with biggest distance to Pareto set
        lmax=Lh[0];
        imax=0; 
		for (i=1; i<=HH-1; i++) {
			if (Lh[i]>lmax) {
				lmax=Lh[i];
            imax=i;
			}
		}  
		  
        //if distance to Pareto set 1st son is less than biggest distance
        //...in population then make substitution
		if (L1<lmax) {
			for (i=0; i<=lchr-1; i++) { PopChrStr[imax][i]=Son1s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax][i]=Son1p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax][i]=Fu1[i]; }
		}
		
        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
        //mutation for 2nd son	
		if (rand()<pmut) {
		Son2p[rand()%(p*(c+d))+1]=rand()%2+1;
          GenVar(Son2s[rand()%lchr+1]);
		}
        //functional for 2nd son
        SetPsi(Psi0);		  
		for (j=0; lchr-1; j++) { Variations(Son2s[j]); }
        GreytoVector(Son2p);
        Func0(Fu2);
        //Distance for 2nd son
        L2=Rast(Fu2);
        //Chromosome with biggest distance to Pareto set
        lmax=Lh[0];
        imax=0;
		  
		for (i=1; i<=HH-1; i++) {
			if (Lh[i]>lmax) {
				lmax=Lh[i];
            imax=i;
			}
		}

        //if distance to Pareto set 2nd son is less than biggest distance
        //...in population then make substitution
		if (L2<lmax) {
			for (i=0; i<=(lchr-1); i++) { PopChrStr[imax][i]=Son2s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax][i]=Son2p[i]; }
			for (i=0; i<= nfu-1; i++) { Fuh[imax][i]=Fu2[i]; }
		}

        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
        //mutation for 3rd son
		if (rand()<pmut) { 
		Son3p[rand()%(p*(c+d))+1]=rand()%2+1;
        GenVar(Son3s[rand()%lchr+1]);
		}
        //functional for 3rd son
        SetPsi(Psi0);	  
		for (j=0; j<=(lchr-1); j++) { Variations(Son3s[j]); }
        GreytoVector(Son3p);
        Func0(Fu3);
        //Distance for 3rd son
        L3=Rast(Fu3);
        //Chromosome with biggest distance to Pareto set
        lmax=Lh[0];
        imax=0;
		for (i=1; i<=HH-1; i++) {
			lmax=Lh[i];
            imax=i;
		}

        //if distance to Pareto set 3rd son is less than biggest distance
        //...in population then make substitution
		if (L3<lmax) {
			for (i=0; i<=lchr-1; i++) { PopChrStr[imax][i]=Son3s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax][i]=Son3p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax][i]=Fu3[i]; }
		}
        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
		
        //mutation for 4th son		
		if (rand()<pmut) { 
		Son4p[rand()%(p*(c+d))+1]=rand()%2+1;
        GenVar(Son4s[rand()%lchr+1]);
		}
		
        //functional for 4th son
        SetPsi(Psi0);		  
		for (j=0; j<=lchr-1; j++) { Variations(Son4s[j]); }
        GreytoVector(Son4p);
        Func0(Fu4);
        //Distance for 4th son
        L4=Rast(Fu4);
        //Chromosome with biggest distance to Pareto set
        lmax=Lh[0];
        imax=0;
		  
		for (i=1;i<=HH-1; i++) {
			if (Lh[i]>lmax) {
			lmax=Lh[i];
            imax=i;
			}
		}
		
        //if distance to Pareto set 4th son is less than biggest distance
        //...in population then make substitution
		if (L4<lmax) {
			for (i=0; i<=lchr-1; i++) {  PopChrStr[imax][i]=Son4s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax][i]=Son4p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax][i]=Fu4[i]; }
		}
		
        //calculating all distances for population	
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
		 }
      rt++;
      //End of cycle for crossoving
	}
    while (!(rt>RR));
	
	// generating new chromosomes
    // Checking Epoch
    pt++;
	
	    //if epoch is over then changing basic
    
		if (pt%Epo==0) 
		{
		  //... на наиболее близкую хромосому к утопической
		  // хромосоме в пространстве нормированных критериев
		  for (i=0; i<=nfu-1; i++) {
			  Fumax=Fuh[0][i];
			Fumin=Fuh[0][i];
			 // ищем максимальное и минимальное значения по каждому функционалу
			 for (k=0; k<=HH-1; k++) {
				 if (Fuh[k][i]>Fumax) { Fumax=Fuh[k][i]; } else { if (Fuh[k][i]< Fumin) { Fumin=Fuh[k][i]; } }
			 }
			 // нормируем критерии, поделив каждое значение на разность между
			// максимумом и минимумом
			if (Fumax!=Fumin) {
				for (k=0; k<=HH-1; k++) { FuhNorm[k][i]=Fuh[k][i]/(Fumax-Fumin); }
			}
		  }
		  
		  // находим хромосому с наименьшей величиной нормы нормированных критериев
		  k=0;
		  suGA=0;	
		for (i=0; i<=nfu-1; i++) { suGA=suGA+Shtraf[i]*sqrt(FuhNorm[0][i]); }
		  suGA=sqrt(suGA);
	//      su=FuhNorm[0,1];
		  for (i=1; i<=HH-1; i++) {
			  su1GA=0;
			  for (j=0; j<=nfu-1; j++) { su1GA=su1GA+Shtraf[j]*sqrt(FuhNorm[i][j]); }
			  su1GA=FuhNorm[i][1];
			  if (su1GA<suGA) {
				  suGA=su1GA;
			  k=i;
			  }
		  }
		  
		  // заменяем базис
		  // строим матрицу для найденной хромосомы
		  SetPsi(Psi0);
		for (j=0; j<=lchr-1; j++) { Variations(PopChrStr[k][j]); }
		  // меняем базисную матрицу на новую
		  SetPsiBas(Psi);
		  //генерируем тождественную хромосому
		for (i=0; i<=lchr-1; i++) { for (j=0; j<=3; j++) { PopChrStr[0][i][j]=0; } }
		for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[0][i]=PopChrPar[k][i]; }
		  //вычисляем все фунционалы для всей популяции	  
		  for (i=0; i<=HH-1; i++) { 
		  SetPsi(Psi0);
		  for (j=0; j<=lchr-1; j++) { Variations(PopChrStr[i][j]); }
		  GreytoVector(PopChrPar[i]);
			Func0(Fuh[i]);
		  }
		  
		  // формируем элиту
		  for (i=0; i<=kel-1; i++) {
			j=rand()% (HH-1)+1+1;
			GreytoVector(PopChrPar[j]);
			ImproveChrom(q,PopChrStr[j]); 
		  }
		  
		  //вычисляем новые расстояния
		for (i=0; i<=HH-1; i++) {  Lh[i]=Rast(Fuh[i]); }
		}
    //конец цикла поколений
	}

while (!(pt>PP));
  ChoosePareto();
  //Сохраним множество Парето в файл
  ofstream fout("output.txt");
  kol=Pareto.size();
  //сортируем решения по возрастанию
	  
for (i=0; i<=kol-2; i++) {
	for (j=i+1; j<=kol-1; j++) {
		if (Fuh[Pareto[i],0]>Fuh[Pareto[j],0]) {
			k=Pareto[i];
        Pareto[i]=Pareto[j];
        Pareto[j]=k;
		}
	}
}	  

  
  for (i=0; i<=kol-1; i++) {
	  //fout << "solution No " << doubeVecToStr(Pareto[i]) << ", ";
	  fout << "solution No " << endl;
	  for (j=0; j<=nfu-1; j++) {
		//  fout << " F" << j << "=" << doubeVecToStr(Fuh[Pareto[i],j]) << ", ";
		fout << " F" << endl;
	  }
  }
  
  fout.close();

		 

///////////////
// Удаление массивов
		 
	 */	 
		
		//cout << "End of program\n";
//конец
return 0;
}

