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
#include "functions.h"



int main(int argc, const char* argv[])
{
//************************  2.  *******************************
//Определяем размерности и задаем значения параметров алгоритма
//*************************************************************

//номера унарных операций
SetLength(O1s,kW);	
for (i=0; i<=kW-1; i++) { O1s[i]=O1sc[i]; }
//номера бинарных операций
SetLength(O2s,kV);
for (i=0; i<=kV-1;i++) { O2s[i]:=O2sc[i]; }	
//базисная матрица
SetLength(Psi,kL,L,L);
SetLength(Psi0,kL,L,L);
for (k=0; k<=kL-1; k++) {
	for (i=0; i<=L-1; i++) {
		for (j=0; j<=(L-1); j++) {
			Psi[k,i,j]=PsiMultBasc[k,i,j];
            Psi0[k,i,j]=PsiMultBasc[k,i,j];
		}
	}
}

//матрица входов
SetLength(M_Entr,kL,2*kInp);
for (k=0; k<=kL-1; k++) {
	for (i=0; i<=(2*kinp - 1); i++) {  M_Entr[k,i]=M_Entrc[k,i]; }
}
//матрица выходов
SetLength(M_Out,Mout,2);
for (i=0; i<=(Mout - 1); i++) {
	for (j=0; j<=1; j++) { M_Out[i,j]=M_Outc[i,j]; }
}
//параметры модели
dt=dtc;
tf=tfc;
SetLength(qmax,p);
SetLength(qmin,p);
SetLength(q,p);
for (i=0; i<=(p-1); i++) {
	qmax[i]=qmaxc[i];
    qmin[i]=qminc[i];
    q[i]=qc[i];
}
//неопределенные параметры
Setlength(qymax,ny);
Setlength(qymin,ny);
Setlength(stepsqy,ny);
SetLength(qy,ny);
SetLength(EAixmax,ny);
SetLength(EAix,ny);
  
for (i=0; i<=ny-1; i++) {
	qymin[i]=qyminc[i];
    qymax[i]=qymaxc[i];
    stepsqy[i]=stepsqyc[i];
    qy[i]=qymin[i];
}
for (i=0; i<=ny-1; i++) { EAix[i]=trunc((qymax[i]-qymin[i])/stepsqy[i]); }
SetEAixmax(EAix);

//вектор состояния объекта
  Setlength(x0,n);
  SetLength(x,n);
  SetLength(xs,n);
  SetLength(fa,n);
  SetLength(fb,n);
  Setlength(xf,high(xfc)+1);  
for (i=0; i<=n-1; i++) { x0[i]=x0c[i]; }
for (i=0; i<=high(xf); i++) { xf[i]=xfc[i]; }
//ограничения на управление
  SetLength(u,m);
  Setlength(umin,m);
  Setlength(umax,m);  
for (i=0; i<=m-1; i++) { umin[i]=uminc[i]; }
for (i=0; i<=m-1; i++) { umax[i]:=umaxc[i]; }
//вектор наблюдения
  SetLength(y,lv);
  //хранение промежуточных значений функционалов
  SetLength(su,nfu);
  SetLength(su1,nfu);
  
  //0-ой слой (переменных и параметров)
  SetLength(V_Entr,kP+kR);
  SetLength(Pnum,kP);
  SetLength(Rnum,kR);
  SetLength(Vs,kP);
  SetLength(Cs,kR);
	  
for (i=0; i<=kP-1; i++) { Pnum[i]:=Pnumc[i]; }	
for (i=0; i<=kR-1; i++) { Rnum[i]=Rnumc[i]; } 
for (i=0; i<=kR-1; i++) { Cs[i]=q[i]; } 

  SetV_Entr;
  SetLength(Prior,kRob);
  SetLength(FlagStop,kRob);
  
 for (i=0; i<=kRob-1; k++) {
	FlagStop[i]=true;
    Prior[i]=kRob-i;
 }
  
  Setlength(PopChrStr,HH,lchr);
  Setlength(PopChrPar,HH,p*(c+d));
  Setlength(Fuh,HH,nfu);
  SetLength(FuhNorm,HH,nfu);
  Setlength(Shtraf,nfu);
	
for (i=0; i<=nfu-1; i++) { Shtraf[i]:=1; }
  Setlength(Lh,HH);
  Setlength(Fu1,nfu);
  Setlength(Fu2,nfu);
  Setlength(Fu3,nfu);
  Setlength(Fu4,nfu);
  Setlength(Son1s,lchr);
  Setlength(Son2s,lchr);
  Setlength(Son3s,lchr);
  Setlength(Son4s,lchr);
  Setlength(Son1p,p*(c+d));
  Setlength(Son2p,p*(c+d));
  Setlength(Son3p,p*(c+d));
  Setlength(Son4p,p*(c+d));
  SetLength(zb,p*(c+d));
  SetLength(z,kL,L);
  SetLength(zs,kL,L);
  
  //************************  2.  *******************************
  //                 ГЕНЕТИЧЕСКИЙ АЛГОРИТМ
  //*************************************************************

  //задаем базисную матрицу сетевого оператора
  SetPsiBas(Psi);
  //генерируем начальную (нулевую) популяцию векторов вариаций
  VectortoGrey(PopChrPar[0]); 
for (i=0; i<= lchr-1; i++) {
	for (j=0; j<=4; j++) { PopChrStr[0,i,j]=0; }
}
//наполняем вектора случайными возможными значениями
for (i=1; i<=HH-1; i++) {
	for (j=0; j<=lchr-1; j++) { GenVar(PopChrStr[i,j]); }
	for (j=0; j<=(p*(c+d)-1); j++) { PopChrPar[i,j]=random(2); }
}
// считаем значения функционалов для каждой хромосомы
for (i=0; i<=(HH-1); i++) {
//для структурной части: берем текущую базисную матрицу
SetPsi(Psi0);
//вносим в нее вариации;
for (j=0; j<=(lchr-1); j++) { Variations(PopChrStr[i,j]); }
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

repeat
    //start of cycle for crossovering
    rt=1;//first couple for crossoving
    repeat
      //select of two parents
      k1=random(HH);
      k2=random(HH);
      ksi=random;
      if ((ksi<(1+alfa*Lh[k1])/(1+Lh[k1])) ||
         (ksi<(1+alfa*Lh[k2])/(1+Lh[k2]))) {
      
        //if true
        ks1=random(lchr);
        ks2=random(p*(c+d));
        //crossover creating four sons
	
		for (i=0; i<=lchr-1; i++) {
		Son1s[i]=PopChrStr[k1,i];
        Son2s[i]=PopChrStr[k2,i];
		}
		
		for (i=0; i<=ks2-1; i++) {
		Son1p[i]=PopChrPar[k1,i];
        Son2p[i]=PopChrPar[k2,i];
        Son3p[i]=PopChrPar[k1,i];
        Son4p[i]=PopChrPar[k2,i];
		}
		
		for (i=ks2; i<=(p*(c+d)-1); i++) {
			Son1p[i]=PopChrPar[k2,i];
          Son2p[i]=PopChrPar[k1,i];
          Son3p[i]=PopChrPar[k2,i];
          Son4p[i]=PopChrPar[k1,i];
		}
				
		for (i=0; i<=(ks1-1); i++) {
			Son3s[i]=PopChrStr[k1,i];
          Son4s[i]=PopChrStr[k2,i];
		}
		
		for (i=ks1; i<=lchr-1; i++) {
		Son3s[i]=PopChrStr[k2,i];
        Son4s[i]=PopChrStr[k1,i];
		}
		
        //mutation for 1st son	
		if (random<pmut) {
		Son1p[random(p*(c+d))]=random(2);
          GenVar(son1s[random(lchr)]);
		}
		
        //functional for 1st son
        SetPsi(Psi0);		  
		for (j=0; j<=lchr-1; j++) { Variations(son1s[j]); }
        GreytoVector(Son1p);
        Func0(Fu1);
        //Distance for 1st son
        L1=Rast(Fu1);
        //Chromosome with biggest distance to Pareto set
        Lmax=Lh[0];
        imax=0; 
		for (i=1; i<=HH-1; i++) {
			if (Lh[i]>Lmax) {
				Lmax=Lh[i];
            imax=i;
			}
		}  
		  
        //if distance to Pareto set 1st son is less than biggest distance
        //...in population then make substitution
		if (L1<Lmax) {
			for (i=0; i<=lchr-1; i++) { PopChrStr[imax,i]=son1s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax,i]=Son1p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax,i]:=Fu1[i]; }
		}
		
		
		
        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
        //mutation for 2nd son	
		if (random<pmut) {
		Son2p[random(p*(c+d))]=random(2);
          GenVar(Son2s[random(lchr)]);
		}
        //functional for 2nd son
        SetPsi(Psi0);		  
		for (j=0; lchr-1; j++) { Variations(Son2s[j]); }
        GreytoVector(Son2p);
        Func0(Fu2);
        //Distance for 2nd son
        L2=Rast(Fu2);
        //Chromosome with biggest distance to Pareto set
        Lmax=Lh[0];
        imax=0;
		  
		for (i=1; i<=HH-1; i++) {
			if (Lh[i]>Lmax) {
				Lmax=Lh[i];
            imax=i;
			}
		}
		

        //if distance to Pareto set 2nd son is less than biggest distance
        //...in population then make substitution
		if (L2<Lmax) {
			for (i=0; i<=(lchr-1); i++) { PopChrStr[imax,i]=Son2s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax,i]=Son2p[i]; }
			for (i=0; i<= nfu-1; i++) { Fuh[imax,i]=Fu2[i]; }
		}
		
		
		
		///////////////////////////////////////////////
		
		
		
        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
        //mutation for 3rd son
		if (random<pmut) { 
		Son3p[random(p*(c+d))]=random(2);
        GenVar(Son3s[random(lchr)]);
		}
        //functional for 3rd son
        SetPsi(Psi0);	  
		for (j=0; j<=(lchr-1); j++) { Variations(Son3s[j]); }
        GreytoVector(Son3p);
        Func0(Fu3);
        //Distance for 3rd son
        L3=Rast(Fu3);
        //Chromosome with biggest distance to Pareto set
        Lmax=Lh[0];
        imax=0;
		for (i=1; i<=HH-1; i++) {
			Lmax=Lh[i];
            imax=i;
		}

        //if distance to Pareto set 3rd son is less than biggest distance
        //...in population then make substitution
		if (L3<Lmax) {
			for (i=0; i<=lchr-1; i++) { PopChrStr[imax,i]=Son3s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax,i]=Son3p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax,i]=Fu3[i]; }
		}
        //calculating all distances for population		  
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
		//////
		
        //mutation for 4th son		
		if (random<pmut) { 
		Son4p[random(p*(c+d))]=random(2);
        GenVar(Son4s[random(lchr)]);
		}
		
        //functional for 4th son
        SetPsi(Psi0);		  
		for (j=0; j<=lchr-1; j++) { Variations(Son4s[j]); }
        GreytoVector(Son4p);
        Func0(Fu4);
        //Distance for 4th son
        L4=Rast(Fu4);
        //Chromosome with biggest distance to Pareto set
        Lmax=Lh[0];
        imax=0;
		  
		for (i=1;i<=HH-1; i++) {
			if (Lh[i]>Lmax) {
			Lmax=Lh[i];
            imax=i;
			}
		}
		
        //if distance to Pareto set 4th son is less than biggest distance
        //...in population then make substitution
		if (L4<Lmax) {
			for (i=0; i<=lchr-1; i++) {  PopChrStr[imax,i]=Son4s[i]; }
			for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[imax,i]=Son4p[i]; }
			for (i=0; i<=nfu-1; i++) { Fuh[imax,i]=Fu4[i]; }
		}
		
        //calculating all distances for population	
		for (i=0; i<=HH-1; i++) { Lh[i]=Rast(Fuh[i]); }
      rt++;
      //End of cycle for crossoving
    until rt>RR;
	
	// generating new chromosomes
    // Checking Epoch
    pt++;
	
	    //if epoch is over then changing basic
    
	if (pt%Epo==0) 
    {
      //... на наиболее близкую хромосому к утопической
      // хромосоме в пространстве нормированных критериев
	  for (i=0; i<=nfu-1; i++) {
		  Fumax=Fuh[0,i];
        Fumin=Fuh[0,i];
		 // ищем максимальное и минимальное значения по каждому функционалу
		 for (k=0; k<=HH-1; k++) {
			 if (Fuh[k,i]>Fumax) { Fumax=Fuh[k,i]; } else { if (Fuh[k,i]< Fumin) { Fumin=Fuh[k,i]; } }
		 }
		 // нормируем критерии, поделив каждое значение на разность между
        // максимумом и минимумом
		if (Fumax!=Fumin) {
			for (k=0; k<=HH-1; k++) { FuhNorm[k,i]=Fuh[k,i]/(Fumax-Fumin); }
		}
	  }
	  
      // находим хромосому с наименьшей величиной нормы нормированных критериев
      k=0;
      suGA=0;	
	for (i=0; i<=nfu-1; i++) { suGA=suGA+Shtraf[i]*sqr(FuhNorm[0,i]); }
      suGA=sqrt(suGA);
//      su=FuhNorm[0,1];
	  for (i=1; i<=HH-1; i++) {
		  su1GA=0;
		  for (j=0; j<=nfu-1; j++) { su1GA=su1GA+Shtraf[j]*sqr(FuhNorm[i,j]); }
		  su1GA=FuhNorm[i,1];
		  if (su1GA<suGA) {
			  suGA=su1GA;
          k=i;
		  }
	  }
	  
      // заменяем базис
      // строим матрицу для найденной хромосомы
      SetPsi(Psi0);
	for (j=0; j<=lchr-1; j++) { Variations(PopChrStr[k,j]); }
      // меняем базисную матрицу на новую
      SetPsiBas(Psi);
      //генерируем тождественную хромосому
	for (i=0; i<=lchr-1; i++) { for (j=0; j<=3; j++) { PopChrStr[0,i,j]=0; } }
	for (i=0; i<=(p*(c+d)-1); i++) { PopChrPar[0,i]=PopChrPar[k,i]; }
      //вычисляем все фунционалы для всей популяции	  
	  for (i=0; i<=HH-1; i++) { 
	  SetPsi(Psi0);
	  for (j=0; j<=lchr-1; j++) { Variations(PopChrStr[i,j]); }
	  GreytoVector(PopChrPar[i]);
        Func0(Fuh[i]);
	  }
	  
      // формируем элиту
	  for (i=0; i<=kel-1; i++) {
		j=random(HH-1)+1;
        GreytoVector(PopChrPar[j]);
        ImproveChrom(q,PopChrStr[j]); 
	  }
	  
      //вычисляем новые расстояния
	for (i=0; i<=HH-1; i++) {  Lh[i]=Rast(Fuh[i]); }
	}
    //конец цикла поколений

until pt>PP;
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
	  for (j=0; j<=nfu-1; j++) {
		//  fout << " F" << j << "=" << doubeVecToStr(Fuh[Pareto[i],j]) << ", ";
		  
	  }
  }
  
  fout.close();

		 }

///////////////
// Удаление массивов
		 
		 
		 
//конец
return 0;
}

