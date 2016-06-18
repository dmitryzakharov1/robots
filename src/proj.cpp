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

    int main(int argc,
        const char * argv[]) {
        unsigned int start_time = clock(); // начальное время
        //************************  2.  *******************************
        //Определяем размерности и задаем значения параметров алгоритма
        //*************************************************************

        //номера унарных операций
        O1s.resize(kW);
        for (i = 0; i < kW; i++) {
            O1s[i] = O1sc[i];
        }
        //номера бинарных операций
        O2s.resize(kV);
        for (i = 0; i < kV; i++) {
            O2s[i] = O2sc[i];
        }
        //базисная матрица
        Psi.resize(kL);
        for (int ii = 0; ii < kL; ii++) {
            Psi[ii].resize(L);
            for (int jj = 0; jj < L; jj++) {
                Psi[ii][jj].resize(L);
            }
        }

        Psi0.resize(kL);
        for (int ii = 0; ii < kL; ii++) {
            Psi0[ii].resize(L);
            for (int jj = 0; jj < L; jj++) {
                Psi0[ii][jj].resize(L);
            }
        }

        for (k = 0; k < kL; k++) {
            for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                    Psi[k][i][j] = PsiMultBasc[k][i][j];
                    Psi0[k][i][j] = PsiMultBasc[k][i][j];
                }
            }
        }

        //матрица входов
        M_Entr.resize(kL);
        for (int ii = 0; ii < kL; ii++) {
            M_Entr[ii].resize(2 * kInP);
        }

        for (k = 0; k < kL; k++) {
            for (i = 0; i < (2 * kInP); i++) {
                M_Entr[k][i] = M_Entrc[k][i];
            }
        }

        //матрица выходов
        M_Out.resize(Mout);
        for (int ii = 0; ii < Mout; ii++) {
            M_Out[ii].resize(2);
        }

        for (i = 0; i < Mout; i++) {
            for (j = 0; j < 2; j++) {
                M_Out[i][j] = M_Outc[i][j];
            }
        }

        //параметры модели
        dt = dtc;
        tf = tfc;
        qmax.resize(p);
        qmin.resize(p);
        q.resize(p);
        for (i = 0; i < p; i++) {
            qmax[i] = qmaxc[i];
            qmin[i] = qminc[i];
            q[i] = qc[i];
        }

        //неопределенные параметры
        qymax.resize(ny);
        qymin.resize(ny);
        stepsqy.resize(ny);
        qy.resize(ny);
        EAixmax.resize(ny);
        EAix.resize(ny);

        for (i = 0; i < ny; i++) {
            qymin[i] = qyminc[i];
            qymax[i] = qymaxc[i];
            stepsqy[i] = stepsqyc[i];
            qy[i] = qyminc[i];
        }

        for (i = 0; i < ny; i++) {
            EAix[i] = trunc((qymax[i] - qymin[i]) / stepsqy[i]);
        }
        //cout << "Start SetEAixmax(EAix)" << endl;
        SetEAixmax(EAix);
        //cout << "End of SetEAixmax(EAix)" << endl;

        //вектор состояния объекта
        x0.resize(n);
        x.resize(n);
        xs.resize(n);
        fa.resize(n);
        fb.resize(n);
        xf.resize(9);
        for (i = 0; i < n; i++) {
            x0[i] = x0c[i];
        }
        for (i = 0; i < xf.size(); i++) {
            xf[i] = xfc[i];
        }

        //ограничения на управление
        u.resize(m);
        //Setlength(umin,m);
        umin.resize(m);
        //Setlength(umax,m);
        umax.resize(m);
        for (i = 0; i < m; i++) {
            umin[i] = uminc[i];
        }
        for (i = 0; i < m; i++) {
            umax[i] = umaxc[i];
        }
        //вектор наблюдения
        y.resize(lv);
        //хранение промежуточных значений функционалов
        su.resize(nfu);
        su1.resize(nfu);

        //0-ой слой (переменных и параметров)
        V_Entr.resize(kP + kR);
        Pnum.resize(kP);
        Rnum.resize(kR);
        Vs.resize(kP);
        Cs.resize(kR);

        for (i = 0; i < kP; i++) {
            Pnum[i] = Pnumc[i];
        }
        for (i = 0; i < kR; i++) {
            Rnum[i] = Rnumc[i];
        }
        for (i = 0; i < kR; i++) {
            Cs[i] = q[i];
        }
        SetV_Entr();
        Prior.resize(kRob);
        FlagStop.resize(kRob);
        for (i = 0; i < kRob; i++) {
            FlagStop[i] = true;
            Prior[i] = kRob - i;
        }
        PopChrStr.resize(HH);
        for (int ii = 0; ii < HH; ii++) {
            PopChrStr[ii].resize(lchr);
            for (int jj = 0; jj < lchr; jj++) {
                PopChrStr[ii][jj].resize(5);
            }
        }

        PopChrPar.resize(HH);
        for (int ii = 0; ii < HH; ii++) {
            PopChrPar[ii].resize(p * (c + d));
        }

        Fuh.resize(HH);
        for (int ii = 0; ii < HH; ii++) {
            Fuh[ii].resize(nfu);
        }

        FuhNorm.resize(HH);
        for (int ii = 0; ii < HH; ii++) {
            FuhNorm[ii].resize(nfu);
        }

        Shtraf.resize(nfu);

        for (i = 0; i < nfu; i++) {
            Shtraf[i] = 1;
        }
        Lh.resize(HH);
        Fu1.resize(nfu);
        Fu2.resize(nfu);
        Fu3.resize(nfu);
        Fu4.resize(nfu);
        Son1s.resize(lchr);
        for (int ii = 0; ii < lchr; ii++) {
            Son1s[ii].resize(5);
        }
        Son2s.resize(lchr);
        for (int ii = 0; ii < lchr; ii++) {
            Son2s[ii].resize(5);
        }
        Son3s.resize(lchr);
        for (int ii = 0; ii < lchr; ii++) {
            Son3s[ii].resize(5);
        }
        Son4s.resize(lchr);
        for (int ii = 0; ii < lchr; ii++) {
            Son4s[ii].resize(5);
        }
        Son1p.resize(p * (c + d));
        Son2p.resize(p * (c + d));
        Son3p.resize(p * (c + d));
        Son4p.resize(p * (c + d));
        zb.resize(p * (c + d));
        z.resize(kL);
        for (int ii = 0; ii < kL; ii++) {
            z[ii].resize(L);
        }

        zs.resize(kL);
        for (int ii = 0; ii < kL; ii++) {
            zs[ii].resize(L);
        }

        //************************  2.  *******************************
        //                 ГЕНЕТИЧЕСКИЙ АЛГОРИТМ
        //*************************************************************
        //задаем базисную матрицу сетевого оператора
        //функция копирует Psi1 d Psi, это ненедо
        SetPsiBas(Psi);
        //генерируем начальную (нулевую) популяцию векторов вариаций 
        VectortoGrey(PopChrPar[0]);
        for (i = 0; i < lchr; i++) {
            for (j = 0; j <= 4; j++) {
                PopChrStr[0][i][j] = 0;
            }
        }

        //наполняем вектора случайными возможными значениями
        //PopChrStr нормально заполняется случайными значениями
        for (i = 1; i <= HH - 1; i++) {
            for (j = 0; j <= lchr - 1; j++) {
                GenVar(PopChrStr[i][j]);
            }
            for (j = 0; j <= (p * (c + d) - 1); j++) {
                PopChrPar[i][j] = rand() % 2;
            }
        }

        // считаем значения функционалов для каждой хромосомы
        for (i = 0; i <= (HH - 1); i++) {
            //для структурной части: берем текущую базисную матрицу
            //cout << "Start SetPsi(Psi0)" << endl;
            //функция также скопирует Psi0 в Psi это ненадо
            SetPsi(Psi0);
            //cout << "End of SetPsi(Psi0)" << endl; 
            //вносим в нее вариации;
            for (j = 0; j <= (lchr - 1); j++) {
                Variations(PopChrStr[i][j]);
            }
            //для параметрической:
            //текущие значения параметров в коде Грея переводим в вектор параметров
            //cout << "Start GreytoVector" << endl;
            GreytoVector(PopChrPar[i]);
            //cout << "End of GreytoVector" << endl;

            // для каждой полученной матрицы сетевого оператора и вектора параметров,
            //определяющих искомое математическое выражение, оцениваем
            //функции приспособленности

            //тяжелая функция, надо ускорять
            //cout << "Start Func0(Fuh[i])" << endl;
            Func0(Fuh[i]);
            //cout << "End of Func0(Fuh[i])" << endl;

        }

        //caculating distances to Pareto set

        for (i = 0; i <= HH - 1; i++) {
            Lh[i] = Rast(Fuh[i]);
        }

        //Start of cycle for generations
        pt = 1; // first current generation

        do {
            //start of cycle for crossovering
            rt = 1; //first couple for crossoving
            do {
                //select of two parents
                k1 = rand() % HH;
                k2 = rand() % HH;
                ksi = rand();
                if ((ksi < (1 + alfa * Lh[k1]) / (1 + Lh[k1])) ||
                    (ksi < (1 + alfa * Lh[k2]) / (1 + Lh[k2]))) {

                    //if true
                    ks1 = rand() % lchr;
                    ks2 = rand() % (p * (c + d));
                    //crossover creating four sons

                    for (i = 0; i <= lchr - 1; i++) {
                        Son1s[i] = PopChrStr[k1][i];
                        Son2s[i] = PopChrStr[k2][i];
                    }

                    for (i = 0; i <= ks2 - 1; i++) {
                        Son1p[i] = PopChrPar[k1][i];
                        Son2p[i] = PopChrPar[k2][i];
                        Son3p[i] = PopChrPar[k1][i];
                        Son4p[i] = PopChrPar[k2][i];
                    }

                    for (i = ks2; i <= (p * (c + d) - 1); i++) {
                        Son1p[i] = PopChrPar[k2][i];
                        Son2p[i] = PopChrPar[k1][i];
                        Son3p[i] = PopChrPar[k2][i];
                        Son4p[i] = PopChrPar[k1][i];
                    }

                    for (i = 0; i <= (ks1 - 1); i++) {
                        Son3s[i] = PopChrStr[k1][i];
                        Son4s[i] = PopChrStr[k2][i];
                    }

                    for (i = ks1; i <= lchr - 1; i++) {
                        Son3s[i] = PopChrStr[k2][i];
                        Son4s[i] = PopChrStr[k1][i];
                    }

                    //mutation for 1st son	
                    if (rand() < pmut) {
                        Son1p[rand() % (p * (c + d))] = rand() % 2;
                        GenVar(Son1s[rand() % lchr]);
                    }

                    //functional for 1st son
                    SetPsi(Psi0);
                    //вместо этой функции

                    for (j = 0; j <= lchr - 1; j++) {
                        Variations(Son1s[j]);
                    }
                    GreytoVector(Son1p);
                    Func0(Fu1);
                    //Distance for 1st son
                    L1 = Rast(Fu1);
                    //Chromosome with biggest distance to Pareto set
                    lmax = Lh[0];
                    imax = 0;
                    for (i = 1; i <= HH - 1; i++) {
                        if (Lh[i] > lmax) {
                            lmax = Lh[i];
                            imax = i;
                        }
                    }

                    //if distance to Pareto set 1st son is less than biggest distance
                    //...in population then make substitution
                    if (L1 < lmax) {
                        for (i = 0; i <= lchr - 1; i++) {
                            PopChrStr[imax][i] = Son1s[i];
                        }
                        for (i = 0; i <= (p * (c + d) - 1); i++) {
                            PopChrPar[imax][i] = Son1p[i];
                        }
                        for (i = 0; i <= nfu - 1; i++) {
                            Fuh[imax][i] = Fu1[i];
                        }
                    }

                    //calculating all distances for population		  
                    for (i = 0; i <= HH - 1; i++) {
                        Lh[i] = Rast(Fuh[i]);
                    }
                    //mutation for 2nd son	
                    if (rand() < pmut) {
                        Son2p[rand() % (p * (c + d))] = rand() % 2;
                        GenVar(Son2s[rand() % lchr]);
                    }
                    //functional for 2nd son
                    SetPsi(Psi0);

                    for (j = 0; j <= lchr - 1; j++) {
                        Variations(Son2s[j]);
                    }
                    GreytoVector(Son2p);
                    Func0(Fu2);
                    //Distance for 2nd son
                    L2 = Rast(Fu2);
                    //Chromosome with biggest distance to Pareto set
                    lmax = Lh[0];
                    imax = 0;

                    for (i = 1; i <= HH - 1; i++) {
                        if (Lh[i] > lmax) {
                            lmax = Lh[i];
                            imax = i;
                        }
                    }

                    //if distance to Pareto set 2nd son is less than biggest distance
                    //...in population then make substitution
                    if (L2 < lmax) {
                        for (i = 0; i <= (lchr - 1); i++) {
                            PopChrStr[imax][i] = Son2s[i];
                        }
                        for (i = 0; i <= (p * (c + d) - 1); i++) {
                            PopChrPar[imax][i] = Son2p[i];
                        }
                        for (i = 0; i <= nfu - 1; i++) {
                            Fuh[imax][i] = Fu2[i];
                        }
                    }

                    //calculating all distances for population		  
                    for (i = 0; i <= HH - 1; i++) {
                        Lh[i] = Rast(Fuh[i]);
                    }
                    //mutation for 3rd son
                    if (rand() < pmut) {
                        Son3p[rand() % (p * (c + d)) + 1] = rand() % 2;
                        GenVar(Son3s[rand() % lchr]);
                    }
                    //functional for 3rd son
                    SetPsi(Psi0);

                    for (j = 0; j <= (lchr - 1); j++) {
                        Variations(Son3s[j]);
                    }
                    GreytoVector(Son3p);
                    Func0(Fu3);
                    //Distance for 3rd son
                    L3 = Rast(Fu3);
                    //Chromosome with biggest distance to Pareto set
                    lmax = Lh[0];
                    imax = 0;
                    for (i = 1; i <= HH - 1; i++) {
                        lmax = Lh[i];
                        imax = i;
                    }

                    //if distance to Pareto set 3rd son is less than biggest distance
                    //...in population then make substitution
                    if (L3 < lmax) {
                        for (i = 0; i <= lchr - 1; i++) {
                            PopChrStr[imax][i] = Son3s[i];
                        }
                        for (i = 0; i <= (p * (c + d) - 1); i++) {
                            PopChrPar[imax][i] = Son3p[i];
                        }
                        for (i = 0; i <= nfu - 1; i++) {
                            Fuh[imax][i] = Fu3[i];
                        }
                    }
                    //calculating all distances for population		  
                    for (i = 0; i <= HH - 1; i++) {
                        Lh[i] = Rast(Fuh[i]);
                    }

                    //mutation for 4th son		
                    if (rand() < pmut) {
                        Son4p[rand() % (p * (c + d)) + 1] = rand() % 2;
                        GenVar(Son4s[rand() % lchr]);
                    }

                    //functional for 4th son
                    SetPsi(Psi0);

                    for (j = 0; j <= lchr - 1; j++) {
                        Variations(Son4s[j]);
                    }
                    GreytoVector(Son4p);
                    Func0(Fu4);
                    //Distance for 4th son
                    L4 = Rast(Fu4);
                    //Chromosome with biggest distance to Pareto set
                    lmax = Lh[0];
                    imax = 0;

                    for (i = 1; i <= HH - 1; i++) {
                        if (Lh[i] > lmax) {
                            lmax = Lh[i];
                            imax = i;
                        }
                    }

                    //if distance to Pareto set 4th son is less than biggest distance
                    //...in population then make substitution
                    if (L4 < lmax) {
                        for (i = 0; i <= lchr - 1; i++) {
                            PopChrStr[imax][i] = Son4s[i];
                        }
                        for (i = 0; i <= (p * (c + d) - 1); i++) {
                            PopChrPar[imax][i] = Son4p[i];
                        }
                        for (i = 0; i <= nfu - 1; i++) {
                            Fuh[imax][i] = Fu4[i];
                        }
                    }

                    //calculating all distances for population	
                    for (i = 0; i <= HH - 1; i++) {
                        Lh[i] = Rast(Fuh[i]);
                    }
                }
                rt++;
                //End of cycle for crossoving
            }
            while (!(rt > RR));

            // generating new chromosomes
            // Checking Epoch
            pt++;

            //if epoch is over then changing basic

            if (pt % Epo == 0) {
                //... на наиболее близкую хромосому к утопической
                // хромосоме в пространстве нормированных критериев
                for (i = 0; i <= nfu - 1; i++) {
                    Fumax = Fuh[0][i];
                    Fumin = Fuh[0][i];
                    // ищем максимальное и минимальное значения по каждому функционалу
                    for (k = 0; k <= HH - 1; k++) {
                        if (Fuh[k][i] > Fumax) {
                            Fumax = Fuh[k][i];
                        } else {
                            if (Fuh[k][i] < Fumin) {
                                Fumin = Fuh[k][i];
                            }
                        }
                    }
                    // нормируем критерии, поделив каждое значение на разность между
                    // максимумом и минимумом
                    if (Fumax != Fumin) {
                        for (k = 0; k <= HH - 1; k++) {
                            FuhNorm[k][i] = Fuh[k][i] / (Fumax - Fumin);
                        }
                    }
                }

                // находим хромосому с наименьшей величиной нормы нормированных критериев
                k = 0;
                suGA = 0;
                for (i = 0; i <= nfu - 1; i++) {
                    suGA = suGA + Shtraf[i] * pow((FuhNorm[0][i]), 2);
                }
                suGA = sqrt(suGA);
                //      su=FuhNorm[0,1];
                for (i = 1; i <= HH - 1; i++) {
                    su1GA = 0;
                    for (j = 0; j <= nfu - 1; j++) {
                        su1GA = su1GA + Shtraf[j] * pow((FuhNorm[i][j]), 2);
                    }
                    su1GA = FuhNorm[i][1];
                    if (su1GA < suGA) {
                        suGA = su1GA;
                        k = i;
                    }
                }

                // заменяем базис
                // строим матрицу для найденной хромосомы
                SetPsi(Psi0);

                for (j = 0; j <= lchr - 1; j++) {
                    Variations(PopChrStr[k][j]);
                }
                // меняем базисную матрицу на новую
                //функция присваивает Psi Psi
                //SetPsiBas(Psi);

                //генерируем тождественную хромосому
                for (i = 0; i <= lchr - 1; i++) {
                    for (j = 0; j <= 3; j++) {
                        PopChrStr[0][i][j] = 0;
                    }
                }
                for (i = 0; i <= (p * (c + d) - 1); i++) {
                    PopChrPar[0][i] = PopChrPar[k][i];
                }
                //вычисляем все фунционалы для всей популяции	  
                for (i = 0; i <= HH - 1; i++) {
                    SetPsi(Psi0);

                    for (j = 0; j <= lchr - 1; j++) {
                        Variations(PopChrStr[i][j]);
                    }
                    GreytoVector(PopChrPar[i]);
                    Func0(Fuh[i]);
                }

                // формируем элиту
                for (i = 0; i <= kel - 1; i++) {
                    j = rand() % (HH - 1) + 1;
                    GreytoVector(PopChrPar[j]);
                    ImproveChrom(PopChrStr[j]);
                }

                //вычисляем новые расстояния
                for (i = 0; i <= HH - 1; i++) {
                    Lh[i] = Rast(Fuh[i]);
                }
            }
            //конец цикла поколений
        }

        while (!(pt > PP));
        ChoosePareto();
        //Сохраним множество Парето в файл
        ofstream fout("output.txt");
        kol = Pareto.size();
        //сортируем решения по возрастанию

        for (i = 0; i <= kol - 2; i++) {
            for (j = i + 1; j <= kol - 1; j++) {
                if (Fuh[Pareto[i]][0] > Fuh[Pareto[j]][0]) {
                    k = Pareto[i];
                    Pareto[i] = Pareto[j];
                    Pareto[j] = k;
                }
            }
        }

        for (i = 0; i <= kol - 1; i++) {
            fout << "solution No " << Pareto[i] << ", ";
            for (j = 0; j <= nfu - 1; j++) {
                fout << " F" << j << "=" << Fuh[Pareto[i]][j] << ", ";
            }
            fout << endl;
        }

        fout.close();

        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "Execution time:" << search_time / double(CLOCKS_PER_SEC) << " seconds" << endl; //в секунды переводим и выводим

        //конец
        return 0;
    }