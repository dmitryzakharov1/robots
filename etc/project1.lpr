program project1;
//types
type
TArrInt=array of integer;
TArrArrInt=array of TArrInt;
TArrArrArrInt=array of TArrArrInt;
TArr5Int=array [0..4]of integer;
TArrArr5Int=array of TArr5Int;
TArrArrArr5int=array of TArrArr5Int;
TArrReal=array of real;
TArrArrReal=array of TArrReal;

const
   deltt=2.5;
   infinity=1e8;
  eps=1e-8;
  eps1=1e-2;
  pokmax=16;
  epstermc=1e-2;
  eps2c=0.5e-2;
  x0c:array [0..17] of double=(0,2,0, 0,4,0, 0,6,0, 0,0,0, 0,0,0, 0,0,0);
  xfc:array[0..8] of real=(-5,0,0, 0,0,0, 5,0,0);
  uminc:array[0..5] of double=(-5,-1, -5,-1, -5,-1);
  umaxc:array[0..5] of double=(5,1, 5,1, 5,1);
  O1sc:array [0..19] of integer=(1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16, 17,18,19,23);
  O2sc:array [0..1] of integer=(1,2);
  qc:array[0..8] of real=(1,1,1, 1,1,1, 1,1,1);
  qminc:array[0..8]of real=(-10,-10,-10, -10,-10,-10, -10,-10,-10);
  qmaxc:array[0..8]of real=(10,10,10, 10,10,10, 10,10,10);
  qyminc:array [0..2]of real =(-10,-10,-10);
  qymaxc:array [0..2]of real =(10,10,10);
  stepsqyc:array[0..2]of real=(20,20,20);
  Pnumc:array [0..8] of integer=(0,1,2, 3,4,5, 6,7,8);
  Rnumc:array [0..8] of integer=(9,10,11, 12,13,14, 15,16,17);
  prepc:array[0..1,0..4,0..1]of real=
      (((-8,1),(-20,1),(-20,-1),(-8,-1),(-8,1)),
       ((20,1),(8,1),(8,-1),(20,-1),(20,1)));
  tfc:real=8;//6;// время окончания
  dtc:real=0.01;//шаг интегрирования

  PsiMultBasc:array [0..3,0..15,0..15] of integer=
       (((0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),   {1}
         (0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0),

         (0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)),


        ((0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0), {2}
         (0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0),

         (0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)),


        ((0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),  {3}
         (0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,2,0, 0,1,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,2, 0,0,1,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 2,0,0,1, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0),

         (0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)),


        ((0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0), {4}
         (0,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,0,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,0),
         (0,0,0,0, 0,0,0,1, 0,0,1,0, 0,0,0,0),

         (0,0,0,0, 0,0,0,0, 1,0,0,1, 0,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,1,0,0, 1,0,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,1,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,1,0),

         (0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,1),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,1,0,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0),
         (0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)));

   M_Entrc: array [0..3,0..11] of integer = ((0,0,0,9, 0,1,0,10, 0,2,0,11),
                                             (0,3,0,12, 0,4,0,13, 0,5,0,14),
                                             (0,6,0,15, 0,7,0,16, 0,8,0,17),
                                             (1,14,1,15, 2,14,2,15, 3,14,3,15));

   M_Outc: array[0..8,0..1] of integer = ((1,15),
                                          (1,13),
                                          (4,15),
                                          (2,15),
                                          (2,13),
                                          (4,14),
                                          (3,15),
                                          (3,13),
                                          (4,13));

 var
   // Parameters of NOP
  kL:integer=4; // number of network operators (layers) in MNOP
  kInp:integer=6;// number of inputs in each layer
  L:integer=16;   // //dimension of network operator matrices
  kp:integer=9;   //size of the set of variables
  kr:integer=9;   //size of the set of parameters
  kw:integer=20;   //size of the set of unary operations
  kv:integer=2;   //size of the set of binary operations
  Mout:integer=9; // number of outputs
  Vs:array of real;//set of variables
  Cs:array of real;//set of parameters
  O1s:array of integer;//set of unary operations
  O2s:array of integer;//set of binary operations
  kOut:integer;// quantity of exit in each matrix
  M_Entr:array of array of integer; // matrix of connects
  V_Entr:array of real; //vector of inputs
  Pnum:array of integer;//vector of number of positions in input vector for variables
  Rnum:array of integer;//vector of number of positions in input vector for parameters
  M_Out:array of array of integer;//matrix of outputs
  z:array of array of real;//vector of nodes
  zs:array of array of string;//string for mathematical expression
  Psi,Psi0:array of array of array of integer;//Network  operator matrices

  // Parameters of GA
  HH:integer=1024; //number of chromosomes in the initial population
  PP:integer=64;  // number of generations (loops of GA)
  RR:integer=64; // number of couples in one generation
  lchr:integer=8; //length of structural part of chromosome
  Epo:integer=24;  // number of generations between epochs
  kel:integer=8;  // number of elitaring chromosomes
  alfa:real=0.4;  // parameter to calculate probability of crossover
  pmut:real=0.7;  // probability of mutation
  PopChrStr:array of array of array [0..4] of integer;//array for structural parts of chromosomes
  PopChrPar:array of array of integer;//array for parametrical parts of chromosomes

  nfu:integer=2;  // number of functionals
  Fuh:array of array of real;// values of functionals for each chromosome
  FuhNorm:array of array of real;
  Shtraf:array of real;
  Lh:array of integer;// values distance to Pareto set
  Pareto:array of integer;// Pareto set

  Son1s,Son2s,Son3s,Son4s:array of array [0..4] of integer;//structural part of sons
  Son1p,Son2p,Son3p,Son4p:array of integer;//parametrical part of sons
  L1,L2,L3,L4:integer;//values distance to Pareto set for sons
  Fu1,Fu2,Fu3,Fu4:array of real;// values of functionals for sons

  p:integer=9;    //number of serching parameters
  c:integer=4;    //number of bits for integer part
  d:integer=12;    //number of bits for fractional part
  q:array of real;//vector of parameters
  qmax,qmin:array of real;// vectors of constraints for parameters
  zb:array of integer;//additional vector

  // Undefined parameters
  ny:integer=3;// dimention of vector of undefine parameters
  qy:array of real;//vextor of undefined parameters
  qymax,qymin:array of real;
  stepsqy:array of real;//vector of steps of undefine parameters

  // Parameters of model
  kRob:integer=3;// number of robots;
  n:integer=18;//dimention of system;
  m:integer=6;//dimention of control;
  ll:integer=42;//dimention of viewing;

  x:array of real;// vector of state
  xs:TArrReal;
  x0:array of real;// vector of initial condition
  xf:array of real;
  fb:array of real;
  fa:array of real;
  u:array of real;// vector of control
  umin:array of real;// vector of control
  umax:array of real;// vector of control
  lv:integer;
  t:real;
  dt:real=0.01; //step of integration;
  tf:real=8;//terminal time;
  dtp:real=0.1;//step of print;

  kChoose:integer;// number of choose chromosome
  f0,f1:real;
  xf1:array of real;
  Prior:array of integer;
  FlagStop:array of boolean;
  epsterm:real=0.01;
  ix,iy:integer;
  EAix,EAixmax:TArrInt;
  sumixmax:real;

  i,j,k, pt,rt,k1,k2,lmax,imax,ks1,ks2:integer;
  ksi,su, su1, Fumax,Fumin:real;

  //************************  1.   *******************************
//ОПИСАНИЕ ВЛОЖЕННЫХ ПРОЦЕДУР И ФУНКЦИЙ
//**************************************************************
Procedure LexPM(var EAix: TArrInt; var flag: boolean);
var
  i,j:integer;
Begin
  i:=ny-1;
  while (i>=0)and(EAix[i]=EAixmax[i]) do i:=i-1;
  if i>=0 then
  begin
    EAix[i]:=EAix[i]+1;
    for j:=i+1 to ny-1 do
      EAix[j]:=0;
    flag:=true;
  end
  else
    flag:=false;
End;
//**************************************************************
Procedure SetEAixmax(EAix1: TArrInt);
var
  i:integer;
  flag:boolean;
Begin
  sumixmax:=0;
  for i:=0 to ny-1 do
    EAixmax[i]:=EAix1[i];
  for i:= 0 to ny-1 do
    EAix[i]:=0;
  repeat
    sumixmax:=sumixmax+1;
    for i:= 0 to ny-1 do
      sumixmax:=sumixmax+EAix[i];
    LexPM(EAix,flag);
  until not flag;
End;
//*************************************************************
Procedure SetPsi(Psi1: TArrArrArrInt);
var
  i,j,k:integer;
Begin
  for k:=0 to kL-1 do
    for i:=0 to L-1 do
      for j:= 0 to L-1 do
        Psi[k,i,j]:=Psi1[k,i,j];
End;
//*************************************************************
  Procedure SetPsiBas(Psi1: TArrArrArrInt);
  var
    k,i,j:integer;
  Begin
    for k:=0 to kL-1 do
      for i:=0 to L-1 do
        for j:= 0 to L-1 do
          Psi0[k,i,j]:=Psi1[k,i,j];
  End;
//*************************************************************
Procedure VectortoGrey(var y: array of integer);
  var
    x,i,j,k:integer;
    r,g1:real;
  Begin
    g1:=1;
    for i := 0 to c - 1 do
      g1:=g1*2;
    for i:=0 to p-1 do
      q[i]:=(q[i]-qmin[i])*g1/(qmax[i]-qmin[i]);
    for i:=0 to p*(c+d)-1 do
      zb[i]:=0;
    for j:=0 to p-1 do
    begin
      x:=trunc(q[j]);
      r:=q[j]-x;
      k:=c+j*(c+d)-1;
      while k>=j*(c+d) do
      begin
        zb[k]:=x mod 2;
        x:=x div 2;
        k:=k-1;
      end;
      k:=c+j*(c+d);
      while k<(c+d)*(j+1) do
      begin
        r:=2*r;
        x:=trunc(r);
        zb[k]:=x;
        r:=r-x;
        k:=k+1;
      end;
      y[j*(c+d)]:=zb[j*(c+d)];
      for i:=j*(c+d)+1 to (j+1)*(c+d)-1 do
        y[i]:=zb[i] xor zb[i-1];
    end;
End;
//*************************************************************
Procedure GenVar(var w:TArr5Int);
  // Генерация элементарной операции
  Function TestSource(j:integer):boolean;
  // если j-номер узла источника, то возвращает false
  var
    i:integer;
    flag:boolean;
  Begin
    flag:=true;
    i:=0;
    while(i<=high(Pnum)) and (j<>Pnum[i]) do i:=i+1;
    if i<=high(Pnum) then flag:=false
    else
    begin
      i:=0;
      while(i<=high(Rnum)) and (j<>Rnum[i]) do i:=i+1;
      if i<=high(Rnum) then flag:=false;
    end;
    result:=flag;
  End;
Begin
  w[0]:=random(kL); //номер слоя
//  w[0]:=random(2)+4;
  w[1]:=random(4);
  case w[1] of
    0,2,3: // замена недиагонального элемента, добавление и удаление дуги
    begin
     w[2]:=random(L-1);
     w[3]:=random(L-w[2]-1)+w[2]+1;
     w[4]:=O1s[random(kW)];
    end;
    1: // замена диагонального элемента
    begin
      w[2]:=random(L);
      w[3]:=w[2];
      w[4]:=O2s[random(kV)];
    end;
  end;
End;
//*************************************************************
Procedure Variations(w:TArr5Int);
// Элементарные операции
// 0 - замена недиагонального элемента
// 1 - замена диагонального элемента
// 2 - добавление дуги
// 3 - удаление дуги
var
  i,j,s1,s2:integer;
Begin
  if (w[1]<>0)or(w[2]<>0)or(w[3]<>0) then
    case w[1] of
      0: if Psi[w[0],w[2],w[3]]<>0 then Psi[w[0],w[2],w[3]]:=w[4];
      1: if Psi[w[0],w[2],w[2]]<>0 then Psi[w[0],w[2],w[2]]:=w[4];
      2: if Psi[w[0],w[2],w[3]]=0 then
           if Psi[w[0],w[3],w[3]]<>0 then Psi[w[0],w[2],w[3]]:=w[4];
      3:
      begin
        s1:=0;
        for i:=0 to w[3]-1 do
          if Psi[w[0],i,w[3]]<>0 then s1:=s1+1;
        s2:=0;
        for j:=w[2]+1 to L-1 do
          if (Psi[w[0],w[2],j]<>0)then s2:=s2+1;
        if s1>1 then
          if s2>1 then
            Psi[w[0],w[1],w[2]]:=0;
      end;
    end;
End;
//*************************************************************
Procedure SetV_Entr;
  var
    i:integer;
  Begin
    for i := 0 to kP-1 do
      V_Entr[Pnum[i]]:=Vs[i];
    for i := 0 to kR-1 do
      V_Entr[Rnum[i]]:=Cs[i];
  End;
//*************************************************************
Function Ro_1(z:real):real;
Begin
  result:=z;
End;
//*************************************************************
Function Ro_2(z:real):real;
Begin
  if abs(z)>sqrt(infinity) then result:=infinity
  else result:=sqr(z);
End;
//*************************************************************
Function Ro_3(z:real):real;
Begin
  result:=-z;
End;
//*************************************************************
Function Ro_10(z:real):real;
Begin
   if z>=0 then result:=1
   else
      result:=-1;
End;
//*************************************************************
Function Ro_4(z:real):real;
Begin
  result:=Ro_10(z)*sqrt(abs(z));
End;
//*************************************************************
Function Ro_5(z:real):real;
Begin
  if abs(z)>eps then result:=1/z
    else result:=Ro_10(z)/eps;
End;
//*************************************************************
Function Ro_6(z:real):real;
Begin
  if z>-ln(eps) then result:=-ln(eps)
  else result:=exp(z);
End;
//*************************************************************
Function Ro_7(z:real):real;
Begin
  if abs(z)<exp(-pokmax) then result:=ln(eps)
    else result:=ln(abs(z));
End;
//*************************************************************
Function Ro_8(z:real):real;
Begin
  if abs(z)>-ln(eps) then
    result:=Ro_10(z)
  else
    result:=(1-exp(-z))/(1+exp(-z));
End;
//*************************************************************
Function Ro_9(z:real):real;
Begin
  if z>=0 then result:=1
    else result:=0;
End;
//*************************************************************
Function Ro_11(z:real):real;
Begin
  result:=cos(z);
End;
//*************************************************************
Function Ro_12(z:real):real;
Begin
  result:=sin(z);
End;
//*************************************************************
Function Ro_13(z:real):real;
Begin
  result:=arctan(z);
End;
//*************************************************************
Function Ro_15(z:real):real;
Begin
  if abs(z)<eps then result:=Ro_10(z)*eps
  else
    result:=Ro_10(z)*exp(ln(abs(z))/3);
End;
//*************************************************************
Function Ro_14(z:real):real;
Begin
  if abs(z)>Ro_15(infinity) then result:=Ro_10(z)*infinity
  else result:=sqr(z)*z;
End;

//*************************************************************
Function Ro_16(z:real):real;
Begin
  if abs(z)<1 then result:=z
  else result:=Ro_10(z);
End;
//*************************************************************
Function Ro_17(z:real):real;
Begin
  result:=Ro_10(z)*ln(abs(z)+1);
End;
//*************************************************************
Function Ro_18(z:real):real;
Begin
  if abs(z)>-ln(eps) then
    result:=Ro_10(z)*infinity
  else
    result:=Ro_10(z)*(exp(abs(z))-1);
End;
//*************************************************************
Function Ro_19(z:real):real;
Begin
  if abs(z)>1/eps then result:=Ro_10(z)*eps
  else result:=Ro_10(z)*exp(-abs(z));
End;
//*************************************************************
Function Ro_20(z:real):real;
Begin
  Result:=z/2;
End;
//*************************************************************
Function Ro_21(z:real):real;
Begin
  Result:=2*z;
End;
//*************************************************************
Function Ro_22(z:real):real;
Begin
  if abs(z)>-ln(eps) then result:=0
  else  result:=exp(abs(z));
End;
//*************************************************************
Function Ro_23(z:real):real;
Begin
  if abs(z)>1/eps then result:=-Ro_10(z)/eps
  else  result:=z-z*sqr(z);
End;
//*************************************************************
Function Ro_24(z:real):real;
Begin
  if z>-ln(eps) then result:=eps/(1+eps)
  else  result:=1/(1+exp(-z));
End;
//*************************************************************
Function Ro_25(z:real):real;
Begin
  if z>0 then result:=1
    else result:=0;
End;
//*************************************************************
Function Ro_26(z:real):real;
Begin
  if abs(z)<eps1 then
    result:=0
  else
    result:=Ro_10(z);
End;
//*************************************************************
Function Ro_27(z:real):real;
Begin
  if abs(z)>1 then
    result:=Ro_10(z)
  else
    result:=Ro_10(z)*(1-sqrt(1-sqr(z)));
End;
//*************************************************************
Function Ro_28(z:real):real;
Begin
  if z*z>ln(infinity) then
    result:=z*(1-eps)
  else
    result:=z*(1-exp(-sqr(z)));
End;
//*************************************************************
Function Xi_0(z1,z2:real):real;
Begin
  result:=z1+z2;
End;
//*************************************************************
Function Xi_1(z1,z2:real):real;
Begin
  result:=z1*z2;
End;
//*************************************************************
Function Xi_2(z1,z2:real):real;
Begin
  if z1>=z2 then result:=z1
  else result:=z2;
End;
//*************************************************************
Function Xi_3(z1,z2:real):real;
Begin
  if z1<z2 then result:=z1
  else result:=z2;
End;
//*************************************************************
Function Xi_4(z1,z2:real):real;
Begin
  result:=z1+z2-z1*z2;
End;
//*************************************************************
Function Xi_5(z1,z2:real):real;
Begin
  result:=Ro_10(z1+z2)*sqrt(sqr(z1)+sqr(z2));
End;
//*************************************************************
Function Xi_6(z1,z2:real):real;
Begin
  result:=Ro_10(z1+z2)*(abs(z1)+abs(z2));
End;
//*************************************************************
Function Xi_7(z1,z2:real):real;
Begin
  result:=Ro_10(z1+z2)*Xi_2(abs(z1),abs(z2));
End;
//*************************************************************
  Procedure RPControl;
  var
    k,i,j:integer;
    zz:real;
  Begin
    SetV_Entr;
    for k := 0 to kL - 1 do
    begin
      for i:=0 to L-1 do
        case psi[k,i,i] of
          2: z[k,i]:=1;
          3: z[k,i]:=-infinity;
          4: z[k,i]:=infinity
          else
          z[k,i]:=0;
        end;
      for i:=0 to kInP-1 do
        if M_Entr[k,2*i]=0 then
          z[k,i]:=V_Entr[M_Entr[k,2*i+1]]
        else
          z[k,i]:=z[M_Entr[k,2*i]-1,M_Entr[k,2*i+1]];
      for i:=0 to L-2 do
        for j:=i+1 to L-1 do
          if Psi[k,i,j]<>0 then
            if Psi[k,j,j]<>0 then
            begin
              case Psi[k,i,j] of
                1: zz:=Ro_1(z[k,i]);
                2: zz:=Ro_2(z[k,i]);
                3: zz:=Ro_3(z[k,i]);
                4: zz:=Ro_4(z[k,i]);
                5: zz:=Ro_5(z[k,i]);
                6: zz:=Ro_6(z[k,i]);
                7: zz:=Ro_7(z[k,i]);
                8: zz:=Ro_8(z[k,i]);
                9: zz:=Ro_9(z[k,i]);
                10: zz:=Ro_10(z[k,i]);
                11: zz:=Ro_11(z[k,i]);
                12: zz:=Ro_12(z[k,i]);
                13: zz:=Ro_13(z[k,i]);
                14: zz:=Ro_14(z[k,i]);
                15: zz:=Ro_15(z[k,i]);
                16: zz:=Ro_16(z[k,i]);
                17: zz:=Ro_17(z[k,i]);
                18: zz:=Ro_18(z[k,i]);
                19: zz:=Ro_19(z[k,i]);
                20: zz:=Ro_20(z[k,i]);
                21: zz:=Ro_21(z[k,i]);
                22: zz:=Ro_22(z[k,i]);
                23: zz:=Ro_23(z[k,i]);
                24: zz:=Ro_24(z[k,i]);
                25: zz:=Ro_25(z[k,i]);
                26: zz:=Ro_26(z[k,i]);
                27: zz:=Ro_27(z[k,i]);
                28: zz:=Ro_28(z[k,i]);
              end;
              case Psi[k,j,j] of
                1: z[k,j]:=Xi_0(z[k,j],zz);
                2: z[k,j]:=Xi_1(z[k,j],zz);
                3: z[k,j]:=Xi_2(z[k,j],zz);
                4: z[k,j]:=Xi_3(z[k,j],zz);
                5: z[k,j]:=Xi_4(z[k,j],zz);
                6: z[k,j]:=Xi_5(z[k,j],zz);
                7: z[k,j]:=Xi_6(z[k,j],zz);
                8: z[k,j]:=Xi_7(z[k,j],zz);
              end;
            end;
    end;
  End;
//*************************************************************
  Procedure Func0(var Fu:array of real);
  var
    i:integer;
  Begin
    RPControl;
    for i:=0 to nfu-1 do
      Fu[i]:=i;
  End;
//*************************************************************
Function Rast(Fu: array of real): integer;
var i,j,k,count:integer;
Begin
  count:=0;
  for i:=0 to HH-1 do
  begin
    j:=0;
    while (j<nfu) and (Fu[j]>=Fuh[i,j]) do j:=j+1;
    if j>=nfu then
    begin
      k:=0;
      while (k<nfu) and (Fu[k]=Fuh[i,k]) do k:=k+1;
      if k<nfu then count:=count+1;
    end;
  end;
  result:=count;
End;
//*************************************************************
Procedure ChoosePareto;
  var
    i,j:integer;
  Begin
    j:=0;
    for i:=0 to HH-1 do
      if Lh[i]=0 then
      begin
        j:=j+1;
        setlength(Pareto,j);
        Pareto[j-1]:=i;
      end;
  End;
//*************************************************************
Procedure GreytoVector(y: TArrInt);
var
  i,j,l1,l:integer;
  g,g1:real;
Begin
  l:=c+d;
  l1:=high(y)+1;
  for i:=0 to l1-1 do
    if i mod l=0 then
      zb[i]:=y[i]
    else
      zb[i]:=zb[i-1] xor y[i];
  j:=-1;
  g1:=1;
  g:=1;
  for i:=0 to c-2 do
    g1:=g1*2;
  for i:=0 to l1-1 do
  begin
    if i mod l=0 then
    begin
      j:=j+1;
      q[j]:=0;
      g:=g1;
    end;
    q[j]:=q[j]+g*zb[i];
    g:=g/2;
  end;
  g1:=g1*2;
  for j := 0 to p - 1 do
    q[j]:=(qmax[j]-qmin[j])*q[j]/g1+qmin[j];
End;
//*************************************************************
Procedure ImproveChrom(q:TArrReal;var StrChrom: TArrArr5Int);
var
  i,j,k:integer;
  flag:boolean;
Begin
  SetPsi(Psi0);
  Func0(Fu1);
  k:=-1;
  for i:=0 to lchr-1 do
  begin
    Variations(StrChrom[i]);
    Func0(Fu2);
    flag:=true;
    for j:=0 to nfu-1 do
      if Fu2[j]>Fu1[j] then flag:=false;
    if flag then
    begin
      for j:=0 to nfu-1 do
        Fu1[j]:=Fu2[j];
      k:=i;
    end;
  end;
  for i:=k+1 to lchr-1 do
    for j:=0 to 3 do
      StrChrom[i,j]:=0;
End;

//*************************************************************


BEGIN
  //************************  2.  *******************************
  //Определяем размерности и задаем значения параметров алгоритма
  //*************************************************************

  //номера унарных операций
  //здесь начинается выполнение
  SetLength(O1s,kW);
  for i:=0 to kW-1 do
      O1s[i]:=O1sc[i];
  //номера бинарных операций
  SetLength(O2s,kV);
  for i:=0 to kV-1 do
      O2s[i]:=O2sc[i];
  //базисная матрица
  SetLength(Psi,kL,L,L);
  SetLength(Psi0,kL,L,L);
    for k:=0 to kL-1 do
          for i := 0 to L - 1 do
            for j := 0 to L-1 do
                begin
                     Psi[k,i,j]:=PsiMultBasc[k,i,j];
                     Psi0[k,i,j]:=PsiMultBasc[k,i,j];
                end;
  //матрица входов
  SetLength(M_Entr,kL,2*kInp);
  for k := 0 to kL - 1 do
        for i:=0 to 2*kinp - 1 do
          M_Entr[k,i]:=M_Entrc[k,i];
  //матрица выходов
  SetLength(M_Out,Mout,2);
  for i := 0 to Mout - 1 do
        for j:=0 to 1 do
           M_Out[i,j]:=M_Outc[i,j];
  //параметры модели
  dt:=dtc;
  tf:=tfc;
  SetLength(qmax,p);
  SetLength(qmin,p);
  SetLength(q,p);
  for i:=0 to p-1 do
  begin
    qmax[i]:=qmaxc[i];
    qmin[i]:=qminc[i];
    q[i]:=qc[i];
  end;
  //неопределенные параметры
  Setlength(qymax,ny);
  Setlength(qymin,ny);
  Setlength(stepsqy,ny);
  SetLength(qy,ny);
  SetLength(EAixmax,ny);
  SetLength(EAix,ny);
  for i:= 0 to ny-1 do
  begin
    qymin[i]:=qyminc[i];
    qymax[i]:=qymaxc[i];
    stepsqy[i]:=stepsqyc[i];
    qy[i]:=qymin[i];
  end;
  for i:=0 to ny-1 do
      EAix[i]:=trunc((qymax[i]-qymin[i])/stepsqy[i]);
  SetEAixmax(EAix);
  //вектор состояния объекта
  Setlength(x0,n);
  SetLength(x,n);
  SetLength(xs,n);
  SetLength(fa,n);
  SetLength(fb,n);
  Setlength(xf,high(xfc)+1);
  for i:=0 to n-1 do
      x0[i]:=x0c[i];
  for i:=0 to high(xf) do
      xf[i]:=xfc[i];

  //ограничения на управление
  SetLength(u,m);
  Setlength(umin,m);
  Setlength(umax,m);
  for i:=0 to m-1 do
      umin[i]:=uminc[i];
  for i:=0 to m-1 do
      umax[i]:=umaxc[i];

  //0-ой слой (переменных и параметров)
  SetLength(V_Entr,kP+kR);
  SetLength(Pnum,kP);
  SetLength(Rnum,kR);
  SetLength(Vs,kP);
  SetLength(Cs,kR);
  for i:=0 to kP-1 do
      Pnum[i]:=Pnumc[i];
  for i:=0 to kR-1 do
      Rnum[i]:=Rnumc[i];
  for i:=0 to kR-1 do
      Cs[i]:=q[i];
  SetV_Entr;
  SetLength(Prior,kRob);
  SetLength(FlagStop,kRob);
  for i := 0 to kRob-1 do
  begin
    FlagStop[i]:=true;
    Prior[i]:=kRob-i;
  end;
  Setlength(PopChrStr,HH,lchr);
  Setlength(PopChrPar,HH,p*(c+d));
  Setlength(Fuh,HH,nfu);
  SetLength(FuhNorm,HH,nfu);
  Setlength(Shtraf,nfu);
  for i := 0 to nfu - 1 do
    Shtraf[i]:=1;
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

  SetPsiBas(Psi);
  //generating population
  VectortoGrey(PopChrPar[0]);
  for i:=0 to lchr-1 do
    for j:=0 to 4 do
      PopChrStr[0,i,j]:=0;
  for i:=1 to HH-1 do
  begin
    for j:=0 to lchr-1 do
      GenVar(PopChrStr[i,j]);
    for j:=0 to p*(c+d)-1 do
      PopChrPar[i,j]:=random(2);
  end;
  // calculating values of functionals
  for i:=0 to HH-1 do
  begin
    SetPsi(Psi0);
    for j:=0 to lchr-1 do
    Variations(PopChrStr[i,j]);
    GreytoVector(PopChrPar[i]);
    Func0(Fuh[i]);
  end;
  //caculating distances to Pareto set
  for i:=0 to HH-1 do
    Lh[i]:=Rast(Fuh[i]);
  //Start of cycle for generations
  pt:=1;  // first current generation
  repeat
    //start of cycle for crossovering
    rt:=1;//first couple for crossoving
    repeat
      //select of two parents
      k1:=random(HH);
      k2:=random(HH);
      ksi:=random;
      if (ksi<(1+alfa*Lh[k1])/(1+Lh[k1])) or
         (ksi<(1+alfa*Lh[k2])/(1+Lh[k2])) then
      begin
        //if true
        ks1:=random(lchr);
        ks2:=random(p*(c+d));
        //crossoving? creating four sons
        for i:=0 to lchr-1 do
        begin
          Son1s[i]:=PopChrStr[k1,i];
          Son2s[i]:=PopChrStr[k2,i];
        end;
        for i:=0 to ks2-1 do
        begin
          Son1p[i]:=PopChrPar[k1,i];
          Son2p[i]:=PopChrPar[k2,i];
          Son3p[i]:=PopChrPar[k1,i];
          Son4p[i]:=PopChrPar[k2,i];
        end;
        for i:=ks2 to p*(c+d)-1 do
        begin
          Son1p[i]:=PopChrPar[k2,i];
          Son2p[i]:=PopChrPar[k1,i];
          Son3p[i]:=PopChrPar[k2,i];
          Son4p[i]:=PopChrPar[k1,i];
        end;
        for i:=0 to ks1-1 do
        begin
          Son3s[i]:=PopChrStr[k1,i];
          Son4s[i]:=PopChrStr[k2,i];
        end;
        for i:=ks1 to lchr-1 do
        begin
          Son3s[i]:=PopChrStr[k2,i];
          Son4s[i]:=PopChrStr[k1,i];
        end;
        //mutation for 1st son
        if random<pmut then
        begin
          son1p[random(p*(c+d))]:=random(2);
          GenVar(son1s[random(lchr)]);
        end;
        //functional for 1st son
        SetPsi(Psi0);
        for j:=0 to lchr-1 do
          Variations(son1s[j]);
        GreytoVector(son1p);
        Func0(Fu1);
        //Distance for 1st son
        L1:=Rast(Fu1);
        //Chromosome with biggest distance to Pareto set
        Lmax:=Lh[0];
        imax:=0;
        for i:=1 to HH-1 do
          if Lh[i]>Lmax then
          begin
            Lmax:=Lh[i];
            imax:=i;
          end;
        if L1<Lmax then
        //if distance to Pareto set 1st son is less than biggest distance
        //...in population then make substitution
        begin
          for i:=0 to lchr-1 do
            PopChrStr[imax,i]:=son1s[i];
          for i:=0 to p*(c+d)-1 do
            PopChrPar[imax,i]:=son1p[i];
          for i:=0 to nfu-1 do
            Fuh[imax,i]:=Fu1[i];
        end;
        //calculating all distances for population
        for i:=0 to HH-1 do
          Lh[i]:=Rast(Fuh[i]);
        //mutation for 2nd son
        if random<pmut then
        begin
          son2p[random(p*(c+d))]:=random(2);
          GenVar(son2s[random(lchr)]);
        end;
        //functional for 2nd son
        SetPsi(Psi0);
        for j:=0 to lchr-1 do
          Variations(son2s[j]);
        GreytoVector(son2p);
        Func0(Fu2);
        //Distance for 2nd son
        L2:=Rast(Fu2);
        //Chromosome with biggest distance to Pareto set
        Lmax:=Lh[0];
        imax:=0;
        for i:=1 to HH-1 do
          if Lh[i]>Lmax then
          begin
            Lmax:=Lh[i];
            imax:=i;
          end;
        if L2<Lmax then
        //if distance to Pareto set 2nd son is less than biggest distance
        //...in population then make substitution
        begin
          for i:=0 to lchr-1 do
            PopChrStr[imax,i]:=son2s[i];
          for i:=0 to p*(c+d)-1 do
            PopChrPar[imax,i]:=son2p[i];
          for i:=0 to nfu-1 do
            Fuh[imax,i]:=Fu2[i];
        end;
        //calculating all distances for population
        for i:=0 to HH-1 do
          Lh[i]:=Rast(Fuh[i]);
        //mutation for 3rd son
        if random<pmut then
        begin
          son3p[random(p*(c+d))]:=random(2);
          GenVar(son3s[random(lchr)]);
        end;
        //functional for 3rd son
        SetPsi(Psi0);
        for j:=0 to lchr-1 do
          Variations(son3s[j]);
        GreytoVector(son3p);
        Func0(Fu3);
        //Distance for 3rd son
        L3:=Rast(Fu3);
        //Chromosome with biggest distance to Pareto set
        Lmax:=Lh[0];
        imax:=0;
        for i:=1 to HH-1 do
          if Lh[i]>Lmax then
          begin
            Lmax:=Lh[i];
            imax:=i;
          end;
        if L3<Lmax then
        //if distance to Pareto set 3rd son is less than biggest distance
        //...in population then make substitution
        begin
          for i:=0 to lchr-1 do
            PopChrStr[imax,i]:=son3s[i];
          for i:=0 to p*(c+d)-1 do
            PopChrPar[imax,i]:=son3p[i];
          for i:=0 to nfu-1 do
            Fuh[imax,i]:=Fu3[i];
        end;
        //calculating all distances for population
        for i:=0 to HH-1 do
          Lh[i]:=Rast(Fuh[i]);
        //mutation for 4th son
        if random<pmut then
        begin
          son4p[random(p*(c+d))]:=random(2);
          GenVar(son4s[random(lchr)]);
        end;
        //functional for 4th son
        SetPsi(Psi0);
        for j:=0 to lchr-1 do
          Variations(son4s[j]);
        GreytoVector(son4p);
        Func0(Fu4);
        //Distance for 4th son
        L4:=Rast(Fu4);
        //Chromosome with biggest distance to Pareto set
        Lmax:=Lh[0];
        imax:=0;
        for i:=1 to HH-1 do
          if Lh[i]>Lmax then
          begin
            Lmax:=Lh[i];
            imax:=i;
          end;
        if L4<Lmax then
        //if distance to Pareto set 4th son is less than biggest distance
        //...in population then make substitution
        begin
          for i:=0 to lchr-1 do
            PopChrStr[imax,i]:=son4s[i];
          for i:=0 to p*(c+d)-1 do
            PopChrPar[imax,i]:=son4p[i];
          for i:=0 to nfu-1 do
            Fuh[imax,i]:=Fu4[i];
        end;
        //calculating all distances for population
        for i:=0 to HH-1 do
          Lh[i]:=Rast(Fuh[i]);
      end;
      rt:=rt+1;
      //End of cycle for crossoving
    until rt>RR;
    // generating new chromosomes
    // Checking Epoch
    pt:=pt+1;
    //if epoch is over then changing basic
    if pt mod Epo=0 then
    begin
      //... на наиболее близкую хромосому к утопической
      // хромосоме в пространстве нормированных критериев
      for i:=0 to nfu-1 do
      begin
        Fumax:=Fuh[0,i];
        Fumin:=Fuh[0,i];
        // ищем максимальное и минимальное значения по каждому функционалу
        for k:=0 to HH-1 do
          if Fuh[k,i]>Fumax then
            Fumax:=Fuh[k,i]
          else
            if Fuh[k,i]< Fumin then
              Fumin:=Fuh[k,i];
        // нормируем критерии, поделив каждое значение на разность между
        // максимумом и минимумом
        if Fumax<>Fumin then
          for k:=0 to HH-1 do
            FuhNorm[k,i]:=Fuh[k,i]/(Fumax-Fumin);
      end;
      // находим хромосому с наименьшей величиной нормы нормированных критериев
      k:=0;
      su:=0;
      for i:=0 to nfu-1 do
        su:=su+Shtraf[i]*sqr(FuhNorm[0,i]);
      su:=sqrt(su);
//      su:=FuhNorm[0,1];
      for i:=1 to HH-1 do
      begin
        su1:=0;
        for j:=0 to nfu-1 do
          su1:=su1+Shtraf[j]*sqr(FuhNorm[i,j]);
//        su1:=sqrt(su1);
        su1:=FuhNorm[i,1];
        if su1<su then
        begin
          su:=su1;
          k:=i;
        end;
      end;
      // заменяем базис
      // строим матрицу для найденной хромосомы
      SetPsi(Psi0);
      for j:=0 to lchr-1 do
        Variations(PopChrStr[k,j]);
      // меняем базисную матрицу на новую
      SetPsiBas(Psi);
      //генерируем тождественную хромосому
      for i:=0 to lchr-1 do
        for j:=0 to 3 do
          PopChrStr[0,i,j]:=0;
      for i:=0 to p*(c+d)-1 do
        PopChrPar[0,i]:=PopChrPar[k,i];
      //вычисляем все фунционалы для всей популяции
      for i:=0 to HH-1 do
      begin
           SetPsi(Psi0);
        for j:=0 to lchr-1 do
          Variations(PopChrStr[i,j]);
        GreytoVector(PopChrPar[i]);
        Func0(Fuh[i]);
      end;
      // формируем элиту
      for i:=0 to kel-1 do
      begin
        j:=random(HH-1)+1;
        GreytoVector(PopChrPar[j]);
        ImproveChrom(q,PopChrStr[j]);
      end;
      //вычисляем новые расстояния
      for i:=0 to HH-1 do
        Lh[i]:=Rast(Fuh[i]);
    end;
    //конец цикла поколений

  until pt>PP;
  ChoosePareto;
 //строим множество Парето
END.



