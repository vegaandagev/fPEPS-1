#include"Header.h"
void GSF(cx_mat * TI,cx_mat * TIp,ofstream & Fidelity,ofstream &Nonlocal,ofstream &OrderFile,ofstream &Spectrum
,ofstream&VNE,double&J,ofstream&Nonlocal1)
{
int Xl=Xi[1];
cx_mat Help0;
cx_mat Help1;
cx_mat Help2;
Help0=trans(TIp[0])*TI[0];
complex<double> s;
s=0;
s=Help0(0,0);
Fidelity<<setprecision (15)<<" "<<abs(s)<<endl;
double Myorder;
/*******************************************************/
complex<double> img(0,1);
cx_mat PZ(2,2); PZ.zeros();        PZ(0,0)=1; PZ(1,1)=-1;
cx_mat PI(2,2);  PI.eye(2,2);
cx_mat PX(2,2); PX.zeros();    PX(0,1)=1; PX(1,0)=1;
cx_mat PY(2,2); PY.zeros();    PY(0,1)=-img; PY(1,0)=img;
PX=0.5*PX;
PZ=0.5*PZ;
PY=0.5*PY;
/***************************************************************/
cx_mat * o1=new cx_mat[1];
cx_mat * o2=new cx_mat[1];
cx_mat * o3=new cx_mat[1];
cx_mat * o1X=new cx_mat[1];
cx_mat * o2X=new cx_mat[1];
cx_mat * o3X=new cx_mat[1];

cx_mat *IHelp=new cx_mat[1];

o1[0]=kron(kron(PZ,PI),kron(PI,PI))-kron(kron(PI,PZ),kron(PI,PI))+
kron(kron(PI,PI),kron(PZ,PI))-kron(kron(PI,PI),kron(PI,PZ));
o2[0]=kron(kron(PZ,PI),kron(PI,PI))+kron(kron(PI,PZ),kron(PI,PI))+
kron(kron(PI,PI),kron(PZ,PI))+kron(kron(PI,PI),kron(PI,PZ));
o3[0]=kron(kron(PI,PI),kron(PI,PI));

o1X[0]=kron(kron(PX,PI),kron(PI,PI))+kron(kron(PI,PX),kron(PI,PI))+
kron(kron(PI,PI),kron(PX,PI))+kron(kron(PI,PI),kron(PI,PX));
o3X[0]=kron(kron(PI,PI),kron(PI,PI));
/****************************************************/
/*******************************Magnetization**************************/
int iS=0;
int jS=1;
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
iS=1;
jS=0;
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0]+Help0;
iS=2;
jS=0;
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0]+Help0;
iS=3;
jS=0;
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0]+Help0;

VNE<<"   "<<real(Help0(0,0)/16.00);
cout<< "Magnetization=" <<real(Help0(0,0)/16.00)<<"   ";

double DerivativeX=real(Help0(0,0)/16.00);
/**********************************************************************/
/*******************************MagnetizationDerivative**************************/
 iS=0;
 jS=1;
IHelp[0]=TIp[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TIp[0])*IHelp[0];
iS=1;
jS=0;
IHelp[0]=TIp[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TIp[0])*IHelp[0]+Help0;
iS=2;
jS=0;
IHelp[0]=TIp[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TIp[0])*IHelp[0]+Help0;
iS=3;
jS=0;
IHelp[0]=TIp[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TIp[0])*IHelp[0]+Help0;

DerivativeX=(abs(DerivativeX)-abs(real(Help0(0,0)/16.00)))*(1.00/DeltaFidelity);

VNE<<setprecision (15)<<" "<<abs(DerivativeX);
cout<< "DerivativeMagnetization=" <<abs(DerivativeX)<<"   ";



/**********************CorrelationZZ************************************************/
o1X[0]=kron(kron(PZ,PZ),kron(PI,PI));
for(int i=0;i<4;i++){
iS=i;
jS=(i+1)%(4);
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));
}





o1X[0]=kron(kron(PI,PI),kron(PZ,PZ));
for(int i=0;i<4;i++){
iS=i;
jS=(i+1)%(4);
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));
}


o1[0]=kron(kron(PI,PZ),kron(PI,PI));
o2[0]=kron(kron(PZ,PI),kron(PI,PI));
iS=0;
jS=1;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));

o1[0]=kron(kron(PI,PI),kron(PI,PZ));
o2[0]=kron(kron(PI,PI),kron(PZ,PI));
iS=0;
jS=1;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));


o1[0]=kron(kron(PI,PZ),kron(PI,PI));
o2[0]=kron(kron(PZ,PI),kron(PI,PI));
iS=2;
jS=3;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));

o1[0]=kron(kron(PI,PI),kron(PI,PZ));
o2[0]=kron(kron(PI,PI),kron(PZ,PI));
iS=2;
jS=3;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));



o1[0]=kron(kron(PI,PI),kron(PZ,PI));
o2[0]=kron(kron(PZ,PI),kron(PI,PI));
iS=0;
jS=2;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));

o1[0]=kron(kron(PI,PI),kron(PI,PZ));
o2[0]=kron(kron(PI,PZ),kron(PI,PI));
iS=0;
jS=2;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));


o1[0]=kron(kron(PI,PI),kron(PZ,PI));
o2[0]=kron(kron(PZ,PI),kron(PI,PI));
iS=1;
jS=3;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));

o1[0]=kron(kron(PI,PI),kron(PI,PZ));
o2[0]=kron(kron(PI,PZ),kron(PI,PI));
iS=1;
jS=3;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));



o1X[0]=kron(kron(PZ,PI),kron(PZ,PI));
for(int i=0;i<4;i++){
iS=i;
jS=(i+1)%(4);
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));
}

o1X[0]=kron(kron(PI,PZ),kron(PI,PZ));
for(int i=0;i<4;i++){
iS=i;
jS=(i+1)%(4);
IHelp[0]=TI[0];
RefreshT2site( o1X,o3X, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];
VNE<<"   "<<real(Help0(0,0));
}

/*************************************g*g*******************************/
o1[0]=kron(kron(PX,PX),kron(PX,PX))*(16.00);
o2[0]=o1[0];
iS=0;
jS=1;
IHelp[0]=TI[0];
RefreshT2site( o1,o2, IHelp, iS,jS);
iS=2;
jS=3;
RefreshT2site( o1,o2, IHelp, iS,jS);
Help0=trans(TI[0])*IHelp[0];

VNE<<"   "<<real(Help0(0,0)/1.00);
cout<< "SymmetryX="<<real(Help0(0,0)/1.00)<<endl;
/*********************************************************************/
if(Xl<=70){
cx_mat Spectrum1(Xl*Xl,Xl*Xl);

for(int i=0;i<Xl;i++)
    for(int j=0;j<Xl;j++)
        for(int m=0;m<Xl;m++)
            for(int n=0;n<Xl;n++)
Spectrum1(i*Xl+j,m*Xl+n)=TI[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0);

cx_mat U;
cx_mat V;
vec s0;
svd_econ(U, s0, V, Spectrum1,'b',"std");
vec NVecH=pow(s0,2);
if(Xl>=15){
for(int i=0;i<100;i++)
Spectrum<<J<<"  "<<setprecision(15)<<NVecH(i)<<endl;
Spectrum<<endl<<endl;}
double nsun=0;
for(int i=0;i<Xl*Xl;i++){
if(NVecH(i)>1.0e-15)
nsun=(-NVecH(i)*((log(NVecH(i)))/log(2)))+nsun;
}
VNE<<setprecision(15)<<"  "<<nsun<<endl<<endl;


}



delete []  o1;
delete []  o2;
delete []  o3;
delete []  o1X;
delete []  o2X;
delete []  o3X;

delete [] IHelp;

}
