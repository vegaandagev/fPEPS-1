#include"Header.h"
using namespace std;
int main()
{
    int DoGap;
    int DoFidelity;
O * OperatorH=new O[1];
O * OperatorV=new O[1];
cx_mat * Itop=new cx_mat[1];
cx_mat * TI=new cx_mat[1];
cx_mat * TIp=new cx_mat[1];
vec * Checker=new vec[1];
Checker[0].zeros(3);
Checker[0](0)=0;
Checker[0](1)=0;
Checker[0](2)=0;
cx_mat Help;
int Xl=Xi[1];
int iS;
int jS;
ofstream OrderFile;
OrderFile.open("Files/OrderFile.text");
ofstream Fidelity;
Fidelity.open("Files/Fidelity.text");
ofstream GapEnergy;
GapEnergy.open("Files/Gap.text");
ofstream Nonlocal;
Nonlocal.open("Files/Nonlocal.text");
ofstream Nonlocal1;
Nonlocal1.open("Files/Nonlocal1.text");
ofstream Spectrum;
Spectrum.open("Files/Spectrum.text");
ofstream VNE;
VNE.open("Files/VNE.text");
ofstream Energyf;
Energyf.open("Files/Energy.text");
int x,w,p;
double Energy;
double J,delta,landa, grid;
IInitial(Itop);

J=CouplingGamma;
grid=Grid;
delta=DeltaFidelity;
if(DeltaFidelity<0.0000001)
    DoFidelity=1;
    else
  DoFidelity=2;

for(x=0;x<Points;x++){
if(x==0)
;
else
J=J+(CSign*grid);

OrderFile<<J;
Nonlocal1<<J;
Fidelity<<J;
Nonlocal<<endl;
Spectrum<<endl;
VNE<<J;
Energyf<<J;
GapEnergy<<J;



for(p=0;p<DoFidelity;p++){
landa=J-delta*p;
OInitial(OperatorH,OperatorV,landa,Energy);
w=1;
for(int q=0;q<=w;q++){
cout<<q<<endl;
DoGap=0;
OptimizeTop(Itop, OperatorH,OperatorV,Energy,Checker,landa,DoGap);
if(abs(Checker[0](0))>0.000001){
if( (((abs(Checker[0](0)-Checker[0](1)))/(abs(Checker[0](0))))< Accuracy) && (q==w) )
;
else if(q==w)
w=w+1;
}else{
if( (abs(Checker[0](0)-Checker[0](1))< 0.00000000001) && (q==w) )
;
else if(q==w)
w=w+1;
}
if((q==w)&&(p==0)){
//DoGap=1;
//OptimizeTop( Itop, OperatorH,OperatorV,Energy,Checker,landa,DoGap);
Storing(Itop,TI);
Energyf<<setprecision(10)<<"  "<<Checker[0](1)<<endl;
GapEnergy<<setprecision(10)<<"  "<<Checker[0](2)<<endl;
}



if((q==w)&&(p==int(DoFidelity-1)))
Storing(Itop,TIp);

}
}

GSF(TI,TIp,Fidelity,Nonlocal,OrderFile,Spectrum,VNE,J,Nonlocal1);
}
delete [] OperatorV;
delete [] Itop;
;delete [] Checker;delete [] TI;
delete [] TIp;delete [] OperatorH;
return 0;
}





