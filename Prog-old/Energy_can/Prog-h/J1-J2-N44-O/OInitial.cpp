#include"Header.h"
void OInitial(O*OperatorH,O*OperatorV,double&J,double &Energy)
{
complex<double> img(0,1);
cx_mat Identity;
Identity.eye(Xi[0],Xi[0]);
vec eigval;
cx_mat h;
cx_mat PZ(2,2); PZ.zeros();        PZ(0,0)=1; PZ(1,1)=-1;
cx_mat PI(2,2);  PI.eye(2,2);
cx_mat I(3,3);  I.eye(3,3);
cx_mat PX(2,2); PX.zeros();    PX(0,1)=1; PX(1,0)=1;
cx_mat PY(2,2); PY.zeros();    PY(0,1)=-img; PY(1,0)=img;
/***********/
PX=PX*(1.00/2.00);
PZ=PZ*(1.00/2.00);
PY=PY*(1.00/2.00);


OperatorH[0].M[0][0]=-J*kron(kron(PX,PI),kron(PI,PI))-J*kron(kron(PI,PX),kron(PI,PI))-J*kron(kron(PI,PI),kron(PX,PI))
-J*kron(kron(PI,PI),kron(PI,PX))+CouplingJ1*kron(kron(PZ,PZ),kron(PI,PI))+CouplingJ1*kron(kron(PZ,PI),kron(PZ,PI))+
CouplingJ1*kron(kron(PI,PI),kron(PZ,PZ))+CouplingJ1*kron(kron(PI,PZ),kron(PI,PZ))-((3*16)/4)*kron(kron(PI,PI),kron(PI,PI))+
CouplingJ2*kron(kron(PZ,PI),kron(PI,PZ))+CouplingJ2*kron(kron(PI,PZ),kron(PZ,PI));
OperatorH[0].M[0][1]=kron(kron(PI,PI),kron(PI,PI));


OperatorH[0].M[1][0]=CouplingJ1*kron(kron(PI,PZ),kron(PI,PI))+CouplingJ2*kron(kron(PI,PI),kron(PI,PZ));
OperatorH[0].M[1][1]=kron(kron(PZ,PI),kron(PI,PI));

OperatorH[0].M[2][0]=CouplingJ2*kron(kron(PI,PZ),kron(PI,PI))+CouplingJ1*kron(kron(PI,PI),kron(PI,PZ));
OperatorH[0].M[2][1]=kron(kron(PI,PI),kron(PZ,PI));

OperatorH[0].M[3][0]=CouplingJ2*kron(kron(PI,PI),kron(PI,PZ));
OperatorH[0].M[3][1]=kron(kron(PZ,PI),kron(PI,PI));

OperatorH[0].M[4][0]=CouplingJ2*kron(kron(PI,PZ),kron(PI,PI));
OperatorH[0].M[4][1]=kron(kron(PI,PI),kron(PZ,PI));


OperatorV[0].M[0][0]=CouplingJ1*kron(kron(PI,PI),kron(PZ,PI))+CouplingJ2*kron(kron(PI,PI),kron(PI,PZ));
OperatorV[0].M[0][1]=kron(kron(PZ,PI),kron(PI,PI));

OperatorV[0].M[1][0]=CouplingJ2*kron(kron(PI,PI),kron(PZ,PI))+CouplingJ1*kron(kron(PI,PI),kron(PI,PZ));
OperatorV[0].M[1][1]=kron(kron(PI,PZ),kron(PI,PI));


OperatorV[0].M[2][0]=0*kron(kron(PI,PI),kron(PI,PI));
OperatorV[0].M[2][1]=0*kron(kron(PI,PI),kron(PI,PI));



//h=-kron(PX,PX)-J*kron(PZ,PI);
//eig_sym(eigval, h);
//OperatorH[0].M[0][0]=(-eigval( (Xi[0]*Xi[0])-1 ))*Identity;
//OperatorH[0].M[0][1].eye(Xi[0],Xi[0]);
//OperatorV[0].M[0][0]=(-eigval( (Xi[0]*Xi[0])-1 ))*Identity;
//OperatorV[0].M[0][1].eye(Xi[0],Xi[0]);
//Energy=0*eigval( (Xi[0]*Xi[0])-1 );
Energy=(3*16);
}



