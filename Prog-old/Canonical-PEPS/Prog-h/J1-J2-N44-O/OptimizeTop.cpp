#include"Header.h"
void OptimizeTop(cx_mat*Itop, O*OperatorH,O*OperatorV,double&Energy,vec*Checker,double&J,int&DoGap )
{
cx_mat * Ou=new cx_mat[2*Int];
cx_mat * Oubar=new cx_mat[2*Int];
cx_mat * Ol=new cx_mat[2*Int];
cx_mat * Olbar=new cx_mat[2*Int];
cx_mat * Or=new cx_mat[2*Int];
cx_mat * Orbar=new cx_mat[2*Int];
cx_mat * Od=new cx_mat[2*Int];
cx_mat * Odbar=new cx_mat[2*Int];
cx_mat * OHelp=new cx_mat[1];
cx_mat * o1=new cx_mat[1];
int iI;
int iS;
for(int q=0; q<2*Int; q++){
Ou[q].zeros(Xi[1],Xi[1]);
Oubar[q].zeros(Xi[1],Xi[1]);
Od[q].zeros(Xi[1],Xi[1]);
Odbar[q].zeros(Xi[1],Xi[1]);
Ol[q].zeros(Xi[1],Xi[1]);
Olbar[q].zeros(Xi[1],Xi[1]);
Or[q].zeros(Xi[1],Xi[1]);
Orbar[q].zeros(Xi[1],Xi[1]);
}
/*********Up**************************************/
for(int q=0; q<Int; q++){
Ou[q]=OperatorH[0].M[q][0];
Oubar[q]=OperatorH[0].M[q][1];
}

for(int q=Int; q<2*Int; q++){
Ou[q]=OperatorH[0].M[q-Int][1];
Oubar[q]=OperatorH[0].M[q-Int][0];
}
/**********Down*****************/
for(int q=0; q<Int; q++){
Od[q]=OperatorH[0].M[q][0];
Odbar[q]=OperatorH[0].M[q][1];
}

for(int q=Int; q<2*Int; q++){
Od[q]=OperatorH[0].M[q-Int][1];
Odbar[q]=OperatorH[0].M[q-Int][0];
}
/**********Right*****************/
for(int q=0; q<Int; q++){
Or[q]=OperatorV[0].M[q][0];
Orbar[q]=OperatorV[0].M[q][1];
}
for(int q=Int; q<2*Int; q++){
Or[q]=OperatorV[0].M[q-Int][1];
Orbar[q]=OperatorV[0].M[q-Int][0];
}
/*************Left*******************/
for(int q=0; q<Int; q++){
Ol[q]=OperatorV[0].M[q][0];
Olbar[q]=OperatorV[0].M[q][1];
}
for(int q=Int; q<2*Int; q++){
Ol[q]=OperatorV[0].M[q-Int][1];
Olbar[q]=OperatorV[0].M[q-Int][0];
}


/*************End---Or---Orbar***********************/
int num=0;
double E1;
int D=Xi[1]*Xi[1]*Xi[1]*Xi[1];
int Xl=Xi[1];
int m=10;
int W=2;
/////////////////////////////////////
cx_vec *Vechelp=new cx_vec[1];
Vechelp[0]=zeros<cx_vec>(D);
cx_vec * q;
cx_vec  a(m+1);
cx_vec  b(m+1);
cx_vec Vec;
cx_vec Vec1;
vec eigval1;
cx_mat eigvec1;
cx_mat V(m,m);
cx_vec* r=new cx_vec[1];
r[0]=randn<cx_vec>(D);
cx_mat Q;
for(int i=0;i<Xl;i++)
    for(int j=0;j<Xl;j++)
        for(int m=0;m<Xl;m++)
            for(int n=0;n<Xl;n++)
r[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n)=Itop[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0);
Vec=r[0]+0.000001*randn<cx_vec>(D);
for(int p=0;p<=W;p++)
{

if(p==0)
q=new cx_vec[m+1];
for(int j=0;j<m+1;j++)
q[j].zeros(D);
a.zeros(m+1);
b.zeros(m+1);
r[0]=Vec;
b[0]=norm(r[0],2);
for(int i=0;i<m;i++){

if(norm(r[0],2)>1.0e-8)
q[i+1]=r[0]/(b(i));
else
break;

Vechelp[0]=q[i+1];
MultiTop(Ou,Oubar,Or,Orbar,Ol,Olbar,Od,Odbar,Vechelp,Xl,OperatorH);

r[0]=(Vechelp[0])-(b(i)*q[i]);


a(i+1)=trace(trans(q[i+1])*r[0]);
r[0]=r[0]-(a(i+1)*q[i+1]);


//orthogonalization
{
for(int n=1;n<=i;n++)//r is determined as to be orthonormal to all other vectors<=i
r[0]=r[0]-(trace(trans(q[n])*r[0])*q[n]);
}
if(i<=m-2)
b(i+1)=norm(r[0],2);
}

V.zeros(m,m);

for(int i=1;i<m-1;i++)
{
V(i,i-1)=b(i);
 V(i,i)=a(i+1);
  V(i,i+1)=b(i+1);
}
V(0,0)=a(1);
V(0,1)=b(1);
V(m-1,m-1)=a(m);
V(m-1,m-2)=b(m-1);

//cout<<V<<endl;
eig_sym(eigval1, eigvec1, V, "standard");


Q.zeros(D,m);
  for(int n=0;n<m;n++)
     Q.col(n)=q[n+1];
       Vec=Q*eigvec1.col(0);




if(p==W && num==0){
p=-1;
m=m+10;
E1=eigval1(0);
num++;
delete [] q;
q=NULL;
}else if(p==W)
{
    num++;
    if(abs(eigval1(0))>0.000001){
if(   (((abs(eigval1(0)-E1))/(abs(eigval1(0))))< 0.000000001)  )
num++;
else if(m<=150){
p=-1;
m=m+10;
E1=eigval1(0);
delete [] q;
q=NULL;
}}else{
if(  ((abs(eigval1(0)-E1))< 0.0000000001)  )
num++;
else if(m<=150){
p=-1;
m=m+10;
E1=eigval1(0);
delete [] q;
q=NULL;
}
}
}

}


cout<< setprecision(15)<<"J="<<J<<"  Energy="<<
((eigval1(0))+Energy)/16<<" m="<<m<<endl;

Checker[0](0)=Checker[0](1);
Checker[0](1)=((eigval1(0))+Energy)/16;
Vec=Vec/norm(Vec,2);
for(int i=0;i<Xl;i++)
    for(int j=0;j<Xl;j++)
        for(int m=0;m<Xl;m++)
            for(int n=0;n<Xl;n++)
Itop[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0)=Vec(i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n);

/**********************************/
if(DoGap==1){
delete [] q;
m=30;
num=0;

if(norm(Vec1,2)>0.0001)
r[0]=Vec1/norm(Vec1,2);
     else{
Vec1=randu<cx_vec>(D);
Vec1=Vec1/norm(Vec1,2);
r[0]=Vec1/norm(Vec1,2);}

r[0]=r[0]-(trace(trans(Vec)*r[0])*Vec);

for(int p=0;p<=W;p++)                            //number of restarted iterations of Lanczos
{
if(p==0)
q=new cx_vec[m+1];
for(int j=0;j<m+1;j++)
q[j].zeros(D);
a.zeros(m+1);
b.zeros(m+1);


	r[0]=Vec1;
	r[0]=r[0]-(trace(trans(Vec)*r[0])*Vec);

      b[0]=norm(r[0],2);

      for(int i=0;i<m;i++){

	if(norm(r[0],2)>1.0e-4)
	  q[i+1]=r[0]/(b(i));
	else
	  break;

q[i+1]=q[i+1]-(trace(trans(Vec)*q[i+1])*Vec);

Vechelp[0]=q[i+1];
MultiTop(Ou,Oubar,Or,Orbar,Ol,Olbar,Od,Odbar,Vechelp,Xl,OperatorH);

r[0]=(Vechelp[0])-(b(i)*q[i]);

	a(i+1)=trace(trans(q[i+1])*r[0]);
	r[0]=r[0]-(a(i+1)*q[i+1]);
	//orthogonalization
	{
	  for(int n=1;n<=i;n++)//r is determined as to be orthonormal to all other vectors<=i
	    r[0]=r[0]-(trace(trans(q[n])*r[0])*q[n]);
	}
	if(i<=m-2)
	  b(i+1)=norm(r[0],2);
      }
V.zeros(m,m);
      for(int i=1;i<m-1;i++)
	{
	  V(i,i-1)=b(i);
	  V(i,i)=a(i+1);
	  V(i,i+1)=b(i+1);
	}
      V(0,0)=a(1);
      V(0,1)=b(1);
      V(m-1,m-1)=a(m);
      V(m-1,m-2)=b(m-1);

eig_sym(eigval1, eigvec1, V, "standard");
Q.zeros(D,m);
for(int n=0;n<m;n++)
	Q.col(n)=q[n+1];
      Vec1=Q*eigvec1.col(0);


if(p==W && num==0){
p=-1;
m=m+10;
E1=eigval1(0);
num++;
delete [] q;
q=NULL;
}else if(p==W)
{
    num++;
    if(abs(eigval1(0))>0.000001){
if(   (((abs(eigval1(0)-E1))/(abs(eigval1(0))))< 0.000000001)  )
num++;
else if(m<=250){
p=-1;
m=m+10;
E1=eigval1(0);
delete [] q;
q=NULL;
}}else{
if(  ((abs(eigval1(0)-E1))< 0.0000000001)  )
num++;
else if(m<=250){
p=-1;
m=m+10;
E1=eigval1(0);
delete [] q;
q=NULL;
}
}
}
}


Checker[0](2)=((eigval1(0)+Energy)/(16))-Checker[0](1);

}

delete [] Ou;delete [] Oubar;delete [] Or;delete [] Orbar;
delete [] q;
delete [] Ol;delete [] Olbar;delete [] Od;delete [] Odbar;
delete [] r;delete [] OHelp;delete [] o1;
delete [] Vechelp;
}



void MultiTop(cx_mat*Ou,cx_mat*Oubar,cx_mat*Or,cx_mat*Orbar,cx_mat*Ol,cx_mat*Olbar,cx_mat*Od,cx_mat*Odbar,
cx_vec*Vec,int&Xl,O*OperatorH){
int iS;
int jS;
cx_mat * o1=new cx_mat[1];
cx_mat * o2=new cx_mat[1];
cx_mat *I=new cx_mat[1];
cx_mat *I1=new cx_mat[1];
cx_mat *Iresult=new cx_mat[1];
I[0].zeros(Xl*Xl*Xl*Xl,2);
Iresult[0].zeros(Xl*Xl*Xl*Xl,2);
for(int i=0;i<Xl;i++)
    for(int j=0;j<Xl;j++)
        for(int m=0;m<Xl;m++)
            for(int n=0;n<Xl;n++)
I[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0)=Vec[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n);
I1[0]=I[0];
/*******************/

for(int q=0; q<=Int; q++){
o1[0]=Ou[q];
o2[0]=Oubar[q];
iS=0;
jS=1;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
}
/********************************/


for(int q=0; q<Int; q++){
o1[0]=Ol[q];
o2[0]=Olbar[q];
iS=0;
jS=2;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
}
/*********/
for(int q=0; q<=Int; q++){
o1[0]=Od[q];
o2[0]=Odbar[q];
iS=2;
jS=3;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
}
/******************/
for(int q=0; q<Int; q++){
o1[0]=Or[q];
o2[0]=Orbar[q];
iS=1;
jS=3;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
}
/**********************************/
o1[0]=OperatorH[0].M[4][0];
o2[0]=OperatorH[0].M[4][1];
iS=2;
jS=1;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
/********************************/

o1[0]=OperatorH[0].M[3][0];
o2[0]=OperatorH[0].M[3][1];
iS=0;
jS=3;
I1[0]=I[0];
RefreshT2site( o1,o2, I1, iS,jS);
Iresult[0]=I1[0]+Iresult[0];
/********************************/



for(int i=0;i<Xl;i++)
    for(int j=0;j<Xl;j++)
        for(int m=0;m<Xl;m++)
            for(int n=0;n<Xl;n++)
Vec[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n)=Iresult[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0);
delete [] o1;
delete [] o2;
delete [] I;
delete [] I1;
delete [] Iresult;
}
