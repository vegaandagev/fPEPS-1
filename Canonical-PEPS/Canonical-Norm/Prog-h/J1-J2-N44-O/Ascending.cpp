#include"Header.h"
void RefreshT2site( cx_mat*oi, cx_mat*oj, cx_mat *I, int&iS, int&jS )
{
int Xl=Xi[1];
int num;
num=0;
cx_mat Help0;
cx_mat Help1;
cx_mat Help3;
complex<double> s;
Help0.zeros(Xl,Xl*Xl*Xl);
Help3.zeros(Xl,Xl*Xl*Xl);
for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help3(i,j*Xl*Xl+m*Xl+n)=I[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0);

if(iS==0 || jS==0){

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help0(i,j*Xl*Xl+m*Xl+n)=Help3(i,j*Xl*Xl+m*Xl+n);

if(iS==0)
Help1=oi[0]*Help0;
else if(jS==0)
Help1=oj[0]*Help0;


for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help3(i,j*Xl*Xl+m*Xl+n)=Help1(i,j*Xl*Xl+m*Xl+n);
num++;
}
if(iS==1 || jS==1){

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help0(j,i*Xl*Xl+m*Xl+n)=Help3(i,j*Xl*Xl+m*Xl+n);


if(iS==1)
Help1=oi[0]*Help0;
else if(jS==1)
Help1=oj[0]*Help0;


for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help3(i,j*Xl*Xl+m*Xl+n)=Help1(j,i*Xl*Xl+m*Xl+n);

if(num==1){
    for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
I[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0)=Help3(i,j*Xl*Xl+m*Xl+n);
}
num++;
}


if(iS==2|| jS==2){

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help0(m,i*Xl*Xl+j*Xl+n)=Help3(i,j*Xl*Xl+m*Xl+n);


if(iS==2)
Help1=oi[0]*Help0;
else if(jS==2)
Help1=oj[0]*Help0;

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help3(i,j*Xl*Xl+m*Xl+n)=Help1(m,i*Xl*Xl+j*Xl+n);

if(num==1){
    for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
I[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0)=Help3(i,j*Xl*Xl+m*Xl+n);
}

num++;
}


if(iS==3|| jS==3){

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help0(n,i*Xl*Xl+j*Xl+m)=Help3(i,j*Xl*Xl+m*Xl+n);

if(iS==3)
Help1=oi[0]*Help0;
else if(jS==3)
Help1=oj[0]*Help0;


for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
Help3(i,j*Xl*Xl+m*Xl+n)=Help1(n,i*Xl*Xl+j*Xl+m);

for(int i=0; i<Xl;i++)
 for(int j=0; j<Xl;j++)
   for(int m=0; m<Xl;m++)
     for(int n=0; n<Xl;n++)
I[0](i*Xl*Xl*Xl+j*Xl*Xl+m*Xl+n,0)=Help3(i,j*Xl*Xl+m*Xl+n);
}

}
/***************************************************************************************************************************/
