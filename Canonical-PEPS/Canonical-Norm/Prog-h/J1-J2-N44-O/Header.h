#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED
#include <iostream>
#include<complex>
#include<math.h>
#include <armadillo>
#include<fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace arma;
const int Lay=2; //Number of layers of tensor Network
const int Int=3; //numbers of interactions+1:for example:See Example in Code to get the line!
const int Xi[Lay]={16,16};  //last term is important
const double CouplingGamma=2.000;
const double CouplingJ2=0.000;  //0.5
const double CouplingJ1=1.000;   //1
const int Points=1;
const int CSign=1;
const double Grid=0.01;
const double DeltaFidelity=0.00;   //if it's equall to zero, fidelity is not calculated.

const double Accuracy=1.0e-8;
const double AccuracyLanczos=1.0e-8;
const int KrylovDimension=30;
const int StartKrylovSpace=15;
const int LanczosRepeat=4;



typedef struct{
cx_mat M[Int+2][2];
}O;
void OInitial(O*,O*,double&,double &);
void IInitial(cx_mat *);
void RefreshT2site(cx_mat*,cx_mat*,cx_mat *,int&,int& );
void OptimizeTop(cx_mat*, O*,O*,double&,vec*,double&,int& );
void Storing(cx_mat*,cx_mat*);
void GSF( cx_mat * , cx_mat * , ofstream &, ofstream &, ofstream &, ofstream &, ofstream &,double&,ofstream &);
void MultiTop(cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_mat*,cx_vec*,int&,O*OperatorH);
#endif // HEADER_H_INCLUDED
