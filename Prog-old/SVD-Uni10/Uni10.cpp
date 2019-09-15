#include <iostream>
#include <assert.h>
#include <map>
#include <time.h>
#include <uni10.hpp>
#include <armadillo>
#include <mkl.h>

using namespace uni10;
using namespace std;
using namespace arma;

int main(){
    wall_clock timer;
    wall_clock timer2;
    double n_secs;
    timer.tic();
    int chi=50;
    int D=25*2;
    Matrix M(D*chi, D*chi);
    M.randomize();
    // carry out SVD
    vector<Matrix> rets = M.svd();
    Matrix U(rets[0].row(), rets[0].col(), rets[0].isDiag());
    Matrix S(rets[1].row(), rets[1].col(), rets[1].isDiag());
    Matrix VT(rets[2].row(), rets[2].col(), rets[2].isDiag());
    // read in the matrice we just write out
    n_secs = timer.toc();
    cout << "took " << n_secs << " seconds" << endl;
    timer.tic();
    mat A;
    mat U1;
    mat V1;
    vec s1;
    A=randn<mat>(D*chi,D*chi);
    //svd_econ(U1, s1, V1, A,"left","std");
    //svd_econ(U1, s1, V1, A,"left","dc");
    svd_econ(U1, s1, V1, A);
    //svd_econ(U1, s1, V1, A,"left");

    n_secs = timer.toc();
    cout << "took " << n_secs << " seconds" << endl;

return 0;
}
