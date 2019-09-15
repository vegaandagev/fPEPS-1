#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED
#include <iostream>
#include <assert.h>
#include <map>
#include <time.h>
#include <iomanip>
#include <limits>
#include <armadillo>
#include "itensor/all.h"
using namespace itensor;
using namespace std;
using namespace arma;


auto   Init_Corner(int &, int &)     -> vector<ITensor>;
auto   Init_Tensors(int &, int &)     -> vector<ITensor>;
auto   Init_Tensors_one(vector<ITensor> & ,int &, int &)  -> vector<ITensor>;
auto   Init_Tensors_fill(vector<ITensor> & ,int &, int &, double &)  -> vector<ITensor>;

auto   Label_TensCTM(vector<ITensor> &, vector<ITensor> & )     -> vector<ITensor>;

auto   norm_CTM( vector<ITensor> &, vector<ITensor> & ) -> double;

auto Left_move( vector<ITensor> &,vector<ITensor> &, int & )  ->  vector<ITensor>;
auto Right_move( vector<ITensor> &,vector<ITensor> &, int & )  ->  vector<ITensor>;
auto Permute_Env(vector<ITensor> & ) -> vector<ITensor>;
auto Permute_Ten(vector<ITensor> & Ten_Iten) -> vector<ITensor>;
auto Permute1_Env(vector<ITensor> & ) -> vector<ITensor>;
auto Permute1_Ten(vector<ITensor> & Ten_Iten) -> vector<ITensor>;


#endif
