#include"Header.h"
void IInitial(cx_mat *Itop) //I= isometry
{
cx_mat d1;
cx_mat C1;
cx_mat V1;
vec s1;
C1=randn<cx_mat>(Xi[Lay-1]*Xi[Lay-1]*Xi[Lay-1]*Xi[Lay-1],1);
svd_econ(d1, s1, V1, C1,'l',"std");
Itop[0]=d1;
}


