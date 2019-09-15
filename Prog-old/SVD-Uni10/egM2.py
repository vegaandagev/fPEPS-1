import pyUni10 as uni10
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab
import random
import copy
import time


def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni


 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.empty((dim0,dim1))
 t0=time.time()
 v=0
 for i in range(dim0):
  for j in range(dim1):
    Mat_np[i,j]=Mat_uni[i*dim1+j]
 print time.time() - t0, "npUni"

 return  Mat_np



t0=time.time()
chi=2
D=2

M = uni10.Matrix( D*chi,D)
#M.save('M')
#dim0=int(M.row())
#dim1=int(M.col())
#Mat_np=np.empty((dim0,dim1))
#Mat_np=np.load('M')


print M
M.randomize()
rets = M.svd();

U = uni10.Matrix(rets[0].row(), rets[0].col(), rets[0].isDiag());
S = uni10.Matrix(rets[1].row(), rets[1].col(), rets[1].isDiag());
VT = uni10.Matrix(rets[2].row(), rets[2].col(), rets[2].isDiag());
U_=(U*S)*VT
print time.time() - t0, "SVD"

print U,S,VT #S[0], S[1], rets[1], rets[1][0], rets[1][1]


#t0=time.time()

#U=U*VT

#print time.time() - t0, "Multi2"

#M = uni10.Matrix(D*chi*D*chi ,1)
#M.randomize()

M1 = uni10.Matrix(D*chi, 1)
M1.randomize()

t0=time.time()

A_np=Mat_uni_to_np(M)

b_np=Mat_uni_to_np(M1)

print time.time() - t0, "MAt to np"

t0=time.time()

x_np=np.linalg.lstsq(A_np, b_np)[0] 
#x_np=sp.linalg.lstsq(A_np, b_np)[0] 
x=Mat_np_to_Uni(x_np)

print time.time() - t0, "Linear"














