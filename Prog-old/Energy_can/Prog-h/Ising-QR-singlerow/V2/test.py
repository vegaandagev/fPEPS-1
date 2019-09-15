import pyuni10 as uni10
import copy
import time
import os, sys
import numpy
np = numpy
import MPSclass2

d = 2
D = 3
L = 4

mps2 = MPSclass2.MPS(d,D,L,"rand",4)
mps1 = MPSclass2.MPS(d*2,D+2,L,"rand",4)

#mps_add=mps2+mps1
#print   mps2.norm()

bdi = uni10.Bond(uni10.BD_IN, D)
bdo = uni10.Bond(uni10.BD_OUT, D)
bdi_pys = uni10.Bond(uni10.BD_IN, d*2)
bdo_pys = uni10.Bond(uni10.BD_OUT, d)
A=uni10.UniTensorR([bdi,bdi_pys,bdo,bdo_pys], "A_middle")
numpy.random.seed(2)
A.Randomize(dn_mu=-1,up_var=1,seed=1)
#A.OrthoRand(dn_mu=-1,up_var=1,seed=1)
A=uni10.Dagger(A)
A=uni10.Transpose(A)
#svd=uni10.Svd(A)
 
#print A, type(A.GetBlock()), A.ElemNum(), Norm(A.GetBlock())
M=A.GetBlock()
#print M, M.shape, M.size, type(M.shape),M.shape[0], M.shape[1]
M_trans=M.T
#print M_trans, M_trans.shape, M_trans.size, type(M_trans.shape),M_trans.shape[0], M_trans.shape[1]

#print "Hi", len(M), M.shape, M.size


x=M.shape[0] 
y=M.shape[1] 

#for i in xrange(x):
# for j in xrange(y):
#  print M[i][j]

#A=uni10.UniTensorR
#A=uni10.Matrix









#print "h"
mps_a=MPSclass2.MPS(2,30,12,"rand",4)   #randuni, rand, ortho, Seed
#print "h"

mps_b=MPSclass2.MPS(3,3,14,"rand",5)   #randuni, rand, ortho
#print "h"

#print mps_a.N, mps_a.D, mps_a.d#, mps_a.tensor   
#print mps_a[0], mps_a[1],mps_a[2],mps_a[3]
#print mps_a[0]
#mps_a[0]=2
#print mps_a[0]

#print  mps_a
#norm=mps_a.norm()
#print norm
#mps_a_=mps_a*(-2.0)     #uniform distribution
#norm=mps_a_.norm()
#print norm/4.0
#print mps_a_[0], mps_a[0]
#print mps_a_[1], mps_a[1]

#mps_a=mps_a.normalize()
#mps_b=mps_b.normalize()
#print mps_a.norm(),mps_b.norm()
#print mps_a.product(mps_b)
mps_ab=(mps_a-mps_b)
#print mps_a.product(mps_b)


#print  mps_ab.norm(), mps_a.norm()+mps_b.norm()+(-2)*mps_a.product(mps_b)


#print mps_a[0].PrintDiagram(), mps_a[1].PrintDiagram()


chi=30
mps_c=mps_a.appSVD(chi)
#print mps_a[0].PrintDiagram(), mps_a[1].PrintDiagram()
#mps_a=mps_a*(1.0e-20)
print mps_c.fidel(mps_a)


#mps_a=mps_a*(1/(mps_a.norm()**(0.5)))
#mps_a=mps_a.normalize()
#chi=2
#print mps_a.norm()


#print mps_a.norm(),mps_b.norm(), mps_a.product(mps_b)


#print mps_a_app.norm(),mps_a.norm(), mps_a_app.product(mps_a),mps_a_app.fidel(mps_a)

#print mps_a.distance(mps_a_app)


#mps_a_app=mps_a.appINV(chi)



















