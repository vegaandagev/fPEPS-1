import pyuni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass2
import UniformDisentangler as UD
import time
import os, sys


N_x=4
N_y=4
D=3
d=2
chi_boundry=10
chi_p=3*D

N_iter=[12,10]
#Opt=["Non","contiguous","cg"]           #  [root & None     and     None & contiguous   and   cg&SVD] 
Opt=["Non","None","SVD", "N",8]           #  [root & None     and     None & contiguous   and   cg&SVD] 
#Opt=["Non","moreAcc","cg", "N",2]           #  [root & None     and     None & contiguous   and   cg&SVD] 

#Opt=["Non","None"]           #  [root & None     and     None & contiguous] 

PEPS_mps=[None]*N_x
PEPS_mps_left=[None]*N_x
PEPS_mps_right=[None]*N_x
PEPS_listten=[None]*N_x

for i in xrange(N_x):
 PEPS_mps[i]=UD.Init_PEPS( N_y, N_x, D, d, i, Opt[3], Opt[4])


for i in xrange(N_x):
 PEPS_listten[i]=UD.mps_to_tensor_left(PEPS_mps[i], N_x, N_y, D, d, i)




MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
#print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
#print  Env_left[0].norm()#, PEPS_listmps[0].norm()

#print  Env_left[N_x-1].norm()**(0.5), Env_right[0].norm()**(0.5)
print "Norm:boundryMethod=",  Env_left[0].product(Env_right[1])
print  "Norm:boundryMethod=",  Env_left[1].product(Env_right[2])


for i in xrange(N_x):
 PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)

PEPS_mps_update=0
################Left-Move######################
for q in xrange(0,N_x-1):
 print "q_left", q
 if q==0:
  Dp=d
  PEPS_mps_update=PEPS_mps_left[0]
 else:
  Dp=D*d

 MPS_R, mps_Q, Uni_list = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)

 if q==0:Q_peps=UD.Qmps_to_Qtensor_left(mps_Q, 1, D, d)
 else: Q_peps=UD.Qmps_to_Qtensor_left(mps_Q, D, D, d)

 PEPS_mps_update=UD.absorption_left(PEPS_mps_left[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p)

##############################################
print "Norm:QRMethod_left",  PEPS_mps_update.norm()
######################################

for i in xrange(N_x):
 PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)


PEPS_mps_update=0
################Right-Move######################
for q in reversed(xrange(1,N_x)):
 print "q_right",q
 if q==N_x-1:
  Dp=d
  PEPS_mps_update=PEPS_mps_right[q]
 else:
  Dp=D*d

 MPS_R, mps_Q, Uni_list = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)

 if q==N_x-1:Q_peps=UD.Qmps_to_Qtensor_right(mps_Q, 1, D, d)
 else: Q_peps=UD.Qmps_to_Qtensor_right(mps_Q, D, D, d)

 PEPS_mps_update=UD.absorption_right(PEPS_mps_right[q-1],q-1, MPS_R, D, d, N_x, N_y,chi_p)

##############################################
print "Norm:QRMethod_right",  PEPS_mps_update.norm()
######################################








