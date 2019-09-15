import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import UniformDisentangler as UD

def   QR_function(a_u,b_u,c_u,d_u):

 N_x=8
 N_y=10
 D=4
 d=2
 chi_boundry=20
 chi_p=4*D

 N_iter=12
 #Opt=["Non","contiguous","cg", 10]       #  [root & None     and     None & contiguous   and   cg&SVD      chi_root]
 #Opt=["Non","Single","cg",20]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD]
 Opt=["Non","root","cg",2]     #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD]
 #Opt=["Non","Ucontiguous","SVD",6]      #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD]



 PEPS_mps=[None]*N_x
 PEPS_mps_left=[None]*N_x
 PEPS_mps_right=[None]*N_x
 PEPS_listten=[None]*N_x

 count_list=[]
 Norm_list=[]
 Norm_list_boundry=[]




 for i in xrange(N_x):
   PEPS_mps[i]=UD.Init_PEPS( N_y, N_x, D, d,i,a_u,b_u,c_u,d_u)


 for i in xrange(N_x):
  PEPS_listten[i]=UD.mps_to_tensor_left(PEPS_mps[i], N_x, N_y, D, d, i)





# MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
# #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
# Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
# #print  Env_left[0].norm()#, PEPS_listmps[0].norm()

# #print  Env_left[N_x-1].norm()**(0.5), Env_right[0].norm()**(0.5)
# N_boundry=Env_left[0].product(Env_right[1])
# #print "Norm:boundryMethod=",  N_boundry
# #print  "Norm:boundryMethod=",  Env_left[1].product(Env_right[2])

# for i in xrange(N_x):
#  PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)

# MPS_R_left=[]
# MPS_R_right=[]
# PEPS_mps_update=0
# ######################################
# for q in xrange(0,N_x-1):
#  if q==0:
#   Dp=d
#   PEPS_mps_update=PEPS_mps[0]
#  else:
#   Dp=D*d

#  MPS_R, mps_Q, Uni_list,Fidel_QR = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)
#  PEPS_mps_update=UD.absorption_left(PEPS_mps_left[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p)
#  MPS_R_left.append(MPS_R)


# N_QR=PEPS_mps_update.norm()
# ##############################################
# print "Norm:QRMethod_left",  N_QR, N_boundry
# ######################################
# return   (N_QR-N_boundry)/N_boundry

 #count_list.append(N_x)
 #Norm_list.append(PEPS_mps_update.norm())
 #Norm_list_boundry.append(Env_left[N_x-2].product(Env_right[N_x-1]))



# for i in xrange(N_x):
#  PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)


# PEPS_mps_update=0
# ################Right-Move######################
# for q in reversed(xrange(1,N_x)):
#  #print "q_right",q
#  if q==N_x-1:
#   Dp=d
#   PEPS_mps_update=PEPS_mps_right[q]
#  else:
#   Dp=D*d

#  MPS_R, mps_Q, Uni_list,Fidel_QR = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)
#  PEPS_mps_update=UD.absorption_right(PEPS_mps_right[q-1],q-1, MPS_R, D, d, N_x, N_y,chi_p)
#  MPS_R_right.append(MPS_R)

#  Norm_val=UD.Norm_based_on_LR(MPS_R,MPS_R_left[q-1], D)
#  Norm_val_boundry=Env_left[q-1].product(Env_right[q])
#  print q, Norm_val, Norm_val_boundry
#  count_list.append(q)
#  Norm_list.append(Norm_val)
#  Norm_list_boundry.append(Norm_val_boundry)
# ##############################################
# print "Norm:QRMethod_right",  PEPS_mps_update.norm()
# ######################################


# file = open("Norm.txt", "w")
# for index in range(len(Norm_list)):
#   file.write(str(count_list[index]) + " " + str(Norm_list[index])+" "+str(Norm_list_boundry[index]) + "\n")
# file.close()






 #########################Ising###################################################

 a_ising=a_u
 a_ising1=c_u



 a_ising.setLabel([0,1,2,3,4])
 a_ising.permute([1,2,0,3,4],3)
 a_ising.setLabel([0,1,2,3,4])

 a_ising1.setLabel([0,1,2,3,4])
 a_ising1.permute([1,2,0,3,4],3)
 a_ising1.setLabel([0,1,2,3,4])


 #print   a_ising.printDiagram()
 B_list=[None]*N_y
 for i in xrange(N_y):
  if i%2==0: 
    B_list[i]=copy.copy(a_ising1) 
  else: 
    B_list[i]=copy.copy(a_ising)

 bdi = uni10.Bond(uni10.BD_IN, 1)
 bdo = uni10.Bond(uni10.BD_IN, D)
 Tem0=uni10.UniTensor([bdi, bdo])
 Tem0.identity()
 Tem0.setLabel([1,-1])
 B_list[0].setLabel([0,-1,2,3,4])
 B_list[0]=B_list[0]*Tem0
 B_list[0].permute([0,1,2,3,4],3)
 #print  B_list[0].printDiagram()

 Tem0.setLabel([4,-4])
 B_list[N_y-1].setLabel([0,1,2,3,-4])
 B_list[N_y-1]=B_list[N_y-1]*Tem0
 B_list[N_y-1].permute([0,1,2,3,4],3)
 #print  B_list[N_y-1].printDiagram()


 PEPS_mps_update=UD.tensor_to_mps_left(B_list, N_x,N_y)
 PEPS_mps_update=PEPS_mps_update.normalize()

 Dp=D*d
 MPS_R, mps_Q, Uni_list, Fidel_QR= UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)


 return   Fidel_QR







 ###############################################################################
 #a_ising1=uni10.UniTensor("test_tensors/heisD4Tens2.dat")
 #a_ising=uni10.UniTensor("test_tensors/heisD4Tens1.dat")

 ##a_ising1=uni10.UniTensor("test_tensors/heisD3Tens2.dat")
 ##a_ising=uni10.UniTensor("test_tensors/heisD3Tens1.dat")

 #a_ising.setLabel([0,1,2,3,4])
 #a_ising.permute([0,2,1,3,4],3)
 #a_ising.setLabel([0,1,2,3,4])

 #a_ising1.setLabel([0,1,2,3,4])
 #a_ising1.permute([0,2,1,3,4],3)
 #a_ising1.setLabel([0,1,2,3,4])


 ##print   a_ising.printDiagram()
 #B_list=[None]*N_y

 #for i in xrange(N_y):
 # if i%2==0: 
 #   B_list[i]=copy.copy(a_ising) 
 # else: 
 #   B_list[i]=copy.copy(a_ising1)


 #bdi = uni10.Bond(uni10.BD_IN, 1)
 #bdo = uni10.Bond(uni10.BD_IN, D)
 #Tem0=uni10.UniTensor([bdi, bdo])
 #Tem0.identity()
 #Tem0.setLabel([1,-1])
 #B_list[0].setLabel([0,-1,2,3,4])
 #B_list[0]=B_list[0]*Tem0
 #B_list[0].permute([0,1,2,3,4],3)
 ##print  B_list[0].printDiagram()

 #Tem0.setLabel([4,-4])
 #B_list[N_y-1].setLabel([0,1,2,3,-4])
 #B_list[N_y-1]=B_list[N_y-1]*Tem0
 #B_list[N_y-1].permute([0,1,2,3,4],3)
 ##print  B_list[N_y-1].printDiagram()


 #PEPS_mps_update=UD.tensor_to_mps(B_list, N_x,N_y)
 ##print  PEPS_mps_update.norm(), PEPS_mps_update[0].printDiagram()
 #PEPS_mps_update=PEPS_mps_update.normalize()

 #Dp=D*d
 #MPS_R, mps_Q, Uni_list = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)



