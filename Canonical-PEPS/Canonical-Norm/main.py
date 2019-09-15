import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import UniformDisentangler as UD

N_x=6
N_y=6
D=2
d=2
chi_boundry=20
chi_single=20
chi_p=6
chi_express=6
Model="ITF"                       #ITF , Heis#
h_coupling=[0.0, 1.0]
N_iter=20
#Opt=["Non","contiguous","cg", 10]           #  [root & None     and     None & contiguous   and   cg&SVD      chi_root] 
Opt=["Non","Single","cg",10,5]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
#Opt=["Non","ACCcontiguous","cg",6,5]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
#Opt=["Non","Ucontiguous","SVD",6]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 

PEPS_mps=[None]*N_x
PEPS_mps_left=[None]*N_x
PEPS_mps_right=[None]*N_x
PEPS_listten=[None]*N_x

PEPS_mps_Qr=[None]*N_x
PEPS_mps_Ql=[None]*N_x
PEPS_Q=[None]*N_x

mps_boundry_left=[None]*N_x
mps_boundry_right=[None]*N_x
mps_boundry_tem=[None]*N_x

count_list=[]
Norm_list=[]
Norm_list_boundry=[]

for i in xrange(N_x):
  #PEPS_mps[i]=UD.Init_PEPS( N_y, N_x, D, d,i)
  PEPS_mps[i]=UD.Init_PEPS( N_y, N_x, D, d,i)



 
for i in xrange(N_x):
  PEPS_listten[i]=UD.mps_to_tensor_left(PEPS_mps[i], N_x, N_y, D, d, i)

#UD.Store_f(PEPS_listten, N_x, N_y)
UD.Reload_f(PEPS_listten, N_x, N_y)




for i in xrange(N_x):
 PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)
for i in xrange(N_x):
 PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)

for i in xrange(N_x):
 PEPS_mps_Ql[i]=copy.copy(PEPS_mps_left[i])
for i in xrange(N_x):
 PEPS_mps_Qr[i]=copy.copy(PEPS_mps_right[i])
for i in xrange(N_x):
 PEPS_Q[i]=copy.copy(PEPS_listten[i])





MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
 #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
#print  Env_left[0].norm(), Env_left[1].norm(), Env_left[2].norm(), Env_left[3].norm()
#print  Env_right[3].norm(), Env_right[2].norm(), Env_right[1].norm(), Env_right[0].norm()
Norm_val_boundry=Env_left[0].product(Env_right[1])
print Env_left[0].product(Env_right[1])
#for i in xrange(N_x):
# print "NormLeft =", i,Env_left[i].norm(), mps_boundry_left[i].norm()
# print "Normright=", i,Env_right[i].norm(), mps_boundry_right[i].norm()

#for i in xrange(N_x-1):
 #print "Norm:boundryMethod=", i,Env_left[i].product(Env_right[i+1]), mps_boundry_left[i].product(mps_boundry_right[i+1])




MPO_Ten=UD.make_boundry_MPO(PEPS_Q, N_y, N_x)
 #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
#print  Env_left[0].norm(), Env_left[1].norm(), Env_left[2].norm(), Env_left[3].norm()
#print  Env_right[3].norm(), Env_right[2].norm(), Env_right[1].norm(), Env_right[0].norm()
Norm_val_boundry=Env_left[0].product(Env_right[1])
print Env_left[0].product(Env_right[1]),Env_left[1].product(Env_right[2])
#for i in xrange(N_x):
# print "NormLeft =", i,Env_left[i].norm(), mps_boundry_left[i].norm()
# print "Normright=", i,Env_right[i].norm(), mps_boundry_right[i].norm()

#for i in xrange(N_x-1):
 #print "Norm:boundryMethod=", i,Env_left[i].product(Env_right[i+1]), mps_boundry_left[i].product(mps_boundry_right[i+1])






Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
Q_norm_list=[]
N_iter_list=[]
Trun_list=[]

Fidel_val=1
E_coulmn=[None]*N_x
E_row=[None]*(N_x-1)
E_mag_coulmn=[None]*N_x
E_mag_row=[None]*(N_x-1)


N_iter_QR=1
N_iter_QR_first=1

j_coupling=3.05
z_coupling=1.0
Grid__coupling=0.00
N_coupling_probe=1

for x_coupling in xrange(N_coupling_probe):

 if x_coupling==0:
  j_coupling=j_coupling
 else:
  j_coupling=j_coupling+Grid__coupling

 h_coupling=[z_coupling, j_coupling]

 #if x_coupling==0:
  #h_coupling=[0, 1]

 if x_coupling==0:
  List_delN=UD.Short_TrotterSteps_start(N_iter_QR_first)
 else:
  List_delN=UD.Short_TrotterSteps(N_iter_QR)

 print h_coupling, List_delN

 E_min=1
 E_0=1
 E_1=1000
 count_iter=0
 for delta, N_iter in List_delN:

  H=UD.Heisenberg(h_coupling,d,Model)
  U_ham = uni10.UniTensor(H.bond(), "U_ham");
  U_ham.putBlock(uni10.takeExp(-delta, H.getBlock()))

  H0=UD.Heisenberg0(h_coupling,d,Model)
  U_ham0 = uni10.UniTensor(H0.bond(), "U_ham0");
  U_ham0.putBlock(uni10.takeExp(-delta, H0.getBlock()))

  HN=UD.HeisenbergN(h_coupling,d,Model)
  U_hamN = uni10.UniTensor(HN.bond(), "U_hamN");
  U_hamN.putBlock(uni10.takeExp(-delta, HN.getBlock()))

  m_mag=UD.Mag(d)
  m_mag0=UD.Mag0(d)
  m_magN=UD.MagN(d)

#  for i in xrange(N_x):
#   PEPS_mps_left[i]=copy.copy(PEPS_mps_leftU[i])
#  for i in xrange(N_x):
#   PEPS_mps_right[i]=copy.copy(PEPS_mps_rightU[i])
#  for i in xrange(N_x):
#   PEPS_listten[i]=copy.copy(PEPS_listtenU[i])


  for i_iter in xrange(N_iter):
   print "step=", i_iter, "delta=", delta
   E_0=E_1*1.0
##############QR-Left#####################################################
   for q in xrange(0,N_x-1):
    if q==0:
     Dp=d
    else:
     Dp=D*d
    MPS_R, mps_Q, Uni_list , Fidel_val= UD.QR_canon(PEPS_mps_Ql[q], N_iter, Dp, D,Opt)
    PEPS_mps_Ql[q]=mps_Q*1.0
    PEPS_Q[q]=UD.mps_to_tensor_left(PEPS_mps_Ql[q], N_x, N_y, D, d, q)
    PEPS_mps_Qr[q]=UD.tensor_to_mps_right(PEPS_Q[q], N_x, N_y)

    PEPS_mps_Ql[q+1]=UD.absorption_left(PEPS_mps_Ql[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p, Trun_list)
    PEPS_Q[q+1]=UD.mps_to_tensor_left(PEPS_mps_Ql[q+1], N_x, N_y, D, d, q+1)
    PEPS_mps_Qr[q+1]=UD.tensor_to_mps_right(PEPS_Q[q+1], N_x, N_y)
    print q,  Norm_val_boundry,  MPS_R.norm(),  (Norm_val_boundry-MPS_R.norm())/(Norm_val_boundry*MPS_R.norm())**(0.5), Fidel_val
    N_iter_list.append((Norm_val_boundry-MPS_R.norm())/Norm_val_boundry)

   print PEPS_mps_Ql[N_x-1].norm()

file = open("norm.txt", "w")
for index in range(len(N_iter_list)):
  file.write(str(index) + " " + str(N_iter_list[index])+" "+ "\n")
file.close()



