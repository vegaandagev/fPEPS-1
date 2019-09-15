import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import UniformDisentangler as UD

N_x=8
N_y=8
D=3
d=2
chi_boundry=20
chi_single=20
chi_p=3*D
chi_express=3*D
Model="ITF"                       #ITF , Heis#
h_coupling=[0.0, 1.0]
N_iter=15
#Opt=["Non","contiguous","cg", 10]           #  [root & None     and     None & contiguous   and   cg&SVD      chi_root] 
Opt=["Non","Single","cg",10,5]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
#Opt=["Non","ACCcontiguous","cg",6,5]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
#Opt=["Non","Ucontiguous","SVD",6]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 

PEPS_mps=[None]*N_x
PEPS_mps_left=[None]*N_x
PEPS_mps_right=[None]*N_x
PEPS_listten=[None]*N_x

PEPS_mps_leftU=[None]*N_x
PEPS_mps_rightU=[None]*N_x
PEPS_listtenU=[None]*N_x

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

#for Location in xrange(N_x):
# #print Location
# mps_boundry_left[Location]=UD.make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single,N_x)
# #print mps_boundry_left[Location].norm()

#for Location in xrange(N_x):
# #print N_x-Location-1
# PEPS_ten=UD.rotate(PEPS_listten[N_x-Location-1])
# mps_boundry_tem[Location]=UD.make_Env_singleLayer(PEPS_ten, Location, mps_boundry_tem[Location-1],d,chi_single,N_x)
# #print mps_boundry_tem[Location].norm()
# mps_boundry_right[N_x-Location-1]=copy.copy(mps_boundry_tem[Location])

#MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
 #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
#Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
#print  Env_left[0].norm(), Env_left[1].norm(), Env_left[2].norm(), Env_left[3].norm()
#print  Env_right[3].norm(), Env_right[2].norm(), Env_right[1].norm(), Env_right[0].norm()
#Norm_val_boundry=Env_left[0].product(Env_right[1])
#print Env_left[0].product(Env_right[1])
#for i in xrange(N_x):
# print "NormLeft =", i,Env_left[i].norm(), mps_boundry_left[i].norm()
# print "Normright=", i,Env_right[i].norm(), mps_boundry_right[i].norm()

#for i in xrange(N_x-1):
 #print "Norm:boundryMethod=", i,Env_left[i].product(Env_right[i+1]), mps_boundry_left[i].product(mps_boundry_right[i+1])

for i in xrange(N_x):
 PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)
for i in xrange(N_x):
 PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)

for i in xrange(N_x):
 PEPS_mps_leftU[i]=copy.copy(PEPS_mps_left[i])
for i in xrange(N_x):
 PEPS_mps_rightU[i]=copy.copy(PEPS_mps_right[i])
for i in xrange(N_x):
 PEPS_listtenU[i]=copy.copy(PEPS_listten[i])


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


N_iter_QR=20
N_iter_QR_first=40

j_coupling=3.5
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
    MPS_R, mps_Q, Uni_list , Fidel_val= UD.QR_canon(PEPS_mps_left[q], N_iter, Dp, D,Opt)
    PEPS_mps_left[q]=mps_Q*1.0
    PEPS_listten[q]=UD.mps_to_tensor_left(PEPS_mps_left[q], N_x, N_y, D, d, q)
    PEPS_mps_right[q]=UD.tensor_to_mps_right(PEPS_listten[q], N_x, N_y)

    PEPS_mps_left[q+1]=UD.absorption_left(PEPS_mps_left[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p, Trun_list)
    PEPS_listten[q+1]=UD.mps_to_tensor_left(PEPS_mps_left[q+1], N_x, N_y, D, d, q+1)
    PEPS_mps_right[q+1]=UD.tensor_to_mps_right(PEPS_listten[q+1], N_x, N_y)
    print "norm_left", q, MPS_R.norm(), Fidel_val
    if q!=0:
     Q_norm_list.append(Fidel_val)



#############QR-right#####################################################
################Right-Move######################
   for q in reversed(xrange(1,N_x)):
    #print "q_right",q
    if q==N_x-1:
     Dp=d
    else:
     Dp=D*d

    PEPS_listten[q], PEPS_listten[q-1], E_coulmn[q], E_row[q-1]=UD.Update_E_H( PEPS_listten[q], PEPS_listten[q-1], U_ham, U_ham0, U_hamN, q, N_x, D, d, N_y, H, H0, HN, chi_express)
    #print E_coulmn[q]
    PEPS_mps_left[q]=UD.tensor_to_mps_left(PEPS_listten[q], N_x, N_y)
    PEPS_mps_right[q]=UD.tensor_to_mps_right(PEPS_listten[q], N_x, N_y)
    PEPS_mps_left[q-1]=UD.tensor_to_mps_left(PEPS_listten[q-1], N_x, N_y)
    PEPS_mps_right[q-1]=UD.tensor_to_mps_right(PEPS_listten[q-1], N_x, N_y)

    MPS_R, mps_Q, Uni_list ,Fidel_val= UD.QR_canon( PEPS_mps_right[q], N_iter, Dp, D, Opt)


    if q!=N_x-1:
     Q_norm_list.append(Fidel_val)

    PEPS_mps_right[q]=mps_Q*1.0
    PEPS_listten[q]=UD.mps_to_tensor_right( PEPS_mps_right[q], N_x, N_y, D, d, q)
    PEPS_mps_left[q]=UD.tensor_to_mps_left( PEPS_listten[q], N_x, N_y)

    print "norm_right", q, MPS_R.norm(), Fidel_val

    #print "WOOO"
    PEPS_mps_right[q-1]=UD.absorption_right(PEPS_mps_right[q-1],q-1, MPS_R, D, d, N_x, N_y,chi_p,Trun_list)
    #print "Done"

    PEPS_listten[q-1]=UD.mps_to_tensor_right(PEPS_mps_right[q-1], N_x, N_y, D, d, q-1)
    PEPS_mps_left[q-1]=UD.tensor_to_mps_left(PEPS_listten[q-1], N_x, N_y)
    #print "norm_right", q, MPS_R.norm()#, (MPS_R.norm()- Norm_val_boundry)/Norm_val_boundry

   #print "q_right", 0

   PEPS_mps_right[0]=UD.update_energy(PEPS_mps_right[0], U_ham,U_ham0,U_hamN, d, D, 0, N_x)
   PEPS_listten[0]=UD.mps_to_tensor_right(PEPS_mps_right[0], N_x, N_y, D, d, 0)
   PEPS_mps_left[0]=UD.tensor_to_mps_left(PEPS_listten[0], N_x, N_y)
   E_val=UD.energy_coulmn_val(PEPS_mps_right[0], H, H0, HN, d, D, 0, N_x)
   E_coulmn[0]=E_val
   #print  "E_coulmn",  E_val

   count=0
   for i in xrange(N_x):
    #print "i", i, E_coulmn[i], sum(E_coulmn[i])
    count=sum(E_coulmn[i])+count

   count1=0
   for i in xrange(N_x-1):
    #print "i", i, E_row[i], sum(E_row[i])
    count1=sum(E_row[i])+count1
   
   #print "E_t", count1, count,(count1+count)/(N_x*N_y)



#   #print "Energy_cal"
###############QR-Left#####################################################
#   for q in xrange(0,N_x-1):
#    if q==0:
#     Dp=d
#    else:
#     Dp=D*d
#    MPS_R, mps_Q, Uni_list , Fidel_val= UD.QR_canon(PEPS_mps_left[q], N_iter, Dp, D,Opt)
#    PEPS_mps_left[q]=mps_Q*1.0
#    PEPS_listten[q]=UD.mps_to_tensor_left(PEPS_mps_left[q], N_x, N_y, D, d, q)
#    PEPS_mps_right[q]=UD.tensor_to_mps_right(PEPS_listten[q], N_x, N_y)

#    PEPS_mps_left[q+1]=UD.absorption_left(PEPS_mps_left[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p)
#    PEPS_listten[q+1]=UD.mps_to_tensor_left(PEPS_mps_left[q+1], N_x, N_y, D, d, q+1)
#    PEPS_mps_right[q+1]=UD.tensor_to_mps_right(PEPS_listten[q+1], N_x, N_y)
#    #print "norm_left", q, MPS_R.norm(), Fidel_val
#    if q!=0:
#     Q_norm_list.append(Fidel_val)


#   #####################   QR-right   ###########################
#   ################   Right-Move   ######################
#   for q in reversed(xrange(1,N_x)):
#    #print "q_right",q
#    if q==N_x-1:
#     Dp=d
#    else:
#     Dp=D*d

#    PEPS_listten[q], PEPS_listten[q-1], E_coulmn[q], E_row[q-1], E_mag_coulmn[q], E_mag_row[q-1]=UD.Energy_E_H( PEPS_listten[q], PEPS_listten[q-1], U_ham, U_ham0, U_hamN, q, N_x, D, d, N_y, H, H0, HN, m_mag, m_mag0, m_magN)

#    PEPS_mps_left[q]=UD.tensor_to_mps_left(PEPS_listten[q], N_x, N_y)
#    PEPS_mps_right[q]=UD.tensor_to_mps_right(PEPS_listten[q], N_x, N_y)
#    PEPS_mps_left[q-1]=UD.tensor_to_mps_left(PEPS_listten[q-1], N_x, N_y)
#    PEPS_mps_right[q-1]=UD.tensor_to_mps_right(PEPS_listten[q-1], N_x, N_y)

#    MPS_R, mps_Q, Uni_list ,Fidel_val= UD.QR_canon( PEPS_mps_right[q], N_iter, Dp, D, Opt)

#    if q!=N_x-1:
#     Q_norm_list.append(Fidel_val)

#    PEPS_mps_right[q]=mps_Q*1.0
#    PEPS_listten[q]=UD.mps_to_tensor_right( PEPS_mps_right[q], N_x, N_y, D, d, q)
#    PEPS_mps_left[q]=UD.tensor_to_mps_left( PEPS_listten[q], N_x, N_y)

#    #print "norm_right", q, MPS_R.norm(), Fidel_val

#    PEPS_mps_right[q-1]=UD.absorption_right(PEPS_mps_right[q-1],q-1, MPS_R, D, d, N_x, N_y,chi_p)
#    PEPS_listten[q-1]=UD.mps_to_tensor_right(PEPS_mps_right[q-1], N_x, N_y, D, d, q-1)
#    PEPS_mps_left[q-1]=UD.tensor_to_mps_left(PEPS_listten[q-1], N_x, N_y)

#   #print "q_right", 0


#   E_val=UD.energy_coulmn_val(PEPS_mps_right[0], H, H0, HN, d, D, 0, N_x)
#   E_coulmn[0]=E_val
#   E_val=UD.energy_coulmn_val(PEPS_mps_right[0], m_mag,m_mag0,m_magN, d, D, 0, N_x)
#   E_mag_coulmn[0]=E_val


#   count=0
#   for i in xrange(N_x):
#    #print "i", i, E_coulmn[i], sum(E_coulmn[i])
#    count=sum(E_coulmn[i])+count

#   count1=0
#   for i in xrange(N_x-1):
#    #print "i", i, E_row[i], sum(E_row[i])
#    count1=sum(E_row[i])+count1
#   #E_0=(count1+count)/(N_x*N_y)
   E_1=(count1+count)/(N_x*N_y)
   print "E_total", count, count1, E_1

   if i_iter==0:
    E_min=E_1*1.0


   count_iter=count_iter+1.0
   E_iter_list.append(E_1)
   N_iter_list.append(count_iter)

#   count=0
#   for i in xrange(N_x):
#    count=sum(E_mag_coulmn[i])+count

#   count1=0
#   for i in xrange(N_x-1):
#    count1=sum(E_mag_row[i])+count1

#   Mag_f=(count1+count)/(N_x*N_y)
#   print "Mag_f", count1, count, Mag_f


#   if E_1>E_0:
#    print "Enegy_break", E_1, E_0, i_iter
##    for i in xrange(N_x):
##     PEPS_mps_left[i]=copy.copy(PEPS_mps_leftU[i])
##    for i in xrange(N_x):
##     PEPS_mps_right[i]=copy.copy(PEPS_mps_rightU[i])
##    for i in xrange(N_x):
##     PEPS_listten[i]=copy.copy(PEPS_listtenU[i])
#    E_1=E_0*1.0
#    #break
#   elif i_iter==0 or E_1 < E_min:
#    E_min=E_1*1.0
#    for i in xrange(N_x):
#     PEPS_mps_leftU[i]=copy.copy(PEPS_mps_left[i])
#    for i in xrange(N_x):
#     PEPS_mps_rightU[i]=copy.copy(PEPS_mps_right[i])
#    for i in xrange(N_x):
#     PEPS_listtenU[i]=copy.copy(PEPS_listten[i])




 h_list.append(h_coupling[1])
 #Mag_f_list.append(Mag_f)
 E_f_list.append(E_1)



file = open("Q_norm.txt", "w")
for index in range(len(N_iter_list)):
  file.write(str(index) + " " + str(Q_norm_list[index])+" "+ "\n")
file.close()


file = open("EnergyIter.txt", "w")
for index in range(len(N_iter_list)):
  file.write(str(N_iter_list[index]) + " " + str(E_iter_list[index])+" "+ "\n")
file.close()



file = open("Trun.txt", "w")
for index in range(len(Trun_list)):
  file.write(str(index) + " " + str(Trun_list[index])+" "+ "\n")
file.close()


#file = open("Mag.txt", "w")
#for index in range(len(h_list)):
#  file.write(str(h_list[index]) + " " + str(Mag_f_list[index])+" "+ "\n")
#file.close()



file = open("Energy.txt", "w")
for index in range(len(h_list)):
  file.write(str(h_list[index]) + " " + str(E_f_list[index])+" "+ "\n")
file.close()



