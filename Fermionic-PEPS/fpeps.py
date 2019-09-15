import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import fUD as UD
from scipy.linalg import block_diag

N_x=4
N_y=4

#No-symmetry
# D=[2]
# chi_boundry=[8]
# chi_single=[8]
# chi_try=[2]
# d_phys=[2]

#Z2-symmetry
D=[3,3]
chi_boundry=[20]
chi_single=[200]
chi_try=[20]
d_phys=[1,1]


#U1-symmetry
#D=[1,2,1,1]
#chi_boundry=[200]
#chi_single=[200]
#chi_try=[20]
#d_phys=[1,1]

interval=+1.0e-2
threshold=+1.0e+1
accuracy=+1.0e-8
Model=["Fer_Z2","off"]               #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1#
N_tebd=100
h_coupling=[1.0, 0.10]
Sys=['Fer','double' ]                 #Fer, Bos, double, single
start_itebd=1.0
division_itebd=5.0
N_iter_SU=[20, "full"]     #"full" or "QR"


print Model, Sys, "h=", h_coupling, "D=",D, "chi_single", chi_single,"chi_boundry", chi_boundry 

PEPS_mps=[None]*N_x
PEPS_mps_left=[None]*N_x
PEPS_mps_right=[None]*N_x
PEPS_listten=[None]*N_x

PEPS_mps_leftU=[None]*N_x
PEPS_mps_rightU=[None]*N_x
PEPS_listtenU=[None]*N_x


mps_boundry_left=[None]*N_x
mps_boundry_right=[None]*(N_x+1)
mps_boundry_temp=[None]*N_x
mps_boundry_up=[None]*(N_x+1)
mps_boundry_down=[None]*N_x


q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d=UD.full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_phys)


#print q_D#, q_chi_boundry, q_chi_single, q_chi_try, q_d

for i in xrange(N_x):
 PEPS_listten[i]=UD.Init_PEPS( N_x, q_D, q_d, i)
 #PEPS_listtenU[i]=UD.Init_PEPS( N_x, q_D, q_d, i)


Landa_col=UD.Landa_f_col( q_D, N_x)
Landa_row=UD.Landa_f_row( q_D, N_x)


#UD.Reload_Landa_row(Landa_row, N_x)
#UD.Reload_Landa_col(Landa_col, N_x)
#UD.Reload_Gamma(PEPS_listten, N_x)
#UD.Reload_f(PEPS_listten, N_x)


#PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, h_coupling, q_d, Model, N_x, q_D,q_chi_try, threshold, interval,Sys)

UD.Store_Gamma( PEPS_listten, N_x)
UD.Store_Landa_row( Landa_row, N_x)
UD.Store_Landa_col( Landa_col, N_x)

#PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)


UD.Reload_f(PEPS_listten, N_x)
PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d, threshold, interval,Sys)
#UD.Store_f(PEPS_listten, N_x)


#PEPS_listten=UD.increase_bond(PEPS_listten,D, N_x)

#PEPS_listten=UD.All_dist(PEPS_listten,N_x, q_D)

# for i in xrange(N_x+1):
#  for j in xrange(N_y):
#    print "row", i, j , Landa_row[i][j]
# 
# for i in xrange(N_x):
#  for j in xrange(N_y+1):
#   print "col", i, j , Landa_col[i][j]





if Sys[1] is "double":

 #MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_x,Sys)
 #Env_left, Env_right=UD.make_ENVMPS( MPO_Ten, N_x, q_chi_boundry)
 #Env_up, Env_down=UD.make_ENV_updown( MPO_Ten, N_x, q_chi_boundry)

 for i in xrange(N_x-1):
  print "Norm:double_Layer=", i, Env_left[i].product_nonsymm(Env_right[i+1]), chi_boundry

 for i in xrange(N_x-1):
  print "Norm:double_Layer=", i, Env_down[i].product_nonsymm(Env_up[i+1]),chi_boundry

 E_double=UD.cal_energy_double(PEPS_listten, N_x, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
 print "E_double", E_double



if Sys[1] is "single":
 for Location in xrange(N_x):
   peps_l=[]
   for i in xrange(N_x):
    peps_l.append(PEPS_listten[i][Location])
   mps_boundry_down[Location]=UD.make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], q_d, q_chi_single, N_x,Sys)

 for Location in reversed(xrange(1,N_x)):
   peps_l=[]
   for i in xrange(N_x):
    peps_l.append(PEPS_listten[i][Location])
   mps_boundry_up[Location]=UD.make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], q_d, q_chi_single, N_x ,Sys)

 for Location in xrange(N_x):
   mps_boundry_left[Location]=UD.make_Env_singleLayer( PEPS_listten[Location], Location, mps_boundry_left[Location-1], q_d, q_chi_single, N_x,Sys)
 for Location in reversed(xrange(1,N_x)):
   mps_boundry_right[Location]=UD.make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], q_d, q_chi_single, N_x,Sys)
 for i in xrange(N_x-1):
  print "Norm:single_Layer=", i, mps_boundry_left[i].product_nonsymm(mps_boundry_right[i+1]), chi_single

 for i in xrange(N_x-1):
  print "Norm:single_Layer=", i, mps_boundry_down[i].product_nonsymm(mps_boundry_up[i+1]), chi_single

 Energy_val=UD.Energy_cal(PEPS_listten, q_d, q_chi_single, N_x, q_D, Model, h_coupling,Sys)
 print "E_single", Energy_val



Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
count_list=[]
Cond_list=[]
Norm_list_boundry=[]


if Sys[1] is "single":

 UD.TEBD_Full(N_x, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_single, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys)


if Sys[1] is "double":

 UD.TEBD_Full_double( N_x, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_boundry, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys)


Energy_val=UD.Energy_cal(PEPS_listten, q_d, q_chi_single, N_x, q_D, Model, h_coupling,Sys)
print "E_single", Energy_val


file = open("EnergyIter1.txt", "w")
for index in range(len(E_iter_list1)):
  file.write(str(index) + " " + str(E_iter_list1[index])+" "+ "\n")
file.close()


#for i in xrange(N_x):
# for j in xrange(N_x+1):
#   print "col", i,j, Landa_col[i][j]

#for i in xrange(N_x+1):
#  for j in xrange(N_x):
#   print "row", i,j, Landa_row[i][j]


#for i in xrange(N_x):
# for j in xrange(N_x):
#  print i, j , PEPS_listten[i][j].printDiagram()
