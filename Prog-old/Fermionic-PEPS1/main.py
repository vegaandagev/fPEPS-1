import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import UniformDisentangler as UD
from scipy.linalg import block_diag
N_x=4
N_y=4

#D=[2]
#chi_boundry=[6]
#chi_single=[6]
#chi_try=[2]
#d_phys=[2]

# D=[2,2]
# chi_boundry=[40,40]
# chi_single=[40,40]
# chi_try=[20,20]
# d_phys=[1,1]

D=[2,2,2]
chi_boundry=[10,10,10]
chi_single=[10,10,10]
chi_try=[20,20,20]
d_phys=[1,1]

interval=+1.0e-2
threshold=+1.0e+1
accuracy=+1.0e-8
Model="Fer_U1"               #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1#
N_tebd=1
z_coupling=1.0
j_coupling=.100
h_coupling=[z_coupling, j_coupling]
start_itebd=1.0
division_itebd=5.0
N_iter_SF=40

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


#print q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d

for i in xrange(N_x):
 PEPS_listten[i]=UD.Init_PEPS( N_x, q_D, q_d, i)
 #PEPS_listtenU[i]=UD.Init_PEPS( N_x, q_D, q_d, i)

#for i in xrange(N_x):
# for j in xrange(N_x):
#  print i, j , PEPS_listten[i][j].printDiagram()

Landa_col=UD.Landa_f_col( q_D, N_x)
Landa_row=UD.Landa_f_row( q_D, N_x)

# UD.Reload_Landa_row(Landa_row, N_x)
# UD.Reload_Landa_col(Landa_col, N_x)
# UD.Reload_Gamma(PEPS_listten, N_x)

#PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)
#PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d, threshold, interval)
#Landa_col_inv,Landa_row_inv=UD.inv_landa_col_row( Landa_col, Landa_row, N_x)
#PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row_inv, Landa_col_inv,N_x)


PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SF, h_coupling, q_d, Model, N_x, q_D,q_chi_try, threshold, interval)


UD.Store_Gamma( PEPS_listten, N_x)
UD.Store_Landa_row( Landa_row, N_x)
UD.Store_Landa_col( Landa_col, N_x)


PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)


for i in xrange(N_x):
   PEPS_listtenU[i]=[  PEPS_listten[i][j]*1.0  for  j  in  xrange(N_x)  ]


#UD.Reload_f(PEPS_listten, N_x)
#PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d, threshold, interval)
#UD.Store_f(PEPS_listten, N_x)


#PEPS_listten=UD.increase_bond(PEPS_listten,D, N_x)

#PEPS_listten=UD.All_dist(PEPS_listten,N_x, q_D)

PEPS_listtenU=UD.Symmetric_non(PEPS_listten, PEPS_listtenU, N_x)





#for Location in xrange(N_x):
#  #print "Location", Location
#  peps_l=[]
#  for i in xrange(N_x):
#   peps_l.append(PEPS_listtenU[i][Location])
#  mps_boundry_down[Location]=UD.make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], q_d, q_chi_single, N_x)




#for Location in reversed(xrange(1,N_x)):
#  #print "Location", Location
#  peps_l=[]
#  for i in xrange(N_x):
#   peps_l.append(PEPS_listtenU[i][Location])
#  mps_boundry_up[Location]=UD.make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], q_d, q_chi_single, N_x )



#for Location in xrange(N_x):
#  mps_boundry_left[Location]=UD.make_Env_singleLayer( PEPS_listtenU[Location], Location, mps_boundry_left[Location-1], q_d, q_chi_single, N_x)
##   print mps_boundry_left[Location].norm()

#for Location in reversed(xrange(1,N_x)):
#  mps_boundry_right[Location]=UD.make_Env_singleLayer_right( PEPS_listtenU[Location], Location, mps_boundry_right[Location+1], q_d, q_chi_single, N_x)
##   print mps_boundry_right[Location].norm()
 


#for i in xrange(N_x-1):
# print "Norm:single_Layer=", i, mps_boundry_left[i].product_nonsymm(mps_boundry_right[i+1])


#for i in xrange(N_x-1):
# print "Norm:single_Layer=", i, mps_boundry_down[i].product_nonsymm(mps_boundry_up[i+1])

MPO_Ten=UD.make_boundry_MPO(PEPS_listtenU, N_x)
Env_left, Env_right=UD.make_ENVMPS( MPO_Ten, N_x, q_chi_boundry)

for i in xrange(N_x-1):
 print "Norm:double_Layer=", i, Env_left[i].product_nonsymm(Env_right[i+1])




# for Location in xrange(N_x):
#   #print "Location", Location
#   peps_l=[]
#   for i in xrange(N_x):
#    peps_l.append(PEPS_listten[i][Location])
#   mps_boundry_down[Location]=UD.make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], q_d, q_chi_single, N_x)
# 
# 
# 
# 
# for Location in reversed(xrange(1,N_x)):
#   #print "Location", Location
#   peps_l=[]
#   for i in xrange(N_x):
#    peps_l.append(PEPS_listten[i][Location])
#   mps_boundry_up[Location]=UD.make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], q_d, q_chi_single, N_x )
# 
# 
# 
# for Location in xrange(N_x):
#   mps_boundry_left[Location]=UD.make_Env_singleLayer( PEPS_listten[Location], Location, mps_boundry_left[Location-1], q_d, q_chi_single, N_x)
# #   print mps_boundry_left[Location].norm()
# 
# for Location in reversed(xrange(1,N_x)):
#   mps_boundry_right[Location]=UD.make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], q_d, q_chi_single, N_x)
# #   print mps_boundry_right[Location].norm()
#  
# 
# 
# for i in xrange(N_x-1):
#  print "Norm:single_Layer=", i, mps_boundry_left[i].product_nonsymm(mps_boundry_right[i+1])
# 
# 
# for i in xrange(N_x-1):
#  print "Norm:single_Layer=", i, mps_boundry_down[i].product_nonsymm(mps_boundry_up[i+1])

MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_x)
Env_left, Env_right=UD.make_ENVMPS( MPO_Ten, N_x, q_chi_boundry)

for i in xrange(N_x-1):
 print "Norm:double_Layer=", i, Env_left[i].product_nonsymm(Env_right[i+1])



E_double=UD.cal_energy_double(PEPS_listten, N_x, h_coupling, q_d, Model,q_chi_boundry,q_D)
print "E_double", E_double


# Energy_val=UD.Energy_cal(PEPS_listten, q_d, q_chi_single, N_x, q_D, Model, h_coupling)
# print "E_single", Energy_val



Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
count_list=[]
Cond_list=[]
Norm_list_boundry=[]



#UD.TEBD_Full(N_x, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_single, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1)

#UD.TEBD_Full_double( N_x, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_boundry, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1)




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


