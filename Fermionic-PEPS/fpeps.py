import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import fUD as UD
from scipy.linalg import block_diag
import math

L=6.0
N_x=4
N_y=4
#La_S=1.0/((N_x-1)*(N_x-1))
#La_S=1.0/((N_x)*(N_x))
La_S=(L*L)/((N_x+1)*(N_x+1))
#La_S=1.0

#print La_S 
#No-symmetry
D=[2]
chi_boundry=[4]
chi_single=[40]
chi_try=[6]
d_phys=[2]

#Z2-symmetry
D=[2,2]
chi_boundry=[60]
chi_single=[120]
chi_try=[60]
d_phys=[2,2]


#U1-symmetry
# D=[1,2,1]
# chi_boundry=[30]
# chi_single=[50]
# chi_try=[40]
# d_phys=[2,2]

interval=+1.0e-2
threshold=[1.0,1.0e+1]
accuracy=+1.0e-8
Model=["Fer_Z2","off", "EVEN"]               #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1, FFI_Z2 # ODD
N_tebd=[20, "on"]
N_e=20.0
RG_Int_Coupl_N=(-4.0*math.pi)
RG_Int_Coupl_D=0.50-math.log(0.49758*((N_e/(N_x*N_x))**0.5))
RG_Int_Coupl=RG_Int_Coupl_N/RG_Int_Coupl_D
#La_S=1.0
#print  "Atractive_interaction", RG_Int_Coupl
RG_Int_Coupl_final=RG_Int_Coupl*0.0
#RG_Int_Coupl_final=RG_Int_Coupl*1.0

#print  "test",  -0.75+(math.log(2.0)/2.0)

RG_Int_Coupl_final=0
#La_S=0.0
N_particle=4.0
N_dens=N_particle/(N_x*N_x)
E_fermi=((2.0*math.pi*N_particle)**0.5)/N_x
E_unit=0.5*N_particle*E_fermi
print  "Atractive_interaction",  RG_Int_Coupl_final
print  "Fermi_Energy",  E_fermi
print  "E_unit",  E_fermi



#RG_Int_Coupl_final=-5.61570361221
h_coupling=[ -1.0/La_S, RG_Int_Coupl_final, (+4.0/La_S)-4.50, 0.0 ]   #J (kinetic energy), g, J1 (chemical potential), h
#h_coupling=[ -1.0/La_S, -5.6157036122, (+2.0/La_S)+1.40, -0.0 ]   #J (kinetic energy), g, J1 (chemical potential), h
#h_coupling=[-1.0, +2.0, -0.0, 0.0]



Sys=['Fer','single', "QR", 15, "Z2", "TEBD_SVD"]                 #Fer, Bos, double, single, "QR, Inv", "U1, Z2", "TEBD_SVD, Previous"
#La_S
start_itebd=1.0
#start_itebd=La_S
division_itebd=5.0
N_iter_SU=[40, "QR"]     # "full" or "QR"


print Model, Sys, "h=", h_coupling, "D=",D, "chi_single", chi_single,"chi_boundry", chi_boundry, "N_x=", N_x


PEPS_mps=[None]*N_x
PEPS_mps_left=[None]*N_x
PEPS_mps_right=[None]*N_x
PEPS_listten=[None]*N_x

PEPS_mps_leftU=[None]*N_x
PEPS_mps_rightU=[None]*N_x
PEPS_listtenU=[None]*N_x


mps_boundry_left=[None]*N_x
mps_boundry_right=[None]*(N_x+1)
mps_boundry_up=[None]*(N_y+1)
mps_boundry_down=[None]*N_y


q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d, q_d_out=UD.full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_phys, d_phys)

#q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d=UD.full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_phys)


#print q_d
#print q_D

for i in xrange(N_x):
 PEPS_listten[i]=UD.Init_PEPS( N_x,N_y, q_D, q_d, i, Model)
 #PEPS_listtenU[i]=UD.Init_PEPS( N_x, q_D, q_d, i)


#for i in xrange(N_x):
# for j in xrange(N_y):
#  print i, j , PEPS_listten[i][j].printDiagram()



Landa_col=UD.Landa_f_col( q_D, N_x,N_y)
Landa_row=UD.Landa_f_row( q_D, N_x,N_y)


#UD.Reload_Landa_row(Landa_row, N_x, N_y)
#UD.Reload_Landa_col(Landa_col, N_x, N_y)
#UD.Reload_Gamma(PEPS_listten, N_x, N_y)
#UD.Reload_f(PEPS_listten, N_x, N_y)


PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, h_coupling, q_d, Model, N_x, N_y, q_D,q_chi_try, threshold, interval,Sys)

UD.Store_Gamma( PEPS_listten, N_x, N_y)
UD.Store_Landa_row( Landa_row, N_x, N_y)
UD.Store_Landa_col( Landa_col, N_x, N_y)

PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x, N_y)


#UD.Store_f(PEPS_listten, N_x, N_y)
#UD.Reload_f(PEPS_listten, N_x, N_y)
#UD.Reload_fRG(PEPS_listten, N_x, N_y)



#for i in xrange(N_x):
# for j in xrange(N_y):
#   print "peps", i, j , PEPS_listten[i][j].printDiagram()


PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, N_y, q_D, q_chi_try, q_d, threshold, interval,Sys)
print "norm_val", norm_val
#UD.Store_f(PEPS_listten, N_x, N_y)


#PEPS_listten=UD.increase_bond(PEPS_listten,D, N_x, N_y)

#PEPS_listten=UD.All_dist(PEPS_listten,N_x,N_y, q_D)

# for i in xrange(N_x+1):
#  for j in xrange(N_y):
#    print "row", i, j , Landa_row[i][j]
# 
# for i in xrange(N_x):
#  for j in xrange(N_y+1):
#   print "col", i, j , Landa_col[i][j]





if Sys[1] is "double":

# MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_x,N_y,Sys)
# Env_left, Env_right=UD.make_ENVMPS( MPO_Ten, N_x,N_y, q_chi_boundry)
# Env_up, Env_down=UD.make_ENV_updown( MPO_Ten, N_x, N_y, q_chi_boundry)

# for i in xrange(N_x-1):
#  print "Norm:double_Layer=", i, Env_left[i].product_nonsymm(Env_right[i+1]), chi_boundry

# for i in xrange(N_y-1):
#  print "Norm:double_Layer=", i, Env_down[i].product_nonsymm(Env_up[i+1]),chi_boundry

 E_double=UD.cal_energy_double(PEPS_listten, N_x, N_y, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
 print "E_double", E_double, E_double/(N_x*N_y)
 N_double=UD.cal_particle_double(PEPS_listten, N_x, N_y, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
 print "N_double", N_double
 E_final=E_double-(((h_coupling[2]-2.0/La_S))*N_double)
 print "E_final", E_final
 E_final=E_double-(((h_coupling[2]-2.0/La_S))*N_e)
 print "E_final", E_final



if Sys[1] is "single":

# for Location in xrange(N_x):
#   mps_boundry_left[Location]=UD.make_Env_singleLayer( PEPS_listten[Location], Location, mps_boundry_left[Location-1], q_d, q_chi_single, N_y,Sys)
# for Location in reversed(xrange(1,N_x)):
#   mps_boundry_right[Location]=UD.make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], q_d, q_chi_single, N_y,Sys, N_x)
# for i in xrange(N_x-1):
#  print "Norm:single_Layer=", i, mps_boundry_left[i].product_nonsymm(mps_boundry_right[i+1]), chi_single

# for Location in xrange(N_y):
#   peps_l=[]
#   for i in xrange(N_x):
#    peps_l.append(PEPS_listten[i][Location])
#   mps_boundry_down[Location]=UD.make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], q_d, q_chi_single, N_x,Sys)

# for Location in reversed(xrange(1,N_y)):
#   peps_l=[]
#   for i in xrange(N_x):
#    peps_l.append(PEPS_listten[i][Location])
#   mps_boundry_up[Location]=UD.make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], q_d, q_chi_single, N_x ,Sys, N_y)

# for i in xrange(N_y-1):
#  print "Norm:single_Layer=", i, mps_boundry_down[i].product_nonsymm(mps_boundry_up[i+1]), chi_single






 Energy_val=UD.Energy_cal(PEPS_listten, q_d, q_chi_single, N_x, N_y, q_D, Model, h_coupling,Sys)
 print "E_single", Energy_val, Energy_val/(N_x*N_y)
 Particle_val=UD.Particle_cal(PEPS_listten, q_d, q_chi_single, N_x, N_y, q_D, Model, h_coupling,Sys)
 print "N_single", Particle_val

 Sz_val=UD.Sz_cal(PEPS_listten, q_d, q_chi_single, N_x, N_y, q_D, Model, h_coupling,Sys)
 print "Sz_single", Sz_val


 E_final=Energy_val-(((h_coupling[2]-4.0/La_S))*Particle_val)
 print "E_final", E_final
 E_final=Energy_val-(((h_coupling[2]-4.0/La_S))*N_e)
 print "E_final", E_final
 print "Eta",  0.50 * math.log(2*abs(E_fermi/E_final))




Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
count_list=[]
Cond_list=[]
Norm_list_boundry=[]



#quit()

if Sys[1] is "single":
 UD.TEBD_Full_Single(N_x,N_y, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_single, q_chi_try, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys)
 E_single=UD.Energy_cal(PEPS_listten, q_d, q_chi_single, N_x,N_y, q_D, Model, h_coupling,Sys)
 print "E_single", E_single
 Sz_val=UD.Sz_cal(PEPS_listten, q_d, q_chi_single, N_x, N_y, q_D, Model, h_coupling,Sys)
 print "Sz_single", Sz_val
 N_double=UD.Particle_cal(PEPS_listten, q_d, q_chi_single, N_x, N_y, q_D, Model, h_coupling,Sys)
 print "N_single", N_double
 E_final=E_single-((h_coupling[2]-4.0/La_S)*N_double)
 print "E_final", E_final
 E_final=E_single-((h_coupling[2]-4.0/La_S)*N_e)
 print "E_final", E_final


if Sys[1] is "double":
 PEPS_listten=UD.TEBD_Full_double( N_x,N_y, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d, q_chi_boundry, q_chi_try, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys)



 E_double=UD.cal_energy_double(PEPS_listten, N_x, N_y, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
 print "E_double", E_double
 N_double=UD.cal_particle_double(PEPS_listten, N_x, N_y, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
 print "N_double", N_double
 E_final=E_double-(((h_coupling[2]-2.0/La_S))*N_double)
 print "E_final", E_final
 E_final=E_double-(((h_coupling[2]-2.0/La_S))*N_e)
 print "E_final", E_final



#  N_double=UD.cal_particle_double_H(PEPS_listten, N_x, N_y, h_coupling, q_d, Model,q_chi_boundry,q_D,Sys)
#  print "NU", N_double
#  E_final=E_double-(h_coupling[2]*N_double)
#  print "E_final", E_final

 
 
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



