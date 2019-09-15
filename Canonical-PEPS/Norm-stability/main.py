import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 
import UniformDisentangler as UD
import time

def QR_function(a_u,b_u,c_u,d_u):

 N_x=8
 N_y=8
 D=2
 d=2
 chi_boundry=10
 chi_p=2*D

 chi_try=6
 interval=+1.0e-2
 threshold=+0.1050e+1
 #accuracy=+1.0e-8

 N_iter=14
 #Opt=["Non","contiguous","cg", 10]           #  [root & None     and     None & contiguous   and   cg&SVD      chi_root] 
 Opt=["Non","Single","cg",20]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
 #Opt=["Non","ACCcontiguous","cg",6]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 
 #Opt=["Non","Ucontiguous","SVD",6]           #  [root & None     and     Single & contiguous, ACCcontiguous   and   cg&SVD] 

 PEPS_mps=[None]*N_x
 PEPS_mps_left=[None]*N_x
 PEPS_mps_right=[None]*N_x
 PEPS_listten=[None]*N_x

 count_list=[]
 Norm_list=[]
 Norm_list_boundry=[]
 MPS_R_right=[]
 MPS_R_left=[]

 PEPS_mps_Qr=[None]*N_x
 PEPS_mps_Ql=[None]*N_x
 PEPS_Q=[None]*N_x

 for i in xrange(N_x):
   PEPS_mps_left[i]=UD.Init_PEPS( N_y, N_x, D, d,i,a_u,b_u,c_u,d_u)

 for i in xrange(N_x):
  PEPS_listten[i]=UD.mps_to_tensor_left(PEPS_mps_left[i], N_x, N_y, D, d, i)

 for i in xrange(N_x):
  PEPS_Q[i]=UD.mps_to_tensor_left(PEPS_mps_left[i], N_x, N_y, D, d, i)

 for i in xrange(N_x):
    PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)


#  for i in xrange(N_x):
#   PEPS_mps_Ql[i]=copy.copy(PEPS_mps_left[i])
#  for i in xrange(N_x):
#   PEPS_mps_Qr[i]=copy.copy(PEPS_mps_right[i])
#  for i in xrange(N_x):
#   PEPS_Q[i]=copy.copy(PEPS_listten[i])

# for i in xrange(N_x):
#  for j in xrange(N_x):
#    PEPS_listten[i][j].randomize()

# UD.Reload_f(PEPS_listten, N_x, N_y)                              

# PEPS_listten, norm_val, count_val=UD.Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval)    #un
# for i in xrange(N_x):                                                                                         #un
#     PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)
# for i in xrange(N_x):
#   PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)

# UD.Store_f(PEPS_listten, N_x, N_y)                                #un

 #UD.Reload_f(PEPS_listten, N_x, N_y)                              
 
 






 #UD.Store_Q(PEPS_Q, N_x, N_y)
 #UD.Reload_Q(PEPS_Q, N_x, N_y)


 for i in xrange(N_x):
  PEPS_listten[i]=UD.mps_to_tensor_left(PEPS_mps_left[i], N_x, N_y, D, d, i)
 for i in xrange(N_x):
    PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_listten[i], N_x, N_y)



 MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
 #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
 Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
 #print  Env_left[0].norm()#, PEPS_listmps[0].norm()

 #print  Env_left[N_x-1].norm()**(0.5), Env_right[0].norm()**(0.5)
 N_boundry=Env_left[0].product(Env_right[1])
 #print "Norm:boundryMethod=",  N_boundry
 #print  "Norm:boundryMethod=",  Env_left[1].product(Env_right[2])


 alpha=0.0
 N_alpha=30
 N_sample=25
 norm_list=[0]*N_alpha
 Qnorm_list=[0]*N_alpha

 norm_list_f=[0]*N_alpha
 Qnorm_list_f=[0]*N_alpha

 norm_list_ff=[0]*N_alpha
 Qnorm_list_ff=[0]*N_alpha


 alpha_list=[None]*N_alpha
 for iter1 in xrange(N_sample):

  for i1 in xrange(N_alpha):
   norm_list_f[i1]=norm_list[i1]+norm_list_f[i1]  
   Qnorm_list_f[i1]=Qnorm_list[i1]+Qnorm_list_f[i1]  


  alpha=0
  for iter in xrange(N_alpha):
   #alpha+=4
   if iter==0:   alpha+=0.0
   else:   alpha+=0.001
   alpha_list[iter]=alpha
   print "iter", iter, "alpha", alpha


   UD.Reload_f(PEPS_listten, N_x, N_y)                  #comm
   UD.Reload_Q(PEPS_Q, N_x, N_y)                        #comm
   PEPS_listten,PEPS_Q=UD.perturb_norm(PEPS_listten,PEPS_Q,alpha,N_x)           #comm


   #print PEPS_listten[0][2].printDiagram()
   MPO_Ten=UD.make_boundry_MPO(PEPS_listten, N_y, N_x)
   #print boundryMPS[1][0].printDiagram(), boundryMPS[1][1].printDiagram()
   Env_left, Env_right=UD.make_ENVMPS(MPO_Ten, N_y, N_x, chi_boundry)
   #print  Env_left[0].norm()#, PEPS_listmps[0].norm()

   #print  Env_left[N_x-1].norm()**(0.5), Env_right[0].norm()**(0.5)
   N_boundry=Env_left[0].product(Env_right[1])
   print "Norm:boundryMethod=",  N_boundry



   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_Q[i][j]                                                     #comm


   for i in xrange(N_x):
    PEPS_mps_left[i]=UD.tensor_to_mps_left(PEPS_Q[i], N_x, N_y)                          #comm


   for i in xrange(N_x):                                                                 #comm
    PEPS_mps_right[i]=UD.tensor_to_mps_right(PEPS_listten[i], N_x, N_y)                       


   PEPS_mps_update=0                                                                     #comm
    ###############Right-Move######################
   for q in reversed(xrange(1,N_x)):

    if q==N_x-1:
     Dp=d
     PEPS_mps_update=PEPS_mps_right[q]
    else:
     Dp=D*d

    MPS_R, mps_Q, Uni_list,Fidel_QR = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)
    PEPS_mps_update=UD.absorption_right(PEPS_mps_right[q-1],q-1, MPS_R, D, d, N_x, N_y,chi_p)
    MPS_R_right.append(MPS_R)



   N_QR=PEPS_mps_update.norm()
#############################################
   print  "Norm:QRMethod_right",  N_QR,  N_boundry





#   PEPS_mps_update=0                                                        #un
#   for q in xrange(0,N_x-1):
#    if q==0:
#     Dp=d
#     PEPS_mps_update=PEPS_mps_left[0]
#    else:
#     Dp=D*d

#    MPS_R, mps_Q, Uni_list,Fidel_QR = UD.QR_canon(PEPS_mps_update, N_iter, Dp, D,Opt)
#    PEPS_Q[q]=UD.mps_to_tensor_left(mps_Q, N_x, N_y, D, d, q)

#    PEPS_mps_update=UD.absorption_left(PEPS_mps_left[q+1],q+1, MPS_R, D, d, N_x, N_y,chi_p)

#   PEPS_Q[N_x-1]=PEPS_mps_update*1.0
#   PEPS_Q[N_x-1]=UD.mps_to_tensor_left(PEPS_mps_update, N_x, N_y, D, d, N_x-1)
#   N_QR=PEPS_mps_update.norm()
#   print  "Norm:QRMethod_left",  N_QR,  N_boundry

#   UD.Store_Q(PEPS_Q, N_x, N_y)                                             #un


   ##############################################
   ######################################
   norm_list[iter]=N_boundry
   Qnorm_list[iter]=N_QR

  if iter1>0:
   norm_list_ff=[ norm_list_f[i]/(iter1)    for i in xrange(len(norm_list_f))   ]  
   Qnorm_list_ff=[ Qnorm_list_f[i]/(iter1)    for i in xrange(len(Qnorm_list_f))   ]  

  file = open("Norm.txt", "w")
  for index in range(N_alpha):
   file.write(str(alpha_list[index]) + " " + str(Qnorm_list_ff[index])+" "+str(norm_list_ff[index])+" "+ "\n")
  file.close()




