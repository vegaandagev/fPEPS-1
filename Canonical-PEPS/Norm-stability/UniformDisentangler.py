import pyUni10 as uni10
import copy
import numpy as np
import scipy as sp
#from numpy import linalg as LA
from scipy import linalg as LA
import MPSclass 
import root
import time


def Store_f(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j].save("Store/Pten/a" + str(i)+str(j))

def Reload_f(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j]=uni10.UniTensor("Store/Pten/a" + str(i)+str(j))

def Store_Q(PEPS_Q, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_Q[i][j].save("Store/Qten/Q" + str(i)+str(j))

def Reload_Q(PEPS_Q, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_Q[i][j]=uni10.UniTensor("Store/Qten/Q" + str(i)+str(j))


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
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np



def perturb_norm(peps,cpeps,alpha,N_x):

 for i in xrange(N_x):
  for j in xrange(N_x):
    peps_mat=peps[i][j].getBlock()
    peps_np=Mat_uni_to_np(peps_mat)

    A_rand=np.random.rand(peps_np.shape[0],peps_np.shape[1])
    norm_val=LA.norm(A_rand)
    A_rand=A_rand*(1.0/abs(norm_val))

    peps_np += alpha*A_rand
    peps_mat=Mat_np_to_Uni(peps_np)
    peps[i][j].putBlock(peps_mat)
    cpeps_mat=cpeps[i][j].getBlock()

    cpeps_np=Mat_uni_to_np(cpeps_mat)

    A_rand=np.random.rand(cpeps_np.shape[0],cpeps_np.shape[1])
    norm_val=LA.norm(A_rand)
    A_rand=A_rand*(+1.0/abs(norm_val))

    cpeps_np += alpha*A_rand
    cpeps_mat=Mat_np_to_Uni(cpeps_np)
    cpeps[i][j].putBlock(cpeps_mat)

 return peps, cpeps



def  Norm_based_on_LR(MPS_R_right, MPS_R_left, D):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 E_a=0
 A_l=[]
 for i in xrange(MPS_R_left.N):

  Ten_L=uni10.UniTensor([MPS_R_left[i].bond(0), bdi, bdi, MPS_R_left[i].bond(2)], "Ten_R")
  Ten_L.putBlock(MPS_R_left[i].getBlock())
  Ten_L.setLabel([-1,-2,0,-3])

  Ten_R=uni10.UniTensor([MPS_R_right[i].bond(0), bdi, bdi, MPS_R_right[i].bond(2)], "Ten_R")
  Ten_R.putBlock(MPS_R_right[i].getBlock())
  Ten_R.setLabel([1,2,0,3])

  A_ten=Ten_R*Ten_L
  A_ten.permute([1,-1,2,-2,3,-3],4)
  A_ten.combineBond([1,-1])
  A_ten.combineBond([2,-2])
  A_ten.combineBond([3,-3])
  A_ten.permute([1,2,3],2)
  
  if i == 0:
    #print A_l[i]
    A_ten.setLabel([-1,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-1,-2,2])
    E_a=A_l_dag*A_ten
    E_a.permute([1,2],1)
    E_a.setLabel([-3,-4])
  elif i == (MPS_R_left.N-1):
    A_ten.setLabel([-3,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-4,-2,1])
    #print A_l[i].printDiagram()
    E_a=A_l_dag*(E_a*A_ten)
  else:
    A_ten.setLabel([-3,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-4,-2,2])
    E_a=A_l_dag*(E_a*A_ten)
    E_a.permute([1,2],1)
    E_a.setLabel([-3,-4])

# list_bond=[]
# for q in xrange(MPS_R_left.N):
#   list_bond.append(A_list[q].bond(2).dim())
# 
# print "hi", max(list_bond)
# 
# mps_R=MPSclass.MPS(A_list[1].bond(1).dim(),max(list_bond),MPS_R_left.N)

# for i in xrange(MPS_R_left.N):
#   mps_R[i]=copy.copy(A_list[i])
    
 return   E_a[0]





def  shrink_mps(mps, Dp, D):
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 A_list=[]
 for i in xrange(0,mps.N,2):
  #print i
  Ten_1=uni10.UniTensor([mps[i].bond(0), bdip, bdi, mps[i].bond(2)], "Ten_R")
  #print mps[i].printDiagram()
  Ten_1.setLabel([1,2,3,4])
  Ten_1.putBlock(mps[i].getBlock())
  Ten_2=uni10.UniTensor([mps[i+1].bond(0), bdip, bdi, mps[i+1].bond(2)], "Ten_R")
  Ten_2.putBlock(mps[i+1].getBlock())
  Ten_2.setLabel([4,-2,-3,-4])
  Results=Ten_1*Ten_2  
  Results.permute([1,2,-2,3,-3,-4],5)  
  Results.combineBond([2,-2])
  Results.combineBond([3,-3])
  Results.combineBond([2,3])
  Results.permute([1,2,-4],2)  
  A_list.append(Results)
    
 list_bond=[]
 for q in xrange(mps.N/2):
   list_bond.append(A_list[q].bond(2).dim())

 mps_R=MPSclass.MPS(A_list[1].bond(1).dim(),max(list_bond),mps.N/2,)

 for i in xrange(mps.N/2):
   mps_R[i]=copy.copy(A_list[i])
    
 return   mps_R


def  absorption_left( PEPS_colmps, q, MPS_R, D, d, N_x, N_y, chi_p):

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 PEPS_ten=mps_to_tensor_left(PEPS_colmps, N_x, N_y, D, d, q)
 A_list=[]

 for i in xrange(len(PEPS_ten)):

   PEPS_ten[i].setLabel([0,1,2,3,4])

   Ten_R=uni10.UniTensor([MPS_R[i].bond(0), bdi, bdi, MPS_R[i].bond(2)], "Ten_R")
   Ten_R.putBlock(MPS_R[i].getBlock())
   Ten_R.setLabel([-1,-2,0,-3])
   #Ten_R.permute([-1,0,-2,-3],3)
   #Ten_R.setLabel([-1,-2,0,-3])

   result=Ten_R*PEPS_ten[i]
   result.permute([-2,-1,1,2,3,-3,4],4)
   result.combineBond([-1,1])
   result.combineBond([-3,4])
   result.permute([-2,-1,2,3,-3],3)
   result.permute([-1,-2,2,3,-3],3)
   result.combineBond([-2,2])
   result.combineBond([-2,3])
   result.permute([-1,-2,-3],2)
   A_list.append(result)

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(A_list[i].bond(2).dim())

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y, )
 for i in xrange(N_y):
     mps_A[i]=copy.copy(A_list[i])
 #print "withouth Truncation", mps_A.norm()
 if mps_A.D>chi_p:
   mps_A=mps_A.appSVD(chi_p)
 #print "After Truncation", mps_A.norm()
 return  mps_A



def  absorption_right( PEPS_colmps, q, MPS_R, D, d, N_x, N_y, chi_p):

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 PEPS_ten=mps_to_tensor_right(PEPS_colmps, N_x, N_y, D, d, q)
 A_list=[]

 for i in xrange(len(PEPS_ten)):

   PEPS_ten[i].setLabel([0,1,2,3,4])

   Ten_R=uni10.UniTensor([MPS_R[i].bond(0), bdi, bdi, MPS_R[i].bond(2)], "Ten_R")
   Ten_R.putBlock(MPS_R[i].getBlock())
   Ten_R.setLabel([-1,-2,3,-3])

   result=Ten_R*PEPS_ten[i]
   result.permute([0,-1,1,2,-2,-3,4],4)
   result.combineBond([-1,1])
   result.combineBond([-3,4])
   result.permute([0,-1,2,-2,-3],3)
   result.permute([-1,-2,2,0,-3],3)
   result.combineBond([-2,2])
   result.combineBond([-2,0])
   result.permute([-1,-2,-3],2)
   A_list.append(result)

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(A_list[i].bond(2).dim())

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y, )
 for i in xrange(N_y):
     mps_A[i]=copy.copy(A_list[i])
 #print "withouth Truncation", mps_A.norm()
 if mps_A.D>chi_p:
   mps_A=mps_A.appSVD(chi_p)
 #print "After Truncation", mps_A.norm()
 return  mps_A


def make_ENVMPS(MPO_Ten, N_y, N_x, chi):

 Env_left=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
  A=copy.copy(MPO_Ten[0][j])
  A.setLabel([1,2,3,4])
  A.permute([2,3,1,4],2)
  A.combineBond([3,1])
  A.permute([2,3,4],2)
  cont_list[j]=copy.copy(A)

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(cont_list[q].bond(2).dim())
 mps_A=MPSclass.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y, )
 for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

 if mps_A.D>chi:
   mps_A=mps_A.appSVD(chi)

 Env_left[0]=mps_A

 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=copy.copy(Env_left[q][j])
   A.setLabel([2,3,4])
   B=copy.copy(MPO_Ten[q+1][j])
   B.setLabel([3,5,6,7])
   Result=A*B
   Result.permute([2,5,6,4,7],3)
   Result.combineBond([2,5])
   Result.combineBond([4,7])
   Result.permute([2,6,4],2)
   Result.setLabel([2,3,4])
   Result.permute([2,3,4],2)
   cont_list[j]=copy.copy(Result)

  list_bond=[]
  for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

  mps_A=MPSclass.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y, )

  for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

  if mps_A.D>chi and q!=(N_x-2):
    #print q, "q"
    mps_A=mps_A.appSVD(chi)

  Env_left[q+1]=mps_A

###############################

 Env_right=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
   A=copy.copy(MPO_Ten[N_x-1][j])
   A.setLabel([1,2,3,4])
   A.permute([2,1,3,4],2)
   A.combineBond([1,3])
   A.permute([2,1,4],2)
   cont_list[j]=copy.copy(A)

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

 mps_A=MPSclass.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y, )
 for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

 if mps_A.D>chi:
   mps_A=mps_A.appSVD(chi)


 Env_right[N_x-1]=mps_A


 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=copy.copy(Env_right[N_x-1-q][j])
   A.setLabel([2,3,4])
   B=copy.copy(MPO_Ten[N_x-2-q][j])
   B.setLabel([5,6,3,7])
   Result=A*B
   Result.permute([6,2,5,7,4],3)
   Result.combineBond([6,2])
   Result.combineBond([7,4])
   Result.permute([6,5,7],2)
   cont_list[j]=copy.copy(Result)
  #print  cont_list[0].bond(1).dim(), cont_list[0].bond(2).dim()
  list_bond=[]
  for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

  mps_A=MPSclass.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y, )
  for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

  if mps_A.D>chi and q != (N_x-2) :
    mps_A=mps_A.appSVD(chi)

  Env_right[N_x-2-q]=mps_A


 return Env_left, Env_right





def  make_boundry_MPO(PEPS_listten, N_y, N_x):

 MPO=[None]*N_x
 for i in xrange(N_x):
   MPO[i]=[None]*N_y

 for i in xrange(N_x):
  for j in xrange(N_y):
   #print i, j, PEPS_listten[i][j].printDiagram()
   
   A=copy.copy(PEPS_listten[i][j])
   A.setLabel([1,2,3,4,5])
   A_t=copy.copy(A)
   A_t.setLabel([-1,-2,3,-4,-5])
   A=A*A_t
   A.permute([1,-1,2,-2,4,-4,5,-5],4)
   A.combineBond([1,-1])
   A.combineBond([2,-2])
   A.combineBond([4,-4])
   A.combineBond([5,-5])
   A.permute([1,2,4,5],2)
   MPO[i][j]=copy.copy(A)

 return   MPO


def mps_to_tensor_left(PEPS_listmps, N_x, N_y, D, d,q):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 B_list=[None]*N_y

 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 elif q==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 return   B_list






def   mps_to_tensor_right( PEPS_listmps, N_x, N_y, D, d, q):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 B_list=[None]*N_y

 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 elif q==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     A.permute([0,1,2,3,4],3)
     B_list[i]=A
 return   B_list




def  tensor_to_mps_left(B_list, N_x,N_y):

 A_list=[None]*N_y
 for i in xrange(N_y):
  A=copy.copy(B_list[i])
  A.setLabel([0,1,2,3,4])
  A.permute([1,0,2,3,4],4)
  A.combineBond([0,2])
  A.combineBond([0,3])
  A.permute([1,0,4],2)
  A_list[i]=A

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
  mps_A[i]=copy.copy(A_list[i])
 return   mps_A



def  tensor_to_mps_right(B_list, N_x,N_y):

 A_list=[None]*N_y
 for i in xrange(N_y):
  A=copy.copy(B_list[i])
  A.setLabel([0,1,2,3,4])
  A.permute([1,3,2,0,4],4)
  A.combineBond([3,2])
  A.combineBond([3,0])
  A.permute([1,3,4],2)
  A_list[i]=A

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
  mps_A[i]=copy.copy(A_list[i])
 return   mps_A




def Init_PEPS( N_y, N_x, D, d, q, a_u,b_u,c_u,d_u):
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 A_list=[None]*N_y
 B_list=[None]*N_y


 a_ising =copy.copy(a_u)#uni10.UniTensor("test_tensors/D3a_u")
 a_ising1=copy.copy(c_u)#uni10.UniTensor("test_tensors/D3c_u")
 a_ising2=copy.copy(b_u)#uni10.UniTensor("test_tensors/D3b_u")
 a_ising3=copy.copy(d_u)#uni10.UniTensor("test_tensors/D3d_u")

# a_ising=uni10.UniTensor("test_tensors/D2a_u")
# a_ising1=uni10.UniTensor("test_tensors/D2c_u")
# a_ising2=uni10.UniTensor("test_tensors/D2b_u")
# a_ising3=uni10.UniTensor("test_tensors/D2d_u")

 a_ising.setLabel([0,1,2,3,4])
 a_ising.permute([1,2,0,3,4],3)
 a_ising.setLabel([0,1,2,3,4])

 a_ising1.setLabel([0,1,2,3,4])
 a_ising1.permute([1,2,0,3,4],3)
 a_ising1.setLabel([0,1,2,3,4])


 a_ising2.setLabel([0,1,2,3,4])
 a_ising2.permute([1,2,0,3,4],3)
 a_ising2.setLabel([0,1,2,3,4])

 a_ising3.setLabel([0,1,2,3,4])
 a_ising3.permute([1,2,0,3,4],3)
 a_ising3.setLabel([0,1,2,3,4])

 #print a_ising3.printDiagram()


 bdi_1 = uni10.Bond(uni10.BD_IN, 1)
 bdo_1 = uni10.Bond(uni10.BD_IN, D)
 Tem1=uni10.UniTensor([bdi_1, bdo_1])
 Tem1.identity()
 Tem1.setLabel([1,-1])

 Tem10=uni10.UniTensor([bdi_1, bdo_1])
 Tem10.identity()
 Tem10.setLabel([0,-10])

 Tem4=uni10.UniTensor([bdi_1, bdo_1])
 Tem4.identity()
 Tem4.setLabel([4,-4])


 Tem3=uni10.UniTensor([bdi_1, bdo_1])
 Tem3.identity()
 Tem3.setLabel([3,-3])


 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")
     A.randomize()

     A=copy.copy(a_ising1)
     A.setLabel([-10,-1,2,3,4])
     A=A*Tem1
     A=A*Tem10
     A.permute([0,1,2,3,4],3)

     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
     #print A_list[i].printDiagram()
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.randomize()
     #A.orthoRand()

     A=copy.copy(a_ising)
     A.setLabel([-10,1,2,3,-4])
     A=A*Tem4
     A=A*Tem10
     A.permute([0,1,2,3,4],3)



     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
     A.randomize()
     #A.orthoRand()

     if q % 2==0:
      A=copy.copy(a_ising1)
     else:
      A=copy.copy(a_ising)

     A.setLabel([-10,1,2,3,4])
     A=A*Tem10
     A.permute([0,1,2,3,4],3)

     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A

 elif q==N_x-1:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
     A.randomize()
     #A.orthoRand()


     A=copy.copy(a_ising3)
     A.setLabel([0,-1,2,-3,4])
     A=A*Tem1
     A=A*Tem3
     A.permute([0,1,2,3,4],3)


     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
     #print A_list[i].printDiagram()
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.randomize()

     A=copy.copy(a_ising2)

     A.setLabel([0,1,2,-3,-4])
     A=A*Tem4
     A=A*Tem3
     A.permute([0,1,2,3,4],3)



     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   else:
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
     A.randomize()
     #A.orthoRand()

     if q % 2==0:
      A=copy.copy(a_ising3)
     else:
      A=copy.copy(a_ising2)

     A.setLabel([0,1,2,-3,4])
     A=A*Tem3
     A.permute([0,1,2,3,4],3)


     B_list[i]=copy.copy(A)
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
 else:
  if q%2==0:
   A_list=[None]*N_y
   for i in xrange(N_y):
    if i == 0:
      A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
      A.randomize()
      #A.orthoRand()


      A=copy.copy(a_ising1)
      A.setLabel([0,-1,2,3,4])
      A=A*Tem1
      A.permute([0,1,2,3,4],3)


      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A
      #print A_list[i].printDiagram()
    elif i ==(N_y-1):
      A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
      A.randomize()
      #A.orthoRand()


      A=copy.copy(a_ising)
      A.setLabel([0,1,2,3,-4])
      A=A*Tem4
      A.permute([0,1,2,3,4],3)



      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A
    else:
      A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
      A.randomize()
      #A.orthoRand()

      if q % 2==0:
       A=copy.copy(a_ising1)
      else:
       A=copy.copy(a_ising)

      #print q, "Hi", A.printDiagram(), a_ising.printDiagram()

      A.setLabel([0,1,2,3,4])

      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A
  else :
   A_list=[None]*N_y
   for i in xrange(N_y):
    if i == 0:
      A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
      A.randomize()
      #A.orthoRand()


      A=copy.copy(a_ising3)
      A.setLabel([0,-1,2,3,4])
      A=A*Tem1
      A.permute([0,1,2,3,4],3)


      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A
      #print A_list[i].printDiagram()
    elif i ==(N_y-1):
      A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
      A.randomize()
      #A.orthoRand()


      A=copy.copy(a_ising2)
      A.setLabel([0,1,2,3,-4])
      A=A*Tem4
      A.permute([0,1,2,3,4],3)



      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A
    else:
      A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
      A.randomize()
      #A.orthoRand()

      if q % 2==0:
       A=copy.copy(a_ising3)
      else:
       A=copy.copy(a_ising2)

      #print q, "Hi", A.printDiagram(), a_ising.printDiagram()

      A.setLabel([0,1,2,3,4])

      B_list[i]=copy.copy(A)
      A.setLabel([0,1,2,3,4])
      A.permute([1,0,2,3,4],2)
      A.combineBond([0,2])
      A.combineBond([0,3])
      A.permute([1,0,4],2)
      A_list[i]=A


 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
     mps_A[i]=copy.copy(A_list[i])

 norm=0.09*mps_A.norm()
 for q in xrange(len(B_list)):
   B_list[q]=B_list[q]*(1/(norm**(0.5/N_y)))






 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   else:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A

 elif q==N_x-1:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.randomize()
     #A.orthoRand()
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   else:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A
   else:
     A=copy.copy(B_list[i])
     A.setLabel([0,1,2,3,4])
     A.permute([1,0,2,3,4],2)
     A.combineBond([0,2])
     A.combineBond([0,3])
     A.permute([1,0,4],2)
     A_list[i]=A

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
     A_list[i].randomize()
     mps_A[i]=copy.copy(A_list[i])

 return mps_A



def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)
  #print invLt[0], invLt[1], invLt[2], invLt[3]
  for i in xrange(D):
      invL2[i] = 0 if ((invLt[i].real) < 1.0e-10) else (1.00 / (invLt[i].real))

  invLanda2.putBlock(qnum,invL2)
 return invLanda2


#########  prerequisite functions  #############
def   setTruncation1(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1),theta.bond(2),bdo_mid])
    GB.assign([bdi_mid,theta.bond(3),theta.bond(4),theta.bond(5)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA

def   sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs



    
    
    
  
  
 
 


# 
# def Init_A_PEPS(N,bdi,bdi1,bdip,bdo,bdo1,bdop):
#  A_list=[None]*N
#  for i in xrange(N):
#   if i == 0:
#     A=uni10.UniTensor([bdi1,bdip,bdi,bdo], "A_first")
#     A.randomize()
#     #A.orthoRand()
#     A.setLabel([0,1,2,3])
#     A.combineBond([1,2])
#     #A.permute([1,0,4],2)
#     A_list[i]=A
#     #print A_list[i].printDiagram()
#   elif i ==(N-1):
#     A=uni10.UniTensor([bdi,bdip,bdi,bdo1], "A_end")
#     A.randomize()
#     #A.orthoRand()
#     A.setLabel([0,1,2,3])
#     A.combineBond([1,2])
#     A_list[i]=A
#     #print A_list[i].printDiagram()
#   else:
#     A=uni10.UniTensor([bdi,bdip,bdi,bdo], "A_middle")
#     A.randomize()
#     #A.orthoRand()
#     A.setLabel([0,1,2,3])
#     A.combineBond([1,2])
#     A_list[i]=A
#     #print A_list[i].printDiagram()
# 
#  mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N, )
#  for i in xrange(N):
#    mps_A[i]=copy.copy(A_list[i])
#  return mps_A


def svd_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1)])
    LA=uni10.UniTensor([bd1,theta.bond(1)])
    GB=uni10.UniTensor([bd1,theta.bond(1)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity1(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())
    #print bdi1
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo1])
    LA=uni10.UniTensor([bdi1,bdo1])
    GB=uni10.UniTensor([bdi1,theta.bond(2),theta.bond(3)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity2(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())

    if bdi.dim()<=bdi1.dim():
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB




def svd_parity3(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())

    if bdi.dim()<=bdi1.dim():
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB



# def svd_parity3(theta):
# 
#  bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
#  bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
#  bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
#  #bd4=uni10.Bond(uni10.BD_IN,theta.bond(7).Qlist())
# 
#  GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
#  LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
#  GB=uni10.UniTensor([bd1,bd2,bd3,bd4,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
# 
#  svds = {}
#  blk_qnums = theta.blockQnum()
#  dim_svd=[]
#  for qnum in blk_qnums:
#      svds[qnum] = theta.getBlock(qnum).svd()
#      GA.putBlock(qnum, svds[qnum][0])
#      LA.putBlock(qnum, svds[qnum][1])
#      GB.putBlock(qnum, svds[qnum][2])
# 
# #    print LA
#  return GA, LA, GB


def Cost_Function(T0, U, U_transpose):

 #print  U_transpose*U
 T0p=copy.copy(T0)
 T0.setLabel([0,1,2,3])
 T0p.setLabel([0,1,2,3])
 Norm0=(T0p*T0)

 T0.setLabel([0,1,2,3])
 T0p.setLabel([0,4,2,3])
 U.setLabel([1,-3])
 U_transpose.setLabel([4,-3])

 Norm1=(T0*(U_transpose*U))*(T0p)
 #print Norm1[0], Norm0[0]
 return Norm0[0]-Norm1[0]


def Env_U(T0, U, U_transpose):
 T0p=copy.copy(T0)
 T0.setLabel([0,1,2,3])
 T0p.setLabel([0,4,2,3])
 U.setLabel([1,-3])
 U_transpose.setLabel([4,-3])
 #print T0.printDiagram(), T0p.printDiagram(), U_transpose.printDiagram()
 Y=(T0*U_transpose)*(T0p)
 Y.permute([1,-3],1)
 return Y



def Env_U(T0, U, U_transpose):
 T0p=copy.copy(T0)
 T0.setLabel([0,1,2,3])
 T0p.setLabel([0,4,2,3])
 U.setLabel([1,-3])
 U_transpose.setLabel([4,-3])
 #print T0.printDiagram(), T0p.printDiagram(), U_transpose.printDiagram()
 Y=(T0*U_transpose)*(T0p)
 Y.permute([1,-3],1)
 return Y


def Cost_Function_root(T0,T1,t0,t1,Uni):

  T0.setLabel([1,5,7,-2])
  T1.setLabel([-2,6,8,4])

  t0.setLabel([1,2,7,-3])
  t1.setLabel([-3,3,8,4])
  Uni.setLabel([5,6,2,3])

  Norm=((t0*Uni)*t1)*(T0*T1)
  #print Norm
  return  Norm[0]


def Env_U_root(T0,T1,t0,t1,Uni):
  T0.setLabel([1,5,7,-2])
  T1.setLabel([-2,6,8,4])

  t0.setLabel([1,2,7,-3])
  t1.setLabel([-3,3,8,4])
  Uni.setLabel([5,6,2,3])

  Norm=((t0)*t1)*(T0*T1)
  Norm.permute([5,6,2,3],2)
  return  Norm


def optimized_root(T0,T1,t0,t1, Uni, location):
  U_update=copy.copy(Uni)
  U=Uni
  E2=10
  fidel=0
  for i in xrange(2):
    E1=Cost_Function_root(T0,T1,t0,t1,U)
    #print "E1", E1
    if E1< 1.0e-14:
      #print "E1< 1.0e-14"
      break
    #print i, E1, E2, abs((E2-E1)/E1)
    fidel=E1
    if E1<E2 or i is 0:
     U_update=copy.copy(U)
     if abs((E2-E1)/E1) < 1.0e-7:
      #print E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
    Y=Env_U_root( T0, T1, t0, t1, U)
    #print Y.printDiagram() 
    svd=Y.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    U.putBlock(temporary_matrix)
    E2=copy.copy(E1)

  E1=Cost_Function_root(T0,T1,t0,t1,U_update)
  #print "simpleRoot-Dis=", fidel
  return U_update



def optimized_U(T0, U, U_transpose):
  U_update=copy.copy(U)
  E2=10
  fidel=0
  for i in xrange(30):
    E1=Cost_Function(T0, U, U_transpose)
    if E1< 1.0e-14:
      #print "E1< 1.0e-14"
      break
    #print i, E1, E2, abs((E2-E1)/E1)
    fidel=E1
    if E1<E2 or i is 0:
     U_update=copy.copy(U)
     if abs((E2-E1)/E1) < 1.0e-7:
      #print E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
    Y=Env_U(T0, U, U_transpose)
    #print Y.printDiagram() 
    svd=Y.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    U.putBlock(temporary_matrix)
    U.setLabel([1,-3])
    U_transpose=copy.copy(U)
    U_transpose.setLabel([4,-3])
    #print   U_transpose*U
    E2=copy.copy(E1)

  U_update.setLabel([1,-3])
  U_update_transpose=copy.copy(U_update)
  U_update_transpose.setLabel([4,-3])
  E1=Cost_Function(T0, U, U_transpose)
  #print "simple-Dis=", fidel
  return U_update

def  Cost_FunctionAll_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])


 Iso0.setLabel([2,10])
 Iso1.setLabel([5,11])
 #Iso0t.setLabel([10,8])
 #Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0))*((Iso1)))*(U))*(t0*t1)

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0))*((Iso1))))*(t0*t1)

 Norm2=(t0*t1)*(t0*t1)

 return Norm2[0]+Norm1[0]-2*Norm0[0]

def  Env_Iso0_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])

 Iso0.setLabel([2,10])
 Iso1.setLabel([5,11])
 #Iso0t.setLabel([10,8])
 #Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0))*((Iso1)))*(U))*(t0*t1)
 Norm0.permute([2,10],1)
 #print Norm0

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=(((T1*T0))*Iso1)*(t0*t1)
 Norm1.permute([2,10],1)

 #print Norm1
 #Norm2=(t0*t1)*(t0*t1)
 return Norm1+(-2)*Norm0

def  Env_Iso1_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])


 Iso0.setLabel([2,10])
 Iso1.setLabel([5,11])
 #Iso0t.setLabel([10,8])
 #Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0)))*(U))*(t0*t1)
 Norm0.permute([5,11],1)

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0))))*(t0*t1)
 Norm1.permute([5,11],1)

 return Norm1+(-2)*Norm0

def  Env_uall_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])

 Iso0.setLabel([2,10])
 Iso1.setLabel([5,11])
 #Iso0t.setLabel([10,8])
 #Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0))*((Iso1))))*(t0*t1)
 Norm0.permute([12,13,10,11],2)

 return Norm0


def  Cost_FunctionAll(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])


 Iso0.setLabel([2,8])
 Iso1.setLabel([5,9])
 Iso0t.setLabel([10,8])
 Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t)))*(U))*(t0*t1)

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t))))*(t0*t1)

 Norm2=(t0*t1)*(t0*t1)

 return Norm2[0]+Norm1[0]-2*Norm0[0]



def  Env_Iso0(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])

 Iso0.setLabel([2,8])
 Iso1.setLabel([5,9])
 Iso0t.setLabel([10,8])
 Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0t))*((Iso1*Iso1t)))*(U))*(t0*t1)
 Norm0.permute([2,8],1)
 #print Norm0

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0t))*((Iso1*Iso1t))))*(t0*t1)
 Norm1.permute([2,8],1)

 #print Norm1
 #Norm2=(t0*t1)*(t0*t1)
 return Norm1+(-2)*Norm0

def  Env_Iso1(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])


 Iso0.setLabel([2,8])
 Iso1.setLabel([5,9])
 Iso0t.setLabel([10,8])
 Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1t)))*(U))*(t0*t1)
 Norm0.permute([5,9],1)

 t0.setLabel([1,10,4,14])
 t1.setLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0*Iso0t))*((Iso1t))))*(t0*t1)
 Norm1.permute([5,9],1)

 return Norm1+(-2)*Norm0

def  Env_uall(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.setLabel([1,2,4,3])
 T1.setLabel([3,5,7,6])

 t0.setLabel([1,12,4,14])
 t1.setLabel([14,13,7,6])

 Iso0.setLabel([2,8])
 Iso1.setLabel([5,9])
 Iso0t.setLabel([10,8])
 Iso1t.setLabel([11,9])
 U.setLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t))))*(t0*t1)
 Norm0.permute([12,13,10,11],2)

 return Norm0


def optimized_All_chi(T0,T1,U,Iso0,Iso1):

 E2=10
 #print "simple-Dis=", fidel
 for i in xrange(10):
   E1=Cost_FunctionAll_chi(T0,T1,U,Iso0,Iso1)
   #print 
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Iso0", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso0_update=copy.copy(Iso0)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso0_chi(T0,T1,U,Iso0,Iso1)
   #print Y.printDiagram() 
   svd=Y.getBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso0.putBlock(temporary_matrix)
   Iso0.setLabel([1,2])
   E2=copy.copy(E1)

 E2=10
 #print "simple-Dis=", fidel
 for i in xrange(10):
   E1=Cost_FunctionAll_chi(T0,T1,U,Iso0,Iso1)
   #print 
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Iso1", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso1_update=copy.copy(Iso1)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso1_chi(T0,T1,U,Iso0,Iso1)
   #print Y.printDiagram() 
   svd=Y.getBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso1.putBlock(temporary_matrix)
   Iso1.setLabel([1,2])
   E2=copy.copy(E1)


 E2=10
 fidel=0
 for i in xrange(3):
   E1=Cost_FunctionAll_chi(T0,T1,U,Iso0,Iso1)
   #print 
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Uni", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    U_update=copy.copy(U)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
    #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
    fidel=E2
   Y=Env_uall_chi(T0,T1,U,Iso0,Iso1)
   #print Y.printDiagram() 
   svd=Y.getBlock().svd()
   temporary_matrix=svd[0]*svd[2]
   U.putBlock(temporary_matrix)
   U.setLabel([1,2,3,4])
   E2=copy.copy(E1)

 return Iso0,Iso1,U








def optimized_All(T0,T1,U,Iso0,Iso1):

 E2=10
 for i in xrange(10):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Iso0", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso0_update=copy.copy(Iso0)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso0(T0,T1,U,Iso0,Iso1)
   #print Y.printDiagram() 
   svd=Y.getBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso0.putBlock(temporary_matrix)
   Iso0.setLabel([1,2])
   E2=copy.copy(E1)

 E2=10
 for i in xrange(10):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   if E1< 1.0e-14:
     break
   #print "Iso1", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso1_update=copy.copy(Iso1)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso1(T0,T1,U,Iso0,Iso1) 
   svd=Y.getBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso1.putBlock(temporary_matrix)
   Iso1.setLabel([1,2])
   E2=copy.copy(E1)


 E2=10
 fidel=0
 for i in xrange(3):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   #print 
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Uni", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    U_update=copy.copy(U)
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
    #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
    fidel=E2
   Y=Env_uall(T0,T1,U,Iso0,Iso1)
   svd=Y.getBlock().svd()
   temporary_matrix=svd[0]*svd[2]
   U.putBlock(temporary_matrix)
   U.setLabel([1,2,3,4])
   E2=copy.copy(E1)

 return Iso0,Iso1,U




def Simple_update_root(mps_A, mps_R, N, bdi, bdip, bdo, bdop, Uni_list):

 for i in xrange(N/2):

  t0=uni10.UniTensor([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
  t0.putBlock(mps_R[2*i].getBlock())

  t1=uni10.UniTensor([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
  t1.putBlock(mps_R[2*i+1].getBlock())


  T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.putBlock(mps_A[2*i].getBlock())

  T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.putBlock(mps_A[2*i+1].getBlock())

  U=optimized_root(T0,T1,t0,t1, Uni_list[i*2],2*i)
  Uni_list[i*2].putBlock(U.getBlock())

 return  Uni_list


####@profile
def Simple_update(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list,count_list):

 #bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 #bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 Iso_list=[None]*N
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensor([bdip,bdo])
  Iso_list[i].identity()

 fidel_val=1
 #fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 #print   "fidel_val", fidel_val
 count_list.append(0)
 Fidel_list.append(fidel_val)
 for i in xrange(N):
  #print "Step", i
  T0=uni10.UniTensor([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
  T0.putBlock(mps_A[i].getBlock())
  T0.setLabel([0,1,2,3])
  U=copy.copy(Iso_list[i])
  #U, s, V=svd_parity(Tempo)
  U.setLabel([1,3])
  U_transpose=copy.copy(U)
  U_transpose.setLabel([4,3])
  U=optimized_U(T0, U, U_transpose)
  Iso_list[i]=U

  #fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #print   "fidel_valIso", fidel_val
  count_list.append(count_list[-1]+1)
  Fidel_list.append(fidel_val)

 #print  count_list, Fidel_list
 
 
 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
  t1.randomize()
  svd=t1.getBlock().svd()
  t1.putBlock(svd[0])
  t1.identity()
  Uni_list.append(copy.copy(t1))

 Uni_list_iden=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
  t1.identity()
  Uni_list_iden.append(copy.copy(t1))


 for q in xrange(1):
  for i in xrange(N/2):
   #print "Step", i
   T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T0.putBlock(mps_A[2*i].getBlock())
   T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T1.putBlock(mps_A[2*i+1].getBlock())

   Iso0,Iso1,U=optimized_All(T0,T1,Uni_list[i*2],Iso_list[2*i],Iso_list[2*i+1])
   Iso_list[i*2].putBlock(Iso0.getBlock())
   Iso_list[i*2+1].putBlock(Iso1.getBlock())
   Uni_list[i*2].putBlock(U.getBlock())

#   Iso_list_copy=[None]*(N)
#   for  i  in  xrange(N):
#    Iso_list_copy[i]=copy.copy(Iso_list[i])
#   Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
#   Env_left[0].setLabel([1,2,3,4])
#   Env_right[1].setLabel([1,2,3,4])
#   N_val=Env_left[0]*Env_right[1]
# 
#   Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
#   Env_left[0].setLabel([1,2,3,4])
#   Env_right[1].setLabel([1,2,3,4])
#   D_val=Env_left[0]*Env_right[1]
# 
#   #print   "fidel_valUni", abs(N_val[0])/(abs(D_val[0])**(0.5))
#   count_list.append(count_list[-1]+1)
#   Fidel_list.append(fidel_val)

 return  Iso_list, Uni_list

####@profile
def   Simple_trivial(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list):

 #bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 #bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 Iso_list=[None]*(N)
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensor([bdip,bdo])
  Iso_list[i].identity()
# Iso_list[0].setLabel([1,2])
# A=copy.copy(Iso_list[0])
# A.setLabel([-1,2])
# Result=A*Iso_list[0]
# Result.permute([1,-1],1)
 
 #print Iso_list[0].printDiagram(), Result

 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
  t1.identity()
  Uni_list.append(copy.copy(t1))

 return  Iso_list, Uni_list










def Simple_update_chi(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, chi_canon, Fidel_list, count_list):

 bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 Iso_list=[None]*(N)
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensor([bdip,bdoChi])
  Iso_list[i].identity()
 
 fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 
 
 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdiChi,bdiChi,bdoChi,bdoChi])
  t1.randomize()
  svd=t1.getBlock().svd()
  t1.putBlock(svd[0])
  t1.identity()
  Uni_list.append(copy.copy(t1))

 Uni_list_iden=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdiChi,bdiChi,bdoChi,bdoChi])
  t1.identity()
  Uni_list_iden.append(copy.copy(t1))



 for i in xrange(N/2):
  #print "Step", i
  T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.putBlock(mps_A[2*i].getBlock())
  #T0.setLabel([0,1,2,3])
  T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.putBlock(mps_A[2*i+1].getBlock())
  #T1.setLabel([0,1,2,3])

  Iso0,Iso1,U=optimized_All_chi(T0,T1,Uni_list[i*2],Iso_list[2*i],Iso_list[2*i+1])
  Iso_list[i*2].putBlock(Iso0.getBlock())
  Iso_list[i*2+1].putBlock(Iso1.getBlock())
  Uni_list[i*2].putBlock(U.getBlock())

  Iso_list_copy=[None]*(N)
  for  i  in  xrange(N):
   Iso_list_copy[i]=copy.copy(Iso_list[i])
  Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_left[0].setLabel([1,2,3,4])
  Env_right[1].setLabel([1,2,3,4])
  N_val=Env_left[0]*Env_right[1]

  Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_left[0].setLabel([1,2,3,4])
  Env_right[1].setLabel([1,2,3,4])
  D_val=Env_left[0]*Env_right[1]

  #print   "fidel_valUni", abs(N_val[0])/(abs(D_val[0])**(0.5))
  count_list.append(count_list[-1]+1)
  Fidel_list.append(fidel_val)




 return  Iso_list, Uni_list


def Fidel_Iso(mps_A, Iso_list, N,bdi,bdi1,bdip,bdo,bdo1,bdop):
 Env_left=[None]*(N)
 for i in xrange( N):
   #print i, 2*i+1
   if i == 0:
      #print "i=inside", i, bdip,bdi,mps_A[i].printDiagram()
      t0=uni10.UniTensor([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
      t0.putBlock(mps_A[i].getBlock())
      t0.setLabel([0,1,2,-6])
      T0=copy.copy(t0)
      T0.setLabel([0,4,2,-7])
      Iso=Iso_list[i]
      Iso.setLabel([1,-3])
      Iso_tran=copy.copy(Iso)
      Iso_tran.setLabel([4,-3])
      Result=(Iso*Iso_tran)*(t0*T0)
      Result.permute([-6,-7],2)
      Env_left[i]=Result
   else:
      #print "i=inside", i
      t0=uni10.UniTensor([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
      t0.putBlock(mps_A[i].getBlock())
      t0.setLabel([-8,1,2,-6])
      T0=copy.copy(t0)
      T0.setLabel([-9,4,2,-7])
      Iso=Iso_list[i]
      Iso.setLabel([1,-3])
      Iso_tran=copy.copy(Iso)
      Iso_tran.setLabel([4,-3])
      Env_left[i-1].setLabel([-8,-9])
      Result=((t0*Env_left[i-1])*T0)*(Iso*Iso_tran)
      Result.permute([-6,-7],2)
      Env_left[i]=Result

 #print Env_left[(N/2)-1][0]

 return ((Env_left[N-1][0])**(0.5))



def Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop):
 Env_right=[None]*(N/2)
 for i in xrange((N/2)-1, 0, -1):
   if i == ((N/2)-1):
      t0=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.putBlock(mps_A[2*i+1].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.putBlock(mps_A[2*i].getBlock())
      t1.setLabel([7,5,6,1])
      T0=copy.copy(t0)
      T0.setLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.setLabel([17,14,6,18])

      Iso=Iso_list[2*i+1]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i+1])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i])
      IsoP_t.setLabel([10,9])

      U=Uni_list[2*i]
      U.setLabel([12,13,10,11])

      Up=Uni_list[2*i-1]
      Up.setLabel([15,14,16,12])

      Result=(((((t0*(Iso*Iso_t))*U)*T0)*Up)*(t1*(IsoP*IsoP_t)))*T1
      Result.permute([7,16,15,17],4)
      Env_right[i]=Result
   else:
      #print "i=insideRight", i
      #print "i=inside", i
      t0=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.putBlock(mps_A[2*i+1].getBlock())
      t0.setLabel([1,2,4,-7])
      t1=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.putBlock(mps_A[2*i].getBlock())
      t1.setLabel([7,5,6,1])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([18,-15,4,-17])
      T1=copy.copy(t1)
      T1.setLabel([17,14,6,18])
      Iso=Iso_list[2*i+1]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i+1])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i])
      IsoP_t.setLabel([10,9])

      U=Uni_list[2*i]
      U.setLabel([12,-16,10,11])

      Up=Uni_list[2*i-1]
      Up.setLabel([15,14,16,12])

      Env_right[i+1].setLabel([-7,-16,-15,-17])
      Result=(((((Env_right[i+1]*(t0*(Iso*Iso_t)))*U)*T0)*Up)*(t1*(IsoP*IsoP_t)))*T1
      
      
      #((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Up))*(T1*T0))*(Env_right[i+1])


      Result.permute([7,16,15,17],4)
      #print Result.printDiagram()
      Env_right[i]=Result

 #print "\n"
 Env_left=[None]*(N/2)
 for i in xrange( 0, (N/2)-1):
   #print i, 2*i+1
   if i == 0:
      #print "i=initLeft", i
      t0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t0.putBlock(mps_A[2*i].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t1.putBlock(mps_A[2*i+1].getBlock())
      t1.setLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.setLabel([18,19,6,20])

      Iso=Iso_list[2*i]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*i+1]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i+1])
      IsoP_t.setLabel([10,9])

      U=Uni_list[2*i]
      U.setLabel([12,13,11,10])

      Result=(((t0*(Iso*Iso_t)*T0)*U)*(t1*(IsoP*IsoP_t)))*T1
      Result.permute([7,13,19,20],0)
      #print Result.printDiagram()
      Env_left[i]=Result
   else:
      #print "i=insideLeft", i
      t0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t0.putBlock(mps_A[2*i].getBlock())
      t0.setLabel([-7,2,4,3])
      t1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t1.putBlock(mps_A[2*i+1].getBlock())
      t1.setLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([-20,12,4,18])
      T1=copy.copy(t1)
      T1.setLabel([18,19,6,20])
      Iso=Iso_list[2*i]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*i+1]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i+1])
      IsoP_t.setLabel([10,9])

      U=Uni_list[2*i]
      U.setLabel([21,13,11,10])

      Up=Uni_list[2*i-1]
      Up.setLabel([-19,12,-13,21])


      Env_left[i-1].setLabel([-7,-13,-19,-20])

      Result=(((((Env_left[i-1]*(t0*(Iso*Iso_t)))*T0)*U)*Up)*(t1*(IsoP*IsoP_t)))*T1
      
      
      #((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*((T1*T0)*Up))*(Env_left[i-1])
      Result.permute([7,13,19,20],4)
      #print Result.printDiagram()
      Env_left[i]=Result

 return Env_left, Env_right

def  Obtain_Env_U(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_u=0

 if Location==0:
      #print "Location=UniFirst", Location
      t0=uni10.UniTensor([mps_A[2*Location].bond(0),bdip,bdi,mps_A[2*Location].bond(2)])
      t0.putBlock(mps_A[2*Location].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[2*Location+1].bond(0),bdip,bdi,mps_A[2*Location+1].bond(2)])
      t1.putBlock(mps_A[2*Location+1].getBlock())
      t1.setLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.setLabel([18,19,6,20])

      Iso=Iso_list[2*Location]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*Location])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*Location+1]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*Location+1])
      IsoP_t.setLabel([10,9])

      #U=Uni_list[2*Location]
      #U.setLabel([12,13,11,10])


      Env_right[Location+1].setLabel([7,13,19,20])
      #print Result.printDiagram()
      Env_u=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Env_right[Location+1]))*(T1*T0)
      #print Env_u.printDiagram()
      Env_u.permute([12,13,11,10],2)
 elif Location==(N-2) or Location==(N-3):
      #print "Location=UniEnd", Location
      if Location%2==0:
       i=Location/2
      else:
       i=(Location/2)+1
      t0=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.putBlock(mps_A[2*i+1].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.putBlock(mps_A[2*i].getBlock())
      t1.setLabel([7,5,6,1])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.setLabel([17,14,6,18])

      Iso=Iso_list[2*i+1]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i+1])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i])
      IsoP_t.setLabel([10,9])

      U=Uni_list[2*i]
      U.setLabel([12,13,10,11])

      Up=Uni_list[2*i-1]
      Up.setLabel([15,14,16,12])

      Env_left[i-1].setLabel([7,16,15,17])
      if Location%2==0:
       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Up))*(Env_left[i-1]))*(T1*T0)
       Env_u.permute([12,13,10,11],2)
       #print Result.printDiagram()
      else:
       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*(Env_left[i-1]))*(T1*T0)
       Env_u.permute([15,14,16,12],2)
       #print Result.printDiagram()
 else:
      if (Location%2)==0:
       #print "Location=uniMidel", Location
       i=(Location/2)
       t0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
       t0.putBlock(mps_A[2*i].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
       t1.putBlock(mps_A[2*i+1].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])
       Iso=Iso_list[2*i]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[2*i])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[2*i+1]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[2*i+1])
       IsoP_t.setLabel([10,9])

       #U=Uni_list[2*i]
       #U.setLabel([21,13,11,10])

       Up=Uni_list[2*i-1]
       Up.setLabel([-19,12,-13,21])


       Env_left[i-1].setLabel([-7,-13,-19,-20])
       Env_right[i+1].setLabel([7,13,19,20])

       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Env_right[i+1]))*(((T1*T0)*Env_left[i-1])*Up))
       Env_u.permute([21,13,11,10],2)
       #print Result.printDiagram()
       #Env_u=Result
      else:
       #print "Location=uniMidelodd", Location
       i=((Location+1)/2)
       t0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
       t0.putBlock(mps_A[2*i].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
       t1.putBlock(mps_A[2*i+1].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])
       Iso=Iso_list[2*i]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[2*i])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[2*i+1]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[2*i+1])
       IsoP_t.setLabel([10,9])

       U=Uni_list[2*i]
       U.setLabel([21,13,11,10])

       #Up=Uni_list[2*i-1]
       #Up.setLabel([-19,12,-13,21])


       Env_left[i-1].setLabel([-7,-13,-19,-20])
       Env_right[i+1].setLabel([7,13,19,20])

       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Env_right[i+1]))*((T1*T0)*Env_left[i-1]))
       Env_u.permute([-19,12,-13,21],2)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_u


def  Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_Iso=0

 if Location==0 or Location==1:
      t0=uni10.UniTensor([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
      t0.putBlock(mps_A[0].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
      t1.putBlock(mps_A[1].getBlock())
      t1.setLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.setLabel([18,19,6,20])

      Iso=Iso_list[0]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[0])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[1]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[1])
      IsoP_t.setLabel([10,9])

      U=Uni_list[0]
      U.setLabel([12,13,11,10])
      Env_right[1].setLabel([7,13,19,20])
      if Location%2==0:
        Env_Iso=(((t1*(IsoP*IsoP_t))*((t0*(Iso_t)))*U)*(Env_right[1]))*(T1*T0)
        #print Env_u.printDiagram()
        Env_Iso.permute([2,8],1)
      else:
        Env_Iso=(((t1*(IsoP_t))*((t0*(Iso*Iso_t)))*U)*(Env_right[1]))*(T1*T0)
        #print Env_u.printDiagram()
        Env_Iso.permute([5,9],1)


 elif Location==(N-1) or Location==(N-2):

      #print "Location=IsoEnd", Location
      t0=uni10.UniTensor([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
      t0.putBlock(mps_A[N-1].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
      t1.putBlock(mps_A[N-2].getBlock())
      t1.setLabel([7,5,6,1])

      T0=copy.copy(t0)
      T0.setLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.setLabel([17,14,6,18])

      Iso=Iso_list[N-1]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[N-1])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[N-2]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[N-2])
      IsoP_t.setLabel([10,9])

      U=Uni_list[N-2]
      U.setLabel([12,13,10,11])

      Up=Uni_list[N-3]
      Up.setLabel([15,14,16,12])

      Env_left[((N-2)/2)-1].setLabel([7,16,15,17])
      if Location==(N-1):
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso_t))))*(U*Up))*(Env_left[((N-2)/2)-1]))*(T1*T0)
       Env_Iso.permute([2,8],1)
      else:
       Env_Iso=((((t1*(IsoP_t))*((t0*(Iso*Iso_t))))*(Up*U))*(Env_left[((N-2)/2)-1]))*(T1*T0)
       Env_Iso.permute([5,9],1)
 else:
      if (Location%2)==0:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensor([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t0.putBlock(mps_A[Location].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[Location+1].bond(0),bdip,bdi,mps_A[Location+1].bond(2)])
       t1.putBlock(mps_A[Location+1].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])

       Iso=Iso_list[Location]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[Location+1]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location+1])
       IsoP_t.setLabel([10,9])

       U=Uni_list[Location]
       U.setLabel([21,13,11,10])

       Up=Uni_list[Location-1]
       Up.setLabel([-19,12,-13,21])

       #print (Location/2)-1, (Location/2)+1
       Env_left[(Location/2)-1].setLabel([-7,-13,-19,-20])
       Env_right[(Location/2)+1].setLabel([7,13,19,20])
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso_t)))*(U*Up))*(Env_right[(Location/2)+1]))*(((T1*T0)*Env_left[(Location/2)-1])))
       Env_Iso.permute([2,8],1)


      else:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensor([mps_A[Location-1].bond(0),bdip,bdi,mps_A[Location-1].bond(2)])
       t0.putBlock(mps_A[Location-1].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t1.putBlock(mps_A[Location].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])
       Iso=Iso_list[Location-1]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location-1])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[Location]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location])
       IsoP_t.setLabel([10,9])

       U=Uni_list[Location-1]
       U.setLabel([21,13,11,10])

       Up=Uni_list[Location-2]
       Up.setLabel([-19,12,-13,21])


       Env_left[((Location-1)/2)-1].setLabel([-7,-13,-19,-20])
       Env_right[(Location+1)/2].setLabel([7,13,19,20])

       Env_Iso=((((t1*(IsoP_t))*((t0*(Iso*Iso_t)))*(U*Up))*(Env_right[(Location+1)/2]))*((T1*T0)*Env_left[((Location-1)/2)-1]))
       Env_Iso.permute([5,9],1)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_Iso





def  Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_Iso=0

 if Location==0 or Location==1:
      t0=uni10.UniTensor([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
      t0.putBlock(mps_A[0].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
      t1.putBlock(mps_A[1].getBlock())
      t1.setLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.setLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.setLabel([18,19,6,20])

      Iso=Iso_list[0]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[0])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[1]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[1])
      IsoP_t.setLabel([10,9])

      U=Uni_list[0]
      U.setLabel([12,13,11,10])
      Env_right[1].setLabel([7,13,19,20])
      if Location%2==0:
        Env_Iso=(((t1*(IsoP*IsoP_t))*((t0*(Iso)))*U)*(T1*T0))*(Env_right[1])
        #print Env_u.printDiagram()
        Env_Iso.permute([11,8],1)
      else:
        Env_Iso=(((t1*(IsoP))*((t0*(Iso*Iso_t)))*U)*(T1*T0))*(Env_right[1])
        #print Env_u.printDiagram()
        Env_Iso.permute([10,9],1)


 elif Location==(N-1) or Location==(N-2):

      #print "Location=IsoEnd", Location
      t0=uni10.UniTensor([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
      t0.putBlock(mps_A[N-1].getBlock())
      t0.setLabel([1,2,4,3])
      t1=uni10.UniTensor([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
      t1.putBlock(mps_A[N-2].getBlock())
      t1.setLabel([7,5,6,1])

      T0=copy.copy(t0)
      T0.setLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.setLabel([17,14,6,18])

      Iso=Iso_list[N-1]
      Iso.setLabel([2,8])
      Iso_t=copy.copy(Isop_list[N-1])
      Iso_t.setLabel([11,8])

      IsoP=Iso_list[N-2]
      IsoP.setLabel([5,9])

      IsoP_t=copy.copy(Isop_list[N-2])
      IsoP_t.setLabel([10,9])

      U=Uni_list[N-2]
      U.setLabel([12,13,10,11])

      Up=Uni_list[N-3]
      Up.setLabel([15,14,16,12])

      Env_left[((N-2)/2)-1].setLabel([7,16,15,17])
      if Location==(N-1):
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso))))*(U*Up))*(T1*T0))*(Env_left[((N-2)/2)-1])
       Env_Iso.permute([11,8],1)
      else:
       Env_Iso=((((t1*(IsoP))*((t0*(Iso*Iso_t))))*(Up*U))*(T1*T0))*(Env_left[((N-2)/2)-1])
       Env_Iso.permute([10,9],1)
 else:
      if (Location%2)==0:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensor([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t0.putBlock(mps_A[Location].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[Location+1].bond(0),bdip,bdi,mps_A[Location+1].bond(2)])
       t1.putBlock(mps_A[Location+1].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])

       Iso=Iso_list[Location]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[Location+1]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location+1])
       IsoP_t.setLabel([10,9])

       U=Uni_list[Location]
       U.setLabel([21,13,11,10])

       Up=Uni_list[Location-1]
       Up.setLabel([-19,12,-13,21])

       #print (Location/2)-1, (Location/2)+1
       Env_left[(Location/2)-1].setLabel([-7,-13,-19,-20])
       Env_right[(Location/2)+1].setLabel([7,13,19,20])
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso)))*(U*Up))*(T1*T0))*(((Env_right[(Location/2)+1])*Env_left[(Location/2)-1])))
       Env_Iso.permute([11,8],1)


      else:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensor([mps_A[Location-1].bond(0),bdip,bdi,mps_A[Location-1].bond(2)])
       t0.putBlock(mps_A[Location-1].getBlock())
       t0.setLabel([-7,2,4,3])
       t1=uni10.UniTensor([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t1.putBlock(mps_A[Location].getBlock())
       t1.setLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.setLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.setLabel([18,19,6,20])
       Iso=Iso_list[Location-1]
       Iso.setLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location-1])
       Iso_t.setLabel([11,8])

       IsoP=Iso_list[Location]
       IsoP.setLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location])
       IsoP_t.setLabel([10,9])

       U=Uni_list[Location-1]
       U.setLabel([21,13,11,10])

       Up=Uni_list[Location-2]
       Up.setLabel([-19,12,-13,21])


       Env_left[((Location-1)/2)-1].setLabel([-7,-13,-19,-20])
       Env_right[(Location+1)/2].setLabel([7,13,19,20])

       Env_Iso=((((t1*(IsoP))*((t0*(Iso*Iso_t)))*(U*Up))*(T1*T0))*((Env_right[(Location+1)/2])*Env_left[((Location-1)/2)-1]))
       Env_Iso.permute([10,9],1)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_Iso


def optimize_uni_locally(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
  #print Uni_list[0]
  E2=0
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].identity()   

  Iso_list_copy=[  copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]


  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)

  Env_left1[0].setLabel([1,2,3,4])
  Env_right1[1].setLabel([1,2,3,4])
  Dominator=Env_left1[0]*Env_right1[1]

  Dominator=(abs(Dominator[0]))**(0.5)
  #print Uni_list[0]
  #print  "Dominator", Dominator 
  Fidel=0
  
  for i in xrange(4):
    Env_u=Obtain_Env_U(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_u.setLabel([1,2,3,4])
    Uni_list[Location].setLabel([1,2,3,4])
    R=Uni_list[Location]*Env_u
    E1=R[0]/(Dominator)
    Fidel=E1
    #print "Uni, step", i,  E1, E2, abs((E2-E1)/E1)
    if E1>E2 or i is 0:
     U_update=copy.copy(Uni_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      #print E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'NotoptimizedUni=i, E1, E2=', i,'  ', E1, '   ', E2
     Uni_list[Location].putBlock(U_update.getBlock())
     Fidel=E2
     break
    #Y=Env_U(T0, T1, U, U_transpose)
    #print Y.printDiagram() 
    svd=Env_u.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    Uni_list[Location].putBlock(temporary_matrix)
    E2=copy.copy(E1)
  #print "Location", Location, "FidelUni", Fidel
  return Fidel

def optimize_uni_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right):
  #print "\n", "Location", Location
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].identity()   
  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)

  Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
  Env_u.setLabel([1,2,3,4,5,6,7,8])
  Uni_list[Location].setLabel([1,2,3,4,5,6,7,8])
  R=Uni_list[Location]*Env_u

  Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
  Env_iso1.setLabel([1,2,3,4,5,6])
  Iso_list[Location].setLabel([1,2,3,4,5,6])
  Dominator=Iso_list[Location]*Env_iso1
  #Dominator=(Dominator[0])**(0.5)

  E=0
  Fidel=0
  #print 'E=', E
  
  E2=0.0
  U_update=copy.copy(Uni_list[Location])
  U_first=copy.copy(Uni_list[Location])
  E_previous=0
  for i in xrange(300):
   Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
   Env_u.setLabel([1,2,3,4,5,6,7,8])
   Uni_list[Location].setLabel([1,2,3,4,5,6,7,8])
   R=Uni_list[Location]*Env_u
   #E1=R[0]/Dominator
   E1=1+(-2)*R[0]+Dominator[0]
   Fidel=R[0]/(Dominator[0]**(0.5))

   #print 'E=', E1,  Fidel

   #E_list.append(E1)
   #Count_list.append(count)
   #print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock()*U_update).trace()
   D_u=(-2.0)*Env_u
   D_u_mat=copy.copy(D_u.getBlock())
   D_u_trans=copy.copy(D_u.getBlock())
   D_u_trans.transpose()
   Uni=copy.copy(Uni_list[Location].getBlock())
   #Iso_t.transpose()
   Uni_t=copy.copy(Uni)
   Uni_t.transpose()
   Z_decent_mat=Uni*D_u_trans*Uni+(-1.00)*D_u_mat
   Z_decent=copy.copy(D_u)
   Z_decent.putBlock(Z_decent_mat)
   Z_decent_trans=copy.copy(Z_decent_mat)
   Z_decent_trans.transpose()
   Norm_Z=((Z_decent_trans*Z_decent_mat)) + ((Z_decent_trans*Uni*Uni_t*Z_decent_mat) * (-0.500))
   #print Iso_t*Iso
   Norm_Z=abs(Norm_Z.trace())
    
   if E1<E_previous or i is 0:
    if abs((E_previous-E1)/E1) < 1.0e-8:
     print E_previous, E1, abs((E_previous-E1)/E1), i
     break
   elif E1>E_previous:
     #print "Notoptimized Unitary"
     break
   E_previous=copy.copy(E1)
   if Norm_Z < 1.0e-7:
    #print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    Temporary=Uni+(2.00)*Gamma*Z_decent_mat
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    Uni_list[Location].putBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_u.setLabel([1,2,3,4,5,6,7,8])
    Uni_list[Location].setLabel([1,2,3,4,5,6,7,8])
    R=Uni_list[Location]*Env_u
    #E1=R[0]/Dominator
    E2=1+(-2)*R[0]+Dominator[0]
    #print 'E2=', E2, E1-E2, -Norm_Z*Gamma 
    if E1-E2 >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    #print 'Numbers=', count 
    Temporary=Uni+Gamma*Z_decent_mat
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    Uni_list[Location].putBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_u.setLabel([1,2,3,4,5,6,7,8])
    Uni_list[Location].setLabel([1,2,3,4,5,6,7,8])
    R=Uni_list[Location]*Env_u
    #E1=R[0]/Dominator
    E2=1+(-2)*R[0]+Dominator[0]
    #print 'E2s=', E2, 'Gamma=', Gamma
    #print 'abs=', abs((0.5)*Norm_Z*Gamma) 
    if abs((0.5)*Norm_Z*Gamma) <1.0e-11 or  abs(E1-E2)<1.0e-11 :
     #print 'break, Gamma is too small', 'E1-E2=', abs(E1-E2)
     Temporary=Uni+Gamma*Z_decent_mat
     svd=Temporary.svd()
     Temporary=svd[0]*svd[2]
     break
    #print E1-E2, (-0.5)*Norm_Z*Gamma  
    if E1-E2 < (0.5)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   Temporary=Uni+Gamma*Z_decent_mat
   svd=Temporary.svd()
   Temporary=svd[0]*svd[2]
   Uni_list[Location].putBlock(Temporary)

  #print "Location", Location, "FidelUni", Fidel






def optimize_iso_locally_SD(mps_A, Iso_list,Isop_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right):
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
     Uni_list_iden[i].identity()

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]

  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)

  Env_iso=Obtain_Env_Iso(mps_A,Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso

  Env_iso1=Obtain_Env_Iso(mps_A,Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso1.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso1
  Dominator=(abs(Dominator[0]))**(0.5)
  E=R[0]/Dominator
  #print 'E=', E
  E1=10
  Fidel=0
  E2=0.0
  U_update=copy.copy(Iso_list[Location])
  U_first=copy.copy(Iso_list[Location])
  E_previous=0
  for i in xrange(100):
   Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
   Env_iso.setLabel([1,2])
   Iso_list[Location].setLabel([1,2])
   R=Iso_list[Location]*Env_iso

   Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
   Env_iso1.setLabel([1,2])
   Iso_list[Location].setLabel([1,2])
   Dominator=Iso_list[Location]*Env_iso1
   Dominator=(abs(Dominator[0]))**(0.5)
   #E1=R[0]/Dominator
   E1=1+(-2)*R[0]+Dominator[0]
   Fidel=R[0]/(Dominator[0]**(0.5))
   #print 'E=', E1 , Fidel 
   #E_list.append(E1)
   #Count_list.append(count)
   #print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().getBlock()*U_update).trace()
   D_u=(-2.0)*Env_iso+Env_iso1
   D_u_mat=copy.copy(D_u.getBlock())
   D_u_trans=copy.copy(D_u.getBlock())
   D_u_trans.transpose()
   Iso=copy.copy(Iso_list[Location].getBlock())
   #Iso_t.transpose()
   Iso_t=copy.copy(Iso)
   Iso_t.transpose()
   Z_decent_mat=Iso*D_u_trans*Iso+(-1.00)*D_u_mat
   Z_decent=copy.copy(D_u)
   Z_decent.putBlock(Z_decent_mat)
   Z_decent_trans=copy.copy(Z_decent_mat)
   Z_decent_trans.transpose()
   Norm_Z=((Z_decent_trans*Z_decent_mat)) + ((Z_decent_trans*Iso*Iso_t*Z_decent_mat) * (-0.500))
   #print Iso_t*Iso
   Norm_Z=abs(Norm_Z.trace())
    
   if E1<E_previous or i is 0:
    if abs((E_previous-E1)/E1) < 1.0e-8:
     print E_previous, E1, abs((E_previous-E1)/E1), i
     break
   elif E1>E_previous:
     #print "Notoptimized Iso"
     break
   E_previous=copy.copy(E1)
   if Norm_Z < 1.0e-7:
    #print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   Gamma=1.0
   while Break_loop is 1:
    Temporary=Iso+(2.00)*Gamma*Z_decent_mat
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    Iso_list[Location].putBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_iso.setLabel([1,2,3,4,5,6])
    Iso_list[Location].setLabel([1,2,3,4,5,6])
    R=Iso_list[Location]*Env_iso
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
    Env_iso1.setLabel([1,2,3,4,5,6])
    Iso_list[Location].setLabel([1,2,3,4,5,6])
    Dominator=Iso_list[Location]*Env_iso1
    #E1=R[0]/Dominator
    E2=1+(-2)*R[0]+Dominator[0]
    #print 'E2=', E2, E1-E2, -Norm_Z*Gamma 
    if E1-E2 >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0
   
   Break_loop=1
   while Break_loop is 1:
    #print 'Numbers=', count 
    Temporary=Iso+Gamma*Z_decent_mat
    svd=Temporary.svd()
    Temporary=svd[0]*svd[2]
    Iso_list[Location].putBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_iso.setLabel([1,2,3,4,5,6])
    Iso_list[Location].setLabel([1,2,3,4,5,6])
    R=Iso_list[Location]*Env_iso
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
    Env_iso1.setLabel([1,2,3,4,5,6])
    Iso_list[Location].setLabel([1,2,3,4,5,6])
    Dominator=Iso_list[Location]*Env_iso1
    #E1=R[0]/Dominator
    E2=1+(-2)*R[0]+Dominator[0]
    #print 'E2s=', E2, 'Gamma=', Gamma
    #print 'abs=', abs((0.5)*Norm_Z*Gamma) 
    if abs((0.5)*Norm_Z*Gamma) <1.0e-11 or  abs(E1-E2)<1.0e-11 :
     #print 'break, Gamma is too small', 'E1-E2=', abs(E1-E2)
     Temporary=Iso+Gamma*Z_decent_mat
     svd=Temporary.svd()
     Temporary=svd[0]*svd[2]
     break
    #print E1-E2, (-0.5)*Norm_Z*Gamma  
    if E1-E2 < (0.5)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   Temporary=Iso+Gamma*Z_decent_mat
   svd=Temporary.svd()
   Temporary=svd[0]*svd[2]
   Iso_list[Location].putBlock(Temporary)
  #print "Location", Location, "FidelIso", Fidel
 
 
 
 
def optimize_iso_locallyUniform(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):


  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]

  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  #Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)
  Fidel=0
  E2=0
  for i in xrange(25):
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_iso.setLabel([1,2])
    Iso_list[Location].setLabel([1,2])
    R=Iso_list[Location]*Env_iso

    E1=abs(R[0])**(0.5)
    Fidel=abs(E1)
    #print "Iso, step", i, E1, E2, R[0],Dominator
    if E1>E2 or i is 0:
     U_update=copy.copy(Iso_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      #print "BreakIso", E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'NotoptimizedIso=i, E1, E2=', i,'  ', E1, '   ', E2
     Iso_list[Location].putBlock(U_update.getBlock())
     Iso_list_copy[Location].putBlock(copy.copy(U_update.getBlock()))
     Fidel=E2
     break
    Env_iso=Env_iso
    svd=Env_iso.getBlock().svd()
    temporary_matrix=(1.0)*svd[0]*svd[2]
    Iso_list[Location].putBlock(temporary_matrix)
    Iso_list_copy[Location].putBlock(temporary_matrix)

    E2=copy.copy(E1)

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso

  E1=abs(R[0])**(0.5)
  Fidel=abs(E1)
  #print "Location", Location, "FidelIso", Fidel#, R[0], Dominator
  return Fidel


def optimize_iso_locally(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):

  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
     Uni_list_iden[i].identity()

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]

  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)
  Fidel=0
  E2=0
  for i in xrange(25):
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_iso.setLabel([1,2])
    Iso_list[Location].setLabel([1,2])
    R=Iso_list[Location]*Env_iso

    #Iso_list_copy=[  copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
    Env_iso1.setLabel([1,2])
    Iso_list[Location].setLabel([1,2])
    Dominator=Iso_list[Location]*Env_iso1
    Dominator=(abs(Dominator[0]))**(0.5)

    E1=abs(R[0]/Dominator)
    Fidel=abs(E1)
    #print "Iso, step", i, E1, E2, R[0],Dominator
    if E1>E2 or i is 0:
     U_update=copy.copy(Iso_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      #print "BreakIso", E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'NotoptimizedIso=i, E1, E2=', i,'  ', E1, '   ', E2
     Iso_list[Location].putBlock(U_update.getBlock())
     Iso_list_copy[Location].putBlock(copy.copy(U_update.getBlock()))
     Fidel=E2
     break
    Env_iso=(-2.0)*Env_iso+Env_iso1
    svd=Env_iso.getBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Iso_list[Location].putBlock(temporary_matrix)
    Iso_list_copy[Location].putBlock(temporary_matrix)

    E2=copy.copy(E1)

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso
  Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso1.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso1
  Dominator=(abs(Dominator[0]))**(0.5)

  E1=abs(R[0]/Dominator)
  Fidel=abs(E1)
  #print "Location", Location, "FidelIso", Fidel#, R[0], Dominator
  return  Fidel

def optimize_isop_locally(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].identity()

  Iso_list_copy=[ copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]

  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)
  Fidel=0
  E2=0

  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso
  Dominator=(abs(Dominator[0]))**(0.5)


  for i in xrange(5):
    Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_isop.setLabel([1,2])
    Isop_list[Location].setLabel([1,2])
    R=Isop_list[Location]*Env_isop


    E1=abs(R[0]/Dominator)
    Fidel=E1
    #print "Isop, step", i,  E1,E2, R[0], Dominator
    if E1>E2 or i is 0:
     U_update=copy.copy(Isop_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      #print "BreakP", E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'NotoptimizedIso=i, E1, E2=', i,'  ', E1, '   ', E2
     Isop_list[Location].putBlock(U_update.getBlock())
     Fidel=E2
     break
    Env_isop=Env_isop
    svd=Env_isop.getBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Isop_list[Location].putBlock(temporary_matrix)
    E2=copy.copy(E1)

  Iso_list_copy=[ copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso
  Dominator=(abs(Dominator[0]))**(0.5)
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_isop.setLabel([1,2])
  Isop_list[Location].setLabel([1,2])
  R=Isop_list[Location]*Env_isop

  E1=abs(R[0]/Dominator)
  Fidel=E1

  #print "Location", Location, "FidelIsop", Fidel#, R[0], Dominator





def Optimize_uni_Iso(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right,N_iter,Fidel_list, count_list):



 for q in xrange(N_iter):

  for i in xrange(N-1):
    #print ""i
    Location=i
    fidel_val=optimize_uni_locally(mps_A, Iso_list,Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
    Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
    #optimize_uni_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    #Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi)
    #print "U", fidel_val
    count_list.append(count_list[-1]+1)
    Fidel_list.append(fidel_val)
    
    
    
    
  for i in xrange(N):
  #for i in xrange(0):
   #print i
   Location=i
   fidel_val=optimize_iso_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   optimize_isop_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #optimize_iso_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
   #Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi)
   #print "T", fidel_val
   count_list.append(count_list[-1]+1)
   Fidel_list.append(fidel_val)




def Full_update(mps_A,Iso_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, N_iter, Fidel_list, count_list):

 #Uni_list=[]
 Isop_list=[]
#  for i in xrange(len(Iso_list)-1):
#   t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
#   t1.randomize()
#   svd=t1.getBlock().svd()
#   t1.putBlock(svd[0])
#   t1.identity()
#   Uni_list.append(copy.copy(t1))

 for i in xrange(len(Iso_list)):
  Isop_list.append(copy.copy(Iso_list[i]))

#   for i in xrange(len(Uni_list)):
#    Uni_list[i].randomize()
#   for i in xrange(len(Isop_list)):
#    Isop_list[i].randomize()
#   for i in xrange(len(Isop_list)):
#    Iso_list[i].randomize()

 Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 for i in xrange(len(Env_left)-1):
  Env_left[i].setLabel([1,2,3,4])
  Env_right[i+1].setLabel([1,2,3,4])
  R=Env_left[i]*Env_right[i+1]
  #print "Env", i, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N-1):
  Location=i
  Env_u=Obtain_Env_U(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_u.setLabel([1,2,3,4])
  Uni_list[Location].setLabel([1,2,3,4])
  R=Uni_list[Location]*Env_u
  #print "Env_Uni", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso
  #print "Env_Iso", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_isop.setLabel([1,2])
  Isop_list[Location].setLabel([1,2])
  R=Isop_list[Location]*Env_isop
  #print "Env_Isop", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N_iter):
   Optimize_uni_Iso(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right,N_iter,Fidel_list, count_list)

 return  Iso_list, Isop_list, Uni_list


def  Full_updateUniform(mps_A,Iso_list,N,bdi,bdi1,bdip,bdo,bdo1,bdop,N_iter,Fidel_list, count_list):


 Uni_list=[]
 Isop_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
  t1.randomize()
  svd=t1.getBlock().svd()
  t1.putBlock(svd[0])
  t1.identity()
  Uni_list.append(copy.copy(t1))

 for i in xrange(len(Iso_list)):
  Isop_list.append(copy.copy(Iso_list[i]))



 Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 for i in xrange(len(Env_left)-1):
  Env_left[i].setLabel([1,2,3,4])
  Env_right[i+1].setLabel([1,2,3,4])
  R=Env_left[i]*Env_right[i+1]
  #print "Env", i, (abs(R[0]))**(0.5)

 #print "\n"


 for i in xrange(N-1):
  Location=i
  Env_u=Obtain_Env_U(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_u.setLabel([1,2,3,4])
  Uni_list[Location].setLabel([1,2,3,4])
  R=Uni_list[Location]*Env_u
  #print "Env_Uni", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso
  #print "Env_Iso", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_isop.setLabel([1,2])
  Isop_list[Location].setLabel([1,2])
  R=Isop_list[Location]*Env_isop
  #print "Env_Isop", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  fidel_val=optimize_iso_locallyUniform(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  count_list.append(count_list[-1]+1)
  Fidel_list.append(fidel_val)
  #print "FidelUniform", fidel_val
  
  
 return  Iso_list

def make_MPS_R(mps_A, Iso_list,N_y, bdi, bdi1, bdip, bdo, bdo1, bdop):

 MPS_list=[None]*(N_y)
 for i in xrange(N_y):
  t0=uni10.UniTensor([ mps_A[i].bond(0), bdip, bdi, mps_A[i].bond(2)])
  t0.putBlock( mps_A[i].getBlock() )
  t0.setLabel([1,2,4,3])
  Iso=Iso_list[i]
  Iso.setLabel([2,8])
  Result=Iso*t0
  Result.permute([1,8,4,3],3)
  Result.combineBond([8,4])
  Result.permute([1,8,3],2)
  MPS_list[i]=Result

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(MPS_list[q].bond(2).dim())

 mps_R=MPSclass.MPS(MPS_list[1].bond(1).dim(),max(list_bond),N_y)

 for i in xrange(N_y):
   mps_R[i]=copy.copy(MPS_list[i])

 return   mps_R

##@profile
def  make_MPS_Q(Uni_list):

 U_even=[]
 V_even=[]

 U_odd=[]
 V_odd=[]
 
 #print len(Uni_list)-1
 for i in xrange(((len(Uni_list)-1)/2)+1):
  #print "StepsEven", i, 2*i
  #print Uni_list[0].printDiagram()
  Uni_list[2*i].setLabel([1,2,3,4])
  Uni_0=copy.copy(Uni_list[2*i])
  Uni_0.permute([1,3,2,4],2)
  U, s, V=svd_parity1(Uni_0)
  U.setLabel([1,5,-1])
  s.setLabel([-1,-2])
  V.setLabel([-2,2,6])
  V=V*s
  V.permute([-1,2,6],1)
  
  U_even.append(U)
  U_even.append(V)
  
 for i in xrange((len(Uni_list)-1)/2):
  #print "StepsOdd", 2*i+1
  Uni_list[2*i+1].setLabel([1,2,3,4])
  Uni_0=copy.copy(Uni_list[2*i+1])
  Uni_0.permute([1,3,2,4],2)
  U, s, V=svd_parity1(Uni_0)
  U.setLabel([1,3,-1])
  s.setLabel([-1,-2])
  V.setLabel([-2,2,4])
  V=V*s
  V.permute([-1,2,4],1)
  #print "Steps",  2*i+1, U.printDiagram()
  U_odd.append(U)
  U_odd.append(V)


 A_list=[None]*(len(Uni_list)+1)
 bdi=uni10.Bond(uni10.BD_IN, 1)
 A_temp=uni10.UniTensor([bdi])
 A_temp.identity()
 A_temp.setLabel([-3])
 #print A_temp
 for i in xrange(len(Uni_list)+1):
   if i==0:
    #print "First", i
    U_even[i].setLabel([1,5,-1])
    A=U_even[i]*A_temp
    A.permute([-3,1,5,-1],3)
    A.combineBond([1,5])
    A.permute([-3,1,-1],2)
    A_list[i]=copy.copy(A)
   elif i == (len(Uni_list)):
    #print "Last", i#, len(U_even)
    A=U_even[i]*A_temp
    A.permute([-1,2,6,-3],3)
    A.combineBond([2,6])
    A.permute([-1,2,-3],2)
    A_list[i]=copy.copy(A)
   else:
    if (i%2) != 0:
     #print "odd", i
     U_even[i].setLabel([-1,2,3])
     U_odd[i-1].setLabel([4,2,5])
     #print U_even[i].printDiagram(),U_odd[i-1].printDiagram() 
     A=U_even[i]*U_odd[i-1]
     A.permute([-1,4,3,5],3)
     A.combineBond([4,3])
     A.permute([-1,4,5],2)
     A_list[i]=copy.copy(A)
    else:
     #print "even", i
     U_even[i].setLabel([1,3,2])
     U_odd[i-1].setLabel([4,5,1])
     A=U_even[i]*U_odd[i-1]
     A.permute([4,5,3,2],3)
     A.combineBond([5,3])
     A.permute([4,5,2],2)
     A_list[i]=copy.copy(A)

 list_bond=[]
 for q in xrange(len(A_list)):
   list_bond.append(A_list[q].bond(2).dim())
 #print max(list_bond)
 mps_Q=MPSclass.MPS( A_list[1].bond(1).dim(), max(list_bond), len(A_list))
 for i in xrange(len(A_list)):
   mps_Q[i]=copy.copy(A_list[i])
 
 
 return   mps_Q


def  Uni_update(Isop_list, Uni_list):

 for i in xrange(((len(Uni_list)-1)/2)+1):
  #print "StepsEven", 2*i
  #print Uni_list[0].printDiagram(), Isop_list[0].printDiagram()
  Uni_list[2*i].setLabel([1,2,3,4])
  Isop_list[2*i].setLabel([3,5])
  Isop_list[2*i+1].setLabel([4,6])
  Uni_0=(Uni_list[2*i]*Isop_list[2*i])*Isop_list[2*i+1]
  Uni_0.permute([1,2,5,6],2)
  Uni_list[2*i]=copy.copy(Uni_0)

 return   Uni_list


def Env_right_R(mps_R):
 E_right=[None]*mps_R.N
 for i in reversed(xrange(mps_R.N)):
  if i == 0:
   mps_R[i].setLabel([-1,-2,-3])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([-1,-2,-4])
   E_right[0]=mps_R_dag*mps_R[i]*E_right[1]
  elif i == (mps_R.N-1):
   mps_R[i].setLabel([-3,-2,1])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([-4,-2,1])
   E_right[i]=mps_R_dag*(mps_R[i])
   E_right[i].setLabel([-3,-4])
  else:
   mps_R[i].setLabel([1,-2,-3])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([2,-2,-4])
   E_right[i]=mps_R_dag*(E_right[i+1]*mps_R[i])
   E_right[i].permute([1,2],1)
   E_right[i].setLabel([-3,-4])
 return   E_right


def Env_left_R(mps_R):
 E_left=[None]*mps_R.N
 for i in xrange(mps_R.N):
  if i == 0:
   mps_R[i].setLabel([-1,-2,1])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([-1,-2,2])
   E_left[0]=mps_R_dag*mps_R[i]
   E_left[0].permute([1,2],1)
   E_left[0].setLabel([-3,-4])
  elif i == (mps_R.N-1):
   mps_R[i].setLabel([-3,-2,1])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([-4,-2,1])
   E_left[i]=mps_R_dag*(E_left[i-1]*mps_R[i])
  else:
   mps_R[i].setLabel([-3,-2,1])
   mps_R_dag=copy.copy(mps_R[i])
   mps_R_dag.setLabel([-4,-2,2])
   E_left[i]=mps_R_dag*(E_left[i-1]*mps_R[i])
   E_left[i].permute([1,2],1)
   E_left[i].setLabel([-3,-4])
 return E_left





def Update_Env_left(mps_R,E_left,location):
 if location==0:
  mps_R[location].setLabel([-1,-2,1])
  mps_R_dag=copy.copy(mps_R[location])
  mps_R_dag.setLabel([-1,-2,2])
  E_left[0]=mps_R_dag*mps_R[location]
  E_left[0].permute([1,2],1)
  E_left[0].setLabel([-3,-4])
 elif location==(mps_R.N-1):
  mps_R[location].setLabel([-3,-2,1])
  mps_R_dag=copy.copy(mps_R[location])
  mps_R_dag.setLabel([-4,-2,1])
  E_left[location-1].setLabel([-3,-4])
  E_left[location]=mps_R_dag*(E_left[location-1]*mps_R[location])
 else:
  mps_R[location].setLabel([-3,-2,1])
  mps_R_dag=copy.copy(mps_R[location])
  mps_R_dag.setLabel([-4,-2,2])
  E_left[location-1].setLabel([-3,-4])
  E_left[location]=mps_R_dag*(E_left[location-1]*mps_R[location])
  E_left[location].permute([1,2],1)
  E_left[location].setLabel([-3,-4])
 return E_left




def Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi):
 N=mps_R.N
 Env_right=[None]*(N/2)
 for i in xrange((N/2)-1, 0, -1):
  if i == ((N/2)-1):
   #print "i=initRight", i
   t0=uni10.UniTensor([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
   #print mps_R[2*i+1].printDiagram(), bdi
   t0.putBlock(mps_R[2*i+1].getBlock())
   t0.setLabel([1,11,4,3])
   t1=uni10.UniTensor([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
   t1.putBlock(mps_R[2*i].getBlock())
   t1.setLabel([7,10,6,1])

   T0=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T0.putBlock(mps_A[2*i+1].getBlock())
   T0.setLabel([18,13,4,3])
   T1=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T1.putBlock(mps_A[2*i].getBlock())
   T1.setLabel([17,14,6,18])

   U=Uni_list[2*i]
   U.setLabel([12,13,10,11])
   Up=Uni_list[2*i-1]
   Up.setLabel([15,14,16,12])
   Result=((((T0*t0)*U)*Up)*(t1))*T1
   Result.permute([7,16,15,17],4)
   Env_right[i]=Result
  else:
   t0=uni10.UniTensor([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
   t0.putBlock(mps_R[2*i+1].getBlock())
   t0.setLabel([1,11,4,-7])
   t1=uni10.UniTensor([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
   t1.putBlock(mps_R[2*i].getBlock())
   t1.setLabel([7,10,6,1])

   T0=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T0.putBlock(mps_A[2*i+1].getBlock())
   T0.setLabel([18,-15,4,-17])
   T1=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T1.putBlock(mps_A[2*i].getBlock())
   T1.setLabel([17,14,6,18])

   U=Uni_list[2*i]
   U.setLabel([12,-16,10,11])

   Up=Uni_list[2*i-1]
   Up.setLabel([15,14,16,12])

   Env_right[i+1].setLabel([-7,-16,-15,-17])
#   Result=((((t1)*(t0))*(U*Up))*(T1*T0))*(Env_right[i+1])
   Result=(((Env_right[i+1]*(T0*t0))*U)*Up)*(t1*T1)

   Result.permute([7,16,15,17],4)
   Env_right[i]=Result
 return   Env_right




def Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi):
 N=mps_A.N
 Env_left=[None]*(N/2)
 for i in xrange( 0, (N/2)-1):
   #print i, 2*i+1
   if i == 0:
      #print "i=initLeft", i
      t0=uni10.UniTensor([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
      t0.putBlock(mps_R[2*i].getBlock())
      t0.setLabel([1,11,4,3])
      t1=uni10.UniTensor([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
      t1.putBlock(mps_R[2*i+1].getBlock())
      t1.setLabel([3,10,6,7])
      #print mps_R[2*i+1].printDiagram()

      T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      T0.putBlock(mps_A[2*i].getBlock())
      T0.setLabel([1,12,4,18])
      T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      T1.putBlock(mps_A[2*i+1].getBlock())
      T1.setLabel([18,19,6,20])

      U=Uni_list[2*i]
      U.setLabel([12,13,11,10])

      Result=(((T0)*(t0))*(U))*(T1*t1)
      Result.permute([7,13,19,20],0)
      Env_left[i]=Result
   else:
      #print "i=insideLeft", i
      t0=uni10.UniTensor([mps_R[2*i].bond(0), bdi, bdi, mps_R[2*i].bond(2)])
      t0.putBlock(mps_R[2*i].getBlock())
      t0.setLabel([-7,11,4,3])
      t1=uni10.UniTensor([mps_R[2*i+1].bond(0), bdi, bdi, mps_R[2*i+1].bond(2)])
      t1.putBlock(mps_R[2*i+1].getBlock())
      t1.setLabel([3,10,6,7])
      #print mps_R[2*i+1].printDiagram()

      T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      T0.putBlock(mps_A[2*i].getBlock())
      T0.setLabel([-20,12,4,18])
      T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      T1.putBlock(mps_A[2*i+1].getBlock())
      T1.setLabel([18,19,6,20])

      U=Uni_list[2*i]
      U.setLabel([21,13,11,10])

      Up=Uni_list[2*i-1]
      Up.setLabel([-19,12,-13,21])

      Env_left[i-1].setLabel([-7,-13,-19,-20])

#      Result=((((t1)*(t0))*(U))*((T1*T0)*Up))*(Env_left[i-1])
      Result=(((Env_left[i-1]*Up)*(t0*T0))*U)*(t1*T1)

      Result.permute([7,13,19,20],4)
      Env_left[i]=Result

 return  Env_left

#@profile
def   Env_left_RU_update(mps_A, mps_R, Uni_list,i,bdip,bdi,Env_left):
 #print "\n"
 if  i==0:
  #print "i=initLeft", i
  t0=uni10.UniTensor([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
  t0.putBlock(mps_R[2*i].getBlock())
  t0.setLabel([1,11,4,3])
  t1=uni10.UniTensor([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
  t1.putBlock(mps_R[2*i+1].getBlock())
  t1.setLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()

  T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.putBlock(mps_A[2*i].getBlock())
  T0.setLabel([1,12,4,18])
  T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.putBlock(mps_A[2*i+1].getBlock())
  T1.setLabel([18,19,6,20])

  U=Uni_list[2*i]
  U.setLabel([12,13,11,10])

  #Result=(((t1)*(t0))*(U))*(T1*T0)
  Result=(((T0)*(t0))*(U))*(T1*t1)

  Result.permute([7,13,19,20],0)
  Env_left[i].putBlock(Result.getBlock())
 else:
  t0=uni10.UniTensor([mps_R[2*i].bond(0), bdi, bdi, mps_R[2*i].bond(2)])
  t0.putBlock(mps_R[2*i].getBlock())
  t0.setLabel([-7,11,4,3])
  t1=uni10.UniTensor([mps_R[2*i+1].bond(0), bdi, bdi, mps_R[2*i+1].bond(2)])
  t1.putBlock(mps_R[2*i+1].getBlock())
  t1.setLabel([3,10,6,7])

  T0=uni10.UniTensor([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.putBlock(mps_A[2*i].getBlock())
  T0.setLabel([-20,12,4,18])
  T1=uni10.UniTensor([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.putBlock(mps_A[2*i+1].getBlock())
  T1.setLabel([18,19,6,20])

  U=Uni_list[2*i]
  U.setLabel([21,13,11,10])

  Up=Uni_list[2*i-1]
  Up.setLabel([-19,12,-13,21])

  Env_left[i-1].setLabel([-7,-13,-19,-20])

  #Result=((((t1)*(t0))*(U))*((T1*T0)*Up))*(Env_left[i-1])
  Result=(((Env_left[i-1]*Up)*(t0*T0))*U)*(t1*T1)

  Result.permute([7,13,19,20],4)
  Env_left[i].putBlock(Result.getBlock())

 return  Env_left


####@profile
def Env_NR_f( mps_A, mps_R, Uni_list, bdip, bdi, position,E_right,E_left, Env_right, Env_left):

 N=mps_R.N
 Env_N=0
 if position == 0:
   E_right[1].setLabel([3,-3])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[0].bond(1).dim())
   t0=uni10.UniTensor([mps_R[0].bond(1), bdo])
   t0.identity()
   t0.setLabel([-2,2])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[0].bond(0).dim())
   t1=uni10.UniTensor([mps_R[0].bond(0), bdo])
   t1.identity()
   t1.setLabel([-1,1])
   Env_N=(t0*E_right[1])*(t1)
   Env_N.permute([-1,-2,-3,1,2,3],3)
 elif position == (mps_R.N-1):
   E_left[position-1].setLabel([1,-1])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[1].bond(1).dim())
   t0=uni10.UniTensor([mps_R[position].bond(1), bdo])
   t0.identity()
   t0.setLabel([-2,2])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[position].bond(2).dim())
   t1=uni10.UniTensor([mps_R[position].bond(2), bdo])
   t1.identity()
   t1.setLabel([-3,3])
   Env_N=(t0*E_left[position-1])*(t1)
   Env_N.permute([-1,-2,-3,1,2,3],3)
 else:
   E_right[position+1].setLabel([3,-3])
   E_left[position-1].setLabel([1,-1])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[1].bond(1).dim())
   t0=uni10.UniTensor([mps_R[1].bond(1), bdo])
   t0.identity()
   t0.setLabel([-2,2])
   Env_N=(t0*E_right[position+1])*(E_left[position-1])
   Env_N.permute([-1,-2,-3,1,2,3],3)





 Env_RU=0

 if position==0 or position==1:
  t0=uni10.UniTensor([mps_R[0].bond(0), bdi, bdi, mps_R[0].bond(2)])
  t0.putBlock(mps_R[0].getBlock())
  t0.setLabel([1,11,4,3])
  t1=uni10.UniTensor([mps_R[1].bond(0),bdi,bdi,mps_R[1].bond(2)])
  t1.putBlock(mps_R[1].getBlock())
  t1.setLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()
  T0=uni10.UniTensor([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
  T0.putBlock(mps_A[0].getBlock())
  T0.setLabel([1,12,4,18])
  T1=uni10.UniTensor([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
  T1.putBlock(mps_A[1].getBlock())
  T1.setLabel([18,19,6,20])

  U=Uni_list[0]
  U.setLabel([12,13,11,10])
  Env_right[1].setLabel([7,13,19,20])
  if position%2==0:
#    Env_RU=(((t1)*U)*(Env_right[1]))*(T1*T0)
    Env_RU=((Env_right[1]*(t1*T1))*U)*T0

    Env_RU.permute([1,11,4,3],3)
    Env_RU.permute([1,11,4,3],3)
    Env_RU.combineBond([11,4])
    Env_RU.permute([1,11,3],3)
  else:
#    Env_RU=((t0*U)*Env_right[1])*(T1*T0)
    Env_RU=(((t0*T0)*U)*T1)*Env_right[1]
    Env_RU.permute([3,10,6,7],3)
    Env_RU.combineBond([10,6])
    Env_RU.permute([3,10,7],3)

 elif position==(N-1) or position==(N-2):

  t0=uni10.UniTensor([mps_R[N-1].bond(0),bdi,bdi,mps_R[N-1].bond(2)])
  t0.putBlock(mps_R[N-1].getBlock())
  t0.setLabel([1,11,4,3])
  t1=uni10.UniTensor([mps_R[N-2].bond(0),bdi,bdi,mps_R[N-2].bond(2)])
  t1.putBlock(mps_R[N-2].getBlock())
  t1.setLabel([7,10,6,1])

  T0=uni10.UniTensor([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
  T0.putBlock(mps_A[N-1].getBlock())
  T0.setLabel([18,13,4,3])
  T1=uni10.UniTensor([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
  T1.putBlock(mps_A[N-2].getBlock())
  T1.setLabel([17,14,6,18])

  U=Uni_list[N-2]
  U.setLabel([12,13,10,11])

  Up=Uni_list[N-3]
  Up.setLabel([15,14,16,12])

  Env_left[((N-2)/2)-1].setLabel([7,16,15,17])
  if position==(N-1):
#   Env_RU=(((T1*T0)*(U*Up))*(Env_left[((N-2)/2)-1]))*(t1)
   Env_RU=(((((Env_left[((N-2)/2)-1])*t1)*T1)*Up)*U)*T0

   Env_RU.permute([1,11,4,3],3)
   Env_RU.combineBond([11,4])
   Env_RU.permute([1,11,3],3)

  else:
#   Env_RU=(((T0*T1)*(Up*U))*Env_left[((N-2)/2)-1])*(t0)
   Env_RU=((T0*t0)*U)*(Env_left[((N-2)/2)-1]*(Up*T1))
   Env_RU.permute([7,10,6,1],3)
   Env_RU.combineBond([10,6])
   Env_RU.permute([7,10,1],3)
 else:
  if (position%2)==0:
    #print "Location=IsoMidelodd", Location
    t0=uni10.UniTensor([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t0.putBlock(mps_R[position].getBlock())
    t0.setLabel([-7,11,4,3])
    t1=uni10.UniTensor([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t1.putBlock(mps_R[position+1].getBlock())
    t1.setLabel([3,10,6,7])
    #print mps_R[2*i+1].printDiagram()


    T0=uni10.UniTensor([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T0.putBlock(mps_A[position].getBlock())
    T0.setLabel([-20,12,4,18])
    T1=uni10.UniTensor([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T1.putBlock(mps_A[position+1].getBlock())
    T1.setLabel([18,19,6,20])

    U=Uni_list[position]
    U.setLabel([21,13,11,10])

    Up=Uni_list[position-1]
    Up.setLabel([-19,12,-13,21])

    #print "position", (position/2)-1
    Env_left[(position/2)-1].setLabel([-7,-13,-19,-20])
    Env_right[(position/2)+1].setLabel([7,13,19,20])
    #Env_RU=(((t1)*((U*Up))*(Env_right[(position/2)+1]))*(((T1*T0)*Env_left[(position/2)-1])))
    Env_RU=(((Env_right[(position/2)+1]*t1)*T1)*U)*((Env_left[(position/2)-1]*Up)*T0)

    Env_RU.permute([-7,11,4,3],3)
    Env_RU.combineBond([11,4])
    Env_RU.permute([-7,11,3],3)

  else:
    t0=uni10.UniTensor([mps_R[position-1].bond(0),bdi,bdi,mps_R[position-1].bond(2)])
    t0.putBlock(mps_R[position-1].getBlock())
    t0.setLabel([-7,11,4,3])
    t1=uni10.UniTensor([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t1.putBlock(mps_R[position].getBlock())
    t1.setLabel([3,10,6,7])

    T0=uni10.UniTensor([mps_A[position-1].bond(0),bdip,bdi,mps_A[position-1].bond(2)])
    T0.putBlock(mps_A[position-1].getBlock())
    T0.setLabel([-20,12,4,18])
    T1=uni10.UniTensor([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T1.putBlock(mps_A[position].getBlock())
    T1.setLabel([18,19,6,20])


    U=Uni_list[position-1]
    U.setLabel([21,13,11,10])

    Up=Uni_list[position-2]
    Up.setLabel([-19,12,-13,21])


    Env_left[((position-1)/2)-1].setLabel([-7,-13,-19,-20])
    Env_right[(position+1)/2].setLabel([7,13,19,20])



#    Env_RU=(((t0*(U*Up))*(Env_left[((position-1)/2)-1]))*((T1*T0)*Env_right[(position+1)/2]))
    Env_RU=((((Env_left[((position-1)/2)-1]*Up)*T0)*t0)*U)*(Env_right[(position+1)/2]*T1)

    Env_RU.permute([3,10,6,7],3)
    Env_RU.combineBond([10,6])
    Env_RU.permute([3,10,7],3)

 return   Env_N,  Env_RU





#@profile
def Env_U_f( mps_A, mps_R, Uni_list, bdip, bdi, position, Env_left, Env_right):


 N=mps_R.N

 Env_RU=0

 if position==0:
  t0=uni10.UniTensor([mps_R[0].bond(0), bdi, bdi, mps_R[0].bond(2)])
  t0.putBlock(mps_R[0].getBlock())
  t0.setLabel([1,11,4,3])
  t1=uni10.UniTensor([mps_R[1].bond(0),bdi,bdi,mps_R[1].bond(2)])
  t1.putBlock(mps_R[1].getBlock())
  t1.setLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()
  T0=uni10.UniTensor([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
  T0.putBlock(mps_A[0].getBlock())
  T0.setLabel([1,12,4,18])
  T1=uni10.UniTensor([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
  T1.putBlock(mps_A[1].getBlock())
  T1.setLabel([18,19,6,20])

  U=Uni_list[0]
  U.setLabel([12,13,11,10])
  
  Env_right[1].setLabel([7,13,19,20])
#  Env_RU=(((t1*t0))*(Env_right[1]))*(T1*T0)
  Env_RU=((Env_right[1]*(t1*T1))*t0)*T0

  Env_RU.permute([12,13,11,10],2)

 elif position==(N-2) or position==(N-3):

  t0=uni10.UniTensor([mps_R[N-1].bond(0),bdi,bdi,mps_R[N-1].bond(2)])
  t0.putBlock(mps_R[N-1].getBlock())
  t0.setLabel([1,11,4,3])
  t1=uni10.UniTensor([mps_R[N-2].bond(0),bdi,bdi,mps_R[N-2].bond(2)])
  t1.putBlock(mps_R[N-2].getBlock())
  t1.setLabel([7,10,6,1])

  T0=uni10.UniTensor([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
  T0.putBlock(mps_A[N-1].getBlock())
  T0.setLabel([18,13,4,3])
  T1=uni10.UniTensor([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
  T1.putBlock(mps_A[N-2].getBlock())
  T1.setLabel([17,14,6,18])

  U=Uni_list[N-2]
  U.setLabel([12,13,10,11])

  Up=Uni_list[N-3]
  Up.setLabel([15,14,16,12])

  Env_left[((N-2)/2)-1].setLabel([7,16,15,17])
  if position==(N-2):
#   Env_RU=((T0*(T1*Up))*(Env_left[((N-2)/2)-1]))*(t1*t0)
   Env_RU=(((((Env_left[((N-2)/2)-1])*Up)*T1)*t1)*(T0*t0))

   Env_RU.permute([12,13,10,11],2)
   #Env_RU.combineBond([11,4])
   #Env_RU.permute([1,11,3],2)
  else:
#   Env_RU=(((t0)*(t1*U))*Env_left[((N-2)/2)-1])*(T1*T0)
   Env_RU=((T0*t0)*U)*(Env_left[((N-2)/2)-1]*(t1*T1))

   Env_RU.permute([15,14,16,12],2)
   #Env_RU.combineBond([10,6])
   #Env_RU.permute([7,10,1],2)
 else:
  if (position%2)==0:
    #print "Location=IsoMidelodd", Location
    t0=uni10.UniTensor([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t0.putBlock(mps_R[position].getBlock())
    t0.setLabel([-7,11,4,3])
    t1=uni10.UniTensor([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t1.putBlock(mps_R[position+1].getBlock())
    t1.setLabel([3,10,6,7])
    #print mps_R[2*i+1].printDiagram()


    T0=uni10.UniTensor([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T0.putBlock(mps_A[position].getBlock())
    T0.setLabel([-20,12,4,18])
    T1=uni10.UniTensor([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T1.putBlock(mps_A[position+1].getBlock())
    T1.setLabel([18,19,6,20])

    U=Uni_list[position]
    U.setLabel([21,13,11,10])

    Up=Uni_list[position-1]
    Up.setLabel([-19,12,-13,21])

    #print (Location/2)-1, (Location/2)+1
    Env_left[(position/2)-1].setLabel([-7,-13,-19,-20])
    Env_right[(position/2)+1].setLabel([7,13,19,20])



#    Env_RU=(((T1*(T0*Up))*(Env_left[(position/2)-1]))*(((t1*t0)*Env_right[(position/2)+1])))
    Env_RU=(((Env_right[(position/2)+1]*t1)*T1))*(((Env_left[(position/2)-1]*Up)*T0)*t0)

    Env_RU.permute([21,13,11,10],2)
    #Env_RU.combineBond([11,4])
    #Env_RU.permute([-7,11,3],2)
  else:
    t0=uni10.UniTensor([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t0.putBlock(mps_R[position+1].getBlock())
    t0.setLabel([-7,11,4,3])
    t1=uni10.UniTensor([mps_R[position+2].bond(0),bdi,bdi,mps_R[position+2].bond(2)])
    t1.putBlock(mps_R[position+2].getBlock())
    t1.setLabel([3,10,6,7])

    T0=uni10.UniTensor([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T0.putBlock(mps_A[position+1].getBlock())
    T0.setLabel([-20,12,4,18])
    T1=uni10.UniTensor([mps_A[position+2].bond(0),bdip,bdi,mps_A[position+2].bond(2)])
    T1.putBlock(mps_A[position+2].getBlock())
    T1.setLabel([18,19,6,20])


    U=Uni_list[position+1]
    U.setLabel([21,13,11,10])

    Up=Uni_list[position]
    Up.setLabel([-19,12,-13,21])

    Env_left[((position-1)/2)].setLabel([-7,-13,-19,-20])
    Env_right[((position+1)/2)+1].setLabel([7,13,19,20])

#    Env_RU=(((t0*(t1*U))*(Env_right[((position+1)/2)+1]))*((T1*T0)*Env_left[((position-1)/2)]))
    Env_RU=Env_left[((position-1)/2)]*(((((Env_right[((position+1)/2)+1]*T1)*t1)*U)*t0)*T0)

    Env_RU.permute([-19,12,-13,21],2)
    #Env_RU.combineBond([10,6])
    #Env_RU.permute([3,10,7],2)

 return   Env_RU
####@profile  
def solve_linear_eq(A,Ap):
 #print A
 Ap_h=copy.copy(Ap)
 Ap_h.transpose()
 Result=uni10.UniTensor(Ap.bond())
 blk_qnums = A.blockQnum()
 #print blk_qnums, Ap.printDiagram()
 blk_qnums1 = Ap.blockQnum()
 for qnum in blk_qnums:
   if qnum in blk_qnums1:
    A_mat=A.getBlock(qnum)
    Ap_mat=Ap.getBlock(qnum)
    A_np=Mat_uni_to_np(A_mat)
    b_np=Mat_uni_to_np(Ap_mat)
    #x_np=np.linalg.lstsq(A_np, b_np, rcond=-1)[0] 
    x_np=sp.linalg.lstsq(A_np, b_np,cond=1.e-12,overwrite_a=False,overwrite_b=False,check_finite=True, lapack_driver='gelsy')[0] 
    #x_np=sp.linalg.lstsq(A_np, b_np,cond=1.e-12,overwrite_a=False,overwrite_b=False,check_finite=True)[0] 

    x=Mat_np_to_Uni(x_np)
    Result.putBlock(qnum, x)
 return Result

#@profile
def  Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list):
 Env_R=0
 Result3=mps_A.norm()

 for q in xrange(N_iter):
  t0=time.time()
  #print "Step", q
  E_right=Env_right_R(mps_R)
  E_left=Env_left_R(mps_R)
  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)

  if Opt_r=="R":
   for i in xrange(mps_R.N): 
    if i%2==0 and i>0:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list, (i/2)-1, bdip, bdi, Env_left)

    Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i, E_right, E_left, Env_right, Env_left)
    Env_RU.setLabel([1,2,3])
    Env_N.setLabel([-1,-2,-3,1,2,3])
    mps_R[i].setLabel([1,2,3])
    mps_R_dagg=copy.copy(mps_R[i])
    mps_R_dagg.setLabel([-1,-2,-3])
    Result1=mps_R_dagg*Env_N*mps_R[i]
    Result2=Env_RU*mps_R[i]
    #print "R", i, abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5))#,Result1[0],Result2[0], Result3
    #print Env_N.printDiagram()
    count_list.append(count_list[-1]+1)
    Fidel_list.append(abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5)))
    if Opt[2]=="SVD":
     U, S, V=svd_parity2(Env_N)
     #print "Hiiiiiii"
     U.transpose()
     V.transpose()
     S=inverse(S)
     U.setLabel([5,6,7,8])
     S.setLabel([4,5])
     V.setLabel([1,2,3,4])
     Env_N_inv=(V*S)*U
     Env_N_inv.permute([1,2,3,6,7,8],3)
     Env_N_inv.setLabel([-1,-2,-3,1,2,3])
     R_new=Env_N_inv*Env_RU
     R_new.permute([-1,-2,-3],2)
    elif  Opt[2]=="cg":
     try:
      R_new=solve_linear_eq(Env_N,Env_RU)
      R_new.setLabel([0,1,2])
      R_new.permute([0,1,2],2)   
     except:
      print "error, SVD"
      U, S, V=svd_parity2(Env_N)
      #print "Hiiiiiii"
      U.transpose()
      V.transpose()
      S=inverse(S)
      U.setLabel([5,6,7,8])
      S.setLabel([4,5])
      V.setLabel([1,2,3,4])
      Env_N_inv=(V*S)*U
      Env_N_inv.permute([1,2,3,6,7,8],3)
      Env_N_inv.setLabel([-1,-2,-3,1,2,3])
      R_new=Env_N_inv*Env_RU
      R_new.permute([-1,-2,-3],2)
      
      
    mps_R[i].putBlock(R_new.getBlock())
    E_left=Update_Env_left(mps_R, E_left,i)
    #Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)

  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  if Opt_u=="U":
   for i in xrange(mps_R.N-1): 
    if (i%2)==1:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list,(i-1)/2,bdip,bdi,Env_left)
    Env_U=Env_U_f(mps_A, mps_R, Uni_list, bdip, bdi, i, Env_left, Env_right)

    Env_U.setLabel([1,2,3,4])
    Uni_list[i].setLabel([1,2,3,4])
    Result1=mps_R.norm()
    Result2=Env_U*Uni_list[i]
    #print "U", i, abs(Result2[0])/(abs(Result1*Result3)**(0.5))#,Result1, Result2[0]
    count_list.append(count_list[-1]+1)
    Fidel_list.append(abs(Result2[0])/(abs(Result1*Result3)**(0.5)))
    svd=Env_U.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    Uni_list[i].putBlock(temporary_matrix)
  #print time.time() - t0, "time_needed:Update"


 Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 E_right=Env_right_R(mps_R)
 E_left=Env_left_R(mps_R)
 Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i,E_right,E_left,Env_right,Env_left)
 Env_RU.setLabel([1,2,3])
 Env_N.setLabel([-1,-2,-3,1,2,3])
 mps_R[i].setLabel([1,2,3])
 mps_R_dagg=copy.copy(mps_R[i])
 mps_R_dagg.setLabel([-1,-2,-3])
 Result1=mps_R_dagg*Env_N*mps_R[i]
 Result2=Env_RU*mps_R[i]
 Fidel_final=abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5))
 return  mps_R, Uni_list, Fidel_final


###@profile
def update_Rnew(Env_N, mps_R_origin,i,bdi,Env_RU,Result3):

 D=bdi.dim()
 #bdi=uni10.Bond(uni10.BD_IN, mps_R_origin[0])
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiN=uni10.Bond(uni10.BD_IN, D)
 bdoN=uni10.Bond(uni10.BD_OUT, D)


 bdiMPS=uni10.Bond(uni10.BD_IN, mps_R_origin[2*i+1].bond(2).dim())
 bdoMPS=uni10.Bond(uni10.BD_OUT, mps_R_origin[2*i].bond(0).dim())
 
 N_open=uni10.UniTensor([mps_R_origin[2*i].bond(0), bdiN,bdiN,bdiN,bdiN, bdiMPS,bdoMPS, bdoN,bdoN,bdoN,bdoN, mps_R_origin[2*i+1].bond(2)], "Ten_R")
 #print Env_N.printDiagram(), N_open.printDiagram()
 N_open.putBlock(Env_N.getBlock())
 N_open.setLabel([-1, -2, -5, -3, -6, -7, 1, 2, 5, 3, 6,7])


 Nu_open=uni10.UniTensor([mps_R_origin[2*i].bond(0), bdiN, bdiN,bdiN, bdiN, bdiMPS], "Ten_R")
 Nu_open.putBlock(Env_RU.getBlock())
 Nu_open.setLabel([1, 2, 5, 3, 6,7])


 Ten_1=uni10.UniTensor([mps_R_origin[2*i].bond(0), bdi, bdi, mps_R_origin[2*i].bond(2)], "Ten_1")
 Ten_1.setLabel([1,2,3,4])
 Ten_1.putBlock(mps_R_origin[2*i].getBlock())
 Ten_1_dagg=copy.copy(Ten_1)
 Ten_1_dagg.setLabel([-1,-2,-3,-4])



 Ten_2=uni10.UniTensor([mps_R_origin[2*i+1].bond(0), bdi, bdi, mps_R_origin[2*i+1].bond(2)], "Ten_2")
 Ten_2.putBlock(mps_R_origin[2*i+1].getBlock())
 Ten_2.setLabel([4,5,6,7])
 Ten_2_dagg=copy.copy(Ten_2)
 Ten_2_dagg.setLabel([-4,-5,-6,-7])


 Ten_1_f=copy.copy(Ten_1)
 Ten_2_f=copy.copy(Ten_2)
 Ten_2_f_dagg=copy.copy(Ten_2)
 Ten_1_f_dagg=copy.copy(Ten_1)
 
 norm1=Ten_1_dagg*((Ten_2_dagg*(N_open*Ten_2))*Ten_1)
 norm2=(Nu_open*Ten_2)*Ten_1
 E0=abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))
 E1=10
 for q in xrange(1):

  N_open1=Ten_2_dagg*(N_open*Ten_2)
  N_open1.permute([-1,-2,-3,-4,1,2,3,4],4)
  Nu_open1=(Nu_open*Ten_2)
  Nu_open1.permute([1,2,3,4],4)
  U, S, V=svd_parity3(N_open1)
  U.transpose()
  V.transpose()
  S=inverse(S)
  U.setLabel([6,7,8,9,10])
  S.setLabel([5,6])
  V.setLabel([1,2,3,4,5])
  N_open1_inv=(V*S)*U
  N_open1_inv.permute([1,2,3,4,7,8,9,10],4)
  N_open1_inv.setLabel([-1,-2,-3,-4,1,2,3,4])
  Ten_1=N_open1_inv*Nu_open1
  Ten_1.permute([-1,-2,-3,-4],3)
  Ten_1.setLabel([1,2,3,4])
  Ten_1_dagg=copy.copy(Ten_1)
  Ten_1_dagg.setLabel([-1,-2,-3,-4])

  norm1=Ten_1_dagg*((Ten_2_dagg*(N_open*Ten_2))*Ten_1)
  norm2=(Nu_open*Ten_2)*Ten_1
  
  #print "1", q, norm1[0], norm2[0],abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))
  E1=abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))


  if E1>E0 and E1<1.0:
   E0=E1*1.0
   Ten_1_f=copy.copy(Ten_1)
   Ten_1_f_dagg=copy.copy(Ten_1)
   Ten_1_f_dagg.setLabel([-1,-2,-3,-4])
  else:
   Ten_1=copy.copy(Ten_1_f)
   Ten_1_dagg=copy.copy(Ten_1_f)
   Ten_1_dagg.setLabel([-1,-2,-3,-4])




  N_open1=Ten_1_dagg*(N_open*Ten_1)
  N_open1.permute([-4,-5,-6,-7,4,5,6,7],4)
  Nu_open1=(Nu_open*Ten_1)
  Nu_open1.permute([4,5,6,7],4)
  U, S, V=svd_parity3(N_open1)
  U.transpose()
  V.transpose()
  S=inverse(S)
  U.setLabel([6,7,8,9,10])
  S.setLabel([5,6])
  V.setLabel([1,2,3,4,5])
  N_open1_inv=(V*S)*U
  N_open1_inv.permute([1,2,3,4,7,8,9,10],4)
  N_open1_inv.setLabel([-4,-5,-6,-7,4,5,6,7])
  Ten_2=N_open1_inv*Nu_open1
  Ten_2.permute([-4,-5,-6,-7],4)
  Ten_2.setLabel([4,5,6,7])
  Ten_2_dagg=copy.copy(Ten_2)
  Ten_2_dagg.setLabel([-4,-5,-6,-7])

  norm1=Ten_1_dagg*(Ten_2_dagg*(N_open*Ten_2)*Ten_1)
  norm2=(Nu_open*Ten_2)*Ten_1
  #print "2", q, norm1[0], norm2[0],abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))
  E1=abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))

  if E1>E0 and E1<1.0: 
   E0=E1*1.0
   Ten_2_f=copy.copy(Ten_2)
   Ten_2_f_dagg=copy.copy(Ten_2)
   Ten_2_f_dagg.setLabel([-4,-5,-6,-7])
  else:
   Ten_2=copy.copy(Ten_2_f)
   Ten_2_dagg=copy.copy(Ten_2_f)
   Ten_2_dagg.setLabel([-4,-5,-6,-7])


# norm1=Ten_1_f_dagg*(Ten_2_f_dagg*(N_open*Ten_2_f)*Ten_1_f)
# norm2=(Nu_open*Ten_2_f)*Ten_1_f
# print "final", q, norm1[0], norm2[0],abs(norm2[0])/(abs(norm1[0]*Result3)**(0.5))


 Ten_1_f.setLabel([1,2,3,4])
 Ten_2_f.setLabel([4,5,6,7])
 mps_R_new=Ten_1_f*Ten_2_f
 mps_R_new.permute([1,2,5,3,6,7],3)
 mps_R_new.combineBond([2,5])
 mps_R_new.combineBond([3,6])
 mps_R_new.combineBond([2,3])

 mps_R_new.permute([1,2,7],2)

 Ten_1_f.combineBond([2,3])
 Ten_2_f.combineBond([5,6])
 Ten_1_f.permute([1,2,4],2)
 Ten_2_f.permute([4,5,7],2)

 mps_R_origin[2*i].putBlock(Ten_1_f.getBlock())
 mps_R_origin[2*i+1].putBlock(Ten_2_f.getBlock())

 return   mps_R_new






###@profile
def  Update_RU_eff(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list, mps_R_origin,bdi_origin):
 Env_R=0
 Result3=mps_A.norm()
 for q in xrange(N_iter):
  t0=time.time()

  #print "Step", q
  E_right=Env_right_R(mps_R)
  E_left=Env_left_R(mps_R)
  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  E_val_0=0
  E_val_1=1
  if Opt_r=="R":
   for i in xrange(mps_R.N): 
    if i%2==0 and i>0:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list, (i/2)-1, bdip, bdi, Env_left)

    Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i, E_right, E_left, Env_right, Env_left)
    Env_RU.setLabel([1,2,3])
    Env_N.setLabel([-1,-2,-3,1,2,3])
    mps_R[i].setLabel([1,2,3])
    mps_R_dagg=copy.copy(mps_R[i])
    mps_R_dagg.setLabel([-1,-2,-3])
    Result1=mps_R_dagg*Env_N*mps_R[i]
    Result2=Env_RU*mps_R[i]
    E_val_1=abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5))
    #print "R", i, E_val_1#,Result1[0],Result2[0]#, Result3
    #print Env_N.printDiagram()
    count_list.append(count_list[-1]+1)
    Fidel_list.append(abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5)))
    if Opt[2]=="SVD":
       R_new=update_Rnew(Env_N, mps_R_origin, i, bdi_origin, Env_RU,Result3)
       R_new.setLabel([-1,-2,-3])
#      U, S, V=svd_parity2(Env_N)
#      #print "Hiiiiiii"
#      U.transpose()
#      V.transpose()
#      S=inverse(S)
#      U.setLabel([5,6,7,8])
#      S.setLabel([4,5])
#      V.setLabel([1,2,3,4])
#      Env_N_inv=(V*S)*U
#      Env_N_inv.permute([1,2,3,6,7,8],3)
#      Env_N_inv.setLabel([-1,-2,-3,1,2,3])
#      R_new=Env_N_inv*Env_RU
#      R_new.permute([-1,-2,-3],2)
    if E_val_1>E_val_0 and E_val_1<1.0:
     E_val_0=E_val_1*1.0
     mps_R[i].putBlock(R_new.getBlock())
    else: 1+1#print "breakR",E_val_1,E_val_0 
    E_left=Update_Env_left(mps_R, E_left,i)
    #Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)

  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  if Opt_u=="U":
   for i in xrange(mps_R.N-1): 
    if (i%2)==1:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list,(i-1)/2,bdip,bdi,Env_left)
    Env_U=Env_U_f(mps_A, mps_R, Uni_list, bdip, bdi, i, Env_left, Env_right)

    Env_U.setLabel([1,2,3,4])
    Uni_list[i].setLabel([1,2,3,4])
    Result1=mps_R.norm()
    Result2=Env_U*Uni_list[i]
    #print "U", i, abs(Result2[0])/(abs(Result1*Result3)**(0.5))#,Result1, Result2[0]
    count_list.append(count_list[-1]+1)
    Fidel_list.append(abs(Result2[0])/(abs(Result1*Result3)**(0.5)))
    svd=Env_U.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    Uni_list[i].putBlock(temporary_matrix)
  #print time.time() - t0, "time_needed:Update"


 Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 E_right=Env_right_R(mps_R)
 E_left=Env_left_R(mps_R)
 Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i,E_right,E_left,Env_right,Env_left)
 Env_RU.setLabel([1,2,3])
 Env_N.setLabel([-1,-2,-3,1,2,3])
 mps_R[i].setLabel([1,2,3])
 mps_R_dagg=copy.copy(mps_R[i])
 mps_R_dagg.setLabel([-1,-2,-3])
 Result1=mps_R_dagg*Env_N*mps_R[i]
 Result2=Env_RU*mps_R[i]
 Fidel_final=abs(Result2[0])/(abs(Result1[0]*Result3)**(0.5))

 return  mps_R_origin, Uni_list, Fidel_final











def   Fidel_basedonQR(MPS_A, MPS_R, MPS_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop):

 chi=200
 if MPS_Q.D>chi:
   MPS_Q=MPS_Q.appSVD(chi)

 chi=200
 if MPS_R.D>chi:
   MPS_R=MPS_R.appSVD(chi)

 A_list=[]
 for i in xrange(N_y):
   t0=uni10.UniTensor([MPS_Q[i].bond(0),bdip,bdi,MPS_Q[i].bond(2)])
   t0.putBlock(MPS_Q[i].getBlock())
   t0.setLabel([-1,2,3,4])
   t1=uni10.UniTensor([MPS_R[i].bond(0),bdi,bdi,MPS_R[i].bond(2)])
   t1.putBlock(MPS_R[i].getBlock())
   t1.setLabel([1,3,5,-4])
   t_result=t0*t1
   t_result.permute([1, -1, 2, 5, -4, 4],4)
   t_result.combineBond([1,-1])
   t_result.combineBond([-4,4])
   t_result.combineBond([2,5])
   t_result.permute([1, 2, -4],2)
   A_list.append(t_result)

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(A_list[q].bond(2).dim())
 #print max(list_bond)
 mps_QR=MPSclass.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y)
 #randuni, rand, ortho
 for i in xrange(N_y):
   mps_QR[i]=copy.copy(A_list[i])
 #print mps_QR[0].printDiagram(),mps_QR.norm() 

 #print MPS_R.norm(), MPS_Q.norm(), mps_QR.norm(), mps_QR.product(MPS_A), mps_QR.fidel(MPS_A)

 return    mps_QR.fidel(MPS_A)



def  Root_update(mps_A, N, bdi, bdip, bdo, bdop, N_iter, Opt,Fidel_list, count_list):


 Uni_list=[None]*(N-1)
 for i in xrange(mps_A.N-1):
  if i%2==0:
    t1=uni10.UniTensor([bdip,bdip,bdo,bdo])
    t1.identity()
    #t1.randomize()
    #svd=t1.getBlock().svd()
    #t1.putBlock(svd[0])
    Uni_list[i]=t1
  else:
    t1=uni10.UniTensor([bdip,bdip,bdop,bdop])
    t1.identity()
    Uni_list[i]=t1
 
 chi=Opt[3]
 mps_R=root.make_mps_R_root(mps_A, N, bdi, bdip, bdo, bdop,chi)
 mps_R_t=mps_R.appSVD(bdi.dim())
 Uni_list=Simple_update_root(mps_A, mps_R_t, N, bdi, bdip, bdo, bdop,Uni_list)
 Opt_r="R"
 Opt_u="U"
 mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter, Opt, Fidel_list, count_list)

 return mps_R, Uni_list

def uni_shrink(Uni_list):
#  Uni_list_copy=[ Uni_list[i]  for i in xrange(len(Uni_list))]
#  U_ten=copy.copy(Uni_list[len(Uni_list)-1])
#  U_ten.identity()
#  Uni_list_copy.append(U_ten)
 #print len(Uni_list)-1
 U_l=[]
 for i in xrange(0,len(Uni_list),4):
  #print "StepsEven", i

  Uni_list[i].setLabel([1,2,3,4])
  Uni_list[i+1].setLabel([5,6,2,7])
  Uni_list[i+2].setLabel([7,8,9,10])

  Result=(Uni_list[i]*Uni_list[i+1])*Uni_list[i+2]
  Result.permute([1,5,6,8,3,4,9,10],8)
  Result.combineBond([1,5])
  Result.combineBond([6,8])
  Result.combineBond([3,4])
  Result.combineBond([9,10])
  Result.permute([1,6,3,9],2)
  U_l.append(Result)
  if  (i+3)<= (len(Uni_list)-1):
   t1=uni10.UniTensor([Uni_list[i+3].bond(0),Uni_list[i+3].bond(3)])
   t1.identity()
   t1.setLabel([5,6])
   t2=uni10.UniTensor([Uni_list[i+3].bond(0),Uni_list[i+3].bond(3)])
   t2.identity()
   t2.setLabel([7,8])
   Uni_list[i+3].setLabel([1,2,3,4])
   Result=(Uni_list[i+3]*t1)*t2
   Result.combineBond([5,1])
   Result.combineBond([2,7])
   Result.combineBond([6,3])
   Result.combineBond([4,8])
   Result.permute([5,2,6,4],2)
   U_l.append(Result)
 return  U_l


def   Sqrt(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum)
   Landam=Landa_cp.getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     if Landam[i*D+j] > 1.0e-12:
      Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
     else:
      Landa_cpm[i*D+j]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 
##@profile
def  back_to_normal(mps_R, bdip, bdi, bdop, bdo):

#   Ten_1=uni10.UniTensor([mps[i].bond(0), bdip, bdi, mps[i].bond(2)], "Ten_R")
#   #print mps[i].printDiagram()
#   Ten_1.setLabel([1,2,3,4])
#   Ten_1.putBlock(mps[i].getBlock())
#   Ten_2=uni10.UniTensor([mps[i+1].bond(0), bdip, bdi, mps[i+1].bond(2)], "Ten_R")
#   Ten_2.putBlock(mps[i+1].getBlock())
#   Ten_2.setLabel([4,-2,-3,-4])
#   Results=Ten_1*Ten_2  
#   Results.permute([1,2,-2,3,-3,-4], 5)
#   Results.combineBond([2,-2])
#   Results.combineBond([3,-3])
#   Results.combineBond([2,3])
#   Results.permute([1,2,-4],2)

 A_list=[]
 for i in xrange(mps_R.N):
  #print "i", i
  t0=uni10.UniTensor([mps_R[i].bond(0),bdip,bdip,bdi,bdi,mps_R[i].bond(2)])
  t0.putBlock(mps_R[i].getBlock())
  t0.setLabel([1,2,-2,3,-3,-4])
  t0.permute([1,2,3,-2,-3,-4],3)
  U, S, V=svd_parity2(t0)
  U.setLabel([1,2,3,10])  
  S=Sqrt(S)
  S.setLabel([10, 4])
  U=U*S
  U.combineBond([2,3])
  U.setLabel([1,2,4])  

  S.setLabel([4, -10])
  V.setLabel([-10,-2,-3,-4])  
  V.permute([-10,-2, -3, -4],3)
  V=S*V
  V.permute([4,-2,-3,-4],3)
  V.combineBond([-2,-3])
  V.setLabel([4,-2,-4])  

  A_list.append(U)
  A_list.append(V)

 list_bond=[]
 for q in xrange(len(A_list)):
   list_bond.append(A_list[q].bond(2).dim())
 mps_R=MPSclass.MPS( A_list[1].bond(1).dim(), max(list_bond), len(A_list))
 for i in xrange(len(A_list)):
  mps_R[i]=copy.copy(A_list[i])
 return  mps_R


def Adding_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list):

 #print "A", mps_A[1].printDiagram(), mps_R[1].printDiagram(), bdi.dim(), bdip.dim(),Uni_list[1].printDiagram(),Uni_list[0].printDiagram() 

 #print "norm_A", mps_A.norm(), mps_R.norm(), len(Uni_list)
 mps_A_update=shrink_mps(mps_A, bdip.dim(), bdi.dim())
 mps_R_update=shrink_mps(mps_R, bdi.dim(), bdi.dim())
 Uni_list_update=uni_shrink(Uni_list)


 #print "norm_A_update", mps_A_update.norm(),mps_R_update.norm(), len(Uni_list_update)
 #print mps_A_update[0].printDiagram(),mps_R_update[0].printDiagram() 
 #print len(Uni_list_update), Uni_list_update[0].printDiagram(), Uni_list_update[1].printDiagram()

 D=bdi.dim()
 Dp=bdip.dim()

 bdi=uni10.Bond(uni10.BD_IN, D*D)
 bdo=uni10.Bond(uni10.BD_OUT, D*D)
 

 bdip=uni10.Bond(uni10.BD_IN, Dp*Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp*Dp)

 #print mps_A_update[1].printDiagram(), mps_R_update[1].printDiagram(), D*D, Dp*Dp,Uni_list_update[0].printDiagram() , Uni_list_update[1].printDiagram()
 #print "Hiiii"
 Opt_r="R"
 Opt_u="U"
 N_iter1=N_iter
 mps_R, Uni_list_new, Fidel_RU=Update_RU(mps_A_update, mps_R_update, Uni_list_update, bdip, bdi, Opt_r, Opt_u, N_iter1,Opt,Fidel_list, count_list)
 
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)


 mps_R=back_to_normal(mps_R, bdi, bdi, bdo, bdo)
 #mps_Q=make_MPS_Q(Uni_list_new)
 #mps_Q=back_to_normal(mps_Q, bdip, bdi, bdop, bdo)
 mps_Q="None"
 return  mps_R, mps_Q, Uni_list_new, Fidel_RU



###@profile

def Adding_updateU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list):

 mps_A_update=shrink_mps(mps_A, bdip.dim(), bdi.dim())
 mps_R_update=shrink_mps(mps_R, bdi.dim(), bdi.dim())
 Uni_list_update=uni_shrink(Uni_list)

 bdi_origin=copy.copy(bdi)
 D=bdi.dim()
 Dp=bdip.dim()

 bdi=uni10.Bond(uni10.BD_IN, D*D)
 bdo=uni10.Bond(uni10.BD_OUT, D*D)

 bdip=uni10.Bond(uni10.BD_IN, Dp*Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp*Dp)

 Opt_r="R"
 Opt_u="U"
 N_iter1=N_iter
 mps_R, Uni_list_new, Fidel_RU=Update_RU_eff(mps_A_update, mps_R_update, Uni_list_update, bdip, bdi, Opt_r, Opt_u, N_iter1,Opt,Fidel_list, count_list,mps_R, bdi_origin)
 
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 mps_Q="None"
 #mps_R=back_to_normal(mps_R, bdi, bdi, bdo, bdo)
 #print "hi0"
 #mps_Q=make_MPS_Q(Uni_list_new)
 #print "hi1"
 #mps_Q=back_to_normal(mps_Q, bdip, bdi, bdop, bdo)
 #print "hi2"


 return  mps_R, mps_Q, Uni_list_new, Fidel_RU






















####@profile
def QR_canon(mps_A_origin, N_iter, Dp, D,Opt):

  
 N_y=mps_A_origin.N               #Even
 #chi_canon            #less than DDdd
 #N_iter=20
 #Norm_init_A=1.0
##############################################################

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)

 count_list=[]
 Fidel_list=[]
 count_list.append(0)
################PEPS-tensor##############################
 if Dp<D:
  mps_Q="None"
  Uni_list="None"
  mps_A=copy.copy(mps_A_origin)
  Iso_list, Uni_list=Simple_trivial(mps_A, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list)
  mps_R=make_MPS_R(mps_A, Iso_list,N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Uni_list=Uni_update(Iso_list, Uni_list)
  mps_Q=make_MPS_Q(Uni_list)
  Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
  print "FidelQR_Dp<D_QR=" , abs(Fidel_final)
  return  mps_R,mps_Q, Uni_list, abs(Fidel_final)


 #print "norm-mps_A", mps_A.norm() 
 #mps_A=Init_A_PEPS(N,bdi,bdi1,bdip,bdo,bdo1,bdop)
 #mps_A=mps_A.normalize()
 #mps_A=mps_A*(Norm_init_A)
 #Norm_val=mps_A_origin.norm()
 #mps_A=mps_A_origin*(1/(Norm_val**0.5))

 Fidel_RU=1.0
 if Opt[0]=="root":
   mps_A=copy.copy(mps_A_origin)
   mps_R, Uni_list=Root_update(mps_A,N_y, bdi, bdip, bdo, bdop, N_iter, Opt,Fidel_list, count_list)
   mps_Q=make_MPS_Q(Uni_list)
   Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, bdi, bdo)
   print "FidelQR_root" , abs(Fidel_final)
   return  mps_R, mps_Q, Uni_list, abs(Fidel_final)
 else:
   mps_A=copy.copy(mps_A_origin)
   Iso_list, Uni_list=Simple_update(mps_A, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list)
   Isop_list=[ copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]
   #Iso_list, Isop_list, Uni_list=Full_update(mps_A,Iso_list,Uni_list,N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, 4, Fidel_list, count_list)
   mps_R=make_MPS_R(mps_A, Iso_list, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   Uni_list=Uni_update(Isop_list, Uni_list)
   mps_Q=0
   if Opt[1]=="Single":
    Opt_r="R"
    Opt_u="U"
    mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt, Fidel_list, count_list)
    mps_Q=make_MPS_Q(Uni_list)
    print "Fidel=QR-A", abs(Fidel_RU)
    return  mps_R, mps_Q, Uni_list, abs(Fidel_RU)
   #Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #print "Fidel=QR-A" , abs(Fidel_final), Fidel_RU
   elif Opt[1]=="contiguous":
    #print N_iter
    Opt_r="R"
    Opt_u="U"
    mps_R, mps_Q, Uni_list, Fidel_RU=Adding_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list)
    print "Fidel=QR-A", abs(Fidel_RU)
    return  mps_R, mps_Q, Uni_list, abs(Fidel_RU)

   #Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #print "Fidel=QR-A" , abs(Fidel_final)
   elif Opt[1]=="ACCcontiguous":
    Opt_r="R"
    Opt_u="U"
    mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt, Fidel_list, count_list)
    mps_R, mps_Q, Uni_list, Fidel_RU=Adding_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list)
    print "Fidel=QR-A", abs(Fidel_RU)
    return  mps_R, mps_Q, Uni_list, abs(Fidel_RU)

   elif Opt[1]=="Ucontiguous":
    Opt_r="R"
    Opt_u="U"
    #mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt, Fidel_list, count_list)
    mps_R, mps_Q, Uni_list, Fidel_RU=Adding_updateU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt,Fidel_list, count_list)
    print "Fidel=QR-A", abs(Fidel_RU)
    #Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
    #print "Fidel=QR-complete" , abs(Fidel_final)
    return  mps_R, mps_Q, Uni_list, abs(Fidel_RU)


#   file = open("QR.txt", "w")
#   for index in range(len(Fidel_list)):
#      file.write(str(count_list[index]) + " " + str(Fidel_list[index])+" "+ "\n")
#   file.close()
#   return  mps_R, mps_Q, Uni_list, 






 file = open("QR.txt", "w")
 for index in range(len(Fidel_list)):
   file.write(str(count_list[index]) + " " + str(Fidel_list[index])+" "+ "\n")
 file.close()
 #print mps_R.norm()


 #Uni_list[0].setLabel([1,2,3,4]) 
 #U_f=copy.copy(Uni_list[0]) 
 #U_f.setLabel([1,2,5,6]) 

 #print Uni_list[0].printDiagram(), U_f*Uni_list[0]

def make_Env_singleLayer(PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x):

 Peps_ket=[ copy.copy(PEPS_listten[i]) for i in xrange(len(PEPS_listten))]
 Peps_bra=[ copy.copy(PEPS_listten[i]) for i in xrange(len(PEPS_listten))]

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).dim()
  D=Peps_bra[0].bond(3).dim()
 else:
  DX=Peps_ket[0].bond(0).dim()
  D=Peps_bra[0].bond(0).dim()

 #print "Location, D, DX", Location, DX, D


 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 IdenX=uni10.UniTensor([bdiX, bdoX])
 IdenX.identity()
 IdenX.setLabel([4,5])

 Iden=uni10.UniTensor([bdi, bdo])
 Iden.identity()
 Iden.setLabel([7,8])

 IdenX1=uni10.UniTensor([bdiX, bdoX])
 IdenX1.identity()
 IdenX1.setLabel([3,6])

 #####################################Zero################################################################
 bdi_1 = uni10.Bond(uni10.BD_IN, 1)
 bdo_1 = uni10.Bond(uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0

 Tem1=uni10.UniTensor([bdi_1])
 Tem1.identity()
 Tem1.setLabel([1])
 #print Tem1

 Tem11=uni10.UniTensor([bdi_1])
 Tem11.identity()
 Tem11.setLabel([11])
 #print Tem11


 Tem10=uni10.UniTensor([bdi_1])
 Tem10.identity()
 Tem10.setLabel([10])
 #print Tem10


 if Location==0:
  mps_list=[]

  Peps_ket[0].setLabel([0,1,2,3,4])
  Peps_bra[0].setLabel([10,11,2,-3,-4])

  results=((Peps_ket[0]*Tem0)*Tem1)*((Peps_bra[0])*Tem11)
  results.permute([10, 4, 3, -4, -3], 5)
  results.combineBond([4,3])
  results.combineBond([4,-4])
  results.permute([10,-3,4],2)
  mps_list.append(results)
  ######################################################3

  results=(IdenX1*Iden)*IdenX

  results.permute([4,3,7,5,8,6],6)
  results.combineBond([4,3])
  results.combineBond([4,7])

  results.combineBond([5,8])
  results.permute([4,6,5],2)
  mps_list.append(results)
  ##################################Middle##############
  for  q  in  xrange(1,len(Peps_ket)-1):
   #print q
   Peps_ket[q].setLabel([0,1,2,3,4])
   Peps_bra[q].setLabel([10,11,2,-3,-4])
   Tem0.identity()
   Tem0.setLabel([0])
   Tem10.identity()
   Tem10.setLabel([10])

   results=((Peps_ket[q]*Tem0))*((Peps_bra[q])*Tem10)
   results.permute([1,11,4,3,-4,-3],6)
   results.combineBond([4,3])
   results.combineBond([4,-4])
   results.combineBond([1,11])
   results.permute([1,-3,4],2)
   mps_list.append(results)
   ######################################################3
   Iden.setLabel([7,8])
   IdenX.setLabel([4,5])
   IdenX1.setLabel([3,6])

   results=(IdenX1*Iden)*IdenX

   results.permute([4,3,7,5,8,6],6)
   results.combineBond([4,3])
   results.combineBond([4,7])

   results.combineBond([5,8])
   results.permute([4,6,5],2)
   mps_list.append(results)
 
  ##################################LAssssssstTTT#####################################################
  q=len(Peps_bra)-1

  results=(IdenX1*Iden)*IdenX

  results.permute([4,3,7,5,8,6],6)
  results.combineBond([4,3])
  results.combineBond([4,7])

  results.combineBond([5,8])
  results.permute([5,6,4],2)
  mps_list.append(results)

  Tem4=uni10.UniTensor([bdi_1])
  Tem4.identity()
  Tem4.setLabel([4])
  Tem10.identity()
  Tem10.setLabel([10])
  Tem11.identity()
  Tem11.setLabel([11])

  Peps_ket[q].setLabel([0,1,2,3,4])
  Peps_bra[q].setLabel([10,-1,2,-3,11])
  #print Peps_ket[q].printDiagram()
  results=((Peps_ket[q])*Tem4)*(((Peps_bra[q])*Tem10)*Tem11)
  results.permute([1,3,-1,-3,0],4)
  results.combineBond([1,3])
  results.combineBond([1,-1])
  results.permute([1,-3,0],2)
  mps_list.append(results)

  list_bond=[]
  for q in xrange(len(mps_list)):
    list_bond.append(mps_list[q].bond(2).dim())

  mps_boundry=MPSclass.MPS(mps_list[1].bond(1).dim(),max(list_bond),len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=copy.copy(mps_list[i])

  E_0=mps_boundry.norm()
  mps_bound=mps_boundry.D
  if mps_boundry.D>chi_boundry:
    mps_boundry=mps_boundry.appSVD(chi_boundry)
    #print "Trunc", mps_bound, chi_boundry, (E_0-mps_boundry.norm())/E_0
  return mps_boundry





######################## Add next Layer ###########

 if Location!=0:
  mps_list=[None]*mps_boundry.N

###############################First-Layer###############################
  mps_list[0]=mps_boundry[0]*1.0


  mps_boundry[1].setLabel([5,0,6])
  Peps_ket[0].setLabel([0,1,2,3,4])
  Tem1.setLabel([1])


  result=(mps_boundry[1]*(Tem1*Peps_ket[0]))
  result.permute([5,6,2,3,4],5)
  result.combineBond([6,4])
  result.combineBond([3,2])
  result.permute([5,3,6],2)
  mps_list[1]=result*1.0
  ###########################################
  mps_boundry[2].setLabel([5,0,6])

  Iden.setLabel([0,10])
  IdenX.setLabel([-5,-6])

  #print mps_boundry[2].printDiagram(), Iden.printDiagram(), IdenX.printDiagram()
  result=(mps_boundry[2]*(Iden*IdenX))
  result.permute([5,6,-5,-6,10],5)
  #results.combineBond([5])
  result.combineBond([5,-5])
  result.combineBond([6,-6])
  result.permute([5,10,6],2)
  mps_list[2]=result*1.0

  ###############################################################################

  for q in xrange(3, mps_boundry.N,2):

   ###########################################
   if q<(mps_boundry.N-3):
    #print q, (q-1)/2

    mps_boundry[q].setLabel([5,0,6])
    Peps_ket[(q-1)/2].setLabel([0,1,2,3,4])

    result=(mps_boundry[q]*(Peps_ket[(q-1)/2]))
    result.permute([5,1,6,2,3,4],6)
    result.combineBond([5,1])
    result.combineBond([6,4])
    result.combineBond([3,2])
    result.permute([5,3,6],2)
    mps_list[q]=result*1.0

  ################
    mps_boundry[q+1].setLabel([5,0,6])
    Iden.setLabel([0,10])
    IdenX.setLabel([-5,-6])

    result=(mps_boundry[q+1]*(Iden*IdenX))
    result.permute([5,6,-5,-6,10],5)
    result.combineBond([5,-5])
    result.combineBond([6,-6])
    result.permute([5,10,6],2)
    mps_list[q+1]=result*1.0
   elif q==(mps_boundry.N-3):
    #print "inside", q+1, (q-1)/2,(q/2)+1
    mps_boundry[q].setLabel([5,0,6])
    Peps_ket[(q-1)/2].setLabel([0,1,2,3,4])

    result=mps_boundry[q]*Peps_ket[(q-1)/2]
    result.permute([5,1,6,2,3,4],6)
    result.combineBond([5,1])
    result.combineBond([6,4])
    result.combineBond([3,2])
    result.permute([5,3,6],2)
    mps_list[q]=result*1.0


    mps_boundry[q+1].setLabel([5,0,6])
    Peps_ket[(q/2)+1].setLabel([0,1,2,3,4])

    result=mps_boundry[q+1]*Peps_ket[(q/2)+1]
    result.permute([5,1,6,2,3,4],6)
    result.combineBond([5,1])
    result.combineBond([6,4])
    result.combineBond([3,2])
    result.permute([5,3,6],2)
    mps_list[q+1]=result*1.0
   else:
    #print "last", q
    mps_list[q]=mps_boundry[q]*1.0



    #mps_A=mps_A.appSVD(chi_p)

  list_bond=[]
  for q in xrange(len(mps_list)):
    list_bond.append(mps_list[q].bond(2).dim())

  mps_boundry=MPSclass.MPS(mps_list[1].bond(1).dim(),max(list_bond),len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=copy.copy(mps_list[i])

  E_0=mps_boundry.norm()
  mps_bound=mps_boundry.D
  if mps_boundry.D>chi_boundry:
    mps_boundry=mps_boundry.appSVD(chi_boundry)
    #print "Trunc", mps_bound, chi_boundry, (E_0-mps_boundry.norm())/E_0


  #print "test", mps_boundry_new.norm()

  ######################### Next Absorption ####################################

  #mps_boundry=mps_boundry_ne*1.0

  if Location == (N_x-1):
   DX=1
   bdiX=uni10.Bond(uni10.BD_IN, DX)
   bdoX=uni10.Bond(uni10.BD_OUT, DX)

  for  q  in  xrange(0, mps_boundry.N, 2):
   if q==0:
    mps_boundry[q].setLabel([10,3,2])
    Peps_bra[q].setLabel([3,11,1,4,-2])
    Tem1.setLabel([11])
    result=((mps_boundry[q]*Peps_bra[q])*Tem1)
    result.permute([10,4,2,-2,1],5)
    result.combineBond([2,-2])
    result.combineBond([2,1])
    result.permute([10,4,2],2)
    mps_list[q]=result*1.


    t1=uni10.UniTensor([mps_boundry[q+1].bond(0),bdiX,bdiphy,mps_boundry[q+1].bond(2)])
    t1.setLabel([10,3,1,2])
    #print mps_boundry[q+1].printDiagram()
    t1.putBlock(mps_boundry[q+1].getBlock())

     #Iden.setLabel([0,10])
    Iden.setLabel([-5,-6])

    result=(t1*Iden)
    result.permute([10,-5,1,3,2,-6],5)
    result.combineBond([10,-5])
    result.combineBond([10,1])
    result.combineBond([2,-6])
    result.permute([10,3,2],2)
    mps_list[q+1]=result*1.0
   elif q< mps_boundry.N-2:
    #print "midle", q
    mps_boundry[q].setLabel([10,3,2])
    Peps_bra[q/2].setLabel([3,11,1,4,-2])

    result=((mps_boundry[q]*Peps_bra[q/2]))
    result.permute([10,11,4,2,-2,1],5)
    result.combineBond([10,11])
    result.combineBond([2,-2])
    result.combineBond([2,1])
    result.permute([10,4,2],2)
    mps_list[q]=result*1.0

    t1=uni10.UniTensor([mps_boundry[q+1].bond(0),bdiX,bdiphy,mps_boundry[q+1].bond(2)])
    t1.setLabel([10,3,1,2])
    t1.putBlock(mps_boundry[q+1].getBlock())
    Iden.setLabel([-5,-6])
    result=(t1*Iden)
    result.permute([10,-5,1,3,2,-6],5)
    result.combineBond([10,-5])
    result.combineBond([10,1])
    result.combineBond([2,-6])
    result.permute([10,3,2],2)
    mps_list[q+1]=result*1.0
   elif q==mps_boundry.N-2:
    #print "last", q
    t1=uni10.UniTensor([mps_boundry[q].bond(0),bdiX,bdiphy,mps_boundry[q].bond(2)])
    t1.setLabel([10,3,1,2])
    t1.putBlock(mps_boundry[q].getBlock())

    Iden.setLabel([-5,-6])

    result=(t1*Iden)
    result.permute([10,-5,1,3,2,-6],5)
    result.combineBond([10,-5])
    result.combineBond([2,-6])
    result.combineBond([2,1])
    result.permute([10,3,2],2)
    mps_list[q]=result*1.0

    #t1=uni10.UniTensor([mps_boundry[q].bond(0),bdiphy,bdiX,mps_boundry[q].bond(2)])
    mps_boundry[q+1].setLabel([10,3,2])
    Peps_bra[q/2].setLabel([3,-2,1,4,11])
    Tem1.setLabel([11])
    #print "Hi",(q+1)/2#, Peps_bra[5],Peps_bra[6], Peps_bra[7] 
    result=(mps_boundry[q+1]*(Peps_bra[q/2]*Tem1))
    result.permute([10,-2,1,4,2],5)
    result.combineBond([10,-2])
    result.combineBond([10,1])
    #result.combineBond([2,1])
    result.permute([10,4,2],2)
    mps_list[q+1]=result*1.0

  list_bond=[]
  for q in xrange(len(mps_list)):
    list_bond.append(mps_list[q].bond(2).dim())

  mps_boundry=MPSclass.MPS(mps_list[1].bond(1).dim(),max(list_bond),len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=copy.copy(mps_list[i])

  mps_bound=mps_boundry.D
  E_0=mps_boundry.norm()
  if mps_boundry.D>chi_boundry:
    mps_boundry=mps_boundry.appSVD(chi_boundry)
    #print "Trunc", mps_bound, chi_boundry, (E_0-mps_boundry.norm())/E_0

  return mps_boundry

def rotate(PEPS_list):

  PEPS_list_copy=[copy.copy(PEPS_list[i])   for i in xrange(len(PEPS_list))]

  for i in xrange(len(PEPS_list)):
    PEPS_list_copy[i].setLabel([0,1,2,3,4])
    PEPS_list_copy[i].permute([3,1,2,0,4],3)
    PEPS_list_copy[i].setLabel([0,1,2,3,4])
  return PEPS_list_copy


def  Cal_norm(PEPS_listten, N_x, D, chi_try, d):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*N_x
 mps_boundry_temp=[None]*N_x

 for Location in xrange(N_x-1):
  mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_try, N_x)

 Location=0
 PEPS_ten=rotate(PEPS_listten[N_x-Location-1])
 mps_boundry_temp[Location]=make_Env_singleLayer(PEPS_ten, Location, mps_boundry_temp[Location-1],d,chi_try,N_x)
 mps_boundry_right[N_x-Location-1]=copy.copy(mps_boundry_temp[Location])

 norm_val=mps_boundry_left[N_x-2].product(mps_boundry_right[N_x-1])

 return   norm_val



def  Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval):

 norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

 count=0
 while count<100:
  print "norm_val", norm_val
  if abs(norm_val)<threshold and abs(norm_val)>(1.0/threshold):  break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(1.0/threshold):
   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

 PEPS_listten=All_dist(PEPS_listten,N_x, D)
 #print "Fixed norm", abs(norm_val)
 return PEPS_listten, abs(norm_val), count


def All_dist(PEPS_listten,N_x, D):
 for i in xrange(N_x-1):
  for j in xrange(N_x):
   PEPS_listten[i][j], PEPS_listten[i+1][j]=equall_dis_H(PEPS_listten[i][j], PEPS_listten[i+1][j],D)
   PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   PEPS_listten[i+1][j]=max_ten(PEPS_listten[i+1][j])
 
 for i in xrange(N_x):
  for j in xrange(N_x-1):
   PEPS_listten[i][j], PEPS_listten[i][j+1]=equall_dis_V(PEPS_listten[i][j], PEPS_listten[i][j+1],D)
   PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   PEPS_listten[i][j+1]=max_ten(PEPS_listten[i][j+1])
 return PEPS_listten

def equall_dis_H(PEPS_1, PEPS_2, D):

 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,5,3,4],3)
 
 #q, r= qr_parity1(A) 
 q,s,V=svd_parity5(A)
 V.setLabel([0,3,4])
 s.setLabel([-100,0])
 V=s*V
 V.permute([-100,3,4],2)
 r=V*1.0
 r.setLabel([-100,3,4])
 q.setLabel([1,2,5,-100])

 A=copy.copy(PEPS_2)
 A.setLabel([4,5,6,7,8])
 A.permute([4, 6, 5,8,7],2)
 U,s,qq=svd_parity6(A) 
 U.setLabel([4,6,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([4,6,-200],2)
 l=U*1.0
 qq.setLabel([-200,5,8,7])
 l.setLabel([4,6,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,6],2)
 U, V, s= setTruncation3(Teta, D)

 #U,s,V=svd_parity2(Teta)

 U.setLabel([-100,3,17])
 s.setLabel([17,-17])
 V.setLabel([-17,-200,6])
 s=Sqrt(s)
 s.setLabel([17,-17])
 
 U=U*s
 V=s*V

 U.permute([-100,3,-17],1)
 U.setLabel([-100,3,4])
 V.permute([17,-200,6],1)
 V.setLabel([4,-200,6])
 
 PEPS_1=q*U
 PEPS_2=qq*V

 PEPS_1.permute([1,2,3,4,5],3)
 PEPS_2.permute([4,5,6,7,8],3)

 return PEPS_1, PEPS_2


def equall_dis_V(PEPS_1, PEPS_2, D):

 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,4,3,5],3)
 
 #q, r= qr_parity1(A) 
 q,s,V=svd_parity5(A)
 V.setLabel([0,3,5])
 s.setLabel([-100,0])
 V=s*V
 V.permute([-100,3,5],2)
 r=V*1.0
 r.setLabel([-100,3,5])
 q.setLabel([1,2,4,-100])

 A=copy.copy(PEPS_2)
 A.setLabel([6,5,7,8,9])
 A.permute([5, 7, 6,8,9],2)
 U,s,qq=svd_parity6(A) 
 U.setLabel([5,7,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([5,7,-200],2)
 l=U*1.0
 qq.setLabel([-200,6,8,9])
 l.setLabel([5,7,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,7],2)
 U, V, s= setTruncation3(Teta, D)

 #U,s,V=svd_parity2(Teta)

 U.setLabel([-100,3,17])
 s.setLabel([17,-17])
 V.setLabel([-17,-200,7])
 s=Sqrt(s)
 s.setLabel([17,-17])
 
 U=U*s
 V=s*V

 U.permute([-100,3,-17],1)
 U.setLabel([-100,3,5])
 V.permute([17,-200,7],1)
 V.setLabel([5,-200,7])
 
 PEPS_1=q*U
 PEPS_2=qq*V

 PEPS_1.permute([1,2,3,4,5],3)
 PEPS_2.permute([6,5,7,8,9],3)

 return PEPS_1, PEPS_2


def svd_parity5(theta):

 bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).dim()*theta.bond(4).dim())
 bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(3).dim()*theta.bond(4).dim())

 bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())
 bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())

 if bdi.dim()<=bdi1.dim():
  #print theta.printDiagram()
  GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo])
  LA=uni10.UniTensor([bdi,bdo])
  GB=uni10.UniTensor([bdi,theta.bond(3),theta.bond(4)])

  svds = {}
  blk_qnums = theta.blockQnum()
  dim_svd=[]

  for qnum in blk_qnums:
      svds[qnum] = theta.getBlock(qnum).svd()
      #print   svds[qnum][0].row(), svds[qnum][0].col()
      GA.putBlock(qnum, svds[qnum][0])
      LA.putBlock(qnum, svds[qnum][1])
      GB.putBlock(qnum, svds[qnum][2])
 else:
  #print theta.printDiagram()
  GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo1])
  LA=uni10.UniTensor([bdi1,bdo1])
  GB=uni10.UniTensor([bdi1,theta.bond(3),theta.bond(4)])

  svds = {}
  blk_qnums = theta.blockQnum()
  dim_svd=[]
  for qnum in blk_qnums:
      svds[qnum] = theta.getBlock(qnum).svd()
      #print   svds[qnum][0].row(), svds[qnum][0].col()
      GA.putBlock(qnum, svds[qnum][0])
      LA.putBlock(qnum, svds[qnum][1])
      GB.putBlock(qnum, svds[qnum][2])

 return GA, LA, GB


def svd_parity6(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(2),theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(2),theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB



def     setTruncation3(theta, chi):
 LA=uni10.UniTensor(theta.bond())
 GA=uni10.UniTensor(theta.bond())
 GB=uni10.UniTensor(theta.bond())
 svds = {}
 blk_qnums = theta.blockQnum()
 dim_svd=[]
 for qnum in blk_qnums:
     svds[qnum] = theta.getBlock(qnum).svd()
     dim_svd.append(int(svds[qnum][1].col()))
 svs = []
 bidxs = []
 for bidx in xrange(len(blk_qnums)):
     svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
 dims = [0] * len(blk_qnums)
 for bidx in bidxs:
     dims[bidx] += 1  
 qnums = []
 for bidx in xrange(len(blk_qnums)):
     qnums += [blk_qnums[bidx]] * dims[bidx]
 bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
 bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
 GA.assign([theta.bond(0),theta.bond(1), bdo_mid])
 GB.assign([bdi_mid, theta.bond(2), theta.bond(3)])
 LA.assign([bdi_mid, bdo_mid])
 degs = bdi_mid.degeneracy()
 for qnum, dim in degs.iteritems():
     if qnum not in svds:
         raise Exception("In setTruncaton(): Fatal error.")
     svd = svds[qnum]
     GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
     GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
     LA.putBlock(qnum, svd[1].resize(dim, dim)  )
 return GA, GB, LA

def sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs


def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 max_list_f=[abs(x) for x in max_list]
 return max(max_list_f)

def max_ten(a):
 Max_val=MaxAbs(a)
 if ( Max_val < 0.5e-1) or (Max_val > 0.5e+1)   :

  if Max_val >= 1:
   print ">1",Max_val
   a=a*(1.00/Max_val)
  if Max_val < 1: 
   print "<1",Max_val
   a=a*(1.00/Max_val)

#  if Max_val >= 1:
#   print ">1",Max_val, Max_val**(1./2.)
#   a=a*(1.00/(Max_val**(1./2.)))
#  if Max_val < 1: 
#   print "<1",Max_val
#   a=a*(1.00/Max_val)

 else: a=a;
 return a






















