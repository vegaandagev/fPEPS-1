import pyuni10 as uni10
import copy
import numpy as np
#from numpy import linalg as LA
import root
import MPSclass2
import scipy.linalg as linalg







def  shrink_mps(mps, Dp, D):
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 A_list=[]
 for i in xrange(0,mps.N,2):
  #print i
  Ten_1=uni10.UniTensorR([mps[i].bond(0), bdip, bdi, mps[i].bond(2)], "Ten_R")
  #print mps[i].printDiagram()
  Ten_1.SetLabel([1,2,3,4])
  Ten_1.PutBlock(mps[i].GetBlock())
  Ten_2=uni10.UniTensorR([mps[i+1].bond(0), bdip, bdi, mps[i+1].bond(2)], "Ten_R")
  Ten_2.PutBlock(mps[i+1].GetBlock())
  Ten_2.SetLabel([4,-2,-3,-4])
  Results=uni10.Contract(Ten_1,Ten_2)
  Results=uni10.Permute(Results,[1,2,-2,3,-3,-4],5)  
  Results=Results.CombineBond([2,-2])
  Results=Results.CombineBond([3,-3])
  Results=Results.CombineBond([2,3])
  Results=uni10.Permute(Results,[1,2,-4],2)  
  A_list.append(Results)
    
 list_bond=[]
 for q in xrange(mps.N/2):
   list_bond.append(A_list[q].bond(2).dim())

 mps_R=MPSclass2.MPS(A_list[1].bond(1).dim(),max(list_bond),mps.N/2)

 for i in xrange(mps.N/2):
   mps_R[i]=A_list[i]*1.0

 return   mps_R


def  absorption_left( PEPS_colmps, q, MPS_R, D, d, N_x, N_y, chi_p):

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 PEPS_ten=mps_to_tensor_left(PEPS_colmps, N_x, N_y, D, d, q)
 A_list=[]

 for i in xrange(len(PEPS_ten)):

   PEPS_ten[i].SetLabel([0,1,2,3,4])

   Ten_R=uni10.UniTensorR([MPS_R[i].bond(0), bdi, bdi, MPS_R[i].bond(2)], "Ten_R")
   Ten_R.PutBlock(MPS_R[i].GetBlock())
   Ten_R.SetLabel([-1,-2,0,-3])
   #Ten_R.uni10.Permute([-1,0,-2,-3],3)
   #Ten_R.SetLabel([-1,-2,0,-3])

   result=uni10.Contract(Ten_R, PEPS_ten[i])
   result=uni10.Permute( result,[-2,-1,1,2,3,-3,4],4)
   result=result.CombineBond([-1,1])
   result=result.CombineBond([-3,4])
   result=uni10.Permute(result,[-2,-1,2,3,-3],3)
   result=uni10.Permute(result,[-1,-2,2,3,-3],3)
   result=result.CombineBond([-2,2])
   result=result.CombineBond([-2,3])
   result=uni10.Permute(result,[-1,-2,-3],2)
   A_list.append(result)

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(A_list[i].bond(2).dim())

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y)
 for i in xrange(N_y):
     mps_A[i]=A_list[i]*1.0
 print "withouth Truncation", mps_A.norm()
 if mps_A.D>chi_p:
   mps_A=mps_A.appSVD(chi_p)
 print "After Truncation", mps_A.norm()
 return  mps_A

def  absorption_right( PEPS_colmps, q, MPS_R, D, d, N_x, N_y, chi_p):

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 PEPS_ten=mps_to_tensor_right(PEPS_colmps, N_x, N_y, D, d, q)
 A_list=[]

 for i in xrange(len(PEPS_ten)):

   PEPS_ten[i].SetLabel([0,1,2,3,4])

   Ten_R=uni10.UniTensorR([MPS_R[i].bond(0), bdi, bdi, MPS_R[i].bond(2)], "Ten_R")
   Ten_R.PutBlock(MPS_R[i].GetBlock())
   Ten_R.SetLabel([-1,-2,3,-3])

   result=uni10.Contract(Ten_R, PEPS_ten[i])
   result=uni10.Permute( result,[0,-1,1,2,-2,-3,4],4)
   result=result.CombineBond([-1,1])
   result=result.CombineBond([-3,4])
   result=uni10.Permute(result,[0,-1,2,-2,-3],3)
   result=uni10.Permute(result,[-1,-2,2,0,-3],3)
   result=result.CombineBond([-2,2])
   result=result.CombineBond([-2,0])
   result=uni10.Permute(result,[-1,-2,-3],2)
   A_list.append(result)

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(A_list[i].bond(2).dim())

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y)
 for i in xrange(N_y):
     mps_A[i]=A_list[i]*1.0
 print "withouth Truncation", mps_A.norm()
 if mps_A.D>chi_p:
   mps_A=mps_A.appSVD(chi_p)
 print "After Truncation", mps_A.norm()
 return  mps_A




def   make_ENVMPS(MPO_Ten, N_y, N_x, chi):

 Env_left=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
  A=MPO_Ten[0][j]*1.0
  A.SetLabel([1,2,3,4])
  A=uni10.Permute(A,[2,3,1,4],2)
  A=A.CombineBond([3,1])
  A=uni10.Permute(A,[2,3,4],2)
  cont_list[j]=A*1.0

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(cont_list[q].bond(2).dim())
 mps_A=MPSclass2.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y)
 for i in xrange(N_y):
   mps_A[i]=cont_list[i]*1.0

 if mps_A.D>chi:
   mps_A=mps_A.appSVD(chi)

 Env_left[0]=mps_A

 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=Env_left[q][j]*1.0
   A.SetLabel([2,3,4])
   B=MPO_Ten[q+1][j]*1.0
   B.SetLabel([3,5,6,7])
   Result=uni10.Contract(A, B)
   Result=uni10.Permute(Result,[2,5,6,4,7],3)
   Result=Result.CombineBond([2,5])
   Result=Result.CombineBond([4,7])
   Result=uni10.Permute(Result,[2,6,4],2)
   Result.SetLabel([2,3,4])
   Result=uni10.Permute(Result,[2,3,4],2)
   cont_list[j]=Result*1.0

  list_bond=[]
  for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

  mps_A=MPSclass2.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y)

  for i in xrange(N_y):
    mps_A[i]=cont_list[i]*1.0

  if mps_A.D>chi and q!=(N_x-2):
     mps_A=mps_A.appSVD(chi)

  Env_left[q+1]=mps_A

###############################

 Env_right=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
   A=MPO_Ten[N_x-1][j]*1.0
   A.SetLabel([1,2,3,4])
   A=uni10.Permute(A,[2,1,3,4],2)
   A=A.CombineBond([1,3])
   A=uni10.Permute(A,[2,1,4],2)
   cont_list[j]=A*1.0

 list_bond=[]
 for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

 mps_A=MPSclass2.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y)
 for i in xrange(N_y):
   mps_A[i]=cont_list[i]*1.0

 if mps_A.D>chi:
   mps_A=mps_A.appSVD(chi)


 Env_right[N_x-1]=mps_A


 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=Env_right[N_x-1-q][j]*1.0
   A.SetLabel([2,3,4])
   B=MPO_Ten[N_x-2-q][j]*1.0
   B.SetLabel([5,6,3,7])
   Result=uni10.Contract(A, B)
   Result=uni10.Permute(Result,[6,2,5,7,4],3)
   Result=Result.CombineBond([6,2])
   Result=Result.CombineBond([7,4])
   Result=uni10.Permute(Result,[6,5,7],2)
   cont_list[j]=Result*1.0
  #print  cont_list[0].bond(1).dim(), cont_list[0].bond(2).dim()
  list_bond=[]
  for i in xrange(N_y):
    list_bond.append(cont_list[i].bond(2).dim())

  mps_A=MPSclass2.MPS( cont_list[0].bond(1).dim(), max(list_bond), N_y)
  for i in xrange(N_y):
   mps_A[i]=cont_list[i]*1.0

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
   A=PEPS_listten[i][j]*1.0
   A.SetLabel([1,2,3,4,5])
   A_t=A*1.0
   A_t.SetLabel([-1,-2,3,-4,-5])
   A=uni10.Contract(A, A_t)
   A=uni10.Permute(A,[1,-1,2,-2,4,-4,5,-5],4)
   A=A.CombineBond([1,-1])
   A=A.CombineBond([2,-2])
   A=A.CombineBond([4,-4])
   A=A.CombineBond([5,-5])
   A=uni10.Permute(A, [1,2,4,5], 2)
   MPO[i][j]=A*1.0

 return   MPO



def mps_to_tensor_left(PEPS_listmps, N_x,N_y,D, d,q):

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
     A=uni10.UniTensorR([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 elif q==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 return   B_list

def mps_to_tensor_right(PEPS_listmps, N_x,N_y,D, d,q):

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
     A=uni10.UniTensorR([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 elif q==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,3,2,0,4],4)
     A.PutBlock(PEPS_listmps[i].GetBlock())
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 return   B_list



def  tensor_to_mps_left(B_list, N_x, N_y):

 A_list=[None]*N_y
 for i in xrange(N_y):
  A=B_list[i]*1.0
  A.SetLabel([0,1,2,3,4])
  A=uni10.Permute(A,[1,0,2,3,4],4)
  A=A.CombineBond([0,2])
  A=A.CombineBond([0,3])
  A=uni10.Permute(A,[1,0,4],2)
  A_list[i]=A

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
   mps_A[i]=A_list[i]*1.0
 return  mps_A



def  tensor_to_mps_right(B_list, N_x, N_y):

 A_list=[None]*N_y
 for i in xrange(N_y):
  A=B_list[i]*1.0
  A.SetLabel([0,1,2,3,4])
  A=uni10.Permute(A,[1,3,2,0,4],4)
  A=A.CombineBond([3,2])
  A=A.CombineBond([3,0])
  A=uni10.Permute(A,[1,3,4],2)
  A_list[i]=A

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
   mps_A[i]=A_list[i]*1.0
 return  mps_A



def Init_PEPS( N_y, N_x, D, d, q, distribution, seed_num):
 #print distribution
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 A_list=[None]*N_y
 B_list=[None]*N_y

 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")
     A.Randomize(UorN=distribution, dn_mu=-10, up_var=10, seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
     #print A_list[i].printDiagram()
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=uni10.UniTensorR([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A

 elif q==N_x-1:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
     #print A_list[i].printDiagram()
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensorR([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
     #print A_list[i].printDiagram()
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
     A.Randomize(UorN=distribution,dn_mu=-10,up_var=10,seed=seed_num)
     #A.orthoRand()
     B_list[i]=A*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
     mps_A[i]=A_list[i]*1.0

 norm=0.5*mps_A.norm()
 for q in xrange(len(B_list)):
   B_list[q]=B_list[q]*(1/(norm**(0.5/N_y)))


 if q == 0:
  for i in xrange(N_y):
   if i == 0:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     a=A.uni10.Permute([1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A

 elif q==N_x-1:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensorR([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.Randomize(dn_mu=-1,up_var=1,seed=1)
     #A.orthoRand()
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   elif i ==(N_y-1):
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A
   else:
     A=B_list[i]*1.0
     A.SetLabel([0,1,2,3,4])
     A=uni10.Permute(A,[1,0,2,3,4],2)
     A=A.CombineBond([0,2])
     A=A.CombineBond([0,3])
     A=uni10.Permute(A,[1,0,4],2)
     A_list[i]=A

 mps_A=MPSclass2.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N_y)
 for i in xrange(N_y):
     mps_A[i]=A_list[i]*1.0

 return mps_A


def inverse(Landa2):
 invLanda2=uni10.UniTensorR(Landa2.bond())
 blk_qnums=Landa2.BlocksQnum()
 for qnum in blk_qnums:
  D=int(Landa2.GetBlock(qnum).shape[0])
  invLt=Landa2.GetBlock(qnum)
  invL2=invLt*0
  #invLt.fill(0)
  for i in xrange(D):
      invL2[i][i] = 0 if ((invLt[i][i]) < 1.0e-10) else (1.00 / (invLt[i][i]))

  invLanda2.PutBlock(qnum,invL2)
 return invLanda2


#########  prerequisite functions  #############
def   setTruncation1(theta, chi):
    LA=uni10.UniTensorR(theta.bond())
    GA=uni10.UniTensorR(theta.bond())
    GB=uni10.UniTensorR(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.GetBlock(qnum).svd()
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
        GA.PutBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.PutBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.PutBlock(qnum, svd[1].resize(dim, dim)  )
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
#     A=uni10.UniTensorR([bdi1,bdip,bdi,bdo], "A_first")
#     A.Randomize(dn_mu=-1,up_var=1,seed=1)
#     #A.orthoRand()
#     A.SetLabel([0,1,2,3])
#     A.CombineBond([1,2])
#     #A.uni10.Permute([1,0,4],2)
#     A_list[i]=A
#     #print A_list[i].printDiagram()
#   elif i ==(N-1):
#     A=uni10.UniTensorR([bdi,bdip,bdi,bdo1], "A_end")
#     A.Randomize(dn_mu=-1,up_var=1,seed=1)
#     #A.orthoRand()
#     A.SetLabel([0,1,2,3])
#     A.CombineBond([1,2])
#     A_list[i]=A
#     #print A_list[i].printDiagram()
#   else:
#     A=uni10.UniTensorR([bdi,bdip,bdi,bdo], "A_middle")
#     A.Randomize(dn_mu=-1,up_var=1,seed=1)
#     #A.orthoRand()
#     A.SetLabel([0,1,2,3])
#     A.CombineBond([1,2])
#     A_list[i]=A
#     #print A_list[i].printDiagram()
# 
#  mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N, 'ortho')
#  for i in xrange(N):
#    mps_A[i]=copy.copy(A_list[i])
#  return mps_A


def svd_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())

    GA=uni10.UniTensorR([theta.bond(0),theta.bond(1)])
    LA=uni10.UniTensorR([bd1,theta.bond(1)])
    GB=uni10.UniTensorR([bd1,theta.bond(1)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.GetBlock(qnum).svd()
        GA.PutBlock(qnum, svds[qnum][0])
        LA.PutBlock(qnum, svds[qnum][1])
        GB.PutBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity1(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())
    #print bdi1
    GA=uni10.UniTensorR([theta.bond(0),theta.bond(1), bdo1])
    LA=uni10.UniTensorR([bdi1,bdo1])
    GB=uni10.UniTensorR([bdi1,theta.bond(2),theta.bond(3)])

    svds = {}
    blk_qnums = theta.BlocksQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        #svds[qnum] = theta.GetBlock(qnum).svd()
        svds[qnum] = linalg.svd(theta.GetBlock(qnum), full_matrices=False, lapack_driver='gesvd')
        GA.PutBlock(qnum, svds[qnum][0])
        LA.PutBlock(qnum, np.diag(svds[qnum][1]))
        GB.PutBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity2(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())

    if bdi.dim()<=bdi1.dim():
     #print theta.printDiagram()
     GA=uni10.UniTensorR([theta.bond(0),theta.bond(1), theta.bond(2), bdo])
     LA=uni10.UniTensorR([bdi,bdo])
     GB=uni10.UniTensorR([bdi,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.BlocksQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         #svds[qnum] = theta.GetBlock(qnum).svd()
         svds[qnum] = linalg.svd(theta.GetBlock(qnum), full_matrices=False, lapack_driver='gesvd')

         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.PutBlock(qnum, svds[qnum][0])
         LA.PutBlock(qnum, np.diag(svds[qnum][1]))
         GB.PutBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensorR([theta.bond(0),theta.bond(1), theta.bond(2), bdo1])
     LA=uni10.UniTensorR([bdi1,bdo1])
     GB=uni10.UniTensorR([bdi1,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.BlocksQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         #svds[qnum] = theta.GetBlock(qnum).svd()
         svds[qnum] = linalg.svd(theta.GetBlock(qnum), full_matrices=False, lapack_driver='gesvd')
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.PutBlock(qnum, svds[qnum][0])
         LA.PutBlock(qnum, np.diag(svds[qnum][1]))
         GB.PutBlock(qnum, svds[qnum][2])

    return GA, LA, GB

# def svd_parity3(theta):
# 
#  bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
#  bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
#  bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
#  #bd4=uni10.Bond(uni10.BD_IN,theta.bond(7).Qlist())
# 
#  GA=uni10.UniTensorR([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
#  LA=uni10.UniTensorR([bd1,bd2,bd3,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
#  GB=uni10.UniTensorR([bd1,bd2,bd3,bd4,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
# 
#  svds = {}
#  blk_qnums = theta.blockQnum()
#  dim_svd=[]
#  for qnum in blk_qnums:
#      svds[qnum] = theta.GetBlock(qnum).svd()
#      GA.PutBlock(qnum, svds[qnum][0])
#      LA.PutBlock(qnum, svds[qnum][1])
#      GB.PutBlock(qnum, svds[qnum][2])
# 
# #    print LA
#  return GA, LA, GB


def Cost_Function(T0, U, U_transpose):

 #print  U_transpose*U
 T0p=T0*1.0
 T0.SetLabel([0,1,2,3])
 T0p.SetLabel([0,1,2,3])
 Norm0=uni10.Contract(T0p, T0)

 T0.SetLabel([0,1,2,3])
 T0p.SetLabel([0,4,2,3])
 U.SetLabel([1,-3])
 U_transpose.SetLabel([4,-3])

 A_t=uni10.Contract(U_transpose,U)
 A_t=uni10.Contract(A_t, T0)
 Norm1=uni10.Contract(A_t, T0p)

 #Norm1=(T0*(U_transpose*U))*(T0p)
 #print Norm1[0], Norm0[0] 
 result=Norm0.GetBlock()[0,0]-Norm1.GetBlock()[0,0]

 return result

def Env_U(T0, U, U_transpose):
 T0p=copy.copy(T0)
 T0.SetLabel([0,1,2,3])
 T0p.SetLabel([0,4,2,3])
 U.SetLabel([1,-3])
 U_transpose.SetLabel([4,-3])
 #print T0.printDiagram(), T0p.printDiagram(), U_transpose.printDiagram()
 Y=(T0*U_transpose)*(T0p)
 Y.uni10.Permute([1,-3],1)
 return Y



def Env_U(T0, U, U_transpose):
 T0p=T0*1.0
 T0.SetLabel([0,1,2,3])
 T0p.SetLabel([0,4,2,3])
 U.SetLabel([1,-3])
 U_transpose.SetLabel([4,-3])
 #print T0.printDiagram(), T0p.printDiagram(), U_transpose.printDiagram()
 A_t=uni10.Contract(T0,U_transpose)
 Y=uni10.Contract(A_t, T0p)
 #Y=(T0*U_transpose)*(T0p)
 Y=uni10.Permute(Y,[1,-3],1)
 return Y




def   Cost_Function_root(T0,T1,t0,t1,Uni,location):


  T0.SetLabel([1,5,7,-2])
  T1.SetLabel([-2,6,8,4])

  t0.SetLabel([1,2,7,-3])
  t1.SetLabel([-3,3,8,4])
  Uni.SetLabel([5,6,2,3])

  Norm=uni10.Contract(uni10.Contract(uni10.Contract(t0,Uni),t1),uni10.Contract(T0,T1))
  #print Norm
  return   Norm.GetBlock()[0,0]


def Env_U_root(T0,T1,t0,t1,Uni,location):
  T0.SetLabel([1,5,7,-2])
  T1.SetLabel([-2,6,8,4])

  t0.SetLabel([1,2,7,-3])
  t1.SetLabel([-3,3,8,4])
  Uni.SetLabel([5,6,2,3])

  Norm=uni10.Contract(uni10.Contract(t0,t1),uni10.Contract(T0,T1))
  #print Norm

  Norm=uni10.Permute(Norm,[5,6,2,3],2)

  return  Norm


def optimized_root(T0,T1,t0,t1, Uni, location):
  U_update=Uni*1.0
  U=Uni
  E2=10
  fidel=0
  for i in xrange(3):
    E1=Cost_Function_root(T0,T1,t0,t1,U,location)
    #print "E1", E1
    if E1< 1.0e-14:
      #print "E1< 1.0e-14"
      break
    #print i, E1, E2, abs((E2-E1)/E1)
    fidel=E1
    if E1<E2 or i is 0:
     U_update=U*1.0
     if abs((E2-E1)/E1) < 1.0e-7:
      #print E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
    Y=Env_U_root( T0, T1, t0, t1, U, location)
    #print Y.printDiagram() 

    #svd=Y.GetBlock().svd()
    #temporary_matrix=svd[0]*svd[2]

    svd=linalg.svd(Y.GetBlock(), full_matrices=False, lapack_driver='gesvd')
    temporary_matrix=np.matmul(svd[0], svd[2])


    U.PutBlock(temporary_matrix)
    E2=copy.copy(E1)

  #E1=Cost_Function_root(T0,T1,t0,t1,U_update,location)
  #print "simpleRoot-Dis=", fidel
  return U_update



def optimized_U(T0, U, U_transpose):
  U_update=U*1.0
  E2=10
  fidel=0
  for i in xrange(30):
    E1=Cost_Function(T0, U, U_transpose)
    if E1< 1.0e-14:
      #print "E1< 1.0e-14"
      break
    #print i, E1, E2, abs((E2-E1)/E1)
    fidel=E1
    if E1<E2 or i==0:
     U_update=U*1.0
     if abs((E2-E1)/E1) < 1.0e-7:
      #print E2, E1, abs((E2-E1)/E1), i
      break
    else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
    Y=Env_U(T0, U, U_transpose)
    #print Y.printDiagram() 
    #svd=Y.GetBlock().svd()
    #temporary_matrix=svd[0]*svd[2]
    svd = linalg.svd(Y.GetBlock(), full_matrices=False, lapack_driver='gesvd')
    temporary_matrix=np.matmul(svd[0], svd[2])
    U.PutBlock(temporary_matrix)
    U.SetLabel([1,-3])
    U_transpose=U*1.0
    U_transpose.SetLabel([4,-3])
    E2=copy.copy(E1)

  U_update.SetLabel([1,-3])
  U_update_transpose=U_update*1.0
  U_update_transpose.SetLabel([4,-3])
  E1=Cost_Function(T0, U, U_transpose)
  return U_update

def  Cost_FunctionAll_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])


 Iso0.SetLabel([2,10])
 Iso1.SetLabel([5,11])
 #Iso0t.SetLabel([10,8])
 #Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0))*((Iso1)))*(U))*(t0*t1)

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0))*((Iso1))))*(t0*t1)

 Norm2=(t0*t1)*(t0*t1)

 return Norm2[0]+Norm1[0]-2*Norm0[0]

def  Env_Iso0_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])

 Iso0.SetLabel([2,10])
 Iso1.SetLabel([5,11])
 #Iso0t.SetLabel([10,8])
 #Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])

 Norm0=((((T1*T0))*((Iso1)))*(U))*(t0*t1)
 Norm0.uni10.Permute([2,10],1)
 #print Norm0

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])

 Norm1=(((T1*T0))*Iso1)*(t0*t1)
 Norm1.uni10.Permute([2,10],1)

 #print Norm1
 #Norm2=(t0*t1)*(t0*t1)
 return Norm1+(-2)*Norm0

def  Env_Iso1_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])


 Iso0.SetLabel([2,10])
 Iso1.SetLabel([5,11])
 #Iso0t.SetLabel([10,8])
 #Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0)))*(U))*(t0*t1)
 Norm0.uni10.Permute([5,11],1)

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])

 Norm1=((((T1*T0)*(Iso0))))*(t0*t1)
 Norm1.uni10.Permute([5,11],1)

 return Norm1+(-2)*Norm0

def  Env_uall_chi(T0,T1,U,Iso0,Iso1):

 t0=copy.copy(T0)
 t1=copy.copy(T1)
 Iso0t=copy.copy(Iso0)
 Iso1t=copy.copy(Iso1)

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])

 Iso0.SetLabel([2,10])
 Iso1.SetLabel([5,11])
 #Iso0t.SetLabel([10,8])
 #Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])

 Norm0=((((T1*T0)*(Iso0))*((Iso1))))*(t0*t1)
 Norm0.uni10.Permute([12,13,10,11],2)

 return Norm0


def  Cost_FunctionAll(T0,T1,U,Iso0,Iso1):

 t0=T0*1.0
 t1=T1*1.0
 Iso0t=Iso0*1.0
 Iso1t=Iso1*1.0

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])


 Iso0.SetLabel([2,8])
 Iso1.SetLabel([5,9])
 Iso0t.SetLabel([10,8])
 Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])


 A_t=uni10.Contract(T1,T0)
 B_t=uni10.Contract(Iso0,Iso0t)
 A_t=uni10.Contract(A_t, B_t)
 B_t=uni10.Contract(Iso1,Iso1t)
 A_t=uni10.Contract(A_t, B_t)
 A_t=uni10.Contract(A_t, U)
 A_t=uni10.Contract(A_t, t0)
 Norm0=uni10.Contract(A_t, t1)


 #Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t)))*(U))*(t0*t1)

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])


 A_t=uni10.Contract(T1,T0)
 B_t=uni10.Contract(Iso0,Iso0t)
 A_t=uni10.Contract(A_t, B_t)
 B_t=uni10.Contract(Iso1,Iso1t)
 A_t=uni10.Contract(A_t, B_t)
 A_t=uni10.Contract(A_t, t0)
 Norm1=uni10.Contract(A_t, t1)

 #print Norm1
# Norm1=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t))))*(t0*t1)

 A_t=uni10.Contract(t1,t0)
 B_t=uni10.Contract(t1,t0)
 Norm2=uni10.Contract(A_t, B_t)


 #Norm2=(t0*t1)*(t0*t1)
 result=Norm2.GetBlock()[0,0]+Norm1.GetBlock()[0,0]-2*Norm0.GetBlock()[0,0]
 return result



def  Env_Iso0(T0,T1,U,Iso0,Iso1):

 t0=T0*1.0
 t1=T1*1.0
 Iso0t=Iso0*1.0
 Iso1t=Iso1*1.0

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])

 Iso0.SetLabel([2,8])
 Iso1.SetLabel([5,9])
 Iso0t.SetLabel([10,8])
 Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])


 A_t=uni10.Contract(T1,T0)
 A_t=uni10.Contract(A_t, Iso0t)
 A_t=uni10.Contract(A_t, Iso1)
 A_t=uni10.Contract(A_t, Iso1t)
 A_t=uni10.Contract(A_t, U)
 A_t=uni10.Contract(A_t, t0)
 Norm0=uni10.Contract(A_t, t1)


 #Norm0=((((T1*T0)*(Iso0t))*((Iso1*Iso1t)))*(U))*(t0*t1)
 Norm0=uni10.Permute(Norm0,[2,8],1)
 #print Norm0

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])


 A_t=uni10.Contract(T1,T0)
 A_t=uni10.Contract(A_t, Iso0t)
 A_t=uni10.Contract(A_t, Iso1)
 A_t=uni10.Contract(A_t, Iso1t)
 #A_t=uni10.Contract(A_t, U)
 A_t=uni10.Contract(A_t, t0)
 Norm1=uni10.Contract(A_t, t1)

# Norm1=((((T1*T0)*(Iso0t))*((Iso1*Iso1t))))*(t0*t1)
 Norm1=uni10.Permute(Norm1,[2,8],1)

 #print Norm1
 #Norm2=(t0*t1)*(t0*t1)
 return Norm1+(-2)*Norm0

def  Env_Iso1(T0,T1,U,Iso0,Iso1):

 t0=T0*1.0
 t1=T1*1.0
 Iso0t=Iso0*1.0
 Iso1t=Iso1*1.0

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])


 Iso0.SetLabel([2,8])
 Iso1.SetLabel([5,9])
 Iso0t.SetLabel([10,8])
 Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])

 A_t=uni10.Contract(T1,T0)
 A_t=uni10.Contract(A_t, Iso0)
 A_t=uni10.Contract(A_t, Iso0t)
 A_t=uni10.Contract(A_t, Iso1t)
 A_t=uni10.Contract(A_t, U)
 A_t=uni10.Contract(A_t, t0)
 Norm0=uni10.Contract(A_t, t1)

 #Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1t)))*(U))*(t0*t1)
 Norm0=uni10.Permute(Norm0,[5,9],1)

 t0.SetLabel([1,10,4,14])
 t1.SetLabel([14,11,7,6])

 A_t=uni10.Contract(T1,T0)
 A_t=uni10.Contract(A_t, Iso0)
 A_t=uni10.Contract(A_t, Iso0t)
 A_t=uni10.Contract(A_t, Iso1t)
 A_t=uni10.Contract(A_t, t0)
 Norm1=uni10.Contract(A_t, t1)



# Norm1=((((T1*T0)*(Iso0*Iso0t))*((Iso1t))))*(t0*t1)
 Norm1=uni10.Permute(Norm1,[5,9],1)

 return Norm1+(-2)*Norm0

def  Env_uall(T0,T1,U,Iso0,Iso1):

 t0=T0*1.0
 t1=T1*1.0
 Iso0t=Iso0*1.0
 Iso1t=Iso1*1.0

 T0.SetLabel([1,2,4,3])
 T1.SetLabel([3,5,7,6])

 t0.SetLabel([1,12,4,14])
 t1.SetLabel([14,13,7,6])

 Iso0.SetLabel([2,8])
 Iso1.SetLabel([5,9])
 Iso0t.SetLabel([10,8])
 Iso1t.SetLabel([11,9])
 U.SetLabel([12,13,10,11])


 A_t=uni10.Contract(T1,T0)
 A_t=uni10.Contract(A_t, Iso0t)
 A_t=uni10.Contract(A_t, Iso0)
 A_t=uni10.Contract(A_t, Iso1t)
 A_t=uni10.Contract(A_t, Iso1)
 A_t=uni10.Contract(A_t, t0)
 Norm0=uni10.Contract(A_t, t1)



 #Norm0=((((T1*T0)*(Iso0*Iso0t))*((Iso1*Iso1t))))*(t0*t1)
 Norm0=uni10.Permute(Norm0,[12,13,10,11],2)

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
   svd=Y.GetBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso0.PutBlock(temporary_matrix)
   Iso0.SetLabel([1,2])
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
   svd=Y.GetBlock().svd()
   temporary_matrix=(-1)*svd[0]*svd[2]
   Iso1.PutBlock(temporary_matrix)
   Iso1.SetLabel([1,2])
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
   svd=Y.GetBlock().svd()
   temporary_matrix=svd[0]*svd[2]
   U.PutBlock(temporary_matrix)
   U.SetLabel([1,2,3,4])
   E2=copy.copy(E1)

 return Iso0,Iso1,U








def    optimized_All(T0,T1,U,Iso0,Iso1):

 E2=10
 for i in xrange(10):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   #print "1", i, E1
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Iso0", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso0_update=Iso0*1.0
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso0(T0,T1,U,Iso0,Iso1)
   #print Y.printDiagram() 
   #svd=Y.GetBlock().svd()
   svd=linalg.svd(Y.GetBlock(), full_matrices=False, lapack_driver='gesvd')
   temporary_matrix=np.matmul(-1.0*svd[0], svd[2])



   #temporary_matrix=(-1)*svd[0]*svd[2]
   Iso0.PutBlock(temporary_matrix)
   Iso0.SetLabel([1,2])
   E2=E1*1.0

 E2=10
 for i in xrange(10):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   #print "2", i, E1
   if E1< 1.0e-14:
     break
   #print "Iso1", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    Iso1_update=Iso1*1.0
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
     #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
     fidel=E2
   Y=Env_Iso1(T0,T1,U,Iso0,Iso1) 
   svd=linalg.svd(Y.GetBlock(), full_matrices=False, lapack_driver='gesvd')
   temporary_matrix=np.matmul(-1.0*svd[0], svd[2])

   #svd=Y.GetBlock().svd()
   #temporary_matrix=(-1)*svd[0]*svd[2]
   Iso1.PutBlock(temporary_matrix)
   Iso1.SetLabel([1,2])
   E2=E1*1.0


 E2=10
 fidel=0
 for i in xrange(3):
   E1=Cost_FunctionAll(T0,T1,U,Iso0,Iso1)
   #print "3", i, E1
   if E1< 1.0e-14:
     #print "E1< 1.0e-14"
     break
   #print "Uni", i, E1, E2, abs((E2-E1)/E1)
   fidel=E1
   if E1<E2 or i is 0:
    U_update=U*1.0
    if abs((E2-E1)/E1) < 1.0e-7:
     #print E2, E1, abs((E2-E1)/E1), i
     break
   else:
    #print 'Notoptimized=i, E1, E2=', i,'  ', E1, '   ', E2
    fidel=E2
   Y=Env_uall(T0,T1,U,Iso0,Iso1)
#   svd=Y.GetBlock().svd()
#   temporary_matrix=svd[0]*svd[2]
   svd=linalg.svd(Y.GetBlock(), full_matrices=False, lapack_driver='gesvd')
   temporary_matrix=np.matmul(svd[0], svd[2])
   U.PutBlock(temporary_matrix)
   U.SetLabel([1,2,3,4])
   E2=copy.copy(E1)

 return Iso0,Iso1,U




def Simple_update_root(mps_A, mps_R, N, bdi, bdip, bdo, bdop, Uni_list):

 for i in xrange(N/2):

  t0=uni10.UniTensorR([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
  t0.PutBlock(mps_R[2*i].GetBlock())

  t1=uni10.UniTensorR([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
  t1.PutBlock(mps_R[2*i+1].GetBlock())


  T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.PutBlock(mps_A[2*i].GetBlock())

  T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.PutBlock(mps_A[2*i+1].GetBlock())

  U=optimized_root(T0,T1,t0,t1, Uni_list[i*2],2*i)
  Uni_list[i*2].PutBlock(U.GetBlock())

 return  Uni_list


#@profile
def   Simple_update(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list):

 #bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 #bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 N_A=mps_A.norm()
 Iso_list=[None]*N
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensorR([bdip,bdo])
  Iso_list[i].Identity()
 fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 #print "fidel_val", fidel_val
 count_list.append(0)
 Fidel_list.append(fidel_val)
 for i in xrange(N):
  #print "Step", i
  T0=uni10.UniTensorR([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
  T0.PutBlock(mps_A[i].GetBlock())
  T0.SetLabel([0,1,2,3])
  U=Iso_list[i]*1.0
  #U, s, V=svd_parity(Tempo)
  U.SetLabel([1,3])
  U_transpose=U*1.0
  U_transpose.SetLabel([4,3])
  U=optimized_U(T0, U, U_transpose)
  Iso_list[i].PutBlock(U.GetBlock())

  fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #print   "fidel_valIso", fidel_val
  count_list.append(count_list[-1]+1)
  Fidel_list.append(fidel_val)

 #print  count_list, Fidel_list
 
 
 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
  t1.Randomize(dn_mu=-1,up_var=1,seed=1)
  #svd=t1.GetBlock().svd()
  svd=linalg.svd(t1.GetBlock(), full_matrices=False, lapack_driver='gesvd')
  t1.PutBlock(svd[0])
  t1.Identity()
  Uni_list.append(t1*1.0)

 Uni_list_iden=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
  t1.Identity()
  Uni_list_iden.append(t1*1.0)


 for q in xrange(2):
  for i in xrange(N/2):
   #print "Step", i
   T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T0.PutBlock(mps_A[2*i].GetBlock())
   T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T1.PutBlock(mps_A[2*i+1].GetBlock())

   Iso0,Iso1,U=optimized_All(T0,T1,Uni_list[i*2],Iso_list[2*i],Iso_list[2*i+1])
   Iso_list[i*2].PutBlock(Iso0.GetBlock())
   Iso_list[i*2+1].PutBlock(Iso1.GetBlock())
   Uni_list[i*2].PutBlock(U.GetBlock())

   Iso_list_copy=[None]*(N)
#    for  i  in  xrange(N):
#     Iso_list_copy[i]=Iso_list[i]*1.0
#    Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
#    Env_left[0].SetLabel([1,2,3,4])
#    Env_right[1].SetLabel([1,2,3,4])
#    N_val=uni10.Contract(Env_left[0],Env_right[1])
#  
#    Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
#    Env_left[0].SetLabel([1,2,3,4])
#    Env_right[1].SetLabel([1,2,3,4])
#    D_val=uni10.Contract(Env_left[0],Env_right[1])
#  
#    #print "fidel_valUni", abs(N_val.GetBlock()[0,0])/(abs(D_val.GetBlock()[0,0]*N_A)**(0.5))
#    fidel_val=abs(N_val.GetBlock()[0,0])/(abs(D_val.GetBlock()[0,0]*N_A)**(0.5))
# 
#    count_list.append(count_list[-1]+1)
#    Fidel_list.append(fidel_val)


 return  Iso_list, Uni_list

#@profile
def   Simple_trivial(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list):

 #bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 #bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 Iso_list=[None]*(N)
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensorR([bdip,bdo])
  Iso_list[i].Identity()
# Iso_list[0].SetLabel([1,2])
# A=copy.copy(Iso_list[0])
# A.SetLabel([-1,2])
# Result=A*Iso_list[0]
# Result.uni10.Permute([1,-1],1)
 
 #print Iso_list[0].printDiagram(), Result

 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
  t1.Identity()
  Uni_list.append(t1*1.0)

 return  Iso_list, Uni_list










def Simple_update_chi(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, chi_canon, Fidel_list, count_list):

 bdiChi=uni10.Bond(uni10.BD_IN, chi_canon)
 bdoChi=uni10.Bond(uni10.BD_OUT, chi_canon)

 Iso_list=[None]*(N)
 for i in xrange(N):
  Iso_list[i]=uni10.UniTensorR([bdip,bdoChi])
  Iso_list[i].Identity()
 
 fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 
 
 Uni_list=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensorR([bdiChi,bdiChi,bdoChi,bdoChi])
  t1.Randomize(dn_mu=-1,up_var=1,seed=1)
  svd=t1.GetBlock().svd()
  t1.PutBlock(svd[0])
  t1.Identity()
  Uni_list.append(copy.copy(t1))

 Uni_list_iden=[]
 for i in xrange(len(Iso_list)-1):
  t1=uni10.UniTensorR([bdiChi,bdiChi,bdoChi,bdoChi])
  t1.Identity()
  Uni_list_iden.append(copy.copy(t1))



 for i in xrange(N/2):
  #print "Step", i
  T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.PutBlock(mps_A[2*i].GetBlock())
  #T0.SetLabel([0,1,2,3])
  T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.PutBlock(mps_A[2*i+1].GetBlock())
  #T1.SetLabel([0,1,2,3])

  Iso0,Iso1,U=optimized_All_chi(T0,T1,Uni_list[i*2],Iso_list[2*i],Iso_list[2*i+1])
  Iso_list[i*2].PutBlock(Iso0.GetBlock())
  Iso_list[i*2+1].PutBlock(Iso1.GetBlock())
  Uni_list[i*2].PutBlock(U.GetBlock())

  Iso_list_copy=[None]*(N)
  for  i  in  xrange(N):
   Iso_list_copy[i]=copy.copy(Iso_list[i])
  Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_left[0].SetLabel([1,2,3,4])
  Env_right[1].SetLabel([1,2,3,4])
  N_val=Env_left[0]*Env_right[1]

  Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_left[0].SetLabel([1,2,3,4])
  Env_right[1].SetLabel([1,2,3,4])
  D_val=Env_left[0]*Env_right[1]

  #print   "fidel_valUni", abs(N_val[0])/(abs(D_val[0])**(0.5))
  count_list.append(count_list[-1]+1)
  Fidel_list.append(fidel_val)




 return  Iso_list, Uni_list


def Fidel_Iso(mps_A, Iso_list, N,bdi,bdi1,bdip,bdo,bdo1,bdop):
 N_A=mps_A.norm()
 Env_left=[None]*(N)
 for i in xrange( N):
   #print i, 2*i+1
   if i == 0:
      #print "i=inside", i, bdip,bdi,mps_A[i].printDiagram()
      t0=uni10.UniTensorR([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
      t0.PutBlock(mps_A[i].GetBlock())
      t0.SetLabel([0,1,2,-6])
      T0=t0*1.0
      T0.SetLabel([0,4,2,-7])
      Iso=Iso_list[i]
      Iso.SetLabel([1,-3])
      Iso_tran=Iso*1.0
      Iso_tran.SetLabel([4,-3])
      A_t=uni10.Contract(Iso,Iso_tran)
      B_t=uni10.Contract(t0,T0)
      Result=uni10.Contract(A_t,B_t)
      Result=uni10.Permute(Result,[-6,-7],2)
      Env_left[i]=Result
   else:
      #print "i=inside", i
      t0=uni10.UniTensorR([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
      t0.PutBlock(mps_A[i].GetBlock())
      t0.SetLabel([-8,1,2,-6])
      T0=t0*1.0
      T0.SetLabel([-9,4,2,-7])
      Iso=Iso_list[i]
      Iso.SetLabel([1,-3])
      Iso_tran=Iso*1.0
      Iso_tran.SetLabel([4,-3])
      Env_left[i-1].SetLabel([-8,-9])

      A_t=uni10.Contract(t0,Env_left[i-1])
      B_t=uni10.Contract(A_t,T0)
      A_t=uni10.Contract(Iso,Iso_tran)
      Result=uni10.Contract(A_t,B_t)
      Result=uni10.Permute(Result,[-6,-7],2)
      Env_left[i]=Result

 #print Env_left[(N/2)-1][0]
 ret = Env_left[N-1].GetRawElem()
 assert ret.shape == (1,1)

 return ((abs(ret[0,0]))**(0.5))/ (N_A**(0.5))



def Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop):
 Env_right=[None]*(N/2)
 for i in xrange((N/2)-1, 0, -1):
   #print i, 2*i+1
   if i == ((N/2)-1):
      #print "i=initRight", i
      t0=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.PutBlock(mps_A[2*i+1].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.PutBlock(mps_A[2*i].GetBlock())
      t1.SetLabel([7,5,6,1])
      #print mps_R[2*i+1].printDiagram()
      T0=t0*1.0
      T0.SetLabel([18,13,4,3])
      T1=t1*1.0
      T1.SetLabel([17,14,6,18])

      Iso=Iso_list[2*i+1]
      Iso.SetLabel([2,8])
      Iso_t=Isop_list[2*i+1]*1.0
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.SetLabel([5,9])

      IsoP_t=Isop_list[2*i]*1.0
      IsoP_t.SetLabel([10,9])

      U=Uni_list[2*i]
      U.SetLabel([12,13,10,11])

      Up=Uni_list[2*i-1]
      Up.SetLabel([15,14,16,12])

      A_t=uni10.Contract(IsoP,IsoP_t)
      A_t=uni10.Contract(A_t, t1)
      B_t=uni10.Contract(Iso,Iso_t)
      B_t=uni10.Contract(B_t, t0)
      A_t=uni10.Contract(A_t, B_t)
      A_t=uni10.Contract(A_t, U)
      A_t=uni10.Contract(A_t, Up)
      A_t=uni10.Contract(A_t, T1)
      Result=uni10.Contract(A_t, T0)

      #Result=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Up))*(T1*T0)
      Result=uni10.Permute(Result,[7,16,15,17],4)
      #print Result.printDiagram()
      Env_right[i]=Result
   else:
      #print "i=insideRight", i
      #print "i=inside", i
      t0=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.PutBlock(mps_A[2*i+1].GetBlock())
      t0.SetLabel([1,2,4,-7])
      t1=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.PutBlock(mps_A[2*i].GetBlock())
      t1.SetLabel([7,5,6,1])
      #print mps_R[2*i+1].printDiagram()
      T0=t0*1.0
      T0.SetLabel([18,-15,4,-17])
      T1=t1*1.0
      T1.SetLabel([17,14,6,18])
      Iso=Iso_list[2*i+1]
      Iso.SetLabel([2,8])
      Iso_t=Isop_list[2*i+1]*1.0
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.SetLabel([5,9])

      IsoP_t=Isop_list[2*i]*1.0
      IsoP_t.SetLabel([10,9])

      U=Uni_list[2*i]
      U.SetLabel([12,-16,10,11])

      Up=Uni_list[2*i-1]
      Up.SetLabel([15,14,16,12])

      Env_right[i+1].SetLabel([-7,-16,-15,-17])
      Result=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(t1,uni10.Contract(IsoP,IsoP_t)),(uni10.Contract(t0,uni10.Contract(Iso,Iso_t)))),uni10.Contract(U,Up)),uni10.Contract(T1,T0)),Env_right[i+1])
      Result=uni10.Permute(Result,[7,16,15,17],4)
      #print Result.printDiagram()
      Env_right[i]=Result

 #print "\n"
 Env_left=[None]*(N/2)
 for i in xrange( 0, (N/2)-1):
   #print i, 2*i+1
   if i == 0:
      #print "i=initLeft", i
      t0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t0.PutBlock(mps_A[2*i].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t1.PutBlock(mps_A[2*i+1].GetBlock())
      t1.SetLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=t0*1.0
      T0.SetLabel([1,12,4,18])
      T1=t1*1.0
      T1.SetLabel([18,19,6,20])

      Iso=Iso_list[2*i]
      Iso.SetLabel([2,8])
      Iso_t=Isop_list[2*i]*1.0
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*i+1]
      IsoP.SetLabel([5,9])

      IsoP_t=Isop_list[2*i+1]*1.0
      IsoP_t.SetLabel([10,9])

      U=Uni_list[2*i]
      U.SetLabel([12,13,11,10])

      A_t=uni10.Contract(IsoP,IsoP_t)
      A_t=uni10.Contract(A_t, t1)
      B_t=uni10.Contract(Iso,Iso_t)
      B_t=uni10.Contract(B_t, t0)
      A_t=uni10.Contract(A_t, B_t)
      A_t=uni10.Contract(A_t, U)
      A_t=uni10.Contract(A_t, T1)
      Result=uni10.Contract(A_t, T0)


      #Result=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*(T1*T0)
      Result=uni10.Permute(Result,[7,13,19,20],0)
      #print Result.printDiagram()
      Env_left[i]=Result
   else:
      #print "i=insideLeft", i
      t0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t0.PutBlock(mps_A[2*i].GetBlock())
      t0.SetLabel([-7,2,4,3])
      t1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t1.PutBlock(mps_A[2*i+1].GetBlock())
      t1.SetLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=t0*1.0
      T0.SetLabel([-20,12,4,18])
      T1=t1*1.0
      T1.SetLabel([18,19,6,20])
      Iso=Iso_list[2*i]
      Iso.SetLabel([2,8])
      Iso_t=Isop_list[2*i]*1.0
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*i+1]
      IsoP.SetLabel([5,9])

      IsoP_t=Isop_list[2*i+1]*1.0
      IsoP_t.SetLabel([10,9])

      U=Uni_list[2*i]
      U.SetLabel([21,13,11,10])

      Up=Uni_list[2*i-1]
      Up.SetLabel([-19,12,-13,21])


      Env_left[i-1].SetLabel([-7,-13,-19,-20])

      Result=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(t1,uni10.Contract(IsoP,IsoP_t)),(uni10.Contract(t0,uni10.Contract(Iso,Iso_t)))),U),uni10.Contract(uni10.Contract(T1,T0),Up)),Env_left[i-1])
      Result=uni10.Permute(Result,[7,13,19,20],4)
      #print Result.printDiagram()
      Env_left[i]=Result

 return Env_left, Env_right

def  Obtain_Env_U(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_u=0

 if Location==0:
      #print "Location=UniFirst", Location
      t0=uni10.UniTensorR([mps_A[2*Location].bond(0),bdip,bdi,mps_A[2*Location].bond(2)])
      t0.PutBlock(mps_A[2*Location].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[2*Location+1].bond(0),bdip,bdi,mps_A[2*Location+1].bond(2)])
      t1.PutBlock(mps_A[2*Location+1].GetBlock())
      t1.SetLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.SetLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.SetLabel([18,19,6,20])

      Iso=Iso_list[2*Location]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*Location])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*Location+1]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*Location+1])
      IsoP_t.SetLabel([10,9])

      #U=Uni_list[2*Location]
      #U.SetLabel([12,13,11,10])


      Env_right[Location+1].SetLabel([7,13,19,20])
      #print Result.printDiagram()
      Env_u=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Env_right[Location+1]))*(T1*T0)
      #print Env_u.printDiagram()
      Env_u.uni10.Permute([12,13,11,10],2)
 elif Location==(N-2) or Location==(N-3):
      #print "Location=UniEnd", Location
      if Location%2==0:
       i=Location/2
      else:
       i=(Location/2)+1
      t0=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      t0.PutBlock(mps_A[2*i+1].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      t1.PutBlock(mps_A[2*i].GetBlock())
      t1.SetLabel([7,5,6,1])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.SetLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.SetLabel([17,14,6,18])

      Iso=Iso_list[2*i+1]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[2*i+1])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[2*i]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[2*i])
      IsoP_t.SetLabel([10,9])

      U=Uni_list[2*i]
      U.SetLabel([12,13,10,11])

      Up=Uni_list[2*i-1]
      Up.SetLabel([15,14,16,12])

      Env_left[i-1].SetLabel([7,16,15,17])
      if Location%2==0:
       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Up))*(Env_left[i-1]))*(T1*T0)
       Env_u.uni10.Permute([12,13,10,11],2)
       #print Result.printDiagram()
      else:
       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*(Env_left[i-1]))*(T1*T0)
       Env_u.uni10.Permute([15,14,16,12],2)
       #print Result.printDiagram()
 else:
      if (Location%2)==0:
       #print "Location=uniMidel", Location
       i=(Location/2)
       t0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
       t0.PutBlock(mps_A[2*i].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
       t1.PutBlock(mps_A[2*i+1].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])
       Iso=Iso_list[2*i]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[2*i])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[2*i+1]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[2*i+1])
       IsoP_t.SetLabel([10,9])

       #U=Uni_list[2*i]
       #U.SetLabel([21,13,11,10])

       Up=Uni_list[2*i-1]
       Up.SetLabel([-19,12,-13,21])


       Env_left[i-1].SetLabel([-7,-13,-19,-20])
       Env_right[i+1].SetLabel([7,13,19,20])

       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(Env_right[i+1]))*(((T1*T0)*Env_left[i-1])*Up))
       Env_u.uni10.Permute([21,13,11,10],2)
       #print Result.printDiagram()
       #Env_u=Result
      else:
       #print "Location=uniMidelodd", Location
       i=((Location+1)/2)
       t0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
       t0.PutBlock(mps_A[2*i].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
       t1.PutBlock(mps_A[2*i+1].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])
       Iso=Iso_list[2*i]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[2*i])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[2*i+1]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[2*i+1])
       IsoP_t.SetLabel([10,9])

       U=Uni_list[2*i]
       U.SetLabel([21,13,11,10])

       #Up=Uni_list[2*i-1]
       #Up.SetLabel([-19,12,-13,21])


       Env_left[i-1].SetLabel([-7,-13,-19,-20])
       Env_right[i+1].SetLabel([7,13,19,20])

       Env_u=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Env_right[i+1]))*((T1*T0)*Env_left[i-1]))
       Env_u.uni10.Permute([-19,12,-13,21],2)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_u


def  Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_Iso=0

 if Location==0 or Location==1:
      t0=uni10.UniTensorR([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
      t0.PutBlock(mps_A[0].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
      t1.PutBlock(mps_A[1].GetBlock())
      t1.SetLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.SetLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.SetLabel([18,19,6,20])

      Iso=Iso_list[0]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[0])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[1]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[1])
      IsoP_t.SetLabel([10,9])

      U=Uni_list[0]
      U.SetLabel([12,13,11,10])
      Env_right[1].SetLabel([7,13,19,20])
      if Location%2==0:
        Env_Iso=(((t1*(IsoP*IsoP_t))*((t0*(Iso_t)))*U)*(Env_right[1]))*(T1*T0)
        #print Env_u.printDiagram()
        Env_Iso.uni10.Permute([2,8],1)
      else:
        Env_Iso=(((t1*(IsoP_t))*((t0*(Iso*Iso_t)))*U)*(Env_right[1]))*(T1*T0)
        #print Env_u.printDiagram()
        Env_Iso.uni10.Permute([5,9],1)


 elif Location==(N-1) or Location==(N-2):

      #print "Location=IsoEnd", Location
      t0=uni10.UniTensorR([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
      t0.PutBlock(mps_A[N-1].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
      t1.PutBlock(mps_A[N-2].GetBlock())
      t1.SetLabel([7,5,6,1])

      T0=copy.copy(t0)
      T0.SetLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.SetLabel([17,14,6,18])

      Iso=Iso_list[N-1]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[N-1])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[N-2]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[N-2])
      IsoP_t.SetLabel([10,9])

      U=Uni_list[N-2]
      U.SetLabel([12,13,10,11])

      Up=Uni_list[N-3]
      Up.SetLabel([15,14,16,12])

      Env_left[((N-2)/2)-1].SetLabel([7,16,15,17])
      if Location==(N-1):
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso_t))))*(U*Up))*(Env_left[((N-2)/2)-1]))*(T1*T0)
       Env_Iso.uni10.Permute([2,8],1)
      else:
       Env_Iso=((((t1*(IsoP_t))*((t0*(Iso*Iso_t))))*(Up*U))*(Env_left[((N-2)/2)-1]))*(T1*T0)
       Env_Iso.uni10.Permute([5,9],1)
 else:
      if (Location%2)==0:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensorR([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t0.PutBlock(mps_A[Location].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[Location+1].bond(0),bdip,bdi,mps_A[Location+1].bond(2)])
       t1.PutBlock(mps_A[Location+1].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])

       Iso=Iso_list[Location]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[Location+1]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location+1])
       IsoP_t.SetLabel([10,9])

       U=Uni_list[Location]
       U.SetLabel([21,13,11,10])

       Up=Uni_list[Location-1]
       Up.SetLabel([-19,12,-13,21])

       #print (Location/2)-1, (Location/2)+1
       Env_left[(Location/2)-1].SetLabel([-7,-13,-19,-20])
       Env_right[(Location/2)+1].SetLabel([7,13,19,20])
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso_t)))*(U*Up))*(Env_right[(Location/2)+1]))*(((T1*T0)*Env_left[(Location/2)-1])))
       Env_Iso.uni10.Permute([2,8],1)


      else:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensorR([mps_A[Location-1].bond(0),bdip,bdi,mps_A[Location-1].bond(2)])
       t0.PutBlock(mps_A[Location-1].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t1.PutBlock(mps_A[Location].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])
       Iso=Iso_list[Location-1]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location-1])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[Location]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location])
       IsoP_t.SetLabel([10,9])

       U=Uni_list[Location-1]
       U.SetLabel([21,13,11,10])

       Up=Uni_list[Location-2]
       Up.SetLabel([-19,12,-13,21])


       Env_left[((Location-1)/2)-1].SetLabel([-7,-13,-19,-20])
       Env_right[(Location+1)/2].SetLabel([7,13,19,20])

       Env_Iso=((((t1*(IsoP_t))*((t0*(Iso*Iso_t)))*(U*Up))*(Env_right[(Location+1)/2]))*((T1*T0)*Env_left[((Location-1)/2)-1]))
       Env_Iso.uni10.Permute([5,9],1)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_Iso





def  Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
 Env_Iso=0

 if Location==0 or Location==1:
      t0=uni10.UniTensorR([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
      t0.PutBlock(mps_A[0].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
      t1.PutBlock(mps_A[1].GetBlock())
      t1.SetLabel([3,5,6,7])
      #print mps_R[2*i+1].printDiagram()
      T0=copy.copy(t0)
      T0.SetLabel([1,12,4,18])
      T1=copy.copy(t1)
      T1.SetLabel([18,19,6,20])

      Iso=Iso_list[0]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[0])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[1]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[1])
      IsoP_t.SetLabel([10,9])

      U=Uni_list[0]
      U.SetLabel([12,13,11,10])
      Env_right[1].SetLabel([7,13,19,20])
      if Location%2==0:
        Env_Iso=(((t1*(IsoP*IsoP_t))*((t0*(Iso)))*U)*(T1*T0))*(Env_right[1])
        #print Env_u.printDiagram()
        Env_Iso.uni10.Permute([11,8],1)
      else:
        Env_Iso=(((t1*(IsoP))*((t0*(Iso*Iso_t)))*U)*(T1*T0))*(Env_right[1])
        #print Env_u.printDiagram()
        Env_Iso.uni10.Permute([10,9],1)


 elif Location==(N-1) or Location==(N-2):

      #print "Location=IsoEnd", Location
      t0=uni10.UniTensorR([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
      t0.PutBlock(mps_A[N-1].GetBlock())
      t0.SetLabel([1,2,4,3])
      t1=uni10.UniTensorR([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
      t1.PutBlock(mps_A[N-2].GetBlock())
      t1.SetLabel([7,5,6,1])

      T0=copy.copy(t0)
      T0.SetLabel([18,13,4,3])
      T1=copy.copy(t1)
      T1.SetLabel([17,14,6,18])

      Iso=Iso_list[N-1]
      Iso.SetLabel([2,8])
      Iso_t=copy.copy(Isop_list[N-1])
      Iso_t.SetLabel([11,8])

      IsoP=Iso_list[N-2]
      IsoP.SetLabel([5,9])

      IsoP_t=copy.copy(Isop_list[N-2])
      IsoP_t.SetLabel([10,9])

      U=Uni_list[N-2]
      U.SetLabel([12,13,10,11])

      Up=Uni_list[N-3]
      Up.SetLabel([15,14,16,12])

      Env_left[((N-2)/2)-1].SetLabel([7,16,15,17])
      if Location==(N-1):
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso))))*(U*Up))*(T1*T0))*(Env_left[((N-2)/2)-1])
       Env_Iso.uni10.Permute([11,8],1)
      else:
       Env_Iso=((((t1*(IsoP))*((t0*(Iso*Iso_t))))*(Up*U))*(T1*T0))*(Env_left[((N-2)/2)-1])
       Env_Iso.uni10.Permute([10,9],1)
 else:
      if (Location%2)==0:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensorR([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t0.PutBlock(mps_A[Location].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[Location+1].bond(0),bdip,bdi,mps_A[Location+1].bond(2)])
       t1.PutBlock(mps_A[Location+1].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])

       Iso=Iso_list[Location]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[Location+1]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location+1])
       IsoP_t.SetLabel([10,9])

       U=Uni_list[Location]
       U.SetLabel([21,13,11,10])

       Up=Uni_list[Location-1]
       Up.SetLabel([-19,12,-13,21])

       #print (Location/2)-1, (Location/2)+1
       Env_left[(Location/2)-1].SetLabel([-7,-13,-19,-20])
       Env_right[(Location/2)+1].SetLabel([7,13,19,20])
       Env_Iso=((((t1*(IsoP*IsoP_t))*((t0*(Iso)))*(U*Up))*(T1*T0))*(((Env_right[(Location/2)+1])*Env_left[(Location/2)-1])))
       Env_Iso.uni10.Permute([11,8],1)


      else:
       #print "Location=IsoMidelodd", Location
       t0=uni10.UniTensorR([mps_A[Location-1].bond(0),bdip,bdi,mps_A[Location-1].bond(2)])
       t0.PutBlock(mps_A[Location-1].GetBlock())
       t0.SetLabel([-7,2,4,3])
       t1=uni10.UniTensorR([mps_A[Location].bond(0),bdip,bdi,mps_A[Location].bond(2)])
       t1.PutBlock(mps_A[Location].GetBlock())
       t1.SetLabel([3,5,6,7])
       #print mps_R[2*i+1].printDiagram()
       T0=copy.copy(t0)
       T0.SetLabel([-20,12,4,18])
       T1=copy.copy(t1)
       T1.SetLabel([18,19,6,20])
       Iso=Iso_list[Location-1]
       Iso.SetLabel([2,8])
       Iso_t=copy.copy(Isop_list[Location-1])
       Iso_t.SetLabel([11,8])

       IsoP=Iso_list[Location]
       IsoP.SetLabel([5,9])

       IsoP_t=copy.copy(Isop_list[Location])
       IsoP_t.SetLabel([10,9])

       U=Uni_list[Location-1]
       U.SetLabel([21,13,11,10])

       Up=Uni_list[Location-2]
       Up.SetLabel([-19,12,-13,21])


       Env_left[((Location-1)/2)-1].SetLabel([-7,-13,-19,-20])
       Env_right[(Location+1)/2].SetLabel([7,13,19,20])

       Env_Iso=((((t1*(IsoP))*((t0*(Iso*Iso_t)))*(U*Up))*(T1*T0))*((Env_right[(Location+1)/2])*Env_left[((Location-1)/2)-1]))
       Env_Iso.uni10.Permute([10,9],1)
       #print Result.printDiagram()
       #Env_u=Result

 return Env_Iso


def optimize_uni_locally(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right):
  #print Uni_list[0]
  E2=0
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].Identity()   

  Iso_list_copy=[  copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]


  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)

  Env_left1[0].SetLabel([1,2,3,4])
  Env_right1[1].SetLabel([1,2,3,4])
  Dominator=Env_left1[0]*Env_right1[1]

  Dominator=(abs(Dominator[0]))**(0.5)
  #print Uni_list[0]
  #print  "Dominator", Dominator 
  Fidel=0
  
  for i in xrange(4):
    Env_u=Obtain_Env_U(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_u.SetLabel([1,2,3,4])
    Uni_list[Location].SetLabel([1,2,3,4])
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
     Uni_list[Location].PutBlock(U_update.GetBlock())
     Fidel=E2
     break
    #Y=Env_U(T0, T1, U, U_transpose)
    #print Y.printDiagram() 
    svd=Env_u.GetBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    Uni_list[Location].PutBlock(temporary_matrix)
    E2=copy.copy(E1)
  #print "Location", Location, "FidelUni", Fidel
  return Fidel

def optimize_uni_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right):
  #print "\n", "Location", Location
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].Identity()   
  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)

  Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
  Env_u.SetLabel([1,2,3,4,5,6,7,8])
  Uni_list[Location].SetLabel([1,2,3,4,5,6,7,8])
  R=Uni_list[Location]*Env_u

  Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
  Env_iso1.SetLabel([1,2,3,4,5,6])
  Iso_list[Location].SetLabel([1,2,3,4,5,6])
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
   Env_u.SetLabel([1,2,3,4,5,6,7,8])
   Uni_list[Location].SetLabel([1,2,3,4,5,6,7,8])
   R=Uni_list[Location]*Env_u
   #E1=R[0]/Dominator
   E1=1+(-2)*R[0]+Dominator[0]
   Fidel=R[0]/(Dominator[0]**(0.5))

   #print 'E=', E1,  Fidel

   #E_list.append(E1)
   #Count_list.append(count)
   #print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().GetBlock()*U_update).trace()
   D_u=(-2.0)*Env_u
   D_u_mat=copy.copy(D_u.GetBlock())
   D_u_trans=copy.copy(D_u.GetBlock())
   D_u_trans.transpose()
   Uni=copy.copy(Uni_list[Location].GetBlock())
   #Iso_t.transpose()
   Uni_t=copy.copy(Uni)
   Uni_t.transpose()
   Z_decent_mat=Uni*D_u_trans*Uni+(-1.00)*D_u_mat
   Z_decent=copy.copy(D_u)
   Z_decent.PutBlock(Z_decent_mat)
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
    Uni_list[Location].PutBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_u.SetLabel([1,2,3,4,5,6,7,8])
    Uni_list[Location].SetLabel([1,2,3,4,5,6,7,8])
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
    Uni_list[Location].PutBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_u=Obtain_Env_U(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_u.SetLabel([1,2,3,4,5,6,7,8])
    Uni_list[Location].SetLabel([1,2,3,4,5,6,7,8])
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
   Uni_list[Location].PutBlock(Temporary)

  #print "Location", Location, "FidelUni", Fidel






def optimize_iso_locally_SD(mps_A, Iso_list,Isop_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right):
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
     Uni_list_iden[i].Identity()

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]

  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)

  Env_iso=Obtain_Env_Iso(mps_A,Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  R=Iso_list[Location]*Env_iso

  Env_iso1=Obtain_Env_Iso(mps_A,Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso1.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
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
   Env_iso.SetLabel([1,2])
   Iso_list[Location].SetLabel([1,2])
   R=Iso_list[Location]*Env_iso

   Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
   Env_iso1.SetLabel([1,2])
   Iso_list[Location].SetLabel([1,2])
   Dominator=Iso_list[Location]*Env_iso1
   Dominator=(abs(Dominator[0]))**(0.5)
   #E1=R[0]/Dominator
   E1=1+(-2)*R[0]+Dominator[0]
   Fidel=R[0]/(Dominator[0]**(0.5))
   #print 'E=', E1 , Fidel 
   #E_list.append(E1)
   #Count_list.append(count)
   #print 'test', (Env_Uni_inner[L_position][L_lay_selected].transpose().GetBlock()*U_update).trace()
   D_u=(-2.0)*Env_iso+Env_iso1
   D_u_mat=copy.copy(D_u.GetBlock())
   D_u_trans=copy.copy(D_u.GetBlock())
   D_u_trans.transpose()
   Iso=copy.copy(Iso_list[Location].GetBlock())
   #Iso_t.transpose()
   Iso_t=copy.copy(Iso)
   Iso_t.transpose()
   Z_decent_mat=Iso*D_u_trans*Iso+(-1.00)*D_u_mat
   Z_decent=copy.copy(D_u)
   Z_decent.PutBlock(Z_decent_mat)
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
    Iso_list[Location].PutBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_iso.SetLabel([1,2,3,4,5,6])
    Iso_list[Location].SetLabel([1,2,3,4,5,6])
    R=Iso_list[Location]*Env_iso
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
    Env_iso1.SetLabel([1,2,3,4,5,6])
    Iso_list[Location].SetLabel([1,2,3,4,5,6])
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
    Iso_list[Location].PutBlock(Temporary)
    #print U_list[L_position][L_lay_selected], Temporary
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    Env_iso.SetLabel([1,2,3,4,5,6])
    Iso_list[Location].SetLabel([1,2,3,4,5,6])
    R=Iso_list[Location]*Env_iso
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
    Env_iso1.SetLabel([1,2,3,4,5,6])
    Iso_list[Location].SetLabel([1,2,3,4,5,6])
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
   Iso_list[Location].PutBlock(Temporary)
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
    Env_iso.SetLabel([1,2])
    Iso_list[Location].SetLabel([1,2])
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
     Iso_list[Location].PutBlock(U_update.GetBlock())
     Iso_list_copy[Location].PutBlock(copy.copy(U_update.GetBlock()))
     Fidel=E2
     break
    Env_iso=Env_iso
    svd=Env_iso.GetBlock().svd()
    temporary_matrix=(1.0)*svd[0]*svd[2]
    Iso_list[Location].PutBlock(temporary_matrix)
    Iso_list_copy[Location].PutBlock(temporary_matrix)

    E2=copy.copy(E1)

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  R=Iso_list[Location]*Env_iso

  E1=abs(R[0])**(0.5)
  Fidel=abs(E1)
  #print "Location", Location, "FidelIso", Fidel#, R[0], Dominator
  return Fidel


def optimize_iso_locally(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):

  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
     Uni_list_iden[i].Identity()

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]

  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)
  Fidel=0
  E2=0
  for i in xrange(25):
    Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_iso.SetLabel([1,2])
    Iso_list[Location].SetLabel([1,2])
    R=Iso_list[Location]*Env_iso

    #Iso_list_copy=[  copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]
    Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
    Env_iso1.SetLabel([1,2])
    Iso_list[Location].SetLabel([1,2])
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
     Iso_list[Location].PutBlock(U_update.GetBlock())
     Iso_list_copy[Location].PutBlock(copy.copy(U_update.GetBlock()))
     Fidel=E2
     break
    Env_iso=(-2.0)*Env_iso+Env_iso1
    svd=Env_iso.GetBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Iso_list[Location].PutBlock(temporary_matrix)
    Iso_list_copy[Location].PutBlock(temporary_matrix)

    E2=copy.copy(E1)

  Iso_list_copy=[copy.copy(Iso_list[i])  for  i  in  xrange(len(Iso_list)) ]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Iso_list_copy,Uni_list_iden,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  R=Iso_list[Location]*Env_iso
  Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso1.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso1
  Dominator=(abs(Dominator[0]))**(0.5)

  E1=abs(R[0]/Dominator)
  Fidel=abs(E1)
  #print "Location", Location, "FidelIso", Fidel#, R[0], Dominator
  return  Fidel

def optimize_isop_locally(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].Identity()

  Iso_list_copy=[ copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]

  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list, Iso_list_copy, Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)
  Fidel=0
  E2=0

  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso
  Dominator=(abs(Dominator[0]))**(0.5)


  for i in xrange(5):
    Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_isop.SetLabel([1,2])
    Isop_list[Location].SetLabel([1,2])
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
     Isop_list[Location].PutBlock(U_update.GetBlock())
     Fidel=E2
     break
    Env_isop=Env_isop
    svd=Env_isop.GetBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Isop_list[Location].PutBlock(temporary_matrix)
    E2=copy.copy(E1)

  Iso_list_copy=[ copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  Dominator=Iso_list[Location]*Env_iso
  Dominator=(abs(Dominator[0]))**(0.5)
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_isop.SetLabel([1,2])
  Isop_list[Location].SetLabel([1,2])
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
#   t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
#   t1.Randomize(dn_mu=-1,up_var=1,seed=1)
#   svd=t1.GetBlock().svd()
#   t1.PutBlock(svd[0])
#   t1.Identity()
#   Uni_list.append(copy.copy(t1))

 for i in xrange(len(Iso_list)):
  Isop_list.append(copy.copy(Iso_list[i]))

#   for i in xrange(len(Uni_list)):
#    Uni_list[i].Randomize(dn_mu=-1,up_var=1,seed=1)
#   for i in xrange(len(Isop_list)):
#    Isop_list[i].Randomize(dn_mu=-1,up_var=1,seed=1)
#   for i in xrange(len(Isop_list)):
#    Iso_list[i].Randomize(dn_mu=-1,up_var=1,seed=1)

 Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 for i in xrange(len(Env_left)-1):
  Env_left[i].SetLabel([1,2,3,4])
  Env_right[i+1].SetLabel([1,2,3,4])
  R=Env_left[i]*Env_right[i+1]
  #print "Env", i, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N-1):
  Location=i
  Env_u=Obtain_Env_U(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_u.SetLabel([1,2,3,4])
  Uni_list[Location].SetLabel([1,2,3,4])
  R=Uni_list[Location]*Env_u
  #print "Env_Uni", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  R=Iso_list[Location]*Env_iso
  #print "Env_Iso", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_isop.SetLabel([1,2])
  Isop_list[Location].SetLabel([1,2])
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
  t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
  t1.Randomize(dn_mu=-1,up_var=1,seed=1)
  svd=t1.GetBlock().svd()
  t1.PutBlock(svd[0])
  t1.Identity()
  Uni_list.append(copy.copy(t1))

 for i in xrange(len(Iso_list)):
  Isop_list.append(copy.copy(Iso_list[i]))



 Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 for i in xrange(len(Env_left)-1):
  Env_left[i].SetLabel([1,2,3,4])
  Env_right[i+1].SetLabel([1,2,3,4])
  R=Env_left[i]*Env_right[i+1]
  #print "Env", i, (abs(R[0]))**(0.5)

 #print "\n"


 for i in xrange(N-1):
  Location=i
  Env_u=Obtain_Env_U(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_u.SetLabel([1,2,3,4])
  Uni_list[Location].SetLabel([1,2,3,4])
  R=Uni_list[Location]*Env_u
  #print "Env_Uni", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_iso.SetLabel([1,2])
  Iso_list[Location].SetLabel([1,2])
  R=Iso_list[Location]*Env_iso
  #print "Env_Iso", Location, (abs(R[0]))**(0.5)

 #print "\n"

 for i in xrange(N):
  Location=i
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_isop.SetLabel([1,2])
  Isop_list[Location].SetLabel([1,2])
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
  t0=uni10.UniTensorR([ mps_A[i].bond(0), bdip, bdi, mps_A[i].bond(2)])
  t0.PutBlock( mps_A[i].GetBlock() )
  t0.SetLabel([1,2,4,3])
  Iso=Iso_list[i]
  Iso.SetLabel([2,8])
  Result=uni10.Contract(Iso,t0)
  Result=uni10.Permute(Result,[1,8,4,3],3)
  Result=Result.CombineBond([8,4])
  Result=uni10.Permute(Result,[1,8,3],2)
  MPS_list[i]=Result

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(MPS_list[q].bond(2).dim())

 mps_R=MPSclass2.MPS(MPS_list[1].bond(1).dim(),max(list_bond),N_y)

 for i in xrange(N_y):
   mps_R[i]=MPS_list[i]*1.0

 return   mps_R

def  make_MPS_Q(Uni_list):

 U_even=[]
 V_even=[]

 U_odd=[]
 V_odd=[]
 
 #print len(Uni_list)-1
 for i in xrange(((len(Uni_list)-1)/2)+1):
  #print "StepsEven", i, 2*i
  #print Uni_list[0].printDiagram()
  Uni_list[2*i].SetLabel([1,2,3,4])
  Uni_0=Uni_list[2*i]*1.0
  Uni_0=uni10.Permute(Uni_0,[1,3,2,4],2)
  U, s, V=svd_parity1(Uni_0)
  U.SetLabel([1,5,-1])
  s.SetLabel([-1,-2])
  V.SetLabel([-2,2,6])
  V=uni10.Contract(V,s)
  V=uni10.Permute(V,[-1,2,6],1)
  
  U_even.append(U)
  U_even.append(V)
  
 for i in xrange((len(Uni_list)-1)/2):
  #print "StepsOdd", 2*i+1
  Uni_list[2*i+1].SetLabel([1,2,3,4])
  Uni_0=Uni_list[2*i+1]*1.0
  Uni_0=uni10.Permute(Uni_0, [1,3,2,4],2)
  U, s, V=svd_parity1(Uni_0)
  U.SetLabel([1,3,-1])
  s.SetLabel([-1,-2])
  V.SetLabel([-2,2,4])
  V=uni10.Contract(V,s)
  V=uni10.Permute(V,[-1,2,4],1)
  #print "Steps",  2*i+1, U.printDiagram()
  U_odd.append(U)
  U_odd.append(V)


 A_list=[None]*(len(Uni_list)+1)
 bdi=uni10.Bond(uni10.BD_IN, 1)
 A_temp=uni10.UniTensorR([bdi])
 A_temp.Identity()
 A_temp.SetLabel([-3])
 #print A_temp
 for i in xrange(len(Uni_list)+1):
   if i==0:
    #print "First", i
    U_even[i].SetLabel([1,5,-1])
    A=uni10.Contract(U_even[i],A_temp)
    A=uni10.Permute(A,[-3,1,5,-1],3)
    A=A.CombineBond([1,5])
    A=uni10.Permute(A,[-3,1,-1],2)
    A_list[i]=A*1.0
   elif i == (len(Uni_list)):
    #print "Last", i#, len(U_even)
    A=uni10.Contract(U_even[i],A_temp)
    A=uni10.Permute(A,[-1,2,6,-3],3)
    A=A.CombineBond([2,6])
    A=uni10.Permute(A,[-1,2,-3],2)
    A_list[i]=A*1.0
   else:
    if (i%2) != 0:
     #print "odd", i
     U_even[i].SetLabel([-1,2,3])
     U_odd[i-1].SetLabel([4,2,5])
     #print U_even[i].printDiagram(),U_odd[i-1].printDiagram() 
     A=uni10.Contract(U_even[i],U_odd[i-1])
     A=uni10.Permute(A,[-1,4,3,5],3)
     A=A.CombineBond([4,3])
     A=uni10.Permute(A,[-1,4,5],2)
     A_list[i]=A*1.0
    else:
     #print "even", i
     U_even[i].SetLabel([1,3,2])
     U_odd[i-1].SetLabel([4,5,1])
     A=uni10.Contract(U_even[i],U_odd[i-1])
     A=uni10.Permute(A,[4,5,3,2],3)
     A=A.CombineBond([5,3])
     A=uni10.Permute(A,[4,5,2],2)
     A_list[i]=A*1.0

 list_bond=[]
 for q in xrange(len(A_list)):
   list_bond.append(A_list[q].bond(2).dim())
 mps_Q=MPSclass2.MPS( A_list[1].bond(1).dim(), max(list_bond), len(A_list))
 for i in xrange(len(A_list)):
   mps_Q[i]=A_list[i]*1.0
 #print mps_Q[0].printDiagram(), mps_Q[mps_Q.N-1].printDiagram(),mps_Q.norm() 
 
 
 
#  bdip=uni10.Bond(uni10.BD_IN, 2)
#  bdiChi=uni10.Bond(uni10.BD_IN, 2)
#  t0=uni10.UniTensorR([mps_Q[0].bond(0),bdip,bdiChi,mps_Q[0].bond(2)])
#  t0.PutBlock(mps_Q[0].GetBlock())
#  t1=uni10.UniTensorR([mps_Q[1].bond(0),bdip,bdiChi,mps_Q[1].bond(2)])
#  t1.PutBlock(mps_Q[1].GetBlock())
#  t2=uni10.UniTensorR([mps_Q[2].bond(0),bdip,bdiChi,mps_Q[2].bond(2)])
#  t2.PutBlock(mps_Q[2].GetBlock())
#  t3=uni10.UniTensorR([mps_Q[3].bond(0),bdip,bdiChi,mps_Q[3].bond(2)])
#  t3.PutBlock(mps_Q[3].GetBlock())
#  
#  t0.SetLabel([0,1,-2,20])
#  t1.SetLabel([20,4,-5,21])
#  t2.SetLabel([21,7,-8,22])
#  t3.SetLabel([22,10,-11,12])
#  Resutl=t0*t1*t2*t3
# 
# 
#  t0.SetLabel([0,1,2,3])
#  t1.SetLabel([3,4,5,6])
#  t2.SetLabel([6,7,8,9])
#  t3.SetLabel([9,10,11,12])
#  Resutl1=t0*t1*t2*t3
#  A=Resutl1*Resutl
#  A.uni10.Permute([2,5,8,11,-2,-5,-8,-11],4)
#  #print A
# 
# 
# 
# 
#  t0.SetLabel([0,1,2,20])
#  t1.SetLabel([20,4,5,21])
#  t2.SetLabel([21,7,8,22])
#  t3.SetLabel([22,10,11,12])
#  Resutl=t0*t1*t2*t3
# 
# 
#  t0.SetLabel([0,-1,2,3])
#  t1.SetLabel([3,-4,5,6])
#  t2.SetLabel([6,-7,8,9])
#  t3.SetLabel([9,-10,11,12])
#  Resutl1=t0*t1*t2*t3
#  A=Resutl1*Resutl
#  A.uni10.Permute([1,4,7,10,-1,-4,-7,-10],4)
#  #print A

 return   mps_Q


def  Uni_update(Isop_list, Uni_list):

 for i in xrange(((len(Uni_list)-1)/2)+1):
  #print "StepsEven", 2*i
  #print Uni_list[0].printDiagram(), Isop_list[0].printDiagram()
  Uni_list[2*i].SetLabel([1,2,3,4])
  Isop_list[2*i].SetLabel([3,5])
  Isop_list[2*i+1].SetLabel([4,6])
  A=uni10.Contract(Uni_list[2*i], Isop_list[2*i])
  Uni_0=uni10.Contract(A,Isop_list[2*i+1])
  Uni_0=uni10.Permute(Uni_0,[1,2,5,6],2)
  Uni_list[2*i]=Uni_0*1.0

 return   Uni_list


def Env_right_R(mps_R):
 E_right=[None]*mps_R.N
 for i in reversed(xrange(mps_R.N)):
  if i == 0:
   mps_R[i].SetLabel([-1,-2,-3])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([-1,-2,-4])
   E_right[0]=uni10.Contract(uni10.Contract(mps_R_dag,mps_R[i]),E_right[1])
  elif i == (mps_R.N-1):
   mps_R[i].SetLabel([-3,-2,1])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([-4,-2,1])
   E_right[i]=uni10.Contract(mps_R_dag,mps_R[i])
   E_right[i].SetLabel([-3,-4])
  else:
   mps_R[i].SetLabel([1,-2,-3])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([2,-2,-4])
   E_right[i]=uni10.Contract(uni10.Contract(mps_R_dag,E_right[i+1]),mps_R[i])
   E_right[i]=uni10.Permute(E_right[i],[1,2],1)
   E_right[i].SetLabel([-3,-4])
 return   E_right


def Env_left_R(mps_R):
 E_left=[None]*mps_R.N
 for i in xrange(mps_R.N):
  if i == 0:
   mps_R[i].SetLabel([-1,-2,1])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([-1,-2,2])
   E_left[0]=uni10.Contract(mps_R_dag,mps_R[i])
   E_left[0]=uni10.Permute(E_left[0],[1,2],1)
   E_left[0].SetLabel([-3,-4])
  elif i == (mps_R.N-1):
   mps_R[i].SetLabel([-3,-2,1])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([-4,-2,1])
   E_left[i]=uni10.Contract(uni10.Contract(mps_R_dag,E_left[i-1]),mps_R[i])
  else:
   mps_R[i].SetLabel([-3,-2,1])
   mps_R_dag=mps_R[i]*1.0
   mps_R_dag.SetLabel([-4,-2,2])
   E_left[i]=uni10.Contract(uni10.Contract(mps_R_dag,E_left[i-1]),mps_R[i])
   E_left[i]=uni10.Permute(E_left[i],[1,2],1)
   E_left[i].SetLabel([-3,-4])
 return E_left





def Update_Env_left(mps_R,E_left,location):
 if location==0:
  mps_R[location].SetLabel([-1,-2,1])
  mps_R_dag=mps_R[location]*1.0
  mps_R_dag.SetLabel([-1,-2,2])
  E_left[0]=uni10.Contract(mps_R_dag,mps_R[location])
  E_left[0]=uni10.Permute(E_left[0],[1,2],1)
  E_left[0].SetLabel([-3,-4])
 elif location==(mps_R.N-1):
  mps_R[location].SetLabel([-3,-2,1])
  mps_R_dag=mps_R[location]*1.0
  mps_R_dag.SetLabel([-4,-2,1])
  E_left[location-1].SetLabel([-3,-4])
  E_left[location]=uni10.Contract(uni10.Contract(mps_R_dag,E_left[location-1]),mps_R[location])
 else:
  mps_R[location].SetLabel([-3,-2,1])
  mps_R_dag=mps_R[location]*1.0
  mps_R_dag.SetLabel([-4,-2,2])
  E_left[location-1].SetLabel([-3,-4])
  E_left[location]=uni10.Contract(uni10.Contract(mps_R_dag,E_left[location-1]),mps_R[location])
  E_left[location]=uni10.Permute(E_left[location],[1,2],1)
  E_left[location].SetLabel([-3,-4])
 return E_left




def Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi):
 N=mps_R.N
 Env_right=[None]*(N/2)
 for i in xrange((N/2)-1, 0, -1):
  if i == ((N/2)-1):
   #print "i=initRight", i
   t0=uni10.UniTensorR([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
   #print mps_R[2*i+1].printDiagram(), bdi
   t0.PutBlock(mps_R[2*i+1].GetBlock())
   t0.SetLabel([1,11,4,3])
   t1=uni10.UniTensorR([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
   t1.PutBlock(mps_R[2*i].GetBlock())
   t1.SetLabel([7,10,6,1])

   T0=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T0.PutBlock(mps_A[2*i+1].GetBlock())
   T0.SetLabel([18,13,4,3])
   T1=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T1.PutBlock(mps_A[2*i].GetBlock())
   T1.SetLabel([17,14,6,18])

   U=Uni_list[2*i]
   U.SetLabel([12,13,10,11])
   Up=Uni_list[2*i-1]
   Up.SetLabel([15,14,16,12])

   A_t=uni10.Contract(T0,t0)
   A_t=uni10.Contract(A_t, U)
   A_t=uni10.Contract(A_t, Up)
   A_t=uni10.Contract(A_t,t1)
   Result=uni10.Contract(A_t,T1)
   #Result=uni10.Contract(A_t, B_t)

   #Result=((((T0*t0)*U)*Up)*(t1))*T1
   Result=uni10.Permute(Result,[7,16,15,17],4)
   Env_right[i]=Result
  else:
   t0=uni10.UniTensorR([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
   t0.PutBlock(mps_R[2*i+1].GetBlock())
   t0.SetLabel([1,11,4,-7])
   t1=uni10.UniTensorR([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
   t1.PutBlock(mps_R[2*i].GetBlock())
   t1.SetLabel([7,10,6,1])

   T0=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
   T0.PutBlock(mps_A[2*i+1].GetBlock())
   T0.SetLabel([18,-15,4,-17])
   T1=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
   T1.PutBlock(mps_A[2*i].GetBlock())
   T1.SetLabel([17,14,6,18])

   U=Uni_list[2*i]
   U.SetLabel([12,-16,10,11])

   Up=Uni_list[2*i-1]
   Up.SetLabel([15,14,16,12])

   Env_right[i+1].SetLabel([-7,-16,-15,-17])
#   Result=((((t1)*(t0))*(U*Up))*(T1*T0))*(Env_right[i+1])

   Result=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_right[i+1],uni10.Contract(T0,t0)),U),Up),uni10.Contract(t1,T1))

   Result=uni10.Permute(Result,[7,16,15,17],4)
   Env_right[i]=Result
 return   Env_right




def Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi):
 N=mps_A.N
 Env_left=[None]*(N/2)
 for i in xrange( 0, (N/2)-1):
   #print i, 2*i+1
   if i == 0:
      #print "i=initLeft", i
      t0=uni10.UniTensorR([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
      t0.PutBlock(mps_R[2*i].GetBlock())
      t0.SetLabel([1,11,4,3])
      t1=uni10.UniTensorR([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
      t1.PutBlock(mps_R[2*i+1].GetBlock())
      t1.SetLabel([3,10,6,7])
      #print mps_R[2*i+1].printDiagram()

      T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      T0.PutBlock(mps_A[2*i].GetBlock())
      T0.SetLabel([1,12,4,18])
      T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      T1.PutBlock(mps_A[2*i+1].GetBlock())
      T1.SetLabel([18,19,6,20])

      U=Uni_list[2*i]
      U.SetLabel([12,13,11,10])


      A_t=uni10.Contract(T0,t0)
      A_t=uni10.Contract(A_t, U)
      B_t=uni10.Contract(T1,t1)
      Result=uni10.Contract(A_t, B_t)

      #Result=(((T0)*(t0))*(U))*(T1*t1)
      Result=uni10.Permute(Result,[7,13,19,20],0)
      Env_left[i]=Result
   else:
      #print "i=insideLeft", i
      t0=uni10.UniTensorR([mps_R[2*i].bond(0), bdi, bdi, mps_R[2*i].bond(2)])
      t0.PutBlock(mps_R[2*i].GetBlock())
      t0.SetLabel([-7,11,4,3])
      t1=uni10.UniTensorR([mps_R[2*i+1].bond(0), bdi, bdi, mps_R[2*i+1].bond(2)])
      t1.PutBlock(mps_R[2*i+1].GetBlock())
      t1.SetLabel([3,10,6,7])
      #print mps_R[2*i+1].printDiagram()

      T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
      T0.PutBlock(mps_A[2*i].GetBlock())
      T0.SetLabel([-20,12,4,18])
      T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
      T1.PutBlock(mps_A[2*i+1].GetBlock())
      T1.SetLabel([18,19,6,20])

      U=Uni_list[2*i]
      U.SetLabel([21,13,11,10])

      Up=Uni_list[2*i-1]
      Up.SetLabel([-19,12,-13,21])

      Env_left[i-1].SetLabel([-7,-13,-19,-20])

#      Result=((((t1)*(t0))*(U))*((T1*T0)*Up))*(Env_left[i-1])
      Result=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_left[i-1],Up),uni10.Contract(t0,T0)),U),uni10.Contract(t1,T1))

      Result=uni10.Permute(Result,[7,13,19,20],4)
      Env_left[i]=Result

 return  Env_left


def   Env_left_RU_update(mps_A, mps_R, Uni_list,i,bdip,bdi,Env_left):
 #print "\n"
 if  i==0:
  #print "i=initLeft", i
  t0=uni10.UniTensorR([mps_R[2*i].bond(0),bdi,bdi,mps_R[2*i].bond(2)])
  t0.PutBlock(mps_R[2*i].GetBlock())
  t0.SetLabel([1,11,4,3])
  t1=uni10.UniTensorR([mps_R[2*i+1].bond(0),bdi,bdi,mps_R[2*i+1].bond(2)])
  t1.PutBlock(mps_R[2*i+1].GetBlock())
  t1.SetLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()

  T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.PutBlock(mps_A[2*i].GetBlock())
  T0.SetLabel([1,12,4,18])
  T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.PutBlock(mps_A[2*i+1].GetBlock())
  T1.SetLabel([18,19,6,20])

  U=Uni_list[2*i]
  U.SetLabel([12,13,11,10])

  #Result=(((t1)*(t0))*(U))*(T1*T0)
  Result=uni10.Contract(uni10.Contract(uni10.Contract(T0,t0),U),uni10.Contract(T1,t1))

  Result=uni10.Permute(Result,[7,13,19,20],0)
  Env_left[i].PutBlock(Result.GetBlock())
 else:
  t0=uni10.UniTensorR([mps_R[2*i].bond(0), bdi, bdi, mps_R[2*i].bond(2)])
  t0.PutBlock(mps_R[2*i].GetBlock())
  t0.SetLabel([-7,11,4,3])
  t1=uni10.UniTensorR([mps_R[2*i+1].bond(0), bdi, bdi, mps_R[2*i+1].bond(2)])
  t1.PutBlock(mps_R[2*i+1].GetBlock())
  t1.SetLabel([3,10,6,7])

  T0=uni10.UniTensorR([mps_A[2*i].bond(0),bdip,bdi,mps_A[2*i].bond(2)])
  T0.PutBlock(mps_A[2*i].GetBlock())
  T0.SetLabel([-20,12,4,18])
  T1=uni10.UniTensorR([mps_A[2*i+1].bond(0),bdip,bdi,mps_A[2*i+1].bond(2)])
  T1.PutBlock(mps_A[2*i+1].GetBlock())
  T1.SetLabel([18,19,6,20])

  U=Uni_list[2*i]
  U.SetLabel([21,13,11,10])

  Up=Uni_list[2*i-1]
  Up.SetLabel([-19,12,-13,21])

  Env_left[i-1].SetLabel([-7,-13,-19,-20])

  #Result=((((t1)*(t0))*(U))*((T1*T0)*Up))*(Env_left[i-1])
  Result=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_left[i-1],Up),uni10.Contract(t0,T0)),U),uni10.Contract(t1,T1))

  Result=uni10.Permute(Result,[7,13,19,20],4)
  Env_left[i].PutBlock(Result.GetBlock())

 return  Env_left


#@profile
def Env_NR_f( mps_A, mps_R, Uni_list, bdip, bdi, position,E_right,E_left, Env_right, Env_left):

 N=mps_R.N
 Env_N=0
 if position == 0:
   E_right[1].SetLabel([3,-3])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[0].bond(1).dim())
   t0=uni10.UniTensorR([mps_R[0].bond(1), bdo])
   t0.Identity()
   t0.SetLabel([-2,2])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[0].bond(0).dim())
   t1=uni10.UniTensorR([mps_R[0].bond(0), bdo])
   t1.Identity()
   t1.SetLabel([-1,1])
   Env_N=uni10.Contract(uni10.Contract(t0,E_right[1]),t1)
   Env_N=uni10.Permute(Env_N,[-1,-2,-3,1,2,3],3)
 elif position == (mps_R.N-1):
   E_left[position-1].SetLabel([1,-1])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[1].bond(1).dim())
   t0=uni10.UniTensorR([mps_R[position].bond(1), bdo])
   t0.Identity()
   t0.SetLabel([-2,2])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[position].bond(2).dim())
   t1=uni10.UniTensorR([mps_R[position].bond(2), bdo])
   t1.Identity()
   t1.SetLabel([-3,3])
   Env_N=uni10.Contract(uni10.Contract(t0,E_left[position-1]),t1)
   Env_N=uni10.Permute(Env_N,[-1,-2,-3,1,2,3],3)
 else:
   E_right[position+1].SetLabel([3,-3])
   E_left[position-1].SetLabel([1,-1])
   bdo = uni10.Bond(uni10.BD_OUT, mps_R[1].bond(1).dim())
   t0=uni10.UniTensorR([mps_R[1].bond(1), bdo])
   t0.Identity()
   t0.SetLabel([-2,2])
   Env_N=uni10.Contract(uni10.Contract(t0,E_right[position+1]),E_left[position-1])
   Env_N=uni10.Permute(Env_N,[-1,-2,-3,1,2,3],3)





 Env_RU=0

 if position==0 or position==1:
  t0=uni10.UniTensorR([mps_R[0].bond(0), bdi, bdi, mps_R[0].bond(2)])
  t0.PutBlock(mps_R[0].GetBlock())
  t0.SetLabel([1,11,4,3])
  t1=uni10.UniTensorR([mps_R[1].bond(0),bdi,bdi,mps_R[1].bond(2)])
  t1.PutBlock(mps_R[1].GetBlock())
  t1.SetLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()
  T0=uni10.UniTensorR([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
  T0.PutBlock(mps_A[0].GetBlock())
  T0.SetLabel([1,12,4,18])
  T1=uni10.UniTensorR([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
  T1.PutBlock(mps_A[1].GetBlock())
  T1.SetLabel([18,19,6,20])

  U=Uni_list[0]
  U.SetLabel([12,13,11,10])
  Env_right[1].SetLabel([7,13,19,20])
  if position%2==0:
#    Env_RU=(((t1)*U)*(Env_right[1]))*(T1*T0)
    Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(Env_right[1],uni10.Contract(t1,T1)),U),T0)

    Env_RU=uni10.Permute(Env_RU,[1,11,4,3],3)
    Env_RU=uni10.Permute(Env_RU,[1,11,4,3],3)
    Env_RU=Env_RU.CombineBond([11,4])
    Env_RU=uni10.Permute(Env_RU,[1,11,3],3)
  else:
#    Env_RU=((t0*U)*Env_right[1])*(T1*T0)



    Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(t0,T0),U),T1),Env_right[1])
    Env_RU=uni10.Permute(Env_RU,[3,10,6,7],3)
    Env_RU=Env_RU.CombineBond([10,6])
    Env_RU=uni10.Permute(Env_RU,[3,10,7],3)

 elif position==(N-1) or position==(N-2):

  t0=uni10.UniTensorR([mps_R[N-1].bond(0),bdi,bdi,mps_R[N-1].bond(2)])
  t0.PutBlock(mps_R[N-1].GetBlock())
  t0.SetLabel([1,11,4,3])
  t1=uni10.UniTensorR([mps_R[N-2].bond(0),bdi,bdi,mps_R[N-2].bond(2)])
  t1.PutBlock(mps_R[N-2].GetBlock())
  t1.SetLabel([7,10,6,1])

  T0=uni10.UniTensorR([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
  T0.PutBlock(mps_A[N-1].GetBlock())
  T0.SetLabel([18,13,4,3])
  T1=uni10.UniTensorR([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
  T1.PutBlock(mps_A[N-2].GetBlock())
  T1.SetLabel([17,14,6,18])

  U=Uni_list[N-2]
  U.SetLabel([12,13,10,11])

  Up=Uni_list[N-3]
  Up.SetLabel([15,14,16,12])

  Env_left[((N-2)/2)-1].SetLabel([7,16,15,17])
  if position==(N-1):
#   Env_RU=(((T1*T0)*(U*Up))*(Env_left[((N-2)/2)-1]))*(t1)
   Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_left[((N-2)/2)-1],t1),T1),Up),U),T0)

   Env_RU=uni10.Permute(Env_RU,[1,11,4,3],3)
   Env_RU=Env_RU.CombineBond([11,4])
   Env_RU=uni10.Permute(Env_RU,[1,11,3],3)

  else:
#   Env_RU=(((T0*T1)*(Up*U))*Env_left[((N-2)/2)-1])*(t0)
   Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(T0,t0),U),uni10.Contract(Env_left[((N-2)/2)-1],uni10.Contract(Up,T1)))
   Env_RU=uni10.Permute(Env_RU,[7,10,6,1],3)
   Env_RU=Env_RU.CombineBond([10,6])
   Env_RU=uni10.Permute(Env_RU,[7,10,1],3)
 else:
  if (position%2)==0:
    #print "Location=IsoMidelodd", Location
    t0=uni10.UniTensorR([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t0.PutBlock(mps_R[position].GetBlock())
    t0.SetLabel([-7,11,4,3])
    t1=uni10.UniTensorR([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t1.PutBlock(mps_R[position+1].GetBlock())
    t1.SetLabel([3,10,6,7])
    #print mps_R[2*i+1].printDiagram()


    T0=uni10.UniTensorR([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T0.PutBlock(mps_A[position].GetBlock())
    T0.SetLabel([-20,12,4,18])
    T1=uni10.UniTensorR([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T1.PutBlock(mps_A[position+1].GetBlock())
    T1.SetLabel([18,19,6,20])

    U=Uni_list[position]
    U.SetLabel([21,13,11,10])

    Up=Uni_list[position-1]
    Up.SetLabel([-19,12,-13,21])

    #print "position", (position/2)-1
    Env_left[(position/2)-1].SetLabel([-7,-13,-19,-20])
    Env_right[(position/2)+1].SetLabel([7,13,19,20])
    #Env_RU=(((t1)*((U*Up))*(Env_right[(position/2)+1]))*(((T1*T0)*Env_left[(position/2)-1])))
    Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_right[(position/2)+1],t1),T1),U),uni10.Contract(uni10.Contract(Env_left[(position/2)-1],Up),T0))

    Env_RU=uni10.Permute(Env_RU,[-7,11,4,3],3)
    Env_RU=Env_RU.CombineBond([11,4])
    Env_RU=uni10.Permute(Env_RU,[-7,11,3],3)

  else:
    t0=uni10.UniTensorR([mps_R[position-1].bond(0),bdi,bdi,mps_R[position-1].bond(2)])
    t0.PutBlock(mps_R[position-1].GetBlock())
    t0.SetLabel([-7,11,4,3])
    t1=uni10.UniTensorR([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t1.PutBlock(mps_R[position].GetBlock())
    t1.SetLabel([3,10,6,7])

    T0=uni10.UniTensorR([mps_A[position-1].bond(0),bdip,bdi,mps_A[position-1].bond(2)])
    T0.PutBlock(mps_A[position-1].GetBlock())
    T0.SetLabel([-20,12,4,18])
    T1=uni10.UniTensorR([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T1.PutBlock(mps_A[position].GetBlock())
    T1.SetLabel([18,19,6,20])


    U=Uni_list[position-1]
    U.SetLabel([21,13,11,10])

    Up=Uni_list[position-2]
    Up.SetLabel([-19,12,-13,21])


    Env_left[((position-1)/2)-1].SetLabel([-7,-13,-19,-20])
    Env_right[(position+1)/2].SetLabel([7,13,19,20])



#    Env_RU=(((t0*(U*Up))*(Env_left[((position-1)/2)-1]))*((T1*T0)*Env_right[(position+1)/2]))
    Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_left[((position-1)/2)-1],Up),T0),t0),U),uni10.Contract(Env_right[(position+1)/2],T1))

    Env_RU=uni10.Permute(Env_RU,[3,10,6,7],3)
    Env_RU=Env_RU.CombineBond([10,6])
    Env_RU=uni10.Permute(Env_RU,[3,10,7],3)

 return   Env_N,  Env_RU





#@profile
def Env_U_f( mps_A, mps_R, Uni_list, bdip, bdi, position, Env_left, Env_right):


 N=mps_R.N

 Env_RU=0

 if position==0:
  t0=uni10.UniTensorR([mps_R[0].bond(0), bdi, bdi, mps_R[0].bond(2)])
  t0.PutBlock(mps_R[0].GetBlock())
  t0.SetLabel([1,11,4,3])
  t1=uni10.UniTensorR([mps_R[1].bond(0),bdi,bdi,mps_R[1].bond(2)])
  t1.PutBlock(mps_R[1].GetBlock())
  t1.SetLabel([3,10,6,7])
  #print mps_R[2*i+1].printDiagram()
  T0=uni10.UniTensorR([mps_A[0].bond(0),bdip,bdi,mps_A[0].bond(2)])
  T0.PutBlock(mps_A[0].GetBlock())
  T0.SetLabel([1,12,4,18])
  T1=uni10.UniTensorR([mps_A[1].bond(0),bdip,bdi,mps_A[1].bond(2)])
  T1.PutBlock(mps_A[1].GetBlock())
  T1.SetLabel([18,19,6,20])

  U=Uni_list[0]
  U.SetLabel([12,13,11,10])
  
  Env_right[1].SetLabel([7,13,19,20])
#  Env_RU=(((t1*t0))*(Env_right[1]))*(T1*T0)
  Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(Env_right[1],uni10.Contract(t1,T1)),t0),T0)

  Env_RU=uni10.Permute(Env_RU,[12,13,11,10],2)

 elif position==(N-2) or position==(N-3):

  t0=uni10.UniTensorR([mps_R[N-1].bond(0),bdi,bdi,mps_R[N-1].bond(2)])
  t0.PutBlock(mps_R[N-1].GetBlock())
  t0.SetLabel([1,11,4,3])
  t1=uni10.UniTensorR([mps_R[N-2].bond(0),bdi,bdi,mps_R[N-2].bond(2)])
  t1.PutBlock(mps_R[N-2].GetBlock())
  t1.SetLabel([7,10,6,1])

  T0=uni10.UniTensorR([mps_A[N-1].bond(0),bdip,bdi,mps_A[N-1].bond(2)])
  T0.PutBlock(mps_A[N-1].GetBlock())
  T0.SetLabel([18,13,4,3])
  T1=uni10.UniTensorR([mps_A[N-2].bond(0),bdip,bdi,mps_A[N-2].bond(2)])
  T1.PutBlock(mps_A[N-2].GetBlock())
  T1.SetLabel([17,14,6,18])

  U=Uni_list[N-2]
  U.SetLabel([12,13,10,11])

  Up=Uni_list[N-3]
  Up.SetLabel([15,14,16,12])

  Env_left[((N-2)/2)-1].SetLabel([7,16,15,17])
  if position==(N-2):
#   Env_RU=((T0*(T1*Up))*(Env_left[((N-2)/2)-1]))*(t1*t0)

   Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_left[((N-2)/2)-1],Up),T1),t1),uni10.Contract(T0,t0))

   Env_RU=uni10.Permute(Env_RU,[12,13,10,11],2)
   #Env_RU.CombineBond([11,4])
   #Env_RU.uni10.Permute([1,11,3],2)
  else:
#   Env_RU=(((t0)*(t1*U))*Env_left[((N-2)/2)-1])*(T1*T0)
   Env_RU=uni10.Contract(uni10.Contract(uni10.Contract(T0,t0),U),uni10.Contract(Env_left[((N-2)/2)-1],uni10.Contract(t1,T1)))

   Env_RU=uni10.Permute(Env_RU,[15,14,16,12],2)
   #Env_RU.CombineBond([10,6])
   #Env_RU.uni10.Permute([7,10,1],2)
 else:
  if (position%2)==0:
    #print "Location=IsoMidelodd", Location
    t0=uni10.UniTensorR([mps_R[position].bond(0),bdi,bdi,mps_R[position].bond(2)])
    t0.PutBlock(mps_R[position].GetBlock())
    t0.SetLabel([-7,11,4,3])
    t1=uni10.UniTensorR([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t1.PutBlock(mps_R[position+1].GetBlock())
    t1.SetLabel([3,10,6,7])
    #print mps_R[2*i+1].printDiagram()


    T0=uni10.UniTensorR([mps_A[position].bond(0),bdip,bdi,mps_A[position].bond(2)])
    T0.PutBlock(mps_A[position].GetBlock())
    T0.SetLabel([-20,12,4,18])
    T1=uni10.UniTensorR([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T1.PutBlock(mps_A[position+1].GetBlock())
    T1.SetLabel([18,19,6,20])

    U=Uni_list[position]
    U.SetLabel([21,13,11,10])

    Up=Uni_list[position-1]
    Up.SetLabel([-19,12,-13,21])

    #print (Location/2)-1, (Location/2)+1
    Env_left[(position/2)-1].SetLabel([-7,-13,-19,-20])
    Env_right[(position/2)+1].SetLabel([7,13,19,20])



#    Env_RU=(((T1*(T0*Up))*(Env_left[(position/2)-1]))*(((t1*t0)*Env_right[(position/2)+1])))
    Env_RU=uni10.Contract((uni10.Contract(uni10.Contract(Env_right[(position/2)+1],t1),T1)),uni10.Contract(uni10.Contract(uni10.Contract(Env_left[(position/2)-1],Up),T0),t0))

    Env_RU=uni10.Permute(Env_RU,[21,13,11,10],2)
    #Env_RU.CombineBond([11,4])
    #Env_RU.uni10.Permute([-7,11,3],2)
  else:
    t0=uni10.UniTensorR([mps_R[position+1].bond(0),bdi,bdi,mps_R[position+1].bond(2)])
    t0.PutBlock(mps_R[position+1].GetBlock())
    t0.SetLabel([-7,11,4,3])
    t1=uni10.UniTensorR([mps_R[position+2].bond(0),bdi,bdi,mps_R[position+2].bond(2)])
    t1.PutBlock(mps_R[position+2].GetBlock())
    t1.SetLabel([3,10,6,7])

    T0=uni10.UniTensorR([mps_A[position+1].bond(0),bdip,bdi,mps_A[position+1].bond(2)])
    T0.PutBlock(mps_A[position+1].GetBlock())
    T0.SetLabel([-20,12,4,18])
    T1=uni10.UniTensorR([mps_A[position+2].bond(0),bdip,bdi,mps_A[position+2].bond(2)])
    T1.PutBlock(mps_A[position+2].GetBlock())
    T1.SetLabel([18,19,6,20])


    U=Uni_list[position+1]
    U.SetLabel([21,13,11,10])

    Up=Uni_list[position]
    Up.SetLabel([-19,12,-13,21])

    Env_left[((position-1)/2)].SetLabel([-7,-13,-19,-20])
    Env_right[((position+1)/2)+1].SetLabel([7,13,19,20])

#    Env_RU=(((t0*(t1*U))*(Env_right[((position+1)/2)+1]))*((T1*T0)*Env_left[((position-1)/2)]))
    Env_RU=uni10.Contract(Env_left[((position-1)/2)],uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(uni10.Contract(Env_right[((position+1)/2)+1],T1),t1),U),t0),T0))

    Env_RU=uni10.Permute(Env_RU,[-19,12,-13,21],2)
    #Env_RU.CombineBond([10,6])
    #Env_RU.uni10.Permute([3,10,7],2)

 return   Env_RU


#@profile  
def solve_linear_eq(A,Ap):
 Ap_h=Ap*1.0
 Ap_h=uni10.Transpose(Ap_h)
 Result=uni10.UniTensorR(Ap.bond())
 blk_qnums = A.BlocksQnum()
 #print blk_qnums, Ap.printDiagram()
 blk_qnums1 = Ap.BlocksQnum()
 for qnum in blk_qnums:
   if qnum in blk_qnums1:
    A_mat=A.GetBlock(qnum)
    Ap_mat=Ap.GetBlock(qnum)
    #A_np=Mat_uni_to_np(A_mat)
    #b_np=Mat_uni_to_np(Ap_mat)
    x_np=np.linalg.lstsq(A_mat, Ap_mat,rcond=-1)[0] 
    #x_np=sp.linalg.lstsq(A_np, b_np)[0] 
    #x=Mat_np_to_Uni(x_np)
    Result.PutBlock(qnum, x_np)
 return Result

#@profile
def Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt):
 Env_R=0
 Result3=mps_A.norm()
 for q in xrange(N_iter[0]):
  E_right=Env_right_R(mps_R)
  E_left=Env_left_R(mps_R)
  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  if Opt_r=="R":
   for i in xrange(mps_R.N): 
    if i%2==0 and i>0:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list, (i/2)-1, bdip, bdi, Env_left)

    Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i, E_right, E_left, Env_right, Env_left)
    Env_RU.SetLabel([1,2,3])
    Env_N.SetLabel([-1,-2,-3,1,2,3])
    mps_R[i].SetLabel([1,2,3])
    mps_R_dagg=mps_R[i]*1.0
    mps_R_dagg.SetLabel([-1,-2,-3])
    Result1=uni10.Contract(uni10.Contract(mps_R_dagg,Env_N),mps_R[i])
    Result2=uni10.Contract(Env_RU,mps_R[i])
    #print "R", i, abs(Result2.GetBlock()[0,0])/(abs(Result1.GetBlock()[0,0]*Result3)**(0.5))
#,Result1[0],Result2[0]
    #print Env_N.printDiagram()
    if Opt[2]=="SVD":
     U, S, V=svd_parity2(Env_N)
     #print "Hiiiiiii"
     U=uni10.Transpose(U)
     V=uni10.Transpose(V)
     #V.transpose()
     #print S
     S=inverse(S)
     #print S
     U.SetLabel([5,6,7,8])
     S.SetLabel([4,5])
     V.SetLabel([1,2,3,4])
     Env_N_inv=uni10.Contract(uni10.Contract(V,S),U)
     Env_N_inv=uni10.Permute(Env_N_inv,[1,2,3,6,7,8],3)
     Env_N_inv.SetLabel([-1,-2,-3,1,2,3])
     R_new=uni10.Contract(Env_N_inv,Env_RU)
     R_new=uni10.Permute(R_new,[-1,-2,-3],2)
    elif  Opt[2]=="cg":
     R_new=solve_linear_eq(Env_N,Env_RU)
     R_new.SetLabel([0,1,2])
     R_new=uni10.Permute(R_new,[0,1,2],2)   

    mps_R[i].PutBlock(R_new.GetBlock())
    E_left=Update_Env_left(mps_R, E_left,i)
    #Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)

  Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
  if Opt_u=="U":
   for i in xrange(mps_R.N-1): 
    if (i%2)==1:
     Env_left=Env_left_RU_update(mps_A, mps_R, Uni_list,(i-1)/2,bdip,bdi,Env_left)
    Env_U=Env_U_f(mps_A, mps_R, Uni_list, bdip, bdi, i, Env_left, Env_right)

    Env_U.SetLabel([1,2,3,4])
    Uni_list[i].SetLabel([1,2,3,4])
    Result1=mps_R.norm()
    Result2=uni10.Contract(Env_U,Uni_list[i])
    #print "U", i, abs(Result2.GetBlock()[0,0])/(abs(Result1*Result3)**(0.5))
#,Result1, Result2[0]
    #svd=Env_U.GetBlock().svd()
    #temporary_matrix=svd[0]*svd[2]
    
    svd=linalg.svd(Env_U.GetBlock(), full_matrices=False, lapack_driver='gesvd')
    temporary_matrix=np.matmul(svd[0], svd[2])

    Uni_list[i].PutBlock(temporary_matrix)

 Env_left=Env_left_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 Env_right=Env_right_RU(mps_A, mps_R, Uni_list,bdip,bdi)
 E_right=Env_right_R(mps_R)
 E_left=Env_left_R(mps_R)
 Env_N, Env_RU=Env_NR_f(mps_A, mps_R, Uni_list, bdip, bdi, i,E_right,E_left,Env_right,Env_left)
 Env_RU.SetLabel([1,2,3])
 Env_N.SetLabel([-1,-2,-3,1,2,3])
 mps_R[i].SetLabel([1,2,3])
 mps_R_dagg=mps_R[i]*1.0
 mps_R_dagg.SetLabel([-1,-2,-3])
 Result1=uni10.Contract(uni10.Contract(mps_R_dagg,Env_N),mps_R[i])
 Result2=uni10.Contract(Env_RU,mps_R[i])
 
 Fidel_final=abs(Result2.GetBlock()[0,0])/(abs(Result1.GetBlock()[0,0]*Result3)**(0.5))


 return  mps_R, Uni_list, Fidel_final





def  Fidel_basedonQR(MPS_A, MPS_R, MPS_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop):

 chi=120
 if MPS_Q.D>chi:
   MPS_Q=MPS_Q.appSVD(chi)

 chi=120
 if MPS_R.D>chi:
   MPS_R=MPS_R.appSVD(chi)


 A_list=[]
 for i in xrange(N_y):
   t0=uni10.UniTensorR([MPS_Q[i].bond(0),bdip,bdi,MPS_Q[i].bond(2)])
   t0.PutBlock(MPS_Q[i].GetBlock())
   t0.SetLabel([-1,2,3,4])
   t1=uni10.UniTensorR([MPS_R[i].bond(0),bdi,bdi,MPS_R[i].bond(2)])
   t1.PutBlock(MPS_R[i].GetBlock())
   t1.SetLabel([1,3,5,-4])
   t_result=uni10.Contract(t0,t1)
   t_result=uni10.Permute(t_result,[1, -1, 2, 5, -4, 4],4)
   t_result=t_result.CombineBond([1,-1])
   t_result=t_result.CombineBond([-4,4])
   t_result=t_result.CombineBond([2,5])
   t_result=uni10.Permute(t_result,[1, 2, -4],2)
   A_list.append(t_result)

 list_bond=[]
 for q in xrange(N_y):
   list_bond.append(A_list[q].bond(2).dim())
 #print max(list_bond)
 mps_QR=MPSclass2.MPS( A_list[1].bond(1).dim(), max(list_bond), N_y)
 #randuni, rand, ortho
 for i in xrange(N_y):
   mps_QR[i]=A_list[i]*1.0
 #print mps_QR[0].printDiagram(),mps_QR.norm() 
 
 #print MPS_R.norm(), MPS_Q.norm(), mps_QR.norm(), mps_QR.product(MPS_A), mps_QR.fidel(MPS_A)

 return    mps_QR.fidel(MPS_A)

def Root_update(mps_A, N, bdi, bdip, bdo, bdop, N_iter,Opt):
 Uni_list=[None]*(N-1)
 for i in xrange(mps_A.N-1):
  if i%2==0:
   t1=uni10.UniTensorR([bdip,bdip,bdo,bdo])
   t1.Identity()
   #svd=t1.GetBlock().svd()
   #t1.PutBlock(svd[0])
   Uni_list[i]=t1
  else:
   t1=uni10.UniTensorR([bdip,bdip,bdop,bdop])
   t1.Identity()
   Uni_list[i]=t1
 
 #mps_R=root.make_mps_R_root(mps_A, N, bdi, bdip, bdo, bdop, bdi.dim())
 #Uni_list=Simple_update_root(mps_A, mps_R, N, bdi, bdip, bdo, bdop,Uni_list)
 chi=0
 if mps_A[1].bond(2).dim()<20:
  chi=mps_A[1].bond(2).dim()+2
 else:
  chi=16
 #print bdi, bdip, mps_A[2].printDiagram()
 mps_R=root.make_mps_R_root(mps_A, N, bdi, bdip, bdo, bdop,chi,N_iter)
 mps_R=mps_R.appSVD(mps_A[1].bond(2).dim())
 #print mps_A[1].PrintDiagram(), mps_R[1].PrintDiagram()
 Uni_list=Simple_update_root(mps_A, mps_R, N, bdi, bdip, bdo, bdop, Uni_list)

 Opt_r="R"
 Opt_u="U"
 mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)

 return mps_R, Uni_list

def uni_shrink(Uni_list):
#  Uni_list_copy=[ Uni_list[i]  for i in xrange(len(Uni_list))]
#  U_ten=copy.copy(Uni_list[len(Uni_list)-1])
#  U_ten.Identity()
#  Uni_list_copy.append(U_ten)
 #print len(Uni_list)-1
 U_l=[]
 for i in xrange(0,len(Uni_list),4):
  #print "StepsEven", i

  Uni_list[i].SetLabel([1,2,3,4])
  Uni_list[i+1].SetLabel([5,6,2,7])
  Uni_list[i+2].SetLabel([7,8,9,10])

  Result=uni10.Contract(uni10.Contract(Uni_list[i],Uni_list[i+1]),Uni_list[i+2])
  Result=uni10.Permute(Result,[1,5,6,8,3,4,9,10],8)
  Result=Result.CombineBond([1,5])
  Result=Result.CombineBond([6,8])
  Result=Result.CombineBond([3,4])
  Result=Result.CombineBond([9,10])
  Result=uni10.Permute(Result,[1,6,3,9],2)
  U_l.append(Result)
  if  (i+3)<= (len(Uni_list)-1):
   #print "Hi"
   t1=uni10.UniTensorR([Uni_list[i+3].bond(0),Uni_list[i+3].bond(3)])
   t1.Identity()
   t1.SetLabel([5,6])
   t2=uni10.UniTensorR([Uni_list[i+3].bond(0),Uni_list[i+3].bond(3)])
   t2.Identity()
   t2.SetLabel([7,8])
   Uni_list[i+3].SetLabel([1,2,3,4])
   Result=uni10.Contract(uni10.Contract(Uni_list[i+3],t1),t2)
   Result=Result.CombineBond([5,1])
   Result=Result.CombineBond([2,7])
   Result=Result.CombineBond([6,3])
   Result=Result.CombineBond([4,8])
   Result=uni10.Permute(Result,[5,2,6,4],2)
   U_l.append(Result)
 return  U_l



def   Sqrt(Landa):
  Landa_cp=Landa*1.0
  blk_qnums=Landa.BlocksQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.GetBlock(qnum).shape[0])
   Landa_cpm=Landa_cp.GetBlock(qnum)
   Landam=Landa_cp.GetBlock(qnum)
   for i in xrange(D):
     if Landam[i][i] > 1.0e-12:
      Landa_cpm[i][i]=Landam[i][i]**(1.00/2.00)
     else:
      Landa_cpm[i][i]=0
   Landa_cp.PutBlock(qnum,Landa_cpm)
  return Landa_cp 

def  back_to_normal(mps_R, bdip, bdi, bdop, bdo):

#   Ten_1=uni10.UniTensorR([mps[i].bond(0), bdip, bdi, mps[i].bond(2)], "Ten_R")
#   #print mps[i].printDiagram()
#   Ten_1.SetLabel([1,2,3,4])
#   Ten_1.PutBlock(mps[i].GetBlock())
#   Ten_2=uni10.UniTensorR([mps[i+1].bond(0), bdip, bdi, mps[i+1].bond(2)], "Ten_R")
#   Ten_2.PutBlock(mps[i+1].GetBlock())
#   Ten_2.SetLabel([4,-2,-3,-4])
#   Results=Ten_1*Ten_2  
#   Results.uni10.Permute([1,2,-2,3,-3,-4], 5)
#   Results.CombineBond([2,-2])
#   Results.CombineBond([3,-3])
#   Results.CombineBond([2,3])
#   Results.uni10.Permute([1,2,-4],2)

 A_list=[]
 for i in xrange(mps_R.N):
  #print "i", i
  t0=uni10.UniTensorR([mps_R[i].bond(0),bdip,bdip,bdi,bdi,mps_R[i].bond(2)])
  t0.PutBlock(mps_R[i].GetBlock())
  t0.SetLabel([1,2,-2,3,-3,-4])
  t0=uni10.Permute(t0,[1,2,3,-2,-3,-4],3)
  U, S, V=svd_parity2(t0)
  U.SetLabel([1,2,3,10])  
  S=Sqrt(S)
  S.SetLabel([10, 4])
  U=uni10.Contract(U,S)
  U=U.CombineBond([2,3])
  U.SetLabel([1,2,4])  

  S.SetLabel([4, -10])
  V.SetLabel([-10,-2,-3,-4])  
  V=uni10.Permute(V,[-10,-2, -3, -4],3)
  V=uni10.Contract(S,V)
  V=uni10.Permute(V,[4,-2,-3,-4],3)
  V=V.CombineBond([-2,-3])
  V.SetLabel([4,-2,-4])  

  A_list.append(U)
  A_list.append(V)

 list_bond=[]
 for q in xrange(len(A_list)):
   list_bond.append(A_list[q].bond(2).dim())
 mps_R=MPSclass2.MPS( A_list[1].bond(1).dim(), max(list_bond), len(A_list))
 for i in xrange(len(A_list)):
  mps_R[i]=A_list[i]*1.0
 return  mps_R


#@profile
def Adding_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt):

 #print "A", mps_A[1].printDiagram(), mps_R[1].printDiagram(), bdi.dim(), bdip.dim(),Uni_list[1].printDiagram(),Uni_list[0].printDiagram() 

 #print "norm_A", mps_A.norm(), mps_R.norm(), len(Uni_list)
 mps_A_update=shrink_mps(mps_A, bdip.dim(), bdi.dim())

 mps_R_update=shrink_mps(mps_R, bdi.dim(), bdi.dim())
 Uni_list_update=uni_shrink(Uni_list)

 D=bdi.dim()
 Dp=bdip.dim()

 bdi=uni10.Bond(uni10.BD_IN, D*D)
 bdo=uni10.Bond(uni10.BD_OUT, D*D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp*Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp*Dp)

 Opt_r="R"
 Opt_u="U"
 mps_R, Uni_list_update, Fidel_RU=Update_RU(mps_A_update, mps_R_update, Uni_list_update, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)
 
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 mps_R=back_to_normal(mps_R, bdi, bdi, bdo, bdo)

 #mps_Q=make_MPS_Q(Uni_list_update)
 #mps_Q=back_to_normal(mps_Q, bdip, bdi, bdop, bdo)
 mps_Q="None"
 return  mps_R, mps_Q, Uni_list_update, Fidel_RU



#@profile
def more_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt):

 mps_A_update=shrink_mps(mps_A, bdip.dim(), bdi.dim())
 mps_R_update=shrink_mps(mps_R, bdi.dim(), bdi.dim())
 Uni_list_update=uni_shrink(Uni_list)

 D=bdi.dim()
 Dp=bdip.dim()

 bdi=uni10.Bond(uni10.BD_IN, D*D)
 bdo=uni10.Bond(uni10.BD_OUT, D*D)

 bdip=uni10.Bond(uni10.BD_IN, Dp*Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp*Dp)
 
 #print "Hi"
 Opt_r="R"
 Opt_u="None"
 N_iter1=[1,10]
 mps_R_update, Uni_list_update, Fidel_RU=Update_RU(mps_A_update, mps_R_update, Uni_list_update, bdip, bdi, Opt_r, Opt_u, N_iter1,Opt)
 #print "Hi1"

 
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 
 bdip=uni10.Bond(uni10.BD_IN, Dp)
 bdop=uni10.Bond(uni10.BD_OUT, Dp)

 mps_R=back_to_normal(mps_R_update, bdi, bdi, bdo, bdo)

 Opt_r="R"
 Opt_u="U"
 mps_R, Uni_list_update, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)


 mps_Q=make_MPS_Q(Uni_list)
 #mps_Q=back_to_normal(mps_Q, bdip, bdi, bdop, bdo)
 #mps_Q="None"
 return  mps_R, mps_Q, Uni_list, Fidel_RU



def Qmps_to_Qtensor_left(Q_mps, D_1, D_2, d):

 bdi1=uni10.Bond(uni10.BD_IN, D_1)
 bdo1=uni10.Bond(uni10.BD_OUT, D_1)

 bdi2=uni10.Bond(uni10.BD_IN, D_2)
 bdo2=uni10.Bond(uni10.BD_OUT, D_2)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 B_list=[None]*Q_mps.N

 for i in xrange(Q_mps.N):
     A=uni10.UniTensorR([Q_mps[i].bond(0),bdi1,bdiphy,bdi2,Q_mps[i].bond(2)])
     A.PutBlock(Q_mps[i].GetBlock())
     A.SetLabel([1,0,2,3,4])
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 return   B_list




def Qmps_to_Qtensor_right(Q_mps, D_1, D_2, d):

 bdi1=uni10.Bond(uni10.BD_IN, D_1)
 bdo1=uni10.Bond(uni10.BD_OUT, D_1)

 bdi2=uni10.Bond(uni10.BD_IN, D_2)
 bdo2=uni10.Bond(uni10.BD_OUT, D_2)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 B_list=[None]*Q_mps.N

 for i in xrange(Q_mps.N):
     A=uni10.UniTensorR([Q_mps[i].bond(0),bdi1,bdiphy,bdi2,Q_mps[i].bond(2)])
     A.PutBlock(Q_mps[i].GetBlock())
     A.SetLabel([1,3,2,0,4])
     A=uni10.Permute(A,[0,1,2,3,4],3)
     B_list[i]=A
 return   B_list



#@profile
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

################PEPS-tensor####################
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

  #print mps_Q[0].PrintDiagram(), mps_A[0].PrintDiagram(), mps_A[1].PrintDiagram()

  return  mps_R,mps_Q, Uni_list



 if Opt[0]=="root":
   mps_A=copy.copy(mps_A_origin)
   mps_R, Uni_list=Root_update(mps_A,N_y, bdi, bdip, bdo, bdop, N_iter,Opt)
   mps_Q=make_MPS_Q(Uni_list)
   Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   print "Fidel=root-A" , abs(Fidel_final)
   return  mps_R, mps_Q, Uni_list
 else:
   mps_A=mps_A_origin*1.0
   Iso_list, Uni_list=Simple_update(mps_A, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, Fidel_list, count_list)
   Isop_list=[ Iso_list[i]*1.0  for i in xrange(len(Iso_list)) ]
   #Iso_list, Isop_list, Uni_list=Full_update(mps_A,Iso_list,Uni_list,N_y, bdi, bdi1, bdip, bdo, bdo1, bdop, 4, Fidel_list, count_list)
   #print "Hi", Uni_list[1].printDiagram()
   mps_R=make_MPS_R(mps_A, Iso_list, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   Uni_list=Uni_update(Isop_list, Uni_list)

   Opt_r="R"
   Opt_u="U"
   mps_R, Uni_list, Fidel_RU=Update_RU(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)
   mps_Q=make_MPS_Q(Uni_list)
   Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   print "Fidel=QR-A" , abs(Fidel_final), Fidel_RU
   #print mps_Q[0].PrintDiagram()
   #print mps_Q[0].PrintDiagram()

   if Opt[1]=="contiguous":
    mps_R, mps_Q, Uni_list, Fidel_RU=Adding_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)
    print "Fidel=QR-A", abs(Fidel_RU)
   elif Opt[1]=="moreAcc":
    #print  "moreAcc"
    mps_R, mps_Q, Uni_list, Fidel_RU=more_update(mps_A, mps_R, Uni_list, bdip, bdi, Opt_r, Opt_u, N_iter,Opt)
    print "Fidel=QR-A", abs(Fidel_RU)



   #Fidel_final=Fidel_basedonQR(mps_A, mps_R, mps_Q, N_y, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #print "Fidel=QR-A" , abs(Fidel_final)
   return  mps_R, mps_Q, Uni_list






 file = open("QR.txt", "w")
 for index in range(len(Fidel_list)):
   file.write(str(count_list[index]) + " " + str(Fidel_list[index])+" "+ "\n")
 file.close()
 #print mps_R.norm()


 #Uni_list[0].SetLabel([1,2,3,4]) 
 #U_f=copy.copy(Uni_list[0]) 
 #U_f.SetLabel([1,2,5,6]) 

 #print Uni_list[0].printDiagram(), U_f*Uni_list[0]



