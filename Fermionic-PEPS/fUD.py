import pyUni10 as uni10
import copy
import numpy as np
import scipy as sp
from numpy import linalg as npLA
from scipy import linalg as LA
import MPSclass 
import math
import TruncateU as TU




def Short_TrotterSteps(N_iterF):
 List_delN=[]


#  Delta_N=(0.05, N_iterF)
#  List_delN.append(Delta_N)
# 
 Delta_N=(0.04, N_iterF)
 List_delN.append(Delta_N)


# Delta_N=(0.005, N_iterF)
# List_delN.append(Delta_N)

#  Delta_N=(0.005, N_iterF)
#  List_delN.append(Delta_N)

 return List_delN




def fermionicOPT(Sys,bdi, bdi1):


 bdo = uni10.Bond(uni10.BD_OUT, bdi.Qlist())
 bdo1 = uni10.Bond(uni10.BD_OUT, bdi1.Qlist())
 bdii = uni10.Bond(uni10.BD_IN, bdi.Qlist())
 bdii1 = uni10.Bond(uni10.BD_IN, bdi1.Qlist())


 T = uni10.UniTensor([bdii,bdii1,bdo,bdo1])
 T.setLabel([1,2,3,4])

 if Sys[0]=="Fer":
  template = np.zeros([bdi.dim(), bdi1.dim(), bdo.dim(), bdo1.dim()])

  for idx in np.ndindex(template.shape):
   if idx[0]==idx[2] and idx[1]==idx[3]:
    if Sys[4]=="Z2":
     if bdi1.Qlist()[idx[1]].prt() == uni10.PRT_ODD and bdi.Qlist()[idx[0]].prt() == uni10.PRT_ODD:
      template[idx] = -1.0
     else:
      template[idx] = +1.0

    if Sys[4]=="U1":
     if bdi.Qlist()[idx[0]].U1() % 2 ==1 and  bdi1.Qlist()[idx[1]].U1() % 2 ==1:
       template[idx] = -1.0
     else:
      template[idx] = +1.0

  T.setRawElem(template.reshape(-1))
 elif Sys[0]=="Bos":
  T.identity()


 T.permute([1,2,3,4],2)

 return T







def   Symmetric_non(PEPS_listten, PEPS_listtenU, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):

   bdi=uni10.Bond(uni10.BD_IN, PEPS_listten[i][j].bond(0).dim())
   bdi1=uni10.Bond(uni10.BD_IN, PEPS_listten[i][j].bond(1).dim())
   bdi2=uni10.Bond(uni10.BD_IN, PEPS_listten[i][j].bond(2).dim())
   bdo=uni10.Bond(uni10.BD_OUT, PEPS_listten[i][j].bond(3).dim())
   bdo1=uni10.Bond(uni10.BD_OUT, PEPS_listten[i][j].bond(4).dim())
   T=uni10.UniTensor([bdi,bdi1,bdi2,bdo,bdo1])
   T.set_zero()
   #blk_qnums = PEPS_listten[i][j].blockQnum()
   #M_tem=[]
   #for qnum in blk_qnums:
    #M_tem.append(UD.Mat_uni_to_np(PEPS_listten[i][j].getBlock(qnum)))
   #Tn=block_diag(  *[ M_tem[i1] for i1 in xrange(len(M_tem)) ]  )
   Tn = get_ndarray(PEPS_listten[i][j])
   T.setRawElem(Tn.reshape(-1))

   PEPS_listtenU[i][j]=T*1.0
   #print  "Hi", i, j, PEPS_listten[i][j], PEPS_listtenU[i][j]
   #print "norm", i, j, PEPS_listtenU[i][j].norm(), PEPS_listten[i][j].norm()

 return  PEPS_listtenU


def get_ndarray(T):
 if type(T) != np.ndarray:
   return ndarray_help(T)
 else:
   assert type(T) == np.ndarray
   nT = np.zeros(T.shape, dtype=object)
   #for i in range(T.shape[0]):
   #  for j in range(T.shape[1]):
   #    nT[i,j] = ndarray_help(T[i,j])
   for idx in np.ndindex(T.shape):
       nT[idx] = ndarray_help(T[idx])
   return nT


def ndarray_help(T):
 nbond = len(T.bond())
 bondD = np.zeros(nbond, dtype=int)
 for i in range(nbond):
   bondD[i] = T.bond(i).dim()
 ret = np.zeros(bondD)
 for idx in np.ndindex(ret.shape):
   ret[idx] = T.at(idx)
 return ret





def full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_in, d_out):

 ######################### No-symmetry #############################################
 if  Model[0] is "Heis"  or  Model[0] is "ITF" or Model[0] is "Fer" or Model[0] is "FFI":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN)
  q_list=[q0_even]
  qchi_list=[q0_even]

  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]

  q_d_in=[q0_even]*d_in[0]
  q_d_out=[q0_even]*d_out[0]

  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out  )

 ###################### Z(2) ######################################
 if Model[0] is "Heis_Z2" or Model[0] is "ITF_Z2" or Model[0] is "Fer_Z2"  or  Model[0] is "FFI_Z2":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN)
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD)
  q_list=[q0_even,q0_odd]
  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]
  q_d_in=[q0_even,q0_odd]
  q_d_out=[q0_even,q0_odd]

  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out  )


 if Model[0] is "Fer_U1": 

  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
  q1_even = uni10.Qnum(1,uni10.PRT_EVEN)
  q2_even = uni10.Qnum(2,uni10.PRT_EVEN)

  q2_odd = uni10.Qnum(2,uni10.PRT_ODD)
  q1_odd = uni10.Qnum(1,uni10.PRT_ODD)

  q_1_odd = uni10.Qnum(-1,uni10.PRT_ODD)
  q_2_odd = uni10.Qnum(-2,uni10.PRT_ODD)

  q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
  q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
  q4_even = uni10.Qnum(4,uni10.PRT_EVEN);
  q5_even = uni10.Qnum(5,uni10.PRT_EVEN);

  q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
  q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
  q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
  q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);
  q_5_even = uni10.Qnum(-5,uni10.PRT_EVEN);

  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  #qchi_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even]
  qchi_list=[q0_even]
  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]

  #qchi_list=[q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even]
  #qchi_list=[q_5_even,q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even,q5_even]
  #qchi_list=[q_1_even,q1_even]

  #q_list=[q_2_even,q_1_even,q0_odd,q1_even,q2_even]
  q_list=[q_1_even,q0_even,q1_even]
  #q_list=[q_2_odd,q_1_odd,q_1_even,q0_odd,q1_even,q1_odd,q2_odd]
  #q_list=[q_1_even,q0_even,q1_even]
  #q_list=[q0_even,q1_even,q2_even]
  #q_list=[q1_even,q3_even,q2_even]
  #q_list=[q_2_even,q_1_even,q0_even,q1_even, q2_even]
  #q_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even, q2_even,q3_even]
  #q_list=[q_1_even,q1_even]

  #q_phys=[q_1_even,q1_even]
  q_phys=[q0_even,q0_even,q1_even,q_1_even]
  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out  )
  q_d_in=q_phys
  q_d_out=q_phys
  
 return q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in, q_d_out


def make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out):
 q_D=[]
 q_chi_boundry=[]
 q_chi_single=[]
 q_chi_try=[]
 q_d_in=[]
 q_d_out=[]

 for i in xrange(len(D)):
  for q in xrange(D[i]):
   q_D.append(q_list[i])


 for i in xrange(len(d_in)):
  for q in xrange(d_in[i]):
   q_d_in.append(q_list[i])

 for i in xrange(len(d_out)):
  for q in xrange(d_out[i]):
   q_d_out.append(q_list[i])



 for i in xrange(len(chi_boundry)):
  for q in xrange(chi_boundry[i]):
   q_chi_boundry.append(qchi_boundry_list[i])


 for i in xrange(len(chi_single)):
  for q in xrange(chi_single[i]):
   q_chi_single.append(qchi_single_list[i])

 for i in xrange(len(chi_try)):
  for q in xrange(chi_try[i]):
   q_chi_try.append(qchi_try_list[i])


 return q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in, q_d_out








def increase_bond( PEPS_listten, D, N_x):

  bdiD=uni10.Bond(uni10.BD_IN, D)
  bdoD=uni10.Bond(uni10.BD_OUT, D)


  T0=uni10.UniTensor([bdiD, PEPS_listten[1][1].bond(3)])
  T0.randomize()
  svd=T0.getBlock().svd()
  T0.putBlock(svd[0])
  #T0.identity()
  T1=copy.copy(T0)
  T0.setLabel([1,0])
  #T1.setLabel([1,2])
  #result=T0*T1
  #result.permute([0,2],1)
  #print result

  for i in xrange(N_x):
   for j in xrange(N_x):
     if i==0 and j==0:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([0,1,2,-3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif i==0 and j==N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([0,-1,2,-3,4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif i==N_x-1 and j==0:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,1,2,3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif i==N_x-1 and j==N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,-1,2,3,4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif j==0 and i>0 and i<N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
#      T0.setLabel([-1,1])
#      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,1,2,-3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif j==N_x-1 and i>0 and i<N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
#      T0.setLabel([-4,4])
#      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,-1,2,-3,4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif i==0 and j>0 and j<N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
#      T0.setLabel([-1,1])
#      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([0,-1,2,-3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     elif i==N_x-1 and j>0 and j<N_x-1:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
#      T0.setLabel([-4,4])
#      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,-1,2,3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
     else:
      #print i, j
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      T0.setLabel([-100,0])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-1,1])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-3,3])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      T0.setLabel([-4,4])
      PEPS_listten[i][j]=PEPS_listten[i][j]*T0
      PEPS_listten[i][j].permute([-100,-1,2,-3,-4],3)
      PEPS_listten[i][j].setLabel([0,1,2,3,4])
      
      
  return  PEPS_listten


def inv_landa_col_row( Landa_col, Landa_row, N_x): 

 Landa_col_inv=[None]*N_x
 for i in xrange(N_x):
  Landa_col_inv[i]=[None]*(N_x+1)

 Landa_row_inv=[None]*(N_x+1)
 for i in xrange(N_x+1):
  Landa_row_inv[i]=[None]*(N_x)

 for i in xrange(N_x):
  for j in xrange(N_x+1):
    Landa_col_inv[i][j]=inverse(Landa_col[i][j])

 for i in xrange(N_x+1):
  for j in xrange(N_x):
    Landa_row_inv[i][j]=inverse(Landa_row[i][j])

 return   Landa_col_inv,   Landa_row_inv



def Landa_f_col(D, N_x, N_y):

 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 bdi = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, D)
 Landa=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa.identity()
 Landa.randomize()
 #Landa=Landa*.0001
 #Landa.orthoRand()
 #Landa=sparce_init(Landa)

 bdi = uni10.Bond(uni10.BD_IN, 1)
 bdo = uni10.Bond(uni10.BD_OUT, 1)
 Landa_iden=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_iden.identity()


 Landa_col=[None]*N_x
 for i in xrange(N_x):
  Landa_col[i]=[None]*(N_y+1)

 for i in xrange(N_x):
  for j in xrange(N_y+1):
   if j==0 or j==N_y:
    Landa_col[i][j]=Landa_iden*1.0
   else:
    Landa_col[i][j]=Landa*1.0

 return   Landa_col


def Landa_f_row(D, N_x,N_y):

 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 bdi = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, D)
 Landa=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa.randomize()
 #Landa.identity()
 #Landa.orthoRand()
 #Landa=sparce_init(Landa)

 bdi = uni10.Bond(uni10.BD_IN, 1)
 bdo = uni10.Bond(uni10.BD_OUT, 1)
 Landa_iden=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_iden.identity()


 Landa_row=[None]*(N_x+1)
 for i in xrange(N_x+1):
  Landa_row[i]=[None]*(N_y)
 
 for i in xrange(N_x+1):
  for j in xrange(N_y):
   if i==0 or i==N_x:
    Landa_row[i][j]=Landa_iden*1.0
   else:
    Landa_row[i][j]=Landa*1.0

 return Landa_row

def sparce_init(Landa):
 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 blk_qnums = Landa.blockQnum()
 for qnum in blk_qnums:
  M=Landa.getBlock(qnum)
  if qnum == q0_even:
   M.randomize()
   M=M
   M[0]=1.00
  else: 
   M.randomize()
   M=M*0.01
   #M[0]=.200
  Landa.putBlock(qnum,M)

 return Landa

def norm_Symmetry(LA):
 norm=0
 blk_qnums = LA.blockQnum()
 for qnum in blk_qnums:
  M=LA.getBlock(qnum)  
  norm=norm+(M.norm()*M.norm())
 norm=norm**(1.00/2.00)
 return norm

def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)

  for i in xrange(D):
   invL2[i] = 0 if ((invLt[i].real) < 1.0e-12) else (1.00 / (invLt[i].real))
  invLanda2.putBlock(qnum,invL2)
 return invLanda2

############@profile
def update_row( PEPS_A, PEPS_B, U_ham, i_x, j_y, Landa_col, Landa_row, D,Sys, N_iterF):

 bdi = uni10.Bond(uni10.BD_IN, D)
 D_dim=[]
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)

 PEPS_A.setLabel([1,2,3,4,5])
 PEPS_B.setLabel([6,7,8,9,10])

 B=PEPS_B*1.0
 B.setLabel([6,-7,8,-9,-10])
 B.permute([6,-7,8,-9,-10],5)
 Swap=fermionicOPT(Sys,B.bond(1), B.bond(2))
 Swap.setLabel([7,-8,-7,8])
 B=Swap*B
 B.permute([6,7,-8,-9,-10],3)
 B.setLabel([6,7,8,9,10])

 U_ham.setLabel([11,12,3,8])

 Landa_row[i_x][j_y].setLabel([-1,1])
 Landa_row[i_x+1][j_y].setLabel([4,6])
 Landa_row[i_x+2][j_y].setLabel([9,-9])
 Landa_col[i_x][j_y].setLabel([-2,2])
 Landa_col[i_x][j_y+1].setLabel([5,-5])
 Landa_col[i_x+1][j_y].setLabel([-7,7])
 Landa_col[i_x+1][j_y+1].setLabel([10,-10])

 PEPS_AA=((((PEPS_A*Landa_row[i_x][j_y])*Landa_row[i_x+1][j_y])*Landa_col[i_x][j_y])*Landa_col[i_x][j_y+1])
 PEPS_BB=(((B*Landa_row[i_x+2][j_y])*Landa_col[i_x+1][j_y])*Landa_col[i_x+1][j_y+1])

 PEPS_AA.permute([-1,-2,3,6,-5],3)
 PEPS_BB.permute([6,-7,8,-9,-10],3)


 if N_iterF[1]=="full":
  Theta=(PEPS_AA*U_ham)*PEPS_BB
  Theta.permute([-1,-2,11,-5,-7,12,-9,-10],4)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ###########################################################################

  blk_qnums = LA.blockQnum()
  Landa_row[i_x+1][j_y].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_row[i_x+1][j_y].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-1,-2,11,-5,6])
  U.permute([-1,-2,11,6,-5],3)

  V.setLabel([6,-7,12,-9,-10])
  V.permute([6,12,-7,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_col[i_x][j_y+1])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-5,5])

  U=((U*invLanda1)*invLanda2)*invLanda3
  U.permute([1,2,11,6,5],3)

  invLanda8=inverse_ten(Landa_row[i_x+2][j_y])
  invLanda9=inverse_ten(Landa_col[i_x+1][j_y])
  invLanda10=inverse_ten(Landa_col[i_x+1][j_y+1])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([7,-7])
  invLanda10.setLabel([-10,10])

  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

  B=V*1.0
  B.setLabel([6,-7,8,-9,-10])
  B.permute([6,-7,8,-9,-10],5)
  Swap=fermionicOPT( Sys, B.bond(1), B.bond(2))
  Swap.setLabel([7,-8,-7,8])
  B=Swap*B
  B.permute([6,7,-8,-9,-10],3)
  B.setLabel([6,7,8,9,10])
  
 elif N_iterF[1]=="QR":

  A=copy.copy(PEPS_AA)
  A.setLabel([-1,-2,3,6,-5])
  A.permute([-1,-2,-5,3,6],3)


  row, colm=cal_rowcol(A)
  if (row<=colm):
   q, V, s=TU.setTruncation(A, row)
  else:
   q, V, s=TU.setTruncation(A, colm)


  s.setLabel([1,0])
  V.setLabel([0,3,6])
  r_u=V*s
  r_u.permute([1,3,6],1)
  q.setLabel([-1,-2,-5,-100])
  r_u.setLabel([-100,3,6])


  A=copy.copy(PEPS_BB)
  A.permute([6,8,-7,-9,-10],2)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   U, qq, s=TU.setTruncation(A, row)
  else:
   U, qq, s=TU.setTruncation(A, colm)

  s.setLabel([0,1])
  U.setLabel([6,8,0])
  l_u=U*s
  l_u.permute([6,8,1],2)

  qq.setLabel([-400,-7,-9,-10])
  l_u.setLabel([6,8,-400])

  Theta=(l_u*U_ham)*r_u
  Theta.permute([-100,11,-400,12],2)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ##########################################################################################

  blk_qnums = LA.blockQnum()
  Landa_row[i_x+1][j_y].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_row[i_x+1][j_y].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-100,11,6])
  U=U*q
  U.permute([-1,-2,11,6,-5],3)

  V.setLabel([6,-400,12])
  V=V*qq
  #V.permute([6,-7,12,-9,-10],3)

  V.permute([6,12,-7,-9,-10],3)

 # Swap=fermionicOPT(Sys,V.bond(1), V.bond(2))
 # Swap.setLabel([12,-7,-12,7])
 # #print Swap.printDiagram(), A.printDiagram(), Swap
 # V=Swap*V
 # V.permute([6,7,-12,-9,-10],3)
 # V.setLabel([6,-7,12,-9,-10])


  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_col[i_x][j_y+1])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-5,5])

  U=((U*invLanda1)*invLanda2)*invLanda3
  U.permute([1,2,11,6,5],3)

  invLanda8=inverse_ten(Landa_row[i_x+2][j_y])
  invLanda9=inverse_ten(Landa_col[i_x+1][j_y])
  invLanda10=inverse_ten(Landa_col[i_x+1][j_y+1])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([7,-7])
  invLanda10.setLabel([-10,10])

  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

  B=V*1.0
  B.setLabel([6,-7,8,-9,-10])
  B.permute([6,-7,8,-9,-10],5)
  Swap=fermionicOPT(Sys,B.bond(1), B.bond(2))
  Swap.setLabel([7,-8,-7,8])
  B=Swap*B
  B.permute([6,7,-8,-9,-10],3)
  B.setLabel([6,7,8,9,10])


 return U, B



############@profile
def update_col(PEPS_A, PEPS_B, U_ham, i_x, j_y, Landa_col, Landa_row, D,Sys, N_iterF):
 bdi = uni10.Bond(uni10.BD_IN, D)
 D_dim=[]
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)

 PEPS_A.setLabel([1,2,3,4,5])
 A=PEPS_A*1.0
 A.setLabel([ 1, 2, 3, 4, 5])
 A.permute([ 1, 2, 3, 4, 5], 6)
 Swap=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap.setLabel([-3,-4,3,4])
 #print A.printDiagram(), Swap.printDiagram()
 A=Swap*A
 A.permute([ 1, 2, -3, -4, 5], 3)
 A.setLabel([ 1, 2, 3, 4, 5])


 PEPS_B.setLabel([6,7,8,9,10])
 U_ham.setLabel([11,12,3,8])

 Landa_row[i_x][j_y].setLabel([-1,1])
 Landa_row[i_x+1][j_y].setLabel([4,-4])
 Landa_row[i_x][j_y+1].setLabel([-6,6])
 Landa_row[i_x+1][j_y+1].setLabel([9,-9])

 Landa_col[i_x][j_y].setLabel([-2,2])
 Landa_col[i_x][j_y+1].setLabel([5,7])
 Landa_col[i_x][j_y+2].setLabel([10,-10])

 PEPS_AA=((((A*Landa_row[i_x][j_y])*Landa_row[i_x+1][j_y])*Landa_col[i_x][j_y])*Landa_col[i_x][j_y+1])
 PEPS_BB=(((PEPS_B*Landa_row[i_x][j_y+1])*Landa_row[i_x+1][j_y+1])*Landa_col[i_x][j_y+2])

 PEPS_AA.permute([-1,-2,3,-4,7],3)
 PEPS_BB.permute([-6,7,8,-9,-10],3)

 if N_iterF[1]=="full":
  Theta=(PEPS_AA*U_ham)*PEPS_BB
  Theta.permute([-1,-2,11,-4,-6,12,-9,-10],4)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ##############################################################################

  blk_qnums = LA.blockQnum()
  Landa_col[i_x][j_y+1].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_col[i_x][j_y+1].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-1,-2,11,-4,7])
  U.permute([-1,-2,11,-4,7],3)

  V.setLabel([7,-6,12,-9,-10])
  V.permute([-6,7,12,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_row[i_x+1][j_y])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-4,4])

  U=((U*invLanda1)*invLanda2)*invLanda3

  U.permute([1,2,11,4,7],3)
  U.setLabel([1,2,3,4,5])
  U=U*1.0
  U.setLabel([ 1, 2, 3, 4, 5])
  U.permute([ 1, 2, 3, 4, 5], 6)
  Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
  Swap.setLabel([-3,-4,3,4])
  U=Swap*U
  U.permute([ 1, 2, -3, -4, 5], 3)
  U.setLabel([ 1, 2, 3, 4, 5])


  invLanda8=inverse_ten(Landa_row[i_x+1][j_y+1])
  invLanda9=inverse_ten(Landa_row[i_x][j_y+1])
  invLanda10=inverse_ten(Landa_col[i_x][j_y+2])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([6,-6])
  invLanda10.setLabel([-10,10])
  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

 elif N_iterF[1]=="QR":

  A=PEPS_AA*1.0
  A.permute([ -1, -2, -4, 3, 7], 3)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   q, V, s=TU.setTruncation(A, row)
  else:
   q, V, s=TU.setTruncation(A, colm)

  s.setLabel([1,0])
  V.setLabel([ 0, 3, 7])
  r_u=V*s
  r_u.permute([ 1, 3, 7],1)
  q.setLabel([ -1, -2, -4, -100])
  r_u.setLabel([ -100, 3, 7])

  A=copy.copy(PEPS_BB)
  A.setLabel([-6,7,8,-9,-10])
  A.permute([7,8,-6,-9,-10],2)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   U, qq, s=TU.setTruncation(A, row)
  else:
   U, qq, s=TU.setTruncation(A, colm)

  s.setLabel([0,1])
  U.setLabel([7,8,0])
  l_u=U*s
  l_u.permute([7,8,1],2)

  qq.setLabel([-400,-6,-9,-10])
  l_u.setLabel([7,8,-400])

  Theta=(l_u*U_ham)*r_u
  Theta.permute([-100,11,-400,12],2)

  #print Theta.printDiagram(), sum(D_dim)
  U, V, LA=TU.setTruncation(Theta, sum(D_dim))

  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
  #print LA

 ##########################################################################################

  blk_qnums = LA.blockQnum()
  Landa_col[i_x][j_y+1].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_col[i_x][j_y+1].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-100,11,7])
  U=U*q
  U.permute([-1,-2,11,-4,7],3)

  V.setLabel([7,-400,12])
  V=V*qq
  V.permute([-6,7,12,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_row[i_x+1][j_y])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-4,4])

  U=((U*invLanda1)*invLanda2)*invLanda3

  U.permute([1,2,11,4,7],3)
  U.setLabel([1,2,3,4,5])
  U=U*1.0
  U.setLabel([ 1, 2, 3, 4, 5])
  U.permute([ 1, 2, 3, 4, 5], 6)
  Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
  Swap.setLabel([-3,-4,3,4])
  U=Swap*U
  U.permute([ 1, 2, -3, -4, 5], 3)
  U.setLabel([ 1, 2, 3, 4, 5])

  invLanda8=inverse_ten(Landa_row[i_x+1][j_y+1])
  invLanda9=inverse_ten(Landa_row[i_x][j_y+1])
  invLanda10=inverse_ten(Landa_col[i_x][j_y+2])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([6,-6])
  invLanda10.setLabel([-10,10])
  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

 return   U,  V



############@profile########################
def simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iterF, h_coupling, d, Model, N_x, N_y, D,chi_try, threshold, interval, Sys):

 H=Heisenberg( h_coupling, d, Model)
 H0=Heisenberg0( h_coupling, d, Model)
 HN=HeisenbergN( h_coupling, d, Model)

 #print H, H0, HN

 for i in xrange(1,600):

  delta=start_itebd/pow(division_itebd,i) 

  if delta>1.0e-1:
   N_iter=N_iterF
  if delta<1.0e-1 and delta>1.0e-3:
   N_iter=N_iterF
  if delta<1.0e-3  and delta>1.0e-5:
   N_iter=N_iterF
  if delta<1.0e-10:
   break

  U_ham = uni10.UniTensor( H.bond(), "U")
  blk_qnums = H.blockQnum()
  for qnum in blk_qnums:
   U_ham.putBlock(qnum, uni10.takeExp(-delta, H.getBlock(qnum)))

  U_ham0 = uni10.UniTensor( H0.bond(), "U");
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
   U_ham0.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))

  U_hamN = uni10.UniTensor( HN.bond(), "U");
  blk_qnums = HN.blockQnum()
  for qnum in blk_qnums:
   U_hamN.putBlock(qnum, uni10.takeExp(-delta, HN.getBlock(qnum)))



  U_evolv=1.0
  print 'delta =', delta
  print "N_iterF=", N_iterF



#  for i_x in xrange(N_x):
#   for j_y in xrange(N_x):
#    PEPS_listten[i_x][j_y]=max_ten(PEPS_listten[i_x][j_y])
    #norm=norm_Symmetry(PEPS_listten[i_x][j_y])
    #PEPS_listten[i_x][j_y]=PEPS_listten[i_x][j_y]*(1.00/norm)

#  PEPS_listten=make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)
#  PEPS_listten, norm_val, count_val=Normalize_PEPS( PEPS_listten, N_x, D, chi_try, d, threshold, interval)
#  Landa_col_inv,Landa_row_inv=inv_landa_col_row( Landa_col, Landa_row, N_x)
#  PEPS_listten=make_PEPS_tensors(PEPS_listten, Landa_row_inv, Landa_col_inv,N_x)


  for q in xrange(N_iterF[0]):

   for i_x in xrange(N_x-1):
    for j_y in xrange(N_y):

     if i_x==0:
      U_evolv=U_ham0*1.0
     elif i_x==N_x-2:
      U_evolv=U_hamN*1.0
     else:
      U_evolv=U_ham*1.0

     Peps_A, Peps_B=update_row(PEPS_listten[i_x][j_y], PEPS_listten[i_x+1][j_y], U_evolv, i_x, j_y, Landa_col, Landa_row, D,Sys,N_iterF)
     PEPS_listten[i_x][j_y]=Peps_A*1.0
     PEPS_listten[i_x+1][j_y]=Peps_B*1.0

   for j_y in xrange(N_y-1):
    for i_x in xrange(N_x):
     if j_y==0:
      U_evolv=U_ham0*1.0
     elif j_y==N_y-2:
      U_evolv=U_hamN*1.0
     else:
      U_evolv=U_ham*1.0
     Peps_A, Peps_B=update_col(PEPS_listten[i_x][j_y], PEPS_listten[i_x][j_y+1], U_evolv, i_x, j_y, Landa_col, Landa_row, D,Sys,N_iterF)
     PEPS_listten[i_x][j_y]=Peps_A*1.0
     PEPS_listten[i_x][j_y+1]=Peps_B*1.0

# for i_x in xrange(N_x):
#  for j_y in xrange(N_x):
#   #PEPS_listten[i_x][j_y]=max_ten(PEPS_listten[i_x][j_y])
#   norm=norm_Symmetry(PEPS_listten[i_x][j_y])
#   PEPS_listten[i_x][j_y]=PEPS_listten[i_x][j_y]*(1.00/norm)


 return PEPS_listten, Landa_col, Landa_row




def make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x, N_y): 

 for i_x in xrange(N_x):
  for j_y in xrange(N_y):
   PEPS_listten[i_x][j_y].setLabel([1,2,3,4,5]) 
   Landa_col[i_x][j_y].setLabel([-2,2])
   Landa_row[i_x][j_y].setLabel([-1,1])
   result=(PEPS_listten[i_x][j_y]*Landa_row[i_x][j_y])*Landa_col[i_x][j_y]
   result.permute([-1,-2,3,4,5],3)
   PEPS_listten[i_x][j_y]=result*1.0

 return PEPS_listten







def  Store_Landa_row(Landa_row, N_x, N_y):
 for j in xrange(N_y):
  for i in xrange(N_x+1):
   Landa_row[i][j].save("StoreGamma/row" + str(i)+"P"+str(j))



def  Store_Landa_col(Landa_col, N_x, N_y):
 for j in xrange(N_y+1):
  for i in xrange(N_x):
   Landa_col[i][j].save("StoreGamma/col" + str(i)+"P"+str(j))



def  Reload_Landa_row(Landa_row, N_x, N_y):
 for j in xrange(N_y):
  for i in xrange(N_x+1):
   Landa_row[i][j]=uni10.UniTensor("StoreGamma/row" + str(i)+"P"+str(j))

def  Reload_Landa_col(Landa_col, N_x, N_y):
 for j in xrange(N_y+1):
  for i in xrange(N_x):
   Landa_col[i][j]=uni10.UniTensor("StoreGamma/col" + str(i)+"P"+str(j))







def  Store_Gamma(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j].save("StoreGamma/Gamma" + str(i)+"P"+str(j))


def  Reload_Gamma(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j]=uni10.UniTensor("StoreGamma/Gamma" + str(i)+"P"+str(j))






def Store_f(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j].save("Store/a" + str(i)+"P"+str(j))

def Reload_f(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j]=uni10.UniTensor("Store/a" + str(i)+"P"+str(j))




def Store_fRG(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j].save("StoreRG/a" + str(i)+"P"+str(j))

def Reload_fRG(PEPS_listten, N_x, N_y):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listten[i][j]=uni10.UniTensor("StoreRG/a" + str(i)+"P"+str(j))









def Store_mps(MPS_A, N_x):
 for i in xrange(N_x):
   MPS_A[i].save("StoreMPS/mps" + str(i))

def Reload_mps(MPS_A, N_x):
 for i in xrange(N_x):
   MPS_A[i]=uni10.UniTensor("StoreMPS/mps" + str(i))

def Store_Q_list(Q_list, N_x):
 for i in xrange(N_x):
   Q_list[i].save("StoreQ/Q" + str(i))

def Reload_Q_list(Q_list, N_x):
 for i in xrange(N_x):
   Q_list[i]=uni10.UniTensor("StoreQ/Q" + str(i))






def   copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU):
 for i in xrange(N_x):
  for j in xrange(N_y):
   PEPS_listtenU[i][j]=PEPS_listten[i][j]*1.0

 return  PEPS_listtenU



####@profile
def Energy_cal( PEPS_listten, d, chi_single, N_x, N_y, D, Model, h_coupling, Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_up=[None]*(N_y+1)
 mps_boundry_down=[None]*N_y

 H=Heisenberg(h_coupling, d, Model)
 H0=Heisenberg0(h_coupling, d, Model)
 HN=HeisenbergN(h_coupling, d, Model)

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y*2)
 for i_ind in xrange(N_y*2):
  T.identity()
  mps_I[i_ind]=T*1.0
 
 
 
##################   Col   #####################

 for Location in reversed(xrange(1,N_x)):
  mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_y,Sys, N_x)

 for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_y,Sys)

 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_col(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_col(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_col(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

  if  i_ind==(N_x-1):break


 E_c=sum(E_coulmn_t)
 E_coulmn=[   E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))   ]


 mps_II=MPSclass.MPS(1,1,N_x*2)
 for i_ind in xrange(N_x*2):
  T.identity()
  mps_II[i_ind]=T*1.0


 #print "\n"
#################   Row   #################
 for Location in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_down[Location]=make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], d, chi_single, N_x,Sys)


 for Location in reversed(xrange(1,N_y)):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys, N_y)


 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_row(mps_II, mps_boundry_up[i_ind+1], peps_l,N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_row(mps_boundry_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)


  if  i_ind==(N_y-1):break

 E_r=sum(E_row_t)
 #print "E=", E_r, E_c, (E_c+E_r)
 E_1=(E_c+E_r)
 
 return E_1







def Particle_cal( PEPS_listten, d, chi_single, N_x, N_y, D, Model, h_coupling, Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_up=[None]*(N_y+1)
 mps_boundry_down=[None]*N_y

 H=N_total(h_coupling, d, Model)
 H0=N_total0(h_coupling, d, Model)
 HN=N_totalN(h_coupling, d, Model)

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y*2)
 for i_ind in xrange(N_y*2):
  T.identity()
  mps_I[i_ind]=T*1.0
 
 
 
##################   Col   #####################

 for Location in reversed(xrange(1,N_x)):
  mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_y,Sys, N_x)

 for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_y,Sys)

 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_col(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_col(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_col(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

  if  i_ind==(N_x-1):break


 E_c=sum(E_coulmn_t)
 E_coulmn=[   E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))   ]


 mps_II=MPSclass.MPS(1,1,N_x*2)
 for i_ind in xrange(N_x*2):
  T.identity()
  mps_II[i_ind]=T*1.0


 #print "\n"
#################   Row   #################
 for Location in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_down[Location]=make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], d, chi_single, N_x,Sys)


 for Location in reversed(xrange(1,N_y)):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys, N_y)


 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_row(mps_II, mps_boundry_up[i_ind+1], peps_l,N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_row(mps_boundry_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)


  if  i_ind==(N_y-1):break

 E_r=sum(E_row_t)
 #print E_row_t, E_coulmn_t
 E_1=(E_c+E_r)
 
 return E_1









def Sz_cal( PEPS_listten, d, chi_single, N_x, N_y, D, Model, h_coupling, Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_up=[None]*(N_y+1)
 mps_boundry_down=[None]*N_y

 H=Sz_total(h_coupling, d, Model)
 H0=Sz_total0(h_coupling, d, Model)
 HN=Sz_totalN(h_coupling, d, Model)

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y*2)
 for i_ind in xrange(N_y*2):
  T.identity()
  mps_I[i_ind]=T*1.0
 
 
 
##################   Col   #####################

 for Location in reversed(xrange(1,N_x)):
  mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_y,Sys, N_x)

 for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_y,Sys)

 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_col(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_col(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_col(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

  if  i_ind==(N_x-1):break


 E_c=sum(E_coulmn_t)
 E_coulmn=[   E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))   ]


 mps_II=MPSclass.MPS(1,1,N_x*2)
 for i_ind in xrange(N_x*2):
  T.identity()
  mps_II[i_ind]=T*1.0


 #print "\n"
#################   Row   #################
 for Location in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_down[Location]=make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], d, chi_single, N_x,Sys)


 for Location in reversed(xrange(1,N_y)):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys, N_y)


 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_row(mps_II, mps_boundry_up[i_ind+1], peps_l,N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_row(mps_boundry_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)


  if  i_ind==(N_y-1):break

 E_r=sum(E_row_t)
 #print E_row_t, E_coulmn_t
 E_1=(E_c+E_r)
 
 return E_1


















#######@profile
def  energy_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, H0, H, HN , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[2*i+1]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[2*i+1]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=((mps_boundry_up[2*i+1]*Swap1))*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Yoho1", A[0]

 for i in xrange(N_x-1):

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  E_val=local_energy_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, H_orig, D,Sys)
  E_coulmn.append(E_val)



#######@profile
def  local_energy_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, H_orig, D,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],2)
 r_u.setLabel([-100,3,9])
 r_u.permute([-100,3,9],2)

 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 #####################################################
 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(4))
 Swap1.setLabel([3,10,-3,-10])
 A_conj=Swap1*A_conj
 A_conj.permute([-6,-7,3,-9,10],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-10,-3,-9],3)

 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-9])
 r_d=V*s
 r_d.permute([1,-3,-9],2)

 q_d.setLabel([-6,-7,-10,-200])
 q_d.permute([-6,-7,-10,-200],4)

 r_d.setLabel([-200,-3,-9])
 r_d.permute([-200,-3,-9],2)



 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],2)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 l_u.setLabel([9,13,-300])
 l_u.permute([9,13,-300],2)


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-9,-10,-13,-14,-15])
 A_conj.permute([-9,-13,-10,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-9,-13,0])
 l_d=U*s
 l_d.permute([-9,-13,1],2)

 qq_d.setLabel([-400,-10,-14,-15])
 qq_d.permute([-400,-10,-14,-15],4)


 l_d.setLabel([-9,-13,-400])
 l_d.permute([-9,-13,-400],2)


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 #print Swap4.printDiagram()
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,26])
 mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location*2+2].setLabel([27,15,28])
 mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,31])
 mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location*2+2].setLabel([32,23,33])
 mps_boundry_down[Location*2+3].setLabel([33,-10,19])


# mps_boundry_left[Location].setLabel([16,6,-60,18])
# mps_boundry_left[Location+1].setLabel([18,11,-110,25])


# mps_boundry_right[Location].setLabel([17,90,-9,19])
# mps_boundry_right[Location+1].setLabel([19,140,-14,24])

######################################################

 A=E_left*mps_boundry_up[Location*2]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_up[2*Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[2*Location+1]


 B=E_right*mps_boundry_up[2*Location+3]
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_down[2*Location+3]

 B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################



######################################################

# l_up=copy.copy(l_u)
# r_up=copy.copy(r_u)
# l_dp=copy.copy(l_up)
# r_dp=copy.copy(r_up)
# l_dp.transpose()
# r_dp.transpose()
# l_dp.setLabel([-400,-10,-13])
# r_dp.setLabel([-10,-200,-3 ])

 
 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig

 #print "Fnit_norm", Norm_h

 #print  Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]




#######@profile
def  energy_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, H0, H, HN , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[2*i+1]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[2*i+1]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])


   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*mps_boundry_left[2*i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*mps_boundry_right[2*i]
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])
   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_left[2*i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=mps_boundry_right[2*i]*E_list_up[i]
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_down[i].setLabel([1,2,3,4])
#  E_list_up[i+1].setLabel([1,2,3,4])
#  A=E_list_down[i]*E_list_up[i+1]
#  print "Yoho", A[0]



 for i in xrange(N_x-1):

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  E_val=local_energy( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H_orig,D,Sys)
  E_coulmn.append(E_val)

#######@profile
def  local_energy( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, D,Sys):

 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)


 r_u.setLabel([-100,3,10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])


######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=A*(Swap1*q_d)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap3*q)
 A=A*mps_boundry_right[2*Location+1]


 B=E_right*mps_boundry_left[2*Location+3]
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[2*Location+3]

 B=B*mps_boundry_left[2*Location+2]
 B=B*(Swap4*qq_d)
 B=B*mps_boundry_right[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################
# l_up=copy.copy(l_u)
# r_up=copy.copy(r_u)
# l_dp=copy.copy(l_up)
# r_dp=copy.copy(r_up)
# l_dp.transpose()
# r_dp.transpose()
# l_dp.setLabel([-400,-10,-13])
# r_dp.setLabel([-10,-200,-3 ])

 
 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig

 #print "Fnit_norm", Norm_h

 #print   Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]



def  update_env_LR(PEPS_listten, N_x, d, chi_single,mps_boundry_left, mps_boundry_right, mps_boundry_temp):
 for Location in xrange(N_x-1):
  PEPS_ten=rotate(PEPS_listten[N_x-Location-1])
  mps_boundry_temp[Location]=make_Env_singleLayer(PEPS_ten, Location, mps_boundry_temp[Location-1],d,chi_single,N_x)
  mps_boundry_right[N_x-Location-1]=copy.copy(mps_boundry_temp[Location])
  for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_x)


 return mps_boundry_right, mps_boundry_left



###########@profile
def Normalize_PEPS(PEPS_listten, N_x,N_y, D, chi_try, d, threshold, interval,Sys):

 norm_val=Cal_norm( PEPS_listten, N_x, N_y,D, chi_try, d,Sys)
 print "Zero_order", norm_val

 count=0

 while count<100:
  #print norm_val, count
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   for i in xrange(N_x):
    for j in xrange(N_y):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten,N_x,N_y,D,chi_try,d,Sys)

  alpha=1.0+interval
  if abs(norm_val)<threshold[0]:
   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten,N_x,N_y,D,chi_try,d,Sys)

 #PEPS_listten=All_dist(PEPS_listten,N_x, D)
 norm_val=Cal_norm(PEPS_listten,N_x,N_y,D,chi_try,d,Sys)
 print "Fixed norm", abs(norm_val)
 return PEPS_listten, abs(norm_val), count


def Cal_norm( PEPS_listten, N_x,N_y, D, chi_try, d,Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*N_x

 for Location in xrange(N_x-1):
  mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_try, N_y,Sys)

 Location=N_x-1
 mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_try, N_y,Sys, N_x)


 norm_val=mps_boundry_left[N_x-2].product_nonsymm(mps_boundry_right[N_x-1])

 return  norm_val






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
 if ( Max_val < 0.2e-1) or (Max_val > 0.2e+1)   :

  if Max_val >= 1:
   a=a*(1.00/Max_val)
  if Max_val < 1: 
   a=a*(1.00/Max_val)
 else: a=a;
 return a

def All_dist(PEPS_listten,N_x, D):

 for i in xrange(N_x-1):
  for j in xrange(N_x):
   PEPS_listten[i][j], PEPS_listten[i+1][j]=equall_dis_H(PEPS_listten[i][j], PEPS_listten[i+1][j],D)
   #PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   #PEPS_listten[i+1][j]=max_ten(PEPS_listten[i+1][j])
 
 for i in xrange(N_x):
  for j in xrange(N_x-1):
   PEPS_listten[i][j], PEPS_listten[i][j+1]=equall_dis_V(PEPS_listten[i][j], PEPS_listten[i][j+1],D)
   #PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   #PEPS_listten[i][j+1]=max_ten(PEPS_listten[i][j+1])
 
 return PEPS_listten



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


def condition_number(N):
 N_mat=N.getBlock()
 A_np=Mat_uni_to_np(N_mat)
 norm_val=npLA.norm(A_np)
 #print norm_val
 A_np=A_np*(1.0/norm_val)
 val=npLA.cond(A_np) 
 return val



def equall_dis_H(PEPS_1, PEPS_2, D):

 D_dim=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)


 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,5,3,4],3)
 
 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)

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

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 U.setLabel([4,6,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([4,6,-200],2)
 l=U*1.0
 qq.setLabel([-200,5,8,7])
 l.setLabel([4,6,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,6],2)
 U, V, s= TU.setTruncation(Teta, sum(D_dim))

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

 D_dim=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)


 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,4,3,5],3)
 
 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)

 
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


 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)



 U.setLabel([5,7,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([5,7,-200],2)
 l=U*1.0
 qq.setLabel([-200,6,8,9])
 l.setLabel([5,7,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,7],2)
 U, V, s= setTruncation3(Teta, sum(D_dim))

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

def  norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):
 val1=((((r_up*r_dp)*iden_h)*N_ten)*l_up)*l_dp
 #print val1
 #val0=((((r_u*r_d)*H)*N_ten)*l_u)*l_d


 val2=((((r_up*r_d)*Ham)*N_ten)*l_up)*l_d
 val3=((((r_u*r_dp)*Ham)*N_ten)*l_u)*l_dp

 #print "2, 3", val2, val3
 
 return val1[0]-val3[0]-val2[0]#+val0[0]






########@profile
def  optimum_0(N_ten, l_u, r_u, r_d, l_d , l_up, r_up,l_dp, r_dp, Ham,iden_h,H):
  

 #l_up.setLabel()
 
 Env_r=((l_dp)*(l_up*iden_h))*N_ten
 Env_r.permute([ -200,-3,-10, -100,3,10],3)

 Env_r1=copy.copy(Env_r)
 Env_r1.transpose()
 Env_r=Env_r+Env_r1



 Env_s=(r_u*N_ten)*((l_u*Ham)*l_dp)
 Env_s.permute([-200,-3,-10],3)

 Env_s1=((((r_d)*Ham)*N_ten)*l_up)*l_d

 Env_s1.permute([-100,3,10],0)
 Env_s1.transpose()
 Env_s=Env_s+Env_s1




 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([-100,3,10,-200,-3,-10])

 
 r_up=A2_inv*Env_s
 r_up.permute([-100,3,10],3)
 r_up.setLabel([-100,3,10])

 r_dp=r_up*1.0
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])
 return r_up, r_dp


########@profile
def  optimum_1(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):
  
 Env_r=((N_ten)*(r_up*r_dp))*iden_h
 Env_r.permute([-10,-13,-400,10,13,-300],3)

 Env_r1=copy.copy(Env_r)
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=(((r_dp*r_u)*N_ten)*(l_u))*Ham
 Env_s.permute([-10,-13,-400],3)

 Env_s1=((((r_up*r_d)*Ham)*N_ten))*l_d
 Env_s1.permute([10,13,-300],0)
 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)

 U.transpose()
 V.transpose()
 S=inverse(S)


 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([10,13,-300,-10,-13,-400])


 l_up=A2_inv*Env_s
 l_up.permute([10,13,-300],3)

 l_dp=l_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])

 return l_up, l_dp


def  optimum_11(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):

 Env_r=((N_ten)*(r_up*r_dp))*iden_h
 Env_r.permute([-10,13,-300,10,-13,-400],3)

 Env_s=(((r_dp*r_u)*N_ten)*(l_u))*Ham
 Env_s.permute([-10,13,-300],2)

 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)

# U, S, V=svd_parityrl(Env_r)

 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([10,-13,-400,-10,13,-300])

 l_up=A2_inv*Env_s
 l_up.permute([10,-13,-400],3)
 l_dp=copy.copy(l_up)
 l_dp.transpose()
 l_dp.setLabel([-10,13,-300])

 return l_up, l_dp


######@profile
def  update_twotensor( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, U_ham0, U_ham, U_hamN, H0, H, HN , D, E_coulmn,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[2*i+1]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[2*i+1]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])

   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*mps_boundry_left[2*i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*mps_boundry_right[2*i]
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])
   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_left[2*i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=mps_boundry_right[2*i]*E_list_up[i]
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)


# for i in xrange(len(PEPS_listten)-1):
#  E_list_down[i].setLabel([1,2,3,4])
#  E_list_up[i+1].setLabel([1,2,3,4])
#  A=E_list_down[i]*E_list_up[i+1]
#  print "Inside", A[0]


 for i in xrange(N_x-1):
  #print "i", i
  Ham_u=1
  if i==0:
    Ham_u= U_ham0
  elif i==(N_x-2):
    Ham_u= U_hamN
  else:
   Ham_u= U_ham

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, Ham_u, H_orig,D,threshold, interval,Sys)
  if Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, Ham_u, H_orig,D,threshold, interval,Sys)




  E_coulmn.append(E_val)

  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])

   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])

   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)





######@profile
def  Update_twotensor_local(PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],3)


 r_u.setLabel([-100,3,10])

 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,9,200])
 q_d.permute([6,7,9,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-9,-200,9,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-9,-200],4)

#####################################################

 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])


 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,11,14,15])
 qq_d.permute([400,11,14,15],4)
 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel( [-400, -11, 400, 11] )
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-11,-14,-15],4)



######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([90,200,9,-200])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])

######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=(A*Swap1)*q_d
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap3*q)
 A=A*mps_boundry_right[2*Location+1]


 B=E_right*mps_boundry_left[2*Location+3]
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[2*Location+3]

 B=B*mps_boundry_left[2*Location+2]
 B=B*(Swap4*qq_d)
 B=B*mps_boundry_right[2*Location+2]
 
 N_ten=A*B
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 N_ten=N_Positiv(N_ten)

######################################################

###############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 D_dim=l_u.bond(0).dim()

 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)


 
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
###########################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])
 Ham.setLabel([-3,-13,3,13])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count


 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
# print "E_2", h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
  r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)



 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)

 A=PEPS_1*1.0
 A.permute([6,7,3,9,10],5)

 Swap1=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap1.setLabel([-3,-9,3,9])
 A=Swap1*A
 A.permute([6,7,-3,-9,10],3)
 PEPS_1=A*1.0


 PEPS_2=l_up*qq
 PEPS_2.permute([11,10,13,14,15],3)


 

 return  PEPS_1,  PEPS_2,  h_h[0]/Norm_h[0]




def norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

#  A=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])
#  A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
#  A=A*mps_boundry_right[2*Location]
#  A=A*mps_boundry_right[2*Location+1]
# 
#  A=((A*Swap4)*Swap5)*(Peps_1*Swap6)
#  A=A*H
# 
#  B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right
# 
#  B=(B*Swap12)*((Peps_2*Swap11)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
#  B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
# 
#  val=A*B

###############
 A=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 T=Peps_1pd*Swap3p


 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]


 A=((A*Swap4p)*Swap5)*(Peps_1p*Swap6)
 A=A*iden_h

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 B=(B*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))

 #print A.printDiagram(), B.printDiagram()
 val1=A*B

 A=E_left*mps_boundry_left[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=((A*Swap4p)*Swap5)*(Peps_1*Swap6)
 A=A*Ham

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right


 #print Swap11.printDiagram(), Swap11p.printDiagram()
 B=B*Swap12
 B1=(Peps_2*Swap11)
 B=B*B1
 B=B*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))

 val2=A*B

 A=E_left*mps_boundry_left[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*Ham

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 B=(B*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2]))
 B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))

 val3=A*B
 #print val, val1, val2, val3
 #if val1[0]<0: print  "Oh, norm<0"
 #return val1[0]-2.0*val2[0]
 #return val[0]+val1[0]-val2[0]-val3[0]
 return val1[0]-val2[0]-val3[0], val1[0]


#@profile

def optimum_0_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

 Env_r=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_r=(Env_r*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))


 Env_r=((Env_r*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_r=Env_r*iden_h

 Env_r=(Env_r*(mps_boundry_right[2*Location+1]*mps_boundry_right[2*Location]))

 Env_r=(Env_r*Swap4p)*(Swap3p)
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap6

 A=(mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1

 A=A*(E_left*Swap2)
 Env_r=Env_r*A
 Env_r.permute([ 53,56,39, 54,57,6,35,37,55,34],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s=(Env_s*Swap12)*(Peps_2p*(Swap11p*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])))

 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
 Env_s=Env_s*Ham
 Env_s=((Env_s))*(Swap6)

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*((Peps_1d*Swap2)*Swap3)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=(A*Swap4)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)



 Env_s1=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s1=(Env_s1*Swap12)*((Peps_2)*((mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])*Swap11))

 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_s1=Env_s1*Ham
 Env_s1=((Env_s1)*(Swap4p*Swap5))

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*(Peps_1*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=Env_s1*((Swap2)*Swap3p)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,35,37,55,34,53,56,39,54,57])

 
 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,35,37,55,34],3)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.transpose()
 Peps_1pd.setLabel([54,57,53,56,39])

 return Peps_1p, Peps_1pd



def optimum_1_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

 Env_r=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_r=Env_r*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_r=Env_r*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])


 Env_r=((Env_r*Swap4p)*Swap5)*(Peps_1p*Swap6)
 Env_r=Env_r*iden_h
 Env_r=((Env_r*Swap10)*Swap11p)
 Env_r=((Env_r*Swap9)*Swap8p)
 Env_r=Env_r*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])

 A=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12
 A=A*Swap7
 A=A*E_right

 Env_r=Env_r*A
 Env_r.permute([ 48,58,46, 49,57,11,31,32,33,15],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_s=Env_s*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 Env_s=((Env_s*Swap4p)*Swap5)*(Peps_1*Swap6)
 Env_s=Env_s*Ham
 Env_s=((Env_s*Swap10)*Swap11)
 Env_s=(Env_s*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))*Peps_2

 Env_s=Env_s*E_right
 Env_s=Env_s*((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)

 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8p)
 Env_s.permute([48,58,46,49,57],5)



 Env_s1=E_left*mps_boundry_left[2*Location]
 Env_s1=Env_s1*mps_boundry_left[2*Location+1]
 Env_s1=Env_s1*(Swap1*((Peps_1d*Swap3)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=((Env_s1*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*Ham
 A=(((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)*(Swap7*(Swap8*Peps_2d)))*E_right

 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11p)
 Env_s1=(Env_s1*mps_boundry_left[2*Location+2])
 Env_s1=(Env_s1*mps_boundry_left[2*Location+3])

 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([11,31,32,33,15,48,58,46, 49,57])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([11,31,32,33,15],3)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.transpose()
 Peps_2pd.setLabel([ 49,57,48,58,46])

 return Peps_2p, Peps_2pd





def Update_twotensor_local_inv(PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, Ham, H_orig, D, threshold, interval, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 Peps_1.setLabel([6,35,37,55,34])
 Peps_1.permute([6,35,37,55,34],3)
 Peps_2.setLabel([11,31,32,33,15])
 Peps_2.permute([11,31,32,33,15],3)

 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([54,57,53,56,39])

 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([49,57,48,58,46])



 Swap1=fermionicOPT(Sys, Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([70,36,35,-60])
 Swap2=Swap1*1.0
 Swap2.setLabel([56,53,-7,36])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([40,-9,57,54])
 Swap4=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap4.setLabel([41,58,40,90])
 Swap6=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(0))
 Swap6.setLabel([42,55,37,38])
 Swap5=fermionicOPT(Sys,Peps_1d.bond(0), Peps_1.bond(2))
 Swap5.setLabel([38,39,58,44])


#####################################################


 Swap7=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap7.setLabel([-14,50,49,57])
 Swap8=fermionicOPT(Sys, Peps_2.bond(0), Peps_2.bond(1))
 Swap8.setLabel([48,58,47,41])
 Swap9=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap9.setLabel([47,46,51,45])
 Swap10=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap10.setLabel([51,43,52,32])
 Swap11=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(1))
 Swap11.setLabel([52,34,-110,31])
 Swap12=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap12.setLabel([33,-150,140,50])



######################################################################



######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])

######################################################
 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1c", Norm_h[0], h_h[0]/Norm_h[0]





 if Sys[5]=="Previous":
  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0
  Peps_1pd=Peps_1d*1.0
  Peps_2pd=Peps_2d*1.0

  Swap8p=Swap8*1.0
  Swap11p=Swap11*1.0
  Swap3p=Swap3*1.0
  Swap4p=Swap4*1.0




 if  Sys[5]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   A.setLabel([ 1, 2, 3, 4, 5])
   A.permute([ 1, 2, 3, 4, 5], 6)
   Swap=fermionicOPT(Sys,A.bond(2), A.bond(3))
   Swap.setLabel([-3,-4,3,4])
   A=Swap*A
   A.permute([ 1, 2, -3, -4, 5], 3)
   A.setLabel([ 1, 2, 3, 4, 5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,4,6,12,9,10],4)

   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,4,-100])
   V.setLabel([-100,6,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   S.setLabel([100,-100])
   V=V*S
   V.permute([6,100,12,9,10],3)
   V.setLabel([11,31,32,33,15])

   U.permute([1,2,11,4,100],4)
   Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
   Swap.setLabel([-11,-4,11,4])
   U=Swap*U
   U.permute([ 1, 2, -11, -4, 100], 3)
   U.setLabel([6,35,37,55,34])

   Peps_1p=U*1.0
   Peps_2p=V*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([54,57,53,56,39])

   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([49,57,48,58,46])

   Swap8p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap8p.setLabel([48,58,47,41])
   Swap11p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap11p.setLabel([52,34,-110,31])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([40,-9,57,54])
   Swap4p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap4p.setLabel([41,58,40,90])







 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([44,45,42,43])
 Ham.setLabel([44,45,42,43])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])

 Peps_1p_init=Peps_1p*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0
 Peps_1pd_init=Peps_1pd*1.0

 valf, E_pos=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(Sys[3]):
  count=count+1
  E_2=E_1*1.0
  val, E_pos=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)

  Peps_2p, Peps_2pd=optimum_1_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count
 #print "valf", valf
 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count








######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=mps_boundry_left[Location*2+1]*A
 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4p)*Swap5)*(Peps_1p*Swap6)


 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right


 B=(B*Swap12)*(Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
 B=((B*Swap10)*Swap11p)*(Swap9*((Peps_2pd*Swap7)*Swap8p))


 N_ten=A*B

 H_orig.setLabel([44,45,42,43])
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2c", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
  Peps_1p, Peps_2p=renormalize_tensor_c_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5,Swap6,Swap7,Swap8p,Swap9,Swap10,Swap11p,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]



def norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3p*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10p
 A=A*iden_h

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2p*Swap11)

 val1=A*B

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10p
 A=A*Ham

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val2=A*B

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_down[2*Location+1]

 A=A*((Swap3p*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2*Swap11)

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0], val1[0]









def optimum_0_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):


 Env_r=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_r=Swap7*Env_r
 Env_r=Env_r*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_r=Env_r*((Peps_2pd*Swap6)*Swap5p)
 Env_r=(Env_r*Swap4p)*Swap12
 Env_r=Env_r*(Peps_2p*Swap11)
 Env_r=Env_r*iden_h

 Env_r=Env_r*(Swap10p*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1]))

 Env_r=(Env_r)*(Swap9)
 Env_r=(Env_r)*(Swap8)
 Env_r=Env_r*Swap3p


 A=E_left*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])
 A=A*Swap1
 A=A*Swap2


 Env_r=(Env_r)*(A)


 Env_r.permute([ 32,33,37, 35,34,6,31,58,59,10],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s=Swap7*Env_s
 Env_s=Env_s*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s=Env_s*((Peps_2d*Swap6)*Swap5)
 Env_s=(Env_s*Swap4)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*Ham


 Env_s=(Env_s*(Swap8*Swap9))*Swap10p


 Env_s=(Env_s*((Peps_1d*Swap3)*Swap2))*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[2*Location]
 Env_s=Env_s*mps_boundry_up[2*Location+1]

 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5p)
 Env_s1=(Env_s1*Swap4p)*Swap12
 Env_s1=Env_s1*(Peps_2*Swap11)
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])*Swap10)

 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3p))
 Env_s1=Env_s1*Peps_1
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*(((mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,31,58,59,10,32,33,37,35,34])

 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,31,58,59,10],3)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.transpose()
 Peps_1pd.setLabel([35,34,32,33,37])

 return Peps_1p, Peps_1pd


#####@profile
def optimum_1_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):

 Env_r=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_r=Env_r*(Peps_1p*Swap1)
 Env_r=Env_r*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_r=Env_r*((Swap3p*Peps_1pd)*Swap2)
 Env_r=(Env_r*(Swap8*Swap9))*Swap10p
 Env_r=Env_r*iden_h

 Env_r=Env_r*Swap4p
 Env_r=Env_r*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_r=Env_r*Swap5p
 Env_r=Env_r*Swap12
 Env_r=Env_r*Swap11

 A=(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])*E_right

 A=Swap7*A
 A=Swap6*A

 Env_r=(Env_r*A)

 Env_r.permute([43,42,47,45,44,54,52,51,53,15],5)
 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[2*Location]
 Env_s=Env_s*mps_boundry_down[2*Location+1]
 Env_s=Env_s*((Swap3*Peps_1d)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10p
 Env_s=Env_s*Ham
 Env_s=((Swap11*Env_s)*Swap12)*Swap4
 Env_s=Env_s*mps_boundry_down[2*Location+2]
 Env_s=Env_s*mps_boundry_down[2*Location+3]
 Env_s=Env_s*((Swap5*Peps_2d)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*(mps_boundry_up[2*Location+3]*mps_boundry_up[2*Location+2])


 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s1=Env_s1*(Peps_1*Swap1)
 Env_s1=Env_s1*mps_boundry_down[2*Location]
 Env_s1=Env_s1*mps_boundry_down[2*Location+1]

 Env_s1=Env_s1*((Swap3p*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((Peps_2*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4p

 Env_s1=Env_s1*mps_boundry_up[2*Location+2]
 Env_s1=Env_s1*mps_boundry_up[2*Location+3]

 Env_s1=Env_s1*mps_boundry_down[2*Location+2]
 Env_s1=Env_s1*mps_boundry_down[2*Location+3]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5p)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([54,52,51,53,15,43,42,47,45,44])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([54,52,51,53,15],3)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.transpose()
 Peps_2pd.setLabel([45,44,43,42,47])

 return Peps_2p, Peps_2pd





def Update_twotensor_local_row_inv(PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


######################################################################


######################################################################

 Peps_1.setLabel([6,31,58,59,10])
 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([35,34,32,33,37])
 Peps_2.setLabel([54,52,51,53,15])
 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([45,44,43,42,47])


 Swap1=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([22,30,31,17])
 Swap2=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap2.setLabel([33,32,-7,30])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([38,36,34,35])
 Swap4=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap4.setLabel([23,41,48,36])
 Swap5=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap5.setLabel([42,43,-10,41])
 Swap6=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap6.setLabel([46,-14,44,45])

 Swap7=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap7.setLabel([24,53,46,20])

 Swap8=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap8.setLabel([37,39,40,38])

 Swap9=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap9.setLabel([56,57,58,39])

 Swap10=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap10.setLabel([25,59,57,54])

 Swap11=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap11.setLabel([55,50,51,52])

 Swap12=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap12.setLabel([47,48,49,50])



######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,266])
 mps_boundry_up[Location*2+1].setLabel([266,25,277])
 mps_boundry_up[Location*2+2].setLabel([277,15,288])
 mps_boundry_up[Location*2+3].setLabel([288,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,311])
 mps_boundry_down[Location*2+1].setLabel([311,-7,322])
 mps_boundry_down[Location*2+2].setLabel([322,23,333])
 mps_boundry_down[Location*2+3].setLabel([333,-10,19])

######################################################

# A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

# A=A*(Peps_1*Swap1)
# A=A*mps_boundry_down[2*Location]
# A=A*mps_boundry_down[2*Location+1]

# A=A*((Swap3*Peps_1d)*Swap2)
# A=(A*(Swap8*Swap9))*Swap10

# B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

# B=Swap7*B
# B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

# B=B*((Peps_2d*Swap6)*Swap5)
# B=(B*Swap4)*Swap12
# B=B*(Peps_2*Swap11)

# N_ten=A*B
# 
# N_ten.permute([40,49,56,55],2)



# H_orig.setLabel([40,49,56,55])
# iden_h=H_orig*1.0
# iden_h.identity()
# iden_h.setLabel([40,49,56,55])
# Norm_h=N_ten*iden_h
# h_h=N_ten*H_orig
# print "E_1r", Norm_h[0], h_h[0]/Norm_h[0]



# Peps_1p=Peps_1*1.0
# Peps_2p=Peps_2*1.0
# 
# Peps_1pd=Peps_1p*1.0
# Peps_1pd.setLabel([32,33,37,35,34])
# Peps_2pd=Peps_2p*1.0
# Peps_2pd.setLabel([43,42,47,45,44])





 if  Sys[5]=="Previous":


  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0

  Peps_1pd=Peps_1p*1.0
  Peps_1pd.transpose()
  Peps_1pd.setLabel([35,34,32,33,37])

  Peps_2pd=Peps_2p*1.0
  Peps_2pd.transpose()
  Peps_2pd.setLabel([45,44,43,42,47])


  Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap3p.setLabel([38,36,34,35])
  Swap4p=fermionicOPT(Sys,Peps_2pd.bond(1), Peps_2pd.bond(0))
  Swap4p.setLabel([23,41,48,36])


  Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
  Swap5p.setLabel([42,43,-10,41])

  Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap10p.setLabel([25,59,57,54])





 if  Sys[5]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   PEPS_B.permute([6,5,8,9,10],5)
   Swap=fermionicOPT(Sys,PEPS_B.bond(1), PEPS_B.bond(2))
   Swap.setLabel([-5,-8,5,8])
   PEPS_B=Swap*PEPS_B
   PEPS_B.permute([6,-5,-8,9,10],3)
   PEPS_B.setLabel([4,7,8,9,10])


   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,5,7,12,9,10],4)
   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,5,-100])
   V.setLabel([-100,7,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   U.permute([1,2,11,100,5],3)
   U.setLabel([6,31,58,59,10])

   S.setLabel([100,-100])
   V=V*S
   V.permute([100,7,12,9,10],5)
   Swap=fermionicOPT(Sys,V.bond(1), V.bond(2))
   Swap.setLabel([-7,-12,7,12])
   V=Swap*V
   V.permute([100,-7,-12,9,10],3)
   V.setLabel([54,52,51,53,15])


   Peps_1p=U*1.0
   Peps_2p=V*1.0
   Peps_1pd=Peps_1p*1.0
   Peps_2pd=Peps_2p*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([35,34,32,33,37])
   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([45,44,43,42,47])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([38,36,34,35])
   Swap4p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap4p.setLabel([23,41,48,36])


   Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap5p.setLabel([42,43,-10,41])

   Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap10p.setLabel([25,59,57,54])




 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([40,49,56,55])

 Ham.setLabel([40,49,56,55])
 
 iden_h=Ham*1.0
 iden_h.setLabel([40,49,56,55])
 iden_h.identity()

 iden_h.setLabel([40,49,56,55])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf, E_pos=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(Sys[3]):
  count=count+1
  E_2=E_1*1.0
  val, E_pos=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)

  Peps_2p, Peps_2pd=optimum_1_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)








 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
  print "warning_norm_in_optimization",  abs(valf), "count", count



 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3p*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10p

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2p*Swap11)

 N_ten=A*B


 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2r", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
  Peps_1p, Peps_2p=renormalize_tensor_r_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5p,Swap6,Swap7,Swap8,Swap9,Swap10p,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]


###@profile

def renormalize_tensor_r_single(Peps_1, Peps_1d,Peps_2, Peps_2d,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval):

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   Peps_1=Peps_1*alpha
   Peps_1d=Peps_1d*alpha
   Peps_2=Peps_2*alpha
   Peps_2d=Peps_2d*alpha
   A=E_left*mps_boundry_up[2*Location]
   A=A*mps_boundry_up[2*Location+1]

   A=A*(Peps_1*Swap1)
   A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

   A=A*((Swap3*Peps_1d)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10

   B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

   B=Swap7*B
   B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

   B=B*((Peps_2d*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2*Swap11)

   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   Peps_1=Peps_1*alpha
   Peps_1d=Peps_1d*alpha
   Peps_2=Peps_2*alpha
   Peps_2d=Peps_2d*alpha
   A=E_left*mps_boundry_up[2*Location]
   A=A*mps_boundry_up[2*Location+1]

   A=A*(Peps_1*Swap1)
   A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

   A=A*((Swap3*Peps_1d)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10

   B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

   B=Swap7*B
   B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

   B=B*((Peps_2d*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2*Swap11)

   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*mps_boundry_up[2*Location]
 A=A*mps_boundry_up[2*Location+1]

 A=A*(Peps_1*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1, Peps_2





###@profile
def renormalize_tensor_c_single(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval):

 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))


 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


   B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

   B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))


   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<threshold[0]:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


   B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

   B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))

   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right
 B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p









#######@profile
def update_twotensor_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, U_ham0, U_ham, U_hamN, H0, H, HN , D, E_coulmn,threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[2*i+1]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[2*i+1]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=((mps_boundry_up[2*i+1]*Swap1))*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):
  #print "i", i
  Ham_u=1
  if i==0:
    Ham_u= U_ham0
  elif i==(N_x-2):
    Ham_u= U_hamN
  else:
   Ham_u= U_ham

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, Ham_u, H_orig,D,threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, Ham_u, H_orig,D,threshold, interval,Sys)




  E_coulmn.append(E_val)

  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)


#######@profile
def  Update_twotensor_local_row(PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,26])
 mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location*2+2].setLabel([27,15,28])
 mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,31])
 mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location*2+2].setLabel([32,23,33])
 mps_boundry_down[Location*2+3].setLabel([33,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location*2]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_up[2*Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[2*Location+1]


 B=E_right*mps_boundry_up[2*Location+3]
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_down[2*Location+3]

 B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 N_ten=N_Positiv(N_ten)
######################################################



######################################################

##############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 
 D_dim=l_u.bond(0).dim()
 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)

 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
##########################################


 
# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])  #new
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])

 Ham.setLabel([-3,-13,3,13])
 
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12 or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
# print "E_2", h_h[0]/Norm_h[0]
# #print Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
   r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 r_up.setLabel([-100,3,9])   #new
 l_up.setLabel([9,13,-300])  #new


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)


 PEPS_2=l_up*qq
 A=PEPS_2*1.0
 A.permute([9,10,13,14,15],5)
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([9,-10,-13,14,15],3)
 PEPS_2=A*1.0
 PEPS_2.permute([9,-10,-13,14,15],3)


 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]




def renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval):
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 norm_val=Norm_h[0]
 #print "Zero_order_local",  norm_val

 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   r_up=r_up*alpha
   r_dp=r_dp*alpha
   l_up=l_up*alpha
   l_dp=l_dp*alpha
   Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
   norm_val=Norm_h[0]
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<threshold[0]:
   r_up=r_up*alpha
   r_dp=r_dp*alpha
   l_up=l_up*alpha
   l_dp=l_dp*alpha
   Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
   norm_val=Norm_h[0]

 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 norm_val=Norm_h[0]
 #print "Fixednorm_local", abs(norm_val)
 return r_up, l_up

 










def   rotate_all(PEPS_listten, N_x):

 PEPS_listten_t=[None]*N_x
 for i in xrange(N_x):
  PEPS_listten_t[i]=[None]*N_x

 for i in xrange(N_x):
  for j in xrange(N_x):
   ten=PEPS_listten[i][j]*1.0
   ten.setLabel([10,1,2,3,4])
   ten.permute([10,1,2,3,4],5)
   Swap1=fermionicOPT(Sys,ten.bond(0), ten.bond(1))
   Swap1.setLabel([-10,-1,10,1])
   #Swap1.identity()

   Swap2=fermionicOPT(Sys,ten.bond(3), ten.bond(4))
   Swap2.setLabel([-3,-4,3,4])
   #Swap2.identity()

   ten=(ten*Swap1)*Swap2
   ten.permute([-1,-10,2,-4,-3],3)
   ten.setLabel([1,2,3,4,5])
   PEPS_listten_t[j][i]=ten*1.0
 return PEPS_listten_t






def sqrt_general(N2):
  N_init=copy.copy(N2)
  blk_qnums = N2.blockQnum()
  for qnum in blk_qnums:
   M=N2.getBlock(qnum)
   eig=M.eigh()
   
   e=Sqrt_mat(eig[0])
   U_trans=copy.copy(eig[1])
   U_trans.transpose()
   M=U_trans*e*eig[1]
   N_init.putBlock(qnum,M)
  return N_init
def Sqrt_mat(e):
 d=int(e.row())
 
 for q in xrange(d):
   ##print e[q] 
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e  
 
def N_Positiv(N):
 N.setLabel([-200,-400,-100,-300])
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 N1=copy.copy(N)
 N1.setLabel([2,3,-100,-300])
 N.setLabel([-200,-400,2,3])

 N=N*N1
 N.permute([-200,-400,-100,-300],2)
 N_final=sqrt_general(N)
 N_final.setLabel([-200,-400,-100,-300])
 N_final.permute([-200,-400,-100,-300], 2)
 return N_final              










############@profile
def update_energy_eff(PEPS_listmps, U_eval, U_eval0, U_evalN, d, D, Location, N_x, chi_express):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 N_y=PEPS_listmps.N
 B_list=[None]*N_y
 if Location == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
 elif Location==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A


 B_list_u=[None]*N_y
 for i in xrange(N_y):
  B_list[i].setLabel([1,3,2,0,4])
  B_list_u[i]=copy.copy(B_list[i])
  B_list_u[i].permute([1,3,2,0,4],6)
  B_list_u[i].combineBond([3,2])
  B_list_u[i].combineBond([3,0])
  B_list_u[i].permute([1,3,4],2)

 list_bond=[]
 for q in xrange(len(B_list_u)):
   list_bond.append(B_list_u[q].bond(2).dim())

 mps_b=MPSclass.MPS(B_list_u[1].bond(1).dim(),max(list_bond)*4,len(B_list_u))

 for i in xrange(len(B_list_u)):
   mps_b[i]=copy.copy(B_list_u[i])


 Norm_init=mps_b.norm()
 #print Norm_init

 for i in xrange(N_y-1):
  row=B_list_u[i].bond(0).dim()*B_list_u[i].bond(1).dim()
  colm=B_list_u[i].bond(2).dim()
  if (row<=colm):
   U,V,s=MPSclass.setTruncation1(B_list_u[i],row)
  else:
   U,V,s=MPSclass.setTruncation1(B_list_u[i],colm)
  B_list_u[i]=copy.copy(U)
  s.setLabel([1,2])
  V.setLabel([2,3])
  V=s*V
  V.permute([1,3],1)
  B_list_u[i+1].setLabel([3,2,4])
  B_list_u[i+1]=V*B_list_u[i+1]
  B_list_u[i+1].permute([1,2,4],2)
  #A_l_can[i+1]=copy.copy(A_l[i+1])

 test_norm=B_list_u[N_y-1]*B_list_u[N_y-1]
 test_norm1=B_list_u[N_y-2]*B_list_u[N_y-2]

 #print test_norm[0]#,test_norm1[0]

#U_eval, U_eval0, U_evalN

 for i in xrange(N_y-1):
  A=uni10.UniTensor([ B_list_u[N_y-1-i].bond(0), B_list[N_y-1-i].bond(1), B_list[N_y-1-i].bond(2), B_list[N_y-1-i].bond(3), B_list_u[N_y-1-i].bond(2)])
  A.putBlock(B_list_u[N_y-1-i].getBlock())
  A.setLabel([1,3,2,0,4])

  B=uni10.UniTensor([ B_list_u[N_y-2-i].bond(0), B_list[N_y-2-i].bond(1), B_list[N_y-2-i].bond(2), B_list[N_y-2-i].bond(3), B_list_u[N_y-2-i].bond(2)])
  B.putBlock(B_list_u[N_y-2-i].getBlock())
  B.setLabel([-11,-3,-2,-10,1])

  if i==0:
    H_eff=U_evalN
    H_eff.setLabel([20,30,-2,2])
  elif i==N_y-2:
    H_eff=U_eval0
    H_eff.setLabel([20,30,-2,2])
  else:
    H_eff=U_eval
    H_eff.setLabel([20,30,-2,2])

  B.setLabel([-11,-3,-2,-10,1])
  B.permute([-11,-3,-10,-2,1],3)
  q,s,V=svd_parity5(B)
  s.setLabel([0,10])
  V.setLabel([0,2,3])
  r_u=V*s
  r_u.permute([10,2,3],2)
  r_u.setLabel([10,-2,1])
  
  q.setLabel([-11,-3,-10,10])
  r_u.setLabel([10,-2,1])

  A.setLabel([1,3,2,0,4])
  A.permute([1,2,3,0,4],2)
  U,s,qq=svd_parity6(A)
  s.setLabel([0,20])
  U.setLabel([1,2,0])
  l_u=U*s
  l_u.permute([1,2,20],2)

  qq.setLabel([20,3,0,4])
  l_u.setLabel([1,2,20])

 ##############QR_update###########################
  A=copy.copy(r_u)
  A.setLabel([-10,2,3])
  B=copy.copy(l_u)
  B.setLabel([3,6,-30])
  H_eff.setLabel([-2,-6,2,6])

  Teta=(A*B)*H_eff
  Teta.permute([-10,-2,-6,-30],2)
  Teta_mat=Teta.getBlock()
  col1=Teta_mat.col()
  row1=Teta_mat.row()
  dim_mat=1
  if col1<row1:
   dim_mat=col1
  else:
   dim_mat=row1
  #print "Info", dim_mat, chi_express
  if chi_express<dim_mat:
   U, V, S= setTruncation3(Teta, chi_express)
  else:
   U, V, S= setTruncation3(Teta, dim_mat)
  U.setLabel([-10,-2,-3])
  V.setLabel([-3,-6,-30])
  S=Sqrt(S)
  S.setLabel([3,-3])
  r_up=U*S
  r_up.permute([-10,-2,3],2)
  r_up.setLabel([10,-2,1])
  l_up=V*S
  l_up.permute([3,-6,-30],2)
  l_up.setLabel([1,2,20])

  U=q*r_up
  U.permute([-11,-3,-2,-10,1],4)
  V=qq*l_up
  V.permute([1,3,2,0,4],4)

  U.setLabel([-11, -3, 20, -10, 1])
  U.permute([-11, -3, 20, -10, 1],4)
  U.combineBond([-3,20])
  U.combineBond([-3,-10])
  U.permute([-11, -3, 1],2)

  V.setLabel([1, 3, 30, 0, 4])
  V.combineBond([3,30])
  V.combineBond([3,0])
  V.permute([1, 3, 4],2)



#  Theta=(A*H_eff)*B
#  #print Theta.printDiagram()
#  Theta.permute([-11, -3, 20, -10, 3, 30, 0, 4],4)
#  U, V, S=setTruncation2(Theta, D)
#  U.setLabel([-11, -3, 20, -10, -1])
#  S.setLabel([-1,1])
#  U=U*S
#  V.setLabel([1, 3, 30, 0, 4])

##################
#  Theta1=A*B
#  Theta2=U*V
#  Theta2.permute([-11, -3, 20, -10, 3, 30, 0, 4],4)
#  Theta2.setLabel([-11, -3, -2, -10, 3, 2, 0, 4])
##################

#  U.permute([-11, -3, 20, -10, 1],4)
#  U.combineBond([-3,20])
#  U.combineBond([-3,-10])
#  U.permute([-11, -3, 1],2)

#  V.combineBond([3,30])
#  V.combineBond([3,0])
#  V.permute([1, 3, 4],2)

#################################
#  norm_test1=Theta1*Theta1
#  norm_test2=Theta2*Theta2
#  norm_testf=Theta1*Theta2
#  print "i", i, norm_test1[0], norm_test2[0], norm_testf[0]/((norm_test2[0]*norm_test1[0])**(0.5))
##################################################
  B_list_u[N_y-1-i]=copy.copy(V)
  B_list_u[N_y-2-i]=copy.copy(U)

  B_list_u[N_y-1-i].setLabel([1,2,3])
  B_list_u[N_y-1-i].permute([2,3,1],2)


  row=B_list_u[N_y-1-i].bond(0).dim()*B_list_u[N_y-1-i].bond(1).dim()
  colm=B_list_u[N_y-1-i].bond(2).dim()
 
  if (row<=colm):
   U,V,s=MPSclass.setTruncation1(B_list_u[N_y-1-i],row)
  else:
   U,V,s=MPSclass.setTruncation1(B_list_u[N_y-1-i],colm)


  U.setLabel([2,3,1])
  U.permute([1,2,3],2)
  B_list_u[N_y-1-i]=copy.copy(U)
  s.setLabel([3,4])
  V.setLabel([4,5])
  V_f=(s*V)
  V_f.permute([3,5],1)
  B_list_u[N_y-2-i].setLabel([-1,2,5])
  B_list_u[N_y-2-i]=V_f*B_list_u[N_y-2-i]
  B_list_u[N_y-2-i].permute([-1,2,3],2)


 list_bond=[]
 for  q  in  xrange(len(B_list_u)):
   list_bond.append(B_list_u[q].bond(2).dim())

 mps_b=MPSclass.MPS( B_list_u[1].bond(1).dim(), max(list_bond), len(B_list_u))

 for  i  in  xrange(len(B_list_u)):
   mps_b[i]=copy.copy(B_list_u[i])


 #mps_b=mps_b.normalize()

 return mps_b


def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 


def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat



def C_i(dim):
 Mat=uni10.Matrix(dim, dim)
 Mat.set_zero()
 for i in xrange(dim):
  for j in xrange(dim):
   if j==i+1:
    Mat[i*dim+j]=j**(0.5) 
 iden_m=Mat*1.0
 iden_m.identity()
 Mat2=Mat*1.0
 Mat2.transpose()
 #commutator=Mat2*Mat+(-1.0)*Mat*Mat2
 #print Mat, Mat2, Mat1 , commutator
 return Mat , Mat2, iden_m


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




def Ham_mat( J, g1, g2, J1, J2, h1, h2):
  dim=4

  Mat=uni10.Matrix(dim*dim, dim*dim)
  Mat.set_zero()
  for i in xrange(dim):
   for j in xrange(dim):
    for m in xrange(dim):
     for n in xrange(dim):
 
      if i==0 and j==1  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==2 and j==3  and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==0 and j==1  and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==2 and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J


      if i==1 and j==0  and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==2  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==1 and j==0  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==2 and j==3  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==0 and j==2 and m==2 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =+J
      if i==2 and j==0  and m==0 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==0 and j==3  and m==3 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==3 and j==0  and m==0 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==1 and j==2  and m==2 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==2 and j==1  and m==1 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J

      if i==1 and j==3  and m==3 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==1  and m==1 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J


      if i==0 and j==1  and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+(2*J2)
      if i==0 and j==2 and m==0 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =+J2+h2
      if i==0 and j==3  and m==0 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J2-h2
      if i==1 and j==0  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)
      if i==1 and j==1  and m==1 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+g2+(2.0*J1)+(2.0*J2)
      if i==1 and j==2  and m==1 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)+J2+h2
      if i==1 and j==3  and m==1 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)+J2-h2
      if i==2 and j==0  and m==2 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J1+h1
      if i==2 and j==1  and m==2 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+J1+(2*J2)+h1
      if i==2 and j==2  and m==2 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2+h1+h2
      if i==2 and j==3  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2+h1-h2
      if i==3 and j==0  and m==3 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J1-h1
      if i==3 and j==1  and m==3 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+J1+(2*J2)-h1
      if i==3 and j==2 and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2-h1+h2
      if i==3 and j==3 and m==3 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2-h1-h2

  #print Mat
  #Mat_np=Mat_uni_to_np(Mat)
  #print Mat_np
  return Mat



# def Ham_mat( J, g1, g2, J1, J2, h1, h2):
#  dim=2
#  c_iu, c_iu_dag, iden=C_i(dim)
#  ham=J*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
#  #print c_iu, c_iu_dag
#  ham=g1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
#  ham=g2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=-1.0*h2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=-1.0*h1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
#  ham=1.0*J2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=1.0*J1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
# 
#  return ham



def Heisenberg(h, d_phys, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 H_c = uni10.UniTensor([bdi, bdo], "Heisenberg")

 sx = matSx()
 sy = matSy()
 sz = matSz()
 idenS = matIden()


# h_ham=uni10.otimes(sx,sx)*(1)+uni10.otimes(sz,sz)*(1)+uni10.otimes(sy,sy)*(-1)
# h_ham=0.25*h_ham
# D=h_ham.eigh()
# print D[0][0]


# h_ham=uni10.otimes(uni10.otimes(sx,sx),iden)+uni10.otimes(uni10.otimes(sz,sz),iden)+(uni10.otimes(uni10.otimes(sy,sy),iden)*(-1))+uni10.otimes(iden,uni10.otimes(sx,sx))+uni10.otimes(iden,uni10.otimes(sz,sz))+(uni10.otimes(iden,uni10.otimes(sy,sy))*(-1))

# h_ham=0.25*h_ham
# D=h_ham.eigh()
# print D[0][0]


# h_ham=uni10.otimes(uni10.otimes(sx,sx),uni10.otimes(iden,iden))+uni10.otimes(uni10.otimes(sz,sz),uni10.otimes(iden,iden))+(uni10.otimes(uni10.otimes(sy,sy),uni10.otimes(iden,iden))*(-1))+uni10.otimes(uni10.otimes(iden,sx),uni10.otimes(sx,iden))+uni10.otimes(uni10.otimes(iden,sz),uni10.otimes(sz,iden))+(uni10.otimes(uni10.otimes(iden,sy),uni10.otimes(sy,iden))*(-1))+uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(sx,sx))+uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(sz,sz))+(uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(sy,sy))*(-1))

# h_ham=0.25*h_ham
# D=h_ham.eigh()
# print D[0][0]


# h_ham=(uni10.otimes(uni10.otimes(sx,sx),uni10.otimes(iden,iden))*(-1))+(uni10.otimes(uni10.otimes(iden,sx),uni10.otimes(sx,iden))*(-1))+(uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(sx,sx))*(-1))+(uni10.otimes(uni10.otimes(sz,iden),uni10.otimes(iden,iden))*(-1))+(uni10.otimes(uni10.otimes(iden,sz),uni10.otimes(iden,iden))*(-1))+(uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(sz,iden))*(-1))+(uni10.otimes(uni10.otimes(iden,iden),uni10.otimes(iden,sz))*(-1))

 #h_ham=0.25*h_ham
# D=h_ham.eigh()
# print D[0][0]

 # print d_phys
 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag

 #########

 #print uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu)
 #H.setRawElem(uni10.otimes(c_iu,c_iu_dag))
 #print H

# print c_id, c_id_dag, uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id)
# H.setRawElem(uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id))
# print H


 if Model[0]=="Fer_Z2" or Model[0]=="Fer_U1" or Model[0]=="Fer" :

#  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
#  ham=ham+h[0]*(uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id))
#  ham=ham+h[1]*0.25*(uni10.otimes(c_iu_dag*c_iu*c_id_dag*c_id,iden)+uni10.otimes(iden,c_iu_dag*c_iu*c_id_dag*c_id))
#  ham=ham+h[2]*0.25*(uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
#  ham=ham+h[2]*0.25*(uni10.otimes(c_id_dag*c_id,iden)+uni10.otimes(iden,c_id_dag*c_id))


  ham=Ham_mat( h[0], 0.25*h[1], 0.25*h[1], 0.25*h[2], 0.25*h[2],  0.25*h[3], 0.25*h[3])
  #ham=Ham_mat( h[0], 0.5*h[1], 0.5*h[1], 0.5*h[2], 0.5*h[2],  0, 0)

#2*2
  #ham=Ham_mat( h[0], 2.0*0.25*h[1], 0.25*2.0*h[1], 2.0*0.25*h[2], 0.25*2.0*h[2], 0.25*h[3], 0.5*h[3] )


  
  H.setRawElem(ham)
  #print "H", H
  if Model[1] is "on":
   H=switch_H(H)



 if Model[0]=="ITF_Z2"  or Model[0]=="ITF":
  ham =h[0]*uni10.otimes(sx,sx)*(-1)+(-0.25)*h[1]*(uni10.otimes(idenS,sz)+uni10.otimes(sz,idenS))
  H.setRawElem(ham)
 if Model[0]=="Heis_Z2" or Model[0]=="Heis_U1" or Model[0]=="Heis":
  ham =(h[0]*uni10.otimes(sz,sz)+h[1]*uni10.otimes(sx,sx)+(-1.0*h[1])*uni10.otimes(sy,sy))*(0.25)
  H.setRawElem(ham)

 if Model[0]=="FFI_Z2" or Model[0]=="FFI":
  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
  ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
  H.setRawElem(ham)


 return H



def Heisenberg0(h, d_phys, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 sx = matSx()
 sy = matSy()
 sz = matSz()
 idenS = matIden()


 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag

#########

 if Model[0]=="Fer_Z2" or  Model[0]=="Fer_U1" or Model[0]=="Fer":

#  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
#  ham=ham+h[0]*(uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id))
#  ham=ham+h[1]*0.25*(2.0*uni10.otimes(c_iu_dag*c_iu*c_id_dag*c_id,iden)+uni10.otimes(iden,c_iu_dag*c_iu*c_id_dag*c_id))
#  ham=ham+h[2]*0.25*(2.0*uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
#  ham=ham+h[2]*0.25*(2.0*uni10.otimes(c_id_dag*c_id,iden)+uni10.otimes(iden,c_id_dag*c_id))
  
  
#  ham=h[2]*0.5*(uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
#  ham=ham+h[2]*0.5*(uni10.otimes(c_id_dag*c_id,iden)+uni10.otimes(iden,c_id_dag*c_id))
  
  ham=Ham_mat( h[0], 0.25*2.0*h[1], 0.25*h[1], 0.25*2.0*h[2], 0.25*h[2],0.5*h[3], 0.25*h[3])
  #ham=Ham_mat( h[0], 0.5*2.0*h[1], 0.5*h[1], 0.5*2.0*h[2], 0.5*h[2],0, 0)

#2*2
  #ham=Ham_mat( h[0], 2.0*0.25*h[1], 0.25*2.0*h[1], 2.0*0.25*h[2], 0.25*2.0*h[2], 0.25*h[3], 0.5*h[3] )
  
  
  H.setRawElem(ham)
  #print "H", H

  if Model[1] is "on":
   H=switch_H(H)

 if Model[0]=="ITF_Z2" or Model[0]=="ITF" :
  ham =h[0]*uni10.otimes(sx,sx)*(-1)+(-0.25)*h[1]*(2.0*uni10.otimes(sz,idenS)+uni10.otimes(idenS,sz))
  H.setRawElem(ham)
 if Model[0]=="Heis_Z2" or Model[0]=="Heis_U1" or Model[0]=="Heis":
  ham =(h[0]*uni10.otimes(sz,sz)+h[1]*uni10.otimes(sx,sx)+(-1.0*h[1])*uni10.otimes(sy,sy))*(0.25)
  H.setRawElem(ham)

 if Model[0]=="FFI_Z2" or Model[0]=="FFI":
  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
  ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
  H.setRawElem(ham)

 return H



def HeisenbergN(h, d_phys, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 sx = matSx()
 sy = matSy()
 sz = matSz()
 idenS = matIden()

 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag

#########


 if Model[0]=="Fer_Z2" or Model[0]=="Fer_U1" or Model[0]=="Fer":

#  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
#  ham=ham+h[0]*(uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id))
#  ham=ham+h[1]*0.25*(uni10.otimes(c_iu_dag*c_iu*c_id_dag*c_id,iden)+2.0*uni10.otimes(iden,c_iu_dag*c_iu*c_id_dag*c_id))
#  ham=ham+h[2]*0.25*(uni10.otimes(c_iu_dag*c_iu,iden)+2.0*uni10.otimes(iden,c_iu_dag*c_iu))
#  ham=ham+h[2]*0.25*(uni10.otimes(c_id_dag*c_id,iden)+2.0*uni10.otimes(iden,c_id_dag*c_id))

#  ham=h[2]*0.0*(uni10.otimes(c_iu_dag*c_iu,iden)+2.0*uni10.otimes(iden,c_iu_dag*c_iu))
#  ham=ham+h[2]*0.0*(uni10.otimes(c_id_dag*c_id,iden)+2.0*uni10.otimes(iden,c_id_dag*c_id))

  ham=Ham_mat( h[0], 0.25*h[1], 0.25*2.0*h[1], 0.25*h[2], 0.25*2.0*h[2], 0.25*h[3], 0.5*h[3] )

#2*2
  #ham=Ham_mat( h[0], 2.0*0.25*h[1], 0.25*2.0*h[1], 2.0*0.25*h[2], 0.25*2.0*h[2], 0.25*h[3], 0.5*h[3] )


  #ham=Ham_mat( h[0], 0.5*h[1], 0.5*2.0*h[1], 0.5*h[2], 0.5*2.0*h[2], 0, 0 )

  H.setRawElem(ham)
  #print "H", H

  if Model[1] is "on":
   H=switch_H(H)


 if Model[0]=="ITF_Z2" or Model[0]=="ITF":
  ham =h[0]*uni10.otimes(sx,sx)*(-1)+(-0.25)*h[1]*(2.0*uni10.otimes(idenS,sz)+uni10.otimes(sz,idenS))
  H.setRawElem(ham)
  #print "N", H
 if Model[0]=="Heis":
  ham =(h[0]*uni10.otimes(sz,sz)+h[1]*uni10.otimes(sx,sx)+(-1.0*h[1])*uni10.otimes(sy,sy))*(0.25)
  H.putBlock(ham)
 if Model[0]=="Heis_Z2" or Model[0]=="Heis_U1" or Model[0]=="Heis":
  ham =(h[0]*uni10.otimes(sz,sz)+h[1]*uni10.otimes(sx,sx)+(-1.0*h[1])*uni10.otimes(sy,sy))*(0.25)
  H.setRawElem(ham)

 if Model[0]=="FFI_Z2"or Model[0]=="FFI":
  ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
  ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
  H.setRawElem(ham)

 return H



def switch_H(H_uni):

 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 H=H_uni*1.0
 M=H_uni.getBlock(q0_odd)
 H.putBlock(q0_even,M)
 M=H_uni.getBlock(q0_even)
 H.putBlock(q0_odd,M)

 return H



def N_total(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0.0, 0.0, 0.25*1.0, 0.25*1.0,0,0)

#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0.5*1.0, 0.5*1.0,0,0)

 H.setRawElem(ham)
 return H
def N_total0(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0.0, 0.0, 0.25*2.0*1.0, 0.25*1.0,0,0)


#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0.25*2.0, 0.25*2.0,0,0)


 #ham=Ham_mat( 0.0, 0.0, 0.0, 0.5*2.0*1.0, 0.5*1.0,0,0)
 H.setRawElem(ham)
 return H
def N_totalN(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0.0, 0.0, 0.25*1.0, 0.25*2.0*1.0,0,0)

#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0.25*2.0, 0.25*2.0,0,0)


 #ham=Ham_mat( 0.0, 0.0, 0.0, 0.5*1.0, 0.5*2.0*1.0,0,0)
 H.setRawElem(ham)
 return H



def Sz_total(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0.0, 0.0, 0,0,0.25*1.0, 0.25*1.0)

#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0,0,0.5, 0.5)

 H.setRawElem(ham)
 return H
def Sz_total0(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0.0, 0.0, 0,0,.25*2.*1., .25*1.)

#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0,0,.5, .5)

 H.setRawElem(ham)
 return H
def Sz_totalN(h, d_phys, Model):
 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
 ham=Ham_mat( 0.0, 0, 0,0,0,.25*1.0, .25*2.*1.)
 #ham=Ham_mat( 0.0, 0, 0,0,0,.5*1.0, .5*2.*1.)

#2*2
 #ham=Ham_mat( 0.0, 0.0, 0.0, 0,0,.5, .5)

 H.setRawElem(ham)
 return H



#def N_total_H(h, d_phys, Model):

# bdi = uni10.Bond(uni10.BD_IN, d_phys)
# bdo = uni10.Bond(uni10.BD_OUT, d_phys)
# H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
# sx = matSx()
# sy = matSy()
# sz = matSz()
# iden = matIden()
# c_i, c_i_dag, iden=C_i(len(d_phys))
# ham =0.5*(uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))

## ham =0.25*(uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))

# H.setRawElem(ham)
# #if Model[1] is "on":
#  #H=switch_H(H)

# return H



#def N_total0_H(h, d_phys, Model):

# bdi = uni10.Bond(uni10.BD_IN, d_phys)
# bdo = uni10.Bond(uni10.BD_OUT, d_phys)
# H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
# sx = matSx()
# sy = matSy()
# sz = matSz()
# iden = matIden()
# c_i, c_i_dag, iden=C_i(len(d_phys))
# ham =0.50*(uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))

## ham =0.25*(2.0*uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))

# H.setRawElem(ham)
# #if Model[1] is "on":
#  #H=switch_H(H)
# return H



#def N_totalN_H(h, d_phys, Model):
# bdi = uni10.Bond(uni10.BD_IN, d_phys)
# bdo = uni10.Bond(uni10.BD_OUT, d_phys)
# H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heisenberg")
# sx = matSx()
# sy = matSy()
# sz = matSz()
# iden = matIden()
# c_i, c_i_dag, iden=C_i(len(d_phys))

# ham =0.5*(uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))

## ham =0.25*(uni10.otimes(c_i_dag*c_i,iden)+2.0*uni10.otimes(iden,c_i_dag*c_i))

# H.setRawElem(ham)
# #if Model[1] is "on":
#  #H=switch_H(H)

# return H




def rotate(PEPS_list):
 PEPS_list_copy=[copy.copy(PEPS_list[i]) for i in xrange(len(PEPS_list))]
 for i in xrange(len(PEPS_list)):
  PEPS_list_copy[i].setLabel([0,1,2,3,4])

  ten=PEPS_list_copy[i]*1.0
  ten.setLabel([0,1,2,3,4])
  ten.permute([0,1,2,3,4],5)

  Swap1=fermionicOPT(Sys,ten.bond(0), ten.bond(1))
  Swap1.setLabel([0,-1,5,1])

  Swap2=fermionicOPT(Sys,ten.bond(0), ten.bond(2))
  Swap2.setLabel([5,-2,6,2])

  Swap3=fermionicOPT(Sys,ten.bond(0), ten.bond(3))
  Swap3.setLabel([6,3,7,-3])

  ten1=((ten*Swap1)*Swap2)*Swap3


  Swap1=fermionicOPT(Sys,ten.bond(3), ten.bond(2))
  Swap1.setLabel([-3,9,11,-2])

  Swap2=fermionicOPT(Sys,ten.bond(3), ten.bond(1))
  Swap2.setLabel([11,8,12,-1])
  #print t1.printDiagram()
  ten=(ten1*Swap1)*Swap2
  ten.permute([12,8,9,7,4],3)
  ten.setLabel([1,2,3,4,5])
  PEPS_list_copy[i]=ten*1.0

  #PEPS_list_copy[i].permute([3,1,2,0,4],3)
  #PEPS_list_copy[i].setLabel([0,1,2,3,4])
 return PEPS_list_copy
#####@profile
def make_Env_singleLayer( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_y,Sys):

 Peps_ket=[]
 for i in xrange( N_y ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange( N_y ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)


 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)


 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

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
  mps_list=[None]*(N_y*2)

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])

  Peps_ket[0].setLabel([0,-1,2,3,4])
  Peps_bra[0].setLabel([-10,11,2,-3,-4])

  Tem10.identity()
  Tem1.identity()
  Tem11.identity()

  mps_list[0]=(((Peps_ket[0]*Tem0)*Tem1)*((Peps_bra[0])*Tem11))*SwapX
  mps_list[0].permute([10, -3, 4, 3, -4], 2)
  ##########################################


  bdi = uni10.Bond(uni10.BD_IN, Peps_ket[0].bond(4).Qlist())
  bdo = uni10.Bond(uni10.BD_OUT, Peps_ket[0].bond(4).Qlist())
  IdenX=uni10.UniTensor([bdi, bdo])
  IdenX.setLabel([1,4])
  IdenX.identity()
  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])

  mps_list[1]=SwapX*IdenX
  mps_list[1].permute([4,3,-4,6,1,11],4)
  #mps_list.append(results)

  #print mps_list[0].printDiagram(), mps_list[1].printDiagram() 

  ####################Middle##############
  for  q  in  xrange(1,N_y):
   Peps_ket[q].setLabel([0,1,2,3,4])
   Peps_bra[q].setLabel([10,11,2,-3,-4])
   Tem0.identity()
   Tem0.setLabel([0])
   Tem10.identity()
   Tem10.setLabel([-10])
   #print "Hi4"

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([-1,-10,1,10])

   mps_list[2*q]=(((Peps_ket[q]*Tem0))*((Peps_bra[q])*Tem10))*SwapX
   mps_list[2*q].permute([-1,11,-3,4,3,-4],3)

   ######################################################


   bdi = uni10.Bond( uni10.BD_IN, Peps_ket[q].bond(4).Qlist())
   bdo = uni10.Bond( uni10.BD_OUT, Peps_ket[q].bond(4).Qlist())
   IdenX=uni10.UniTensor([bdi, bdo])
   IdenX.setLabel([1,4])
   IdenX.identity()

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])

   mps_list[2*q+1]=IdenX*SwapX
   mps_list[2*q+1].permute([4,3,-4,6,1,11],4)

  mps_boundry=MPSclass.MPS(2,2,len(mps_list))


  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0



  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()


  return mps_boundry

  ######################## Add next Layer ###########

 if Location!=0:
  mps_list=[None]*(N_y*2)

  ###############################First-Layer###############################

  mps_boundry[0].setLabel([1,2,3])
  Tem1.setLabel([4])
  SwapX=fermionicOPT(Sys,mps_boundry[0].bond(1), Peps_ket[0].bond(1))
  SwapX.setLabel([5,6,2,4])
  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([1,5,3,6],2)

  mps_boundry[1].setLabel([5,0,6])
  Peps_ket[0].setLabel([0,1,2,3,4])
  Tem1.setLabel([1])
  mps_list[1]=(mps_boundry[1]*(Peps_ket[0]))
  mps_list[1].permute([5,1,2,3,6,4],4)

  ###########################################
 
  #print mps_list[0].printDiagram(), mps_list[1].printDiagram()
   ###############################################################################

  for q in xrange(1,N_y):
   ###########################################
    SwapX=fermionicOPT(Sys,mps_boundry[2*q].bond(1), Peps_ket[q-1].bond(4))
    SwapX.setLabel([10,-4,0,4])
    mps_boundry[2*q].setLabel([6,0,-6])

    mps_list[2*q]=mps_boundry[2*q]*SwapX
    mps_list[2*q].permute([6,4,10,-6,-4],3)



    mps_boundry[2*q+1].setLabel([5,0,6])
    Peps_ket[q].setLabel([0,1,2,3,4])

    mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_ket[q]
    mps_list[2*q+1].permute([5,1,2,3,6,4],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
   mps_boundry[i]=mps_list[i]*1.0

  #print mps_list[2].printDiagram(), mps_list[3].printDiagram()

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD( sum(chi_boundry_list) )
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t


######################### Next Absorption ####################################
  for  q  in  xrange(0, mps_boundry.N, 2):
    mps_boundry[q].setLabel([10,3,2])
    Peps_bra[q/2].setLabel([3,11,1,4,-2])
    mps_list[q]=((mps_boundry[q]*Peps_bra[q/2]))
    mps_list[q].permute([10,11,4,2,1,-2],3)




    mps_boundry[q+1].setLabel([2,1,3,5])
    SwapX=fermionicOPT(Sys, Peps_bra[q/2].bond(4), mps_boundry[q+1].bond(2) )
    SwapX.setLabel([ -5, -3, -2, 3 ])
    mps_list[q+1]=mps_boundry[q+1]*SwapX
    mps_list[q+1].permute([2,1,-2,-3,5,-5],4)




  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
 
  return mps_boundry



#####@profile
def make_Env_singleLayer_right( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_y,Sys, N_Position):

 Peps_ket=[]
 for i in xrange( N_y ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)




 Peps_bra=[]
 for i in xrange(N_y):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)


 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)


 #############   Zero   ##############################################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

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

 if Location==N_Position-1:
  mps_list=[None]*(2*N_y)

  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])

  SwapXX=fermionicOPT(Sys, Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapXX.setLabel([5,6,3,-4])

  Peps_ket[0].setLabel([0,-1,2,3,4])
  Peps_bra[0].setLabel([-10,11,2,-3,-4])

  Tem10.identity()
  Tem1.identity()
  Tem11.identity()
  Tem10.setLabel([10])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  Temr=Tem0*1.0
  Temr.setLabel([-3])
  
  mps_list[0]=(Peps_bra[0]*Temr)*(SwapX*Tem1)
  mps_list[0].permute([11, 10, -4, 2, -1], 2)


  Temr.setLabel([5])
  mps_list[1]=(((Peps_ket[0]*Temr)*SwapXX))
  mps_list[1].permute([ -4, 2, -1,0,6,4], 4)


  ####################   Middle   ##########################################
  for  q  in  xrange(1,len(Peps_ket)):
   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])

   SwapXX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapXX.setLabel([5,6,3,-4])

   Peps_ket[q].setLabel([0,-1,2,3,4])
   Peps_bra[q].setLabel([-10,11,2,-3,-4])

   Tem10.identity()
   Tem1.identity()
   Tem11.identity()
   Tem10.setLabel([10])
   Tem1.setLabel([1])
   Tem0.setLabel([0])
   Temr=Tem0*1.0
   Temr.setLabel([-3])
   
   mps_list[2*q]=(Peps_bra[q]*Temr)*(SwapX)
   mps_list[2*q].permute([11,1, 10, -4, 2, -1], 3)

   Temr.setLabel([5])
   mps_list[2*q+1]=(((Peps_ket[q]*Temr)*SwapXX))
   mps_list[2*q+1].permute([ -4, 2, -1,0,6,4], 4)

  mps_boundry=MPSclass.MPS(2,2,len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=copy.copy(mps_list[i])

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  return mps_boundry

  ########################   Add next Layer   ################################

 if Location!=(N_Position-1):
  mps_list=[None]*(N_y*2)

  ###############################   First-Layer   ###############################

  mps_boundry[0].setLabel([5,-3,6])

  Peps_bra[0].setLabel([-10,11,2,-3,-4])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  Temr=Tem0*1.0
  Temr.setLabel([5])

  mps_list[0]=(Temr*mps_boundry[0])*(Peps_bra[0])
  mps_list[0].permute([11,-10,2,6,-4],3)


  mps_boundry[1].setLabel([7,5,8])

  SwapXX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapXX.setLabel([5,6,3,-4])

  mps_list[1]=SwapXX*mps_boundry[1]
  mps_list[1].permute([7,-4,3,8,6],3)

###########################################

  for q in xrange(1, len(Peps_ket)):
###########################################
   mps_boundry[2*q].setLabel([5,-3,6])

   Peps_bra[q].setLabel([-10,11,2,-3,-4])


   mps_list[2*q]=(mps_boundry[2*q])*((Peps_bra[q]))

   mps_list[2*q].permute([5,11,-10,2,6,-4],4)


   mps_boundry[2*q+1].setLabel([7,5,8])
   SwapXX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapXX.setLabel([5,6,3,-4])


   mps_list[2*q+1]=SwapXX*mps_boundry[2*q+1]
   mps_list[2*q+1].permute([7,-4,3,8,6],3)

  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t



  mps_boundry[0].setLabel([11,-10,2,6])
  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])
  Tem1.setLabel([1])

  mps_list[0]=(SwapX*Tem1)*mps_boundry[0]
  mps_list[0].permute([11,10,6,2,-1],2)

  mps_boundry[1].setLabel([7,3,8])
  Peps_ket[0].setLabel([0,-1,2,3,4])
  mps_list[1]=mps_boundry[1]*Peps_ket[0]
  mps_list[1].permute([7,2,-1,0,8,4],4)

######################### Next Absorption ####################################
  for  q  in  xrange(1, len(Peps_ket)):
   mps_boundry[2*q].setLabel([11,-10,2,6])
   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])
   Tem1.setLabel([1])

   mps_list[2*q]=SwapX*mps_boundry[2*q]
   mps_list[2*q].permute([11,1,10,6,2,-1],3)

   mps_boundry[2*q+1].setLabel([7,3,8])
   Peps_ket[q].setLabel([0,-1,2,3,4])
   mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_ket[q]
   mps_list[2*q+1].permute([7,2,-1,0,8,4],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
 
  return mps_boundry


#####@profile
def make_Env_singleLayer_down( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)




 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).Qlist()
  D=Peps_bra[0].bond(3).Qlist()
 else:
  DX=Peps_ket[0].bond(0).Qlist()
  D=Peps_bra[0].bond(0).Qlist()

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

 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

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

 Temr=uni10.UniTensor([bdi_1])
 Temr.identity()
 Temr.setLabel([10])

 Temrr=uni10.UniTensor([bdi_1])
 Temrr.identity()
 Temrr.setLabel([10])

 if Location==0:

  mps_list=[None]*(N_x*2)

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])
  Peps_ket[0].setLabel([0,-1,2,3,4])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  
  mps_list[0]=(Peps_ket[0]*Tem0)*(SwapX*Tem1)
  mps_list[0].permute([10, 4, -10,2,3], 2)
  ##########################################

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])
  Peps_bra[0].setLabel([-10,-11,2,-3,-4])
  Tem11.setLabel([-11])

  mps_list[1]=SwapX*(Peps_bra[0]*Tem11)
  mps_list[1].permute([-10,2,3,11,-3,6],4)
  #mps_list.append(results)

  ####################Middle##############
  for  q  in  xrange(1,len(Peps_ket)):

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])
   Peps_ket[q].setLabel([6,-1,2,3,4])
   Tem1.setLabel([1])
   
   mps_list[2*q]=((Peps_ket[q]))*(SwapX*Tem1)
   mps_list[2*q].permute([10,6, 4, -10,2,3], 3)
   ##########################################

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])
   Peps_bra[q].setLabel([-10,-11,2,-3,-4])
   Tem11.setLabel([-11])

   mps_list[2*q+1]=SwapX*(Peps_bra[q]*Tem11)
   mps_list[2*q+1].permute([-10,2,3,11,-3,6],4)


  mps_boundry=MPSclass.MPS(2,2,len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()
  return mps_boundry

  ########################  Add next Layer  ###########

 if Location!=0:
  mps_list=[None]*(N_x*2)
  ###############################First-Layer###############################
  #print  mps_boundry[0].printDiagram()
  mps_boundry[0].setLabel([1,2,3])
  Tem1.setLabel([4])
  SwapX=fermionicOPT(Sys,mps_boundry[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([5,6,2,4])
  Tem1.setLabel([6])


  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([1,5,3,4],2)

  mps_boundry[1].setLabel([3,0,6])
  Peps_bra[0].setLabel([4,0,2,7,8])
  mps_list[1]=(mps_boundry[1]*Peps_bra[0])
  mps_list[1].permute([3,4,2,8,6,7],4)

  ###########################################
 

   ###############################################################################

  for q in xrange(1,len(Peps_ket)):
   ###########################################

    SwapX=fermionicOPT(Sys,mps_boundry[2*q].bond(1), Peps_bra[q].bond(0))
    SwapX.setLabel([5,6,2,4])
    mps_boundry[2*q].setLabel([1,2,3])

    mps_list[2*q]=mps_boundry[2*q]*SwapX
    mps_list[2*q].permute([1,6,5,3,4],3)



    mps_boundry[2*q+1].setLabel([3,0,6])
    Peps_bra[q].setLabel([4,0,2,7,8])

    mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_bra[q]
    mps_list[2*q+1].permute([3,4,2,8,6,7],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()


######################### Next Absorption ####################################

  mps_boundry[0].setLabel([6,5,4])
  #print "Peps_ket", Peps_ket[0].printDiagram()
  Peps_ket[0].setLabel([ 3, 5, 2, -4, -2 ])
  Tem1.setLabel([3])
  mps_list[0]=mps_boundry[0]*(Peps_ket[0]*Tem1)
  mps_list[0].permute([6,-2,4,2,-4],2)

  
  mps_boundry[1].setLabel([4,2,8,6])
  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(3), mps_boundry[1].bond(2) )
  SwapX.setLabel([ -4, -2, 3, 8 ])
  mps_list[1]=mps_boundry[1]*SwapX
  mps_list[1].permute([4,2,3,-2,6,-4],4)

  #print "Hi", mps_list[0].printDiagram(), mps_list[1].printDiagram()

  for  q  in  xrange(1, len(Peps_ket)):

   mps_boundry[2*q].setLabel([6,5,4])
   Peps_ket[q].setLabel([ 3, 5, 2, -4, -2 ])
   mps_list[2*q]=mps_boundry[2*q]*Peps_ket[q]
   mps_list[2*q].permute([6,3,-2,4,2,-4],3)

   mps_boundry[2*q+1].setLabel([4,2,8,6])
   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(3), mps_boundry[2*q+1].bond(2) )
   SwapX.setLabel( [ -4, -2, 3, 8 ] )
   mps_list[2*q+1]=mps_boundry[2*q+1]*SwapX
   mps_list[2*q+1].permute([4,2,3,-2,6,-4],4)




  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()

  return mps_boundry




#####@profile
def make_Env_singleLayer_up( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys, N_y):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):

  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])
  Peps_bra.append(A_conj*1.0)

 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

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

 Temr=uni10.UniTensor([bdi_1])
 Temr.identity()
 Temr.setLabel([10])

 Temrr=uni10.UniTensor([bdi_1])
 Temrr.identity()
 Temrr.setLabel([10])

 if Location==N_y-1:

  mps_list=[None]*(N_x*2)

  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(1), Peps_bra[0].bond(0) )
  SwapX.setLabel( [1,10,-1,-10] )
  Peps_ket[0].setLabel( [0,-1,2,3,4] )
  Tem1.setLabel( [10] )
  Temr.setLabel( [4] )

  mps_list[0]=(Peps_ket[0]*Temr)*(SwapX*Tem1)
  mps_list[0].permute( [0, 1, -10, 2 , 3] , 2 )

  ##########################################

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])
  Peps_bra[0].setLabel([-10,-11,2,-3,-4])
  Tem11.setLabel([11])

  mps_list[1]=(SwapX*Tem11)*Peps_bra[0]
  mps_list[1].permute([-10,2,3,-11,-3,6],4)
  #mps_list.append(results)

  ####################Middle##############
  for  q  in  xrange(1,len(Peps_ket)):

   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(1), Peps_bra[q].bond(0) )
   SwapX.setLabel([1,10,-1,-10])
   Peps_ket[q].setLabel([6,-1,2,3,4])
   Tem1.setLabel([4])
   
   mps_list[2*q]=((Peps_ket[q]))*(SwapX*Tem1)
   mps_list[2*q].permute([ 10, 6, 1, -10,2, 3], 3)
   ##########################################

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])
   Peps_bra[q].setLabel([-10,-11,2,-3,-4])
   Tem11.setLabel([11])

   mps_list[2*q+1]=(SwapX*Tem11)*(Peps_bra[q])
   mps_list[2*q+1].permute([-10,2,3,-11,-3,6],4)


  mps_boundry=MPSclass.MPS(2,2,len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()
  return mps_boundry

 ######################### Next Absorption ####################################

 if Location!=(N_y-1):

  mps_list=[None]*(N_x*2)
  mps_boundry[0].setLabel([6,5,4])
  Peps_ket[0].setLabel([ 3, -2, 2, -4, 5 ])
  Tem1.setLabel([3])
  mps_list[0]=mps_boundry[0]*(Peps_ket[0]*Tem1)
  mps_list[0].permute([6,-2,2,-4,4],3)

  mps_boundry[1].setLabel([4,-2,6])
  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(3), mps_boundry[1].bond(1) )
  SwapX.setLabel([ -4, 8, 3, -2 ])

  mps_list[1]=mps_boundry[1]*SwapX
  mps_list[1].permute([3,4,8,-4,6],3)



  for  q  in  xrange(1, len(Peps_ket)):

   mps_boundry[2*q].setLabel([6,5,4])
   Peps_ket[q].setLabel( [ 3, -2, 2, -4, 5 ] )
   mps_list[2*q]=mps_boundry[2*q]*Peps_ket[q]
   mps_list[2*q].permute( [3,6,-2,2,-4,4], 4 )

   mps_boundry[2*q+1].setLabel( [4,-2,6] )
   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(3), mps_boundry[2*q+1].bond(1) )
   SwapX.setLabel( [ -4, 8, 3, -2 ] )
   mps_list[2*q+1]=mps_boundry[2*q+1]*SwapX
   mps_list[2*q+1].permute([3,4,8,-4,6],3)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

#   A_t=mps_boundry.norm()
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))

########################  Add next Layer  ###########

###############################First-Layer###############################
  mps_boundry[0].setLabel([3,-2,2,-4])
  SwapX=fermionicOPT(Sys, mps_boundry[0].bond(1), Peps_bra[0].bond(0) )
  SwapX.setLabel([20,6,-2,4])
  Tem1.setLabel([6])
  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([3,20,4,2,-4],2)

  mps_boundry[1].setLabel([3,8,-4])
  Peps_bra[0].setLabel([4,0,2,7,8])
  mps_list[1]=(mps_boundry[1]*Peps_bra[0])
  mps_list[1].permute([4,2,3,0,7,-4],4)
###########################################
###########################################



  for q in xrange(1,len(Peps_ket)):
   ###########################################
   mps_boundry[2*q].setLabel([3,-2,2,-4])
   SwapX=fermionicOPT(Sys, mps_boundry[2*q].bond(1), Peps_bra[q].bond(0) )
   SwapX.setLabel([20,6,-2,4])
   mps_list[2*q]=(mps_boundry[2*q]*SwapX)
   mps_list[2*q].permute([6,3,20,4,2,-4],3)


   mps_boundry[2*q+1].setLabel([3,8,-4])
   Peps_bra[q].setLabel([4,0,2,7,8])
   mps_list[2*q+1]=(mps_boundry[2*q+1]*Peps_bra[q])
   mps_list[2*q+1].permute([4,2,3,0,7,-4],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t

  return mps_boundry





























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

def  bond_list_dist(T):
 bond_list=list(T.bond())
 bond_Inlist=[]
 bond_OUTlist=[]
 for i in xrange(len(bond_list)):
  if bond_list[i].type()==1:
   bond_Inlist.append(bond_list[i])
  else: 
   bond_OUTlist.append(bond_list[i])
 return bond_Inlist, bond_OUTlist

def  cal_rowcol(T):
 blk_qnums = T.blockQnum()
 row_list=[]
 col_list=[]
 for qnum in blk_qnums:
  M_tem=T.getBlock(qnum)
  col_list.append(M_tem.col())
  row_list.append(M_tem.row())
 return  sum(row_list),  sum(col_list)


def make_ENV_updown(MPO_Ten, N_x,N_y, chi):

 chi_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_list.append(dim)

 chi_dim=sum(chi_list)
 #print "chi_dim", chi_dim
 bdi1 = uni10.Bond(uni10.BD_IN, 1)
 Ten1=uni10.UniTensor([bdi1])
 Ten2=uni10.UniTensor([bdi1])
 Ten1.identity()
 Ten2.identity()

###################################################################
 Env_up=[None]*N_y
 cont_list=[None]*N_x
 for j in xrange(N_x):
  MPO_Ten[j][N_y-1].setLabel([1,-1,2,-2,4,-4,5,-5])
  Ten1.setLabel([5])
  Ten2.setLabel([-5])

  A=(MPO_Ten[j][N_y-1]*Ten1)*Ten2
  A.permute([1,-1,2,-2,4,-4],4)
  cont_list[j]=copy.copy(A)

 mps_A=MPSclass.MPS( 2, 2, N_x )

 for i in xrange(N_x):
  mps_A[i]=cont_list[i]*1.0

 Env_up[N_y-1]=mps_A.appSVD(chi_dim)


 for q in xrange(N_y-1):
  for j in xrange(N_x):
   A=Env_up[N_y-1-q][j]*1.0
   A.setLabel([10,5,-5,20])
   B=MPO_Ten[j][N_y-2-q]*1.0
   B.setLabel([1,-1,2,-2,4,-4,5,-5])
   Result=A*B
   Result.permute([10,1,-1,2,-2,20,4,-4],5)
   cont_list[j]=Result*1.0


  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=cont_list[i]*1.0


  Env_up[N_y-2-q]=mps_A.appSVD(sum(chi_list))
#############################################################


###################################################################
 Env_down=[None]*N_y
 cont_list=[None]*N_x
 for j in xrange(N_x):
  MPO_Ten[j][0].setLabel([1,-1,2,-2,4,-4,5,-5])
  Ten1.setLabel([2])
  Ten2.setLabel([-2])

  A=(MPO_Ten[j][0]*Ten1)*Ten2
  A.permute([1,-1,5,-5,4,-4],4)
  cont_list[j]=A*1.0

 mps_A=MPSclass.MPS( 2, 2, N_x )

 for i in xrange(N_x):
  mps_A[i]=cont_list[i]*1.0

 Env_down[0]=mps_A.appSVD(chi_dim)


 for q in xrange(N_y-1):
  for j in xrange(N_x):
   A=Env_down[q][j]*1.0
   A.setLabel([10,2,-2,20])
   B=MPO_Ten[j][q+1]*1.0
   B.setLabel([1,-1,2,-2,4,-4,5,-5])
   Result=A*B
   Result.permute([10,1,-1,5,-5,20,4,-4],5)
   cont_list[j]=Result*1.0


  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=cont_list[i]*1.0


  Env_down[q+1]=mps_A.appSVD(sum(chi_list))


 if N_y==1:
  Env_down=[None]*3
  Env_up=[None]*3
  bdi=uni10.Bond(uni10.BD_IN,1)
  bdo=uni10.Bond(uni10.BD_OUT,1) 
  T=uni10.UniTensor([bdi,bdi,bdi,bdo])
  mps_I=MPSclass.MPS(1,1,N_x)
  for i_ind in xrange(mps_I.N):
   T.identity()
   mps_I[i_ind]=T*1.0
  for q in xrange(N_y+1):
   Env_down[q]=mps_I*1.0
   Env_up[q]=mps_I*1.0
  
 return Env_up, Env_down



def make_ENVMPS( MPO_Ten, N_x, N_y, chi):

 chi_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_list.append(dim)

 chi_dim=sum(chi_list)
 #print "chi_dim", chi_dim

 bdi1 = uni10.Bond(uni10.BD_IN, 1)

 #bdi1 = uni10.Bond(uni10.BD_IN, 1)

 Ten1=uni10.UniTensor([bdi1])
 Ten2=uni10.UniTensor([bdi1])
 Ten1.identity()
 Ten2.identity()

 Env_left=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
  MPO_Ten[0][j].setLabel([1,-1,2,-2,4,-4,5,-5])

  #print MPO_Ten[0][j].bond(0)
  #Ten1=uni10.UniTensor([MPO_Ten[0][j].bond(0)])
  Ten2=uni10.UniTensor([MPO_Ten[0][j].bond(1)])  
  Ten1.setLabel([1])
  Ten2.setLabel([-1])
  Ten1.identity()
  Ten2.identity()


  A=(MPO_Ten[0][j]*Ten1)*Ten2
  A.permute([2,-2,4,-4,5,-5],4)
  cont_list[j]=copy.copy(A)

 mps_A=MPSclass.MPS( 2, 2, N_y)

 for i in xrange(N_y):
  mps_A[i]=copy.copy(cont_list[i])


 Env_left[0]=mps_A.appSVD(chi_dim)

 #print mps_A[0].printDiagram(), Env_left[0][0].printDiagram()

# Env_left[0]=mps_A.appSVD(chi_dim)

 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=copy.copy(Env_left[q][j])
   A.setLabel([10,1,-1,20])
   B=MPO_Ten[q+1][j]*1.0
   B.setLabel([1,-1,2,-2,3,-3,4,-4])
   Result=A*B
   Result.permute([10,2,-2,3,-3,20,4,-4],5)
   cont_list[j]=Result*1.0

  mps_A=MPSclass.MPS( 2, 2, N_y)
  for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

  Env_left[q+1]=mps_A.appSVD(sum(chi_list))


##################################################

###############################

 Env_right=[None]*N_x
 cont_list=[None]*N_y
 for j in xrange(N_y):
   MPO_Ten[N_x-1][j].setLabel([1,-1,2,-2,3,-3,4,-4])
   Ten1.setLabel([3])
   Ten2.setLabel([-3])

   A=(MPO_Ten[N_x-1][j]*Ten1)*Ten2
   A.permute([2,-2,1,-1,4,-4],4)
   cont_list[j]=copy.copy(A)


 mps_A=MPSclass.MPS( 2, 2, N_y)

 for i in xrange(N_y):
  mps_A[i]=cont_list[i]*1.0

 #print mps_A.norm()
 Env_right[N_x-1]=mps_A.appSVD(sum(chi_list))
#  print N_x-1,  Env_right[N_x-1].norm()#, ( Env_right[N_x-1].norm()-mps_A.norm()) / mps_A.norm()


 for q in xrange(N_x-1):
  for j in xrange(N_y):
   A=copy.copy(Env_right[N_x-1-q][j])
   A.setLabel([10,3,-3,20])
   B=copy.copy(MPO_Ten[N_x-2-q][j])
   B.setLabel([1,-1,2,-2,3,-3,4,-4])
   Result=A*B
   Result.permute([10,2,-2,1,-1,20,4,-4],5)
   cont_list[j]=copy.copy(Result)


  mps_A=MPSclass.MPS( 2, 2, N_y)
  for i in xrange(N_y):
   mps_A[i]=copy.copy(cont_list[i])

  Env_right[N_x-2-q]=mps_A.appSVD(sum(chi_list))
#   print N_x-2-q, Env_right[N_x-2-q].norm()#, mps_A.norm(), ( Env_right[N_x-2-q].norm()-mps_A.norm() ) / mps_A.norm()


 return Env_left, Env_right

def make_boundry_MPO( PEPS_listten, N_x,N_y,Sys):

 MPO=[None]*N_x
 for i in xrange(N_x):
   MPO[i]=[None]*N_y

 for i in xrange(N_x):
  for j in xrange(N_y):
   A=copy.copy(PEPS_listten[i][j])
   A.setLabel([1,2,3,4,5])
   A.permute([1,2,3,4,5],3)
   A_t=A*1.0
   A_t.transpose()
   A_t.setLabel([-4,-5,-1,-2,3])
   A_t.permute([-1,-2,3,-4,-5],5)
   A.permute([1,2,3,4,5],5)

   #print A.printDiagram(), Swap1.printDiagram(), Swap2.printDiagram() , A_t.printDiagram(), Swap4.printDiagram(), Swap3.printDiagram()

   Swap2=fermionicOPT(Sys,A_t.bond(3), A_t.bond(4))
   Swap2.setLabel([-6,7,-4,-5])

   Swap3=fermionicOPT(Sys,A_t.bond(0), A_t.bond(1))
   Swap3.setLabel([8,9,-1,-2])

   Swap4=fermionicOPT(Sys,A_t.bond(0), A.bond(1))
   Swap4.setLabel([-8,10,8,2])

   Swap1=fermionicOPT(Sys,A.bond(3), A_t.bond(4))
   Swap1.setLabel([6,11,4,7])

   
   A=(A*Swap1)*(((A_t*Swap2)*Swap3)*Swap4)
   A.permute([1,-8,10,9,6,-6,5,11],4)
   #A.combineBond([1,-1])
   #A.combineBond([2,-2])
   #A.combineBond([4,-4])
   #A.combineBond([5,-5])
   #A.permute([1,2,4,5],2)
   MPO[i][j]=A*1.0

 return   MPO


def Init_PEPS( N_x,N_y, D, d, q, Model):
  bdi=uni10.Bond(uni10.BD_IN, D)
  bdo=uni10.Bond(uni10.BD_OUT, D)

  bdiphy=uni10.Bond(uni10.BD_IN, d)
  bdophy=uni10.Bond(uni10.BD_OUT, d)
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN)
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD)

  q_list=[q0_even]
  bdi1=uni10.Bond(uni10.BD_IN, q_list)
  bdo1=uni10.Bond(uni10.BD_OUT, q_list)


  q_list=[q0_odd]
  bdi11=uni10.Bond(uni10.BD_IN, q_list)
  bdo11=uni10.Bond(uni10.BD_OUT, q_list)


  A_list=[None]*N_y
# B_list=[None]*N_x

  if q == 0:
   if N_y==1:
    if Model[2]=="ODD":
     A=uni10.UniTensor([bdi11,bdi1,bdiphy,bdo,bdo1], "A_first")
    else:
     A=uni10.UniTensor([bdi1,bdi1,bdiphy,bdo,bdo1], "A_first")

    A.randomize()
    A.setLabel([0,1,2,3,4])
    A.permute([0,1,2,3,4],3)
    A_list[0]=A
   else:
    for i in xrange(N_y):
     if i == 0:
       if Model[2]=="ODD":
        A=uni10.UniTensor([bdi11,bdi1,bdiphy,bdo,bdo], "A_first")
       else:
        A=uni10.UniTensor([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")

       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     elif i ==(N_y-1):
       A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     else:
       A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
  elif q==N_x-1:
   if N_y==1:
    A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo1,bdo1], "A_first")
    A.randomize()
    A.setLabel([0,1,2,3,4])
    A.permute([0,1,2,3,4],3)
    A_list[0]=A
   else:
    for i in xrange(N_y):
     if i == 0:
       A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     elif i ==(N_y-1):
       A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     else:
       A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
  else:
   if N_y==1:
    A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo,bdo1], "A_first")
    A.randomize()
    A.setLabel([0,1,2,3,4])
    A.permute([0,1,2,3,4],3)
    A_list[0]=A
   else:
    for i in xrange(N_y):
     if i == 0:
       A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     elif i ==(N_y-1):
       A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A
     else:
       A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
       A.randomize()
       A.setLabel([0,1,2,3,4])
       A.permute([0,1,2,3,4],3)
       A_list[i]=A




  return  A_list





def inverse_ten(Landa2):
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
def   setTruncationMPS(theta, chi):
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
        svs, bidxs = sv_mergemps(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),bdo_mid])
    GB.assign([bdi_mid,  theta.bond(2)])
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

def   sv_mergemps(svs, bidxs, bidx, sv_mat, chi, len_qn):
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



    
    
    
  
  
 
 

def svd_parityrl(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])
    GB=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB

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




def svd_parity5(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bdi1B=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bdi1.combine(bdi1B)
    bdo1=uni10.Bond(uni10.BD_OUT,bdi1.Qlist())


    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).Qlist())
    bdiB=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())
    bdiC=uni10.Bond(uni10.BD_IN,theta.bond(2).Qlist())
    
    bdi.combine(bdiB)
    bdi.combine(bdiC)
    bdo=uni10.Bond(uni10.BD_OUT,bdi.Qlist()) 
    #print "hi" , bdi.dim(), bdi1.dim(), theta.printDiagram()
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

 bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).Qlist())
 bdi1A=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())
 bdi.combine(bdi1A)
 bdo=uni10.Bond(uni10.BD_OUT,bdi.Qlist())

 bdi1=uni10.Bond(uni10.BD_IN,theta.bond(2).Qlist())
 bdiA=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
 bdiB=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())

 bdi1.combine(bdiA)
 bdi1.combine(bdiB)
 bdo1=uni10.Bond(uni10.BD_OUT,bdi1.Qlist())
 #print bdi.dim(), bdi1.dim()
# bdi1=uni10.Bond(uni10.BD_IN,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())
# bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())

# bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
# bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())

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







def setTruncation2(theta, chi):
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
    GA.assign([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo_mid])
    GB.assign([bdi_mid, theta.bond(4), theta.bond(5), theta.bond(6), theta.bond(7)])
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


def setTruncation3(theta, chi):
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


def svd_parity3(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB

def svd_parity4(theta):
    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim()*theta.bond(8).dim()*theta.bond(9).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim()*theta.bond(8).dim()*theta.bond(9).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3),theta.bond(4), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(5),theta.bond(6), theta.bond(7),theta.bond(8),theta.bond(9)])

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
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3),theta.bond(4), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(5),theta.bond(6), theta.bond(7),theta.bond(8),theta.bond(9)])

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

def Sqrt(Landa):
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


def Sqrt_minor(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum)
   Landam=Landa_cp.getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     if Landam[i*D+j] > 1.0e-10:
      Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
     else:
      Landa_cpm[i*D+j]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 


####@profile
def  cal_particle_double(PEPS_listten, N_x, N_y, h_coupling, d, Model, chi_boundry, D,Sys):

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 H=N_total(h_coupling, d, Model)
 H0=N_total0(h_coupling, d, Model)
 HN=N_totalN(h_coupling, d, Model)


 MPO_Ten=make_boundry_MPO(PEPS_listten, N_x, N_y,Sys)
 Env_left, Env_right=make_ENVMPS(MPO_Ten, N_x,N_y, chi_boundry)


 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_double_col(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_double_col(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_double_col(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

 E_c=sum(E_coulmn_t)
 E_coulmn=[ E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))  ]

#################   Row   #################

 mps_II=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_II.N):
  T.identity()
  mps_II[i_ind]=T*1.0


 MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,N_y,Sys)
 Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x,N_y, chi_boundry)

 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_double_row( mps_II, Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_double_row( Env_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)



 E_r=sum(E_row_t)
 #print "E=", E_r, E_c, (E_c+E_r)/(N_x*N_y)
 E_row=[ E_row_t[i]*1.0   for i in xrange(len(E_row_t))  ]
 #E_1=(E_c+E_r)/(N_x*N_y)
 E_1=(E_c+E_r)
 print E_row_t, E_coulmn_t
 return E_1





def  cal_particle_double_H(PEPS_listten, N_x, N_y, h_coupling, d, Model, chi_boundry, D,Sys):

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 H=N_total_H(h_coupling, d, Model)
 H0=N_total0_H(h_coupling, d, Model)
 HN=N_totalN_H(h_coupling, d, Model)


 MPO_Ten=make_boundry_MPO(PEPS_listten, N_x, N_y,Sys)
 Env_left, Env_right=make_ENVMPS(MPO_Ten, N_x,N_y, chi_boundry)


 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_double_col(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_double_col(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_double_col(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

 E_c=sum(E_coulmn_t)
 E_coulmn=[ E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))  ]

#################   Row   #################

 mps_II=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_II.N):
  T.identity()
  mps_II[i_ind]=T*1.0


 MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,N_y,Sys)
 Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x,N_y, chi_boundry)

 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_double_row( mps_II, Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_double_row( Env_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)



 E_r=sum(E_row_t)
 #print "E=", E_r, E_c, (E_c+E_r)/(N_x*N_y)
 E_row=[ E_row_t[i]*1.0   for i in xrange(len(E_row_t))  ]
 #E_1=(E_c+E_r)/(N_x*N_y)
 E_1=(E_c+E_r)
 print E_row_t, E_coulmn_t
 return E_1





def  cal_energy_double(PEPS_listten, N_x, N_y, h_coupling, d, Model, chi_boundry, D,Sys):

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 H=Heisenberg(h_coupling, d, Model)
 H0=Heisenberg0(h_coupling, d, Model)
 HN=HeisenbergN(h_coupling, d, Model)


 MPO_Ten=make_boundry_MPO(PEPS_listten, N_x, N_y,Sys)
 Env_left, Env_right=make_ENVMPS(MPO_Ten, N_x,N_y, chi_boundry)


 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_double_col(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_double_col(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)
  else:
   energy_double_col(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_y, H0, H, HN, D, E_coulmn_t,Sys)

 E_c=sum(E_coulmn_t)
 E_coulmn=[ E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))  ]

#################   Row   #################

 mps_II=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_II.N):
  T.identity()
  mps_II[i_ind]=T*1.0


 MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,N_y,Sys)
 Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x,N_y, chi_boundry)

 E_row_t=[]
 for i_ind in xrange(N_y):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   energy_double_row( mps_II, Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  elif i_ind==N_y-1:
   energy_double_row( Env_down[i_ind-1], mps_II, peps_l, N_x, H0, H, HN, D, E_row_t,Sys)
  else:
   energy_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, H0, H, HN, D, E_row_t,Sys)



 E_r=sum(E_row_t)
 #print "E=", E_r, E_c, (E_c+E_r)/(N_x*N_y)
 E_row=[ E_row_t[i]*1.0   for i in xrange(len(E_row_t))  ]
 #E_0=E_1*1.0
 #E_1=(E_c+E_r)/(N_x*N_y)
 #print E_row_t, E_coulmn_t
 E_1=(E_c+E_r)

 return E_1















def  energy_double_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, H0, H, HN , D, E_coulmn,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[i]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[i]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   #E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=Swap1*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H
  #print "Hi2", i
  E_val=energy_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, H_orig,D,Sys)
  E_coulmn.append(E_val)



#######@profile
def  energy_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, H_orig, D,Sys ):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 #mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])
 #mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,32])
 #mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location+1].setLabel([32,23,-10,19])
 #mps_boundry_down[Location*2+3].setLabel([33,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[Location]
 A=A*mps_boundry_up[Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[Location+1]


 B=E_right#*mps_boundry_up[2*Location+3]
 B=(Swap2*qq_d)*B
 #B=B*mps_boundry_down[2*Location+3]

 #B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 #B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0], H_orig



 return h_h[0]/Norm_h[0]








#######@profile
def energy_double_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, H0, H, HN , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   iden.setLabel([5])
   Peps_list=Peps_list*iden
   iden.setLabel([11])
   Peps_list_conj=Peps_list_conj*iden

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)
   
   iden.setLabel([13])
   mps_boundry_leftA=mps_boundry_left[i]*iden
   iden.setLabel([15])
   mps_boundry_rightA=mps_boundry_right[i]*iden

   result=(PEP_com*mps_boundry_rightA)*mps_boundry_leftA

   result.permute([12,10,9,14],4)
   E_list_up[i]=result
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)


   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_up[i+1].setLabel([13,5,11,15])
   result=((((PEP_com*E_list_up[i+1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([12,10,9,14],4)
   E_list_up[i]=result



 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   iden.setLabel([10])
   PEP_com=PEP_com*iden
   iden.setLabel([9])
   PEP_com=PEP_com*iden

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   result=((PEP_com*mps_leftA)*mps_rightA)

   result.permute([13,5,11,15],4)
   E_list_down[i]=result

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*((Peps_list_conj*Swap3)*Swap2))*Swap4

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_down[i-1].setLabel([12,10,9,14])

   result=((((PEP_com*E_list_down[i-1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([13,5,11,15],4)
   E_list_down[i]=result

#  for i in xrange(len(PEPS_listten)-1):
#   E_list_down[i].setLabel([1,2,3,4])
#   E_list_up[i+1].setLabel([1,2,3,4])
#   A=E_list_down[i]*E_list_up[i+1]
#   print "double" , A[0]


 for i in xrange(N_x-1):

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  E_val=double_local_energy_col( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H_orig,D,Sys)
  E_coulmn.append(E_val)



#######@profile
def  double_local_energy_col(PEPS_listten,E_list_down,E_list_up,mps_boundry_left,mps_boundry_right, Location, H_orig,D,Sys):

 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)


 r_u.setLabel([-100,3,10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])

######################################################################
 mps_boundry_left[Location].setLabel([16,6,-60,18])
 mps_boundry_left[Location+1].setLabel([18,11,-110,25])


 mps_boundry_right[Location].setLabel([17,90,-9,19])
 mps_boundry_right[Location+1].setLabel([19,140,-14,24])

######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*q)
 A=A*(Swap3*q_d)
 A=A*mps_boundry_right[Location]

 B=E_right*mps_boundry_left[Location+1]
 B=B*(Swap4*qq)
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_right[Location+1]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################


 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "Fnit_norm", Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]


######@profile
def  TEBD_Full_Single( N_x, N_y,PEPS_listten, E_iter_list, D, accuracy, N_tebd, i, d, chi_single, chi_try, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys):


 PEPS_listtenU=[None]*N_x
 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_up=[None]*(N_y+1)
 mps_boundry_down=[None]*N_y



 for i in xrange(N_x):
  PEPS_listtenU[i]=Init_PEPS( N_x,N_y, D, d, i,Model)


 Fidel_val=1
 E_coulmn=[]
 E_row=[]
 E_mag_coulmn=[]
 E_mag_row=[]
 E_0=1.0
 E_1=1.0

 Error_try=[]
 E_tebd=[]

 Energy_list=[]
 Sz_list=[]
 Nval_list=[]



 E_00=1.0
 E_11=1.0
 Energy_val=2.0

 mps_I=MPSclass.MPS(1,1,2*N_y)
 for i_ind in xrange(mps_I.N):
  mps_I[i_ind].identity()

 mps_II=MPSclass.MPS(1,1,2*N_x)
 for i_ind in xrange(mps_II.N):
  mps_II[i_ind].identity()


 List_delN=Short_TrotterSteps(N_tebd)

 print  h_coupling,  List_delN
 E_00_f=10.0
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 for delta, N_iter in List_delN:
  #delta=0
  print delta, N_iter

  H=Heisenberg(h_coupling, d, Model)
  H0=Heisenberg0(h_coupling, d, Model)
  HN=HeisenbergN(h_coupling, d, Model)

  U_ham = uni10.UniTensor( H.bond(), "U")
  blk_qnums = H.blockQnum()
  for qnum in blk_qnums:
   U_ham.putBlock(qnum, uni10.takeExp(-delta, H.getBlock(qnum)))

  U_ham0 = uni10.UniTensor( H0.bond(), "U");
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
   U_ham0.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))

  U_hamN = uni10.UniTensor( HN.bond(), "U");
  blk_qnums = HN.blockQnum()
  for qnum in blk_qnums:
   U_hamN.putBlock(qnum, uni10.takeExp(-delta, HN.getBlock(qnum)))

  PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x,N_y, D, chi_try, d, threshold, interval,Sys)
  Energy_val=Energy_cal(PEPS_listten, d, chi_single, N_x,N_y, D, Model, h_coupling,Sys)
  E_0=Energy_val
  #E_00=E_00_f*1.0
  E_11=Energy_val*1.0
  print "E_0", Energy_val

  PEPS_listtenU=copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU)
  PEPS_listten=copy_f( PEPS_listtenU, N_x,N_y, PEPS_listten)

  for q_iter in xrange(N_iter[0]):
   PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, N_y, D, chi_try, d, threshold, interval,Sys)
   Energy_val=Energy_cal(PEPS_listten, d, chi_single, N_x,N_y, D, Model, h_coupling,Sys)
   Sz_val   =     Sz_cal(PEPS_listten, d, chi_single, N_x, N_y, D, Model, h_coupling,Sys)
   N_val   =     Particle_cal(PEPS_listten, d, chi_single, N_x, N_y, D, Model, h_coupling,Sys)

   print "E_0", Energy_val
   print "Sz_val", Sz_val
   print "N_val", N_val

   Energy_list.append(Energy_val)
   Sz_list.append(Sz_val)
   Nval_list.append(N_val)


   file = open("E_peps0.txt", "w")
   for index in range(len(Energy_list)):
    file.write(str(index) + " " + str(Energy_list[index])+" " + str(Sz_list[index])+" " + str(Nval_list[index])+" "+ "\n")
   file.close()



##########################   Col   #############################
   for Location in reversed(xrange(1,N_x)):
    mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_y,Sys,N_x)
#############
   E_coulmn_t=[]
   for i_ind in xrange(N_x):
    if i_ind==0:
     update_twotensor(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_twotensor(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y,U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)
    else:
     update_twotensor(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_y,U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)

    if  i_ind==(N_x-1):break
    mps_boundry_left[i_ind]=make_Env_singleLayer(PEPS_listten[i_ind], i_ind, mps_boundry_left[i_ind-1], d, chi_single, N_y,Sys)


###########################   Row   ###############################

   for Location in reversed(xrange(1,N_y)):
    peps_l=[]
    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][Location])
    mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys,N_y)

   E_row_t=[]

   for i_ind in xrange(N_y):
    peps_l=[]
    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][i_ind])

    if i_ind==0:
     update_twotensor_row(mps_II, mps_boundry_up[i_ind+1], peps_l, N_x, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t,threshold, interval,Sys)
    elif i_ind==N_y-1:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_II, peps_l, N_x,U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t, threshold, interval,Sys)
    else:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t,threshold, interval,Sys)

    for i in xrange(N_x):
     PEPS_listten[i][i_ind]=peps_l[i]

    if  i_ind==(N_y-1):break
    mps_boundry_down[i_ind]=make_Env_singleLayer_down( peps_l, i_ind, mps_boundry_down[i_ind-1], d, chi_single, N_x,Sys)




   E_00=E_11*1.0
   E_11=(sum(E_coulmn_t)+sum(E_row_t))
   #E_11=sum(E_coulmn_t)
   #E_11=sum(E_row_t)
   #print  sum(E_coulmn_t),  sum(E_row_t)
   print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t),  sum(E_row_t), (E_11-E_00) / E_11 
   E_iter_list1.append(E_11)


   Error_try.append((E_11-E_00) / E_11)
   E_tebd.append(E_11)

   file = open("TEBD.txt", "w")
   for index in range(len(Error_try)):
    file.write(str(index) + " " + str(Error_try[index])+" " + str(E_tebd[index])+" "+ "\n")
   file.close()

   if N_iter[1]=="off": 
    Store_f(PEPS_listten, N_x, N_y)


   break_criteria="off"
   if N_iter[1]=="on":
    if E_11 < E_00 and q_iter>0:
     Store_f(PEPS_listten, N_x, N_y)
     PEPS_listtenU=copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU)


    if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
     if abs((E_11-E_00)/E_11)<accuracy:
      print "loop_finished:reached_levelofaccuracy"
      E_iter_list1.append("loop_finished:reached_levelofaccuracy")
     else: 
      print "Didnot_get_lowered";
      E_iter_list1.append("Didnot_get_lowered")
     PEPS_listten=copy_f( PEPS_listtenU, N_x,N_y, PEPS_listten)
     break_criteria="on"
     break

   if break_criteria=="on": break


# 
# 
# 
#    E_00=E_11*1.0
#    E_11=(sum(E_coulmn_t)+sum(E_row_t))
#    #E_11=sum(E_coulmn_t)
#    #E_11=sum(E_row_t)
#    #print  sum(E_coulmn_t),  sum(E_row_t)
#    print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t),  sum(E_row_t), (E_11-E_00) / E_11 
#    E_iter_list1.append(E_11)
# 
#    if E_11 < E_00 and q_iter>0:
#     Store_f(PEPS_listten, N_x, N_y)
#     PEPS_listtenU=copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU)
#     #E_00_f=E_11*1.0
# 
# 
#    Error_try.append((E_11-E_00) / E_11)
#    E_tebd.append(E_11)
# 
#    file = open("TEBD.txt", "w")
#    for index in range(len(Error_try)):
#     file.write(str(index) + " " + str(Error_try[index])+" " + str(E_tebd[index])+" "+ "\n")
#    file.close()
# 
# 
# 
#    if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
#     if abs((E_11-E_00)/E_11)<accuracy:print "loop_finished:reached_levelofaccuracy"
#     else: print "Didnot_get_lowered";
#     PEPS_listten=copy_f( PEPS_listtenU, N_x,N_y, PEPS_listten)
#     
#     break


####@profile
def  TEBD_Full_double( N_x,N_y, PEPS_listten, E_iter_list, D, accuracy, N_tebd, i, d, chi_boundry, chi_try, mps_boundry_left, mps_boundry_right, threshold, interval, Model,h_coupling,count_list,E_iter_list1,Sys):

 PEPS_mps_leftU=[None]*N_x
 PEPS_mps_rightU=[None]*N_x
 PEPS_listtenU=[None]*N_x
 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*N_x
 mps_boundry_temp=[None]*N_x

 for i in xrange(N_x):
  PEPS_listtenU[i]=Init_PEPS( N_x,N_y, D, d, i,Model)

 Fidel_val=1
 E_coulmn=[]
 E_row=[]
 E_mag_coulmn=[]
 E_mag_row=[]
 E_0=1.0
 E_1=1.0

 E_00=1.0
 E_11=1.0


 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_y)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 mps_II=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_II.N):
  T.identity()
  mps_II[i_ind]=T*1.0



 List_delN=Short_TrotterSteps(N_tebd)

 print h_coupling, List_delN
 E_00_f=10.0
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 for delta, N_iter in List_delN:
  #delta=0
  print delta, N_iter

  H=Heisenberg(h_coupling, d, Model)
  H0=Heisenberg0(h_coupling, d, Model)
  HN=HeisenbergN(h_coupling, d, Model)

  U_ham = uni10.UniTensor( H.bond(), "U")
  blk_qnums = H.blockQnum()
  for qnum in blk_qnums:
   U_ham.putBlock(qnum, uni10.takeExp(-delta, H.getBlock(qnum)))

  U_ham0 = uni10.UniTensor( H0.bond(), "U");
  blk_qnums = H0.blockQnum()
  for qnum in blk_qnums:
   U_ham0.putBlock(qnum, uni10.takeExp(-delta, H0.getBlock(qnum)))

  U_hamN = uni10.UniTensor( HN.bond(), "U");
  blk_qnums = HN.blockQnum()
  for qnum in blk_qnums:
   U_hamN.putBlock(qnum, uni10.takeExp(-delta, HN.getBlock(qnum)))


  PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x,N_y, D, chi_try, d, threshold, interval,Sys)
  Energy_val=cal_energy_double(PEPS_listten, N_x,N_y, h_coupling, d, Model, chi_boundry, D,Sys)
  E_0=Energy_val
  E_11=1.0*Energy_val
  print "E_0", Energy_val

  PEPS_listtenU=copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU)
  PEPS_listten=copy_f( PEPS_listtenU, N_x,N_y, PEPS_listten)

  for q_iter in xrange(N_iter):
   #Energy_val=cal_energy_double(PEPS_listten, N_x,N_y, h_coupling, d, Model, chi_boundry, D,Sys)
   #print q_iter, "E_0", Energy_val

###########################   Col   #############################
###   make_right_Env   ######
   MPO_Ten=make_boundry_MPO( PEPS_listten, N_x, N_y,Sys)
   Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, N_y, chi_boundry)
#############
   E_coulmn_t=[]
   for i_ind in xrange(N_x):
    if i_ind==0:
     update_local_double(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_y, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_local_double(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_y, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)
    else:
     update_local_double(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_y, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_coulmn_t,threshold, interval,Sys)


 
    
    if  i_ind==(N_x-1):break
    MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,N_y,Sys)
    Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, N_y,chi_boundry)


#    mps_boundry_left[i_ind]=make_Env_singleLayer(PEPS_listten[i_ind], i_ind, mps_boundry_left[i_ind-1], d, chi_single, N_x)
#    norm_val=mps_boundry_left[i_ind].product(mps_boundry_right[i_ind+1])
#    if norm_val>threshold or norm_val<(1.0/threshold):
#     #print "Warning", norm_val
#     PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval)
#     #print "Fixed", norm_val, "Num_iter", count_val
#     mps_boundry_right, mps_boundry_left=update_env_LR(PEPS_listten, N_x, d, chi_single, mps_boundry_left, mps_boundry_right, mps_boundry_temp)


##########################   Row   ###############################
   #PEPS_listten=rotate_all(PEPS_listten, N_x)


   MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,N_y,Sys)
   Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x,N_y, chi_boundry)

   #Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, chi_boundry)


   E_row_t=[]
   for i_ind in xrange(N_y):

    peps_l=[]

    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][i_ind])

    if i_ind==0:
     update_local_double_row( mps_II, Env_up[i_ind+1], peps_l, N_x, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t, threshold, interval,Sys)
    elif i_ind==N_y-1:
     update_local_double_row( Env_down[i_ind-1], mps_II, peps_l, N_x, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t, threshold, interval,Sys)
    else:
     update_local_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, U_ham0, U_ham, U_hamN, H0, H, HN, D, E_row_t, threshold, interval,Sys)

    for i in xrange(N_x):
     PEPS_listten[i][i_ind]=peps_l[i]

    if i_ind==(N_y-1):break
    MPO_Ten=make_boundry_MPO( PEPS_listten,  N_x,N_y,Sys)
    Env_up, Env_down=make_ENV_updown( MPO_Ten,  N_x,N_y, chi_boundry)



   E_00=E_11*1.0
   E_11=(sum(E_coulmn_t)+sum(E_row_t))#/(N_x*N_y)
   print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t), sum(E_row_t),(E_11-E_00) / E_11 , sum(E_row_t)
   E_iter_list1.append(E_11)

   if E_11 < E_00 and q_iter>0:
    Store_f(PEPS_listten, N_x, N_y)
    PEPS_listtenU=copy_f(PEPS_listten,N_x,N_y, PEPS_listtenU)
    #E_00_f=E_11*1.0


   if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
    if abs((E_11-E_00)/E_11)<accuracy:print "loop_finished:reached_levelofaccuracy"
    else: print "Didnot_get_lowered";
    PEPS_listten=copy_f( PEPS_listtenU, N_x,N_y, PEPS_listten)
    break


 return PEPS_listten





####@profile
def update_local_double( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, U_ham0, U_ham, U_hamN,H0, H, HN , D, E_coulmn,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[i]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[i]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])

   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[i]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[i]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)


#  for i in xrange(len(PEPS_listten)-1):
#   E_list_down[i].setLabel([1,2,3,4])
#   E_list_up[i+1].setLabel([1,2,3,4])
#   A=E_list_down[i]*E_list_up[i+1]
#   print "Inside", A[0]



 for i in xrange(N_x-1):

  Ham_u=1
  if i==0:
    Ham_u= U_ham0
  elif i==(N_x-2):
    Ham_u= U_hamN
  else:
   Ham_u= U_ham

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=double_local_update( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H_orig, Ham_u, D, threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=double_local_update_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H_orig, Ham_u, D, threshold, interval,Sys)


  E_coulmn.append(E_val)
  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[i]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)


###@profile
def norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*iden_h

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))

 val1=A*B

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)
 A=A*Ham

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))

 val2=A*B

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*Ham

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0]






###@profile
def norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*iden_h

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val1=A*B

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val2=A*B

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0]




###@profile
def optimum_0_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 Env_r=mps_boundry_right[Location+1]*E_right
 Env_r=(Env_r*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 Env_r=((Env_r*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))
 Env_r=Env_r*iden_h

 Env_r=(Env_r*mps_boundry_right[Location])
 Env_r=(Env_r*Swap4)*(Swap3)
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap6

 A=(mps_boundry_left[Location]*Swap1)
 A=A*(E_left*Swap2)
 Env_r=Env_r*A
 Env_r.permute([ 53,56,39, 54,57,6,35,37,55,34],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=mps_boundry_right[Location+1]*E_right
 Env_s=(Env_s*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
 Env_s=Env_s*Ham
 Env_s=((Env_s))*(Swap6)

 A=E_left*(mps_boundry_left[Location]*Swap1)
 A=A*((Peps_1d*Swap2)*Swap3)
 A=A*mps_boundry_right[Location]
 A=(A*Swap4)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)



 Env_s1=mps_boundry_right[Location+1]*E_right
 Env_s1=(Env_s1*Swap12)*((Peps_2*Swap11)*mps_boundry_left[Location+1])
 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))
 Env_s1=Env_s1*Ham
 Env_s1=((Env_s1)*(Swap4*Swap5))

 A=E_left*(mps_boundry_left[Location]*Swap1)
 A=A*(Peps_1*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[Location]

 Env_s1=Env_s1*((Swap2)*Swap3)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,35,37,55,34,53,56,39,54,57])

 
 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,35,37,55,34],5)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([53,56,39,54,57])

 return Peps_1p, Peps_1pd





###@profile
def optimum_0_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):


 Env_r=E_right*mps_boundry_up[Location+1]
 Env_r=Swap7*Env_r
 Env_r=Env_r*mps_boundry_down[Location+1]

 Env_r=Env_r*((Peps_2pd*Swap6)*Swap5)
 Env_r=(Env_r*Swap4)*Swap12
 Env_r=Env_r*(Peps_2p*Swap11)
 Env_r=Env_r*iden_h

 Env_r=Env_r*(Swap10*mps_boundry_up[Location])

 Env_r=(Env_r)*(Swap9)
 Env_r=(Env_r)*(Swap8)
 Env_r=Env_r*Swap3


 A=E_left*mps_boundry_down[Location]
 A=A*Swap1
 A=A*Swap2


 Env_r=(Env_r)*(A)


 Env_r.permute([ 32,33,37, 35,34,6,31,58,59,10],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=E_right*mps_boundry_up[Location+1]
 Env_s=Swap7*Env_s
 Env_s=Env_s*mps_boundry_down[Location+1]

 Env_s=Env_s*((Peps_2d*Swap6)*Swap5)
 Env_s=(Env_s*Swap4)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*Ham


 Env_s=(Env_s*(Swap8*Swap9))*Swap10


 Env_s=(Env_s*((Peps_1d*Swap3)*Swap2))*mps_boundry_down[Location]
 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[Location]

 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*mps_boundry_up[Location+1]
 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*mps_boundry_down[Location+1]

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5)
 Env_s1=(Env_s1*Swap4)*Swap12
 Env_s1=Env_s1*(Peps_2*Swap11)
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*(mps_boundry_up[Location]*Swap10)
 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3))
 Env_s1=Env_s1*Peps_1
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*((mps_boundry_down[Location]*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,31,58,59,10,32,33,37,35,34])

 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,31,58,59,10],5)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([32,33,37,35,34])

 return Peps_1p, Peps_1pd


###@profile
def optimum_1_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):

 Env_r=E_left*mps_boundry_up[Location]
 Env_r=Env_r*(Peps_1p*Swap1)
 Env_r=Env_r*mps_boundry_down[Location]
 Env_r=Env_r*((Swap3*Peps_1pd)*Swap2)
 Env_r=(Env_r*(Swap8*Swap9))*Swap10
 Env_r=Env_r*iden_h

 Env_r=Env_r*Swap4
 Env_r=Env_r*mps_boundry_down[Location+1]
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap12
 Env_r=Env_r*Swap11

 A=mps_boundry_up[Location+1]*E_right
 A=Swap7*A
 A=Swap6*A

 Env_r=(Env_r*A)

 Env_r.permute([43,42,47,45,44,54,52,51,53,15],5)
 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*mps_boundry_up[Location]
 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[Location]
 Env_s=Env_s*((Swap3*Peps_1d)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10
 Env_s=Env_s*Ham
 Env_s=((Swap11*Env_s)*Swap12)*Swap4
 Env_s=Env_s*mps_boundry_down[Location+1]
 Env_s=Env_s*((Swap5*Peps_2d)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*mps_boundry_up[Location+1]

 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*mps_boundry_up[Location]
 Env_s1=Env_s1*(Peps_1*Swap1)
 Env_s1=Env_s1*mps_boundry_down[Location]
 Env_s1=Env_s1*((Swap3*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((Peps_2*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4

 Env_s1=Env_s1*mps_boundry_up[Location+1]

 Env_s1=Env_s1*mps_boundry_down[Location+1]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([54,52,51,53,15,43,42,47,45,44])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([54,52,51,53,15],5)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([43,42,47,45,44])

 return Peps_2p, Peps_2pd



###@profile
def optimum_1_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 Env_r=E_left*mps_boundry_left[Location]
 Env_r=Env_r*(Swap1*((Peps_1pd*Swap3)*Swap2))
 Env_r=Env_r*mps_boundry_right[Location]
 Env_r=((Env_r*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_r=Env_r*iden_h
 Env_r=((Env_r*Swap10)*Swap11)
 Env_r=((Env_r*Swap9)*Swap8)
 Env_r=Env_r*mps_boundry_left[Location+1]
 A=mps_boundry_right[Location+1]*Swap12
 A=A*Swap7
 A=A*E_right
 Env_r=Env_r*A

 Env_r.permute([ 48,58,46, 49,57,11,31,32,33,15],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*mps_boundry_left[Location]
 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3)*Swap2))
 Env_s=Env_s*mps_boundry_right[Location]
 Env_s=((Env_s*Swap4)*Swap5)*(Peps_1*Swap6)
 Env_s=Env_s*Ham
 Env_s=((Env_s*Swap10)*Swap11)
 Env_s=(Env_s*mps_boundry_left[Location+1])*Peps_2
 Env_s=Env_s*E_right
 Env_s=Env_s*(mps_boundry_right[Location+1]*Swap12)
 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8)
 Env_s.permute([48,58,46,49,57],5)



 Env_s1=E_left*mps_boundry_left[Location]
 Env_s1=Env_s1*(Swap1*((Peps_1d*Swap3)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[Location]
 Env_s1=((Env_s1*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*Ham
 A=((mps_boundry_right[Location+1]*Swap12)*(Swap7*(Swap8*Peps_2d)))*E_right
 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11)
 Env_s1=(Env_s1*mps_boundry_left[Location+1])
 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([11,31,32,33,15,48,58,46, 49,57])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([11,31,32,33,15],5)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([48,58,46, 49,57])

 return Peps_2p, Peps_2pd






###@profile
def double_local_update_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, Ham ,D, threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 Peps_1.setLabel([6,35,37,55,34])
 Peps_1.permute([6,35,37,55,34],5)
 A.setLabel([6,35,37,55,34])
 A.permute([6,35,37,55,34],5)

 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(0))
 Swap1.setLabel([70,-60,35,36])
 Swap2=Swap1*1.0
 Swap2.setLabel([-7,36,56,53])
 Swap3=fermionicOPT(Sys,A.bond(4), A.bond(3))
 Swap3.setLabel([57,54,40,-9])
 Swap4=Swap3*1.0
 Swap4.setLabel([40,58,41,90])
 Swap6=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap6.setLabel([37,55,42,38])
 Swap5=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap5.setLabel([38,39,58,44])


#####################################################

 A1=Peps_2*1.0
 A1.setLabel([11,31,32,33,15])
 A1.permute([11,31,32,33,15],5)
 Peps_2.setLabel([11,31,32,33,15])
 Peps_2.permute([11,31,32,33,15],5)

 Swap7=fermionicOPT(Sys,A1.bond(3), A1.bond(4))
 Swap7.setLabel([49,57,-14,50])
 Swap8=fermionicOPT(Sys,A1.bond(0), A1.bond(1))
 Swap8.setLabel([48,58,47,41])
 Swap9=fermionicOPT(Sys,A1.bond(0), A1.bond(2))
 Swap9.setLabel([51,45,47,46])
 Swap10=fermionicOPT(Sys,A1.bond(0), A1.bond(2))
 Swap10.setLabel([52,32,51,43])
 Swap11=fermionicOPT(Sys,A1.bond(0), A1.bond(1))
 Swap11.setLabel([-110,34,52,31])
 Swap12=fermionicOPT(Sys,A1.bond(3), A1.bond(4))
 Swap12.setLabel([33,50,140,-150])



######################################################################
 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([53,56,39,54,57])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([48,58,46,49,57])

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1c", Norm_h[0], h_h[0]/Norm_h[0]


 #N_ten=N_Positiv(N_ten)

######################################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 Peps_1p=Peps_1*1.0
 Peps_2p=Peps_2*1.0
 Peps_1pd=Peps_1d*1.0
 Peps_2pd=Peps_2d*1.0

 Peps_2d.setLabel([48,58,46,49,57])
 Peps_1d.setLabel([53,56,39,54,57])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([44,45,42,43])
 Ham.setLabel([44,45,42,43])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf=norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)

  Peps_2p, Peps_2pd=optimum_1_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count







 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([53,56,39,54,57])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([48,58,46,49,57])

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))

 N_ten=A*B

 H_orig.setLabel([44,45,42,43])
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2c", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold or Norm_h[0]<(1.0/threshold):
  Peps_1p, Peps_2p=renormalize_tensor_c(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]



###@profile
def renormalize_tensor_c(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval):

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold and abs(norm_val)>(1.0/threshold):break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_left[Location]
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location]
   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
   B=mps_boundry_right[Location+1]*E_right
   B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(1.0/threshold):
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_left[Location]
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location]
   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
   B=mps_boundry_right[Location+1]*E_right
   B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p


###@profile
def  Update_twotensor_double_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D, threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])



 Peps_1.setLabel([6,31,58,59,10])
 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([32,33,37,35,34])
 Peps_2.setLabel([54,52,51,53,15])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([43,42,47,45,44])


 Swap1=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([31,17,22,30])
 Swap2=Swap1*1.0
 Swap2.setLabel([33,32,-7,30])
 Swap3=fermionicOPT(Sys,Peps_1.bond(4), Peps_1.bond(3))
 Swap3.setLabel([34,35,38,36])
 Swap4=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap4.setLabel([48,36,23,41])
 Swap5=Swap4*1.0
 Swap5.setLabel([42,43,-10,41])
 Swap6=fermionicOPT(Sys,Peps_2.bond(4), Peps_2.bond(3))
 Swap6.setLabel([44,45,46,-14])

 Swap7=fermionicOPT(Sys,Peps_2.bond(4), Peps_2.bond(3))
 Swap7.setLabel([46,53,24,20])

 Swap8=fermionicOPT(Sys,Peps_1.bond(2), Peps_1.bond(4))
 Swap8.setLabel([37,38,40,39])

 Swap9=Swap8*1.0
 Swap9.setLabel([58,39,56,57])

 Swap10=Swap3*1.0
 Swap10.setLabel([57,59,25,54])

 Swap11=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap11.setLabel([55,50,51,52])

 Swap12=Swap11*1.0
 Swap12.setLabel([47,48,49,50])



######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,322])
 mps_boundry_down[Location+1].setLabel([322,23,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 
 N_ten.permute([40,49,56,55],2)



 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1r", Norm_h[0], h_h[0]/Norm_h[0]



 Peps_1p=Peps_1*1.0
 Peps_2p=Peps_2*1.0
 
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([32,33,37,35,34])
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([43,42,47,45,44])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([40,49,56,55])

 Ham.setLabel([40,49,56,55])
 
 iden_h=Ham*1.0
 iden_h.setLabel([40,49,56,55])
 iden_h.identity()

 iden_h.setLabel([40,49,56,55])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf=norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)

  Peps_2p, Peps_2pd=optimum_1_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)








 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]
 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)
 N_ten=A*B
 #N_ten.permute([40,49,56,55],2)


 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2r", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold or Norm_h[0]<(1.0/threshold):
  Peps_1p, Peps_2p=renormalize_tensor_r(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]




###@profile
def renormalize_tensor_r(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval):

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]
 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold and abs(norm_val)>(1.0/threshold):break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_up[Location]
   A=A*(Peps_1p*Swap1)
   A=A*mps_boundry_down[Location]
   A=A*((Swap3*Peps_1pd)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10
   B=E_right*mps_boundry_up[Location+1]
   B=Swap7*B
   B=B*mps_boundry_down[Location+1]
   B=B*((Peps_2pd*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2p*Swap11)
   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(1.0/threshold):
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_up[Location]
   A=A*(Peps_1p*Swap1)
   A=A*mps_boundry_down[Location]
   A=A*((Swap3*Peps_1pd)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10
   B=E_right*mps_boundry_up[Location+1]
   B=Swap7*B
   B=B*mps_boundry_down[Location+1]
   B=B*((Peps_2pd*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2p*Swap11)
   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

  A=E_left*mps_boundry_up[Location]
  A=A*(Peps_1p*Swap1)
  A=A*mps_boundry_down[Location]
  A=A*((Swap3*Peps_1pd)*Swap2)
  A=(A*(Swap8*Swap9))*Swap10
  B=E_right*mps_boundry_up[Location+1]
  B=Swap7*B
  B=B*mps_boundry_down[Location+1]
  B=B*((Peps_2pd*Swap6)*Swap5)
  B=(B*Swap4)*Swap12
  B=B*(Peps_2p*Swap11)
  N_ten=A*B
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p






####@profile
def double_local_update( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, Ham ,D, threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],3)


 r_u.setLabel([-100,3,10])

 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,9,200])
 q_d.permute([6,7,9,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-9,-200,9,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-9,-200],4)

#####################################################

 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])


 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,11,14,15])
 qq_d.permute([400,11,14,15],4)
 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel( [-400, -11, 400, 11] )
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-11,-14,-15],4)



######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([90,200,9,-200])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*q_d)
 A=A*mps_boundry_right[Location]
 A=A*(Swap3*q)


 B=mps_boundry_left[Location+1]*E_right
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[Location+1]
 B=B*(Swap4*qq_d)


 N_ten=A*B
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0]


 N_ten=N_Positiv(N_ten)

######################################################

###############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 D_dim=l_u.bond(0).dim()
 #chi_dim=sum(D)
 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)
 
   
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
###########################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])
 Ham.setLabel([-3,-13,3,13])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count


 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
 #print "E_2", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold or Norm_h[0]<(1.0/threshold):
  r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)

 A=PEPS_1*1.0
 A.permute([6,7,3,9,10],5)

 Swap1=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap1.setLabel([-3,-9,3,9])
 A=Swap1*A
 A.permute([6,7,-3,-9,10],3)
 PEPS_1=A*1.0


 PEPS_2=l_up*qq
 PEPS_2.permute([11,10,13,14,15],3)


 

 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]







####@profile
def update_local_double_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, U_ham0, U_ham, U_hamN, H0, H, HN , D, E_coulmn,threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[i]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[i]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   #E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=Swap1*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):
  #print "i", i
  Ham_u=1
  if i==0:
    Ham_u= U_ham0
  elif i==(N_x-2):
    Ham_u= U_hamN
  else:
   Ham_u= U_ham

  H_orig=1
  if i==0:
   H_orig= H0
  elif i==(N_x-2):
   H_orig= HN
  else:
   H_orig= H

  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_double_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, Ham_u, H_orig,D,threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_double_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, Ham_u, H_orig,D,threshold, interval,Sys)


  E_coulmn.append(E_val)

  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)


####@profile
def  Update_twotensor_double_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D, threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,32])
 mps_boundry_down[Location+1].setLabel([32,23,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(q*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*(Swap3*q_d)

 B=E_right*mps_boundry_up[Location+1]
 B=(Swap2*qq_d)*B

 B=B*mps_boundry_down[Location+1]
 B=B*(Swap4*qq)
 N_ten=A*B
 
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0]


 N_ten=N_Positiv(N_ten)
######################################################

#############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])
 
 
 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 
 D_dim=l_u.bond(0).dim()

 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)
 
 
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 
 

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
#########################################


 
# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])  #new
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])

 Ham.setLabel([-3,-13,3,13])
 
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12 or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
 #print "E_2", Norm_h[0], h_h[0]/Norm_h[0]
# #print Norm_h[0]


 if Norm_h[0]> threshold or Norm_h[0]<(1.0/threshold):
   r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 r_up.setLabel([-100,3,9])   #new
 l_up.setLabel([9,13,-300])  #new


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)


 PEPS_2=l_up*qq
 A=PEPS_2*1.0
 A.permute([9,10,13,14,15],5)
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([9,-10,-13,14,15],3)
 PEPS_2=A*1.0
 PEPS_2.permute([9,-10,-13,14,15],3)


 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]


















