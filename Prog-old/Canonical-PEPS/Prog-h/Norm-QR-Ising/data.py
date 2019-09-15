import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass 


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



    
    
    
  
  
 
 



def Init_A_PEPS(N,bdi,bdi1,bdip,bdo,bdo1,bdop):
 A_list=[None]*N
 for i in xrange(N):
  if i == 0:
    A=uni10.UniTensor([bdi1,bdip,bdi,bdo], "A_first")
    A.randomize()
    A.orthoRand()
    A.setLabel([0,1,2,3])
    A.combineBond([1,2])
    #A.permute([1,0,4],2)
    A_list[i]=A
    #print A_list[i].printDiagram()
  elif i ==(N-1):
    A=uni10.UniTensor([bdi,bdip,bdi,bdo1], "A_end")
    A.randomize()
    A.orthoRand()
    A.setLabel([0,1,2,3])
    A.combineBond([1,2])
    A_list[i]=A
    #print A_list[i].printDiagram()
  else:
    A=uni10.UniTensor([bdi,bdip,bdi,bdo], "A_middle")
    A.randomize()
    A.orthoRand()
    A.setLabel([0,1,2,3])
    A.combineBond([1,2])
    A_list[i]=A
    #print A_list[i].printDiagram()

 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N, 'ortho')
 for i in xrange(N):
   mps_A[i]=copy.copy(A_list[i])
 return mps_A


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

def Cost_Function(T0, U, U_transpose):

 #print   U_transpose*U
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
 

def Simple_update(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, chi):

 bdiChi=uni10.Bond(uni10.BD_IN, chi)
 bdoChi=uni10.Bond(uni10.BD_OUT, chi)

 Iso_list=[]*(N)
 for i in xrange(N):
  #print "Step", i
  T0=uni10.UniTensor([mps_A[i].bond(0),bdip,bdi,mps_A[i].bond(2)])
  T0.putBlock(mps_A[i].getBlock())
  T0.setLabel([0,1,2,3])
  #T0.setLabel([0,1,3,2],3)
  Tempo=uni10.UniTensor([bdip,bdoChi])
  Tempo.identity()
  #Tempo.randomize()
  U, s, V=svd_parity(Tempo)
  U.setLabel([1,3])
  U_transpose=copy.copy(U)
  #U_transpose.setLabel([1,2,5,6,90,100 ])
  #print U*U_transpose
  U_transpose.setLabel([4,3])
  #print U*U_transpose
  #cost=Cost_Function(T0, T1, U, U_transpose)
  #Y=Env_U(T0, T1, U, U_transpose)
  #print Y*U
  
  U=optimized_U(T0, U, U_transpose)
  
  Iso_list.append(U)


#  mps_list=[]
#  for i in xrange(0,N,2):
#   #print i, i/2, N
#   T0=uni10.UniTensor([mps_A[i].bond(0),bdi,bdi_phys,bdi,mps_A[i].bond(2)])
#   T0.putBlock(mps_A[i].getBlock())
#   #print T0.printDiagram()
#   T0.permute([0,1,2,3,4],3)
#   
#   T1=uni10.UniTensor([mps_A[i+1].bond(0),bdi,bdi_phys,bdi,mps_A[i+1].bond(2)])
#   T1.putBlock(mps_A[i+1].getBlock())
#   #print T1.printDiagram()
#   #T1.permute([0,1,2,4,3],3)
#   T1.setLabel([4,5,6,8,7])
#   Iso_list[i/2].setLabel([1,2,5,6,9,10])
#   #print "T1", T1.printDiagram()
#   
#   teta=(T0*T1)*Iso_list[i/2]
#   teta.permute([3,0,9,8,7,10],3)
#   #print "teta", teta.printDiagram()
#   chi1=0
#   if i==0:
#     if chi<=(teta.bond(0).dim()*teta.bond(1).dim()*teta.bond(2).dim()):
#        chi1=chi
#     else:
#        chi1=(teta.bond(0).dim()*teta.bond(1).dim()*teta.bond(2).dim())
#   elif i==(N-2):
#     if chi<=(teta.bond(3).dim()*teta.bond(4).dim()*teta.bond(5).dim()):
#       chi1=chi
#     else:
#       chi1=(teta.bond(3).dim()*teta.bond(4).dim()*teta.bond(5).dim())
#   else:
#     if chi<=(teta.bond(0).dim()*teta.bond(1).dim()*teta.bond(2).dim()):
#        chi1=chi
#     else:
#        chi1=(teta.bond(0).dim()*teta.bond(1).dim()*teta.bond(2).dim())
# 
#   #print "i=", i, teta.printDiagram(), "chi1", chi1
#   U, V, s=setTruncation1(teta, chi1)
# 
#   U.setLabel([3,0,9,-1])
#   U.permute([0,9,3,-1],3)
#   U.combineBond([9,3])
#   mps_list.append(U)
# 
#   V.setLabel([-2,8,7,10])
#   s.setLabel([-1,-2])
#   V=s*V
#   V.permute([-1,10,8,7],3)
#   V.combineBond([10,8])
#   mps_list.append(V)
#   #print V.printDiagram()
#  list_bond=[]
#  for q in xrange(N):
#    list_bond.append(mps_list[q].bond(2).dim())
# 
#  mps_R=MPSclass.MPS(mps_list[1].bond(1).dim(),max(list_bond), N,'ortho')
#  #randuni, rand, ortho
#  for i in xrange(N):
#   mps_R[i]=copy.copy(mps_list[i])
#  #print mps_R[9].printDiagram()
#mps_R,
 return  Iso_list




def Fidel_Iso(mps_A, Iso_list, N,bdi,bdi1,bdip,bdo,bdo1,bdop):
 Env_left=[None]*(N)
 for i in xrange( N):
   #print i, 2*i+1
   if i == 0:
      #print "i=inside", i
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
   #print i, 2*i+1
   if i == ((N/2)-1):
      #print "i=initRight", i
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


      Result=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Up))*(T1*T0)
      Result.permute([7,16,15,17],4)
      #print Result.printDiagram()
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
      Result=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U*Up))*(T1*T0))*(Env_right[i+1])
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


      Result=(((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*(T1*T0)
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

      Result=((((t1*(IsoP*IsoP_t))*((t0*(Iso*Iso_t))))*(U))*((T1*T0)*Up))*(Env_left[i-1])
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

  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop)

  Env_left1[0].setLabel([1,2,3,4])
  Env_right1[1].setLabel([1,2,3,4])
  Dominator=Env_left1[0]*Env_right1[1]

  Dominator=(Dominator[0])**(0.5)
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
     print 'NotoptimizedUni=i, E1, E2=', i,'  ', E1, '   ', E2
     Uni_list[Location].putBlock(U_update.getBlock())
     Fidel=E2
     break
    #Y=Env_U(T0, T1, U, U_transpose)
    #print Y.printDiagram() 
    svd=Env_u.getBlock().svd()
    temporary_matrix=svd[0]*svd[2]
    Uni_list[Location].putBlock(temporary_matrix)
    E2=copy.copy(E1)
  print "Location", Location, "FidelUni", Fidel


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

  print "Location", Location, "FidelUni", Fidel






def optimize_iso_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right):
  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].identity()   
  #print "hiiiiiiiiiiiiiii", Uni_list_iden[1]
  Env_left1, Env_right1=Pr_Env_right(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi)
  #Dominator=Env_left1[0]*Env_right1[1]
  #Dominator=(Dominator[0])**(0.5)

  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
  Env_iso.setLabel([1,2,3,4,5,6])
  Iso_list[Location].setLabel([1,2,3,4,5,6])
  R=Iso_list[Location]*Env_iso

  Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
  Env_iso1.setLabel([1,2,3,4,5,6])
  Iso_list[Location].setLabel([1,2,3,4,5,6])
  Dominator=Iso_list[Location]*Env_iso1
  Dominator=(Dominator[0])**(0.5)

  E=R[0]/Dominator
  #print 'E=', E
  E1=10
  Fidel=0
  E2=0.0
  U_update=copy.copy(Iso_list[Location])
  U_first=copy.copy(Iso_list[Location])
  E_previous=0
  for i in xrange(100):
   Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
   Env_iso.setLabel([1,2,3,4,5,6])
   Iso_list[Location].setLabel([1,2,3,4,5,6])
   R=Iso_list[Location]*Env_iso
   Env_iso1=Obtain_Env_Iso(mps_A, Iso_list,Uni_list_iden, N, bdi_phys, bdi, Location,Env_left1, Env_right1)
   Env_iso1.setLabel([1,2,3,4,5,6])
   Iso_list[Location].setLabel([1,2,3,4,5,6])
   Dominator=Iso_list[Location]*Env_iso1

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
  print "Location", Location, "FidelIso", Fidel
 
 
def optimize_iso_locally(mps_A, Iso_list, Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right):

  Uni_list_iden=[  copy.copy(Uni_list[i])  for i in xrange(len(Uni_list)) ]
  for i in xrange(len(Uni_list)):
   Uni_list_iden[i].identity()   

  Iso_list_copy=[  copy.copy(Iso_list[i])  for i in xrange(len(Iso_list)) ]

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
    print "Iso, step", i, E1, E2, abs((E2-E1)/E1)
    if E1>E2 or i is 0:
     U_update=copy.copy(Iso_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      print "Break", E2, E1, abs((E2-E1)/E1), i
      break
    else:
     print 'NotoptimizedIso=i, E1, E2=', i,'  ', E1, '   ', E2
     Iso_list[Location].putBlock(U_update.getBlock())
     Fidel=E2
     break
    Env_iso=(-2.0)*Env_iso+Env_iso1
    svd=Env_iso.getBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Iso_list[Location].putBlock(temporary_matrix)
    Iso_list_copy[Location].putBlock(temporary_matrix)

    E2=copy.copy(E1)

  print "Location", Location, "FidelIso", Fidel


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

  Env_isop1=Obtain_Env_Isop(mps_A, Iso_list,Iso_list_copy,Uni_list_iden, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left1, Env_right1)
  Env_isop1.setLabel([1,2])
  Isop_list[Location].setLabel([1,2])
  Dominator=Isop_list[Location]*Env_isop1
  Dominator=(abs(Dominator[0]))**(0.5)


  for i in xrange(5):
    Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
    Env_isop.setLabel([1,2])
    Isop_list[Location].setLabel([1,2])
    R=Isop_list[Location]*Env_isop


    E1=abs(R[0]/Dominator)
    Fidel=E1
    print "Isop, step", i,  E1,E2, abs((E2-E1)/E1)
    if E1>E2 or i is 0:
     U_update=copy.copy(Isop_list[Location])
     if abs((E2-E1)/E1) < 1.0e-7:
      print "Break", E2, E1, abs((E2-E1)/E1), i
      break
    else:
     print 'NotoptimizedIso=i, E1, E2=', i,'  ', E1, '   ', E2
     Isop_list[Location].putBlock(U_update.getBlock())
     Fidel=E2
     break
    Env_isop=Env_isop
    svd=Env_isop.getBlock().svd()
    temporary_matrix=(-1.0)*svd[0]*svd[2]
    Isop_list[Location].putBlock(temporary_matrix)
    E2=copy.copy(E1)

  print "Location", Location, "FidelIsop", Fidel





def Optimize_uni_Iso(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right,N_iter):


 for q in xrange(2):

  for i in xrange(N):
  #for i in xrange(0):
   #print i
   Location=i
   optimize_iso_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   optimize_isop_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #optimize_iso_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
   #Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi)


 #print Uni_list[0]
 for q in xrange(N_iter):

  for i in xrange(N-1):
    #print ""i
    Location=i
    optimize_uni_locally(mps_A, Iso_list,Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
    Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
    #optimize_uni_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
    #Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi)

  for i in xrange(N):
  #for i in xrange(0):
   #print i
   Location=i
   optimize_iso_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   optimize_isop_locally(mps_A, Iso_list,Isop_list, Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
   Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Isop_list,Uni_list,  N, bdi, bdi1, bdip, bdo, bdo1, bdop)
   #optimize_iso_locally_SD(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi, Location,Env_left, Env_right)
   #Env_left, Env_right=Pr_Env_right(mps_A, Iso_list,Uni_list, N, bdi_phys, bdi)




def Full_update(mps_A,Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop,N_iter):

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

 # for i in xrange(len(Isop_list)):
#   Isop_list[i].randomize()
#  for i in xrange(len(Isop_list)):
#   Iso_list[i].randomize()



 Env_left, Env_right=Pr_Env_right(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)
 for i in xrange(len(Env_left)-1):
  Env_left[i].setLabel([1,2,3,4])
  Env_right[i+1].setLabel([1,2,3,4])
  R=Env_left[i]*Env_right[i+1]
  print "Env", i, (abs(R[0]))**(0.5)

 print "\n"


 for i in xrange(N-1):
  Location=i
  Env_u=Obtain_Env_U(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location,Env_left, Env_right)
  Env_u.setLabel([1,2,3,4])
  Uni_list[Location].setLabel([1,2,3,4])
  R=Uni_list[Location]*Env_u
  print "Env_Uni", Location, (abs(R[0]))**(0.5)


 print "\n"

 for i in xrange(N):
  Location=i
  Env_iso=Obtain_Env_Iso(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_iso.setLabel([1,2])
  Iso_list[Location].setLabel([1,2])
  R=Iso_list[Location]*Env_iso
  print "Env_Iso", Location, (abs(R[0]))**(0.5)

 print "\n"


 for i in xrange(N):
  Location=i
  Env_isop=Obtain_Env_Isop(mps_A, Iso_list,Isop_list,Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right)
  Env_isop.setLabel([1,2])
  Isop_list[Location].setLabel([1,2])
  R=Isop_list[Location]*Env_isop
  print "Env_Isop", Location, (abs(R[0]))**(0.5)

 print "\n"




 Optimize_uni_Iso(mps_A, Iso_list, Isop_list, Uni_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop, Location, Env_left, Env_right,N_iter)


 return  Iso_list, Uni_list 





#################   Canonical Algorithm   #####################
D=3
Dp=6
N=10               #Even
chi=3            #less than DDdd
N_iter=20
Norm_init_A=1.0
##############################################################



bdi=uni10.Bond(uni10.BD_IN, D)
bdo=uni10.Bond(uni10.BD_OUT, D)

bdiChi=uni10.Bond(uni10.BD_IN, chi)
bdoChi=uni10.Bond(uni10.BD_OUT, chi)

bdip=uni10.Bond(uni10.BD_IN, Dp)
bdop=uni10.Bond(uni10.BD_OUT, Dp)

bdi1=uni10.Bond(uni10.BD_IN, 1)
bdo1=uni10.Bond(uni10.BD_OUT, 1)


################PEPS-tensor##############################
mps_A=Init_A_PEPS(N,bdi,bdi1,bdip,bdo,bdo1,bdop)
mps_A=mps_A.normalize()
print "norm-mps_A", mps_A.norm() 
mps_A=mps_A*(Norm_init_A)

Iso_list=Simple_update(mps_A, N, bdi, bdi1, bdip, bdo, bdo1, bdop, chi)


fidel_val=Fidel_Iso(mps_A, Iso_list, N, bdi, bdi1, bdip, bdo, bdo1, bdop)

#print "Fidel", fidel_val

Iso_list, Uni_list=Full_update(mps_A,Iso_list,N, bdi, bdi1, bdip, bdo, bdo1, bdop,N_iter)






