# Author: Reza Haghshenas rezahaghshenass@gmail.com
# Edit to pyuni10 v2.0: Matt O'Rourke mattorourke41@gmail.com
import pyuni10 as uni10
ut = uni10
import copy
import time
import numpy
np = numpy
import scipy.linalg as linalg
#########  prerequisite functions  #############
def   setTruncation(theta, chi):
    LA=uni10.UniTensorR(theta.bond())
    GA=uni10.UniTensorR(theta.bond())
    GB=uni10.UniTensorR(theta.bond())
    svds = {}
    blk_qnums = theta.BlocksQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = linalg.svd(theta.GetBlock(qnum), full_matrices=False, lapack_driver='gesvd')
        #svds[qnum] = theta.GetBlock(qnum).svd()
        dim_svd.append(int(len(svds[qnum][1])))
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
    GA.Assign([theta.bond(0), theta.bond(1),bdo_mid])
    GB.Assign([bdi_mid,  theta.bond(2)])
    LA.Assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.PutBlock(qnum, svd[0][:,:dim])
        GB.PutBlock(qnum, svd[2][:dim,:])
        LA.PutBlock(qnum, numpy.diag(svd[1][:dim])  )
    return GA, GB, LA

def setTruncationL(theta,chi):
    U=uni10.UniTensorR(theta.bond())
    S=uni10.UniTensorR(theta.bond())
    V=uni10.UniTensorR(theta.bond())
    svds = {}
    blk_qnums = theta.BlocksQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = linalg.svd(theta.GetBlock(qnum))
        #svds[qnum] = theta.GetBlock(qnum).svd()
        dim_svd.append(int(len(svds[qnum][1])))
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
    U.Assign([theta.bond(0),bdo_mid])
    V.Assign([bdi_mid, theta.bond(1), theta.bond(2)])
    S.Assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        U.PutBlock(qnum, svd[0][:,:dim])
        V.PutBlock(qnum, svd[2][:dim,:])
        S.PutBlock(qnum, numpy.diag(svd[1][:dim])  )
    return U, V, S

def   sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + len(sv_mat)
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
            if(cur1 < len(ori_svs)) and cur2 < len(sv_mat):
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
            elif cur2 < len(sv_mat) :
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
        bidxs = [bidx] * len(sv_mat)
        svs = [sv_mat[i] for i in xrange(len(sv_mat))]
       else: bidxs = [bidx];  svs = [sv_mat[0]];
    return svs, bidxs



############### MPS---Square-root #############################
class  MPS:

#################################################################################
 def __init__(self, physical=2, Dimension=2, Number=2, rand_fuc='rand',seed_val=1):
   self.N = Number
   self.D = Dimension
   self.d = physical
   self.tensor=numpy.zeros(Number, dtype=object)
   bdi = uni10.Bond(uni10.BD_IN, self.D)
   bdo = uni10.Bond(uni10.BD_OUT, self.D)
   bdi1 = uni10.Bond(uni10.BD_IN, 1)
   bdo1 = uni10.Bond(uni10.BD_OUT, 1)
   bdi_pys = uni10.Bond(uni10.BD_IN, self.d)
   A_fixed=uni10.UniTensorR([bdi,bdi_pys,bdo], "A_middle")
   A_fixed.Randomize(dn_mu=-1,up_var=1)
   #print seed_val 
   for i in xrange(self.N):
    if i == 0:
     A=uni10.UniTensorR([bdi1,bdi_pys,bdo], "A_0")
     if rand_fuc is 'rand' or 'randuni':
      A.Randomize(dn_mu=-1,up_var=1,seed=seed_val)
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.OrthoRand()
      self.tensor[i]=A
    elif i ==((self.N)-1):
     A=uni10.UniTensorR([bdi,bdi_pys,bdo1], "A_N")
     if rand_fuc is 'rand' or 'randuni':
      A.Randomize(dn_mu=-1,up_var=1,seed=seed_val)
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.OrthoRand()
      self.tensor[i]=A
    else:
     A=uni10.UniTensorR([bdi,bdi_pys,bdo], "A_middle")
     if rand_fuc is 'rand':
      A.Randomize(dn_mu=-1,up_var=1,seed=seed_val)
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.OrthoRand()
      self.tensor[i]=A
     elif rand_fuc is 'randuni':
      self.tensor[i]=ut.UniTensorR(A_fixed)


#################################################################################
 def __getitem__(self,i):
        return self.tensor[i]
#################################################################################
 def __setitem__(self, i, A):
        self.tensor[i]=A

#################################################################################
 def norm(self):
   A_l=[ self.tensor[i]*1.0  for i in xrange(len(self.tensor))  ]
   E_a=0

   for i in xrange(len(A_l)):
     if i == 0:
       A_l[i].SetLabel([-1,-2,1])
       A_l_dag=A_l[i]*1.0
       A_l_dag.SetLabel([-1,-2,2])
       E_a = ut.Contract(A_l_dag, A_l[i])
       E_a = ut.Permute(E_a, [1,2],1)
       E_a.SetLabel([-3,-4])
     elif i == (len(A_l)-1):
       A_l[i].SetLabel([-3,-2,1])
       A_l_dag=A_l[i]*1.0
       A_l_dag.SetLabel([-4,-2,1])
       E_a = ut.Contract(E_a, A_l[i])
       E_a = ut.Contract(A_l_dag, E_a)
     else:
       A_l[i].SetLabel([-3,-2,1])
       A_l_dag=A_l[i]*1.0
       A_l_dag.SetLabel([-4,-2,2])
       E_a = ut.Contract(E_a, A_l[i])
       E_a = ut.Contract(A_l_dag, E_a)
       E_a = ut.Permute(E_a, [1,2],1)
       E_a.SetLabel([-3,-4])

   ret = E_a.GetRawElem()
   assert ret.shape == (1,1)
   return ret[0,0]


##################################################################################
 def __mul__(self, m):
  mps_t=MPS(self.d,self.D,self.N,'randuni')
  for i in xrange(self.N):
    mps_t[i]=self[i]*1.0
  for i in xrange(mps_t.N):
   if m<0:
    mps_t[i]=mps_t[i]*(abs(m)**(1.0/mps_t.N))
    if i==0:
     mps_t[0]=mps_t[0]*(-1.0)
   else:
    mps_t[i]=mps_t[i]*(abs(m)**(1.0/mps_t.N))
  return mps_t

#################################################################################
 def __add__(self, mps):
  list_bond=[]
  list_bond1=[]

  for q in xrange(self.N):
   list_bond.append(self[q].bond(2).dim())

  for q in xrange(self.N):
   list_bond1.append(mps[q].bond(2).dim())

  mps_f=MPS(self.d,max(list_bond)+max(list_bond1),self.N,'randuni')

  for q in xrange(self.N):
   if q ==0:
    D=self[q].bond(2).dim()
    D1=mps[q].bond(2).dim()
    self_mat=self[q].GetBlock()
    mps_mat=mps[q].GetBlock()
    bdi_t = uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo_t = uni10.Bond(uni10.BD_OUT, D+D1)
    mps_f_t=uni10.UniTensorR([bdi_t,bdo_t])
    mps_f_mat=mps_f_t.GetBlock()
    mps_f_mat.fill(0)
    for i in xrange(mps_f.d):
     for j in xrange(D+D1):
       if j<D:
        mps_f_mat[i][j]=self_mat[i][j]
       elif j>=D and (j-D)<D1:
        mps_f_mat[i][j]=mps_mat[i][j-D]

    bdi1=uni10.Bond(uni10.BD_IN, 1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, D+D1)

    mps_f[q]=uni10.UniTensorR([bdi1,bdi,bdo], "A_0")
    mps_f[q].PutBlock(mps_f_mat)

   elif q ==(self.N-1):
    D=self[q].bond(0).dim()
    D1=mps[q].bond(0).dim()
    self_mat=self[q].GetBlock()
    mps_mat=mps[q].GetBlock()

    bdi_t = uni10.Bond(uni10.BD_IN, mps_f.d*(D+D1))
    bdo_t = uni10.Bond(uni10.BD_OUT, 1)
    mps_f_t=uni10.UniTensorR([bdi_t,bdo_t])
    mps_f_mat=mps_f_t.GetBlock()
    mps_f_mat.fill(0)

    for i in xrange(mps_f.d*(D+D1)):
     for j in xrange(1):
       if i<(D*mps_f.d):
        mps_f_mat[i][j]=self_mat[i][j]
       elif i>=(D*mps_f.d) and (i-(D*mps_f.d))<(D1*mps_f.d):
        mps_f_mat[i][j]=mps_mat[i-(D*mps_f.d)][j]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, 1)

    mps_f[q]=uni10.UniTensorR([bdi1,bdi,bdo], "A_N")
    mps_f[q].PutBlock(mps_f_mat)

   else:
    D=self[q].bond(0).dim()
    Dy=self[q].bond(2).dim()
    D1=mps[q].bond(0).dim()
    Dy1=mps[q].bond(2).dim()
    self_mat=self[q].GetBlock()
    mps_mat=mps[q].GetBlock()

    bdi_t = uni10.Bond(uni10.BD_IN, mps_f.d*(D+D1))
    bdo_t = uni10.Bond(uni10.BD_OUT, Dy1+Dy)
    mps_f_t=uni10.UniTensorR([bdi_t,bdo_t])
    mps_f_mat=mps_f_t.GetBlock()
    mps_f_mat.fill(0)

    for i in xrange(mps_f.d*(D+D1)):
     for j in xrange(Dy1+Dy):
       if i<mps_f.d*D and j<Dy:
        mps_f_mat[i][j]=self_mat[i][j]
       elif i>=mps_f.d*D and j>=Dy and (i-mps_f.d*D)<mps_f.d*D1 and (j-Dy)<Dy1:
        mps_f_mat[i][j]=mps_mat[i-mps_f.d*D][j-Dy]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, Dy1+Dy)

    mps_f[q]=uni10.UniTensorR([bdi1,bdi,bdo], "A_midle")
    mps_f[q].PutBlock(mps_f_mat)

  return mps_f
###################################################################################
 def __sub__(self, other):
  return self+other*(-1.0)
#################################################################################


#########################################################################
 def fidel(self,mps):
   mps_norm=mps.norm()
   self_norm=self.norm()
   assert ( mps_norm>1.e-14), "MPS_norm is too small"           
   assert ( self_norm>1.e-14), "MPS_norm is too small"         

   E_a=0
   for i in xrange(mps.N):
     if i == 0:
       self[i].SetLabel([-1,-2,1])
       mps[i].SetLabel([-1,-2,2])
       E_a= ut.Contract(mps[i], self[i])
       E_a = ut.Permute(E_a, [1,2],1)
       E_a.SetLabel([-3,-4])
     elif i == (mps.N-1):
       self[i].SetLabel([-3,-2,1])
       mps[i].SetLabel([-4,-2,1])
       E_a = ut.Contract(E_a, self[i])
       E_a = ut.Contract(mps[i], E_a)
     else:
       self[i].SetLabel([-3,-2,1])
       mps[i].SetLabel([-4,-2,2])
       E_a = ut.Contract(E_a, self[i])
       E_a = ut.Contract(mps[i], E_a)
       E_a = ut.Permute(E_a, [1,2],1)
       E_a.SetLabel([-3,-4])

   return E_a.GetBlock()[0,0] / (numpy.sqrt(mps_norm) * numpy.sqrt(self_norm))
##################################################################################


##############################Product###########################################
 def product(self,mps):

   E_a=0
   for i in xrange(self.N):
     if i == 0:
       self[i].SetLabel([-1,-2,1])
       mps[i].SetLabel([-1,-2,2])
       E_a= ut.Contract(mps[i], self[i])
       E_a = ut.Permute(E_a,[1,2],1)
       E_a.SetLabel([-3,-4])
     elif i == (self.N-1):
       self[i].SetLabel([-3,-2,1])
       mps[i].SetLabel([-4,-2,1])
       E_a = ut.Contract(E_a,self[i])
       E_a = ut.Contract(mps[i],E_a)
     else:
       self[i].SetLabel([-3,-2,1])
       mps[i].SetLabel([-4,-2,2])
       E_a = ut.Contract(E_a, self[i])
       E_a = ut.Contract(mps[i],E_a)
       E_a = ut.Permute(E_a,[1,2],1)
       E_a.SetLabel([-3,-4])

   return E_a.GetBlock()[0,0]


 def distance(self,mps):
   return self.product(self)+mps.product(mps)-2.0*(mps.product(self))

####################################################################################
 def normalize(self):
  mps_t=MPS(self.d,self.D,self.N,'randuni')
  for i in xrange(self.N):
    mps_t[i]=self[i]*1.0
  mps_t=mps_t*(1/(mps_t.norm()**(0.5)))
  return mps_t




#####################################################################################
 def appSVD(self,chi):
   assert (self.D>=chi)
   A_l_can= numpy.zeros(self.N, dtype=object)
   A_l=[  self.tensor[i]*1.0  for i in xrange(len(self.tensor))  ]
   A_l_can=[  self.tensor[i]*1.0  for i in xrange(len(self.tensor))  ]

   for i in xrange(self.N-1):
           row=A_l[i].bond(0).dim()*A_l[i].bond(1).dim()
           colm=A_l[i].bond(2).dim()
            
           if (row<=colm):
            U,V,s=setTruncation(A_l[i],row)
           else:
            U,V,s=setTruncation(A_l[i],colm)
           A_l_can[i]=U*1.0
           s.SetLabel([1,2])
           V.SetLabel([2,3])
           V=ut.Contract(s , V)
           V=ut.Permute(V,[1,3],1)
           A_l[i+1].SetLabel([3,2,4])
           A_l[i+1]=ut.Contract(V , A_l[i+1])
           A_l[i+1]=ut.Permute(A_l[i+1],[1,2,4],2)
           A_l_can[i+1]=A_l[i+1]*1.0


   for i in xrange((self.N-1)):
      A_l_can[self.N-1-i].SetLabel([1,2,3])
      A_l_can[self.N-1-i] = ut.Permute(A_l_can[self.N-1-i], [1,2,3],1)

      row=A_l_can[self.N-1-i].bond(0).dim()
      colm=A_l_can[self.N-1-i].bond(2).dim()*A_l_can[self.N-1-i].bond(1).dim()

      if (row<=colm and row<=chi):
           U,V,s=setTruncationL(A_l_can[self.N-1-i],row)
      elif (row<=colm and row>chi):
           U,V,s=setTruncationL(A_l_can[self.N-1-i],chi)
      elif (row>colm and colm<=chi):
           U,V,s=setTruncationL(A_l_can[self.N-1-i],colm)
      elif (row>colm and colm>chi):
           U,V,s=setTruncationL(A_l_can[self.N-1-i],chi)

      V.SetLabel([1,2,3])
      V = ut.Permute(V, [1,2,3],2)
      A_l_can[self.N-1-i]= ut.UniTensorR(V)
      s.SetLabel([2,3])
      U.SetLabel([1,2])
      U = ut.Contract( U, s)
      A_l_can[self.N-i-2].SetLabel([-2,-1,1])
      A_l_can[self.N-i-2] = ut.Contract(A_l_can[self.N-i-2] , U)
      A_l_can[self.N-i-2] = ut.Permute(A_l_can[self.N-i-2], [-2,-1,3],2)

   mps_app=MPS( self.d, chi, self.N)   #randuni, rand, ortho
   for i in xrange(self.N):
     mps_app[i]= A_l_can[i]*1.0

   return mps_app
#####################################################################################
