import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
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

def   svd_parity(theta):
    #print theta,theta.getBlock().svd()
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(2).Qlist())

    #A_{m<n}=U_{mm}S_{mm}V_{mn}

    LA=uni10.UniTensor([theta.bond(1), bdo])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1),bdo])
    GB=uni10.UniTensor([theta.bond(1), bdo1])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
    for qnum in blk_qnums:
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0])
        GB.putBlock(qnum, svd[2])
        LA.putBlock(qnum, svd[1])
#    print LA
    return GA, LA,GB

def   qr_parity(theta):

        #bd1=copy.copy(theta.bond(3))
        #bd1.change(uni10.BD_IN)
        bd1=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())

        GA=uni10.UniTensor(uni10.CTYPE,[theta.bond(0),theta.bond(1)])
        LA=uni10.UniTensor(uni10.CTYPE,[bd1, theta.bond(1)])

        svds = {}
        blk_qnums = theta.blockQnum()
        dim_svd=[]
        for qnum in blk_qnums:
                svds[qnum] = theta.getBlock(qnum).qr()
                GA.putBlock(qnum, svds[qnum][0])
                LA.putBlock(qnum, svds[qnum][1])

        #    print LA
        return GA, LA
##########################################################
def    qr_parity1(theta):

        #bd1=copy.copy(theta.bond(3))
        #bd1.change(uni10.BD_IN)
        bd1=uni10.Bond(uni10.BD_IN,theta.bond(2).Qlist())

        GA=uni10.UniTensor(uni10.CTYPE,[theta.bond(0),theta.bond(1),theta.bond(2)])
        LA=uni10.UniTensor(uni10.CTYPE,[bd1, theta.bond(2)])

        svds = {}
        blk_qnums = theta.blockQnum()
        dim_svd=[]
        for qnum in blk_qnums:
                svds[qnum] = theta.getBlock(qnum).qr()
                GA.putBlock(qnum, svds[qnum][0])
                LA.putBlock(qnum, svds[qnum][1])

        #    print LA
        return GA, LA






############### MPS---Square-root #############################
class    MPS:

#################################################################################
#Use the __init__() function to assign values to object properties
#The self parameter is a reference to the class instance itself, and is used to access variables that belongs to the class.
 #@profile
 def __init__(self, physical=2, Dimension=2, Number=2, rand_fuc='rand'):
   self.N = Number
   self.D = Dimension
   self.d = physical
   self.tensor=[None]*Number
   bdi = uni10.Bond(uni10.BD_IN, self.D)
   bdo = uni10.Bond(uni10.BD_OUT, self.D)
   bdi1 = uni10.Bond(uni10.BD_IN, 1)
   bdo1 = uni10.Bond(uni10.BD_OUT, 1)
   bdi_pys = uni10.Bond(uni10.BD_IN, self.d)
   A_fixed=uni10.UniTensor([bdi,bdi_pys,bdo], "A_middle")
   A_fixed.randomize()

   for i in xrange(self.N):
    if i == 0:
     A=uni10.UniTensor([bdi1,bdi_pys,bdo], "A_0")
     if rand_fuc is 'rand' or 'randuni':
      A.randomize()
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.orthoRand()
      self.tensor[i]=A
    elif i ==((self.N)-1):
     A=uni10.UniTensor([bdi,bdi_pys,bdo1], "A_N")
     if rand_fuc is 'rand' or 'randuni':
      A.randomize()
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.orthoRand()
      self.tensor[i]=A
    else:
     A=uni10.UniTensor([bdi,bdi_pys,bdo], "A_middle")
     if rand_fuc is 'rand':
      A.randomize()
      self.tensor[i]=A
     elif rand_fuc is 'ortho':
      A.orthoRand()
      self.tensor[i]=A
     elif rand_fuc is 'randuni':
      self.tensor[i]=copy.copy(A_fixed)


#################################################################################
 def __getitem__(self,i):
        return self.tensor[i]
#################################################################################
 def __setitem__(self, i, A):
        self.tensor[i]=A

#################################################################################
 def norm(self):
   #print len(self.tensor), self.tensor
   A_l=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   E_a=0
   for i in xrange(len(A_l)):
     if i == 0:
       #print A_l[i]
       A_l[i].setLabel([-1,-2,1])
       A_l_dag=copy.copy(A_l[i])
       A_l_dag.setLabel([-1,-2,2])
       E_a=A_l_dag*A_l[i]
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])
     elif i == (len(A_l)-1):
       A_l[i].setLabel([-3,-2,1])
       A_l_dag=copy.copy(A_l[i])
       A_l_dag.setLabel([-4,-2,1])
       #print A_l[i].printDiagram()
       E_a=A_l_dag*(E_a*A_l[i])
     else:
       A_l[i].setLabel([-3,-2,1])
       A_l_dag=copy.copy(A_l[i])
       A_l_dag.setLabel([-4,-2,2])
       E_a=A_l_dag*(E_a*A_l[i])
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])
   

###test# only for N=5###
#   A_l[0].setLabel([1,-1,2])
#   A_l[1].setLabel([2,-2,3])
#   A_l[2].setLabel([3,-3,4])
#   A_l[3].setLabel([4,-4,5])
#   A_l[4].setLabel([5,-5,6])

#   A_l0_dag=copy.copy(A_l[0])
#   A_l1_dag=copy.copy(A_l[1])
#   A_l2_dag=copy.copy(A_l[2])
#   A_l3_dag=copy.copy(A_l[3])
#   A_l4_dag=copy.copy(A_l[4])
  
#   A_l0_dag.setLabel([1,-1,7])
#   A_l1_dag.setLabel([7,-2,8])
#   A_l2_dag.setLabel([8,-3,9])
#   A_l3_dag.setLabel([9,-4,10])
#   A_l4_dag.setLabel([10,-5,6])
   #norm=(A_l[0]*A_l0_dag)*(A_l[1]*A_l1_dag)*(A_l[2]*A_l2_dag)*(A_l[3]*A_l3_dag)*(A_l[4]*A_l4_dag)
#   print norm[0]
   return E_a[0]


##################################################################################
 def __mul__(self, m):
  mps_t=MPS(self.d,self.D,self.N,'randuni')
  for i in xrange(self.N):
    mps_t[i]=copy.copy(self[i])
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

  #print "Max", max(list_bond), max(list_bond1)
  mps_f=MPS(self.d,max(list_bond)+max(list_bond1),self.N,'randuni')

  for q in xrange(self.N):
   if q ==0:
    D=self[q].bond(2).dim()
    D1=mps[q].bond(2).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d,D+D1)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(D+D1):
       if j<D:
        mps_f_mat[i*(D+D1)+j]=self_mat[i*D+j]
       elif j>=D and (j-D)<D1:
        mps_f_mat[i*(D+D1)+j]=mps_mat[i*D1+(j-D)]

    bdi1=uni10.Bond(uni10.BD_IN, 1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, D+D1)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_0")
    mps_f[q].putBlock(mps_f_mat)

   elif q ==(self.N-1):
    D=self[q].bond(0).dim()
    D1=mps[q].bond(0).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d*(D+D1),1)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(D+D1):
       if j<D:
        mps_f_mat[j*mps_f.d+i]=self_mat[j*self.d+i]
       elif j>=D and (j-D)<D1:
        mps_f_mat[j*mps_f.d+i]=mps_mat[(j-D)*mps.d+i]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, 1)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_N")
    mps_f[q].putBlock(mps_f_mat)

   else:
    D=self[q].bond(0).dim()
    Dy=self[q].bond(2).dim()
    D1=mps[q].bond(0).dim()
    Dy1=mps[q].bond(2).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d*(D+D1),Dy1+Dy)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(Dy1+Dy):
      for m in xrange(D+D1):
       if j<Dy and m<D:
        mps_f_mat[m*(Dy1+Dy)*mps_f.d+i*(Dy1+Dy)+j]=self_mat[m*Dy*self.d+i*Dy+j]
       elif j>=Dy and m>=D and (j-Dy)<Dy1 and (m-D)<D1:
        mps_f_mat[m*(Dy1+Dy)*mps_f.d+i*(Dy1+Dy)+j]=mps_mat[(m-D)*Dy1*mps.d+i*Dy1+(j-Dy)]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, Dy1+Dy)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_midle")
    mps_f[q].putBlock(mps_f_mat)

  #print mps_f[0].norm(),self[0].norm()+mps[0].norm()
  #print "Start", mps_f[0].printDiagram(),mps_f[1].printDiagram(),mps_f[2].printDiagram(),mps_f[3].printDiagram()
  return mps_f
###################################################################################
 def __sub__(self, other):
  return self+other*(-1.0)
#################################################################################


#########################################################################
 def fidel(self,mps):
   #print len(self.tensor), self.tensor
   mps_1=MPS(mps.d,mps.D,mps.N,'randuni')
   for i in xrange(self.N):
    mps_1[i]=copy.copy(mps[i])
   self_1=MPS(self.d,self.D,self.N,'randuni')
   for i in xrange(self.N):
    self_1[i]=copy.copy(self[i])

   mps_1=mps_1.normalize()
   self_1=self_1.normalize()

   A_l=[  copy.copy(self_1.tensor[i])  for i in xrange(len(self_1.tensor))  ]
   B_l=[  copy.copy(mps_1.tensor[i])  for i in xrange(len(mps_1.tensor))  ]
   E_a=0
   for i in xrange(len(A_l)):
     if i == 0:
       #print A_l[i]
       A_l[i].setLabel([-1,-2,1])
       #A_l_dag=copy.copy(A_l[i])
       B_l[i].setLabel([-1,-2,2])
       E_a=B_l[i]*A_l[i]
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])
     elif i == (len(A_l)-1):
       A_l[i].setLabel([-3,-2,1])
       #B_l[i]=copy.copy(A_l[i])
       B_l[i].setLabel([-4,-2,1])
       E_a=B_l[i]*(E_a*A_l[i])
     else:
       A_l[i].setLabel([-3,-2,1])
       #A_l_dag=copy.copy(A_l[i])
       B_l[i].setLabel([-4,-2,2])
       E_a=B_l[i]*(E_a*A_l[i])
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])
   

###test# only for N=5###
#   A_l[0].setLabel([1,-1,2])
#   A_l[1].setLabel([2,-2,3])
#   A_l[2].setLabel([3,-3,4])
#   A_l[3].setLabel([4,-4,5])
#   A_l[4].setLabel([5,-5,6])

#   A_l0_dag=copy.copy(A_l[0])
#   A_l1_dag=copy.copy(A_l[1])
#   A_l2_dag=copy.copy(A_l[2])
#   A_l3_dag=copy.copy(A_l[3])
#   A_l4_dag=copy.copy(A_l[4])
  
#   A_l0_dag.setLabel([1,-1,7])
#   A_l1_dag.setLabel([7,-2,8])
#   A_l2_dag.setLabel([8,-3,9])
#   A_l3_dag.setLabel([9,-4,10])
#   A_l4_dag.setLabel([10,-5,6])
   #norm=(A_l[0]*A_l0_dag)*(A_l[1]*A_l1_dag)*(A_l[2]*A_l2_dag)*(A_l[3]*A_l3_dag)*(A_l[4]*A_l4_dag)
#   print norm[0]
   return E_a[0]
##################################################################################


##############################Product###########################################
 def product(self,mps):
   #print len(self.tensor), self.tensor
   mps_1=MPS(mps.d,mps.D,mps.N,'randuni')
   for i in xrange(self.N):
    mps_1[i]=copy.copy(mps[i])
   self_1=MPS(self.d,self.D,self.N,'randuni')
   for i in xrange(self.N):
    self_1[i]=copy.copy(self[i])
   #mps_1=mps_1.normalize()
   #self_1=self_1.normalize()
   
   A_l=[  copy.copy(self_1.tensor[i])  for i in xrange(len(self_1.tensor))  ]
   B_l=[  copy.copy(mps_1.tensor[i])  for i in xrange(len(mps_1.tensor))  ]

   #print A_l[2],B_l[2], 

   E_a=0
   for i in xrange(len(A_l)):
     if i == 0:
       #print A_l[i]
       A_l[i].setLabel([-1,-2,1])
       #A_l_dag=copy.copy(A_l[i])
       B_l[i].setLabel([-1,-2,2])
       E_a=B_l[i]*A_l[i]
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])
     elif i == (len(A_l)-1):
       A_l[i].setLabel([-3,-2,1])
       #B_l[i]=copy.copy(A_l[i])
       B_l[i].setLabel([-4,-2,1])
       E_a=B_l[i]*(E_a*A_l[i])
     else:
       A_l[i].setLabel([-3,-2,1])
       #A_l_dag=copy.copy(A_l[i])
       B_l[i].setLabel([-4,-2,2])
       E_a=B_l[i]*(E_a*A_l[i])
       E_a.permute([1,2],1)
       E_a.setLabel([-3,-4])

   #print E_a
   return E_a[0]


 def distance(self,mps):
   return self.product(self)+mps.product(mps)-2.0*(mps.product(self))

####################################################################################
 def normalize(self):
  mps_t=MPS(self.d,self.D,self.N,'randuni')
  for i in xrange(self.N):
    mps_t[i]=copy.copy(self[i])
  mps_t=mps_t*(1/(mps_t.norm()**(0.5)))
  return mps_t




#####################################################################################
 def appSVD(self,chi):
   assert (self.D>=chi)
   A_l=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   A_l_can=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]

   #print   A_l
   #U,s,V=svd_parity(A_l[0])
   #print self[0].printDiagram(),self[1].printDiagram(),self[2].printDiagram(),self[3].printDiagram()
   for i in xrange(self.N-1):
           #print i
           #print A_l[0].bond(2).dim(), (A_l[0].bond(0).dim()*A_l[0].bond(1).dim())
           row=A_l[i].bond(0).dim()*A_l[i].bond(1).dim()
           colm=A_l[i].bond(2).dim()
           if (row<=colm):
            U,V,s=setTruncation1(A_l[i],row)
           else:
            U,V,s=setTruncation1(A_l[i],colm)
           #print U.printDiagram(),s.printDiagram(),V.printDiagram()
           A_l_can[i]=copy.copy(U)
           s.setLabel([1,2])
           V.setLabel([2,3])
           V=s*V
           V.permute([1,3],1)
           A_l[i+1].setLabel([3,2,4])
           A_l[i+1]=V*A_l[i+1]
           A_l[i+1].permute([1,2,4],2)
           #print Results.printDiagram()
           #Just important for last Step:
           A_l_can[i+1]=copy.copy(A_l[i+1])


   mps_t=MPS( self.d, self.D, self.N, 'randuni')   #randuni, rand, ortho
   for i in xrange(self.N):
     mps_t[i]=copy.copy(A_l_can[i])
     #print  i, mps_t[i].printDiagram() 
#check unitary
     #print A_l_can[self.N-1]
   #print mps_t.norm(), self.norm(), mps_t.fidel(self),mps_t.norm(), self.norm()
   #A_l_can[2].setLabel([1,2,3]);
   #A_l_can_dagg=copy.copy(A_l_can[2]) 
   #A_l_can_dagg.setLabel([1,2,4]);
   #Results=A_l_can_dagg*A_l_can[2]
   #print Results
   #print mps_t

   for i in xrange((self.N-1)):
      A_l_can[self.N-1-i].setLabel([1,2,3])
      A_l_can[self.N-1-i].permute([2,3,1],2)
      #print "\n,\n,\n,\n, First", i, A_l_can[self.N-1-i].printDiagram()#,  (A_l_can[self.N-1-i].bond(0).dim()*A_l_can[self.N-1-i].bond(1).dim()), '\n'

      row=A_l_can[self.N-1-i].bond(0).dim()*A_l_can[self.N-1-i].bond(1).dim()
      colm=A_l_can[self.N-1-i].bond(2).dim()

      if (row<=colm and row<=chi):
           U,V,s=setTruncation1(A_l_can[self.N-1-i],row)
      elif (row<=colm and row>chi):
           U,V,s=setTruncation1(A_l_can[self.N-1-i],chi)
      elif (row>colm and colm<=chi):
           U,V,s=setTruncation1(A_l_can[self.N-1-i],colm)
      elif (row>colm and colm>chi):
           U,V,s=setTruncation1(A_l_can[self.N-1-i],chi)

      row=U.bond(0).dim()*U.bond(1).dim()
      colm=U.bond(2).dim()
      if (row<=colm):
       U1,V1,s1=setTruncation1(U,row)
      else:
       U1,V1,s1=setTruncation1(U,colm)

      U1.setLabel([2,3,1])
      U1.permute([1,2,3],2)
      A_l_can[self.N-1-i]=copy.copy(U1)
      s1.setLabel([1,2])
      V1.setLabel([2,3])
      s.setLabel([3,4])
      V.setLabel([4,5])
      V_f=(s*V)*(s1*V1)
      V_f.permute([1,5],1)
      A_l_can[self.N-i-2].setLabel([-1,2,5])
      #print "i=",i, V_f.printDiagram(), A_l_can[self.N-i-2].printDiagram()
      A_l_can[self.N-i-2]=V_f*A_l_can[self.N-i-2]
      A_l_can[self.N-i-2].permute([-1,2,1],2)
      #print "update", A_l_can[self.N-i-2].printDiagram()

   mps_app=MPS( self.d, chi, self.N, 'randuni')   #randuni, rand, ortho
   for i in xrange(self.N):
     mps_app[i]=copy.copy(A_l_can[i])
     #print mps_app[i].printDiagram()


   #mps_app=mps_app.normalize()
   #print mps_app.norm(),self.norm(), mps_app.fidel(self),mps_app.product(self)#,self.norm()
   #A_l_can[2].setLabel([1,2,3]);
   #A_l_can_dagg=copy.copy(A_l_can[2]) 
   #A_l_can_dagg.setLabel([1,2,4]);
   #Results=A_l_can_dagg*A_l_can[2]
   #print Results
   #print mps_t

   return mps_app
###################################################################
 def appINV(self,chi):
   assert (self.D>=chi)
   A_l=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   A_l_can=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   return mps_app
#####################################################################################

