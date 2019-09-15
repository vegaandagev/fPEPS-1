import pyuni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass2



#mps_a=MPSclass.MPS(4,64,4,'ortho')   #randuni, rand, ortho
#mps_b=MPSclass.MPS(4,20,4,'rand')   #randuni, rand, ortho

#print mps_a.N, mps_a.D, mps_a.d#, mps_a.tensor   
#print MPS.__doc__,MPS.__name__,MPS.__module__,MPS.__bases__,MPS.__dict__
#print mps_a[0], mps_a[1],mps_a[2],mps_a[3]
#print mps_a[0]
#mps_a[0]=2
#print mps_a[0]

#print  mps_a
# norm=mps_a.norm()
# print norm
#mps_a_=mps_a*-2.0     #uniform distribution
#norm=mps_a_.norm()
# print norm/4.0
# print mps_a_[0], mps_a[0]
# print mps_a_[1], mps_a[1]

#mps_a=mps_a*(1/(mps_a.norm()**(0.5)))
#mps_b=mps_b*(1/(mps_b.norm()**(0.5)))
##print mps_a.norm(),mps_b.norm()
#mps_ab=mps_a+mps_b
#mps_ab=(mps_a*3-mps_b*4)*0.5

#fidel=mps_a.fidel(mps_b)
##fidel=mps_a.fidel(mps_a)
#print mps_a.product(mps_b)
#mps_a=mps_a*(-1)
#print mps_a.product(mps_b)
#print mps_a
#print mps_ab.D, mps_ab.norm(), (mps_a.norm()*9+mps_b.norm()*16-12*2*mps_a.product(mps_b))*0.25
#mps_a=mps_a*(1/(mps_a.norm()**(0.5)))
#mps_a=mps_a.normalize()
#chi=2
#print mps_a.norm()
#mps_a_app=mps_a.appSVD(chi)


#print mps_a.norm(),mps_b.norm(), mps_a.product(mps_b)


#print mps_a_app.norm(),mps_a.norm(), mps_a_app.product(mps_a),mps_a_app.fidel(mps_a)

#print mps_a.distance(mps_a_app)


#mps_a_app=mps_a.appINV(chi)


def Init_A_PEPS( A_list, N, bdi_phys, bdi, bdi1, bdo, bdo1):
 for i in xrange(N):
  if i == 0:
    A=uni10.UniTensorR([bdi_phys,bdi1,bdo,bdo], "A_middle")
    A.randomize()
    #A.orthoRand()
    print A
    A.SetLabel([0,1,2,3])
    A_dagg=copy.copy(A)
    A_dagg.SetLabel([0,-1,-2,-3])
    A=A*A_dagg
    A.uni10.Permute([1,-1,2,-2,3,-3],4)
    A.CombineBond([1,-1])
    A.CombineBond([3,-3])
    A.CombineBond([2,-2])
    A_list[i]=A
  elif i ==(N-1):
    A=uni10.UniTensorR([bdi_phys,bdi,bdo,bdo1], "A_middle")
    A.randomize()
    #A.orthoRand()

    A.SetLabel([0,1,2,3])
    A_dagg=copy.copy(A)
    A_dagg.SetLabel([0,-1,-2,-3])
    A=A*A_dagg
    A.uni10.Permute([1,-1,2,-2,3,-3],4)
    A.CombineBond([1,-1])
    A.CombineBond([3,-3])
    A.CombineBond([2,-2])
    A_list[i]=A
  else:
    A=uni10.UniTensorR([bdi_phys,bdi,bdo,bdo], "A_middle")
    A.randomize()
    #A.orthoRand()

    A.SetLabel([0,1,2,3])
    A_dagg=copy.copy(A)
    A_dagg.SetLabel([0,-1,-2,-3])
    A=A*A_dagg
    A.uni10.Permute([1,-1,2,-2,3,-3],4)
    A.CombineBond([1,-1])
    A.CombineBond([3,-3])
    A.CombineBond([2,-2])
    A_list[i]=A
 mps_A=MPSclass.MPS( A_list[1].bond(1).dim(), A_list[1].bond(0).dim(), N, 'ortho')
 for i in xrange(N):
   mps_A[i]=copy.copy(A_list[i])
 return mps_A

def make_YZY_mps(mps_Y,mps_Z):
 #print bdi
 Y_list=[]
 for i in xrange(N):
  A=copy.copy(mps_Y[i])
  B=copy.copy(mps_Z[i])
  C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  C1.SetLabel([7,6,8,9])

  result=(A1*B1)*C1
  result.uni10.Permute([7,4,1,2,8,9,5,3],2)
  result.CombineBond([7,4])
  result.CombineBond([7,1])
  result.CombineBond([9,5])
  result.CombineBond([9,3])
  result.uni10.Permute([7,2,8,9],2)
  result.CombineBond([2,8])
  result.uni10.Permute([7,2,9],2)
  Y_list.append(result)

 mps_YZY=MPSclass.MPS(Y_list[1].bond(1).dim(),Y_list[1].bond(0).dim(),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_YZY[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 return mps_YZY

def make_ZYZ_mps(mps_Y,mps_Z):
 Y_list=[]
 for i in xrange(N):
  A=copy.copy(mps_Z[i])
  B=copy.copy(mps_Y[i])
  C=copy.copy(mps_Z[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  C1.SetLabel([7,6,8,9])

  result=(A1*B1)*C1
  result.uni10.Permute([7,4,1,2,8,9,5,3],2)
  result.CombineBond([7,4])
  result.CombineBond([7,1])
  result.CombineBond([9,5])
  result.CombineBond([9,3])
  result.uni10.Permute([7,2,8,9],2)
  result.CombineBond([2,8])
  result.uni10.Permute([7,2,9],2)
  Y_list.append(result)
 mps_ZYZ=MPSclass.MPS(Y_list[1].bond(1).dim(),Y_list[1].bond(0).dim(),N,'ortho')
 for i in xrange(N):
   mps_ZYZ[i]=copy.copy(Y_list[i])
 return mps_ZYZ

#####################################################################################
def Appmake_YZY_mps(mps_Y,mps_Z,chi,bdi):
 N=mps_Y.N
 Y_list=[None]*mps_Y.N
 for i in xrange(mps_Y.N):
  A=mps_Y[i]*1.0
  B=mps_Z[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result


 list_bond=[]
 for q in xrange(mps_Y.N):
   list_bond.append(Y_list[q].bond(2).dim())

# print list_bond
 mps_ZY=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),mps_Y.N)
 #randuni, rand, ortho
 for i in xrange(mps_Y.N):
  mps_ZY[i]=Y_list[i]*1.0
 #print mps_YZY[1].printDiagram()
 
 if mps_ZY.D>chi:
   #print "used chi for multiolication=", chi 
   mps_ZY=mps_ZY.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_ZY.D 
   mps_ZY=mps_ZY.appSVD(mps_ZY.D)


 for i in xrange(mps_Y.N):
  A=mps_ZY[i]*1.0
  B=mps_Y[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(mps_Y.N):
   list_bond.append(Y_list[q].bond(2).dim())

 mps_YZY=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),N)
 #randuni, rand, ortho
 for i in xrange(mps_Y.N):
  mps_YZY[i]=Y_list[i]*1.0
 #print mps_YZY[1].printDiagram()
 
 if mps_YZY.D>chi:
   #print "used chi for multiolication=", chi 
   mps_YZY=mps_YZY.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_YZY.D 
   mps_YZY=mps_YZY.appSVD(mps_YZY.D)
 return mps_YZY



def Appmake_ZYZ_mps(mps_Y,mps_Z,chi,bdi):
 N=mps_Y.N
 Y_list=[None]*N
 for i in xrange(N):
  A=mps_Z[i]*1.0
  B=mps_Y[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())


 mps_YZ=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),N)
 #randuni, rand, ortho
 for i in xrange(N):
  mps_YZ[i]=Y_list[i]*1.0
 #print mps_YZY[1].printDiagram()
 
 if mps_YZ.D>chi:
   #print "used chi for multiolication=", chi 
   mps_YZ=mps_YZ.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_YZ.D 
   mps_YZ=mps_YZ.appSVD(mps_YZ.D)

 for i in xrange(N):
  A=mps_YZ[i]*1.0
  B=mps_Z[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())



 mps_ZYZ=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),N)
 #randuni, rand, ortho
 for i in xrange(N):
  mps_ZYZ[i]=Y_list[i]*1.0
 #print mps_YZY[1].printDiagram()
 
 if mps_ZYZ.D>chi:
   #print "used chi for multiolication=", chi 
   mps_ZYZ=mps_ZYZ.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_ZYZ.D 
   mps_ZYZ=mps_ZYZ.appSVD(mps_ZYZ.D)
 
 return mps_ZYZ
###################################################################################


def Appmake_YV_mps(mps_A,mps_V,chi):

 Y_list=[None]*N
 for i in xrange(N):
  A=copy.copy(mps_A[i])
  B=copy.copy(mps_V[i])
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=(A1*B1)#*C1
  result.uni10.Permute([4,1,2,6,5,3],2)
  result.CombineBond([4,1])
  result.CombineBond([5,3])
  result.uni10.Permute([4,2,6,5],2)
  result.CombineBond([2,6])
  result.uni10.Permute([4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())


 mps_AV=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_AV[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 
 if mps_AV.D>chi:
   #print "used chi for multiolication=", chi 
   mps_AV=mps_AV.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_AV.D 
   mps_AV=mps_AV.appSVD(mps_AV.D)


 return mps_AV







#####################################################################################




def test(mps_Z,mps_Y):
 A_list=[None]*N
 for i in xrange(N):
  A_list[i]=uni10.UniTensorR([mps_Z[i].bond(0),bdi,bdi,mps_Z[i].bond(2)])
  A_list[i].PutBlock(copy.copy(mps_Z[i].GetBlock()))

 B_list=[None]*N
 for i in xrange(N):
  B_list[i]=uni10.UniTensorR([mps_Y[i].bond(0),bdi,bdi,mps_Y[i].bond(2)])
  B_list[i].PutBlock(copy.copy(mps_Y[i].GetBlock()))

 A_list[0].SetLabel([-100,0,1,-2])
 A_list[1].SetLabel([-2,3,4,-3])
 A_list[2].SetLabel([-3,5,6,-4])
 A_list[3].SetLabel([-4,7,8,-200])
 Z_mat=(A_list[0]*A_list[1])*(A_list[2]*A_list[3])
 #print Z_mat.printDiagram()
 Z_mat.uni10.Permute([-100,0,3,5,7,1,4,6,8,-200],5)
 Z_mat.CombineBond([0,-100])
 Z_mat.CombineBond([0,3])
 Z_mat.CombineBond([0,5])
 Z_mat.CombineBond([0,7])
 Z_mat.CombineBond([1,4])
 Z_mat.CombineBond([1,6])
 Z_mat.CombineBond([1,8])
 Z_mat.CombineBond([1,-200])
 Z_mat=Z_mat.GetBlock()

 B_list[0].SetLabel([-100,0,1,-2])
 B_list[1].SetLabel([-2,3,4,-3])
 B_list[2].SetLabel([-3,5,6,-4])
 B_list[3].SetLabel([-4,7,8,-200])
 Y_mat=(B_list[0]*B_list[1])*(B_list[2]*B_list[3])
 Y_mat.uni10.Permute([-100,0,3,5,7,1,4,6,8,-200],5)
 Y_mat.CombineBond([0,-100])
 Y_mat.CombineBond([0,3])
 Y_mat.CombineBond([0,5])
 Y_mat.CombineBond([0,7])
 Y_mat.CombineBond([1,4])
 Y_mat.CombineBond([1,6])
 Y_mat.CombineBond([1,8])
 Y_mat.CombineBond([1,-200])
 Y_mat=Y_mat.GetBlock()
 #print Y_mat,
 print Z_mat*Y_mat


def test1(mps_Z,mps_Y):
 A_list=[None]*N
 for i in xrange(N):
  A_list[i]=uni10.UniTensorR([mps_Z[i].bond(0),bdi,bdi,mps_Z[i].bond(2)])
  A_list[i].PutBlock(copy.copy(mps_Z[i].GetBlock()))

 B_list=[None]*N
 for i in xrange(N):
  B_list[i]=uni10.UniTensorR([mps_Y[i].bond(0),bdi,bdi,mps_Y[i].bond(2)])
  B_list[i].PutBlock(copy.copy(mps_Y[i].GetBlock()))

 A_list[0].SetLabel([-100,0,1,-2])
 A_list[1].SetLabel([-2,3,4,-200])
 Z_mat=(A_list[0]*A_list[1])
 #print Z_mat.printDiagram()
 Z_mat.uni10.Permute([-100,0,3,1,4,-200],3)
 Z_mat.CombineBond([0,-100])
 Z_mat.CombineBond([0,3])
 Z_mat.CombineBond([1,4])
 Z_mat.CombineBond([1,-200])
 Z_mat=Z_mat.GetBlock()

 B_list[0].SetLabel([-100,0,1,-2])
 B_list[1].SetLabel([-2,3,4,-200])
 Y_mat=(B_list[0]*B_list[1])
 Y_mat.uni10.Permute([-100,0,3,1,4,-200],3)
 Y_mat.CombineBond([0,-100])
 Y_mat.CombineBond([0,3])
 Y_mat.CombineBond([1,4])
 Y_mat.CombineBond([1,-200])
 Y_mat=Y_mat.GetBlock()
 #print Y_mat,
 print Z_mat*Y_mat

def Fidel(A, B):
 E=A*B
 
 N=E.trace()
 D=(A.norm()*A.norm()*B.norm()*B.norm())**(0.5)
 #print "fidelity:",B.norm(),E.norm(), N/D
 return N/D


def FidelMPSYZ(mps_Z, mps_Y, mps_Iden,bdi):
 N=mps_Y.N
 Y_list=[None]*N
 for i in xrange(N):
  A=copy.copy(mps_Y[i])
  B=copy.copy(mps_Z[i])
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=(A1*B1)#*C1
  result.uni10.Permute([4,1,2,6,5,3],2)
  result.CombineBond([4,1])
  result.CombineBond([5,3])
  result.uni10.Permute([4,2,6,5],2)
  result.CombineBond([2,6])
  result.uni10.Permute([4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())

# print list_bond
 mps_ZY=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_ZY[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 
# if mps_ZY.D>chi:
#   print "used chi for multiolication=", chi 
#   mps_ZY=mps_ZY.appSVD(chi)
# else:
#   print "used chi for multiplication=",  mps_ZY.D 
#   mps_ZY=mps_ZY.appSVD(mps_ZY.D)

 return mps_ZY.fidel(mps_Iden)



def FidelMPSZAZ(mps_A, mps_Z, mps_Iden,chi, bdi):
 N=mps_A.N
 Y_list=[None]*N
 for i in xrange(N):
  A=copy.copy(mps_Z[i])
  B=copy.copy(mps_A[i])
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=(A1*B1)#*C1
  result.uni10.Permute([4,1,2,6,5,3],2)
  result.CombineBond([4,1])
  result.CombineBond([5,3])
  result.uni10.Permute([4,2,6,5],2)
  result.CombineBond([2,6])
  result.uni10.Permute([4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())


 mps_AZ=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_AZ[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 
 if mps_AZ.D>chi:
   #print "used chi for test=", chi 
   mps_AZ=mps_AZ.appSVD(chi)
 else:
   #print "used chi for test=",  mps_AZ.D 
   mps_AZ=mps_AZ.appSVD(mps_AZ.D)

 for i in xrange(N):
  A=copy.copy(mps_AZ[i])
  B=copy.copy(mps_Z[i])
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=(A1*B1)#*C1
  result.uni10.Permute([4,1,2,6,5,3],2)
  result.CombineBond([4,1])
  result.CombineBond([5,3])
  result.uni10.Permute([4,2,6,5],2)
  result.CombineBond([2,6])
  result.uni10.Permute([4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())



 mps_ZAZ=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_ZAZ[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 
 if mps_ZAZ.D>chi:
   #print "used chi for multiolication=", chi 
   mps_ZAZ=mps_ZAZ.appSVD(chi)
 else:
   #print "used chi for multiplication=",  mps_ZAZ.D 
   mps_ZAZ=mps_ZAZ.appSVD(mps_ZAZ.D)
 
 return mps_ZAZ.fidel(mps_Iden)




def FidelMPSAYY( mps_A, mps_Y, mps_Iden, chi, bdi):
 N=mps_A.N
 Y_list=[None]*N
 for i in xrange(N):
  A=mps_Y[i]*1.0
  B=mps_Y[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())


 mps_YY=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),N)
 #randuni, rand, ortho
 for i in xrange(N):
  mps_YY[i]=Y_list[i]*1.0
 #print mps_YZY[1].printDiagram()
 
 if mps_YY.D>chi:
   #print "used chi for test=", chi 
   mps_YY=mps_YY.appSVD(chi)
 else:
   #print "used chi for test=",  mps_YY.D 
   mps_YY=mps_YY.appSVD(mps_YY.D)

 
 return mps_YY.fidel(mps_A)


def FidelMPSYV(mps_Y,mps_V, mps_Iden, chi, bdi):
 N=mps_Y.N
 Y_list=[None]*N
 for i in xrange(N):
  A=copy.copy(mps_Y[i])
  B=copy.copy(mps_V[i])
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=(A1*B1)#*C1
  result.uni10.Permute([4,1,2,6,5,3],2)
  result.CombineBond([4,1])
  result.CombineBond([5,3])
  result.uni10.Permute([4,2,6,5],2)
  result.CombineBond([2,6])
  result.uni10.Permute([4,2,5],2)
  Y_list[i]=result

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())


 mps_YY=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')
 #randuni, rand, ortho
 for i in xrange(N):
  mps_YY[i]=copy.copy(Y_list[i])
 #print mps_YZY[1].printDiagram()
 
 if mps_YY.D>chi:
   #print "used chi for test=", chi 
   mps_YY=mps_YY.appSVD(chi)
 else:
   #print "used chi for test=",  mps_YY.D 
   mps_YY=mps_YY.appSVD(mps_YY.D)


 return mps_YY.fidel(mps_Iden)

def YY_multy(mps_Y_best,bdi):
 N=mps_Y_best.N
 Y_list=[None]*N
 for i in xrange(N):
  A=mps_Y_best[i]*1.0
  B=mps_Y_best[i]*1.0
  #C=copy.copy(mps_Y[i])
  A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
  B1=uni10.UniTensorR([B.bond(0),bdi,bdi,B.bond(2)])
  #C1=uni10.UniTensorR([C.bond(0),bdi,bdi,C.bond(2)])
  A1.PutBlock(A.GetBlock())
  B1.PutBlock(B.GetBlock())
  #C1.PutBlock(C.GetBlock())

  A1.SetLabel([1,2,-2,3])
  B1.SetLabel([4,-2,6,5])
  #C1.SetLabel([7,6,8,9])

  result=uni10.Contract(A1,B1)#*C1
  result=uni10.Permute(result,[4,1,2,6,5,3],2)
  result=result.CombineBond([4,1])
  result=result.CombineBond([5,3])
  result=uni10.Permute(result,[4,2,6,5],2)
  result=result.CombineBond([2,6])
  result=uni10.Permute(result,[4,2,5],2)
  Y_list[i]=result
  #if i==0 or i==1:
    #print "in", Y_list[i].printDiagram()

 list_bond=[]
 for q in xrange(N):
   list_bond.append(Y_list[q].bond(2).dim())

 mps_YY=MPSclass2.MPS(Y_list[1].bond(1).dim(),max(list_bond),N)

 for i in xrange(N):
  mps_YY[i]=Y_list[i]*1.0
 return  mps_YY
#################   Canonical Algorithm   #####################
def  make_mps_R_root(mps_A_init,N_y, bdi, bdip, bdo, bdop,chi,N_iter):
 A_list=[]
 for i in xrange(mps_A_init.N):
  t0=uni10.UniTensorR([mps_A_init[i].bond(0),bdip,bdi,mps_A_init[i].bond(2)])
  t0.PutBlock(mps_A_init[i].GetBlock())
  t1=t0*1.0
  t0.SetLabel([1,2,3,4])
  t1.SetLabel([-1,2,-3,-4])
  Resutl=uni10.Contract(t0,t1)
  Resutl=uni10.Permute(Resutl,[1,-1,3,-3,4,-4],4)
  Resutl=Resutl.CombineBond([1,-1])
  Resutl=Resutl.CombineBond([4,-4])
  Resutl=Resutl.CombineBond([3,-3])
  Resutl=uni10.Permute(Resutl,[1,3,4],2)
  A_list.append(Resutl)
 #print  A_list[0].printDiagram()
 list_bond=[]
 for q in xrange(mps_A_init.N):
   list_bond.append(A_list[q].bond(2).dim())

 mps_A_origin=MPSclass2.MPS(A_list[1].bond(1).dim(),max(list_bond),mps_A_init.N)

 for i in xrange(mps_A_init.N):
  mps_A_origin[i]=A_list[i]*1.0
  
  
  
 D=bdi.dim()
 N=mps_A_init.N
 #chi=2
 #N_iter=20
 Delta=1.0e-6
 Norm_init_A=0.50

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)


 #################### Identity ############################
 Iden=uni10.UniTensorR([bdi,bdi1,bdo,bdo1], "A_middle")
 Iden.Identity()
 Iden=uni10.Permute(Iden,[1,0,2,3],3)
 #print Iden.printDiagram
 Iden_com=Iden*1.0
 Iden_com=Iden_com.CombineBond([0,2])
 mps_Iden=MPSclass2.MPS(D*D,1,N)   #randuni, rand, ortho

 for i in xrange(N):
  mps_Iden[i]=Iden_com*1.0
 ##########################################################

 ################PEPS-tensor##############################
 #A_list=[None]*N
 #mps_A=Init_A_PEPS(A_list,N,bdi_phys,bdi,bdi1,bdo,bdo1 )
 mps_A=mps_A_origin*1.0
#  if mps_A.norm()>1:
#    print "Warning", mps_A.norm()
  
 mps_A=mps_A.normalize()
 mps_A=mps_A*(Norm_init_A)






 mps_AF=copy.copy(mps_A)
 #print mps_AF[0].PrintDiagram(),  mps_Iden[0].PrintDiagram()
 mps_A=mps_AF+mps_Iden*(Delta)
 mps_A_non=mps_AF+mps_Iden*(0.0)


 #print "hiiiii", mps_A.norm(), A_mat.norm()*A_mat.norm(),A_mat[4]


 mps_Z=copy.copy(mps_Iden)
 mps_Y=copy.copy(mps_A)


 Iter=[]
 ZY_acc_mps=[]
 ZY_acc_mat=[]
 ZAZ_acc_mps=[]
 ZAZ_acc_mat=[]
 AYY_acc_mps=[]
###################First-Step#############################
#print mps_Y[0], mps_A[0]
 E_1=0
 E_0=1
 mps_Y_best=copy.copy(mps_Y)
 for i in xrange(N_iter[1]):
  #print "\n"
  #print "Step", i

  #mps_YZY=make_YZY_mps(mps_Y,mps_Z)
  #mps_ZYZ=make_ZYZ_mps(mps_Y,mps_Z)

  mps_YZY=Appmake_YZY_mps(mps_Y,mps_Z,chi,bdi)
  mps_ZYZ=Appmake_ZYZ_mps(mps_Y,mps_Z,chi,bdi)



  mps_Ynew=((mps_Y*3.0)-mps_YZY)*(0.5)
  mps_Znew=((mps_Z*3.0)-mps_ZYZ)*(0.5)


  if mps_Ynew.D>chi:
   #print "used chi for Y=", chi 
   mps_Y=mps_Ynew.appSVD(chi)
  else:
   #print "used chi for Y=",  mps_Ynew.D 
   mps_Y=mps_Ynew.appSVD(mps_Ynew.D)

  if mps_Znew.D>chi:
   #print "used chi for Z=", chi
   mps_Z=mps_Znew.appSVD(chi)
  else:
   #print "used chi for Z=", mps_Znew.D 
   mps_Z=mps_Znew.appSVD(mps_Znew.D)

  #print "Fidel:YZ", FidelMPSYZ(mps_Z, mps_Y, mps_Iden, bdi)
  #print "Fidel:ZAZ", FidelMPSZAZ(mps_A,mps_Z, mps_Iden, chi,  bdi)
  #print "Fidel:ZA_nonZ", FidelMPSZAZ(mps_A_non,mps_Z, mps_Iden, chi, bdi)
  print "Fidel:A_nonYY", FidelMPSAYY(mps_A_non,mps_Y, mps_Iden, chi, bdi)
  E_0=copy.copy(E_1)
  E_1=FidelMPSAYY(mps_A_non, mps_Y, mps_Iden, chi, bdi)
  if E_1<E_0:
   #print "Break"
   break
  else: mps_Y_best=copy.copy(mps_Y)
  #Iter.append(i)
  #ZY_acc_mps.append(FidelMPSYZ(mps_Z, mps_Y, mps_Iden,bdi))
  #ZAZ_acc_mps.append(FidelMPSZAZ(mps_A,mps_Z, mps_Iden, chi,bdi))
  #AYY_acc_mps.append(FidelMPSAYY(mps_A_non,mps_Y, mps_Iden, chi,bdi))





 #mps_YY=YY_multy(mps_Y_best,bdi)






  #print mps_Y_best[0].printDiagram()
  #print mps_Y_best[1].printDiagram()



  #print mps_YY[0].printDiagram()
  #print mps_YY[1].printDiagram()
 #norm1=mps_YY.norm()
 #norm2=mps_A_origin.norm()

 #mps_Y_best=mps_Y_best*((norm2**(0.25))/(norm1**(0.25)))
 #mps_YY=YY_multy(mps_Y_best,bdi)

 #print  mps_Y_best.norm(),  mps_A_origin.norm(), mps_YY.norm()


 print "Fidel_square", FidelMPSAYY(mps_A_non,mps_Y_best, mps_Iden, chi, bdi)


 
#  Y_list=[None]*mps_Y.N
#  for i in xrange(mps_Y.N):
#   A=copy.copy(mps_Y[i])
#   A1=uni10.UniTensorR([A.bond(0),bdi,bdi,A.bond(2)])
#   A1.PutBlock(A.GetBlock())
#   A1.SetLabel([1,2,-2,3])
# 
#   A1.uni10.Permute([1,-2,2,3],3)
#   A1.CombineBond([-2,2])
#   A1.uni10.Permute([1,-2,3],2)
#   Y_list[i]=A1
# 
#  list_bond=[]
#  for q in xrange(N):
#    list_bond.append(Y_list[q].bond(2).dim())
# 
# 
#  mps_Y_best=MPSclass.MPS(Y_list[1].bond(1).dim(),max(list_bond),N,'ortho')



 return mps_Y_best
 


# Gamma_a.save("Store/Gamma_a")
# Gamma_b.save("Store/Gamma_b")
# Gamma_a=uni10.UniTensorR("Store/Gamma_a")
# Gamma_b=uni10.UniTensorR("Store/Gamma_b")

  #ZAZ_acc_mat.append(Fidel(Z_mat*A_mat*Z_mat, Iden_mat))

  #print Y_mat*Z_mat
  #test(mps_Z,mps_Y)
  #test1(mps_Z,mps_Y)
##############################################################

#file = open("Acc.txt", "w")
#for index in range(len(ZY_acc_mps)):
#    file.write(str(Iter[index]) + " " + str(ZY_acc_mps[index])+" " + str(ZAZ_acc_mps[index])+" "+ "\n")
#file.close()



#file = open("Acc.txt", "w")
#for index in range(len(ZY_acc_mps)):
#    file.write(str(Iter[index]) + " " + str(ZY_acc_mps[index])+" " + str(ZY_acc_mat[index])+" " + str(ZAZ_acc_mps[index])+" " + str(ZAZ_acc_mat[index]) + "\n")
#file.close()


############################Inverse###################################

#M_mat=copy.copy(A_mat)
# V_mat=copy.copy(Z_mat)
# #V_mat=copy.copy(Iden_mat)
# 
# mps_V=copy.copy(mps_Z)
# 
# print "Vnorm_0", V_mat.norm()*V_mat.norm(),mps_V.norm()
# 
# 
# for i in xrange(N_iterInv):
# 
#   print "\n"
#   print "Step", i
# ###############################
# #  mps_YV=Appmake_YV_mps(mps_Y,mps_V,chi)
# #  mps_t=((mps_Iden*2.0)-mps_YV)
# #  mps_YV=Appmake_YV_mps(mps_V,mps_t,chi)
# #  if mps_YV.D>chi:
# #   print "used chi for V=", chi 
# #   mps_V=mps_YV.appSVD(chi)
# #  else:
# #   print "used chi for V=",  mps_YV.D 
# #   mps_V=mps_YV.appSVD(mps_YV.D)
# 
# ###############3thrd#############################
# #   mps_YV=Appmake_YV_mps(mps_Y,mps_V,chi)
# #   mps_t=((mps_Iden*3.0)-mps_YV)
# #   if mps_t.D>chi:
# #    mps_t=mps_t.appSVD(chi)
# #   else:
# #    mps_t=mps_t.appSVD(mps_t.D)
# #   mps_t=Appmake_YV_mps(mps_YV,mps_t,chi)
# #   mps_t=((mps_Iden*3.0)-mps_t)
# #   if mps_t.D>chi:
# #    mps_t=mps_t.appSVD(chi)
# #   else:
# #    mps_t=mps_t.appSVD(mps_t.D)
# #   mps_V=Appmake_YV_mps(mps_V,mps_t,chi)
# ###################################################
# 
# 
# ###############6thrd#############################
# 
#   mps_YV=Appmake_YV_mps(mps_Y,mps_V,chi)
#   mps_t=((mps_Iden*1.0)-mps_YV)
#   if mps_t.D>chi:
#    mps_t=mps_t.appSVD(chi)
#   else:
#    mps_t=mps_t.appSVD(mps_t.D)
#   mps_t=Appmake_YV_mps(mps_YV,mps_t,chi)
#   mps_t=((mps_Iden*1.0)-mps_t)
#   if mps_t.D>chi:
#    mps_t=mps_t.appSVD(chi)
#   else:
#    mps_t=mps_t.appSVD(mps_t.D)
#   mps_t1=copy.copy(mps_t)
# #########
#   mps_t=((mps_Iden*3.0)-mps_YV)
#   if mps_t.D>chi:
#    mps_t=mps_t.appSVD(chi)
#   else:
#    mps_t=mps_t.appSVD(mps_t.D)
#   mps_t=Appmake_YV_mps(mps_YV,mps_t,chi)
#   mps_t=((mps_Iden*3.0)-mps_t)
#   if mps_t.D>chi:
#    mps_t=mps_t.appSVD(chi)
#   else:
#    mps_t=mps_t.appSVD(mps_t.D)
#   mps_t2=copy.copy(mps_t)
# ########
#   mps_t=((mps_Iden*2.0)-mps_YV)
#   if mps_t.D>chi:
#    mps_t=mps_t.appSVD(chi)
#   else:
#    mps_t=mps_t.appSVD(mps_t.D)
#   mps_t3=copy.copy(mps_t)
#   mps_t=Appmake_YV_mps(mps_V,mps_t3,chi)
#   mps_t=Appmake_YV_mps(mps_t,mps_t2,chi)
#   mps_V=Appmake_YV_mps(mps_t,mps_t1,chi)
# ###################################################
# 
#   #V_mat=V_mat*(Iden_mat*2.0+(-1.0)*Y_mat*V_mat)
#   #V_mat=V_mat*(Iden_mat*3.0+(-1.0)*(Y_mat*V_mat)*(3*Iden_mat+(-1.0)*Y_mat*V_mat))
#   #V_mat=V_mat*(Iden_mat*2.0+(-1.0)*(Y_mat*V_mat))*(3*Iden_mat+(-1.0)*Y_mat*V_mat*(3*Iden_mat+(-1.0)*Y_mat*V_mat))*(Iden_mat+(-1.0)*Y_mat*V_mat*(Iden_mat+(-1.0)*Y_mat*V_mat))
# 
#   #print "Vnorm", V_mat.norm()*V_mat.norm(),mps_V.norm()
#   print "Vnorm", mps_V.norm()
# 
#   #print "Fidel:AV", Fidel(Y_mat*V_mat,Iden_mat )
#   #print "Fidel:AV", Fidel(Y_mat_non*V_mat,Iden_mat )
#   print "Fidel:YZ", FidelMPSYV(mps_Y,mps_V, mps_Iden, chi)
#   print "Fidel:ZA_nonZ", FidelMPSZAZ(mps_A_non,mps_V, mps_Iden, chi)
#   print "Fidel:ZAZ", FidelMPSZAZ(mps_A,mps_V, mps_Iden, chi)
# 
# 
# 
# 
# 
