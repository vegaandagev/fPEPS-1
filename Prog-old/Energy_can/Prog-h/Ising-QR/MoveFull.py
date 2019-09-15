import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
#import random
import copy
#import line_profiler
import TruncateU
import basicB
import basicA
import basicC
#import time


def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 #sv_mat = uni10.Matrix( len(max_list), len(max_list), max_list, True)
 #return sv_mat.absMax()
 max_list_f=[abs(x) for x in max_list]
 #print max_list_f, max(max_list_f)
 return max(max_list_f)

def norm_CTM(c):
 Max_val=abs(MaxAbs(c))
 if (( Max_val < 0.50e-1) or (Max_val > 0.50e+1)) and abs(Max_val)>1.0e-12:
  #print  "Max_val", Max_val
  c*=(1.00/Max_val)
 return c


def distance(theta,A):
   blk_qnums = theta.blockQnum()
   val=0
   for qnum in blk_qnums:
    T1=theta.getBlock(qnum)
    T2=A.getBlock(qnum)
    print "col", theta.getBlock(qnum).col(), qnum 
    for  i  in xrange( int(theta.getBlock(qnum).col()) ):
     if abs(T1[i]) > 1.00e-11:  
      val=val+abs((T1[i]-T2[i]) / T1[i])
      #if abs((T1[i]-T2[i]) / T1[i]) > 1.00e+1: print "hi", T1[i], T2[i], i 
     else: val=val+(T1[i]-T2[i]); #print T1[i]-T2[i]
   return val 


def check_eigenvalues(s):
 M_s=s.getBlock()
 p=0
 for i in xrange(int(M_s.row())):
    for j in xrange(int(M_s.col())):
     if i==j:
      if M_s[i*int(M_s.col())+j]<1.e-8:
       p=i
#       print p, M_s[i*int(M_s.col())+j]
       break
    else:
        continue  # executed if the loop ended normally (no break)
    break  # executed if 'continue' was skipped (break)

 return p


def pick_vec(p,U):
 x=int(U.row())
 y=int(U.col())
 Vec_up=uni10.Matrix(x,1)
 for i in xrange(x):
      Vec_up[i]=U[i*y+p]

 return Vec_up


def add_vec_to_mat(U,vec):
 x=int(U.row())
 y=int(U.col())
 U_new=copy.copy(U)
 U_new.resize(x, y+1)
 p=y
 for i in xrange(x):
   U_new[i*(y+1)+p]=vec[i*(1)+0]
 return U_new


def Add_onevector(Rb_mat, R_mat, p, V, U):
 U_mat=U.getBlock()
 V_mat=V.getBlock()
 U_mat.transpose()

 Vec_up=pick_vec(p,U_mat)
 Vec_vp=pick_vec(p,V_mat)

 Vec_up.randomize()
 Vec_vp.randomize()

 Vec_up=R_mat*Vec_up
 Rb_mat.transpose()
 Vec_vp=Rb_mat*Vec_vp
 R_mat=add_vec_to_mat(R_mat,Vec_vp)
 Rb_mat=add_vec_to_mat(Rb_mat,Vec_up)
 Rb_mat.transpose()

 A=Rb_mat*R_mat   
 return A, Rb_mat,  R_mat

#@profile 
def produce_projectives(theta,theta1,chi_dim):
 theta=copy.copy(theta)
 theta1=copy.copy(theta1)
 
 theta.setLabel([1,2,20,3,4,40])
 theta.permute([1,2,20,3,4,40],3) 

 
 U, s, V=TruncateU.svd_parity1(theta)
 
 if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
  s=s*(1.00/MaxAbs(s))

 U.setLabel([1,2,20,-1,-2,-3])
 s.setLabel([-1,-2,-3,3,4,5])
 R=U*s
 R.permute([1,2,20,3,4,5],3)


 theta1.setLabel([1,2,20,3,4,40])
 theta1.permute([1,2,20,3,4,40],3) 

 
 U, s, V=TruncateU.svd_parity1(theta1)

 if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
  s=s*(1.00/MaxAbs(s))

 U.setLabel([1,2,20,-1,-2,-3])
 s.setLabel([-1,-2,-3,6,7,8])
 Rb=U*s
 Rb.permute([1,2,20,6,7,8],3)
 Rb.permute([6,7,8,1,2,20],3)
 
 
 A=R*Rb
 A.permute([6,7,8,3,4,5],3)


 V, U, s=TruncateU.setTruncation(A, chi_dim) 

 if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
  s=s*(1.00/MaxAbs(s))


 U.setLabel([-1,3,4,5])
 V.setLabel([6,7,8,-1])

# if MaxAbs(s) > 1.0e-8:
#  s=s*(1.00/MaxAbs(s)) 

 s=TruncateU.inverse(s)
 s=TruncateU.Sqrt(s)
 
 s.setLabel([6,-1])
 U=s*U
 U.permute([6,3,4,5],1)

 s.setLabel([-1,9])
 V=V*s
 V.permute([6,7,8,9],3) 


 R.permute([1,2,20,3,4,5],3)
 U.transpose()
 U1x=R*U
 U1x.permute([1,2,20,6],3)

 Rb.permute([6,7,8,1,2,20],3)
 V.transpose()
 U1x_trans=Rb*V
 U1x_trans.permute([9,1,2,20],1)

 return U1x, U1x_trans 









#def produce_projectives(theta,theta1,chi_dim):
# theta=copy.copy(theta)
# theta1=copy.copy(theta1)
# 
# theta.setLabel([1,2,20,3,4,40])
# theta.permute([1,2,20,3,4,40],3) 

# 
# U, s, V=TruncateU.svd_parity1(theta)


# if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
#  s=s*(1.00/MaxAbs(s))

#  
# U.setLabel([1,2,20,-1,-2,-3])
# s.setLabel([-1,-2,-3,3,4,5])
# R=U*s
# R.permute([1,2,20,3,4,5],3)


#################################################
## R, q=TruncateU.lq_parity1(theta)
## R.setLabel([1,2,20,3,4,5])
################################################


# theta1.setLabel([1,2,20,3,4,40])
# theta1.permute([1,2,20,3,4,40],3) 

# 
# U, s, V=TruncateU.svd_parity1(theta1)

# if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
#  s=s*(1.00/MaxAbs(s))


# U.setLabel([1,2,20,-1,-2,-3])
# s.setLabel([-1,-2,-3,6,7,8])
# Rb=U*s
# Rb.permute([1,2,20,6,7,8],3)
# Rb.permute([6,7,8,1,2,20],3)
# 
#################################################
## Rb, q=TruncateU.lq_parity1(theta1)
## Rb.setLabel([1,2,20,6,7,8])
## Rb.permute([6,7,8,1,2,20],3)
################################################
# 
# A=Rb*R
# A.permute([6,7,8,3,4,5],3)

# V, U, s=TruncateU.setTruncation(A, chi_dim) 

## M_s=s.getBlock()
## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print "1", M_s[i*int(M_s.col())+j]


# U.setLabel([-1,3,4,5])
# V.setLabel([6,7,8,-1])

# if MaxAbs(s) > 1.0e+5 or MaxAbs(s) < 1.0e-2:
#  #print "MaxAbs(s)1", MaxAbs(s)
#  s=s*(1.00/MaxAbs(s))
#   
# s=TruncateU.inverse(s)
# s=TruncateU.Sqrt(s)
# 
# s.setLabel([6,-1])
# U=s*U
# U.permute([6,3,4,5],1)

# s.setLabel([-1,9])
# V=V*s
# V.permute([6,7,8,9],3) 


# R.permute([1,2,20,3,4,5],3)
# U.transpose()
# U1x=R*U
# U1x.permute([1,2,20,6],3)



# Rb.permute([6,7,8,1,2,20],3)
# V.transpose()
# U1x_trans=V*Rb
# U1x_trans.permute([9,1,2,20],1)
################################################################
# return U1x, U1x_trans 



#def produce_projectives1(theta,theta1,chi_dim):
# theta=copy.copy(theta)
# theta1=copy.copy(theta1)
# 
# theta.setLabel([1,2,20,3,4,40])
# theta.permute([1,2,20,3,4,40],3) 

# 
# U, s, V=TruncateU.svd_parity1(theta)
## M_s=s.getBlock()
# #print M_s

# if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
#  #print "MaxAbs(s)1", MaxAbs(s)
#  s=s*(1.00/MaxAbs(s))

## M_s=s.getBlock()

## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print "1", M_s[i*int(M_s.col())+j]
#  
# U.setLabel([1,2,20,-1,-2,-3])
# s.setLabel([-1,-2,-3,3,4,5])
# R=U*s
# R.permute([1,2,20,3,4,5],3)



# theta1.setLabel([1,2,20,3,4,40])
# theta1.permute([1,2,20,3,4,40],3) 

# 
# U, s, V=TruncateU.svd_parity1(theta1)

# if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
#  #print "MaxAbs(s)1", MaxAbs(s)
#  s=s*(1.00/MaxAbs(s))


# U.setLabel([1,2,20,-1,-2,-3])
# s.setLabel([-1,-2,-3,6,7,8])
# Rb=U*s
# Rb.permute([1,2,20,6,7,8],3)
# Rb.permute([6,7,8,1,2,20],3)
# 
# 
#  
# A=R*Rb
# A.permute([6,7,8,3,4,5],3)

## print A
# V, U, s=TruncateU.setTruncation(A, chi_dim) 


## M_s=s.getBlock()

## for i in xrange(int(M_s.row())):
##  for j in xrange(int(M_s.col())):
##   if i==j:
##    print "1", M_s[i*int(M_s.col())+j]

# U.setLabel([-1,3,4,5])
# V.setLabel([6,7,8,-1])

# if MaxAbs(s) > 1.0e+7 or MaxAbs(s) < 1.0e-1:
#  #print "MaxAbs(s)1", MaxAbs(s)
#  s=s*(1.00/MaxAbs(s))
#   
# s=TruncateU.inverse(s)
# s=TruncateU.Sqrt(s)
# 
# s.setLabel([6,-1])
# U=s*U
# U.permute([6,3,4,5],1)

# s.setLabel([-1,9])
# V=V*s
# V.permute([6,7,8,9],3) 


# R.permute([1,2,20,3,4,5],3)
# U.transpose()
# U1x=R*U
# U1x.permute([1,2,20,6],3)



# Rb.permute([6,7,8,1,2,20],3)
# V.transpose()
# U1x_trans=Rb*V
# U1x_trans.permute([9,1,2,20],1)

# return U1x, U1x_trans 




#@profile 
def  add_left1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,chi,D):

 
 chi_dim=0
 for i in xrange(len(chi)):
  chi_dim+=chi[i]

# t0=time.time()

 CTM_1 = uni10.Network("Network/CTM1.net")
 CTM_1.putTensor('c1',c1)
 CTM_1.putTensor('c2',c2)
 CTM_1.putTensor('Ta1',Ta1)
 CTM_1.putTensor('Ta2',Ta2)
 CTM_1.putTensor('Tb1',Tb1)
 CTM_1.putTensor('Tb4',Tb4)
 CTM_1.putTensor('a',a)
 CTM_1.putTensor('b',b)
 theta=CTM_1.launch()
 theta.permute([100, 300, -300 , 400, 200 ,-200],3)

# print time.time() - t0, "up"

 
 CTM_2 = uni10.Network("Network/CTM2.net")
 CTM_2.putTensor('c3',c3)
 CTM_2.putTensor('c4',c4)
 CTM_2.putTensor('Ta3',Ta3)
 CTM_2.putTensor('Ta4',Ta4)
 CTM_2.putTensor('Tb2',Tb2)
 CTM_2.putTensor('Tb3',Tb3)
 CTM_2.putTensor('c',c)
 CTM_2.putTensor('d',d)
 theta1=CTM_2.launch()
 theta1.permute([100, 300, -300 , 400, 200 ,-200],3)


 
 U1x, U1x_trans=produce_projectives(theta,theta1, chi_dim)
 #print "norm1", U1x.norm() 

 
 theta.permute([  400, 200 ,-200, 100, 300, -300], 3)
 theta1.permute([  400, 200 ,-200, 100, 300, -300], 3)

 U2x, U2x_trans=produce_projectives(theta,theta1, chi_dim)
 #print "norm2", U1x.norm() 

 
# t0=time.time()

 Ta4p=copy.copy(Ta4)
 Ta4p.setName("Ta4p")
 Tb4p=copy.copy(Tb4)
 Ta2p=copy.copy(Ta2)
 Tb2p=copy.copy(Tb2)
 U1xb=copy.copy(U1x)
 U1xb_trans=copy.copy(U1x_trans)
 U2xb=copy.copy(U2x)
 U2xb_trans=copy.copy(U2x_trans)
 ap=copy.copy(a)
 bp=copy.copy(b)
 cp=copy.copy(c)
 dp=copy.copy(d)

 CTM_3 = uni10.Network("Network/CTM3.net")
 CTM_3.putTensor('c1',c1)
 CTM_3.putTensor('c2',c2)
 CTM_3.putTensor('Ta1',Ta1)
 CTM_3.putTensor('Ta2',Ta2)
 CTM_3.putTensor('Ta4',Ta4)
 CTM_3.putTensor('Tb1',Tb1)
 CTM_3.putTensor('Tb2',Tb2)
 CTM_3.putTensor('Tb4',Tb4)
 CTM_3.putTensor('a',a)
 CTM_3.putTensor('b',b)
 CTM_3.putTensor('c',c)
 CTM_3.putTensor('d',d)
 CTM_3.putTensor('U1x',U1x)
 CTM_3.putTensor('U1x_trans',U1x_trans)
 CTM_3.putTensor('U2x',U2x)
 CTM_3.putTensor('U2x_trans',U2x_trans)
 theta=CTM_3.launch() 
 theta.permute([100, 300, -300 , 400, 200 ,-200],3)
# print time.time() - t0, "upup"

 

 
 CTM_4 = uni10.Network("Network/CTM4.net")
 CTM_4.putTensor('c3',c3)
 CTM_4.putTensor('c4',c4)
 CTM_4.putTensor('Ta2p',Ta2p)
 CTM_4.putTensor('Ta3',Ta3)
 CTM_4.putTensor('Ta4p',Ta4p)
 CTM_4.putTensor('Tb2p',Tb2p)
 CTM_4.putTensor('Tb3',Tb3)
 CTM_4.putTensor('Tb4p',Tb4p)
 CTM_4.putTensor('ap',ap)
 CTM_4.putTensor('bp',bp)
 CTM_4.putTensor('cp',cp)
 CTM_4.putTensor('dp',dp)
 CTM_4.putTensor('U1xb',U1xb)
 CTM_4.putTensor('U1xb_trans',U1xb_trans)
 CTM_4.putTensor('U2xb',U2xb)
 CTM_4.putTensor('U2xb_trans',U2xb_trans)
 theta1=CTM_4.launch() 
 theta1.permute([100, 300, -300 , 400, 200 ,-200],3)
 

 U3x, U3x_trans=produce_projectives(theta,theta1, chi_dim)
 #print "norm3", U3x.norm() 
 theta.permute([  400, 200 ,-200, 100, 300, -300], 3)
 theta1.permute([  400, 200 ,-200, 100, 300, -300], 3)

 U4x, U4x_trans=produce_projectives(theta,theta1, chi_dim)
 #print "norm4", U3x.norm() 


###############################################################################
 c1.setLabel([0,1])
 Tb1.setLabel([1,2,-2,3])
 c1bar=c1*Tb1
 c1bar.permute([0,2,-2,3],3)

 c4.setLabel([0,1])
 Ta3.setLabel([1,2,-2,3])
 c4bar=c4*Ta3
 c4bar.permute([0,2,-2,3],1)

 ##############################
 U3x_trans.setLabel([4,0,2,-2])
 c1bar=c1bar*U3x_trans
 c1bar.permute([4,3],1)
 ############################
 U3x.setLabel([0,2,-2,4])
 c4bar=c4bar*U3x
 c4bar.permute([4,3],0)

 
 #############################
 Tb4.setLabel([0,1,-1,2])
 a.setLabel([1,-1,3,-3,4,-4,5,-5])
 U3x.setLabel([2,5,-5,6])
 U1x_trans.setLabel([7,0,3,-3])
 Tb4bar=((Tb4*U3x)*a)*U1x_trans
 Tb4bar.permute([7,4,-4,6],1)
# absorb_1 = uni10.Network("Network/absorb1.net")
# absorb_1.putTensor('Tb4',Tb4)
# absorb_1.putTensor('a',a)
# absorb_1.putTensor('U3x',U3x)
# absorb_1.putTensor('U1x_trans',U1x_trans)
# Tb4bar=absorb_1.launch() 

 ###########################
 Ta4.setLabel([0,1,-1,2])
 c.setLabel([1,-1,3,-3,4,-4,5,-5])
 U1x.setLabel([2,5,-5,6])
 U3x_trans.setLabel([7,0,3,-3])
 Ta4bar=((Ta4*U3x_trans)*c)*U1x
 Ta4bar.permute([7,4,-4,6],1)
#############################
###############################

 #################3################################
 c2.setLabel([0,1])
 Ta1.setLabel([3,2,-2,0])
 c2bar=c2*Ta1
 c2bar.permute([1,2,-2,3],3)
# c3.setLabel([0,1])
 c3.setLabel([1,0])
 Tb3.setLabel([3,2,-2,1])
 c3bar=c3*Tb3
 c3bar.permute([3,0,2,-2],1)
 ##############################
 U4x_trans.setLabel([4,1,2,-2])
 c2bar=c2bar*U4x_trans
 c2bar.permute([3,4],2)
 ############################
 U4x.setLabel([0,2,-2,4])
 c3bar=c3bar*U4x
# c3bar.permute([4,3],1)
 c3bar.permute([3,4],1)
 #############################
 Ta2.setLabel([0,1,-1,2])
 b.setLabel([4,-4,5,-5,1,-1,3,-3])
 U4x.setLabel([2,3,-3,7])
 U2x_trans.setLabel([6,0,5,-5])
 Ta2bar=((Ta2*U4x)*b)*U2x_trans
 Ta2bar.permute([6,4,-4,7],3)
 ###########################
 Tb2.setLabel([0,1,-1,2])
 d.setLabel([4,-4,5,-5,1,-1,3,-3])
 U2x.setLabel([2,3,-3,7])
 U4x_trans.setLabel([6,0,5,-5])
 Tb2bar=((Tb2*U4x_trans)*d)*U2x
 Tb2bar.permute([6,4,-4,7],3)
 ###########################
 ###########################

 Tb2bar, Ta2bar=equall_dis_qr(Tb2bar, Ta2bar )
 Ta4bar, Tb4bar=equall_dis_qr(Ta4bar, Tb4bar )
 Ta4bar.permute([62,58,57,17],1)
 Tb4bar.permute([17,64,58,57],1)
 

 c3=norm_CTM(c3bar)
 c2=norm_CTM(c2bar)
 Ta2=norm_CTM(Ta2bar)
 Tb2=norm_CTM(Tb2bar)

 c1=norm_CTM(c1bar)
 c4=norm_CTM(c4bar)
 Ta4=norm_CTM(Ta4bar)
 Tb4=norm_CTM(Tb4bar)


 return c1, Ta4, Tb4, c4, c2, Ta2, Tb2, c3
 
 
 
def equall_dis_qr(plist0, plist1 ):
 
 plist0.setLabel([62,58,57,17])
 plist1.setLabel([17,64,58,57])
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64,58,57],1)
 
 q, r=qr_parity2(plist0)
 q.setLabel([62,58,57,1])
 r.setLabel([1,0])
 

 l, qq=lq_parity2(plist1)
 qq.setLabel([-1,64,58,57])
 l.setLabel([0,-1])

 teta=l*r
 teta.permute([1,-1],1)

 U, s, V =svd_parity2(teta)
 s=TruncateU.Sqrt(s)

 s.setLabel([0,17])
 U.setLabel([1,0])
 U=U*s

 s.setLabel([17,0])
 V.setLabel([0,-1])
 V=s*V

 plist0=q*U
 plist1=V*qq
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64,58,57],3)
 return plist0, plist1
 
def equall_dis_qr1(plist0, plist1 ):
 
 plist0.setLabel([62,58,57,17])
 plist1.setLabel([64,17])
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64],1)
 teta=plist0*plist1
 teta.permute([62,58,57,64],3)

 chi_dim=(teta.bond(3).dim())

 U, s, V =TruncateU.svd_parity(teta,chi_dim)
 s=TruncateU.Sqrt(s)

 s.setLabel([0,17])
 U.setLabel([62,58,57,0])
 plist0=U*s

 s.setLabel([17,0])
 V.setLabel([0,64])
 plist1=s*V
 
 plist0.permute([62,58,57,17],3)
 plist1.permute([17,64],1)
 return plist0, plist1
 
 
 
def equall_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):
  
 Norm=magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 if Norm[0] < 0: c1=-1.0*c1;

 while Norm[0]<1.0e-1: 
  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d=checking_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  Norm=magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  #print Norm[0]

 while Norm[0]>1.0e+5: 
  c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d=checking_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  Norm=magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
  #print Norm[0]

 return c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4
 
def checking_norm(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):
 Norm=magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d)
 if Norm[0] < 0: c1=-1.0*c1;
 if Norm[0] < 1.0e-1:
  c1*=1.4
  c2*=1.4
  c3*=1.4
  c4*=1.4
  Ta1*=1.4
  Ta2*=1.4
  Ta3*=1.4
  Ta4*=1.4
  Tb1*=1.4
  Tb2*=1.4
  Tb3*=1.4
  Tb4*=1.4

 elif Norm[0] > 1.e+5:
  c1*=0.92
  c2*=0.92
  c3*=0.92
  c4*=0.92
  Ta1*=0.92
  Ta2*=0.92
  Ta3*=0.92
  Ta4*=0.92
  Tb1*=0.92
  Tb2*=0.92
  Tb3*=0.92
  Tb4*=0.92
 return c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d
 
 
def  magnetization_value(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d):

 CTM_net = uni10.Network("Network/CTM.net")
 CTM_net.putTensor('c1',c1)
 CTM_net.putTensor('c2',c2)
 CTM_net.putTensor('c3',c3)
 CTM_net.putTensor('c4',c4)
 CTM_net.putTensor('Ta1',Ta1)
 CTM_net.putTensor('Ta2',Ta2)
 CTM_net.putTensor('Ta3',Ta3)
 CTM_net.putTensor('Ta4',Ta4)
 CTM_net.putTensor('Tb1',Tb1)
 CTM_net.putTensor('Tb2',Tb2)
 CTM_net.putTensor('Tb3',Tb3)
 CTM_net.putTensor('Tb4',Tb4)
 CTM_net.putTensor('a',a)
 CTM_net.putTensor('b',b)
 CTM_net.putTensor('c',c)
 CTM_net.putTensor('d',d)
 norm=CTM_net.launch()
 
# c1.setLabel([4,1])
# c2.setLabel([3,7])
# c3.setLabel([24,4])
# c4.setLabel([18,22])
# Ta1.setLabel([2,6,-6,3]) 
# Ta2.setLabel([14,10,-10,7]) 
# Ta3.setLabel([22,19,-19,23]) 
# Ta4.setLabel([18,15,-15,11]) 
# Tb1.setLabel([1,5,-5,2]) 
# Tb2.setLabel([4,17,-17,14]) 
# Tb3.setLabel([23,20,-20,24]) 
# Tb4.setLabel([11,8,-8,4])
# a.setLabel([8,-8,12,-12,9,-9,5,-5])
# b.setLabel([9,-9,13,-13,10,-10,6,-6])
# c.setLabel([15,-15,19,-19,16,-16,12,-12])
# d.setLabel([16,-16,20,-20,17,-17,13,-13])
# norm=(((((c3*Tb2)*Tb3)*(d))*(((c4*Ta3)*Ta4)*c))*(((c2*Ta2)*Ta1)*b))*(((c1*Tb1)*Tb4)*a) 
#print 'hi1', '\n',norm,
 return norm


def  Env_energy_h(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1,c2,c3,c4,Ta1,Tb1,Ta2,Tb2,Ta3,Tb3,Ta4, Tb4,D,d_phys)
 
 E_ab=basicB.Energy_ab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0, a_u, b_u)

 return E_ab

def  Env_energy_v(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,c_u,a_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicB.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ca=basicB.Energy_ca(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,c_u,a_u)

 return E_ca

def  Env_energy_D(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,c_u,b_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicA.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ca=basicA.energy_cb(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,c_u,b_u)

 return E_ca

def  Env_energy_D1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,d_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8=basicA.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys)

 E_ca=basicA.energy_ad(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d, H0,a_u,d_u)

 return E_ca

def Env_energy_three(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)

 E_cab=basicC.energy_cab(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H0)
 return E_cab


def Env_energy_three1(c1,c2,c3,c4,Ta1,Ta2,Ta3,Ta4,Tb1,Tb2,Tb3,Tb4,a,b,c,d,a_u,b_u,c_u,d_u,H0,D,d_phys):

 E1, E2, E3, E4, E5, E6, E7, E8, a_u, b_u, c_u, d_u,a, b, c, d=basicC.produce_Env(a,b,c,d,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4,D,d_phys,a_u, b_u,c_u,d_u)

 E_abd=basicC.energy_abd(E1, E2, E3, E4, E5, E6, E7, E8, a, b, c,d,a_u,b_u,c_u,d_u, H0)
 return E_abd





 
def permuteN(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 a.permute([3,-3,2,-2,1,-1,0,10],4)
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 ##print'a', a
 b.setLabel([0,10,1,-1,2,-2,3,-3])
 b.permute([3,-3,2,-2,1,-1,0,10],4)
 b.setLabel([0,10,1,-1,2,-2,3,-3])

 c.setLabel([0,10,1,-1,2,-2,3,-3])
 c.permute([3,-3,2,-2,1,-1,0,10],4)
 c.setLabel([0,10,1,-1,2,-2,3,-3])

 d.setLabel([0,10,1,-1,2,-2,3,-3])
 d.permute([3,-3,2,-2,1,-1,0,10],4)
 d.setLabel([0,10,1,-1,2,-2,3,-3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([0,1],0)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([0,1],2)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,-1,2])
 Ta1.permute([2,1,-1,0],1)
 Ta1.setLabel([0,1,-1,2])
 
 Ta2.setLabel([0,1,-1,2])
 Ta2.permute([2,1,-1,0],1)
 Ta2.setLabel([0,1,-1,2])
 
 Ta3.setLabel([0,1,-1,2])
 Ta3.permute([2,1,-1,0],3)
 Ta3.setLabel([0,1,-1,2])

 Ta4.setLabel([0,1,-1,2])
 Ta4.permute([2,1,-1,0],3)
 Ta4.setLabel([0,1,-1,2])

 Tb1.setLabel([0,1,-1,2])
 Tb1.permute([2,1,-1,0],1)
 Tb1.setLabel([0,1,-1,2])
 
 Tb2.setLabel([0,1,-1,2])
 Tb2.permute([2,1,-1,0],1)
 Tb2.setLabel([0,1,-1,2])
 
 Tb3.setLabel([0,1,-1,2])
 Tb3.permute([2,1,-1,0],3)
 Tb3.setLabel([0,1,-1,2])

 Tb4.setLabel([0,1,-1,2])
 Tb4.permute([2,1,-1,0],3)
 Tb4.setLabel([0,1,-1,2])

 
def permuteN1(a, b,c,d ,c1, c2,c3,c4,Ta1, Tb1,Ta2, Tb2,Ta3, Tb3,Ta4, Tb4):
 
 ##print'a', a
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 a.permute([3,-3,2,-2,1,-1,0,10],4)
 a.setLabel([0,10,1,-1,2,-2,3,-3])
 ##print'a', a
 b.setLabel([0,10,1,-1,2,-2,3,-3])
 b.permute([3,-3,2,-2,1,-1,0,10],4)
 b.setLabel([0,10,1,-1,2,-2,3,-3])

 c.setLabel([0,10,1,-1,2,-2,3,-3])
 c.permute([3,-3,2,-2,1,-1,0,10],4)
 c.setLabel([0,10,1,-1,2,-2,3,-3])

 d.setLabel([0,10,1,-1,2,-2,3,-3])
 d.permute([3,-3,2,-2,1,-1,0,10],4)
 d.setLabel([0,10,1,-1,2,-2,3,-3])

 
 c1.setLabel([0,1])
 c1.permute([1,0],1)
 c1.setLabel([0,1])

 c2.setLabel([0,1])
 c2.permute([0,1],2)
 c2.setLabel([0,1])
 
 c3.setLabel([0,1])
 c3.permute([1,0],1)
 c3.setLabel([0,1])

 c4.setLabel([0,1])
 c4.permute([0,1],0)
 c4.setLabel([0,1])


 Ta1.setLabel([0,1,-1,2])
 Ta1.permute([2,1,-1,0],3)
 Ta1.setLabel([0,1,-1,2])
 
 Ta2.setLabel([0,1,-1,2])
 Ta2.permute([2,1,-1,0],3)
 Ta2.setLabel([0,1,-1,2])
 
 Ta3.setLabel([0,1,-1,2])
 Ta3.permute([2,1,-1,0],1)
 Ta3.setLabel([0,1,-1,2])

 Ta4.setLabel([0,1,-1,2])
 Ta4.permute([2,1,-1,0],1)
 Ta4.setLabel([0,1,-1,2])

 Tb1.setLabel([0,1,-1,2])
 Tb1.permute([2,1,-1,0],3)
 Tb1.setLabel([0,1,-1,2])
 
 Tb2.setLabel([0,1,-1,2])
 Tb2.permute([2,1,-1,0],3)
 Tb2.setLabel([0,1,-1,2])
 
 Tb3.setLabel([0,1,-1,2])
 Tb3.permute([2,1,-1,0],1)
 Tb3.setLabel([0,1,-1,2])

 Tb4.setLabel([0,1,-1,2])
 Tb4.permute([2,1,-1,0],1)
 Tb4.setLabel([0,1,-1,2])
 



 
 
def qr_parity2(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3)])
    LA=uni10.UniTensor([bd1, theta.bond(3)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).qr()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])

    return GA, LA
 
def lq_parity2(theta):
    bd1=uni10.Bond(uni10.BD_OUT,theta.bond(0).Qlist())

    
    LA=uni10.UniTensor([theta.bond(0),bd1])
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3)])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).lq()
        GA.putBlock(qnum, svds[qnum][1])
        LA.putBlock(qnum, svds[qnum][0])

    return  LA, GA
 
def svd_parity2(theta):

    LA=uni10.UniTensor([theta.bond(0), theta.bond(1)])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1)])
    GB=uni10.UniTensor([theta.bond(0), theta.bond(1)])
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
 
 
 
