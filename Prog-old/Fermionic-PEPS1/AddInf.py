def fermionicOPT(bdi, bdi1):

 bdo=copy.copy(bdi)
 bdo1=copy.copy(bdi1)


 bdo.change(uni10.BD_OUT)
 bdo1.change(uni10.BD_OUT)

 T = uni10.UniTensor([bdi,bdi1,bdo,bdo1])

 # identity

 template = np.zeros([bdi.dim(), bdi1.dim(), bdo.dim(), bdo1.dim()])
 for idx in np.ndindex(template.shape):
  if idx[0]==idx[2] and idx[1]==idx[3] : template[idx] = 1

 T.setRawElem(template.reshape(-1))

 Tn = get_ndarray(T)

 for idx in np.ndindex(Tn.shape):
  if template[idx]==1 and Tn[idx] == 0:
   print "Identity: bad at {}".format(idx)
   break
  elif Tn[idx] == 1 and template[idx] == 0:
   print "Identity: bad at {}".format(idx)
   break
 #print "Identity: they are identical if nothing prints before this"


 # swap gate

 for idx in np.ndindex(template.shape):
     if idx[0]==idx[2] and idx[1]==idx[3]:
         if bdi.Qlist()[idx[0]] == bdi1.Qlist()[idx[1]] and bdi.Qlist()[idx[0]].prt() == uni10.PRT_ODD:
             template[idx] = -1
             #print "ODD",  uni10.PRT_EVEN
         elif bdo.Qlist()[idx[2]] == bdo1.Qlist()[idx[3]] and bdo.Qlist()[idx[2]].prt() == uni10.PRT_ODD:
             raise ValueError("""Qnums of incoming bonds don't match but Qnums of
                                     outgoing bonds do""")
         else:
             template[idx] = 1

 T.setRawElem(template.reshape(-1))
 #print T
 Tn = get_ndarray(T)

 for idx in np.ndindex(Tn.shape):
     if template[idx] != Tn[idx]:
         print "Swap gate: bad at {}".format(idx)
         break
 #print "swap gate: they are identical if nothing prints before this"


 return T














 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)
 Swap2=fermionicOPT(A.bond(3), A.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(A.bond(0), A.bond(1))
 Swap3.setLabel([-1,-2,8,9])
 A.permute([1,2,3,4,5],3)
 A.setLabel([1,2,3,4,5])

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
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
 r_d.permute([1,-3,-10],2)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)
 Swap2=fermionicOPT(A.bond(3), A.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(A.bond(0), A.bond(1))
 Swap3.setLabel([-1,-2,8,9])
 A.permute([1,2,3,4,5],3)
 A.setLabel([1,2,3,4,5])


 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([3,8,11,10])
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
 l_d.permute([-10,-13,1],2)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])















 ####################################################
 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)
 Swap2=fermionicOPT(A.bond(3), A.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(A.bond(0), A.bond(1))
 Swap3.setLabel([-1,-2,8,9])
 A.permute([1,2,3,4,5],3)
 A.setLabel([1,2,3,4,5])

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 Swap1=fermionicOPT(A_conj.bond(2), A_conj.bond(4))
 Swap1.setLabel([-3,10,3,-10])
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


 r_d.setLabel([-200,-3,-10])  #new
 r_d.permute([-200,-3,-10],2)










 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)
 Swap2=fermionicOPT(A.bond(3), A.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(A.bond(0), A.bond(1))
 Swap3.setLabel([-1,-2,8,9])
 A.permute([1,2,3,4,5],3)
 A.setLabel([1,2,3,4,5])


 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
# Swap1=fermionicOPT(A_conj.bond(2), A_conj.bond(0))
# Swap1.setLabel([3,8,11,10])
# A_conj=A_conj*Swap1

 #A_conj.permute([10,9,11,-6,7],3)
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

 l_d.setLabel([-10,-13,-400])  #new
 l_d.permute([-10,-13,-400],2)










 mps_boundry_left[Location].setLabel([16,6,-60,18])
 mps_boundry_left[Location+1].setLabel([18,11,-110,25])


 mps_boundry_right[Location].setLabel([17,90,-9,19])
 mps_boundry_right[Location+1].setLabel([19,140,-14,24])












