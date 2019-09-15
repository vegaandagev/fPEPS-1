import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import sys
import time
import main 
import ipeps 

h_list=[]
acc_list=[]
h=3.5

for q in xrange(1):
 t0=time.time()

 h=h+(0.0)
#for q in reversed(xrange(20)):
# h=h-(0.20)
 print q, "h=", h
 a_u,b_u,c_u,d_u=ipeps.ipeps_function(h)
 acc=main.QR_function(a_u,b_u,c_u,d_u) 
 print "accuracy", acc
 h_list.append(h)
 acc_list.append(acc)
 print time.time() - t0, "whole-update"

 
file = open("accIsing.txt", "w")
for index in range(len(h_list)):
  file.write(str(h_list[index]) + " " + str(acc_list[index])+" "+ "\n")
file.close()

