import numpy as np
from math import pi
from scipy import linalg
import time, itertools
import matplotlib.pyplot as plt
from svd_dgesvd import svd_dgesvd
from numpy.linalg import svd as svd_dgesvd
import copy
import operator

def spin_operators(S):
	d = int(np.rint(2*S + 1))
	dz = np.zeros(d)
	mp = np.zeros(d-1)

	for n in range(d-1):
		dz[n] = S - n
		mp[n] = np.sqrt((2*S - n)*(n + 1))

	dz[d - 1] = -S
	Sp = np.diag(mp,1)
	Sm = np.diag(mp,-1)
	Sx = 0.5*(Sp + Sm)
	Sy = -0.5j*(Sp - Sm)
	Sz = np.diag(dz)

	return Sp,Sm,Sx,Sy,Sz

def make_H_xxx_mpo(L,Jp,Jz,hz):
	s0 = np.eye(2)
	sp = np.array([[0.,1.],[0.,0.]])
	sm = np.array([[0.,0.],[1.,0.]])
	sz = np.array([[0.5,0.],[0.,-0.5]])
	w_list = []
	w2_list = []
	
	for i in range(L):
		w = np.zeros((5,5,2,2),dtype=np.float)
		w[0,:4] = [s0,sp,sm,sz]
		w[0:,4] = [hz[i]*sz, Jp[i]/2.*sm, Jp[i]/2.*sp, Jz[i]*sz, s0]
		w_list.append(np.real(w))
		w = np.tensordot(w,w,axes=(3,2))
		w = np.transpose(w,[0,3,1,4,2,5])
		w = np.reshape(w,(25,25,2,2))
		w2_list.append(np.real(w))
	return w_list,w2_list
	
def mpo_expectation_value(psi,O_mpo_list,LB = None, RB = None):
	L = psi.L
	
	D = O_mpo_list[L-1].shape[1]
	Rp = np.zeros([1,1,D],dtype=float)
	if np.any(RB) == None:
		Rp[0,0,D-1] = 1.
	else:
		Rp[0,0,:] = RB
	
	D = O_mpo_list[0].shape[0]
	Lp = np.zeros([1,1,D],dtype=float)
	if np.any(LB) == None:
		Lp[0,0,0] = 1.
	else:
		Lp[0,0,:] = LB

	for i in np.arange(L-1,-1,-1):
		Rp = np.tensordot(psi.B[i], Rp, axes=(2,0))
		Rp = np.tensordot(O_mpo_list[i], Rp, axes=([1,2],[3,0]))            
		Rp = np.tensordot(np.conj(psi.B[i]), Rp, axes=([0,2],[1,3]))
		Rp = np.transpose(Rp,(2,0,1))

	return np.tensordot(Lp,Rp,axes = ([0,1,2],[0,1,2]))

def mpo_expectation_value_2(psi,O_mpo_list):
	D = O_mpo_list[0].shape[0]
	L = psi.L
	Rp = np.zeros([1,1,D],dtype=float)
	Rp[0,0,D-1] = 1.
	
	Lp = np.zeros([1,1,D],dtype=float)
	Lp[0,0,0] = 1.
	
	for i in np.arange(L-1,-1,-1):
		Rp = np.tensordot(psi.B[i], Rp, axes=(2,0))
		Rp = np.tensordot(O_mpo_list[i], Rp, axes=([1,2],[3,0]))            
		Rp = np.tensordot(np.conj(psi.B[i]), Rp, axes=([0,2],[1,3]))
		Rp = np.transpose(Rp,(2,0,1))
		
	N = np.ones([1,1])
	for i in np.arange(L):
		N = np.tensordot(N,np.conj(psi.B[i]), axes=(1,1))
		N = np.tensordot(N,psi.B[i], axes=([0,1],[1,0]))
	N = np.trace(N)

	return np.tensordot(Lp,Rp,axes = ([0,1,2],[0,1,2])).item(),N
	
def mpo_to_full(H_mpo_list):
	D = H_mpo_list[0].shape[0]
	vL = np.zeros(D)
	vL[0] = 1.
	vR = np.zeros(D)
	vR[D-1] = 1.
	
	L = len(H_mpo_list)
	d =  H_mpo_list[0].shape[2]

	H_full = np.tensordot(vL,H_mpo_list[0],axes=(0,0))
	for i in range(0,L-1):
		H_full = np.tensordot(H_full,H_mpo_list[i+1],axes=(2*i,0))

	H_full = np.tensordot(H_full,vR,axes=(L*2-2,0))
	H_full=np.transpose(H_full,np.hstack([np.arange(0,L*2,2),np.arange(0,L*2,2)+1]))
	H_full=np.reshape(H_full,[d**(L),d**(L)])
	
	return H_full

