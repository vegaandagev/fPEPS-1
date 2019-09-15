import numpy as np
import time
import scipy.optimize
from scipy.linalg import expm
import mpo as mpo
import itertools
import copy
import sys 
import cPickle
import os
import gzip

np.set_printoptions(precision=2,linewidth=1000)

############ Define the parameter ###########################
N_layer_list = [0,1,2]
L_list = [8]
eta_list = [1,2,3,4,5,6,7,8,9]
seed_list = range(30)

J = 1.
N_steps = 40
met = "BFGS" 
d = 2
sweep_met = "col"

get_exact = True

class mps():
    def __init__(self):
        self.L = None
        self.chi = None
        self.d = None
        self.B = None
        self.s = None  
        self.dtype = None

    @classmethod	
    def initial_state(cls,L,chi,d):
        psi = cls()
        psi.d = d
        psi.L = L
        psi.dtype = float
        psi.chi = chi
        psi.B = []
        psi.s = []
        for i in range(psi.L):
            psi.B.append(np.zeros((d,d,1,1)))
            psi.B[i][0,0,0,0] = 1.
            psi.B[i][1,1,0,0] = 1.
            #psi.B[i] = np.reshape(psi.B[i],[d*d,1,1])
            psi.s.append(np.ones(1))
        psi.s.append(np.ones(1))
        return psi

def initial_p_list(L,N_layer,d):
    if N_layer > 0:
        p_list = []
        for i in range(N_layer):
            p_list.append([])
            for j in range(L):
                H = np.eye(d*d) + (0.5 - np.random.rand(d*d,d*d))/1000    
                print 'H',H      
                p = []
                for a in range(d*d):
                    for b in range(a+1):
                        p.append(H[a,b])                # Takes the lower triangle of H and stores as parameters
                        #print 'p', p                
                p_list[i].append(np.array(p))
        #print 'p_list[0]',p_list[0][1].shape
    else:
        p_list = []
        for i in range(1):
            p_list.append([])
            for j in range(L):
                H = np.eye(d*d)          
                p = []
                for a in range(d*d):
                    for b in range(a+1):
                        p.append(H[a,b])                # Takes the lower triangle of H and stores as parameters
                p_list[i].append(np.array(p))
    return p_list

def update_mpo(psi,U_list,chi_max):
    " Updates the B and s matrices using U_bond and the TEBD protocol "
    d = psi.d
    L = psi.L
    for p in [0,1]:
        for i_bond in np.arange(p,L-1,2): 
            ia=np.mod(i_bond-1,L)
            ib=i_bond
            ic=np.mod(i_bond+1,L)
            chia = psi.B[ib].shape[2]
            chib = psi.B[ib].shape[3]
            chic = psi.B[ic].shape[3]

            # Construct theta matrix #
            C = np.tensordot(np.reshape(psi.B[ib],[d,d,chia,chib]),np.reshape(psi.B[ic],[d,d,chib,chic]),axes=(3,2))
            if np.any(U_list):
                C = np.tensordot(C,U_list[ib],axes=([1,4],[0,1]))
            else:	
                C = np.tensordot(C,np.reshape(np.eye(d**2,d**2),[d,d,d,d]),axes=([1,4],[0,1]))
                
            theta = np.tensordot(np.diag(psi.s[ib]),C,axes=(1,1))

            theta = np.transpose(theta,(1,4,0,2,5,3))
            theta = np.reshape(theta,(d*d*chia,d*d*chic))
            C = np.transpose(C,(0,4,1,2,5,3))
            C = np.reshape(C,(d*d*chia,d*d*chic))
            # Schmidt decomposition #
            try:
                X, Y, Z = np.linalg.svd(theta,full_matrices=0)
            except np.linalg.linalg.LinAlgError:
                print 'SVD did not converge, diagonalizing theta_dagger*theta'
                Y, Z = np.linalg.eigh(np.dot(theta.conj().T,theta))
                piv = np.argsort(Y)[::-1]        
                Y = np.sqrt(np.abs(Y[piv]))
                Z = np.conj(Z[:,piv].T)
            Z=Z.T
            W = np.dot(C,Z.conj())                    
            chib = np.min([np.sum(Y>10.**(-14)), chi_max])

            # Obtain the new values for B and l #
            invsq = np.sqrt(sum(Y[:chib]**2))/2.**(L/2.)
            psi.s[ic] = Y[0:chib]/invsq 
            psi.B[ib] = np.reshape(W[:,:chib],(d,d,chia,chib))/invsq
            psi.B[ic] = np.transpose(np.reshape(Z[:,:chib],(d,d,chic,chib)),(0,1,3,2))

def gen_U_list(p_list,L,d):
    N_layer = len(p_list)
    U_list = []
    for i in range(N_layer):
        U_list.append([])
        for j in range(L):
            A = np.zeros([d*d,d*d])
            cnt = 0
            for a in range(d*d):
                for b in range(a+1):
                    A[a,b] = p_list[i][j][cnt]
                    A[b,a] = p_list[i][j][cnt]
                    cnt += 1
            U_list[i].append(np.reshape(expm(1j*A),[d,d,d,d]))
    return U_list

def update_U_list(U_list, p , i, j, d):
        A = np.zeros([d*d,d*d])
        cnt = 0
        for a in range(d*d):
            for b in range(a+1):
                A[a,b] = p[cnt]
                A[b,a] = p[cnt]
                cnt += 1
            U_list[i][j] = np.reshape(expm(1j*A),[d,d,d,d])
        return U_list


def gen_H_2_deleted_cell(psi, U_list, H_mpo_list, Rp, Lp, m,i):
        L = psi.L
        d = U_list[0][0].shape[0]
        D = H_mpo_list[0].shape[0]
        N_layer = len(U_list)

        eta = np.zeros([d,d,d])
        for k in range(d):
                eta[k,k,k] = 1.
        
        E = copy.copy(Lp[L-i-2])

        if np.mod(i, 2) ==1:
                for x in range(2):
                        if m == N_layer-1:
                                offset = 2*x
                        else:
                                offset = 4*x
                        E = np.tensordot(E, psi.B[i+1], axes = ([offset+0,offset+1],[3,1]))
                        for j in range(m):
                                E = np.tensordot(E, U_list[j][i], axes = ([offset+0,offset+1], [1,3]))                        
                        for j in range(m+1, N_layer-1):
                                E = np.tensordot(E, U_list[j][i], axes = ([offset+2,offset+3], [1,3]))
                        if m!= N_layer-1:
                                E = np.tensordot(E, U_list[-1][i], axes = ([offset+2], [1]))
                                E = np.tensordot(E, H_mpo_list[i+1], axes = ([-1,offset+2], [2,1]))
                                E = np.tensordot(E, np.conj(U_list[-1][i].T), axes = ([-1,offset+2], [0,2]))
                                for j in range(N_layer-2, m, -1):
                                        E = np.tensordot(E, np.conj(U_list[j][i].T), axes = ([offset+2,offset+3], [0,2]))
                                for j in range(m-1, -1, -1):
                                        E = np.tensordot(E, np.conj(U_list[j][i].T), axes = ([offset+4,offset+5], [0,2]))
                                E = np.tensordot(E, psi.B[i+1], axes = ([offset+4,offset+5],[0,3]))
                                E = np.tensordot(E, eta, axes = ([4*N_layer+3+4+x,-2], [0,1]))

                        else:

                                E = np.tensordot(E, np.transpose(H_mpo_list[i+1],axes=[2,1,0,3]), axes = ([offset+1], [1]))
                                for j in range(N_layer-2, -1, -1):
                                        E = np.tensordot(E, np.conj(U_list[j][i].T), axes = ([offset+2,offset+3], [0,2]))
                                E = np.tensordot(E, psi.B[i+1], axes = ([offset+2,offset+3],[0,3]))
                                E = np.tensordot(E, eta, axes = ([4*N_layer+3+2+x,-2], [0,1]))

                E = np.trace(E, axis1 = 4*N_layer+3+offset, axis2= -1)

                if m ==N_layer-1:
                        if N_layer == 2:
                                inds = [0, 1,2,5, 8,9,10]
                                i2 = inds + [11+ii for ii in inds]
                                inds2 = [4,5,6,8,10,11, 12]
                                i1 = inds2 + [9 + ii for ii in inds2]
                                H_2_del = np.tensordot(E, Rp[i], axes = (i1, i2))

                        elif N_layer ==1:
                                inds = [0, 3,6, 7, 10, 13]
                                inds2 = [4, 6,8, 9, 11, 13]
                                H_2_del = np.tensordot(E, Rp[i], axes = (inds2, inds))
                        H_2_del = np.transpose(H_2_del, axes = [8,9,0,4,10,11,5,1, 12,13,2,6, 14,15,7,3])
                        

                elif m ==0:
                        
                        inds = [0, 3,4,5,6,7,10]
                        i2 = inds + [11+ii for ii in inds]
                        inds2 = [8,9,10,11,12,13,14]
                        i1 = inds2 + [7 + ii for ii in inds2]
                        H_2_del = np.tensordot(E, Rp[i], axes = (i1, i2))
                        H_2_del = np.transpose(H_2_del, axes = [8,9,0,1,10,11,2,3, 12,13,4,5, 14,15,6,7])


        else:
                if m!=0:
                        for x in range(2):
                                E = np.tensordot(E, psi.B[i+1], axes = (4*x+0,3))
                                E = np.tensordot(E, U_list[0][i], axes = ([-2,4*x+0], [1,3]))
                                for j in range(1, m-1):
                                        E = np.tensordot(E, U_list[j][i], axes = ([4*x+0,4*x+1], [1,3]))
                                for j in range(m+1, N_layer):
                                        E = np.tensordot(E, U_list[j][i], axes = ([4*x+2,4*x+3], [1,3]))
                                E = np.tensordot(E, H_mpo_list[i+1], axes = ([4*x+2,4*x+3,4*x+4], [2,1,3]))
                                for j in range(N_layer-1, m, -1):
                                        E = np.tensordot(E, np.conj(U_list[j][i].T), axes = ([4*x+2,4*x+3], [0,2]))
                                for j in range(m-1, 0, -1):
                                        E = np.tensordot(E, np.conj(U_list[j][i].T), axes = ([4*x+4,4*x+5], [0,2]))
                                E = np.tensordot(E, np.conj(U_list[0][i].T), axes = ([4*x+4], [0]))
                                E = np.tensordot(E, psi.B[i+1], axes = ([4*x+4,-2], [3,0]))
                                E = np.tensordot(E, eta, axes = ([4*N_layer+3+x+4,-2], [0,1]))
                        E = np.trace(E, axis1 = 4*N_layer+3+4, axis2= -1)

                        inds = [0,1,2, 5,8, 9,10]
                        i1 = [8+ ii for ii in range(7)]+ [15+ ii for ii in range(7)]
                        i2 = inds + [ 11 +ii for ii in inds]
                        
                        H_2_del = np.tensordot(E, Rp[i], axes = (i1, i2))
                        H_2_del = np.transpose(H_2_del, axes = [8,9,0,1,10,11, 2,3,12,13, 4,5,14,15, 6,7])

                elif m == 0:
                        for x in range(2):
                                E = np.tensordot(E, np.transpose(psi.B[i+1],axes = [0,2,1,3]),  axes = (2*x+0,3))
                                for j in range(1, N_layer):
                                        E = np.tensordot(E, U_list[j][i], axes = ([2*x+1,2*x+2], [1,3]))
                                E = np.tensordot(E, H_mpo_list[i+1], axes = ([2*x+1,2*x+2,2*x+3], [2,1,3]))
                                for j in range(N_layer-1):
                                        E = np.tensordot(E, np.conj(U_list[N_layer-j-1][i].T), axes = ([2*x+1,2*x+2], [0,2]))
                                E = np.tensordot(E, psi.B[i+1], axes = ([2*x+2], [3]))
                                E = np.tensordot(E, eta, axes = ([4*N_layer+3+2+x,-2], [0,1]))

                        E = np.trace(E, axis1 = 4*N_layer+3+2, axis2= -1)

                        if N_layer == 2:
                                inds = [0,3,4,5,6,7,10]
                                inds2 = [4, 6, 7,8,9,10, 12]
                                i1 =  inds2  + [9 + ii for ii in inds2]
                                i2 =inds + [11+ ii for ii in inds]
                                H_2_del = np.tensordot(E, Rp[i], axes = (i1, i2))
                                
                        elif N_layer == 1:
                                inds = [0, 3, 6, 7, 10, 13]
                                inds2 = [4, 6, 8, 9, 11, 13]
                                H_2_del = np.tensordot(E, Rp[i], axes = (inds2, inds))
                                
                        H_2_del = np.transpose(H_2_del, axes = [8,9,4,0,10,11,1,5, 12,13,6,2, 14,15,3,7])

        return H_2_del
        
def functional(p,H2_del, d):
        A = np.zeros([d*d,d*d])
        cnt = 0
        for a in range(d*d):
                for b in range(a+1):
                        A[a,b] = p[cnt]
                        A[b,a] = p[cnt]
                        cnt += 1
        U = np.reshape(expm(1j*A),[d,d,d,d])
        E = np.tensordot(H2_del, U, axes = ([0,1,2,3], [0,2,1,3]))
        E = np.tensordot(E, np.conj(U.T), axes = ([0,1,2,3], [1,3,0,2]))
        E = np.tensordot(E, U, axes = ([0,1,2,3], [0,2,1,3]))
        E = np.tensordot(E, np.conj(U.T), axes = ([0,1,2,3], [1,3,0,2]))
        return -np.abs(E)
        
def get_psi(p_list,chi,L):
    N_layer = len(p_list)
    psi = mps.initial_state(L,chi,2)

    U_list = gen_U_list(p_list,L,2)
    for i in range(N_layer):
        print i,chi
        update_mpo(psi,U_list[i],chi)
    return psi
    
def update_right_products_H_2(psi, U_list, H_mpo_list, Rp, i):
        L = psi.L
        d = U_list[0][0].shape[0]
        D = H_mpo_list[0].shape[0]
        N_layer = len(U_list)

        eta = np.zeros([d,d,d])
        for j in range(d):
                eta[j,j,j] = 1.

        E = copy.copy(Rp[i])
                
        if np.mod(i, 2) ==1:
                for x in range(2):
                        E = np.tensordot(E, psi.B[i+1], axes = (0,2))
                        for j in range(N_layer):
                                E = np.tensordot(E, U_list[j][i], axes = ([0,1], [0,2]))
                        E = np.tensordot(E, H_mpo_list[i+1], axes = ([0, -1], [0,2]))
                        E = np.tensordot(E, np.conj(U_list[-1][i].T), axes = ([0,1,-1], [1,3,0]))
                        for j in range(1, N_layer):
                                E = np.tensordot(E, np.conj(U_list[N_layer-j-1][i].T), axes = ([0,1], [1,3]))
                        E = np.tensordot(E, psi.B[i+1], axes = (0,2))
                        E = np.tensordot(E, eta, axes = ([4*N_layer+3+x,-2], [0,1]))

                E = np.trace(E, axis1 = 4*N_layer+3, axis2= -1)
                E = np.transpose(E, [1,0]+ range(2, len(E.shape)))                
                E = np.transpose(E, range(4*N_layer+3)+ [4*N_layer+4,4*N_layer+3]+ range(4*N_layer+5, len(E.shape)))

        else:
                for x in range(2):
                        E = np.tensordot(E, psi.B[i+1], axes = (0,2))
                        E = np.tensordot(E, U_list[0][i], axes = ([0,1,-2], [0,2,1]))
                        for j in range(1,N_layer):
                                E = np.tensordot(E, U_list[j][i], axes = ([0,1], [0,2]))
                        E = np.tensordot(E, np.transpose(H_mpo_list[i+1], [0, 2,1,3]), axes = ([0], [0]))
                        for j in range(N_layer):
                                E = np.tensordot(E, np.conj(U_list[N_layer-j-1][i].T), axes = ([0,1], [1,3]))
                        E = np.tensordot(E, psi.B[i+1], axes = ([0,-1],[2,0]))
                        E = np.tensordot(E, eta, axes = ([4*N_layer+3+x,-2], [0,1]))
                E = np.trace(E, axis1 = 4*N_layer+3, axis2= -1)

        if len(Rp) > i+1:
                Rp[i+1] = E
        else:
                Rp.append(E)

        return Rp
        
                                               

def right_products_H_2(psi, U_list, H_mpo_list):
        L = psi.L
        d = U_list[0][0].shape[0]
        D = H_mpo_list[0].shape[0]
        N_layer = len(U_list)

        eta = np.zeros([d,d,d])
        for i in range(d):
                eta[i,i,i] = 1.

        E = np.zeros([1,D,1])
        E[0,0,0]=1
        Rp = np.zeros([d]*(2*(N_layer+1)))
        pairs = itertools.product([0,1], repeat = (N_layer+1))
        for x in pairs:
                Rp[tuple([val for val in x for _ in (0, 1)])]=1
        Rp = np.tensordot( psi.B[0], Rp,  axes = (1, 0))
        E = np.tensordot(E, Rp, axes = (0, 1))
        E = np.tensordot( E, H_mpo_list[0], axes = ([0,2*N_layer+4], [0, 2]))
        E = np.tensordot(E, Rp.T, axes = ([0,-1], [-2,0]))
        E = np.tensordot(E, eta, axes = ([0, -1], [0,1]))
        E = np.tensordot(E, E, axes = ([-1], [-1]))

        Rpp=[E]
                
        for i in range(L-1):
                Rpp = update_right_products_H_2(psi, U_list, H_mpo_list, Rpp, i) 

        return Rpp


def update_left_products_H_2(psi, U_list, H_mpo_list, Lp, i):
        L = psi.L
        d = U_list[0][0].shape[0]
        D = H_mpo_list[0].shape[0]
        N_layer = len(U_list)

        eta = np.zeros([d,d,d])
        for j in range(d):
                eta[j,j,j] = 1.

        E = copy.copy(Lp[L-2-i])

        if np.mod(i, 2) ==1:
                for x in range(2):
                        E = np.tensordot(E, psi.B[i+1], axes = ([0,1],[3,1]))
                        for j in range(N_layer-1):
                                E = np.tensordot(E, U_list[j][i], axes = ([0,1], [1,3]))
                        E = np.tensordot(E, U_list[-1][i], axes = ([0], [1]))
                        E = np.tensordot(E, H_mpo_list[i+1], axes = ([-1,0], [2,1]))
                        E = np.tensordot(E, np.conj(U_list[-1][i].T), axes = ([-1,0], [0,2]))

                        for j in range(1, N_layer):
                                E = np.tensordot(E, np.conj(U_list[N_layer-j-1][i].T), axes = ([0,1], [0,2]))
                        E = np.tensordot(E, psi.B[i+1], axes = ([0,1],[0,3]))
                        E = np.tensordot(E, eta, axes = ([4*N_layer+3+x,-2], [0,1]))
                E = np.trace(E, axis1 = 4*N_layer+3, axis2= -1)

        else:
                for x in range(2):
                        E = np.tensordot(E, psi.B[i+1], axes = (0,3))
                        E = np.tensordot(E, U_list[0][i], axes = ([-2,0], [1,3]))
                        for j in range(1, N_layer):
                                E = np.tensordot(E, U_list[j][i], axes = ([0,1], [1,3]))
                        E = np.tensordot(E, H_mpo_list[i+1], axes = ([0,1,2], [2,1,3]))
                        for j in range(N_layer-1):
                                E = np.tensordot(E, np.conj(U_list[N_layer-j-1][i].T), axes = ([0,1], [0,2]))
                        E = np.tensordot(E, np.conj(U_list[0][i].T), axes = ([0], [0]))
                        E = np.tensordot(E, psi.B[i+1], axes = ([0,-2], [3,0]))
                        E = np.tensordot(E, eta, axes = ([4*N_layer+3+x,-2], [0,1]))

                E = np.trace(E, axis1 = 4*N_layer+3, axis2= -1)

        if len(Lp) > L-1 - i:
                Lp[L-1-i] = E
        else:
                Lp.append(E)
                

        return Lp
        


        
def left_products_H_2(psi, U_list, H_mpo_list):
        L = psi.L
        d = U_list[0][0].shape[0]
        D = H_mpo_list[0].shape[0]
        N_layer = len(U_list)

        eta = np.zeros([d,d,d])
        for i in range(d):
                eta[i,i,i] = 1.

        Lp = np.zeros([1]+ [d]*2*N_layer + [D] + [d]*2*N_layer+[1]+ [1]+ [d]*2*N_layer + [D] + [d]*2*N_layer+[1])
        pairs = itertools.product([0,1], repeat = (4*N_layer))
        for x in pairs:
                ind = [val for val in x for _ in (0, 1)]
                ll = len(ind)
                ind = [0] + ind[0:ll/4]+[D-1]+ ind[ll/4 : ll/2]+[0] + [0] + ind[ll/2 : 3*ll/4]+[D-1]+ ind[3*ll/4 : ] + [0]
                Lp[tuple(ind)]=1
        E = Lp
        Lpp = [E]

        for i in range(L-2,-1,-1):
                Lpp = update_left_products_H_2(psi, U_list, H_mpo_list, Lpp, i) 

        return Lpp


	
def dec2bin(d,n):
    a =  [0]*(n - len(list(bin(d)[2:]))) + list(bin(d)[2:])
    return np.array(a,dtype=int)

def exp_H_2(psi,H_mpo_list):
	d = psi.B[0].shape[0]
	D = H_mpo_list[0].shape[0]
	chi1 = psi.B[0].shape[2]
	chi2 = psi.B[0].shape[3]
	
	eta = np.zeros([d,d,d])
	for i in range(d):
		eta[i,i,i] = 1.
		
	E = np.zeros([chi1,D,chi1,chi1,D,chi1]) #a,x,ap,b,y,bp
	E[0,0,0,0,0,0] = 1.
	for i in range(psi.L):
		chi1 = psi.B[i].shape[2]
		chi2 = psi.B[i].shape[3]
		B = psi.B[i]
		
		E = np.tensordot(E,B,axes=[0,2])
		E = np.tensordot(E,H_mpo_list[i],axes=[[0,6],[0,2]])
		E = np.tensordot(E,np.conj(B),axes=[[0,7],[2,1]])
		E = np.tensordot(E,eta,axes=[[3,6],[0,1]])
		
		E = np.tensordot(E,B,axes=[0,2])
		E = np.tensordot(E,H_mpo_list[i],axes=[[0,7],[0,2]])
		E = np.tensordot(E,np.conj(B),axes=[[0,8],[2,1]])
		E = np.tensordot(E,eta,axes=[[4,7],[0,1]])
		
		E = np.trace(E,axis1=3,axis2=7)
	return E[0,D-1,0,0,D-1,0]
    
def mpo_expectation_value(B_list,O_mpo_list):
    D = O_mpo_list[0].shape[0]
    L = len(B_list)
    Rp = np.zeros([1,1,D],dtype=float)
    Rp[0,0,D-1] = 1.
	
    Lp = np.zeros([1,1,D],dtype=float)
    Lp[0,0,0] = 1.

    for i in np.arange(L-1,-1,-1):
        Rp = np.tensordot(B_list[i], Rp, axes=(2,0))
        Rp = np.tensordot(O_mpo_list[i], Rp, axes=([1,2],[3,0]))
        Rp = np.tensordot(np.conj(B_list[i]), Rp, axes=([0,2],[1,3]))
        Rp = np.transpose(Rp,(2,0,1))
    return np.tensordot(Lp,Rp,axes = ([0,1,2],[0,1,2])).item()
	
def get_energies(psi,H_mpo_list):
    L = psi.L
    E_list = []
    d,d,chi1,chi2 = psi.B[0].shape
		
    for j in range(d**L):
        alpha = dec2bin(j,L)
        B_list = []
        for l in range(L):
            chi1,chi2 = psi.B[l].shape[2],psi.B[l].shape[3]
            M = np.reshape(psi.B[l],[d,d,chi1,chi2])
            B_list.append(M[alpha[l],:,:,:])
	
        N = np.ones([1,1])
        for i in np.arange(L):
            N = np.tensordot(N,np.conj(B_list[i]), axes=(1,1))
            N = np.tensordot(N,B_list[i], axes=([0,1],[1,0]))
        N = np.trace(N)
        E_list.append(np.real(mpo_expectation_value(B_list,H_mpo_list)))
    return E_list

def trace_H(H2_mpo_list):
    D = H2_mpo_list[0].shape[0]
    L = len(H_mpo_list)
    E = np.zeros([D])
    E[0] = 1.
    for i in range(L):
        M = np.trace(H2_mpo_list[i],axis1=2,axis2=3)
        E = np.tensordot(E,M,axes=(0,0))
    return E[D-1]
    
def output_it():
    psi = get_psi(p_list,chi,L)
    t_all.append(time.time() - t0)
    
    if N_layer == 0:
        v_all.append(H2 - exp_H_2(psi,H_mpo_list))
    else:
        v_all.append(H2 + var_func(p))
    
    data = {}
    if get_exact:
        E_list = np.sort(get_energies(psi,H_mpo_list))
        E_all.append(E_list)
        print L,len(p_list),"%06.1f"%(time.time() - t0), "%.3f"%(np.min(E_list)),"%.3f"%(np.max(E_list)),"%02.3f"%(v_all[-1])," : ",psi.s[L/2]
        np.savetxt(fn, E_all)
        data['E'] = E_all
    else:
        print L,len(p_list),"%06.1f"%(time.time() - t0),"%02.3f"%(v_all[-1])," : ",psi.s[L/2]
        
    data['t'] = t_all
    data['v'] = v_all
    data['L'] = L
    data['eta'] = eta
    data['hz'] = hz
    data['p_list'] = p_list
    f = gzip.open(fo,'wb')
    cPickle.dump(data,f)
    f.close()
    

for seed in seed_list:
    np.random.seed(seed)
    print 'seed', seed
    
    for L in L_list:
        Jz = L*[J]
        Jp = L*[J]
        hz0 = 2*(0.5-np.random.rand(L))
        
        for eta in eta_list:
            hz = eta*hz0
            print 'hi', hz
            
            for N_layer in N_layer_list:
                chi = d**(2*N_layer)
                
                fn = "E_all_sf_trunc_L_" + str(L) + "_N_layer_" + str(N_layer) + "_eta_%02.2f"%eta +"_Jp_%02.2f"%Jp[0] + "_" + met + "_seed_" + str(seed)
                fo = "data_sf_trunc_L_" + str(L) + "_N_layer_" + str(N_layer) + "_eta_%02.2f"%eta +"_Jp_%02.2f"%Jp[0] + "_" + met  + "_seed_" + str(seed) + ".pklz"
    
                if os.path.isfile(fo)==False:
                    open(fo, 'w').close() 
                    H_mpo_list,H2_mpo_list = mpo.make_H_xxx_mpo(L,Jp,Jz,hz)

                    ########### Get the Full Hamiltonian ##########################
                    if get_exact == True:
                        H_full = mpo.mpo_to_full(H_mpo_list)
                        e_exact = np.sort(np.linalg.eigvals(H_full))
                        print L,0,"%06.1f"%0, "%.3f"%(np.min(e_exact)),"%.3f"%(np.max(e_exact)),"%02.3f"%0
                        E_all = [e_exact]
            
                    t_all = [0]
                    v_all = [0]

                    ########### Now find the approximation ########################
                    psi = mps.initial_state(L,1,2)
                    H2 = trace_H(H2_mpo_list)
      
                    i=j=0
                    var_func = lambda p: functional(p, T, d)
                    p_list = initial_p_list(L,N_layer,d)
                    U_list = gen_U_list(p_list, L, d)

                    Lp= left_products_H_2(psi, U_list, H_mpo_list)
                    Rp= right_products_H_2(psi, U_list, H_mpo_list)
        
                    t0 = time.time()
                    if N_layer > 0:
                        if sweep_met == "col":
                            for k in range(N_steps):       
                                if np.mod(k, 2) == 0:
                                    jRange = range(L-1)
                                else:
                                    jRange = range(L-2, -1, -1)           
                                for j in jRange:
                                    for i in range(N_layer):
                                        T = gen_H_2_deleted_cell(psi, U_list, H_mpo_list, Rp, Lp, i, j)
                                        p0 = p_list[i][j]
                                        p = scipy.optimize.minimize(var_func,x0=p0,method=met)['x']

                                        p_list[i][j] = p
                                        U_list = update_U_list(U_list,p, i,j,d)

                                        if np.mod(k, 2) == 0:
                                            Rp= update_right_products_H_2(psi, U_list, H_mpo_list, Rp, j)
                                        else: 
                                            Lp= update_left_products_H_2(psi, U_list, H_mpo_list, Lp, j)
                           
                                        sys.stdout.flush()
                                output_it()

                        else:
                            Nsweep = 0
                            for k in range(N_steps):
                                for i in range(N_layer):
                        
                                    if np.mod(Nsweep, 2) == 0:
                                        jRange = range(L-1)
                                    else:
                                        jRange = range(L-2, -1, -1)
                            
                                    for j in jRange: 
                                        T = gen_H_2_deleted_cell(psi, U_list, H_mpo_list, Rp, Lp, i, j)
                                        p0 = p_list[i][j]
                                        p = scipy.optimize.minimize(var_func,x0=p0,method=met)['x']

                                        p_list[i][j] = p
                                        U_list = update_U_list(U_list,p, i,j,d)

                                        if np.mod(Nsweep, 2) == 0:
                                            Rp= update_right_products_H_2(psi, U_list, H_mpo_list, Rp, j)
                                        else: 
                                            Lp= update_left_products_H_2(psi, U_list, H_mpo_list, Lp, j)
                           
                                        sys.stdout.flush()
                                    Nsweep = Nsweep+1 
                                output_it()
                    else:
                        output_it()
                        
                                    
                        
