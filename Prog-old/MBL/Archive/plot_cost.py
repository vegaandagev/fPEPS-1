import numpy as np
import pylab as pl
import cPickle
from scipy.stats import sem
import os
import gzip

L = 8
N_seed = 30
Jp = 1.
eta_list = [0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,12,14,16,24,48]
pl.figure( figsize=(5, 3.2))
lp = []
n_max = 40
met = "BFGS" 
n = 1

for N_layer in [0,1,2]:
    C = []
    for seed in range(N_seed):
        v = []
        for eta in eta_list:
            fo = "data_sf_trunc_L_" + str(L) + "_N_layer_" + str(N_layer) + "_eta_%02.2f"%eta + "_Jp_%02.2f"%Jp + "_" + met + "_seed_" + str(seed) + ".pklz"
            f = gzip.open(fo,'rb')
            data = cPickle.load(f)
            W = np.max(data['E']-np.min(data['E']))
            W = 1.
            f.close()
            if len(data['v']) >= n:
                v.append(data['v'][1].real/2.**L/W**2)
        if len(v) == len(eta_list):
            C.append(v)
    if N_layer==0:
        pl.errorbar(eta_list,np.mean(C,axis=0),sem(C,axis=0),color='b')
        p, = pl.semilogx(eta_list,np.mean(C,axis=0),'-^b')
    if N_layer==1:
        pl.errorbar(eta_list,np.mean(C,axis=0),sem(C,axis=0),color='g')
        p, = pl.semilogx(eta_list,np.mean(C,axis=0),'-vg')
    if N_layer==2:
        pl.errorbar(eta_list,np.mean(C,axis=0),sem(C,axis=0),color='r')
        p, = pl.semilogx(eta_list,np.mean(C,axis=0),'-sr')
    lp.append(p)

pl.xlabel('$W$')
pl.title('$L = %.0f$'%(L))
pl.yticks([0,0.5,1.0])
#pl.ylabel('$\overline{\\langle H^2\\rangle_{\\alpha} - \\langle H\\rangle_{\\alpha}^2}$')
#pl.legend(lp,['$N_{\\mathrm{layer}} = 0$','$N_{\\mathrm{layer}} = 1$','$N_{\\mathrm{layer}} = 2$'],loc='upper right')
pl.subplots_adjust(bottom=0.15, left = 0.2, right=0.95, top=0.9)
pl.savefig('fig_cost_L_%.0f.pdf'%L)
pl.show()