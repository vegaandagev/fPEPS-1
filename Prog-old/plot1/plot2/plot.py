from numpy import loadtxt
import matplotlib.pyplot as plt
import numpy as np


R=loadtxt("Z.txt")

X=R[:,0]
Y=R[:,1]
print X, Y

X1=R[:,0]
Y1=R[:,2]
print X1, Y1


#X2=R[:,0]
#Y2=R[:,3]


#X3=R[:,0]
#Y3=R[:,4]

#X4=R[:,0]
#Y4=R[:,5]

#X5=R[:,0]
#Y5=R[:,6]





plt.plot( X, Y,'bs',label='D=4',markersize=np.sqrt(200.) )
plt.plot( X, Y1,'g.',label='D=5',markersize=np.sqrt(200.) )
#plt.plot( X, Y2,'r>',label='simple 9-PEPS',markersize=np.sqrt(200.) )
#plt.plot( X, Y3,'m*',label='Full Z_{2}-symmetric 5-PEPS',markersize=np.sqrt(200.) )
#plt.axhline(y=-0.4375, label='best value', color='r', linestyle='-')
#plt.xlim([2.5,8])
#plt.ylim([-0.415,-0.439])
plt.xlabel('$D_{cut}$', fontsize=15)
plt.ylabel('E', fontsize=15)
plt.legend(loc='center')
plt.savefig('Z.pdf')
plt.clf()

