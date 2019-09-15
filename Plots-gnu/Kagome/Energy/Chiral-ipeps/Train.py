import numpy as np
import matplotlib.pyplot as plt
import codecs

filecp = codecs.open('chiral.txt')
file_data = np.loadtxt(filecp, usecols=(0,1),skiprows=1)
print(file_data)
x = file_data[:,0]
print x, x[3]

y = file_data[:,1]
print y, y[2]

Dfy=[]
Dsy=[]
print len(x)
for i in xrange(len(x)):
 if i < int(len(x)-1):
  #print  i, y[i]-y[i+1], x[i]-x[i+1]
  Dfy.append((y[i]-y[i+1])/(x[i]-x[i+1]))
print Dfy 


print len(x)
for i in xrange(len(x)):
 if i < int(len(x)-2):
  #print  i, Dfy[i]-Dfy[i+1], x[i]-x[i+1]
  Dsy.append((Dfy[i]-Dfy[i+1])/(x[i]-x[i+1]))
print Dsy 


file = open("Der.txt", "w")
for index in range(len(x)-2):
    file.write(str(x[index]) + "    " + str(Dfy[index])+"    " + str(Dsy[index]) + "\n")
file.flush()
file.close()
