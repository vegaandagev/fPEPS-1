import numpy as np
import matplotlib.pyplot as plt
import codecs

filecp = codecs.open('CorrelationHX.txt')
file_data = np.loadtxt(filecp, usecols=(0,1),skiprows=1)
print(file_data)
x = file_data[:,0]
print x
y = file_data[:,1]
print "y", y



filecp = codecs.open('CorrelationHZ.txt')
file_data = np.loadtxt(filecp, usecols=(0,1),skiprows=1)
y1 = file_data[:,1]
print "y1", y1



filecp = codecs.open('CorrelationHY.txt')
file_data = np.loadtxt(filecp, usecols=(0,1),skiprows=1)
y2 = file_data[:,1]
print "y2", y2



Dfy=[]
for i in xrange(len(y)):
  Dfy.append(y[i]+y1[i]+y2[i])
print Dfy 



file = open("Corr.txt", "w")
for index in range(len(x)):
    file.write(str(x[index]) + "    " + str(Dfy[index]) + "\n")
file.flush()
file.close()
