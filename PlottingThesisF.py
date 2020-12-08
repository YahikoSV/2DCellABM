# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 18:21:02 2020

@author: Potatsu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#UPLOAD THE FFING MATRIXES

time = np.loadtxt("time.txt")
A = np.loadtxt("grow10_matrix.txt")
#N = np.loadtxt("deform80_3_graph.txt")

#N[0] = 1


#Matrixlist
y = 9
x = int(len(A)/y)
Matrixlist = A

Matrixfixed = np.zeros((x,y+1))
for i in range (1,y+1):
    Matrixfixed[:,i] = Matrixlist[(i-1)*x:i*x]
 
font = {'size'   : 30}
matplotlib.rc('font', **font)

plt.figure(figsize=(200,200))

#plt.plot(time,N)
#plt.title('No of cells vs time Plot')
#plt.ylabel('No. of cells')
#plt.xlabel('Accepted time steps')

plt.show()

d = 4





fig, ax = plt.subplots(figsize=(15,15)) # note we must use plt.subplots, not plt.subplot

# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

for j in range (0,len(Matrixfixed)):
    if Matrixfixed[j,5] == 2:
        circle1 = plt.Circle((Matrixfixed[j,1],Matrixfixed[j,2]), Matrixfixed[j,4], color = 'blue')
        circle2 = plt.Circle((Matrixfixed[j,6],Matrixfixed[j,7]), Matrixfixed[j,4], color = 'blue')
        ax.add_artist(circle1)
        ax.add_artist(circle2)
    elif Matrixfixed[j,5] == 1:
        circle1 = plt.Circle((Matrixfixed[j,1],Matrixfixed[j,2]), Matrixfixed[j,4], color = 'red')
        ax.add_artist(circle1)
    else:
        dummy = 1
plt.title('r_g = 0.5, r_d = 0.5, \u03B5 = -20') 
plt.xlim(-90,90)
plt.ylim(-90,90)
plt.xlabel('x')
plt.ylabel('y')
plt.show()       

'''
#Graphing Gowth Rates:

time2 = np.copy(time)
Control  = np.loadtxt("Control_Graph.txt")
Growth1 = np.loadtxt("Growth1_Graph.txt")
Growth2 = np.loadtxt("Growth2_Graph.txt")
Deform1 = np.loadtxt("Deform1_Graph.txt")
Deform2 = np.loadtxt("Deform2_Graph.txt")
Attract1 = np.loadtxt("Attract1_Graph.txt")
Attract2 = np.loadtxt("Attract2_Graph.txt")

Control[0] = 1
Growth1[0] = 1
Growth2[0] = 1
Deform1[0] = 1
Deform2[0] = 1
Attract1[0] = 1
Attract2[0] = 1

plt.plot(time,Growth1)
plt.plot(time,Control)
plt.plot(time,Growth2)
plt.title('Cell Population with varying rates of growth')
plt.legend(['r_g = .2','r_g = .5','r_g = .8'])
plt.show()

plt.plot(time,Deform1)
plt.plot(time,Control)
plt.plot(time,Deform2)
plt.legend(['r_d = .2','r_d = .5','r_d = .8'])
plt.title('Cell Population with varying rates of defromation')
plt.show()

plt.plot(time,Attract1)
plt.plot(time,Attract2)
plt.plot(time,Control)
plt.legend(['v_pot = 0.0','v_pot = -5.0','v_pot = -20.0'])
plt.title('Cell Population with varying attraction strengh')
plt.show()


Control  = np.loadtxt("Control_Graph.txt")
Control[0] = 1
plt.plot(time,Control)
plt.title('Control Cell Population (r_g = .5, r_d = .5, v_pot = -20.0)')
plt.show()
'''


