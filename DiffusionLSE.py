import numpy as np
import scipy as sp
from random import *
from matplotlib import pyplot as plt 
import time
start_time = time.time()


Tau=100                 
#r=(random())              # Random number generator

# Define storage arrays
N=[]
T=[]

#Monte -Carlo


Simul=10000                   
for j in range(Simul):

   #initial Conditions since photon  ejected from the centre
   #x,y,z=position
   #R=Radius
   #ns=number of scatters
   #t=time 
        
        
        
    R,x,y,z,ns,t=(0,0,0,0,0,0)

 
    
   
    

    for i in range(10000):
        phi=(uniform(0, 360))          #random value of cos Theta between -1 and 1
        theta=(uniform(-1, 1))         #azimuthal angle between 0 and 2 pi
 
    #Positions friom random scattering                                
        x=x + r*np.cos(phi)
        y=y + r*np.sin(phi)
        z=z + r*theta
        R=(x**2 + y**2 + z**2)**(1/2)     #Escape radius
        ns=ns+1                            
        t=t+r                             #Time to photon escape
        
        N_Tau=r/R*Tau                     #New tau dependent on R
        r=(1-np.exp(-N_Tau))/N_Tau        # N_tau step length
        if R>=2:
            break
            
    T.append(t)
    N.append(ns)
#Plotting the Distributions
fig=plt.figure(figsize=(15,10))
plt.hist(T,bins=150, label='Tau = 100')


plt.legend()
plt.xlabel('T')                 
plt.ylabel('Distrib')
plt.show()
plt.savefig('distributiontime100random')

print (time.time() - start_time)
