
# coding: utf-8

# In[103]:

import scipy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit  
#import pyneb as pn
import h5py
#get_ipython().magic('matplotlib nbagg')

T_eff=5770
sigma=5.67e-5

def getS(tau):
    T4 = (3.0/4.0) * T_eff**4 * (tau+ 2.0/3.0)
    S = (sigma/np.pi)*T4
    return S
  
tau = np.logspace(-4,1,51)
#print(tau)
#print(len(tau))

S=getS(tau)

'''plt.plot(tau, S)  
plt.ylabel('source function')
plt.xlabel('tau')
plt.savefig('S(tau)')'''


# In[104]:

plt.close("all")
#get_ipython().magic('matplotlib nbagg')

y = np.linspace(-5,5,21)          #y=j
#print(y)

def getPhi(y):
    phi = np.exp(-(y**2))
    return phi                    #phi=k

phi = getPhi(y)
'''plt.figure()
plt.plot(y, phi)
plt.xlabel('doppler width from line center')
plt.ylabel('doppler line profile')
plt.savefig('Dopplerlineprofile')'''


# In[105]:

tau2=np.outer(tau, 10*phi+1)            #makes Tau2(jk)=tau(j) phi(k) ##frequency and depth dependent..

'''fig, ax = plt.subplots()
ax.imshow(tau2)
plt.savefig('tau2')

#print(tau2.shape)

fig2, ax2 = plt.subplots()
ax2.plot(tau2[0,:])
ax2.plot(tau2[-1,:])
ax2.set_yscale("log")
plt.ylabel('tau')
plt.xlabel('tau2')
plt.savefig('tautau2')'''


# In[106]:

#def transporteq(tau2,y)
#j=frequency index, k=depth index

j=0
I_nu=sum(S*np.exp(-tau2[:,j])*np.log(10)*tau2[:,j]*0.1)
print('I_nu=',I_nu)

def function_K(tau2):
    K=np.transpose(np.exp(-tau2)*np.log(10)*tau2*0.1)
    return K

K = function_K(tau2) 
'''plt.close("all")
plt.imshow(np.transpose(K))
plt.xlabel('tau2')
plt.ylabel('K')
plt.show()
plt.savefig('K')'''


# In[107]:

I=np.matmul(K,S)
'''plt.close("all")
plt.plot(y,I)
plt.ylabel('Intensity')
plt.xlabel('Frequancy variable y')
plt.show()
plt.savefig('yI')'''


# In[108]:

taugray_wheretauisone=np.zeros((21))    
print(taugray_wheretauisone.shape)         # 21
print(tau.shape)                           # 51
print(tau2.shape)                          # 51,21

for i in range(21):
    taugray_wheretauisone[i]=np.interp(1,tau2[:,i],tau)



'''plt.close("all")
plt.plot(taugray_wheretauisone,I, color='r', label='Eddington-Barbier')
plt.plot(tau, S, color='g', label='Source function') 
plt.xlabel('gray optical depth')

plt.ylabel('Intensity')
plt.legend()
plt.show()
plt.savefig('eddingtonbarbierapprox')'''


# In[109]:

#1.0.4 Inversion
S_inv=np.zeros((21))
KT=np.transpose(K)
KTK=np.dot(KT,K)
KTK_inv=np.linalg.inv(KTK)
KTK_invKT=np.dot(KTK_inv,KT)
S_inv=np.dot(KTK_invKT,I)

plt.close("all")
plt.plot(tau, S_inv) #color='r' #label='inverted S function')
plt.plot(tau, S) #color='g' #label='Source function') 
plt.xscale('log')
plt.xlabel(r'$\tau$',size=18)
plt.ylabel('S',size=18)
plt.legend()
plt.show()
#plt.savefig('invertedSfunction')
#print(S_inv-S)


# In[42]:

#help(np.interp)


# In[110]:

H=np.zeros([51,51])                            #1.17
nk=47
H[0,0:3]=(2,-4,2)
H[1,0:4]=(-4,10,-8,2)

for i in range(nk):
    H[i+2,i:5+i]=(2,-8,12,-8,2)

H[49,47:51]=(2,-8,10,-4)
H[50,48:51]=(2,-4,2)
print(H[40:51,40:51])


# In[111]:

print(np.sum(H,axis=0))


# In[112]:

lam=np.linspace(8e-5, 8e-4, 51)             #1.18
X2=np.zeros(51)
KI=np.matmul(KT,I)

for i in range(51):
    S_reg=np.matmul(np.linalg.inv(KTK+lam[i]*H),np.matmul(KT,I))
    X2[i]=np.sum((S-S_reg)**2/S**2)

#print(lam)
#print(X2)
plt.close("all")    
plt.plot(lam,X2)
plt.xlabel('lambda')
plt.ylabel('chi^2')
plt.show()
plt.savefig('lambdamin')


# In[23]:

np.where(np.min(X2)==X2)


# In[113]:

print(lam[7])


# In[114]:

for i in range(51):
    
    S_reg=np.matmul(np.linalg.inv((np.matmul(KT,K))+lam[7]*H), np.matmul(KT,I))
    X2[i]=np.sum((S-S_reg)**2/S**2)

plt.close("all") 
plt.plot(tau,S_reg, color='b',label='S_regularization')
plt.plot(tau,S,color='g',label='S')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log tau')
plt.ylabel('log intensity')
plt.legend()
plt.savefig('Sregularization')


# In[115]:

#noise = np.random.normal(0,1,100)

# 0 is the mean of the normal distribution you are choosing from
# 1 is the standard deviation of the normal distribution
# 100 is the number of elements you get in array noise

I_1=I*(1+np.random.uniform(-0.01,0.01,21))
I_10=I*(1+np.random.uniform(-0.1,0.1,21))


S_reg_1=np.matmul(np.linalg.inv((np.matmul(KT,K))+8e-4*H), np.matmul(KT,I_1))    #changed for largest lambda
S_reg_10=np.matmul(np.linalg.inv((np.matmul(KT,K))+8e-4*H), np.matmul(KT,I_10))


plt.close("all") 
#plt.plot(y,I_1)
#plt.plot(y,I_10)
plt.plot(tau, S_reg_1, color='b',label='S_reg, 1% noise')
plt.plot(tau, S_reg_10, color='r',label='S_reg, 10% noise')
plt.plot(tau, S,color='g',label='S')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log tau')
plt.ylabel('log I')
plt.legend()
plt.savefig('Snoise')


# In[116]:

def Quadrature():
    mu=np.linspace(0.05,0.95,5)
    h=0.25                                              #h=0.225 normalized = 0.25
    a=np.array([h*14/45,h*64/45,h*24/45,h*64/45,h*14/45])
    J_bode=h*np.sum(a*mu)
    return J_bode,a,mu

J_bode,a,mu = Quadrature()
print(J_bode,a,mu)


# In[117]:

print((0.95-0.05)/4)
#print(tau)


# In[118]:

def lmbd_matrix(tau, mu):
  '''
  Input:
    tau -- 1D optical depth scale
    mu  -- positional angle
  Output:
    matrix for lambda transformation
  '''
  tau_l = tau / mu
  r0 = 0.0
  rn = 1.0

  nd   = len(tau_l)
  dtau = np.roll(tau_l, -1) - tau_l # Don't need dtau[-1]
  dts  = dtau + np.roll(dtau, 1)    # Don't need dts[0]

  a = 2.0 / (dts * np.roll(dtau, 1))
  c = 2.0 / (dts * dtau)
  abc = np.ones(nd, dtype = float)

  # upper boundary condition
  #   I^- = r0 * I^+
  ff0  = (1.0 - r0) / (1.0 + r0)
  abc0 = 1.0 + 2.0 * ff0 / dtau[0]
  c0 = 2.0 / dtau[0]**2

  # lower boundary condition
  #   I^+ = rn * I^-
  ffn  = (1.0 - rn) / (1.0 + rn)
  abcn = 1.0 + 2.0 * ffn / dtau[-2]
  an = 2.0 / dtau[-2]**2

  # elimination
  f     = np.zeros(nd, dtype = float)
  z_mat = np.zeros((nd, nd), dtype = float)
  f[0]  = abc0 / c0
  z_mat[0, 0] = 1.0 / (abc0 + c0)
  for k in range(1, nd - 1):
    f[k] = (abc[k] + a[k] * f[k - 1] / (1.0 + f[k - 1])) / c[k]
    z_mat[0:k + 1, k] = 1.0 / (c[k] * (1.0 + f[k])) * np.append(a[k] * z_mat[0:k, k - 1], 1)

  mat = np.zeros((nd, nd), dtype = float)
  mat[:, -1] = 1.0 / ( abcn + an * (f[-2] / (1.0 + f[-2]))) * np.append(an * z_mat[0:-1, -2], 1)

  # backsubstitute
  for k in range(nd - 2, -1, -1):
    mat[:, k] = 1.0 / (1.0 + f[k]) * mat[:, k + 1] + z_mat[:, k]

  # this routine has been converted from IDL, where col-row order is different,
  # so that this matrix has to be transposed
  return mat.T


# In[125]:

def lmbd(tau):
    Lambda=np.zeros((101,101))
    for i in range (5):
        Lambda+= a[i]*lmbd_matrix(tau,mu[i])
    return Lambda
        

##3.8
tau=np.logspace(-4,6,num=101)
L=lmbd(tau)
S=1-np.exp(-tau)
J=np.matmul(L,S)

'''plt.close("all")
plt.plot(tau,S, color='g')
plt.plot(tau,J)
plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.show()'''


# In[120]:

plt.close("all")
plt.imshow(L)
plt.show()


# In[148]:

#inversion
epsilon=0.1

B=S
S_inversion=np.matmul((np.linalg.inv(np.identity(101)-(1-epsilon)*L)),epsilon*B)

plt.close("all")
plt.plot(tau,B,  label='B')
plt.plot(tau,J, label='J')
plt.plot(tau,S_inversion, label='S')
plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.show()


# In[164]:

#!/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import sys
#IMPORT YOUR LMBD FUNCTION HERE

def lambda_iterate(tau, B, epsilon, plot = False):
  '''
  Input:
    tau -- optical depth scale
    B -- Planck function as a function of tau
    epsilon -- destruction probability per extinction
  Output:
    J -- mean intensity
    S -- source function after iterations
    step -- number of iterations needed
  '''
  # 1st guess
  step = 0
  S = B
  Lambda = lmbd(tau) #NOTE, DONT FORGET TO IMPORT THIS FUNCTION
  if plot:
    J = Lambda.dot(S)
    fig = plt.figure()
    main_ax = fig.add_axes([0.1, 0.09, 0.87, 0.87])
    main_ax.plot(tau, B, '.k', linewidth = 2, markersize = 11, label = 'B')
    main_ax.plot(tau, J, '.-r', linewidth = 0.5, markersize = 3, alpha = 0.75, label = 'J')
    main_ax.plot(tau, S, '.-g', linewidth = 0.5, markersize = 3, alpha = 0.75, label = 'S')

  # iteration
  start_time = time.clock()
  while True:
    S_prev = S
    S = (1.0 - epsilon) * Lambda.dot(S) + epsilon * B
    dev = np.amax( abs(S_prev - S) / S )
    if dev < 1e-2 * epsilon:
      break
    step += 1
    if plot and (step % 10 == 0):
    #if plot and (step  % 10 == 1):
      J = Lambda.dot(S)
      main_ax.plot(tau, J, '.-r', linewidth = 0.5, markersize = 3, alpha = 0.75)
      main_ax.plot(tau, S, '.-g', linewidth = 0.5, markersize = 3, alpha = 0.75)

  end_time = time.clock()
  dt = end_time - start_time
  
  J = Lambda.dot(S)
  if plot:
    main_ax.set_xscale(u'log')
    main_ax.set_yscale(u'log')
    main_ax.set_xlim((7e-5, 1e+2)) # 1.3e+6
    main_ax.set_ylim((1.5e-6, 1.3e0))
    main_ax.set_xlabel('Continuum optical depth')
    main_ax.set_ylabel('Source function, mean intensity')
    main_ax.legend(loc = 'best', fancybox = True)
    plt.savefig('lambda_liter2.pdf', dpi = 600, facecolor = 'w', edgecolor = 'w',
      orientation = 'portrait', format = 'pdf', bboxinches = 'tight')
    #plt.savefig('SN314', dpi = 600, facecolor = 'w', edgecolor = 'w',
      #orientation = 'portrait', format = 'pdf', bboxinches = 'tight')
    plt.show()
    plt.close(fig)
  return (J, S, step, dt)


# In[165]:

#epsilon=0.01                                               #3.14
#B=1-np.exp(-tau)


#lambda_iterate(tau, B, epsilon, plot= True)


# In[166]:

epsilon=0.009                                             #3.13

lambda_iterate(tau, B, epsilon, plot= True)


# In[178]:

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
import sys
#IMPORT YOUR LMBD FUNCTION HERE

def lambda_ap_iterate(tau, B, epsilon, plot = False):
  '''
  Input:
    tau -- optical depth scale
    B -- Planck function as a function of tau
    epsilon -- destruction probability per extinction
  Output:
    J -- mean intensity
    S -- source function after iterations
    step -- number of iterations needed
  '''
  # 1st guess
  step = 0
  
  Lambda = lmbd(tau) #NOTE, DONT FORGET TO IMPORT THIS FUNCTION
   
  S = np.matmul(Lstar,B)
  if plot:
    J = Lambda.dot(S)
    fig = plt.figure()
    main_ax = fig.add_axes([0.1, 0.09, 0.87, 0.87])
    main_ax.plot(tau, B, '.k', linewidth = 2, markersize = 11, label = 'B')
    main_ax.plot(tau, J, '.-r', linewidth = 0.5, markersize = 3, alpha = 0.75, label = 'J')
    main_ax.plot(tau, S, '.-g', linewidth = 0.5, markersize = 3, alpha = 0.75, label = 'S')

  # iteration
  start_time = time.clock()
  while True:
    S_prev = S
    S =np.matmul(np.linalg.inv(np.identity(101)-(1-epsilon)*Lstar), ((1.0 - epsilon) * (Lambda-Lstar)).dot(S) + epsilon * B)
    dev = np.amax( abs(S_prev - S) / S )
    if dev < 1e-2 * epsilon:
      break
    step += 1
    if plot and (step % 10 == 0):
    #if plot and (step  % 10 == 1):
      J = Lambda.dot(S)
      main_ax.plot(tau, J, '.-r', linewidth = 0.5, markersize = 3, alpha = 0.75)
      main_ax.plot(tau, S, '.-g', linewidth = 0.5, markersize = 3, alpha = 0.75)

  end_time = time.clock()
  dt = end_time - start_time
  
  J = Lambda.dot(S)
  if plot:
    main_ax.set_xscale(u'log')
    main_ax.set_yscale(u'log')
    main_ax.set_xlim((7e-5, 1e+2)) # 1.3e+6
    main_ax.set_ylim((1.5e-6, 1.3e0))
    main_ax.set_xlabel('Continuum optical depth')
    main_ax.set_ylabel('Source function, mean intensity')
    main_ax.legend(loc = 'best', fancybox = True)
    plt.savefig('lambda_liter2.pdf', dpi = 600, facecolor = 'w', edgecolor = 'w',
      orientation = 'portrait', format = 'pdf', bboxinches = 'tight')
    #plt.savefig('SN314', dpi = 600, facecolor = 'w', edgecolor = 'w',
      #orientation = 'portrait', format = 'pdf', bboxinches = 'tight')
    plt.show()
    plt.close(fig)
  return (J, S, step, dt)


# In[181]:

def Lambstar():
    Lstar=np.zeros((101,101))
    for i in range (101):
        Lstar[i][i]=L[i][i] 
    return Lstar
        
Lstar=Lambstar()
epsilon=0.001                                               #3.5
B=1-np.exp(-tau)


lambda_ap_iterate(tau, B, epsilon, plot= True)


# In[185]:

epsilon=np.linspace(0.001,0.9,20)
yL=np.zeros(20)
yLacc=np.zeros(20)

for i in range(20):
    JL,SL,syL,dtL=lambda_iterate(tau, B, epsilon[i], plot= False)
    JLacc, SLacc, syLacc,dtLacc=lambda_ap_iterate(tau, B, epsilon[i], plot= False)
    yL[i] += syL
    yLacc[i] += syLacc

plt.close("all")
plt.plot(epsilon,yL, label='Lambda', color='r')
plt.plot(epsilon,yLacc, label='ALI', color='b')
plt.xlabel(r'$\tau$',size=18)
plt.ylabel('Number of Steps',size=18)
plt.legend(fontsize=20)
plt.show()


# In[ ]:



