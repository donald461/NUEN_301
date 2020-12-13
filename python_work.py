# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:21:54 2020

@author: Donald
"""

import numpy as np

import matplotlib.pyplot as plt

 
# decay constant for DNP groups, 1/s
lambda_6 = np.array([0.0124 , 0.0305  , 0.111   , 0.301   , 1.14    , 3.01   ])
# delayed neutron fractions per DNP froups
beta_6   = np.array([0.00021, 0.00142 , 0.00127 , 0.00258 , 0.00075 , 0.00027])
# mean generation time, s
MGT = 1e-3
# number of delayed neutron presursor groups
nb_DNP = len(beta_6) 

def inhour(s, beta_, lama_, MGT):
    out = MGT * s
    for i in range(len(beta_)):
        out += beta_[i]*s/(s+lambda_6[i])
    return out

# creating arrays for s values in chunks to avoid division by 0 if s=-lambda[i]
eps = 1e-8
npts = 10000
# empty list 
s = []
# first portion, from -Big to -lambda(end)
s.append( np.linspace(-7, -lambda_6[-1]-eps,npts))
# other portions, from -lambda(k) to -lambda(k-1)
for k in range(nb_DNP-1,0,-1):
    s.append(np.linspace(-lambda_6[k]+eps,-lambda_6[k-1]-eps,npts))
# last portion, from -lambds(1st) to +Big
s.append(np.linspace(-lambda_6[0]+eps,3,npts))

plt.close("all")
plt.figure()
for k in range(nb_DNP):
        plt.plot(s[k],inhour(s[k],beta_6,lambda_6,MGT))
# set axis limits
plt.axis([s[0][0],s[-1][-1],-0.1,0.1])
plt.grid()
plt.show()


    