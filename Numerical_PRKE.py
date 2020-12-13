# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:33:42 2020

@author: Donald
"""

import numpy as np
import matplotlib.pyplot as plt


# decay constant for DNP groups, 1/s
lambda_i = np.array([0.0124 , 0.0305  , 0.111   , 0.301   , 1.14    , 3.01   ])
# delayed neutron fractions per DNP froups
beta_i   = np.array([0.00021, 0.00142 , 0.00127 , 0.00258 , 0.00075 , 0.00027])
beta_tot = sum(beta_i)
n_dnp = 6
MGT = 2e-4
rho = beta_tot/11
# initial values for n and c
X0 = np.array([1., beta_i[0]/(lambda_i[0]*MGT), beta_i[1]/(lambda_i[1]*MGT)\
               , beta_i[2]/(lambda_i[2]*MGT), beta_i[3]/(lambda_i[3]*MGT)\
               , beta_i[4]/(lambda_i[4]*MGT), beta_i[5]/(lambda_i[5]*MGT)])
# end time for simulation
Tend = 30
# number of time steps
n_steps = 30
# time step size 
dt = Tend / n_steps


# identity matrix
I = np.eye(n_dnp+1)
A = np.zeros((n_dnp+1,n_dnp+1))
A[0,0] = (rho-beta_tot)/MGT
for i in range(n_dnp):
    A[0,i+1] = lambda_i[i]
    A[i+1,0] = beta_i[i]/MGT
    A[i+1,i+1] = -lambda_i[i]
# final form of the linear system matrix
M = I - dt*A

# storage place for plotting solution later
sol = np.zeros((7,n_steps+1))
sol[:,0] = X0

# loop through time steps
for i in range(n_steps):
        # end time steps values
        X1 = np.linalg.solve(M,X0)
        # store
        sol[:,i+1] = X1
        # X1 becomes initial value for next time step
        X0 = np.copy(X1)    

# Make plot with 1 neutron group
beta_ = 650e-5
lambda_ = 0.3
# initial values for n and c
X0_1 = np.array([1., beta_/(lambda_*MGT)])
# end time for simulation
Tend_1 = 30
# number of time steps
n_steps_1 = 3000
# time step size 
dt_1 = Tend_1 / n_steps_1


# identity matrix
I_1 = np.eye(2)
A_1 = np.array([[(rho-beta_)/MGT, lambda_], [beta_/MGT, -lambda_]])
# final form of the linear system matrix
M_1 = I_1 - dt_1*A_1

# storage place for plotting solution later
sol_1 = np.zeros((2,n_steps_1+1))
sol_1[:,0] = X0_1

# loop through time steps
for i in range(n_steps_1):
        # end time steps values
        X1_1 = np.linalg.solve(M_1,X0_1)
        # store
        sol_1[:,i+1] = X1_1
        # X1 becomes initial value for next time step
        X0_1 = np.copy(X1_1)
        
time_30 = np.linspace(0,Tend,n_steps+1)
time_3000 = np.linspace(0,Tend,n_steps_1+1)
plt.figure(2, dpi=250)
plt.plot(time_30,sol[0,:]/sol[0,0],label='n with 6 groups')
plt.plot(time_30,sol[1,:]/sol[1,0],label='c with 6 groups')
plt.plot(time_3000,sol_1[0,:]/sol_1[0,0],label='n with 1 groups')
plt.plot(time_3000,sol_1[1,:]/sol_1[1,0],label='c with 1 groups')
plt.xlabel('time [s]')
plt.ylabel('Neutron population')
plt.title('Comparing ^ and 1 neutron group solution')
plt.grid()
plt.legend()
plt.show()


