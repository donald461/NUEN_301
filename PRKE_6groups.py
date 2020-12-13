
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 10:30:03 2020

@author: Donald
"""
import numpy as np
import matplotlib.pyplot as plt

print('Case number = 7')
# Value for external source
S = 2
# decay constant
lambda_ = 0.3
# DNP fraction
beta_ = 650e-5
# mean generation time [s] (big lambda)
MGT = 2e-4
# number of delayed neutron precursor groups
nb_DNP = 1
# given initial neutron population
n_init = 1
# define initial precursor
c_init = beta_/(lambda_ * MGT)
print('n(0) = ', n_init)
print('c(0) = ', c_init)

# Define initial and final rho values
rho_init = -beta_/9
rho = 0
# Write down quadratic and give s1, s2
disc = (MGT*lambda_+beta_-rho)**2 + 4*rho*lambda_*MGT
s1 = (rho-beta_-lambda_*MGT + np.sqrt(disc)) / (2*MGT)
s2 = (rho-beta_-lambda_*MGT - np.sqrt(disc)) / (2*MGT)
print('rooot s1: ', s1)
print('rooot s2: ', s2)

# Check if particular solution is needed and if so solve for it
if S == 0:
    n_part = 0
    c_part = 0
    print('Neutron population particular solution: ', n_part)
if S != 0:
    # check if final reactivity is 0 if not constant solutions if so linear solitions
    if rho != 0:
        n_part = -S * MGT / rho
        c_part = beta_ * (lambda_ * MGT) * n_part
        print('Neutron population particular solution: ', n_part)
    else:
        a = (S * lambda_) / ((-beta_ / MGT) + lambda_)
        b = a - S
        d = b / -lambda_  
        print('Neutron population particular solution = ', a, 't')

# solve for the coefficients A1 and A2
mat = np.array([[1,1],[lambda_/(lambda_+s1),lambda_/(lambda_+s2)]])
rhs = np.ones(2)
coef = np.linalg.solve(mat, rhs)
# difine the A_1 and A_2 coeffeicints
A_1 = coef[0]
A_2 = coef[1]
# difine the B_1 and B_2 coeffeicints
B_1 = beta_/MGT * (coef[0]/(lambda_*s1))
B_2 = beta_/MGT * (coef[1]/(lambda_*s2))
print('Amplitude A1 = ', A_1)
print('Amplitude A2 = ', A_2)
print('Amplitude B1 = ', B_1)
print('Amplitude B2 = ', B_2)

# Caluculate the initial promt jump
P_J = (rho_init - beta_) / (rho - beta_)
print('Value of prompt jump = ', P_J)

# Calcualte n(5) and n(30)
if S != 0:
    if rho != 0:
        n_5 = A_1*np.exp(s1*5) + A_2*np.exp(s2*5) + n_part
        n_30 = A_1*np.exp(s1*30) + A_2*np.exp(s2*30) + n_part
        n_500 = A_1*np.exp(s1*500) + A_2*np.exp(s2*500) + n_part
    else:
        n_5 =  A_1*np.exp(s1*5) + A_2*np.exp(s2*5) + a * 5
        n_30 = A_1*np.exp(s1*30) + A_2*np.exp(s2*30) + a * 30
        n_500 =  A_1*np.exp(s1*500) + A_2*np.exp(s2*500) + a * 500
else:
    n_5 =  A_1*np.exp(s1*5) + A_2*np.exp(s2*5)
    n_30 = A_1*np.exp(s1*30) + A_2*np.exp(s2*30)
    n_500 =  A_1*np.exp(s1*500) + A_2*np.exp(s2*500)  
    
print('n(5 sec) = ', n_5)
print('n(30 sec) = ', n_30)
print('n(500 sec) = ', n_500)

# create function for ploting
if S != 0:
    if rho != 0:
        n = lambda time: A_1*np.exp(s1*time) + A_2*np.exp(s2*time) + n_part
        c = lambda time: B_1*np.exp(s1*time) + B_2*np.exp(s2*time) + c_part
    else:
        n = lambda time: A_1*np.exp(s1*time) + A_2*np.exp(s2*time) + a * time
        c = lambda time: B_1*np.exp(s1*time) + B_2*np.exp(s2*time) + b * time + d
else:
    n = lambda time: A_1*np.exp(s1*time) + A_2*np.exp(s2*time)
    c = lambda time: B_1*np.exp(s1*time) + B_2*np.exp(s2*time)
    
t = np.linspace(0,500,10000)
#plt.close('all')
plt.figure(dpi=250)
plt.semilogy(t,n(t)/n(0),label='n')
plt.semilogy(t,c(t)/c(0),label='c')
plt.grid()
plt.legend()
plt.show()