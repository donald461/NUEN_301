# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:33:42 2020

@author: Donald
"""

import numpy as np
import matplotlib.pyplot as plt

def numericalSolution(rho_init, alex_variables):
    printValue = alex_variables *10
    print("Alex variables: "+ str(printValue))
    #rho_init = 
    
    # decay constant
    lambda_ = 0.3
    # DNP fraction
    beta_ = 650e-5
    # mean generation time [s]
    MGT = 2e-4
    
    # diffine intial and final reactivity
    rho_init = 0
    rho = 0.5*beta_
    # Define the source term
    S = 0
    
    # initial values for n and c
    X0 = np.array([1., beta_/(lambda_*MGT)])
    # end time for simulation
    Tend = 30
    # number of time steps
    n_steps = 30
    # time step size 
    dt = Tend / n_steps
    
    
    # identity matrix
    I = np.eye(2)
    A = np.array([[(rho-beta_)/MGT, lambda_], [beta_/MGT, -lambda_]])
    # final form of the linear system matrix
    M = I - dt*A
    
    # storage place for plotting solution later
    sol = np.zeros((2,n_steps+1))
    sol[:,0] = X0
    
    # loop through time steps
    for i in range(n_steps):
            # end time steps values
            X1 = np.linalg.solve(M,X0)
            # store
            sol[:,i+1] = X1
            # X1 becomes initial value for next time step
            X0 = np.copy(X1)
         
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
    plt.close('all')
    plt.figure(2, dpi=250)
    plt.plot(time_30,sol[0,:]/sol[0,0],label='n dt=1')
    plt.plot(time_30,sol[1,:]/sol[1,0],label='cd t=1')
    plt.plot(time_3000,sol_1[0,:]/sol_1[0,0],label='n dt=0.01')
    plt.plot(time_3000,sol_1[1,:]/sol_1[1,0],label='cd t=0.01')
    plt.grid()
    plt.legend()
    plt.show()
    
if __name__ == '__main__':
    # here is solution for X
    print("===========================================================")
    print("Solution for Homework X: Question Y: Part A.i")
    var = 30
    aVar = 7
    numericalSolution_1(var, aVar)
    
    
    #print(table)
    # here is solution for y
    print("\n===========================================================")
    print("Solution for Homework X: Question Y: Part A.i")
    var = 20
    aVar = 20
    numericalSolution(var, aVar)
    
    
    


