import math
import scipy
import scipy.optimize as optimize
from scipy.optimize import minimize

import numpy as np
import pandas as pd
from numpy import linalg as LNG 
import HMM

N=2

# forward-algo
def alpha_scaled_opt(a,b, o, p, cut):
    seq_length = len(o)*cut
    
    c = np.zeros(int(seq_length/cut)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, int(seq_length/cut)))
    alpha[:, 0] = b[:, o[0][0],o[0][1]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    
    
    

    for t in range(1, int(seq_length/cut)):   
        
        for i in range(0, N):
            alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt(a,b, o, scaling_factors, cut):
    seq_length = len(o)*cut
    
    beta = np.zeros((N, int(seq_length/cut)))
    
    length = int(seq_length/cut)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(int(seq_length/cut)-2,-1,-1):             
        for i in range(0, N):             
            for l in range(0, N):
                beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma(alpha, beta):

        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi( a, b, o, alpha, beta):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        
        denom = 0
        for i in range(0, N):
            for j in range(0, N):
                denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]
                
        
        for i in range(0, N):
            for j in range(0, N):
                ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]) / denom
    
    return ksi









def opt_params(MU,RR,seq, n_st, lmbd_initial, epsilon, n_neanderthals, bnds, ind_num, cut, n_step, nd_portion):
    lmbd_0 = lmbd_initial
    d = MU * cut
    seq_length = cut*len(seq[ind_num])
    if n_neanderthals != 0:

        for iii in range(n_step):           
       
            A_0 = HMM.initA(cut,RR, lmbd_0[0]/d, lmbd_0[3])###ошибка в порядке аргументов
            b_0 = HMM.initB(MU, cut, lmbd_0[0:3], n_st+1)
            P=[0.95, 0.05]
        
            alpha, sc_factors = alpha_scaled_opt(A_0,b_0, seq [ind_num], P, cut)
            beta = beta_scaled_opt(A_0, b_0, seq [ind_num], sc_factors, cut)        
            gamma = def_gamma(alpha, beta)
            ks = def_ksi( A_0, b_0, seq [ind_num], alpha, beta)
            coeff_a = np.zeros((N, N))

            for i in range(N):
                for j in range(N):                
                    for t in range(int(seq_length/cut)-1):
                        coeff_a[i, j] += ks[i, j, t] 

            def multi_Q(lmbd):

                lmbd=np.array(lmbd)
                Q=0
                A = HMM.initA(cut,RR,lmbd[0]/d, lmbd[3])
                b = HMM.initB(MU, cut, lmbd[0:3],  n_st+1)


                for ii in range(N):
                    for jj in range(N):
                        Q += math.log(A[ii,jj]) * coeff_a[ii,jj]

                for ii in range(N):
                    for t in range(int(seq_length/cut)):
                        Q += math.log ( b[ii, seq [ind_num][t][0], seq [ind_num][t][1]]) * gamma[ii,t]

                return -Q

            def gradient_respecting_bounds(bounds, fun, eps=epsilon):
                """bounds: list of tuples (lower, upper)"""
                def gradient(x):
                    fx = fun(x)
                    grad = np.zeros(len(x))
                    for k in range(len(x)):
                        d = np.zeros(len(x))
                        d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
                        grad[k] = (fun(x + d) - fx) / d[k]
                    return grad
                return gradient

            opt_result = scipy.optimize.minimize(multi_Q, 
                                                 lmbd_0,                                                 
                                                 bounds=bnds,                                                  
                                                 method="L-BFGS-B" )


            if LNG.norm(lmbd_0-opt_result.x) < epsilon:
                break
            else:
                lmbd_0[0:4] = opt_result.x[0:4]
                

    return lmbd_0
