import HMM
import numpy as np
from numpy import linalg as LNG 
import math
import scipy
import scipy.optimize as optimize
from scipy.optimize import minimize
import math
from numpy import linalg as LNG 

import useful as usfl

N=2
# forward-algo
def alpha_scaled_opt(a,b, o, p, cut):

    
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    N=2
    alpha = np.zeros((N, len(o)))
    alpha[:, 0] = b[:, o[0][0],o[0][1]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])  
    
    

    for t in range(1, len(o)):   
        
        for i in range(0, N):
            alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt(a,b, o, scaling_factors, cut):

    N=2
    
    beta = np.zeros((N,len(o)))
    beta[:, len(o) - 1] = np.ones(N)*scaling_factors[len(o)-1] 
    

    for t in range(len(o)-2,-1,-1):             
        for i in range(0, N):             
            for l in range(0, N):
                beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma(alpha, beta):
    N=2

        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi( a, b, o, alpha, beta):
    N=2
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


    
    




def new_lambda_i_coal2(O, Gamma):
    nom, denom = 0, 0
    i=0
    
    for o in O:
        gamma=Gamma[i]
        i+=1
        
        for t in range(0, len(o), 1):   
            nom += o[t, 1] *  gamma[1, t]
            denom += gamma[1,t]
    return nom/denom
    
    
    
def new_lambda_n2(O, Gamma):
    nom, denom = 0, 0
    i=0
    for o in O:    
        gamma=Gamma[i]
        i+=1        

        for t in range(0, len(o), 1):   
            nom += o[t, 1] *  gamma[0, t]  + o[t, 0]  * gamma[1, t]
            denom += gamma[0, t] + gamma[1, t]    
    
    return nom/ denom
          
          
def new_lambda_af2(O, Gamma):
    nom, denom = 0, 0
    i=0
    
    for o in O:    
        gamma=Gamma[i]
        i+=1   
        for t in range(0, len(o), 1):   
            nom += o[t, 0] * gamma[0, t]
            denom += gamma[0, t]  
    
    return nom/ denom
          











    


def new_a_ij2(Gamma, Ksi):
    a = np.zeros((5, 5 ))
    nom = 0
    denom = 0
    
    
    for i in range(N):
        for j in range(N):
            for ksi, gamma in zip(Ksi, Gamma):
                for t in range(len(gamma[0])-1):
                    nom += ksi[i,j,t]
                    denom += gamma[i, t]
            a[i,j] = nom/denom
            nom=0
            denom=0
            
    return a








def E_step2(cut,  p, O, n_states, mu,rr, lambda_old):

   
    b = HMM.initB(mu, cut, lambda_old[0:3], n_states+1)
    a = HMM.initA(cut,rr, lambda_old[4]/(mu*cut), lambda_old[3])
    
    GAMMA=[]
    KSI=[]
    for o in O:
        alpha, sc_factors = alpha_scaled_opt(a,b, o, p, cut)
        beta = beta_scaled_opt(a, b, o, sc_factors, cut)    
        gamma = def_gamma(alpha, beta)    
        ksi = def_ksi( a, b, o, alpha, beta)
        GAMMA.append(gamma)
        KSI.append(ksi)
        
    

    
    a_new=new_a_ij2(GAMMA, KSI)
    
    lambda_new4=(a_new[0][1]+a_new[1][0])*mu/rr
    lambda_new3=a_new[0][1]/(a_new[0][1]+a_new[1][0])
    


    return new_lambda_i_coal2(O, GAMMA), new_lambda_n2(O, GAMMA), new_lambda_af2(O, GAMMA), lambda_new3, lambda_new4






def EM_algorithm2(p,o, n_states, mu, rr, lambda_0, epsilon, cut ):
    
    lmbd = np.array(lambda_0)
    
    em_steps = 0

    for i in range(100):


        lmbd_new = np.array(E_step2(cut,  p, o, n_states, mu, rr, lmbd))
        if lmbd_new[4]>0.1:
            print('Oops. something went wrong with sample')
            return(lambda_0)
        em_steps += 1
        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new
        

#    print('Число шагов в EM -алгоритме', em_steps )
    return lmbd_new










def E_step3(cut,  p, O, n_states, mu,rr, lambda_old):

   
    b = HMM.initB(mu, cut, lambda_old[0:3], n_states+1)
    a = HMM.initA(cut,rr, lambda_old[4]/(mu*cut), lambda_old[3])
    
    GAMMA=[]
    for o in O:
        alpha, sc_factors = alpha_scaled_opt(a,b, o, p, cut)
        beta = beta_scaled_opt(a, b, o, sc_factors, cut)    
        gamma = def_gamma(alpha, beta)    

        GAMMA.append(gamma)       

    return new_lambda_i_coal2(O, GAMMA), new_lambda_n2(O, GAMMA), new_lambda_af2(O, GAMMA), lambda_old[3], lambda_old[4]











def EM_algorithm3(p,o, n_states, mu, rr, lambda_0, epsilon, cut ):
    
    lmbd = np.array(lambda_0)
    
    em_steps = 0

    for i in range(100):


        lmbd_new = np.array(E_step3(cut,  p, o, n_states, mu, rr, lmbd))
        em_steps += 1
        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new


        
#    print('Число шагов в EM -алгоритме', em_steps )
    return lmbd_new
    
    
    
    
    
    
    
    
    
    

N=2   
a_gaps=np.identity(N)
b_gaps=np.ones((N,1))    
    
#EM-common+gaps
 # forward-algo
def alpha_scaled_opt_gaps(a,b, o, p, gaps):   
    c = np.zeros(len(o)) #scaling factors, которые как раз позволяют не обнулиться
    
    alpha = np.zeros((N, len(o)))
    alpha[:, 0] = b[:, o[0][0],o[0][1]] * p
    
    c[0] = 1 / sum(alpha[:, 0])
    alpha[:, 0] = alpha[:, 0] / sum(alpha[:, 0])
    for t in range(1, len(o)):   
        if usfl.point_in_set(t, gaps)==True:

            for i in range(0, N):
                
                alpha[i, t] = np.dot(alpha[:, t-1],a_gaps[:,i]) * b_gaps[i][0]     

        else:
        
            for i in range(0, N):
            
                alpha[i, t] = np.dot(alpha[:, t-1],a[:,i]) * b[i, o[t][0],o[t][1]] 
            
        c[t] = 1 / sum(alpha[:,t]) #сохраняем множители        
        alpha[:, t] = alpha[:, t] / sum(alpha[:,t])     
        
    return alpha, c

# Backward procedure. Scaled case.
def beta_scaled_opt_gaps(a,b, o, scaling_factors, gaps):
    
    beta = np.zeros((N, len(o)))
    
    length = len(o)
    beta[:, length - 1] = np.ones(N)*scaling_factors[length-1] 
    

    for t in range(len(o)-2,-1,-1):        
        if usfl.point_in_set(t, gaps)==True:
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a_gaps[i, l] * b_gaps[l][0] * beta[l, t+1]       
        else:     
            for i in range(0, N):             
                for l in range(0, N):
                    beta[i, t] += a[i, l] * b[l, o[t+1][0],o[t+1][1]] * beta[l, t+1]
                
        beta[:, t] = beta[:, t] * scaling_factors[t]

    return beta 

# gamma matrix
def def_gamma_gaps(alpha, beta, gaps):
        
    gamma = np.zeros((N,len(alpha[0])))
    for m in range(0, len(alpha[0])):
        denom = sum(alpha[:, m]*beta[:,m])
        
        for i in range(0,N):
            gamma[i, m] = (alpha[i, m] * beta[i, m]) / denom
    return gamma


# ksi[i, j, t]
def def_ksi_gaps( a, b, o, alpha, beta, gaps):
    
    M = len(o)
    ksi = np.zeros((N, N, M-1))
    
    for t in range(0, M-1):
        if point_in_set(t, gaps)==True:
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a_gaps[i, j] * b_gaps[j][0] * beta[j, t+1]        
            for i in range(0, N):
                for j in range(0, N):       
                    ksi[i, j, t] = (alpha[i, t]*a_gaps[i, j]*b_gaps[j,0]* beta[j, t+1]) / denom
        else:
        
            denom = 0
            for i in range(0, N):
                for j in range(0, N):
                    denom += alpha[i, t] * a[i, j] * b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]
                
        
            for i in range(0, N):
                for j in range(0, N):
                    ksi[i, j, t] = (alpha[i, t]*a[i, j]*b[j, o[t+1][0],o[t+1][1]] * beta[j, t+1]) / denom
    
    return ksi

    
def new_lambda_i_gaps(O, Gamma, gaps):
    nom, denom = 0, 0
    i=0
    
    for o in O:
        gamma=Gamma[i]
        i+=1
        
        for t in range(0, len(o), 1):
            if usfl.point_in_set(t, gaps)==False:   
                nom += o[t, 1] *  gamma[1, t]
                denom += gamma[1,t]
    return nom/denom
   
def new_lambda_n_gaps(O, Gamma, gaps):
    nom, denom = 0, 0
    i=0
    for o in O:    
        gamma=Gamma[i]
        i+=1        

        for t in range(0, len(o), 1):   
            if usfl.point_in_set(t, gaps)==False:
                nom += o[t, 1] *  gamma[0, t]  + o[t, 0]  * gamma[1, t]
                denom += gamma[0, t] + gamma[1, t]    
    
    return nom/ denom
          
def new_lambda_af_gaps(O, Gamma, gaps):
    nom, denom = 0, 0
    i=0
    
    for o in O:    
        gamma=Gamma[i]
        i+=1   
        for t in range(0, len(o), 1):   
            if usfl.point_in_set(t, gaps)==False:
                nom += o[t, 0] * gamma[0, t]
                denom += gamma[0, t]  
    
    return nom/ denom
          
    
    
    
    
def E_step_gaps(cut,  p, O, n_states, mu,rr, lambda_old, gaps):

   
    b = HMM.initB(mu, cut, lambda_old[0:3], n_states+1)
    a = HMM.initA(cut,rr, lambda_old[4]/(mu*cut), lambda_old[3])
    
    GAMMA=[]
    for o in O:
        alpha, sc_factors = alpha_scaled_opt_gaps(a,b, o, p, gaps)
        beta = beta_scaled_opt_gaps(a, b, o, sc_factors, gaps)    
        gamma = def_gamma_gaps(alpha, beta, gaps)    

        GAMMA.append(gamma)       

    return new_lambda_i_gaps(O, GAMMA, gaps), new_lambda_n_gaps(O, GAMMA, gaps), new_lambda_af_gaps(O, GAMMA, gaps), lambda_old[3], lambda_old[4]   
    
    
    
    
def EM_algorithm_gaps(p,o_mas, n_states, mu, rr, lambda_0, epsilon, cut, n_em_steps, gaps ):
    
    lmbd = np.array(lambda_0)
    
    em_steps = 0

    for i in range(n_em_steps):


        lmbd_new = np.array(E_step_gaps(cut,  p, o_mas, n_states, mu, rr, lmbd,gaps))
        if lmbd_new[4]>0.1:
            print('Oops. something went wrong with sample')
            return(lambda_0)
        em_steps += 1
        print(lmbd_new-lmbd)
        if LNG.norm(lmbd_new-lmbd) < epsilon:
            break
        lmbd = lmbd_new
        

#    print('Число шагов в EM -алгоритме', em_steps )
    return lmbd_new
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
