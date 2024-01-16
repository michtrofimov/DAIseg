import argparse

import random
import numpy as np
import random
import sklearn

import math



#Ti: Introgression of Nd
def initA(Ti,cut,a) -> np.array:
    A = np.zeros((2,2))
    
    A[0][1]=Ti*RR*cut*a
    A[0][0]=1-A[0][1]
 
    A[1][0]=Ti*RR*cut*(1-a)
    A[1][1]=1-A[1][0]
    
    return A

#Ti: Introgression of Nd
#Taf: Time out of Africa
#Tn: Time of Split between Nd and Sapiens

def initB(m,cut, lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(2,n_st,n_st))
    meani = lmbd[0]
    meann = lmbd[1]
    meanaf = lmbd[2]
    
    Pi = np.empty(n_st)
    Paf=np.empty(n_st)
    Pn=np.empty(n_st)

    
    Pi[0]=np.exp(-meani)
    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    
    sumi=0
    sumaf=0
    sumn=0
    
    for i in range(1,n_st):
        Pi[i]=Pi[i-1]*meani/i
        Paf[i]=Paf[i-1]*meanaf/i
        Pn[i]=Pn[i-1]*meann/i
        
        sumi=sumi+Pi[i]
        sumaf=sumaf+Paf[i]
        sumn=sumn+Pn[i]

    Pi[0]=1-sumi
    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st): 
        for j in range(n_st):
            B[0][i][j]=Paf[i]*Pn[j]
            B[1][i][j]=Pn[i]*Pi[j]
                  
    return B




# forward-algo
def alpha_scaled_opt(a,b, o, p, cut):
    
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



def viterbi(V, initial_distribution, a, b):
    
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * b[:, V[0][0],V[0][1]])
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(a[:, j]) + np.log(b[j, V[t][0], V[t][1]])
 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)

    return result




def get_HMM_tracts(seq):
    migrating_tracts = []
    for i in range(N):
        migrating_tracts.append([])
    start=0
    for i in range(1,len(seq)):
        if seq[i]!=seq[i-1]:
            migrating_tracts[seq[i-1]].append([start,i-1])
            start=i
    migrating_tracts[seq[len(seq)-1]].append([start,len(seq)-1])
    return migrating_tracts


parser = argparse.ArgumentParser(description='DAIseg')
parser.add_argument('--obs1', type= str, help='File with observations with respect to Outgroup')
parser.add_argument('--obs2', type= str, help='File with observations with respect to Archaic')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o', type= str, help = 'Name of output file' )

args = parser.parse_args()



f = open(args.HMM_par, 'r')
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())
seq_length = int(f.readline())
t_n = float(f.readline())
t_ooa = float(f.readline())
t_i = float(f.readline())
f.close()

print(t_n, t_ooa, t_i)


N = 2 # number of hidden states

seq1, seq2 = [], []
with open(args.obs1, 'r') as f1, open(args.obs2, 'r') as f2:
    for line1, line2 in zip(f1, f2):
        row = line1.replace('\n','').split(' ')
        row = [int(i) for i in row]
        seq1.append(row)
        row = line2.replace('\n','').split(' ')
        row = [int(i) for i in row]
        seq2.append(row)
        
seq1=np.array(seq1)
seq1 = np.transpose(seq1)
seq2=np.array(seq2)
seq2 = np.transpose(seq2)

seq=[]
for i in range(len(seq1)):
    seq.append(np.column_stack((seq1[i], seq2[i])))   

seq=np.array(seq)

print(np.max(seq))


d = 2 * MU * L
n_st = np.max(seq)

lmbd_opt = d* np.array([t_n, t_ooa, t_i])
a = 0.05
A = initA(L,lmbd_opt[0]/d, a)
B = initB(MU, L, lmbd_opt[0:3],   n_st+1)
P=[0.95, 0.05]

tracts_HMM_result = []
for idx in range(0, len(seq)):
    tracts_HMM =  get_HMM_tracts(viterbi(seq [idx], P, A, B))
    for k in range(N):
        for j in range(len(tracts_HMM[k])):
            tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]
            tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]

    tracts_HMM_result.append(tracts_HMM)


with open(args.o, "w") as f:
    for i in tracts_HMM_result:
        f.write(str(i)+'\n')


