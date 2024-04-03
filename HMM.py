import numpy as np
import pandas as pd

import math


#number of states
N=2 


#Ti: Introgression of Nd
def initA(cut,RR,Ti, a) -> np.array:
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


#Ti: Introgression of Nd
#Taf: Time out of Africa
#Tn: Time of Split between Nd and Sapiens

def initBwN(m,cut, lmbd, n_st) -> np.array: 
    
    B = np.empty(shape=(2,n_st))
    meann = lmbd[1]
    meanaf = lmbd[2]


    Paf=np.empty(n_st)
    Pn=np.empty(n_st)

    Paf[0]=np.exp(-meanaf)
    Pn[0]=np.exp(-meann)
    

    sumaf=0
    sumn=0
    
    for i in range(1,n_st):
        Paf[i]=Paf[i-1]*meanaf/i
        Pn[i]=Pn[i-1]*meann/i
        
        sumaf=sumaf+Paf[i]
        sumn=sumn+Pn[i]

    Paf[0]=1-sumaf
    Pn[0]=1-sumn
    
    for i in range(n_st): 
        B[0][i]=Paf[i]
        B[1][i]=Pn[i]
               

    return B



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

def viterbinD(V, initial_distribution, a, b):
    
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * b[:, V[0]])
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(a[:, j]) + np.log(b[j, V[t]])
 
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



#something for classification report

def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]

        ri += 1

    return ranges



def ne_set_interval(set_intervals):
    
    ne_set = []
    for j in range(len(set_intervals)):
        
        if j==0:
            ne_set.append([0,set_intervals[j][0]-1])
        if j== len(set_intervals)-1:
            ne_set.append([set_intervals[j][1]+1, seq_length-1])
        if j!=0 and j!= len(set_intervals)-1:
            ne_set.append([set_intervals[j-1][1]+1, set_intervals[j][0]-1])
    return ne_set

def len_tracts(set_intervals):
    if len(set_intervals)==0:
        return 0
    else:
        s=0
        for j in range(len(set_intervals)):
            s+= set_intervals[j][1]-set_intervals[j][0]+1
        return s
    
def confusion_mtrx(real, res_HMM):
    conf_matrix = np.zeros((2,2))
    for i in range(N):
        for j in range(N):
            conf_matrix[i,j] =int(len_tracts(intersections(real[i], res_HMM[j])))
            

    return conf_matrix

def classification_rpt(conf_matrix):
    clas_report = {}
    for i in range(N):
        dd={}
        dd['precision'] = round(conf_matrix[i,i]/sum(conf_matrix[:,i]),7)
        dd['recall'] = round(conf_matrix[i,i]/sum(conf_matrix[i,:]),7)
        dd['f1-score'] = round(2*dd['recall']*dd['precision']/(dd['recall']+dd['precision']), 7)
        clas_report[str(i)] = dd
    return clas_report




def createDataFrame(seq_length,RR,MU,cut,lmbd_opt, n_st, seq, n_neanderthal,  n_ref_pop, n_obs_seq, nd_true_tracts, nd_portion, N_ND):
    
    
    
    # return European tracts with input=Neanderthal tracts
    def tracts_eu(tr_nd, seq_length):
        result = []

        if tr_nd[0][0] > 0:
            result.append([0,tr_nd[0][0]-1])

        for i in range(len(tr_nd)-1):
            result.append([tr_nd[i][1]+1, tr_nd[i+1][0]-1])

        if tr_nd[-1][1]!=seq_length-1:
            result.append([tr_nd[-1][1]+1,seq_length-1])

        return result   
    df= pd.DataFrame(columns=['State', 'Value', 'Precision/recall/f', 'n_eu',
                                       'n_neand', 'L',  'n_ref_pop','n_e_nd'])
    
    d =  MU * cut
    if n_neanderthal != 0:
        


        a = initA(cut, RR, lmbd_opt[0]/d, lmbd_opt[3])
        b = initB(MU, cut, lmbd_opt[0:3],   n_st+1)
        P=[0.95, 0.05]

        for idx in range(0, n_obs_seq):
            tracts_HMM =  get_HMM_tracts(viterbi(seq [idx], P, a, b))
            for k in range(N):
                for j in range(len(tracts_HMM[k])):
                    tracts_HMM[k][j][0]= cut * tracts_HMM[k][j][0]
                    tracts_HMM[k][j][1]= cut * (tracts_HMM[k][j][1]+1)-1
                    
                    
        

            real_tracts_in_states = []           
            real_tracts_in_states.append(tracts_eu(nd_true_tracts[idx], seq_length))
            real_tracts_in_states.append(nd_true_tracts[idx])
            cl_report = classification_rpt(confusion_mtrx(real_tracts_in_states, tracts_HMM))


            
            for j in range(N):
                df.loc[len(df.index)] = [j, cl_report[str(j)]['precision'], 'precision',idx, n_neanderthal, cut,
                                         n_ref_pop,  N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['recall'], 'recall',idx, n_neanderthal, cut, 
                                         n_ref_pop,  N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['f1-score'], 'f1-score',idx, n_neanderthal, cut,
                                         n_ref_pop, N_ND]
        
    else:
        d = MU * cut


        a = initA(cut,RR,lmbd_opt[0]/d, lmbd_opt[3])
        b = initBwN(MU, cut, lmbd_opt[0:3],   n_st+1)
        P=[0.95, 0.05]
        

        for idx in range(0, n_obs_seq):
            tracts_HMM =  get_HMM_tracts(viterbinD(seq [idx], P, a, b))
            for k in range(N):
                for j in range(len(tracts_HMM[k])):
                    tracts_HMM[k][j][0]= cut * tracts_HMM[k][j][0]
                    tracts_HMM[k][j][1]= cut * tracts_HMM[k][j][1]

            real_tracts_in_states = []
            real_tracts_in_states.append(tracts_eu(nd_true_tracts[idx], seq_length))
            real_tracts_in_states.append(nd_true_tracts[idx])
            cl_report = classification_rpt(confusion_mtrx(real_tracts_in_states, tracts_HMM))
            

            for j in range(N):
                df.loc[len(df.index)] = [j, cl_report[str(j)]['precision'], 'precision',idx, n_neanderthal, cut,n_ref_pop, N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['recall'], 'recall',idx, n_neanderthal, cut, n_ref_pop, N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['f1-score'], 'f1-score',idx, n_neanderthal, cut,n_ref_pop,  N_ND]


    return df


#lmbd_opt is mas[mas] - optimization is for every observable genome
def createDataFrame2(seq_length,RR,MU,cut,lmbd_opt, n_st, seq, n_neanderthal,  n_ref_pop, n_obs_seq, nd_true_tracts, nd_portion, N_ND):
    
    
    
    # return European tracts with input=Neanderthal tracts
    def tracts_eu(tr_nd, seq_length):
        result = []

        if tr_nd[0][0] > 0:
            result.append([0,tr_nd[0][0]-1])

        for i in range(len(tr_nd)-1):
            result.append([tr_nd[i][1]+1, tr_nd[i+1][0]-1])

        if tr_nd[-1][1]!=seq_length-1:
            result.append([tr_nd[-1][1]+1,seq_length-1])

        return result   
    df= pd.DataFrame(columns=['State', 'Value', 'Precision/recall/f', 'n_eu',
                                       'n_neand', 'L',  'n_ref_pop','n_e_nd'])
    
    d =  MU * cut
    if n_neanderthal != 0:
        


        for idx in range(0, n_obs_seq):


            a = initA(cut, RR, lmbd_opt[idx][0]/d, lmbd_opt[idx][3])
            b = initB(MU, cut, lmbd_opt[idx][0:3],   n_st+1)
            P=[0.97, 0.03]
            tracts_HMM =  get_HMM_tracts(viterbi(seq [idx], P, a, b))
            for k in range(N):
                for j in range(len(tracts_HMM[k])):
                    tracts_HMM[k][j][0]= cut * tracts_HMM[k][j][0]
                    tracts_HMM[k][j][1]= cut * (tracts_HMM[k][j][1]+1)-1
                    
                    
        

            real_tracts_in_states = []           
            real_tracts_in_states.append(tracts_eu(nd_true_tracts[idx], seq_length))
            real_tracts_in_states.append(nd_true_tracts[idx])
            cl_report = classification_rpt(confusion_mtrx(real_tracts_in_states, tracts_HMM))

            
            for j in range(N):
                df.loc[len(df.index)] = [j, cl_report[str(j)]['precision'], 'precision',idx, n_neanderthal, cut,
                                         n_ref_pop,  N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['recall'], 'recall',idx, n_neanderthal, cut, 
                                         n_ref_pop,  N_ND]
                df.loc[len(df.index)] = [j, cl_report[str(j)]['f1-score'], 'f1-score',idx, n_neanderthal, cut,
                                         n_ref_pop, N_ND]
        
