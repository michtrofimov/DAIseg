import argparse
import argparse
import numpy as np
import HMM
import sys
import useful as usfl

import EM 
parser = argparse.ArgumentParser(description='DAIseg') 



parser.add_argument('--location', type=str, help='File with first-last positions on chr')
parser.add_argument('--gaps', type=str, help='File with gaps')
parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser.add_argument('--EM_steps', type=str, help='number of EMsteps')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o', type= str, help = 'Name of output file' )
parser.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )
parser.add_argument('--obs_af', type=str, help='File with observations with respect to Africans')
parser.add_argument('--obs_archaic', type=str, help='File with observations with respect to Archaic reference genomes')
parser.add_argument('--obs_samples', type=str, help='File with samples names')


args = parser.parse_args()

with open(args.location,'r') as f1:

    seq_start, seq_end = f1.readline().split(' ')
    seq_start = int(seq_start)
    seq_end = int(seq_end.replace('\n',''))    

N = 2 # number of hidden states


f = open(args.HMM_par, 'r')
GEN_time = float(f.readline())
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())

#seq_start, seq_end = f.readline().split(' ')
#seq_start = int(seq_start)
#seq_end = int(seq_end.replace('\n',''))

Lambda_0=np.zeros(5)
Lambda_0[1] = float(f.readline())/GEN_time*MU*L
Lambda_0[2] = float(f.readline())/GEN_time*MU*L
Lambda_0[0] = float(f.readline())/GEN_time*MU*L
Lambda_0[4] = float(f.readline())/GEN_time*MU* L
Lambda_0[3] = float(f.readline())

f.close()




seq1, seq2 = [], []
with open(args.obs_af, 'r') as f1, open(args.obs_archaic, 'r') as f2:
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

n1=seq1.max()
n2=seq2.max()
seq=[]
for i in range(len(seq1)):
    seq.append(np.column_stack((seq1[i], seq2[i])))   

SEQ=np.array(seq)
N_st=SEQ.max()+1



if args.gaps is not None:
    with open(args.gaps,'r') as f:
        l=f.readline()
    m=l.replace('\n','').replace('[','').replace(']','').split(',')
    m=[[int(m[2*i]), int(m[2*i+1])] for i in range(int(len(m)/2))]   
    domain=usfl.exclude_gaps([[seq_start, seq_end]], m)


    #list of gaps numbers consistent with windows and starting position
    gaps_numbers=[]

    seq_start_mas=[]
    seq_end_mas=[]

    len_mas=[]
    for i in range(len(domain)):
        if (domain[i][0] // L)*L + (domain[0][0]% L) >= domain[i][0]:
            seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L))
        else:
            seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L)+L)
        if (domain[i][1]//L)*L-1+(domain[0][0]% L) <= domain[i][1]:
            seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L))
        else:
            seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L)-L)        
        
        len_mas.append(int((domain[i][1]-domain[i][0])/1000))

        if i!=len(domain)-1:
            gaps_numbers.append([int((domain[i][1]-domain[0][0])/1000),int((domain[i+1][0]-domain[0][0])/1000)] )

    domain=[[seq_start_mas[i], seq_end_mas[i]] for i in range(len(domain))]






    SEQ_mas=[]
    for i in range(len(len_mas)):
        p1=int((seq_start_mas[i]-seq_start_mas[0])/1000)
        p2=int((seq_end_mas[i]-seq_start_mas[0])/1000)
        SEQ_mas.append(SEQ[:,p1:(p2+1)])

else:
    SEQ_mas=[SEQ]
    gaps_numbers=[[]]
    seq_start_mas=[seq_start]
    domain=[[seq_start, seq_end]]





def run_daiseg(lmbd_opt,seq, n_st, idx, start):
    d = MU * L
    A = HMM.initA(L,RR, lmbd_opt[4]/d, lmbd_opt[3])
    B = HMM.initB(MU, L, lmbd_opt[0:3],   n_st)
    P=[0.97, 0.03]

    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi(seq [idx], P, A, B))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+start
           tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]+start-1

    return tracts_HMM

def run_daiseg_all(lmbd_0):
    tracts_HMM_mas=[]

    
    for idx in range(0, len(seq)):    
        tracts_HMM=[[],[]]
        for i in range(len(SEQ_mas)):
            tr=run_daiseg(lmbd_0, SEQ_mas[i], N_st, idx, seq_start_mas[i])
            for j in range(N):   
               for k in tr[j]:             
                   tracts_HMM[j].append( k )
 

        tracts_HMM_mas.append([tracts_HMM[j] for j in range(N)])
    return tracts_HMM_mas





def EM_function2(seq, lambda_0):
    P=[0.95, 0.05]
    n_EM_steps = 10
    epsilon = 1e-6
    d =  L * MU
    N_neanderthal = 6




    Lambda_new=EM.EM_algorithm2(P, seq,  N_st, MU, RR, lambda_0, epsilon, L)
    return Lambda_new
    
    
    
def EM_function3(seq,lambda_0):
    P=[0.95, 0.05]
    n_EM_steps = 10
    epsilon = 1e-6
    d =  L * MU
    N_neanderthal=6
    Lambda_new=EM.EM_algorithm3(P, seq,  N_st, MU, RR, lambda_0, epsilon, L)
    return Lambda_new
    
    
    
P=[0.95, 0.05]
n_EM_steps = 10
epsilon = 1e-6

def EM_gaps(seq, lambda_0, n_st):
    return EM.EM_algorithm_gaps(P, seq, n_st, MU, RR, lambda_0, epsilon, L, int(args.EM_steps), gaps_numbers )

if args.EM=='no': 
    Tracts_HMM_mas = run_daiseg_all(Lambda_0)




        
if args.EM=='yes': 

    Lambda_opt = EM_gaps(SEQ, Lambda_0, N_st)    
    Tracts_HMM_mas = run_daiseg_all(Lambda_opt)
        
      
with open(args.obs_samples,'r') as f:
    names=f.readlines()

names=[str(names[i].replace('\n','')) for i in range(len(names))]
print(names)  

with open(args.o, "w") as f:
   for i in range(len(Tracts_HMM_mas)):
       if i % 2 ==0:
           f.write(names[int(i // 2)]+'\t0\t'+str(Tracts_HMM_mas[i][1])+'\n')
       else:
           f.write(names[int(i // 2)]+'\t1\t'+str(Tracts_HMM_mas[i][1])+'\n')       



    

    
    










