import argparse
import argparse
import numpy as np
import HMM


import EM0104
parser = argparse.ArgumentParser(description='DAIseg') 

parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o', type= str, help = 'Name of output file' )
parser.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )


args = parser.parse_args()


N = 2 # number of hidden states


f = open(args.HMM_par, 'r')
GEN_time = float(f.readline())
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())

seq_start, seq_end = f.readline().split(' ')
seq_start = int(seq_start)
seq_end = int(seq_end.replace('\n',''))

Lambda_0=np.zeros(5)
Lambda_0[1] = float(f.readline())/GEN_time*MU*L
Lambda_0[2] = float(f.readline())/GEN_time*MU*L
Lambda_0[0] = float(f.readline())/GEN_time*MU*L
Lambda_0[4] = float(f.readline())/GEN_time*MU* L
Lambda_0[3] = float(f.readline())


f.close()




seq1, seq2 = [], []
with open('obs.outgroup.txt', 'r') as f1, open('obs.neand.txt', 'r') as f2:
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

print(n1,n2)

seq=[]
for i in range(len(seq1)):
    seq.append(np.column_stack((seq1[i], seq2[i])))   

SEQ=np.array(seq)
N_st=SEQ.max()+1

def run_daiseg(lmbd_opt,seq, n_st, idx):
    d = MU * L
    A = HMM.initA(L,RR, lmbd_opt[4]/d, lmbd_opt[3])
    B = HMM.initB(MU, L, lmbd_opt[0:3],   n_st)
    P=[0.97, 0.03]

    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi(seq [idx], P, A, B))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+seq_start
           tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]+seq_start-1

    return tracts_HMM






def EM_function2(seq, lambda_0):
    P=[0.95, 0.05]
    n_EM_steps = 10
    epsilon = 1e-6
    d =  L * MU
    N_neanderthal = 6




    Lambda_new=EM0104.EM_algorithm2(P, seq,  N_st, MU, RR, lambda_0, epsilon, L)
    return Lambda_new
    
    
    
def EM_function3(seq,lambda_0):
    P=[0.95, 0.05]
    n_EM_steps = 10
    epsilon = 1e-6
    d =  L * MU
    N_neanderthal=6
    Lambda_new=EM0104.EM_algorithm3(P, seq,  N_st, MU, RR, lambda_0, epsilon, L)
    return Lambda_new
    

tracts_HMM_result = []
if args.EM=='no':
    for idx in range(0, len(seq)):
        tracts_HMM_result.append(run_daiseg(Lambda_0, SEQ, N_st, idx))




        
if args.EM=='yes': 
    for idx in range(0, len(seq)):
        if args.EM_est == 'coal':
            
            Lambda_opt = EM_function3(SEQ, Lambda_0)

        
        if args.EM_est == 'all':
            Lambda_opt = EM_function2(SEQ, Lambda_0)
    
    
     
        tracts_HMM_result.append(run_daiseg(Lambda_opt, SEQ, N_st, idx)) 
        
      


with open(args.o, "w") as f:
   for i in tracts_HMM_result:
       f.write(str(i)+'\n') 



    

    
    










