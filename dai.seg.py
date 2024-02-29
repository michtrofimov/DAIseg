import argparse
import argparse
import numpy as np
import HMM
import EM
parser = argparse.ArgumentParser(description='DAIseg') 

parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser.add_argument('--EM_times', type=str, help='Do EM for one or all samples? One or all') 
parser.add_argument('--obs_out', type= str, help='File with observations with respect to Outgroup')
parser.add_argument('--obs_neand', type= str, help='File with observations with respect to Archaic')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o', type= str, help = 'Name of output file' )

args = parser.parse_args()



f = open(args.HMM_par, 'r')
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())

seq_start, seq_end = f.readline().split(' ')
seq_start = int(seq_start)
seq_end = int(seq_end.replace('\n',''))

Lambda_0=np.zeros(4)
Lambda_0[1] = float(f.readline())
Lambda_0[2] = float(f.readline())
Lambda_0[0] = float(f.readline())
a=float(f.readline())
f.close()


GEN_time=29
N = 2 # number of hidden states

seq1, seq2 = [], []
with open(args.obs_out, 'r') as f1, open(args.obs_neand, 'r') as f2:
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

SEQ=np.array(seq)
N_st=SEQ.max()+1

def run_daiseg(lmbd_opt,seq, n_st, idx):
    d = MU * L
    A = HMM.initA(L,RR, lmbd_opt[0]/d, a)
    B = HMM.initB(MU, L, lmbd_opt[0:3],   n_st)
    P=[0.97, 0.03]

    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi(seq [idx], P, A, B))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+seq_start
           tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]+seq_start-1

    return tracts_HMM





def EM_function(ind):

    n_EM_steps = 10
    epsilon = 1e-12
    d =  L * MU
    N_neanderthal=6
    bnds = ((Lambda_0[0], Lambda_0[0]+d*1000000/GEN_time), (Lambda_0[1], Lambda_0[1]+d*100000000/GEN_time), 
           (Lambda_0[2], Lambda_0[2]+d*50000000/GEN_time), (a-0.000001, a+0.000000001))

    
    ind = 0
    Lambda_new=EM.opt_params(MU,RR,SEQ, N_st, Lambda_0, epsilon, N_neanderthal, bnds, ind, L, n_EM_steps, a)
    return Lambda_new
    

tracts_HMM_result = []
if args.EM=='no':
    for idx in range(0, len(seq)):
        tracts_HMM_result.append(run_daiseg(Lambda_0, SEQ, N_st, idx))

        
if args.EM=='yes': 
    if args.EM_times=='one':
        Lambda_opt = EM_function(0)
        for idx in range(0, len(seq)):   
            tracts_HMM_result.append(run_daiseg(Lambda_opt, SEQ, N_st, idx)) 

        
    if args.EM_times=='all':
        for idx in range(0, len(seq)): 
            Lambda_opt = EM_function(idx)
            tracts_HMM_result.append(run_daiseg(Lambda_opt, SEQ, N_st, idx))             


with open(args.o, "w") as f:
   for i in tracts_HMM_result:
       f.write(str(i)+'\n') 



    

    
    










