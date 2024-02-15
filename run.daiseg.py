import argparse
from importlib import reload
reload(argparse)
import argparse
import numpy as np
import HMM


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

seq_start, seq_end = f.readline().split(' ')
seq_start = int(seq_start)
t_n = float(f.readline())
t_ooa = float(f.readline())
t_i = float(f.readline())
a=float(f.readline())
f.close()



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





d = MU * L
n_st = np.max(seq)

lmbd_opt = d* np.array([t_n, t_ooa, t_i])
a = 0.05
A = HMM.initA(L,RR, lmbd_opt[0]/d, a)
B = HMM.initB(MU, L, lmbd_opt[0:3],   n_st+1)
P=[0.95, 0.05]

tracts_HMM_result = []
for idx in range(0, len(seq)):
    print(HMM.viterbi(seq[idx], P, A, B))

    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi(seq [idx], P, A, B))

    for k in range(N):
        for j in range(len(tracts_HMM[k])):
            tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]
            tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]

    tracts_HMM_result.append(tracts_HMM)


with open(args.o, "w") as f:
    for i in tracts_HMM_result:
        f.write(str(i)+'\n')





print(HMM.get_HMM_tracts([0,0,0,1,1]))