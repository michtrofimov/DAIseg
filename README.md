# DAIseg
DAIseg method is created to detect ancient introssed segments using unadmixed outgroup population and several reference archaic genomes(instead of Skov et al'18).

## Data preparation

## Estimation of parameters
There are three main parameters in DAIseg model:
1. mean coalescent time between Neanderthals and AMH t_n,
2. mean coalescent time between Outgroup and Ingroup t_ooa,
3. introgression time of archaic segments into Ingroup t_i.

For example, Outgroup = some African population, Ingroup = Modern Europeans.

The resulting values of parameters are in par.file.txt with the following structure(line by line):
1. Mutation rate1.25e-8
2. Recombination rate 1e-8
3. Window size L
4. Full sequence length in bp 
5. t_n
6. t_ooa
7. t_i


## HMM 
The command line of  HMM.py is 
HMM.py --HMM_par par.file.txt --obs1 obs1.txt --obs2 obs2.txt --o output.txt
