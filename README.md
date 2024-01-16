# DAIseg
DAIseg method is created to detect ancient introssed segments using unadmixed outgroup population and several reference archaic genomes(instead of Skov et al'18).

## Data preparation
we need to prepare to files with observations connected with Outgroup and Neanderthals.

## Estimation of parameters
There are three main parameters in DAIseg model:
1. mean coalescent time between Neanderthals and AMH __t_n__,
2. mean coalescent time between Outgroup and Ingroup __t_ooa__,
3. introgression time of archaic segments into Ingroup __t_i__.

For example, Outgroup = some African population, Ingroup = Modern Europeans.

The resulting values of parameters are in par.file.txt with the following structure(line by line):
* Mutation rate,
* Recombination rate, 
* Window size,
* Full sequence length in bp, 
* t_n,
* t_ooa,
* t_i.


## HMM 
The command line of  HMM.py is 
__HMM.py --HMM_par par.file.txt --obs1 obs1.txt --obs2 obs2.txt --o output.txt__
