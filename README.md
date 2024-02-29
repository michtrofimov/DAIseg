# DAIseg
DAIseg method is created to detect ancient introssed segments using unadmixed outgroup population and several reference archaic genomes.

# Files's summary
*  __outgroup.txt__(Africa), __observations.txt__(European), __archaic.txt__(Neanderthals) are .txt files which consist  of the samples's ids written in a column
```note
NA18484
NA18489
GM19129
```
*  __all.chr22.vcf.gz{.tbi}__ files where  all reference genomes(Outgroup and Archaic) and observable samples simultaneously with snps only (excluding indels, deletions etc.) is in it. The main reason of it is to avoid inconsistencies.
*  __obs.....txt__ is 
```note
0 0 0 0 
0 2 0 1
0 0 1 0
```
Here is the example of four  observations sequences (one for  each column) with respect to one fixed reference populations. Each row corresponds to the number of variants obtained in the window of size L=1000.
*  __par.file.txt__
```note
1.25e-08    #mutation rate μ
1e-08    #recombination rate
1000    #window size
start end    #position of the first SNP in .vcf file
lambda_arch    #the mean value of derived alleles in a window of size L accumulated during time t_arch and mutation rate μ 
lambda_split    #the mean value of derived alleles in a window of size L accumulated during time t_split and mutation rate μ 
lambda_intr    #the mean value of derived alleles in a window of size L accumulated during time t_intr and mutation rate μ 
0.025    #admixture proportion of archaic introgression
```

# Pipeline briefly
0. (optionally) Run __panel.preparation.sh__
1. Run __./script.eu.sh__
2. (optionally)





## Step 0. Data preparation. Merging 1000GP  and Archaic genomes
The link to download 1000GP panel is http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 

The link to download archaic samples http://cdna.eva.mpg.de/neandertal/Vindija/VCF/ 



In the  file __panel.preparation.sh__ change names of CURRENTDIR (path to current directory), NAME1000(name of vcf file of 1000GP),  DIRNEAND(directory with neanderthal genomes) and run it. The resulting vcf.gz file is __all.chr22.vcf.gz{.tbi}__

## Step 1. Data preparation. Make observations
You need in  __all.chr22.vcf.gz{.tbi}__,  __outgroup.txt__, __observations.txt__, __archaic.txt__ to run  __./script.eu.sh__ and to  make observation  files __obs.neand.txt__, __obs.outgroup.txt__ and file with stard-end positions __par.file.txt__



By default  __par.file.txt__ is 





## Step 2. Estimation of parameters (optionally)
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


## Step 3. Run DAI.seg HMM 
The command line of  HMM.py is 

__HMM.py --HMM_par par.file.txt --obs1 obs1.txt --obs2 obs2.txt --o output.txt__



# Working with 1000GP (optionally)
If you want to work with samples from 1000GP there are two convinient files panel.preparation.sh and script.eu.sh

The link to download 1000GP panel is http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz where i is in 1..22.

The link to download archaic samples http://cdna.eva.mpg.de/neandertal/Vindija/VCF/ 


__panel.preparation.sh__ is file which prepare panel (remove indelsand multiallelic snps) and glue it with 3 archaic samples. 


__script.eu.sh__ 
