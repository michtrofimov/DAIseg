# DAIseg
DAIseg method is created to detect ancient introssed segments using unadmixed outgroup population and several reference archaic genomes. 

__Input__: .vcf.gz{.tbi} file where merging all Neanderthal, Outgroupand  ingroup observable samples together and three .txt files with ids of each group.

__Output__: .txt file where each line corresponds to the array of tracts with modern and archaic ancestries.


# Pipeline briefly
0. (optionally) Run __panel.preparation.sh__ with samples' name files to merge 1000GP, neanderthal samples and obtain .vcf.gz file.
1. Using .vcf.gz{.tbi} and files with samples's names to run __./script.eu.sh__ to make observation files.
3. Run __dai.seg.py__ to obtain archaic tracts of samples from  __observations.txt__ with using or no-using  EM algorithm.



# Files's summary
*  __outgroup.txt__(Africa), __archaic.txt__(Neanderthals)  and __obs.samples.txt__(European),are .txt files which consist  of the samples's ids of reference Africans and Neanderthals and observable Europeans written in a column
```note
NA18484
NA18489
GM19129
```


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

By default, the splitting times 550.000, 70000 and 55.000 are used to make lambda_arch, lambda_split and lambda_intr on Step 1 and make initiall guess for EM algorithm on Step 2.


*  __all.chr22.vcf.gz{.tbi}__ files where  all reference genomes(Outgroup and Archaic) and observable samples simultaneously with snps only (excluding indels, deletions etc.) is in it. The main reason of it is to avoid inconsistencies.
  
* __output.txt__ is a  file 
```note
[[[t_1,t_2], [t_3,t_4], [t_5,t_6]], [[t_2+1, t_3-1], [t_4+1, t_5-1]]]
[[[t'_1,t'_2], [t'_3,t'_4], [t'_5,t'_6]], [[t'_2+1, t'_3-1], [t'_4+1, t'_5-1]]]
...
...
```
where each two lines correspond to the one diploid sample from obs.samples.txt.





## Step 0. Merging 1000GP  and Archaic genomes
Download 1000GP panel is 
>http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 

and  archaic samples 
>http://cdna.eva.mpg.de/neandertal/Vindija/VCF/ 

Make .txt files with samples's names  __outgroup.txt__, __obs.samples.txt__, __archaic.txt__

Add full path to files  of 1000GP and three neanderthals to variables __$NAME1000__ and __$n1, $n2, $n3__ in  __panel.preparation.sh__ , change $CHR variable and run script. The resulting vcf.gz file is __all.chr22.vcf.gz{.tbi}__

## Step 1.  Make observations

You need in  __all.chr22.vcf.gz{.tbi}__,  __outgroup.txt__, __observations.txt__, __archaic.txt__ to run  

>__./make.observations.sh__

and to  make observation  files __obs.neand.txt__, __obs.outgroup.txt__ and file with default parameters and start-end positions __par.file.txt__ (see the File's summary paragraph). 





## Step 2.1 Run DAI.seg without EM algorithm
>  python dai.seg.py  --EM no --HMM_par par.file.txt  --o output.tracts.txt

where file par.file.txt was created on the previous Step. 

## Step 2.1 Run DAI.seg using EM algorithm

There are three main parameters in DAIseg model:
1. Mean value of derived alleles in a window of size L accumulated during time t_arch  __lambda_archaic__,
2. Mean value of derived alleles in a window of size L accumulated during time t_split  __lambda_split__,
3. Mean value of derived alleles in a window of size L accumulated during time t_intr  __lambda_intr__

where  t_arch, t_split, t_intr are COALESCENT times 

__par.file.txt__ obtained on the Step 1 could be used as the initial guess for EM algorithm.

There are two possible options to estimate parameters: 
use only __one__ observable sample 
> python dai.seg.py --EM yes --EM_times one --o output.tracts.txt

to obtain single __par.file.0.txt__ file with parameters 
or use   __--EM all__ observable samples
 
> python dai.seg.py --EM yes --EM_times all--o output.tracts.txt

to make estimations of coalescent times for all samples.
