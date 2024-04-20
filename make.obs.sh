#!/bin/bash
CHR=22
panelfinal=all.chr${CHR}.vcf.gz

eu=$3
arch=$2
out=$1

for i in ${out} ${arch} 
do
	bcftools query -S $i  -f '%POS\n'   ${panelfinal} > $i.positions.txt # print number of postions in separate text file
	bcftools query -S $i  -f '[%GT]\n'  ${panelfinal}  |sed 's/[^0]//g'  | awk '{ print length }'> $i.0.txt # считаем количество ноликов в строке
	bcftools query -S $i -f '[%GT]\n'  ${panelfinal} |sed 's/[^1]//g'  | awk '{ print length }'> $i.1.txt # считаем количество единичек в строке
	paste  $i.0.txt $i.1.txt > $i.al.spec.txt #join 

	awk '{
		if ($1>"0"&& $2>"0"){print "1\t1"}
		else 
			if ($1>"0"&& $2=="0"){print "1\t-1"}
			else {print "-1\t1" }
		}' $i.al.spec.txt > $i.spec.txt	
done 

for i in ${out} ${arch} 
do
paste $i.positions.txt $i.spec.txt  > chr${CHR}.$i.reference.txt
done

for i in ${out} ${arch} 
do
rm $i.spec.*
rm $i.al.*
rm $i.0.*
rm $i.1.*
rm $i.positions.*
done


bcftools query -S obs.samples.txt  -f '[%GT ]\n'  ${panelfinal} |sed  's/|/ /g' >  obs.chr${CHR}.ingroup.txt


python3 obs.py chr${CHR}.${out}.reference.txt chr${CHR}.${arch}.reference.txt obs.chr${CHR}.ingroup.txt

for i in ${out} ${arch} 
do
rm chr${CHR}.$i.reference.txt
done
rm obs.chr${CHR}.ingroup.txt



rm positions.chr${CHR}.txt










