# Genetic Diversity and Population Divergence
# calculate the heterozygosity 
vcfutils_pl=/home/software/bcftools-1.3.1/vcfutils.pl
samtools mpileup -t AD,ADF,ADR,DP,SP -Bugf ${ref} ind.sorted.rmdup.bam | bcftools call -c - \
 | ${vcfutils_pl} vcf2fq -d 10 -D 100 > ind.fq

seqtk seq -a -q20 -n N ind.fq > ind.f
sed -i '/>/d' ind.f
cat ind.f | tr a-z A-Z > ind.fa

awk -v RS='A' 'END {print "A "--NR}' ind.fa >> ind.a
awk -v RS='T' 'END {print "T "--NR}' ind.fa >> ind.a
awk -v RS='C' 'END {print "C "--NR}' ind.fa >> ind.a
awk -v RS='G' 'END {print "G "--NR}' ind.fa >> ind.a
awk -v RS='K' 'END {print "K "--NR}' ind.fa >> ind.a
awk -v RS='M' 'END {print "M "--NR}' ind.fa >> ind.a
awk -v RS='R' 'END {print "R "--NR}' ind.fa >> ind.a
awk -v RS='S' 'END {print "S "--NR}' ind.fa >> ind.a
awk -v RS='W' 'END {print "W "--NR}' ind.fa >> ind.a
awk -v RS='Y' 'END {print "Y "--NR}' ind.fa >> ind.a
awk -v RS='N' 'END {print "N "--NR}' ind.fa >> ind.a 

cat ind.a | awk '{print $2}' | head -n 10 | awk '{ sum += $1; } END { print  sum}' >> ind.B
cat ind.a | awk '{print $2}' | sed -n '5,10p' | awk '{ sum += $1; } END { print  sum}' >> ind.B
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' ind.B > ind.B2
awk '(sum=$2/$1){print sum}' ind.B2 >  ${heter_d}${ind}.B3

# calculate the Tajima's D, nucleotide diversity (π), and population genetic differentiation (FST) 
vcftools --vcf filtered_autosomal.vcf --keep pop.list --TajimaD 50000 --out pop_tajima
vcftools --vcf filtered_autosomal.vcf --keep pop.list --window-pi 50000 --out pop_pi
vcftools --vcf filtered_autosomal.vcf --weir-fst-pop pop1.list --weir-fst-pop pop2.list --fst-window-size 50000 --out pop1_2

# calculated the dxy, python scripts were developed by Simon Martin (https://github.com/simonhmartin)
python2.7 parseVCF.py -v filtered_autosomal.vcf -o filtered_autosomal.calls -m geno
python2.7 mergeCalls.py -f filtered_autosomal.calls -m all -I ${ref}.fai -o filtered_autosomal.calls.merge

python2.7 egglib_sliding_window.py \
-i filtered_maf_ld.calls.merge -o filtered_maf_ld.100k.csv \
-w 50000 -m 5 -s 50000 -a dxy \
-p "pop1[ind1,ind2,ind3]; pop2[ind1,ind2,ind3]
