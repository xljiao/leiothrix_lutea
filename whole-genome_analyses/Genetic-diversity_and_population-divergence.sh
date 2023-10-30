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

