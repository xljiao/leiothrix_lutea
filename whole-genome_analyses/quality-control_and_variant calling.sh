# Quality-control 
# FastQC v0.11.8 to assess read quality, Fastp v0.21 to trim poor-quality reads and Illumina adapters with default settings

# FastQC
fastqc -o fastqc/ -t 40 -f fastq *.fastq.gz 

# fastp
fastp -i ind.R1.fastq.gz -o ind.R1.fp.fastq.gz -I ind.R2.fastq.gz -O ind.fp.fastq.gz -q 15 -u 50 -n 5 -l 150 -t 0 -h ind.fastp.html

# Assemble reads with reference genome using bwa
ref=ASM1340044v1.fa
bwa index -a bwtsw ${ref}
  
bwa mem ${ref} ind.R1.fp.fastq.gz ind.R2.fp.fastq.gz | samtools view -b -S > ind.bam

# BAM files sort
samtools sort ind.bam -o ind.sorted.bam 

# BAM files remove duplicates 
samtools rmdup ind.sorted.bam ind.sorted.rmdup.bam 

# Call variants using HaplotypeCaller in GATK v4.0.9.
gatk HaplotypeCaller -R ${ref} -ERC GVCF -I ind.sorted.rmdup.bam -O ind.g.vcf --genotyping-mode DISCOVERY --pcr-indel-model CONSERVATIVE /
--sample-ploidy 2 --min-base-quality-score 5 --kmer-size 10 --kmer-size 25 --native-pair-hmm-threads 40
gatk CombineGVCFs -R ${ref} -V gvcf.list -O com.g.vcf.gz
gatk GenotypeGVCFs -R ${ref} -V com.g.vcf.gz -O com.vcf.gz

# VCFtools v0.1.13 to remove indels and filter VCF according to the following criteria
/home/software/vcftools_0.1.13/bin/vcftools --vcf com.vcf.gz \
--min-alleles 2 --max-alleles 2 --maf 0.05 --max-maf 0.95 \
--minQ 30 --minDP 5 --maxDP 100 --max-missing 0.95 --remove-indels \
--recode --out flitered_vcf

# Retain autosomal SNPs using BLAST+2.2.26
parus_major_genome=parus.major.fa
makeblastdb -in ${parus_major_genome} -dbtype nucl  -title parus_major_fa
blastn -query ${ref} -db ${parus_major_genome} -evalue 1e-40 -outfmt 7 -out pm_ll_blast.out





