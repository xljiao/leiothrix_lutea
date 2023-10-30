# Population Structure and Phylogenetic Analysis

# Prune autosomal SNPs at high linkage disequilibrium across all individuals using Plink v1.9
# First, the vcf file is converted to a ped/map file in plink format
vcftools --vcf filtered_autosomal.vcf --plink --out filtered_autosomal

# We set the independent pair-wise filter at a correlation threshold of 0.1 for a window size of 50 kilobases (kb) and a step size of 1 kb
plink --file filtered_autosomal --indep-pairwise 50kb 1 0.1 --allow-extra-chr --out filtered_autosomal --threads 40
plink --file filtered_autosomal --extract filtered_autosomal.prune.in --out filtered_autosomal_ld --allow-extra-chr --make-bed

# We performe principal component analysis (PCA) using Plink v1.9 
plink --bfile filtered_autosomal_ld --allow-extra-chr --pca 10 --out filtered_autosomal_ld  --threads 20

# ADMIXTURE v1.3.0 is run to analyze population structure with cross-validation and 100 bootstraps, and using K values ranging from 2 to 6. 
plink --bfile filtered_autosomal_ld --recode 12 --out filtered_autosomal_ld_ad --allow-extra-chr
for K in 1 2 3 4 5 6 
do
admixture --cv -B100 -j40 filtered_autosomal_ld_ad.ped $K | tee log${K}.out
done

# We use FastTree2 to reconstruct an ML tree based on the general time reversible (GTR) model of nucleotide substitution using concatenated SNPs of autosomal dataset.
# Convert vcf format to fasta format using vcf2phylip.py (https://github.com/edgardomortiz/vcf2phylip)
python ${vcf2phylip} -i filtered_autosomal.vcf -o filtered_autosomal -f
${fasttree} -nt -gtr < filtered_autosomal.min4.fasta > filtered_autosomal.tree

