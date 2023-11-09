# calculate the pairwise genetic distance matrix
vcftools --vcf filtered_autosomal.vcf --fst-window-size 50000 --weir-fst-pop pop1 --weir-fst-pop pop2 --out pop1_pop2
