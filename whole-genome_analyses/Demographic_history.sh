# Demographic History

# PSMC
vcfutils_pl=~/bcftools-1.3.1/vcfutils.pl

samtools mpileup -C50 -uf ${ref} ind.sotred.rmdup.bam | bcftools call -c -V indels | ${vcfutils_pl} vcf2fq -d 6 -D 36 | gzip > ind.fq.gz
fq2psmcfa} -q20 ind.fq.gz > ind.psmcfa

# using the python scripts remove_scf.py to remove the sex chromosomes.
/usr/bin/python2.7 remove_scf.py chr.list ind.psmcfa

# run psmc
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ind.psmc ind.psmcfa 

# perform 100 bootstrap
seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc ind.psmcfa | sh &

# plot
cat ind_round-*.psmc > ind_combined.psmc
~/psmc/utils/psmc_plot.pl -u 4.6e-9 -g 2.5 -M 'pop1' -p pop1_plot ind_combined.psmc

# smc++
bgzip -c filtered_autosomal.vcf > filtered_autosomal.vcf.gz
tabix -p filtered_autosomal.vcf.gz
vcf=filtered_autosomal.vcf.gz
for chr in `cat chr.lst`
do
smc++ vcf2smc --cores 10 --ignore-missing $vcf ~/pop1_out/$chr.smc.gz $chr pop1:ind1,ind2,ind3
done

smc++ estimate --timepoints 1 1000000 --em-iterations 30 --cores 20 --spline pchip -o ~/pop1_analysis/ 4.6e-9 ~/pop1_out/*.smc.gz

smc++ plot plot.pdf -g 2.5 ~/pop1/model.final.json ~/pop2/model.final.json

# joint-demographic history
# we created the dadi format SFS files by python script 'vcf_to_dadi.py' from Popgen Pipeline Platform (PPP, https://github.com/jaredgk/PPP) and examined all dempgraphic models using dadi_pipeline (https://github.com/dportik/dadi_pipeline).
