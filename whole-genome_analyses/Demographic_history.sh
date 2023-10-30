# Demographic History

# PSMC
vcfutils_pl=~/bcftools-1.3.1/vcfutils.pl

samtools mpileup -C50 -uf ${ref} ind.sotred.rmdup.bam | bcftools call -c -V indels | ${vcfutils_pl} vcf2fq -d 6 -D 36 | gzip > ind.fq.gz
fq2psmcfa} -q20 ind.fq.gz > ind.psmcfa

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

cp ~/pop1_out/model.final.json estimate_iter30_0.1ma/hw_model.final.json
cp bw_analysis/model.final.json estimate_iter30_0.1ma/bw_model.final.json
cp sy_analysis/model.final.json estimate_iter30_0.1ma/sy_model.final.json
run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/hw_g_0.1ma.pdf $PWD/hw_model.final.json docker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/bw_g_0.1ma.pdf $PWD/bw_model.final.json docker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/sy_g_0.1ma.pdf $PWD/sy_model.final.jsondocker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/year_0.1ma.pdf -g 3.783783 $PWD/hw_model.final.json docker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/year_0.1ma.pdf -g 4.171668 $PWD/bw_model.final.json docker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/year_0.1ma.pdf -g 4.906858 $PWD/sy_model.final.jsondocker run --rm -v $PWD:/home/haoyan/Podoces/Corvus_cornix_V5/smc++/estimate_iter30_0.1ma/ terhorst/smcpp:latest plot $PWD/all_g_0.1ma $PWD/hw_model.final.json $PWD/bw_model.final.json $PWD/sy_model.final.json
smc++ plot pop1.pdf -g 2.5 ~/pop1_analysis/model.final.json

