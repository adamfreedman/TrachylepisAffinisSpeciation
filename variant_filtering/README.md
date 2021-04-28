# Variant filtering
## LD filtering
RAD-seq data sets typically produce more than one variant call per RAD locus. Given that these loci are relatively short, single nucleotide variants (SNVs) found on the same locus will be in strong linkage disequilibrium (LD). Inclusion of SNVs in strong LD has the potential to distort evolutionary inferences. Therefore, our first step in filtering the genotypes generated with [STACKS](https://catchenlab.life.illinois.edu/stacks/) is to downsample our variant call format (vcf) genotypes by randomly selecting one polymorphic site per RAD locus. We do this with a custom python script [SampleOneSnpPerRadLocusFromVcf.py](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/blob/master/variant_filtering/scripts/SampleOneSnpPerRadLocusFromVcf.py):

```bash
SampleOneSnpPerRadLocusFromVcf.py trachylepis_RAD.vcf
```
which produces the new, downsampled file *oneperrad_trachylepis_RAD.vcf*.

## Locus and individual filtering
In our paper, we demonstrate through simulations that, due to the demographic history of *Trachylepis affinis*, filtering (particularly on minor allele frequency), will lead to incorrect selection of the best demographic model. Nevertheless, it has been shown elsewhere that F<sub>ST</sub> estimators can be biased by rare alleles, and that ancestry proportions inferred by software such as [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655.full), which we use in our paper, can be strongly influenced by filtering on MAF and other dataset features, we produce a filtered data set to generate estimates of population differentiation and ancestry proportions, to be compared with results on unfiltered data. 

First, we tabulated the fraction of missing genotypes per individual with our script [summarizemissing_by_ind.py](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/blob/master/variant_filtering/scripts/summarizemissing_by_ind.py), which produces [trachylepis_RAD_oneperrad_missingstats_by_ind.txt](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/blob/master/variant_filtering/data/trachylepis_RAD_oneperrad_missingstats_by_ind.txt), a table with three columns: sample id, count of missing genotypes, and frequency of missing genotypes. Then, we extracted individuals with less than 25% of genotypes as follows: 
```bash
grep -v id trachylepis_RAD_oneperrad_missingstats_by_ind.txt |awk '$3<=0.25{print $0}' > ids_maxmis_25pcent_keep.txt
``` 
Then we use [vcftools](https://vcftools.github.io/examples.html) to only retain polymorphic sites with a MAF >= 0.05, and to exclude all individuals with a fraction of missing genotypes greater than 25%, as follows:
```bash
vcftools --vcf oneperrad_trachylepis_RAD.vcf --recode --keep ids_maxmis_25pcent_keep.txt --maf 0.05 --out oneperrad_trachylepis_RAD_maxmis25pcent_maf05
```
This produces a set of filtered genotypes for 194 individuals and 14896 loci.
