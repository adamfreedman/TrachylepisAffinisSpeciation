# Demographic inference

## dadi
To test the fit of a suite of demographic models to patterns of genomic variation in our sampled *T. affinis* populations, we used [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master/), a method that infers demographic history using a diffusion approximation to the allele frequency spectrum. For more details on the method see [Gutenkunst et al. 2009](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695). Dadi can compare up to three populations. We perform demographic inference for three types of models based upon broad patterns of population structure observed in our data:
* two populations: forest vs. ecotone
* two populations: populations in Southwest Province, Cameroon (SWP) vs. all other forest populations (SF)
* three populations: SWP, SF and ecotone
    * first population split between ecotone and (SWP,SF)
    * second split between SWP and SF

### Conversion to dadi format
For our analyses of empirical data, we converted our genotypes (without minor allele frequency or missingness filters; see [variant_filtering](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/tree/master/variant_filtering) for details to dadi input format using the R package radiator version 0.0.4, using the vcf2dadi function. This function takes the vcf file and a tab-separated text "strata" file that indicates the population (stratum) to which an individual sample belongs, eg. for forest vs. ecotone, "ECO" and "FOR", implemented in this fashion:
```bash
vcf2dadi(data="oneperrad_trachylepis_RAD.vcf",strata="taffinis_dadi_strata_forVecotone.tsv",pop.levels = c("ECO","FOR"),common.markers = TRUE)
```
In later analyses we conducted simulations of genotypes for the purpose of examining F<sub>ST</sub> distributions (see below) we discovered that this function had been deprecated, and that its successor function, tidy_vcf generated unresolvable errors. I thus wrote a simple python script, ConvertVcfToDadi.py, for performing this conversion. I performed tests to confirm that, although the output order of genomic positions differs from that of vcf2dadi, the contents of the files with respect to allele counts are identical.

 

 
