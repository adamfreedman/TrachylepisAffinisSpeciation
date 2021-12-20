# Principal components analysis (PCA)
To visualize population structure among our 14 sampled populations of *Trachyleipis affinis*, we performed a PCA on the *T. affinis* genotypes. Results are presented in Supplementary Data Figure S3 of our manuscript.

First, we converted our vcf file of genotypes to [plink](https://www.cog-genomics.org/plink2),version 1.90b3s, format using vcftools as follows:
```bash
vcftools --vcf oneperrad_trachylepis_RAD.vcf --plink --out oneperrad_trachylepis_RAD
```
Next we use *plink* to recode variants from basespace to integer,so that zeroes reflect missing data, while 1 and 2 represent alternative variants.
```bash
plink --file oneperrad_trachylepis_RAD --recode 12 --out oneperrad_taffinis_intrecode
```

Because vcftools doesn't like "Un" chromosomes, it by default converts the chromsome field (field 1) to zero. For this reason, and to make the PCA estimation play nice with our data, we rewrite the map file:
```bash
awk '{print "1","tafsnp"$2,"-9","-9"}' oneperrad_taffinis_intrecode.map > new.map
mv new.map oneperrad_taffinis_intrecode.map
```
Then, we perform PCA with the *smartpca* module of [Eigensoft](https://www.hsph.harvard.edu/alkes-price/software/), version 6.1.4. This requires creating a parameter file that specifies the data files, number of PCs to output, etc. provided here as *taffpcs.par*.
```bash
/PATH/TO/smartpca -p taffpcs.par > taffpca.log
```
Finally, for to facilitate downstream plotting and visualization, we append population ids to the PC loadings for individuals. We first extract the individual-level data:
```bash
tail -208 taffpcs > taff_smartpca_10pcs_Rtable.txt
```
Then, we use a table that maps individual ids to populations to append population ids to the PCA results:
```bash
python AppendPopToPCs.py ids_renumbered.txt taff_smartpca_10pcs_Rtable.txt
```
Finally, we create the pca plot for PCs 1 and 2 in R using the script [pcaplotcode.R]](https://github.com/adamfreedman/TrachylepisAffinisSpeciation/tree/master/pca).
