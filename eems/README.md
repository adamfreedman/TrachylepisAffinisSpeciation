# EEMS
To visualize spatial patterns of gene flow, and geographic regions with less (or more) gene flow than expected under an isolation-by-distance model, we estimated effective migration surfaces with [eems](https://github.com/dipetkov/eems). To do this, we first converted our unfiltered genotypes to binary plink (i.e. bed) format with the program plink. 

Next, we calculated parwise genetic distances between samples using the eems program bed2diffs as follows:
```bash
bed2diffs_v1 --bfile oneperrad_trachylepis --nthreads 1
``` 
where "oneperrad_trachylepis" is the prefix of the plink SNP files.

Next we estimate effective migration surfaces with runeems_snps. As this method is based on MCMC, we perform 5 separate runs with MCMC parameters, the number of demes (200), the number of sites (i.e. SNPs = 56055), the number of individuals (208), and the paths to input and output files in a required configuration file, an example of which is provided as taff_params_1.ini. An eexample command line for one iteration of migration surface estimation is then:
```bash
runeems_snps --params taff_params_1.ini  
```

The datapath in taff_params_1.ini specifies th full path and prefix of three required files: the datapath.diffs file output by bed2diffs_v1, the datapath.coord file, which specifies the geographic coordinates of each sample (longitude and latitude) in the order in which they appear in the plink file, and the datapath.outer file which specifies the coordinates defining the polygon in which eems inferences take place. These coordinates were the same ones used for our gradient forest analysis.

Before combining the outputs of the five eems runs, we first confirmed that the MCMC chains were properly mixed. To do this, we plotted the chains using the R script MixingPlot.R. This script was used to produce Extended Data Figure S11 in our manuscript, and shows that chains are well mixed.

Finally, we integrate the 5 runs to produce a plot of the effective migration surface with the R script EemsPlot.R. Note that for this script to work correctly, eems must be installed (i.e. it is not a stand-alone R package that can be pulled from a remote repository such as CRAN), as well as geos and gdal. We used versions 3.6.2 and 2.3.0 for geos and gdal, respectively. We customized the eems.colors variable in this module so as to produce the colors schema in Figure 1D of our manuscript. 
 
