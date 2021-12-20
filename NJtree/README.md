# Neighbor-joining tree
As a method for representing both population structure and individual-level differentiation, we generated a neighbor-joining tree from a matrix of genetic distances in the form of 1 minus identity-by-state (IBS). This matrix is calculated with plink as follows:
```bash
plink --file oneperrad_trachylepis --distance square 1-ibs
```
where oneperrad_trachylepis is the prefix of the set of plink files generated from the unfiltered genotype set with ped,map,and fam suffixes, i.e. uncompressed plink format.

After generating the distance matrix, we use treeplot.R to generate a neighbor-joining tree, using the distance matrix file plink.mdist and wcols_ids_renumbered.txt, which is a table of sample ids, corresponding populaton ids, and colors to use in coloring sample labels by population id. As population labels confirmed a structuring of the tree largely concordant with results obtained with ADMIXTURE results for unfiltered genotypes for K=7--the number of ancestral populations with the lowest cross-validation error--we redrew the plot without labels, and manually drew colored polygons around each clade. Geographic populations that showed evidence of admixture between ancestral populations were bi-colored according to the colors used for admixture plots. The resulting figure corresponds to Supplementary Data Figure S2 in our manuscript. 

