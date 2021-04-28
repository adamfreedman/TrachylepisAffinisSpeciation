# Variant filtering
## LD filtering
RAD-seq data sets typically produce more than one variant call per RAD locus. Given that these loci are relatively short, single nucleotide variants (SNVs) found on the same locus will be in strong linkage disequilibrium (LD). Inclusion of SNVs in strong LD has the potential to distort evolutionary inferences. Therefore, our first step in filtering the genotypes generated with [STACKS](https://catchenlab.life.illinois.edu/stacks/) is to downsample our variant call format (vcf) genotypes by randomly selecting one polymorphic site per RAD locus. We do this with a custom python script:

```bash
SampleOneSnpPerRadLocusFromVcf.py trachylepis_RAD.vcf
```
which produces the new, downsampled file *oneperrad_trachylepis_RAD.vcf*.

