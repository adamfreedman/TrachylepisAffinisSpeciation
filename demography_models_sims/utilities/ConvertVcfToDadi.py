from collections import defaultdict
import argparse

def BuildStrataDict(stratafile):
    # strata file has a header line
    strata = open(stratafile,'r')
    strata.readline()
    strata_dict = {}
    for line in strata:
        id,population = line.strip().split('\t')
        strata_dict[id] = population
    return strata_dict

def ParseGenotypes(vcflinedict,strata_dict):
    allele1_dict = defaultdict(int)
    allele2_dict = defaultdict(int)
    for sampleid in strata_dict.keys():
        gtype = vcflinedict[sampleid].split(':')[0]
        allele1_dict[strata_dict[sampleid]] += gtype.count('0') 
        allele2_dict[strata_dict[sampleid]] += gtype.count('1')
        
    return allele1_dict,allele2_dict

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert Stacks-generated vcf genotypes to dadi input")
    parser.add_argument('-strata','--strata-file',dest='strata',type=str,help='strata file, col1=id,col2=pop')
    parser.add_argument('-vcf','--vcf-file',dest='vcf',type=str,help='vcf genotypes infile')
    parser.add_argument('-o','--outfile',dest='out',type=str,help='name of dadi input file exported')
    opts = parser.parse_args()
    
    strata_dict = BuildStrataDict(opts.strata)
    pops = list(set(strata_dict.values()))
    pops.sort() 
    fout = open(opts.out,'w')
    out_header = 'IN_GROUP\tOUT_GROUP\tAllele1\t%s\tAllele2\t%s\tMARKERS\n' % ('\t'.join(pops),'\t'.join(pops))
    fout.write(out_header)
    
    fopen = open(opts.vcf,'r')
    for line in fopen:
        if line[:6] == '#CHROM': 
            fields = line[1:].strip().split('\t')
            print(fields)
        elif line[0] =='#':
            pass    
        else:
            linedict = dict(zip(fields,line.strip().split()))
            markerstring = '%s__%s__%s' % (linedict['CHROM'],linedict['ID'],linedict['POS'])    
            allele1_dict,allele2_dict = ParseGenotypes(linedict,strata_dict)
            allele1_counts = [str(allele1_dict[pop]) for pop in pops]
            allele2_counts = [str(allele2_dict[pop]) for pop in pops] 
            fout.write('---\t---\t%s\t%s\t%s\t%s\t%s\n' % (linedict['REF'],'\t'.join(allele1_counts),linedict['ALT'],'\t'.join(allele2_counts),markerstring))  
    
    fout.close()
    

