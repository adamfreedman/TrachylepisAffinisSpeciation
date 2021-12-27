from numpy.random import rand
import argparse
from collections import defaultdict

def ExtractAlleleFrequenciesFromCounts(dadialleles_pop1,dadialleles_pop2,variantcount1,variantcount2):
    pop1 = variantcount1/float(dadialleles_pop1)
    pop2 = variantcount2/float(dadialleles_pop2)
    return pop1,pop2

def GenerateGenotypesFromFrequencies(freq1,freq2,gtypes1,gtypes2):
    gtypes = defaultdict(list)
    for i in range(gtypes1):
       alleles = []
       for j in range(2):
            if freq1 <= rand():
                alleles.append('1')
            else:
                alleles.append('0')
       alleles.sort()
       gtypes['pop1'].append('/'.join(alleles))

    for i in range(gtypes2):
       alleles = []
       for j in range(2):
            if freq2 <= rand():
                alleles.append('1')
            else:
                alleles.append('0')
       alleles.sort()
       gtypes['pop2'].append('/'.join(alleles))
    
    return gtypes        

def CalcMafFromGtypes(gtype_dict):
    alleles = ''
    for pop in gtype_dict:
        for gtype in gtype_dict[pop]:
            alleles+=gtype.replace('/','')
    maf=min(alleles.count('0'),alleles.count('1'))/float(len(alleles))
    return maf

def BuildGtypeOutString(gtype_dict):
    gtype_string_dict = {'0/0':'0/0:56:56,0:.,77.63,.','0/1':'0/1:46:22,24:.,63.77,.','1/1':'1/1:23:0,23:.,31.88,.'}
    gtype_list = []
    for gtype in gtype_dict['pop1']:
        gtype_list.append(gtype_string_dict[gtype])
    for gtype in gtype_dict['pop2']:
        gtype_list.append(gtype_string_dict[gtype])
    return '\t'.join(gtype_list)  

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='takes table of minor allele frequencies sampled from 2DSFS to produce a fake vcf file')
    parser.add_argument('-i','--sfs-infile',dest='sfs',type=str,help='filename of tab-separated allele counts')
    parser.add_argument('-a1','--pop1_downprojected_alleles',dest='a1',type=int,help='number of alleles dadi downprojected population 1')
    parser.add_argument('-a2','--pop2_downprojected_alleles',dest='a2',type=int,help='number of alleles dadi downprojected population 2')
    parser.add_argument('-d1','--pop1_ndiploids',dest='d1',type=int,help='number of diploid genotypes to sample from frequencies')
    parser.add_argument('-d2','--pop2_ndiploids',dest='d2',type=int,help='number of diploid genotypes to sample from frequencies')
    parser.add_argument('-vout','--vcfout',dest='vcfout',type=str,help='name of output vcf file')
    opts = parser.parse_args()
    ids = []
    for i in range(opts.d1):
        ids.append('ECO%s' % str(i+1))
    for j in range(opts.d2):
        ids.append('FOR%s' % str(j+1))
    idstring = '\t'.join(ids)
    header = '##fileformat=VCFv4.0\n##fileDate=20161107\n##source="Stacks v1.32"\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">\n##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihood">\n'

    fout = open(opts.vcfout,'w')
    fout.write(header)
    fout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % idstring)

    sfsin = open(opts.sfs,'r')
    counter = 0
    labels =  sfsin.readline().strip().split('\t')
    for line in sfsin:
        pop1count,pop2count = [int(i) for i in line.strip().split('\t')]
        pop1freq, pop2freq = ExtractAlleleFrequenciesFromCounts(opts.a1,opts.a2,pop1count,pop2count)

        if min((pop1count+pop2count)/float(opts.a1+opts.a2),1-(pop1count+pop2count)/float(opts.a1+opts.a2))>=0.05:
            counter+=1
            print 'counter ==', counter
            mafpass = 0
            while mafpass == 0:
                gtypes = GenerateGenotypesFromFrequencies(pop1freq,pop2freq,opts.d1,opts.d2)
                if CalcMafFromGtypes(gtypes) >= 0.05:
                    mafpass += 1
                    gtypes_formatted = BuildGtypeOutString(gtypes)
                    fout.write('un\t%s\t%s\tC\tA\t.\tPASS\t.\tGT:DP:AD:GL\t%s\n' % (counter,counter,gtypes_formatted))


fout.close()
                
