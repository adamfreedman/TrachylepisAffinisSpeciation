from collections import defaultdict
misdict = defaultdict(int)
fopen=open('oneperrad_trachylepis_RAD.vcf','r')
fout=open('trachylepis_RAD_oneperrad_missingstats_by_ind.txt','w')
fout.write('id\tmiscount\tmisfreq\n')

for line in fopen:
    if '#CHROM' in line:
        fields = line.strip().split('\t')
    elif line[0] == '#' and 'CHROM' not in line:
        pass
    else:
        zerocount = 0
        linedict =  dict(zip(fields,line.strip().split()))
        for sampleid in fields[8:]:
            if linedict[sampleid].split(':')[0] == './.':
                misdict[sampleid]+=1
for sampleid in misdict:
    fout.write('%s\t%s\t%s\n' % (sampleid,misdict[sampleid],misdict[sampleid]/float(56055)))

fout.close()
