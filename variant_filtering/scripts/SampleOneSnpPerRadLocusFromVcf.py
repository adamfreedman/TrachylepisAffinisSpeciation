import sys
from collections import defaultdict
from os.path import basename
genotype_dict=defaultdict(list)
from numpy.random import randint

vcfin=open(sys.argv[1],'r')
vcfout=open('oneperrad_'+ basename(sys.argv[1]),'w')
for line in vcfin:
    if line[0]=='#':
        vcfout.write(line)
    else:
       linelist=line.strip().split()
       genotype_dict[int(linelist[2])].append(line)
          
keys=genotype_dict.keys()
keys.sort()    
print('dictionary build completed')    
   
counter=0    
for i in range(len(keys)):
    counter+=1
    if counter%1000==0:
        print('processing rad locus...%s' % counter)
    randselect=randint(0,len(genotype_dict[keys[i]]))
    vcfout.write(genotype_dict[keys[i]][randselect])
    
vcfout.close() 
