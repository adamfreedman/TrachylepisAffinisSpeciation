import sys
ids=open(sys.argv[1],'r')

iddict={}
for line in ids:
    linelist=line.strip().split()
    iddict[linelist[0]]=linelist[1]
    
pcs=open(sys.argv[2],'r')

fout=open('wpoplabels_'+sys.argv[2],'w')
fout.write('id\tpop\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n')
for line in pcs:
    linelist=line.strip().split()
    id=linelist[0]
    fout.write('%s\t%s\t%s\n' % (id,iddict[id],"\t".join(linelist[1:])))
    
fout.close()  
