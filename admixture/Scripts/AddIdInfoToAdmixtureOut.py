import sys
idsin=open(sys.argv[1],'r')
admixin=open(sys.argv[2],'r')

fout=open('withids_'+sys.argv[2],'w')
numids=int(sys.argv[3])

for i in range(numids):
    idline=idsin.readline().strip()
    props=admixin.readline().strip()
    fout.write('%s %s\n' % (idline,props))


fout.close()
