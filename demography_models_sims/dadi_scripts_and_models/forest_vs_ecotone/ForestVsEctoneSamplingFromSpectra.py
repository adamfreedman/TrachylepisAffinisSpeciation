import dadi
import numpy
from Models_2D import anc_asym_mig_size

def Create2dsfsIndexingArray(sfs):
    rows,columns =  sfs.shape
    indexing_list = []
    for i in range(rows):
        for j in range(columns):
            indexing_list.append([i+1,j+1])
    index_array = numpy.array(indexing_list)
    return index_array

def SampleSfs(indexarray,probs,nsamples=56055):
    samples =  numpy.random.choice(19872,size=nsamples,p=probs)
    allele_combo = []
    for sample in samples:
        allele_combo.append(indexarray[sample,])
    return allele_combo 
### 1. load empirical data ###
snps1 = "/n/holylfs/LABS/informatics/adamf/trachylpeis_rad/dadi/dadi_forVseco_input_20170915_141444.tsv"

### 2.make dadi data dictionary from empirical data file ###
dd1 = dadi.Misc.make_data_dict(snps1)

### 3. provide populations labels ###
pop_ids=["ECO", "FOR"]

### 4. provide downward projection of number of alleles assayed ###
proj_1 = [91,215]

### 5. create 2D SFS from data dictionary, ids, and projection data ###
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)

### 6. normalize 2D observation counts by sum of 2D SFS to get approximate probabilities of occurence
fs_probs = fs_1/numpy.sum(fs_1)
fs_probs_1d =  fs_probs.reshape(19872,)
fs_probs_1d = numpy.ma.filled(fs_probs_1d,0)
index_array = Create2dsfsIndexingArray(fs_1)

### sample sfs
samples = SampleSfs(index_array,fs_probs_1d)
fout = open('sim_allele_counts.txt','w')
fout.write('ECO\tFOR\n')
for sample in samples:
   fout.write('%s\t%s\n' % (sample[0],sample[1]))

fout.close()

### MODEL 
best_model_params = [0.1009, 0.5933, 29.9198, 8.7145, 1.1616, 0.4132, 5.3655, 0.0269]
#pts = [100,110,120]
pts=120
modelsfs = anc_asym_mig_size(best_model_params, proj_1, pts)
modelprobs = modelsfs/numpy.sum(modelsfs)
modelprobs = modelprobs.reshape(19872,)
modelprobs = numpy.ma.filled(modelprobs,0)
model_samples = SampleSfs(index_array,modelprobs)
modelout = open('modelsim_alelle_counts.txt','w')
modelout.write('ECO\tFOR\n')
for sample in model_samples:
   modelout.write('%s\t%s\n' % (sample[0],sample[1]))
modelout.close()
