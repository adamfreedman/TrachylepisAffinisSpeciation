import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
from datetime import datetime
import matplotlib.pyplot as plt

'''
usage: python dadi_2D_02_second_optimizations.py

Requires the Models_2D.py script to be in same working directory. This is where
all the population model functions are stored for this script. 

Script will perform optimizations from multiple starting points using a
1-fold perturbed set of USER SELECTED starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate based on ll or AIC.

The outputs from here will be labeled according to model and a user
selected prefix name, and are written to the working directory:

"Round3_[prefix]_[model_name]_optimized.txt"

example output file from one model:

Model	param_set	Replicate	log-likelihood	theta	AIC	optimized_params
Divergence with no migration	parameter set = [nu1, nu2, T]	1	-443.41	317.15	892.82	0.6795	0.7315	0.141	
Divergence with no migration	parameter set = [nu1, nu2, T]	2	-179.8	141.7	365.6	4.3798	1.6121	0.7995	
Divergence with no migration	parameter set = [nu1, nu2, T]	3	-241.28	116.99	488.56	3.1119	17.7282	1.1496	
Divergence with no migration	parameter set = [nu1, nu2, T]	4	-1763.04	54.2	3532.08	9.0388	0.2684	5.3911	
Divergence with no migration	parameter set = [nu1, nu2, T]	5	-177.49	150.28	360.98	2.7819	2.6884	0.8431	

*Note the params are written to as many columns as there are values, without any column heading.



Requires user to edit sections of code marked with #**************

You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 


############################################
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-Matplotlib
-dadi
############################################

Dan Portik
daniel.portik@uta.edu
April 2017
'''
t_begin = datetime.now()

#===========================================================================
#get snps file 

#**************
snps1 = "dadi_input_20170915_141444.tsv"

#Create python dictionary from snps file
dd1 = dadi.Misc.make_data_dict(snps1)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["ECO", "FOR"]
#projection sizes, in ALLELES not individuals
proj_1 = [91,215]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
#create function to run models with user input parameters (presumably from best runs)
#run optimization 'x' times on perturbed best starting params

def Two_Pop_Models(pts, fs, outfile, reps, y, model_name, params):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round3_{0}_{1}_optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(reps) + int(1)
        
    if model_name == "no_divergence":
        #####################################
        #No Divergence
        #####################################
        print "---------------------------------------------------"
        print "No Divergence",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_divergence
        params = []
        print "parameter set = [none]"
        
        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("No divergence model"+'\t')
            fh_out.write("parameter set = [none]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #simulate the model with the optimized parameters
            sim_model = func_exec(params, fs.sample_sizes, pts)
            
            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))
            
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'

        
    elif model_name == "no_mig":
        fh_out = open(outname, 'a')
        #####################################
        # Divergence with no migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0]
        upper_bound = [30, 30, 10]
        print "parameter set = [nu1, nu2, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*3)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 10]
        print "parameter set = [nu1, nu2, m, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*4)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 20, 20, 10]
        print "parameter set = [nu1, nu2, m12, m21, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_sym_mig":
        fh_out = open(outname, 'a') 
        #####################################
        #Divergence with ancient symmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "sec_contact_sym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 10, 10]
        print "parameter set = [nu1, nu2, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*5)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sec_contact_asym_mig":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1, nu2, m12, m21, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact"+'\t')
            fh_out.write("parameter set = [nu1, nu2, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    ##########################################################################
    #Similar models but with size change allowed
    ##########################################################################
    
    elif model_name == "no_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with no migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with no migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.no_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with no migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
        
    elif model_name == "sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with symmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with symmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with symmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        

    elif model_name == "asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with asymmetric migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with asymmetric migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
 
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with asymmetric migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
 
 
    elif model_name == "anc_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient symmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient symmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient symmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

        
    elif model_name == "anc_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence with ancient asymmetrical migration, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence with ancient asymmetrical migration, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.anc_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence with ancient asymmetrical migration, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'


    elif model_name == "sec_contact_sym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and symmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and symmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_sym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and symmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*7)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'

           
    elif model_name == "sec_contact_asym_mig_size":
        fh_out = open(outname, 'a')
        #####################################
        #Divergence and asymmetrical secondary contact, size change
        #####################################
        print "---------------------------------------------------"
        print "Divergence and asymmetrical secondary contact, size change",'\n','\n'

        #first call a predefined model
        model_call = Models_2D.sec_contact_asym_mig_size

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 10, 10]
        print "parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"
    
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Divergence and asymmetrical secondary contact, size change"+'\t')
            fh_out.write("parameter set = [nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "base parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*8)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()
        print "---------------------------------------------------", '\n'
        
    
        

#======================================================================================
# Finally, execute model with appropriate arguments
# Two_Pop_Models(pts, fs, outfile, reps, maxiter, model_name, params):

# pts:  grid choice (list of three numbers, ex. [20,30,40]

# fs:  spectrum object name

# outfile:  prefix for output naming -> "Round3_{0}_{1}_optimized.txt".format(outfile,model_name)

# reps:  integer to control number of replicates, ex. 10

# maxiter:  max number of iterations per optimization step (not intuitive! see dadi user group)

# model_name:  "no_divergence", "no_mig", "sym_mig", "asym_mig", "anc_sym_mig", "anc_asym_mig",
#        "sec_contact_sym_mig", "sec_contact_asym_mig", "no_mig_size", "sym_mig_size",
#        "asym_mig_size", "anc_sym_mig_size", "anc_asym_mig_size", "sec_contact_sym_mig_size",
#        "sec_contact_asym_mig_size"

# params:  list of best parameters to start optimizations from



#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script

#************** "no_divergence"
#leave blank, no parameters
no_divergence_params = []

#************** "no_mig"
#3 Values
no_mig_params = [0.5012,1.9103,0.4392] #ok

#************** "sym_mig"
#4 Values
sym_mig_params = [1.3293,4.1482,0.0626,2.9249] #ok

#************** "asym_mig"
#5 Values
asym_mig_params = [1.4353,4.0388,0.0345,0.0817,2.9302] #ok

#************** "anc_sym_mig"
#5 Values
anc_sym_mig_params = [1.1453,4.4656,0.0855,2.0577,0.0281] #ok

#************** "anc_asym_mig"
#6 Values
anc_asym_mig_params = [0.9427,3.653,0.107,0.1693,1.2648,0.0391] #ok

#************** "sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [1.2837,3.9,0.0638,0.3196,2.4557] #ok

#************** "sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.9781,2.7379,0.0733,0.2236,1.277,0.1317] #ok

#************** "no_mig_size"
#6 Values
no_mig_size_params = [0.1353,6.3531,2.7628,1.9423,0.1485,0.0637] #ok

#************** "sym_mig_size"
#7 Values
sym_mig_size_params = [0.2826,1.4253,6.7603,2.8593,0.081,0.9944,0.0634] #ok

#************** "asym_mig_size"
#8 Values
asym_mig_size_params = [0.4254,0.0417,1.64,4.4519,0.0747,0.0697,0.8194,2.227] #ok

#************** "anc_sym_mig_size"
#7 Values
anc_sym_mig_size_params = [0.1687,1.2689,12.2287,8.4463,0.3529,0.4814,0.0364] #ok

#************** "anc_asym_mig_size"
#8 Values
anc_asym_mig_size_params = [0.1574,0.6818,2.8193,10.809,0.2986,0.4259,0.4569,0.0237] #ok

#************** "sec_contact_sym_mig_size"
#7 Values
sec_contact_sym_mig_size_params = [0.2581,0.6261,2.748,3.471,0.1797,0.4097,0.055] #ok

#************** "sec_contact_asym_mig_size"
#8 Values
sec_contact_asym_mig_size_params = [0.4982,1.2158,0.7572,2.1861,0.1006,0.1073,0.3594,0.9744] #ok




#===========================================================================
#**************
#Input some of the basic reusable arguments here
pts = [50,60,70]
fs = fs_1
outfile = "EcotoneVsForest_coarsegrid"
reps = int(100)
maxiter = int(50)


#===========================================================================
# Two_Pop_Models(pts, fs, outfile, reps, maxiter, model_name, params):

# Here it is set up to call each model one by one sequentially, which could finish relatively quickly.
# If it takes too long, create multiple verisions of this script, block out some models (use a hash or delete),
# and execute one version for every core you have available. It will greatly speed up these steps,
# and sometimes if extrapolations fail the script will crash too and this could prevent it from
# happening too many times.
# There are 15 models to test here. If you do nothing, they will all run.

Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_divergence", no_divergence_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "no_mig_size", no_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sym_mig_size", sym_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "asym_mig_size", asym_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_sym_mig_size", anc_sym_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "anc_asym_mig_size", anc_asym_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig_size", sec_contact_sym_mig_size_params)
Two_Pop_Models(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig_size", sec_contact_asym_mig_size_params)


#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================

