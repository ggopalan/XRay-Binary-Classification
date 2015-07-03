XRay-Binary-Classification
========================

Contained in this repository are the data files and R code used for the analyses in 
"Classifying X-ray Binaries: A Probabilistic Approach" by G. Gopalan, S. Vrtilek, and 
L. Bornn. Our hope is that one may use this code to both recreate the analysis in this 
paper and apply it to other data sets.

1) Data Sets: There are four main folders of CCI data used in the paper:
	i) The 24 training systems are in the folder "Training".
	ii) The 6 validation systems from Table 1 are in "Validation".
	iii) The 6 burster systems from Table 2 are in "Burster_Test".
	iv) The 5 unclassified systems from Table 3 are in "Unclass".
	
Note on the formatting of input data: the first column is the system name, the second is the second is the system type, the third is the date, the fourth is the intensity, the fifth is color 1, and the sixth is color 2.

2) Main Code: The flow of the code is as follows.
    i) First, "observational_study_design.R" reads in the CCI files in "Training" and 
    subsamples the data to create a training set as  discussed in Section 2. The 
    resultant file is output to "dat.RData".
    ii) Second, "rcpp_bbh_inv_samples.r" is run to generate the samples necessary for the
    MCMC algorithm we employ. The resultant output is stored in "rcpp_bbh_inverse.RData"
    iii) Third, "bbh_ess_sampler.r" is run to sample from the posterior distribution of 
    model parameters and latent variables using elliptical slice sampling (Murray, Adams,
    and MacKay 2010). The resultant output is stored in "bbh_ess_10000reps.RData"
    iv) Fourth, "predictive_distribution_sampler.r" is run to generate the posterior predictive samples
    of compact object type of all the observations in a file specified in the first command line argument (where the third through fifth coumns contain CCI values) and the output file name is specified in the second command line argument. The resultant output is stored in "system_name.RData".
    
3) Output: The ".RData" files in the "Validation_Output" folder store the posterior 
   predictive draws of the compact object type for all observations in the system (denoted as Y_pred in the paper)      and are located in the Validation_Output
   folder. 
    
NOTES:
1) This code was run on the Harvard Odyssey supercomputer and so it will take substantially
longer on a personal computer; if possible the code should be run on a supercomputer.
2) More instructions for running the code are found in the comments within each R file.
