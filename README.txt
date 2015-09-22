An octave function, will fit lineshape data, along with dissociation constant, isomerization rate, and bound chemical shift values to determine microscopic rate constants kab, kba, kbc, kcb, kcd, and kdc as described in Chapter III.
Data are 2D extracted hsqcs comprising lineshapes that will be summed over nitrogen to generate a 1D proton lineshape to fit. Data must be stored in:
  ./dat/MUT/FGPconcentration_CYPconcentration/hsqc_RES#_RESconformation.txt 
e.g.: ./data/WT/1000_50/hsqc_6_trans.txt
Data folder must contain data for residues 6 and 7, cis and trans, and also contain an ‘hsqc_11_ref.txt’ file that contains reference data for GSW peptide reference peak
Output file is in Octave format, with parameters as listed below

function [] = Fit_Lineshapes(MUT, KD, kiso, freq_bound, F_Out, Cyp_Conc, FGP_Conc, iterations, acquisition_time, acquisition_points, zero_fill)

MUT = string protein/mutant name (default, “WT”)
KD = Dissociation constant to fit (M), if set to 0, will not use in fit
kiso = isomerization rate, at 1 mM peptide, 20 uM protein (s^-1), if set to 0, will not use in fit
freq_bound = vector of proton frequencies, in Hz, of bound residues, if set to 0, will not use in fit
F_out = output file (Default “temp.out”)
Cyp_Conc = Cyp concentrations to use, in uM, first value must be free peptide, with Cyp_Conc==0 (default = [0, 5, 10, 20, 5, 10, 20, 50, 100]) 
FGP_Conc = Peptide concentrations to use, in uM (default = [1000, 500, 500, 500, 1000, 1000, 1000, 1000, 1000])
iterations = number of fits to perform (default = 200)
acquisition_time = collected FID time (in s, default 1)
acquisition_points = number of points acquired in FID (default 14045)
zero_fill = if zero filling is applied, number of points (default 2^16)

Output File Parameters (N = number of iterations)
PARS_OUT = fit parameters, a 10xN matrix of parameters, with each column corresponding to [kab; kba; kbc; kcb; kcd; kdc; w_6_bound_trans, w_6_bound_cis, w_7_bound_trans, w_7_bound_cis]
r2_OUT = r^2 values for each fit
f_OUT = XxN matrix of fits, where X is the number of data points fit
full_data_vector = linearized 1D proton data set to which fit is performed
data_freq_linear = frequencies corresponding to each data point in full_data_vector
Populations_OUT = 5xN matrix of the populations occupied given 1 mM peptide and 20 uM Cyp, given the fit parameters, with each column corresponding to [free_trans; bound_trans; bound_cis; free_cis; free_Cyp]

