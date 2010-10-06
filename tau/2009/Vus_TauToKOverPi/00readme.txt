#
# This file is outdated 
# [SwB Tue Oct  5 04:35:16 PDT 2010]  
# Now, just input to vus.input this result:
# G(K- nu(tau))/G(pi- nu(tau)) =   0.064362 +-   0.000940 
# from TauFit/readpdg_constrained_aleph_hcorr.log
#
# From TauTo1Prong/log/average.log
# BR_TauToPiNu =  (10.7878 +- 0.0556) %
# BR_TauToKNu = (0.6969439 +- 0.00975197 )%
# with Correlation Co-efficient = 0.0846088
 
# So, BR_TauToKNu/BR_TauToKNu = 0.0646048 and the error is Delta,

# f = a/b => df/da  = 1/b; df/fb = -a/b**2

# sigma_f**2 = (sigma_a/b)**2 + (a sigma_b/b**2)**2 - 2 (a/b**3) Cov(a,b)

# But Cov(a,b) = sigma_a sigma_b rho_ab

# So, (sigma_f/f)**2 = sigma_f**2 * (b**2/a**2) 
#  =  (b**2/a**2)  * [ (sigma_a/b)**2 + (a sigma_b/b**2)**2 - 2 (a/b**3) Cov(a,b) ]
#  = (sigma_a/a) **2 + (sigma_b/b)**2 - 2 (sigma_a/a) * (sigma_b/b)* rho_ab

# Here sigma_a/a = 0.0556/10.7878 = 0.00515397
# and  sigma_b/b = 0.00975197/0.6969439 = 0.0139925

# So, (sigma_f/f)**2 =  0.00515397**2 +  0.0139925**2 - 2 * 0.00515397 *  0.0139925 * 0.0846088 = 0.00021015

# So, sigma_f = 0.000936547

# =>  RK/Rpi = 0.0646048 +- 0.000936547

###################################
# Now run

../../../Common/epc -s0 vus.input | tee vus.log
