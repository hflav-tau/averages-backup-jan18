# From TauFit/readpdg_constrained_aleph_hcorr.log
# BR_TauToE = (17.831 +- 0.0398) %
# BR_TauToPiNu =(10.8302 +- 0.0512) %
# BR_TauToKNu = (0.6971 +- 0.0096) %
# Corr between tau- --> e- nubar(e) nu(tau) and tau- --> pi- nu(tau) =   0.000961
# Corr between tau- --> e- nubar(e) nu(tau) and tau- --> K- nu(tau) =   0.041556
# Corr between tau- --> pi- nu(tau) and tau- --> K- nu(tau) =  -0.004707

#####################
 
# If the only error contributing to gtaugmu_pi and gtaugmu_e were from BR_TauToPiNu and BR_TauToE,
# then the Correlation co-efficient between gtaugmu_pi and gtaugmu_e would be also  0.000961. But it is not..

# From gtaugmu_pi/log/gtaugmu_pi.log 

# gtaugmu_pi = 0.996564 +- 0.00299738 out of which 0.00236167 is the error-component due to measured BR_TauToPiNu

# From gtaugmu_e/log/gtaugmu_e.log

# gtaugmu_e = 1.00105  +- 0.00207702 out of which 0.00112125 is the error-component due to measured BR_TauToE

# So, the covariance  between gtaugmu_pi and gtaugmu_e =  (0.000961 * 0.00236167 * 0.00112125) = 2.54475e-09
# and correlation is 100 * ( 0.000961 * 0.00236167 * 0.00112125 ) / ( 0.00299738 * 0.00207702 ) = 0.0408755 %
#####################

# If the only error contributing to gtaugmu_K and gtaugmu_e were from BR_TauToKNu and BR_TauToE,
# then the Correlation co-efficient between gtaugmu_K and gtaugmu_e would be also  0.041556. But it is not..

# From gtaugmu_k/log/gtaugmu_k.log

# gtaugmu_k = 0.986033 +- 0.007233 out of which 0.00680697 is the error-component due to measured BR_TauToKNu

# From gtaugmu_e/log/gtaugmu_e.log

# gtaugmu_e = 1.00105  +- 0.00207702 out of which 0.00112125 is the error-component due to measured BR_TauToE

# So, the covariance  between gtaugmu_K and gtaugmu_e =  (0.041556 * 0.00680697 * 0.00112125) = 3.17168e-07
# and the correlation is 100 * ( 0.041556 * 0.00680697 * 0.00112125) / (0.007233 * 0.00207702) =  2.11121 %
#####################

# If the only error contributing to gtaugmu_pi and gtaugmu_K were from BR_TauToPiNu and BR_TauToKNu,
# then the Correlation co-efficient between gtaugmu_pi and gtaugmu_K would be also -0.004707. But it is not..

# From gtaugmu_pi/log/gtaugmu_pi.log 

# gtaugmu_pi = 0.996564 +- 0.00299738 out of which 0.00236167 is the error-component due to measured BR_TauToPiNu

# From gtaugmu_k/log/gtaugmu_k.log

# gtaugmu_k = 0.986033 +- 0.007233 out of which 0.00680697 is the error-component due to measured BR_TauToKNu

# So, the covariance between gtaugmu_pi and gtaugmu_k =  (-0.004707 * 0.00236167 * 0.00680697) = -7.56689e-08
# and the correlation is 100 * (-0.004707 * 0.00236167 * 0.00680697) / (0.00299738 * 0.007233) = -0.349025

