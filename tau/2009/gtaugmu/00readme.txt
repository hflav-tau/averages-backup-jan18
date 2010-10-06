# From TauFit/readpdg_constrained_aleph_hcorr.log
# BR_TauToPiNu =(10.8302 +- 0.0512) %
# BR_TauToKNu = (0.6971 +- 0.0096)%
# Corr between tau- --> pi- nu(tau) and tau- --> K- nu(tau) =  -0.004707
 
# If the only error contributing to gtaugmu_pi and gtaugmu_K were from BR_TauToPiNu and BR_TauToKNu,
# then the Correlation co-efficient between gtaugmu_pi and gtaugmu_K would be also -0.004707. But it is not..

# From gtaugmu_pi/log/gtaugmu_pi.log 

# gtaugmu_pi = 0.996564 +- 0.00299738 out of which 0.00236167 is the error-component due to measured BR_TauToPiNu

# From gtaugmu_k/log/gtaugmu_k.log

# gtaugmu_k = 0.986033 +- 0.007233 out of which 0.00680697 is the error-component due to measured BR_TauToKNu

# So, Correlation Co-efficient between gtaugmu_pi and gtaugmu_k => (-0.004707 * 0.00236167 * 0.00680697) / (0.00299738 * 0.007233) = -0.00349025

# Inputs to 
# ../../../combos/swagato_example/Common/average.perl 
# are
# 0.996564
# 0.00299738*0.00299738 = 8.98429e-06
# 0.986033
# 0.007233*0.007233 = 5.23163e-05
# -0.004707 * 0.00236167 * 0.00680697 = -7.56689e-08

# So, run it by simply typing | source 00readme.txt |

../../../combos/swagato_example/Common/average.perl 0.996564 8.98429e-06 0.986033 5.23163e-05 -7.56689e-08

# x1 = 0.996564 +- 0.00299738052305676 x2 = 0.986033 +- 0.00723300076040367 rho = -0.00349025396370819
# w1 = 0.8525683534374 w2 = 0.1474316465626 w1+w2 = 1
# <x> = 0.995011397330049 +- 0.00276560397410467

