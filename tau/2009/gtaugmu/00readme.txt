# From TauTo1Prong/log/average.log
# BR_TauToPiNu =  (10.7878 +- 0.0556) %
# BR_TauToKNu = (0.6969439 +- 0.00975197 )%
# with Correlation Co-efficient = 0.0846088
 
# If the only error contributing to gtaugmu_pi and gtaugmu_K were from BR_TauToPiNu and BR_TauToKNu,
# then the Correlation co-efficient between gtaugmu_pi and gtaugmu_K would be also 0.0846088. But it is not..

# From gtaugmu_pi/log/gtaugmu_pi.log 

# gtaugmu_pi = 0.99463 +- 0.00316042 out of which 0.0025697 is the error-component due to measured BR_TauToPiNu

# From gtaugmu_k/log/gtaugmu_k.log

# gtaugmu_k = 0.985941 +- 0.00733234 out of which 0.00691564 is the error-component due to measured BR_TauToKNu

# So, Correlation Co-efficient between gtaugmu_pi and gtaugmu_k => (0.0846088 * 0.0025697 * 0.00691564) / (0.00316042 * 0.00733234) = 0.0648848

# Inputs to 
# ../../../combos/swagato_example/Common/average.perl 
# are
# 0.99463
# 0.00316042*0.00316042 = 9.98825e-06
# 0.985941
# 0.00733234*0.00733234 = 5.37632e-05
# 0.0846088 * 0.0025697 * 0.00691564 = 1.50359e-06

# So, run it by simply typing | source 00readme.txt |

../../../combos/swagato_example/Common/average.perl 0.99463 9.98825e-06 0.985941 5.37632e-05 1.50359e-06 

# x1 = 0.99463 +- 0.00316041927598222 x2 = 0.985941 +- 0.00733233932657239 rho = 0.0648846806742611
# w1 = 0.860321640213966 w2 = 0.139678359786034 w1+w2 = 1
# <x> = 0.993416334731819 +- 0.00296700633768414

