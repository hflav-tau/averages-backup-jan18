# From TauTo1Prong/log/average.log
# BR_TauToPiNu =  (10.7878 +- 0.0556) %
# BR_TauToKNu = (0.6969439 +- 0.00975197 )%
# with Correlation Co-efficient = 0.0846088
 
# If the only error contributing to gtaugmu_pi and gtaugmu_K were from BR_TauToPiNu and BR_TauToKNu,
# then the Correlation co-efficient between gtaugmu_pi and gtaugmu_K would be also 0.0846088. But it is not..

# From gtaugmu_pi/log/gtaugmu_pi.log 

# gtaugmu_pi = 0.994612 +- 0.00316036 out of which 0.00256966 is the error-component due to measured BR_TauToPiNu

# From gtaugmu_k/log/gtaugmu_k.log

# gtaugmu_k = 0.985923 +- 0.00733221 out of which 0.00691551 is the error-component due to measured BR_TauToKNu

# So, Correlation Co-efficient between gtaugmu_pi and gtaugmu_k => (0.0846088 * 0.00256966 * 0.00691551) / (0.00316036*0.00733221) = 0.0648849

# Inputs to 
# ../../../combos/swagato_example/Common/average.perl 
# are
# 0.994612
# 0.00316036*0.00316036 = 9.98788e-06
# 0.985923
# 0.00733221*0.00733221 = 5.37613e-05
# 0.0846088 * 0.00256966 * 0.00691551 = 1.50354e-06

# So, run it by simply typing | source 00readme.txt |

../../../combos/swagato_example/Common/average.perl 0.994612 9.98788e-06 0.985923 5.37613e-05 1.50354e-06

# x1 = 0.994612 +- 0.00316036073890308 x2 = 0.985923 +- 0.00733220976241133 rho = 0.0648848713186181
# w1 = 0.86032191840585 w2 = 0.13967808159415 w1+w2 = 1
# <x> = 0.993398337149028 +- 0.00296695191487956
