#!/usr/local/bin/perl

# The BaBar Paper is http://arXiv.org/pdf/0707.2981

# Labels channels in tau- -> h-h-h+nu as 1=pipipi, 2=pikpi, 3=pikk, 4=kkk

# BR_1 = (8.83 +- 0.01 +- 0.13)e-2 
# BR_2 = (0.273 +- 0.002 + 0.009)e-2 
# BR_3 = (0.1346 +- 0.0010 +- 0.0036)e-2 
# BR_4 = (1.58 +- 0.13 +- 0.12)e-5

$BR_1 = 8.8337e-2;
$BR_2 = 0.27257e-2;
$BR_3 = 0.13461e-2;
$BR_4 = 1.5777e-5;

$stat_1 = 0.0074e-2;
$stat_2 = 0.0018e-2;
$stat_3 = 0.001e-2;
$stat_4 = 0.13e-5;

$syst_1 = 0.126724e-2;
$syst_2 = 0.00924406e-2;
$syst_3 = 0.00364131e-2;
$syst_4 = 0.12308e-5;

$sigma_1 = sqrt($stat_1**2 + $syst_1**2);
$sigma_2 = sqrt($stat_2**2 + $syst_2**2);
$sigma_3 = sqrt($stat_3**2 + $syst_3**2);
$sigma_4 = sqrt($stat_4**2 + $syst_4**2);

printf "%s %.4e %s %.4e %s %.4e %s %.4e\n", "BR_1 = ",$BR_1," stat_1 = ",$stat_1," syst_1 = ",$syst_1," sigma_1 = ",$sigma_1;
printf "%s %.4e %s %.4e %s %.4e %s %.4e\n", "BR_2 = ",$BR_2," stat_2 = ",$stat_2," syst_2 = ",$syst_2," sigma_2 = ",$sigma_2;
printf "%s %.4e %s %.4e %s %.4e %s %.4e\n", "BR_3 = ",$BR_3," stat_3 = ",$stat_3," syst_3 = ",$syst_3," sigma_3 = ",$sigma_3;
printf "%s %.4e %s %.4e %s %.4e %s %.4e\n", "BR_4 = ",$BR_4," stat_4 = ",$stat_4," syst_4 = ",$syst_4," sigma_4 = ",$sigma_4;

# The total correlation coefficients between channels are:
$rhotot_12 = 0.543535;
$rhotot_13 = 0.390346;
$rhotot_14 = 0.031469;
$rhotot_23 = 0.177495;
$rhotot_24 = 0.0931907;
$rhotot_34 = 0.0870484;

# The statistical correlation between channels:
$rhostat_12 = -0.129981;
$rhostat_13 =  0.00136523;
$rhostat_14 = -0.00186661;
$rhostat_23 = -0.150561;
$rhostat_24 =  0.00619404;
$rhostat_34 = -0.0948907;

# 
# The usual formula for Covariance_IJ =
#       rhotot_IJ * sigma_I * sigma_J = rhostat_IJ * stat_I * stat_J + rhosyst_IJ * syst_I * syst_J
# 

$rhosyst_12 = ($rhotot_12 * $sigma_1 * $sigma_2 - $rhostat_12 * $stat_1 * $stat_2) / ($syst_1 * $syst_2);
$rhosyst_13 = ($rhotot_13 * $sigma_1 * $sigma_3 - $rhostat_13 * $stat_1 * $stat_3) / ($syst_1 * $syst_3);
$rhosyst_14 = ($rhotot_14 * $sigma_1 * $sigma_4 - $rhostat_14 * $stat_1 * $stat_4) / ($syst_1 * $syst_4);
$rhosyst_23 = ($rhotot_23 * $sigma_2 * $sigma_3 - $rhostat_23 * $stat_2 * $stat_3) / ($syst_2 * $syst_3);
$rhosyst_24 = ($rhotot_24 * $sigma_2 * $sigma_4 - $rhostat_24 * $stat_2 * $stat_4) / ($syst_2 * $syst_4);
$rhosyst_34 = ($rhotot_34 * $sigma_3 * $sigma_4 - $rhostat_34 * $stat_3 * $stat_4) / ($syst_3 * $syst_4);

# Here are the actual numbers from the Supporting Document V19
$rhosyst_12 = 0.556623;
$rhosyst_13 = 0.406991;
$rhosyst_14 = 0.0466315;
$rhosyst_23 = 0.197039;
$rhosyst_24 = 0.138868;
$rhosyst_34 = 0.16347;

printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_12 = ",$rhostat_12," rhosyst_12 = ",$rhosyst_12," rhotot_12 = ",$rhotot_12;
printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_13 = ",$rhostat_13," rhosyst_13 = ",$rhosyst_13," rhotot_13 = ",$rhotot_13;
printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_14 = ",$rhostat_14," rhosyst_14 = ",$rhosyst_14," rhotot_14 = ",$rhotot_14;
printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_23 = ",$rhostat_23," rhosyst_23 = ",$rhosyst_23," rhotot_23 = ",$rhotot_23;
printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_24 = ",$rhostat_24," rhosyst_24 = ",$rhosyst_24," rhotot_24 = ",$rhotot_24;
printf "%s %.4e %s %.4e %s %.4e\n", "rhostat_34 = ",$rhostat_34," rhosyst_34 = ",$rhosyst_34," rhotot_34 = ",$rhotot_34;

printf "%s %.4e\n", "used rhotot_12 = ",($rhostat_12*$stat_1*$stat_2 + $rhosyst_12*$syst_1*$syst_2)/($sigma_1*$sigma_2);
printf "%s %.4e\n", "used rhotot_13 = ",($rhostat_13*$stat_1*$stat_3 + $rhosyst_13*$syst_1*$syst_3)/($sigma_1*$sigma_3);
printf "%s %.4e\n", "used rhotot_14 = ",($rhostat_14*$stat_1*$stat_4 + $rhosyst_14*$syst_1*$syst_4)/($sigma_1*$sigma_4);
printf "%s %.4e\n", "used rhotot_23 = ",($rhostat_23*$stat_2*$stat_3 + $rhosyst_23*$syst_2*$syst_3)/($sigma_2*$sigma_3);
printf "%s %.4e\n", "used rhotot_24 = ",($rhostat_24*$stat_2*$stat_4 + $rhosyst_24*$syst_2*$syst_4)/($sigma_2*$sigma_4);
printf "%s %.4e\n", "used rhotot_34 = ",($rhostat_34*$stat_3*$stat_4 + $rhosyst_34*$syst_3*$syst_4)/($sigma_3*$sigma_4);

# In the language of COMBOS, rhosyst_IJ*syst_I*syst_J is written as = Sum_i Delta_I,i Delta_J,i
# where i = 1..n are the n systematic variations and Delta_I,i = change in BR_I due to i^th systematic error.
# 
# The problem is to find all values of Delta_I,i, which can then be used to obtain
# the un-correlated parts of systematic errors = sqrt(syst_I **2 - Sum_i (Delta_I,i)**2).
# For consistency, the un-correlated parts should be positive numbers.
# 
# From Table III, i=1,2,6 are completely correlated between all the 4 channels.
#
# Thus we can write down all the components of the systematic errors with the notation:
#
# Delta_I_0 = sqrt(0.9**2 + 0.3**2 + 0.1**2) % = 0.953939201416946 % of the BR_I 
#           is the component which is fully correlated between channels 1, 2, 3 and 4

$common_1=0.9;
$Delta_1_0_1 = $common_1 * 1e-2 * $BR_1;
$Delta_2_0_1 = $common_1 * 1e-2 * $BR_2;
$Delta_3_0_1 = $common_1 * 1e-2 * $BR_3;
$Delta_4_0_1 = $common_1 * 1e-2 * $BR_4;

$common_2=0.3;
$Delta_1_0_2 = $common_2 * 1e-2 * $BR_1;
$Delta_2_0_2 = $common_2 * 1e-2 * $BR_2;
$Delta_3_0_2 = $common_2 * 1e-2 * $BR_3;
$Delta_4_0_2 = $common_2 * 1e-2 * $BR_4;

$common_3=0.1;
$Delta_1_0_3 = $common_3 * 1e-2 * $BR_1;
$Delta_2_0_3 = $common_3 * 1e-2 * $BR_2;
$Delta_3_0_3 = $common_3 * 1e-2 * $BR_3;
$Delta_4_0_3 = $common_3 * 1e-2 * $BR_4;

$common = sqrt(0.9**2 + 0.3**2 + 0.1**2);
$Delta_1_0 = $common * 1e-2 * $BR_1;
$Delta_2_0 = $common * 1e-2 * $BR_2;
$Delta_3_0 = $common * 1e-2 * $BR_3;
$Delta_4_0 = $common * 1e-2 * $BR_4;

#
# Delta_I_u is the un-correlated part 
#
# and the correlated parts: Delta_I_a, Delta_I_b, Delta_I_c, Delta_I_d, Delta_I_e, Delta_I_f are such that
# 
# source a induces correlation between channels 1 and 2 only
# source b induces correlation between channels 1 and 3 only
# source c induces correlation between channels 1 and 4 only
# source d induces correlation between channels 2 and 3 only
# source e induces correlation between channels 2 and 4 only
# source f induces correlation between channels 3 and 4 only
#
# We have 6 equations : 

# $rhosyst_12 * $syst_1 * $syst_2 = $Delta_1_0 * $Delta_2_0 + $Delta_1_a * $Delta_2_a
# $rhosyst_13 * $syst_1 * $syst_3 = $Delta_1_0 * $Delta_3_0 + $Delta_1_b * $Delta_3_b
# $rhosyst_14 * $syst_1 * $syst_4 = $Delta_1_0 * $Delta_4_0 + $Delta_1_c * $Delta_4_c
# $rhosyst_23 * $syst_2 * $syst_3 = $Delta_2_0 * $Delta_3_0 + $Delta_2_d * $Delta_3_d
# $rhosyst_24 * $syst_2 * $syst_4 = $Delta_2_0 * $Delta_4_0 + $Delta_2_e * $Delta_4_e
# $rhosyst_34 * $syst_3 * $syst_4 = $Delta_3_0 * $Delta_4_0 + $Delta_3_f * $Delta_4_f

# These can be parameterized as unknown fractions of the part of syst_I remaining after removing the fully correlated part

# $Delta_1_a = $a_1 * $remaining_1;
# $Delta_1_b = $b_1 * $remaining_1;
# $Delta_1_c = $c_1 * $remaining_1;
  
# $Delta_2_a = $a_2 * $remaining_2;
# $Delta_2_d = $d_2 * $remaining_2;
# $Delta_2_e = $e_2 * $remaining_2;
  
# $Delta_3_b = $b_3 * $remaining_3;
# $Delta_3_d = $d_3 * $remaining_3;
# $Delta_3_f = $f_3 * $remaining_3;
  
# $Delta_4_c = $c_4 * $remaining_4;
# $Delta_4_e = $e_4 * $remaining_4;
# $Delta_4_f = $f_4 * $remaining_4;

$remaining_1 = sqrt($syst_1**2 - $Delta_1_0**2);
$remaining_2 = sqrt($syst_2**2 - $Delta_2_0**2);
$remaining_3 = sqrt($syst_3**2 - $Delta_3_0**2);
$remaining_4 = sqrt($syst_4**2 - $Delta_4_0**2);

# Thus we have

$a_1a_2 = ($rhosyst_12 * $syst_1 * $syst_2 - $Delta_1_0 * $Delta_2_0) / ($remaining_1* $remaining_2);
$b_1b_3 = ($rhosyst_13 * $syst_1 * $syst_3 - $Delta_1_0 * $Delta_3_0) / ($remaining_1* $remaining_3);
$c_1c_4 = ($rhosyst_14 * $syst_1 * $syst_4 - $Delta_1_0 * $Delta_4_0) / ($remaining_1* $remaining_4);
$d_2d_3 = ($rhosyst_23 * $syst_2 * $syst_3 - $Delta_2_0 * $Delta_3_0) / ($remaining_2* $remaining_3);
$e_2e_4 = ($rhosyst_24 * $syst_2 * $syst_4 - $Delta_2_0 * $Delta_4_0) / ($remaining_2* $remaining_4);
$f_3f_4 = ($rhosyst_34 * $syst_3 * $syst_4 - $Delta_3_0 * $Delta_4_0) / ($remaining_3* $remaining_4);

# A set of possible solutions are:

if ($a_1a_2 > 0) {$a_1 = sqrt($a_1a_2); $a_2 = $a_1;} else {$a_1 = sqrt(-$a_1a_2); $a_2 = -$a_1;}
if ($b_1b_3 > 0) {$b_1 = sqrt($b_1b_3); $b_3 = $b_1;} else {$b_1 = sqrt(-$b_1b_3); $b_3 = -$b_1;}
if ($c_1c_4 > 0) {$c_1 = sqrt($c_1c_4); $c_4 = $c_1;} else {$c_1 = sqrt(-$c_1c_4); $c_4 = -$c_1;}
if ($d_2d_3 > 0) {$d_2 = sqrt($d_2d_3); $d_3 = $d_2;} else {$d_2 = sqrt(-$d_2d_3); $d_3 = -$d_2;}
if ($e_2e_4 > 0) {$e_2 = sqrt($e_2e_4); $e_4 = $e_2;} else {$e_2 = sqrt(-$e_2e_4); $e_4 = -$e_2;}
if ($f_3f_4 > 0) {$f_3 = sqrt($f_3f_4); $f_4 = $f_3;} else {$f_3 = sqrt(-$f_3f_4); $f_4 = -$f_3;}

# print "a_1a_2 = ",$a_1a_2," => a_1 = ",$a_1," a_2 = ",$a_2,"\n";
# print "b_1b_3 = ",$b_1b_3," => b_1 = ",$b_1," b_3 = ",$b_3,"\n";
# print "c_1c_4 = ",$c_1c_4," => c_1 = ",$c_1," c_4 = ",$c_4,"\n";
# print "d_2d_3 = ",$d_2d_3," => d_2 = ",$d_2," d_3 = ",$d_3,"\n";
# print "e_2e_4 = ",$e_2e_4," => e_2 = ",$e_2," e_4 = ",$e_4,"\n";
# print "f_3f_4 = ",$f_3f_4," => f_3 = ",$f_3," f_4 = ",$f_4,"\n";

# Then finally we have

$Delta_1_a = $a_1 * $remaining_1;
$Delta_1_b = $b_1 * $remaining_1;
$Delta_1_c = $c_1 * $remaining_1;

$Delta_2_a = $a_2 * $remaining_2;
$Delta_2_d = $d_2 * $remaining_2;
$Delta_2_e = $e_2 * $remaining_2;

$Delta_3_b = $b_3 * $remaining_3;
$Delta_3_d = $d_3 * $remaining_3;
$Delta_3_f = $f_3 * $remaining_3;

$Delta_4_c = $c_4 * $remaining_4;
$Delta_4_e = $e_4 * $remaining_4;
$Delta_4_f = $f_4 * $remaining_4;

$Delta_1_u = sqrt($syst_1**2 - ($Delta_1_0**2 + $Delta_1_a**2 + $Delta_1_b**2 + $Delta_1_c**2));
$Delta_2_u = sqrt($syst_2**2 - ($Delta_2_0**2 + $Delta_2_a**2 + $Delta_2_d**2 + $Delta_2_e**2));
$Delta_3_u = sqrt($syst_3**2 - ($Delta_3_0**2 + $Delta_3_b**2 + $Delta_3_d**2 + $Delta_3_f**2));
$Delta_4_u = sqrt($syst_4**2 - ($Delta_4_0**2 + $Delta_4_c**2 + $Delta_4_e**2 + $Delta_4_f**2));

# Printout of the results:

printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_0 = ",$Delta_1_0," = ",$common," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_0_1 = ",$Delta_1_0_1," = ",$common_1," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_0_2 = ",$Delta_1_0_2," = ",$common_2," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_0_2 = ",$Delta_1_0_3," = ",$common_3," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_a = ",$Delta_1_a," = ",100.*$Delta_1_a/$BR_1," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_b = ",$Delta_1_b," = ",100.*$Delta_1_b/$BR_1," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_c = ",$Delta_1_c," = ",100.*$Delta_1_c/$BR_1," % of ",$BR_1;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_1_u = ",$Delta_1_u," = ",100.*$Delta_1_u/$BR_1," % of ",$BR_1;

printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_0 = ",$Delta_2_0," = ",$common," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_0_1 = ",$Delta_2_0_1," = ",$common_1," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_0_2 = ",$Delta_2_0_2," = ",$common_2," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_0_2 = ",$Delta_2_0_3," = ",$common_3," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_a = ",$Delta_2_a," = ",100.*$Delta_2_a/$BR_2," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_d = ",$Delta_2_d," = ",100.*$Delta_2_d/$BR_2," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_e = ",$Delta_2_e," = ",100.*$Delta_2_e/$BR_2," % of ",$BR_2;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_2_u = ",$Delta_2_u," = ",100.*$Delta_2_u/$BR_2," % of ",$BR_2;

printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_0 = ",$Delta_3_0," = ",$common," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_0_1 = ",$Delta_3_0_1," = ",$common_1," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_0_2 = ",$Delta_3_0_2," = ",$common_2," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_0_2 = ",$Delta_3_0_3," = ",$common_3," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_b = ",$Delta_3_b," = ",100.*$Delta_3_b/$BR_3," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_d = ",$Delta_3_d," = ",100.*$Delta_3_d/$BR_3," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_f = ",$Delta_3_f," = ",100.*$Delta_3_f/$BR_3," % of ",$BR_3;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_3_u = ",$Delta_3_u," = ",100.*$Delta_3_u/$BR_3," % of ",$BR_3;

printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_0 = ",$Delta_4_0," = ",$common," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_0_1 = ",$Delta_4_0_1," = ",$common_1," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_0_2 = ",$Delta_4_0_2," = ",$common_2," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_0_2 = ",$Delta_4_0_3," = ",$common_3," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_c = ",$Delta_4_c," = ",100.*$Delta_4_c/$BR_4," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_e = ",$Delta_4_e," = ",100.*$Delta_4_e/$BR_4," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_f = ",$Delta_4_f," = ",100.*$Delta_4_f/$BR_4," % of ",$BR_4;
printf "%s %.4e %s %.4e %s %.4e\n", "Delta_4_u = ",$Delta_4_u," = ",100.*$Delta_4_u/$BR_4," % of ",$BR_4;
