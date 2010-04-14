This directory is to cross-check the Scale Factors obtained in TauTo3Prongs/average_alucomb.log

To produce the sfactors run

root -l -b -q sfactor.cc >!  sfactor.log ; tail -4 sfactor.log

Scale Factors :
PIMPIMPIPNU = 1.61983 PIMKMPIPNU = 2.40957 PIMKMKPNU = 2.32292 KMKMKPNU = 5.43807 HMHMHPNU = 1.42995 

Note that this agrees with output of
grep "^S-factor_1" ../TauTo3Prong/log/average_alucomb.log

S-factor_1 1.6198264507 2.409577e+00 2.322927e+00 5.438124e+00
S-factor_1 1.4299445296

Enter this information in the beginning of babar.pl, belle.pl, pdg.pl, pdg_hhh.pl

$rescale = 1;
$sfac1 = 1.61983 * $rescale;
$sfac2 = 2.40957 * $rescale;
$sfac3 = 2.32292 * $rescale;
$sfac4 = 5.43807 * $rescale;
$sfac5 = 1.42995 * $rescale;

Note that in pdg.pl as well as pdg.input these measurements blocks from BEGIN to END are commented out
so as to exclude measurements with large error from the calculation of the Scale Factors

* BEGIN       OPAL PimKmPipNu published  ABBIENDI  04J
* BEGIN       CLEO PimKmPipNu published RICHICHI   99
* BEGIN       OPAL PimKmKpNu published   ABBIENDI 00D
* BEGIN       CLEO PimKmKpNu published   RICHICHI 99
* BEGIN       ALEPH PimKmKpNu published   BARATE  98

At this point one has to note that the measuerment number for the HHH results have changed.
So, now in average.input the lines have to be changed to read:

PARAMETER CHI2_N_SYM_NSUM  3 0 0
PARAMETER CHI2_N_SYM_01    13 1234 -1234
PARAMETER CHI2_N_SYM_02    14 1234 -1234
PARAMETER CHI2_N_SYM_03    15 1234 -1234

Then run |source todo|

mv average.log average_with_rescale_1.log

grep SCALE log/average_with_rescale_1.log | tail -1
 CHI2_N_SYM: CHI2/NDOF, SCALE FAC, CL =   0.86663166  0.930930534  0.572821736

Then enter this information in the beginning of babar.pl, belle.pl, pdg.pl, pdg_hhh.pl

$rescale = 0.930930534;
$sfac1 = 1.61983 * $rescale;
$sfac2 = 2.40957 * $rescale;
$sfac3 = 2.32292 * $rescale;
$sfac4 = 5.43807 * $rescale;
$sfac5 = 1.42995 * $rescale;

Note that this means

$sfac1 = 1.50795
$sfac2 = 2.24314
$sfac3 = 2.16248
$sfac4 = 5.06247
$sfac5 = 1.33118

Upto numerical precision of comparison this agrees with output of 
grep "^S-factor_2" ../TauTo3Prong/log/average_alucomb.log

S-factor_2 1.5087078077 2.244283e+00 2.163576e+00 5.065073e+00
S-factor_2 1.3318516161

Then run |source todo| again. From average.log now we have

 CHI2_N_SYM      :   M_PIMPIMPIPNU= 0.8940266E-01  +-  0.9775369E-03            CL =  0.4433
                      M_PIMKMPIPNU= 0.2918724E-02  +-  0.1578414E-03            RHO= +0.2284
                       M_PIMKMKPNU= 0.1428260E-02  +-  0.6053603E-04            RHO= +0.1944 +0.0239
                        M_KMKMKPNU= 0.2176073E-04  +-  0.7442379E-05            RHO= -0.0120 +0.0720 +0.0481


which means scale factor are

pipipi: 0.9775369E-03/ 0.7034080E-03 = 1.38972
pikpi :  0.1578414E-03 / 0.7097857E-04 = 2.22379
pikk :  0.6053603E-04 / 0.2791064E-04 = 2.16892
kkk :  0.7442379E-05 /0.1469934E-05 = 5.06307

which is quite close to final scale factors obtained in TauTo3Prongs/log/average_alucomb.log

S-factor   1.3888719483 2.197452e+00 2.147277e+00 5.064842e+00

