#!/usr/local/bin/perl

# print $ARGV[0]."\n";

$rescale = 1;
$rescale = 0.930930534;
$sfac1 = 1.61983 * $rescale;
$sfac2 = 2.40957 * $rescale;
$sfac3 = 2.32292 * $rescale;
$sfac4 = 5.43807 * $rescale;
$sfac5 = 1.42995 * $rescale;

printf "* B(tau- -> pi- pi- pi+ nu) [ex K0]\n";
printf "BEGIN    BaBar PimPimPipNu published PRL100:011801,2008\n";
printf "\n";
printf "MEASUREMENT m_PimPimPipNu statistical systematic\n";
printf "DATA        m_PimPimPipNu statistical systematic\n";
printf "            8.83370e-02   %.5e 0.00000e+00\n", 1.26940e-03*$sfac1;
printf "\n";
printf "STAT_CORR_WITH BaBar PimKmPipNu published 0.543535\n";
printf "STAT_CORR_WITH BaBar PimKmKpNu published 0.390346\n";
printf "STAT_CORR_WITH BaBar KmKmKpNu published 0.031469\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> pi- K- pi+ nu) [ex. K0]\n";
printf "BEGIN    BaBar PimKmPipNu published PRL100:011801,2008\n";
printf "\n";
printf "MEASUREMENT m_PimKmPipNu statistical systematic\n";
printf "DATA        m_PimKmPipNu statistical systematic\n";
printf "            2.72570e-03  %.5e 0.00000e+00\n",9.41768e-05*$sfac2;
printf "\n";
printf "STAT_CORR_WITH BaBar PimPimPipNu published 0.543535\n";
printf "STAT_CORR_WITH BaBar PimKmKpNu published 0.177495\n";
printf "STAT_CORR_WITH BaBar KmKmKpNu published 0.0931907\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> pi- K- K+ nu)\n";
printf "BEGIN    BaBar PimKmKpNu published PRL100:011801,2008\n";
printf "\n";
printf "MEASUREMENT m_PimKmKpNu  statistical systematic\n";
printf "DATA        m_PimKmKpNu  statistical systematic\n";
printf "            1.34610e-03  %.5e 0.00000e+00\n",3.77613e-05*$sfac3;
printf "\n";
printf "STAT_CORR_WITH BaBar PimPimPipNu published 0.390346\n";
printf "STAT_CORR_WITH BaBar PimKmPipNu published 0.177495\n";
printf "STAT_CORR_WITH BaBar KmKmKpNu published 0.0870484\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> K- K- K+ nu)\n";
printf "BEGIN    BaBar KmKmKpNu published PRL100:011801,2008\n";
printf "\n";
printf "MEASUREMENT m_KmKmKpNu   statistical systematic\n";
printf "DATA        m_KmKmKpNu   statistical systematic\n";
printf "            1.57770e-05  %.5e 0.00000e+00\n",1.79021e-06*$sfac4;
printf "\n";
printf "STAT_CORR_WITH BaBar PimPimPipNu published 0.031469\n";
printf "STAT_CORR_WITH BaBar PimKmPipNu published 0.0931907\n";
printf "STAT_CORR_WITH BaBar PimKmKpNu published 0.0870484\n";
printf "\n";
printf "END\n";
