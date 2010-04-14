#!/usr/local/bin/perl

# print $ARGV[0]."\n";

$rescale = 1.211572;
$sfac1 = 1.61983 * $rescale;
$sfac2 = 2.40957 * $rescale;
$sfac3 = 2.32292 * $rescale;
$sfac4 = 5.43807 * $rescale;
$sfac5 = 1.42995 * $rescale;

printf "* B(tau- -> h- h- h+ nu) [ex K0]\n";
printf "BEGIN       DELPHI HmHmHpNu published ABDALLAH	 06A\n";
printf "\n";
printf "MEASUREMENT m_PimPimPipNu statistical systematic\n";
printf "DATA        m_PimPimPipNu statistical systematic\n";
printf "            9.317e-2      %.5e 0.0\n", sqrt( 0.090e-2 ** 2 + 0.082e-2 ** 2)*$sfac5;
printf "\n";
printf "END\n";
printf "\n";
printf "BEGIN       OPAL HmHmHpNu published  AKERS	 95Y\n";
printf "\n";
printf "MEASUREMENT m_PimPimPipNu statistical systematic\n";
printf "DATA        m_PimPimPipNu statistical systematic\n";
printf "	    9.87e-2       %.5e 0.0\n", sqrt( 0.10e-2 ** 2 + 0.24e-2 ** 2) * $sfac5;
printf "\n";
printf "END\n";
printf "\n";
printf "BEGIN       CLEO HmHmHpNu published  BALEST	 95C\n";
printf "\n";
printf "MEASUREMENT m_PimPimPipNu statistical systematic\n";
printf "DATA        m_PimPimPipNu statistical systematic\n";
printf "            9.51e-2       %.5e 0.0\n", sqrt( 0.07e-2 ** 2 + 0.20e-2 ** 2) * $sfac5;
printf "\n";
printf "END\n";


