#!/usr/local/bin/perl

# print $ARGV[0]."\n";

$sfac1 = 1.63576;
$sfac2 = 2.08347;
$sfac3 = 1.55071;
$sfac4 = 5.44599;
$sfac5 = 1.43062;

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


