#!/usr/bin/perl
$x1     = $ARGV[0];
$sigma1 = sqrt($ARGV[1]);
$x2     = $ARGV[2];
$sigma2 = sqrt($ARGV[3]);
$rho    = $ARGV[4]/($sigma1*$sigma2);
print "x1 = ", $x1, " +- ", $sigma1, " x2 = ", $x2, " +- " , $sigma2, " rho = ", $rho, "\n";
$den    = $sigma1*$sigma1 + $sigma2*$sigma2 - 2 * $rho * $sigma1 * $sigma2;
$w1     = ($sigma2*$sigma2 - $rho * $sigma1 * $sigma2)/$den;
$w2     = ($sigma1*$sigma1 - $rho * $sigma1 * $sigma2)/$den;
print "w1 = ", $w1, " w2 = ", $w2, " w1+w2 = ", $w1 + $w2, "\n";
$aver   = $w1 * $x1 + $w2 * $x2;
$sigma  = sqrt((( 1 - $rho * $rho) * $sigma1* $sigma1 * $sigma2 * $sigma2)/$den);
print "<x> = ",$aver, " +- ", $sigma, "\n";
