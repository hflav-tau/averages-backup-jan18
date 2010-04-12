///////////////////////////////////////////////////////////////////////////////

	information on tau/2009/TauToHmHmHpNu-alu2

///////////////////////////////////////////////////////////////////////////////

This directory contains input data files to check the procedure that
has been used by Swagato in tau/2009/TauToHmHmHpNu to do a combined
average of tau -> hhh nu BRs (with h = pi, K).

In the original directory, the published total correlation coefficients
between different channels in the same experiment have been split in
their statistical contribution and their systematic contribution and the
latter ones have been represented by systematic contributions dependant
from common fictitious external parameters. By so doing the
statistical parte of the correlation could be commuinicated to Combos
with the STAT_CORR_WITH cards.

In this directory the total correlation is represented as statistical
correlation and fully communicated to Combos with the STAT_CORR_WITH
cards. Since Combos interprets this latter card as correlation only
regarding the statistical error in the input cards, the statistical
error has been seto to the total error and the systematic error has
been set to zero.  The whole procedure is possible because the BaBar
and Belle resulta do not have any common systematic between them and
all systematic correlations are included in the published total
correlation coefficients.

The directory includes two executable scripts in the R language.
- babar.r outputs babar.input
- belle.r outputs belle.input
The script are left as reference, there is no need to run them because
the .input files are also stored in the subversion repository.
