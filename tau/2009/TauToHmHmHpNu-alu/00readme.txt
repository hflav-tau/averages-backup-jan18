///////////////////////////////////////////////////////////////////////////////

	test averaging procedures (A.Lusiani)

here we test new procedures implemented in Combos to
- average >2 statistically correlated measurements
- include in the average measurements that corrspond to sum of quantities
	
warning: this is a test directory and does not work like the standard ones

///////////////////////////////////////////////////////////////////////////////

babar.r and belle.r are two R language scripts that when executed
generate on their standard output the combos cards babar.input and
belle.input.  In order to run them the R language package must be
installed. The Makefile has been configured to only usr R it it
exists.

checkTauToHmHmHpNu_v1.r is an R script that does the averages just
beteween BaBar and Belle.

checkTauToHmHmHpNu_v2.r reads directly the Combos cards that are
combined in TauToHmHmHpNu/average.input and averages all results by
BaBar, Belle and PDG, including results for indifferentiated hhh nu
final states.

has been checkTauToHmHmHpNu_v2.r requires a non standard R package names
"maxLik", which can be privately installed by doing:

shell> R
> install.packages("maxLik")
> q()

///////////////////////////////////////////////////////////////////////////////

	alucomb.r

a Combos replacement that emulates most of the current functionality
and has some extra features has been written in R and put in
"../../../Common/bin/alucomb.r"
it can by run just like Combos with a .input file as argument

- can read Combos .input files
- uses the card SUMOFMEAS to indicate measurements that correspond to
  sum of quantities we want to average
- uses the card ERROR_CORR_WITH to indicate total (stat. + syst.)
  correlations between measurements

///////////////////////////////////////////////////////////////////////////////
