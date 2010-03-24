This directory contains cards and code to check the averages of the
tau -> hhh nu BRs.

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

checkTauToHmHmHpNu_v2.r requires a non standard R package names
"maxLik", which can be privately installed by doing:

shell> R
> install.packages("maxLik")
> q()
