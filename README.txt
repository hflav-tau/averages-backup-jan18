# To compile
cd combos ; make ; cd -

# To run combos
cd tau/2009/TauToPimKzsNu ; make ; cd -

# To obtain PDG style averages
cd tau/2009/TauToPimKzsNu ; make pdg ; cd -

# Naming Convention adopted so that we can easily assimiliate files into latex: 

For daughters of tau decays:
First negative charge, then positive charge, then neutrals, then neutrino's; within each category : low to high mass.

Eg. tau- -> e- mu- pi- K- rho- K*-    e+ mu+ pi+ K+ rho+ K*+    gamma pi0 K0S K0L eta rho0 omega K*0    eta' phi nuebar numbar nut
    Tau  To Em Mm  Pim Km Rhom Kstarm Ep Mp  Pip Kp Rhop Kstarp Gamma Piz Kzs Kzl Eta Rhoz Omega Kstarz Etap Phi Nueb   Numbar Nu

When writting out decaymodes, no space needs to be inserted between particle names, but first letter is capitalized.

Avoid latex conflicts by replacing numbers [0-9] in names with something appropriate.

We follow a style similar to hepnames/pennames. 
http://mirror.ctan.org/macros/latex/contrib/hepnames/hepnames.pdf



