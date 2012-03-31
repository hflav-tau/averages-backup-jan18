## ////////////////////////////////////////////////////////////////////////////
##
## utility functions for HFAG-tau
##

##
## return numeric id for sorting labels like "Gamma5", "Gamma3by5"
## <n>by<m> are sorted after <n> in ascending order ov <m>
##
alucomb2.gamma.num.id = function(gamma.name) {
  gamma.name = sub("Unitarity", "Gamma1000", gamma.name, fixed=TRUE)
  gamma.name = sub("GammaAll", "Gamma999", gamma.name, fixed=TRUE)
  num1 = as.numeric(str_match(gamma.name, ("^(\\D+)(\\d+[.]?\\d*)"))[,3])
  tmp = str_match(gamma.name, ("^\\D+\\d+(by|)(\\d*)"))[,3]
  num2 = ifelse(tmp=="", 0, as.numeric(tmp))
  return(1000*num1 + num2)
}
