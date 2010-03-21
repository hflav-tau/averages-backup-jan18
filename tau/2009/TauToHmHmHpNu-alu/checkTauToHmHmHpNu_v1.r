#!/usr/bin/env Rscript

##
## test 1st average of tau -> hhh nu by Swagato
## - including statistical correlations between 4 results in Combos
##

library(methods)

rc = setClass("bibitem",
         representation(experiment = "character",
                        method = "character",
                        status = "character",
                        where = "character",
                        results = "list"
                        )
         )

rc = setClass("parameter",
         representation(value = "numeric",
                        error = "numeric"
                        ),
         prototype(value = numeric(0),
                   error = numeric(0)
                   )
         )

rc = setClass("result",
         representation(value = "numeric",
                        stat = "numeric",
                        syst = "numeric",
                        params = "list",
                        syst.terms = "list"
                        ),
         prototype(value = numeric(0),
                   stat = numeric(0),
                   syst = numeric(0)
                   )
         )

babar = new("bibitem",
  experiment="BaBar",
  method="count",
  status="published",
  where="PRL100:011801,08"
  )

babar@results = c(babar@results,
  list(PimPimPipNu = new("result",
         value = 8.83370e-02,
         stat = 1.26940e-03,
         syst = 0.00000e+00
         )))

babar@results = c(babar@results,
  list(PimKmPipNu = new("result",
         value = 2.72570e-03,
         stat = 9.41768e-05,
         syst = 0.00000e+00
         )))

babar@results = c(babar@results,
  list(PimKmKpNu = new("result",
         value = 1.34610e-03,
         stat = 3.77613e-05,
         syst = 0.00000e+00
         )))

babar@results = c(babar@results,
  list(KmKmKpNu = new("result",
         value = 1.57770e-05,
         stat = 1.79021e-06,
         syst = 0.00000e+00
         )))

belle = new("bibitem",
  experiment="Belle",
  method="other",
  status="published",
  where="arXiv:1001.0083"
  )

belle@results = c(belle@results,
  list(PimPimPipNu = new("result",
         value = 8.42000e-02,
         stat = 2.55196e-03,
         syst = 0.00000e+00
         )))

belle@results = c(belle@results,
  list(PimKmPipNu = new("result",
         value = 3.30000e-03,
         stat = 1.65303e-04,
         syst = 0.00000e+00
         )))

belle@results = c(belle@results,
  list(PimKmKpNu = new("result",
         value = 1.55000e-03,
         stat = 5.59017e-05,
         syst = 0.00000e+00
         )))

belle@results = c(belle@results,
  list(KmKmKpNu = new("result",
         value = 3.29000e-05,
         stat = 2.58699e-06,
         syst = 0.00000e+00
         )))

quant.names = names(babar@results)
meas.names = c(
  paste("babar", names(babar@results), sep="."),
  paste("belle", names(belle@results), sep=".")
  )

meas.err = unlist(lapply(c(babar@results,belle@results), function(x) {x@stat}))
names(meas.err) = meas.names

covariance = diag(meas.err)^2
rownames(covariance) = names(meas.err)
colnames(covariance) = names(meas.err)

corr = covariance * 0

corr["babar.PimPimPipNu", "babar.PimKmPipNu"] = 0.543535
corr["babar.PimPimPipNu", "babar.PimKmKpNu"] = 0.390346
corr["babar.PimPimPipNu", "babar.KmKmKpNu"] = 0.031469
corr["babar.PimKmPipNu", "babar.PimKmKpNu"] = 0.177495
corr["babar.PimKmPipNu", "babar.KmKmKpNu"] = 0.0931907
corr["babar.PimKmKpNu", "babar.KmKmKpNu"] = 0.0870484

corr["belle.PimPimPipNu", "belle.PimKmPipNu"] = 0.1749885
corr["belle.PimPimPipNu", "belle.PimKmKpNu"] = 0.04948276
corr["belle.PimPimPipNu", "belle.KmKmKpNu"] = -0.05346557
corr["belle.PimKmPipNu", "belle.PimKmKpNu"] = 0.08026913
corr["belle.PimKmPipNu", "belle.KmKmKpNu"] = 0.03505142
corr["belle.PimKmKpNu", "belle.KmKmKpNu"] = -0.008312885

corr = (corr + t(corr))/2
covariance = covariance + (corr * (meas.err %o% meas.err))

##
## delta is a matrix that multiplied by the measurements vector
## results for each quantity we average the linear combination of
## measurements, with a minus sign
##
## when quantities corrspond to measurements, delta has -1
## for each column/row that matches quantity with measurement
##
delta = rbind(diag(rep(-1,4)), diag(rep(-1,4)))

meas = unlist(lapply(c(babar@results,belle@results), function(x) {x@value}))
names(meas) = meas.names

invcov = solve(covariance)

covariance.quant = solve(t(delta) %*% invcov %*% delta)
rownames(covariance.quant) = quant.names
colnames(covariance.quant) = quant.names
quant.err = sqrt(diag(covariance.quant))

quant = drop(-covariance.quant %*% t(delta) %*% (invcov %*% meas))
names(quant) = quant.names

cat("averages\n")
show(quant)
cat("errors\n")
show(quant.err)
cat("covariance\n")
show(covariance.quant)
