#!/usr/bin/env Rscript

require(methods, quietly=TRUE)
require(docopt, quiet=TRUE)
require(stringr, quiet=TRUE)
require(grid, quiet=TRUE)
require(ggplot2, quiet=TRUE)
require(scales, quiet=TRUE)
require(Matrix, quietly=TRUE)
require(latex2exp, quiet=TRUE)

source("../../../Common/bin/aluelab3.r")
source("../../../Common/bin/alureport.r")

"Usage:
  alucomb2-plot-leptuniv-sm.r [options] <file>

Arguments:
  <file>        .rdata file produced by aluelab3-results.r

Options:
" -> doc

##
## functions
##

load.in.list <- function(.file.name) { load(.file.name); as.list(environment()) }

##--- save plot
save.plot = function(name, width=my.graph.width, height=my.graph.height, horiz.pixels=my.graph.horiz.pixels) {
  dpi = my.graph.horiz.pixels/my.graph.width
  file.png = paste(name, "png", sep=".")
  ggsave(filename=file.png, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", file.png, "produced\n")
  file.pdf = paste(name, "pdf", sep=".")
  ggsave(filename=file.pdf, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", file.pdf, "produced\n")
}

##
## HFAG label
##
hfag.label = function(title="HFAG-Tau", subtitle="Summer 2016", scale=1.2, fsratio=0.78, x=unit(0.845,"npc"), y = unit(0.125,"npc")) {
  ##
  ## draw title at desired position, font size scaled
  ##
  hl.title = textGrob(
    title,
    x = x,
    y = y,
    gp = gpar(fontface="bold.italic", col="white", cex=scale)
  )
  ##
  ## draw subtitle below title, with smaller font than title by fsratio
  ##
  hl.subtitle = textGrob(
    subtitle,
    x = hl.title$x,
    y = hl.title$y - unit(1.3/(fsratio*scale), "grobheight", hl.title),
    gp = gpar(fontface="bold.italic", col="black", cex=fsratio*scale)
  )
  ##
  ## draw black title background
  ##
  hl.title.bkg = rectGrob(
    x = hl.title$x,
    y = hl.title$y,
    width = unit(1/scale, "grobwidth", hl.title) + unit(0.02, "npc"),
    height = unit(1, "grobheight", hl.title) + unit(0.02, "npc"),
    gp=gpar(fill="black", col="black", lwd=1, cex=scale)
  )
  ##--- draw subtitle background
  hl.box = rectGrob(
    x = hl.subtitle$x,
    y = hl.subtitle$y,
    width = unit(1/(fsratio*scale), "grobwidth", hl.title.bkg) - unit(0.003, "npc"),
    ## height = unit(1, "grobheight", hl.subtitle) + 0.4*unit(1/(fsratio), "grobheight", hl.title),
    height = unit(1, "grobheight", hl.subtitle) + unit(0.02, "npc"),
    gp=gpar(fill="white", col="black", lwd=1, cex=fsratio*scale)
  )
  grobTree(hl.title.bkg, hl.title, hl.box, hl.subtitle)
}

##
## main code
##

##--- get options and args
opts <- docopt(doc, quoted_args = TRUE)
## opts = list(file="hfag-2016.rdata")
## opts = list(file="tau-br-fit-hfag-2016-elab-vadirect.rdata")

##--- load alucomb2 fit data
load(opts$file)

rc = quant$quant.expr.add(
  "Be_by_tau_tau",
  1/tau_mu * m_tau^5/m_mu^5 * (phspf_mebymmu * delta_mu_W * delta_mu_gamma)/ (phspf_mebymmu * delta_mu_W * delta_mu_gamma))

##
## plot
##

##--- plot height 6cm
my.graph.height = 8/2.54
##--- plot width = 3/2 plot height
my.graph.width = my.graph.height*1.5
##--- set desired horizontal pixels for a plot
my.graph.horiz.pixels = 560

##
## rescale plot size for PC screen
## not used presently
## unfortunately the plot on X11 screen does not describe well the plot on PDF and PNG
##
my.graph.x11.scale=1

try(dev.off(), silent=TRUE)
dev.new(width=my.graph.width*my.graph.x11.scale, height=my.graph.height*my.graph.x11.scale)

slope = quant$val("Be_by_tau_tau")/1e15
slope.min = (quant$val("Be_by_tau_tau") - quant$err("Be_by_tau_tau"))/1e15
slope.max = (quant$val("Be_by_tau_tau") + quant$err("Be_by_tau_tau"))/1e15

##--- plot limits
xlim.min = 289
xlim.max = 292
ylim.min = 0.177
ylim.max = 0.179

##
## compute band of SM prediction
## - take xmin on the left of plot left limit and xmax on the right of plot right limit
## - for xmin and xmax compute the min and max of the SM prediction for the leptonic BR from the lifetime
##
band.xmin = xlim.min - 0.2*(xlim.max-xlim.min)
band.xmax = xlim.max + 0.2*(xlim.max-xlim.min)
band.xmin.ymin = band.xmin*slope.min
band.xmin.ymax = band.xmin*slope.max
band.xmax.ymin = band.xmax*slope.min
band.xmax.ymax = band.xmax*slope.max

plot.df = data.frame(
  tau = quant$val("tau_tau")*1e15,
  Bl = quant$val("Be_lept"),
  tau.min = (quant$val("tau_tau") - quant$err("tau_tau")) *1e15,
  tau.max = (quant$val("tau_tau") + quant$err("tau_tau")) *1e15,
  Bl.min = quant$val("Be_lept") - quant$err("Be_lept"),
  Bl.max = quant$val("Be_lept") + quant$err("Be_lept"))

rc = ggplot(data = plot.df, aes(x = tau, y = Bl)) +
  theme_bw() +
  coord_cartesian(xlim = c(xlim.min, xlim.max), ylim = c(ylim.min, ylim.max))

##
## compute proper size in the two axis units in order to get the
## errorbar wiskers of the same size on x and y
##
rc2 = ggplot_build(rc)
wisker.size = 0.03 * my.graph.width
wisker.size.x = wisker.size * (rc2$layout$panel_ranges[[1]]$x.range[2] - rc2$layout$panel_ranges[[1]]$x.range[1]) / my.graph.width
wisker.size.y = wisker.size * (rc2$layout$panel_ranges[[1]]$y.range[2] - rc2$layout$panel_ranges[[1]]$y.range[1]) / my.graph.height
rm(rc2)

rc = rc +
  geom_errorbar(aes(ymin=Bl.min,  ymax=Bl.max),  color="Red",
                width=wisker.size.x,
                size=0.6666) +
  geom_errorbarh(aes(xmin=tau.min, xmax=tau.max), color="Red",
                 height=wisker.size.y,
                 size=0.6666) +
  geom_point(color="Red") +
  geom_polygon(
    data =
      data.frame(
        x=c(band.xmin, band.xmin, band.xmax, band.xmax),
        y=c(band.xmin.ymax, band.xmin.ymin, band.xmax.ymin, band.xmax.ymax)),
    aes(x=x, y=y), fill="Yellow") +
  ## labs(title="Lepton Universality Test") +
  labs(x=TeX("$\\tau_{\\tau}$(fs)"), y=TeX("$\\textit{B}_{\\tau e}^{\\prime}$")) +
  annotation_custom(hfag.label())

print(rc)

fname = sub("[.][^.]*$", "", basename(opts$file))
fname = sub("tau-br-fit-", "", fname)
fname = paste0("plot-tau-lept-univ-sm-", fname)

save.plot(fname)
