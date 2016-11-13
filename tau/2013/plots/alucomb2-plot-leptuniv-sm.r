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
  file.png = paste(name, "png", sep=".")
  dpi = my.graph.horiz.pixels/my.graph.width
  ggsave(filename=file.png, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", file.png, "produced\n")
  file.pdf = paste(name, "pdf", sep=".")
  ggsave(filename=file.pdf, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", file.pdf, "produced\n")
}

##
## HFAG label
##
hfag.label = function(title="HFAG-Tau", subtitle="Summer 2016", fsratio=0.78, x=unit(0.81,"npc"), y = unit(0.16,"npc")) {
  hl.title = textGrob(
    title,
    x = x,
    y = y,
    gp = gpar(fontface="bold.italic", col="white")
  )

  hl.subtitle = textGrob(
    subtitle,
    x = hl.title$x,
    y = hl.title$y - unit(1.05,"lines"),
    gp = gpar(fontface="bold.italic", col="black", cex=fsratio)
  )

  hl.title.bkg = rectGrob(
    x = hl.title$x,
    y = hl.title$y,
    width = unit(0.4, "char") + unit(1, "grobwidth", hl.title),
    height = unit(0.4, "char") + unit(1, "grobheight", hl.title),
    gp=gpar(fill="black")
  )
  
  hl.box = rectGrob(
    x = hl.title$x,
    y = hl.title$y - unit(fsratio, "grobheight", hl.subtitle),
    width = unit(0.4, "char") + unit(1, "grobwidth", hl.title),
    height = unit(0.4, "char") + unit(1, "grobheight", hl.title) + unit(fsratio, "grobheight", hl.subtitle) + unit(0.9/2, "lines"),
    gp=gpar(fill="white")
  )

  grobTree(hl.box, hl.title.bkg, hl.title, hl.subtitle)
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

my.graph.height = 6/2.54
my.graph.width= my.graph.height*1.5
my.graph.horiz.pixels = 560
my.graph.x11.scale=1

try(dev.off(), silent=TRUE)
dev.new(width=my.graph.width*my.graph.x11.scale, height=my.graph.height*my.graph.x11.scale)

slope = quant$val("Be_by_tau_tau")/1e15
slope.min = (quant$val("Be_by_tau_tau") - quant$err("Be_by_tau_tau"))/1e15
slope.max = (quant$val("Be_by_tau_tau") + quant$err("Be_by_tau_tau"))/1e15

xlim.min = 289
xlim.max = 292
ylim.min = 0.177
ylim.max = 0.179

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

rc2 = ggplot_build(rc)

wisker.size = 0.03 * my.graph.width
wisker.size.x = wisker.size * (rc2$layout$panel_ranges[[1]]$x.range[2] - rc2$layout$panel_ranges[[1]]$x.range[1]) / my.graph.width
wisker.size.y = wisker.size * (rc2$layout$panel_ranges[[1]]$y.range[2] - rc2$layout$panel_ranges[[1]]$y.range[1]) / my.graph.height

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
  labs(x=TeX("$\\tau_{\\tau}$(fs)"), y=TeX("$B^{\\prime}_{\\tau e}")) +
  annotation_custom(hfag.label())
  
print(rc)

fname = sub("[.][^.]*$", "", basename(opts$file))
fname = sub("tau-br-fit-", "", fname)
fname = paste0("plot-tau-lept-univ-sm-", fname)

save.plot(fname)
