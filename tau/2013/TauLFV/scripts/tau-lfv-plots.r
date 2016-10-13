#!/usr/bin/env Rscript

## /////////////////////////////////////////////////////////////////////////////
##
##	tau-lfv.plots.r
##
## - read Tau LFV data in yaml
## - produce Tau LFV limits plot
## - produce Tau LFV combinations plot
##
## /////////////////////////////////////////////////////////////////////////////

'usage: tau-lfv-plots.r
' -> doc

require(docopt, quiet=TRUE)
require(stringr, quiet=TRUE)
require(yaml, quietly=TRUE)
require(ggplot2, quiet=TRUE)
require(scales, quiet=TRUE)
require(grid, quiet=TRUE)
require(latex2exp, quiet=TRUE)
## require(tikzDevice, quiet=TRUE)
## require(gridExtra, quiet=TRUE)
## require(gtable, quiet=TRUE)

##
## functions
##

##--- save plot
save.plot = function(name, plot=last_plot(), width=dev.size()[1], height=dev.size()[2], dpi=200) {
  ##
  fname.png = paste(name, "png", sep=".")
  ggsave(filename=fname.png, plot=plot, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", fname.png, "produced\n")
  ##
  fname.pdf = paste(name, "pdf", sep=".")
  ggsave(filename=fname.pdf, plot=plot, width=width, height=height)
  cat(file=stderr(), "file", fname.pdf, "produced\n")
  ##
  ##  fname.tex = paste(name, "tex", sep=".")
  ##  tikz(file = fname.tex, sanitize=TRUE)
  ##  print(plot)
  ##  dev.off()
  ##  cat(file=stderr(), "file", fname.tex, "produced\n")
}

##
## HFAG label
##
hfag.label = function(title="HFAG-Tau", subtitle="Summer 2016", fsratio=0.78, x=unit(0.9,"npc"), y = unit(0.95,"npc")) {
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

##--- convert list to data.frame
list.to.df = function(lst, keyname=NA, stringsAsFactors=FALSE) {
  rc = do.call(rbind, lapply(lst, function(x) as.data.frame(x, stringsAsFactors=stringsAsFactors)))
  if (!is.na(keyname)) {
    rc[[keyname]] = as.character(names(lst))
  }
  rc
}

##
## get linear minor breaks for log scale major breacks
##
log10.minor.breaks = function (...) {
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    major_breaks = 10^seq(minx, maxx, by=1)
    minor_sub_breaks = 2:10    
    minor_breaks =  as.vector(minor_sub_breaks %o% major_breaks)
    return(minor_breaks)
  }
}

##
## do tau LFV plot
##
tau.lfv.plot = function(data, name="plot") {
  rc = ggplot(data, aes(descr, limit)) +
    geom_point(aes(shape=factor(exp), color=factor(exp))) +
    ## labs(title="") +
    labs(y="90% CL upper limits") +
    scale_x_discrete(labels = TeX, expand=c(0,1)) +
    scale_y_continuous(
      ## limit=c(0.5e-8, 1e-6),
      trans = log10_trans(),
      breaks = trans_breaks('log10', function(x) 10^x, n=2),
      minor_breaks = log10.minor.breaks(),
      labels = trans_format('log10', math_format(10^.x))
    ) +
    annotation_logticks(sides = "l") +
    annotation_custom(hfag.label()) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 60, hjust=1, vjust=1, size=6,
        margin=unit(c(0.25, 0,    0, 0), "cm")
      ),
      axis.text.y = element_text(
        margin=unit(c(0   , 0.35, 0, 0.1), "cm")
      ),
      axis.ticks.length = unit(-0.2, "cm"),
      axis.ticks = element_line(size = 0.3),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.position = "bottom",
      legend.margin=unit(0, "null"),
      plot.margin=unit(c(0.04,0.02,0.01,0.02), "null"),
      panel.grid.major = element_line(colour = "gray80", size=0.3),
      panel.grid.minor = element_line(colour = "gray90", size=0.3)
      ## panel.grid.major.x = element_blank(),
      ## panel.grid.minor.x = element_blank()
    )
  
  print(rc)
  save.plot(name, width=my.width, height=my.height)
}

##
## main code
##

opts <- docopt(doc)

my.width = 6.5
my.height = 4
try(dev.off(), silent=TRUE)
dev.new(width=my.width, height=my.height)

##--- read Tau LFV data
tau.lfv.data = yaml.load_file("tau-lfv-data.yaml")

##
## limits plot
##
limits.df = list.to.df(tau.lfv.data$limits)
limits.df$descr = paste0("$", limits.df$descr, "$")
limits.df$descr = factor(limits.df$descr, unique(limits.df$descr))
limits.df$ref = NULL

tau.lfv.plot(limits.df, "tau-lfv-limits")

##
## combinations plot
##
combs.df = list.to.df(tau.lfv.data$combs)
combs.df$descr = paste0("$", combs.df$descr, "$")
combs.df$descr = factor(combs.df$descr, unique(combs.df$descr))
combs.df$exp = "HFAG combination"

##
## get limits data from
## - all exp limits for which a combination exists
## - the combined limit from combs.df
## - remove data not used for combinations from CLEO and ATLAS
## - refactor exp info for plotting with the remaining labels for plotting
##
data.df = rbind(
  limits.df[
    limits.df$gamma %in% combs.df$gamma &
    limits.df$exp != "CLEO" &
    limits.df$exp != "ATLAS", ],
  combs.df)
data.df$exp = factor(data.df$exp, c("BaBar", "Belle", "LHCb", "HFAG combination"))

tau.lfv.plot(data.df, "tau-lfv-combs")
