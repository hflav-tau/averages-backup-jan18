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
require(RColorBrewer, quiet=TRUE)
## require(tikzDevice, quiet=TRUE)
## require(gridExtra, quiet=TRUE)
## require(gtable, quiet=TRUE)

##
## functions
##

##
## save last plot
##
save.plot = function(name, plot = last_plot(), width=my.graph.width, height=my.graph.height, horiz.pixels=my.graph.horiz.pixels) {
  dpi = my.graph.horiz.pixels/my.graph.width
  ## png
  file.png = paste(name, "png", sep=".")
  rc = ggsave(filename=file.png, plot=plot, width=width, height=height, dpi=dpi)
  if (FALSE) {
    file.png.conv = paste(name, "-conv.", "png", sep="")
    system2("convert", c("-trim", "-border", "8x8", "-bordercolor", "white", file.png, file.png.conv), stdout=TRUE, stderr=TRUE)
    system2("mv", c("-f", file.png.conv, file.png), stdout=TRUE, stderr=TRUE)
    }
  cat(file=stderr(), "file", file.png, "produced\n")
  ## pdf
  if (FALSE) {
    file.pdf = paste(name, "pdf", sep=".")
    rc = ggsave(filename=file.pdf, plot=plot, width=width, height=height, dpi=dpi)
    cat(file=stderr(), "file", file.pdf, "produced\n")
  }
  ## svg
  if (FALSE) {
    file.svg = paste(name, "svg", sep=".")
    rc = ggsave(filename=file.svg, plot=plot, width=width, height=height, dpi=dpi)
    cat(file=stderr(), "file", file.svg, "produced\n")
  }
  ## tikz
  if (FALSE) {
    file.tex = paste(name, "tex", sep=".")
    tikz(file.tex, standAlone = FALSE, width=width, height=height, sanitize=TRUE)
    grid.draw(plot)
    dev.off()
    cat(file=stderr(), "file", file.tex, "produced\n")
  }
}

##
## HFAG label
##
hfag.label = function(title="HFAG-Tau", subtitle="Summer 2016", fsratio=0.78, x=unit(0.92,"npc"), y = unit(0.95,"npc")) {
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
    labs(title=TeX("90% CL upper limits on $\\tau$ LFV decays")) +
    labs(y="") +
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
      plot.title = element_text(
        hjust = 0.5,
        margin = margin(0.00, 0.00, 0.00, 0.00, "cm")
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
##      axis.title.y = element_text(
##        margin = margin(0.00, 0.00, 0.00, 0.00, "cm")
##      ),
      axis.text.x = element_text(
        angle = 60, hjust=1, vjust=1, size=6,
        margin = margin(0.30, 0.00, -0.40, 0.00, "cm")
      ),
      axis.text.y = element_text(
        margin = margin(0.00, 0.30, 0.00, 0.00, "cm")
      ),
      axis.ticks.length = unit(-0.2, "cm"),
      axis.ticks = element_line(size = 0.3),
      plot.margin = margin(0.02, 0.01, 0.00, 0.01, "null"),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.position = "bottom",
      legend.spacing = unit(0, "null"),
      legend.margin = margin(0, 0, 0, 0, "null"),
      panel.grid.major = element_line(colour = "gray80", size=0.3),
      panel.grid.minor = element_line(colour = "gray90", size=0.3)
      ## panel.grid.major.x = element_blank(),
      ## panel.grid.minor.x = element_blank()
    ) +
    scale_color_brewer(palette="Set1")
  
  print(rc)
  save.plot(name, width=my.graph.width, height=my.graph.height)
}

##
## main code
##

opts <- docopt(doc)

## 
## my.graph.width = 6.5
## my.graph.height = 4
## try(dev.off(), silent=TRUE)
## dev.new(width=my.graph.width, height=my.graph.height)

##--- plot width = 14cm
my.graph.width = 18/2.54
##--- plot height
my.graph.height = my.graph.width * 800 / 1300
##--- desired horizontal pixels for a plot
my.graph.horiz.pixels = 1300

##
## rescale plot size for PC screen
## not used presently
## unfortunately the plot on X11 screen does not describe well the plot on PDF and PNG
## one issue: plot adapts to size of X11 window
##
my.graph.x11.scale=1

#
## open graphic device
##
try(dev.off(), silent=TRUE)
##--- window must be sized for largest ever plot in this code
dev.new(width=my.graph.width*my.graph.x11.scale, height=my.graph.height*my.graph.x11.scale, xpos=400, ypos=540)

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
