#!/usr/bin/env Rscript

## /////////////////////////////////////////////////////////////////////////////
##
##	tau-lfv.plots.r
##
## - test code for plotting Tau LFV upper limits
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
## require(gridExtra, quiet=TRUE)
## require(gtable, quiet=TRUE)

##
## functions
##

##--- save plot
save.plot = function(name, plot=last_plot(), width=dev.size()[1], height=dev.size()[2], dpi=80) {
  file.png = paste(name, "png", sep=".")
  ggsave(filename=file.png, plot=plot, width=width, height=height, dpi=dpi)
  cat(file=stderr(), "file", file.png, "produced\n")
  file.pdf = paste(name, "pdf", sep=".")
  ggsave(filename=file.pdf, plot=plot, width=width, height=height)
  cat(file=stderr(), "file", file.pdf, "produced\n")
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

  hl.title.bkg = rectGrob(
    x = hl.title$x,
    y = hl.title$y,
    width = unit(0.4, "char") + unit(1, "grobwidth", hl.title),
    height = unit(0.4, "char") + unit(1, "grobheight", hl.title),
    gp=gpar(fill="black")
  )
  
  hl.subtitle = textGrob(
    subtitle,
    x = hl.title$x,
    y = hl.title$y - unit(1.05,"lines"),
    gp = gpar(fontface="bold.italic", col="black", cex=fsratio)
  )

  hl.box = rectGrob(
    x = hl.title$x,
    y = hl.title$y - unit(fsratio, "grobheight", hl.subtitle),
    width = unit(0.4, "char") + unit(1, "grobwidth", hl.title),
    height = unit(0.4, "char") + unit(1, "grobheight", hl.title) + unit(fsratio, "grobheight", hl.subtitle) + unit(0.9/2, "lines"),
    gp=gpar(fill="NA")
  )

  grobTree(hl.title.bkg, hl.title, hl.subtitle, hl.box)
}

##--- convert list to data.frame
list.to.df = function(lst, keyname=NA, stringsAsFactors=FALSE) {
  rc = do.call(rbind, lapply(lst, function(x) as.data.frame(x, stringsAsFactors=stringsAsFactors)))
  if (!is.na(keyname)) {
    rc[[keyname]] = as.character(names(lst))
  }
  rc
}

log10_minor_break = function (...) {
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

##
## main code
##

opts <- docopt(doc)

my.width = 6.5
my.height = 4
try(dev.off(), silent=TRUE)
dev.new(width=my.width, height=my.height)

tau.lfv.data.combs = read.csv("tau-lfv-data-combs.csv", stringsAsFactors=FALSE)
tau.lfv.data = yaml.load_file("tau-lfv-data.yaml")

##--- data.frame to plot LVF combinations limits
data.df = data.frame(
  gamma = factor(tau.lfv.data.combs$descr, unique(tau.lfv.data.combs$descr)),
  limit = tau.lfv.data.combs$limit,
  exp = "HFAG"
)

##--- data.frame to plot LVF limits
data.df = list.to.df(tau.lfv.data$limits)
data.df$descr = factor(data.df$descr, unique(data.df$descr))

rc = ggplot(data.df, aes(descr, limit)) +
  geom_point(aes(shape=factor(exp), color=factor(exp))) +
  ## labs(title="") +
  labs(y="90% CL upper limits") +
  scale_x_discrete(labels = TeX, expand=c(0,1)) +
  scale_y_continuous(
    ## limit=c(0.5e-8, 1e-6),
    trans = log10_trans(),
    breaks = trans_breaks('log10', function(x) 10^x, n=2),
    ## minor_breaks = trans_breaks("log10", function(x) 10^x, n=20),
    ## minor_breaks = print(c(sapply(breaks, function(x) seq(0, x, x/10)))),
    minor_breaks = as.vector(2:10 %o% 10^((-9):(-5))),
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
save.plot("plot", width=my.width, height=my.height)
