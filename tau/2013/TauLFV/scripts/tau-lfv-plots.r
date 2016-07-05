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

##--- format ggplot2 axis labels in scientific format
fmt.num.scientific = function(val) {
  ## turn in to character string in scientific notation
  val <- format(val, scientific = TRUE)
  ## quote the part before the exponent to keep all the digits
  val <- gsub("^(.*)e", "'\\1'e", val)
  ## turn the 'e+' into plotmath format
  val <- gsub("e", "%*%10^", val)
  val <- gsub("'1'%[*]%", "", val)
  ## return this as an expression
  parse(text=val)
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

##
## main code
##

opts <- docopt(doc)

my.width = 6.5
my.height = 4
try(dev.off(), silent=TRUE)
dev.new(width=my.width, height=my.height)

num.gamma = 50
num.results = 150
num.exp = 5
labels.exp = c("CLEO", "BaBar", "Belle", "LHCb", "ATLAS")

gamma = paste0("G", ceiling(runif(num.results, 0, num.gamma)))
gamma[1] = "$\\mu\\gamma$"

data.df = data.frame(
  gamma = gamma,
  limit = 10^rnorm(num.results, -8, 1),
  exp = labels.exp[as.integer(ceiling(runif(num.results, 0, num.exp)))]
)

rc = ggplot(data.df, aes(gamma, limit)) +
  geom_point(aes(shape=factor(exp), color=factor(exp))) +
  ## labs(title="") +
  labs(y="90% CL upper limits") +
  scale_x_discrete(labels = TeX) +
  scale_y_continuous(
    ## limit=c(0, 4),
    labels=fmt.num.scientific,
    trans = log10_trans()
  ) +
  annotation_custom(hfag.label()) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1, vjust=1, size=7),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = "bottom",
    legend.margin=unit(0, "null"),
    plot.margin=unit(c(0.04,0.02,0.01,0.02), "null")
  )

print(rc)
save.plot("plot", width=my.width, height=my.height)
