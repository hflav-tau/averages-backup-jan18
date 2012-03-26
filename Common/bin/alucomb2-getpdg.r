## ////////////////////////////////////////////////////////////////////////////
##
## alucomb2.getpdg(node) -- get PDG11 data for a node 
##

require(XML, quietly=TRUE)
require(stringr, quietly=TRUE)

## ////////////////////////////////////////////////////////////////////////////
## functions

get.xmlval = function(tree, xpath) {
  val = paste(xpathSApply(tree, xpath, function(x) xmlValue(x)))
  val = gsub("&nbsp;", " ", val, fixed=TRUE)
  val = gsub("[\\s\ua0]+", " ", val, perl=TRUE)
  val = gsub("(^[\\s\ua0]+|[\\s\ua0]+$)", "", val, perl=TRUE)
  val = gsub("\u2212", "-", val)
  return(val)
}

utf8ToHex = function(x) {
  sprintf("%x", utf8ToInt(x))
}

## ////////////////////////////////////////////////////////////////////////////
## code

alucomb2.getpdg = function(node) {
  ## page.uri = paste("http://pdglive.lbl.gov/popupblockdata.brl?nodein=", node, sep="")
  page.uri = paste("http://pdg8.lbl.gov/rpp2011v1-viewer/popupblockdata.brl?nodein=", node, sep="")
  ## page.tree = htmlParse(page.uri, useInternalNodes = TRUE, encoding = "UTF-8")
  page.tree = htmlTreeParse(page.uri, useInternalNodes = TRUE, encoding = "UTF-8")
  pdg.text = get.xmlval(page.tree, "//div[@id='content']")

  pdginfo = list()

  ##--- get hidden value to get the correct units for PDG fit/average
  rc = str_match(pdg.text,
    "back to[^\ub1]*\\s+([[:digit:].]+)\\s+\ub1\\s+([[:digit:].]+)[^\ub1]*10\\s+([+-]|)\\s+(\\d+)")
  if (is.na(rc[1])) return(NULL)
  hidden.val = as.numeric(rc[2])
  hidden.pow = as.numeric(paste(rc[4], rc[5], sep=""))
  hidden.val = hidden.val * 10^hidden.pow
  
  ##--- get scale of measurement
  units.val = 1
  rc = str_match(pdg.text, "VALUE\\s*\\(10([+-])(\\d+)\\s*\\)")

  if (!is.na(rc[1])) {
    units.exp = -as.numeric(rc[3])
    if (rc[2] == "-") {
      units.exp = -units.exp
    }
    units.val = 10^units.exp
  }
  rc = str_match(pdg.text, "VALUE\\s*\\(%\\s*\\)")
  if (!is.na(rc[1])) {
    units.val = 100
  }

  our.str.match = paste(
    "([[:digit:].]+)\\s+\ub1\\s+([[:digit:].]+)PM[^\ub1]+OUR %s\\s+",
    "(Error includes scale factor of ([[:digit:]]+[.]*[[:digit:]]*)[.]|)", sep="")

  for (our.keyw in c("AVERAGE", "FIT")) {
    rc = str_match(pdg.text, sprintf(our.str.match, our.keyw))
    if (is.na(rc[1])) next
    val.err = try(as.numeric(rc[c(2,3)]), silent=TRUE)
    if (val.err[1] != 0) {
      units.our = 10^round(log(hidden.val/val.err[1])/log(10))
    } else {
      units.our = 1
    }
    if (inherits(val.err, "try-error")) next
    scale = 1
    if (!(rc[4] == "")) {
      scale = as.numeric(rc[5])
    }
    label = paste("PDG", our.keyw)
    pdginfo[[label]][c("val", "stat")] = val.err * units.our
    pdginfo[[label]][["syst"]] = scale
  }

  pdg.text = gsub("Error includes scale factor of [[:digit:]]+([.][[:digit:]]*|)[.]", "", pdg.text)
  
  rc = str_match_all(pdg.text,
    paste("(([[:digit:].]+)\\s+\ub1\\s+([[:digit:].]+)",
          "(\\s+\ub1\\s+([[:digit:].]+)|document.write.*\\)\\};)",
          "(?:\\s+f&a|\\s+avg|)(?:\\s+[[:alnum:].]+|)",
          "\\s+PM[^\ub1*]+[]];([[:alpha:]]+\\s+[[:alnum:]][[:alnum:]][[:alpha:]]*)\\s+",
          "|We do not use the following data)", sep=""))[[1]]

  if (length(rc) > 0) {
    for (i in 1:nrow(rc)) {
      if (rc[i,3] == "") break
      val = suppressWarnings(try(as.numeric(rc[i,c(3,4,6)]), silent=TRUE))
      if (inherits(rc, "try-error")) {
        cat("error", rc, "\n")
      }
      label = rc[i,7]
      pdginfo[[label]][c("val", "stat", "syst")] = val / units.val
    }
  }    

  return(do.call(rbind, lapply(pdginfo, function(x) as.data.frame(t(x)))))
}
