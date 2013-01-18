#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## automatically get the bibtex data of the PDG tau section references
##

require(stringr, quietly=TRUE)
require(XML, quietly=TRUE)
require(RCurl, quietly=TRUE)

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

getUrlRetry = function(uri, curl = getCurlHandle(), delay=2) {
  for (i in 1:2) {
    rc = try(getURLContent(uri, followlocation = TRUE), silent=TRUE)
    if (!inherits(rc, "try-error")) break;
    Sys.sleep(delay)
  }
  if (inherits(rc, "try-error")) {
    warning(paste("error in getURLContent(", uri, ")"))
  }
  return(rc)
}

alucomb2.getpdgrefs = function(node) {
  ## Sys.setlocale('LC_ALL','C')
  curl = getCurlHandle(.encoding = "UTF-8")
  
  page.base = "http://pdg8.lbl.gov/rpp2011v1-viewer"
  page.uri = file.path(page.base, "reference.brl?parcode=S035")
  page.tree = htmlParse(page.uri, useInternalNodes = TRUE, encoding = "UTF-8")

  ref.tag1 = xpathSApply(page.tree, "//form[@name='FORM']/table[2]/tr[position() >= 5]/td[2][following-sibling::td/table]/text()", function(x) gsub("(^[\\s\ua0]+|[\\s\ua0]+$)", "", xmlValue(x), perl=TRUE))
  ref.tag2 = xpathSApply(page.tree, "//form[@name='FORM']/table[2]/tr[position() >= 5]/td[3][following-sibling::td/table]/text()", function(x) gsub("(^[\\s\ua0]+|[\\s\ua0]+$)", "", xmlValue(x), perl=TRUE))
  
  ref.onclick = xpathSApply(page.tree, "//form[@name='FORM']/table[2]/tr[position() >= 5]/td[4]/table[1]/tr[1]/td[1]/a[1]/attribute::onclick")
  ref.onclick = str_match(ref.onclick, "(REFERENCE_info.brl\\?rkeyv=\\d+)\\&")[,2]
  
  rc = mapply(function(ref.onclick, ref.tag1, ref.tag2) {
    pdg.key = paste(ref.tag1, ref.tag2, sep=".")
    page.uri = file.path(page.base, ref.onclick)

    page.html = getUrlRetry(page.uri, curl=curl)
    if (inherits(page.html, "try-error")) {
      cat("error on", pdg.key, "loading", page.uri, "\n", file=stderr())
      return(NULL)
    }
    page.tree = htmlParse(page.html, useInternalNodes = TRUE, asText=TRUE)
    insp.uri = xpathSApply(page.tree, "//frameset/frame[@name = 'head']/attribute::src")

    page.html = getUrlRetry(insp.uri, curl=curl)
    if (inherits(page.html, "try-error")) {
      cat("error on", pdg.key, "loading", insp.uri, "\n", file=stderr())
      return(NULL)
    }
    page.tree = htmlParse(page.html, useInternalNodes = TRUE, encoding = "UTF-8", asText=TRUE)
    bibtex.uri = xpathSApply(page.tree, "//ul[@class='tight_list']//li[1]/a[contains(text(), 'BibTeX')]/attribute::href")
    if (length(bibtex.uri) != 1) {
      cat("error on", pdg.key, "1 bibtex.uri needed,", length(bibtex.uri), "found\n", file=stderr())
      return(NULL)
    }

    page.html = getUrlRetry(paste("http://inspirehep.net/", bibtex.uri, sep=""), curl=curl)
    if (inherits(page.html, "try-error")) {
      cat("error on", pdg.key, "loading", bibtex.uri, "\n", file=stderr())
      return(NULL)
    }
    page.tree = htmlParse(page.html, useInternalNodes = TRUE, encoding = "UTF-8", asText=TRUE)
    bibtex.txt = xpathSApply(page.tree, "//div[@class='pagebodystripemiddle']/pre/text()", xmlValue)
    bibtex.txt = gsub("(^[\\s\ua0]+|[\\s\ua0]+$)", "", bibtex.txt, perl=TRUE)
    bibtex.key = str_match(bibtex.txt, "\\@article\\{([^,]*),")[,2]
    bibtex.txt = sub("\\}\\s*$",
      paste("      pdgkey         = \"", pdg.key, "\",\n}", sep=""),
      bibtex.txt, perl=TRUE)
    bibtex.txt = gsub("-+\\&gt\\;", "-->", bibtex.txt, perl=TRUE)
    bibtex.txt = gsub("&gt\\;", ">", bibtex.txt, perl=TRUE)
    bibtex.txt = gsub("&lt\\;", "<", bibtex.txt, perl=TRUE)
    cat("\\pdgbibalias{", pdg.key, "}{", bibtex.key, "}\n", sep="")
    cat(bibtex.txt, "\n", sep="")
    return(data.frame(bibtex=bibtex.txt, pdkkey=pdg.key))
  }, ref.onclick, ref.tag1, ref.tag2)
  
  return(invisible(NULL))
}

## ////////////////////////////////////////////////////////////////////////////
## code

args = commandArgs(TRUE)
rc = alucomb2.getpdgrefs()
