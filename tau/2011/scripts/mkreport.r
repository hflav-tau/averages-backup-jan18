#!/usr/bin/env Rscript

##
## mkreport.r
##

require("optparse", quietly=TRUE)
require("stringr", quietly=TRUE)
source("../../../Common/bin/alucomb2-hfag-tau.r")
source("../../../Common/bin/alucomb2-fit.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## example function to use vars of an env in the function env
##
mkreport.rc = function(alucomb.env) {
  attach(alucomb.env)
  rc = try({
  })
  detach(alucomb.env)
  return(invisible(NULL))
}

##
## get the inspire bibtex key from the measurement tags
##
get.reference = function(tags) {
  bib.pdg.table = list(
    "DEL-AMO-SANCHEZ.11E" = "delAmoSanchez:2010pc",
    "AUBERT.10B" = "Aubert:2009ag",
    "AUBERT.10F" = "Aubert:2009qj",
    "HAYASAKA.10" = "Hayasaka:2010np",
    "LEE.10" = "Lee:2010tc",
    "LEES.10A" = "not found",
    "MIYAZAKI.10" = "Miyazaki:2009wc",
    "MIYAZAKI.10A" = "Miyazaki:2010qb",
    "AUBERT.09AK" = "Aubert:2009ra",
    "AUBERT.09D" = "Aubert:2009ys",
    "AUBERT.09W" = "Aubert:2009ap",
    "GROZIN.09A" = "Grozin:2008nw",
    "INAMI.09" = "Inami:2008ar",
    "MIYAZAKI.09" = "Miyazaki:2008mw",
    "AUBERT.08" = "Aubert:2007mh",
    "AUBERT.08AE" = "Aubert:2008nj",
    "AUBERT.08K" = "Aubert:2007kx",
    "FUJIKAWA.08" = "Fujikawa:2008ma",
    "HAYASAKA.08" = "Hayasaka:2007vc",
    "MIYAZAKI.08" = "Miyazaki:2007zw",
    "NISHIO.08" = "Nishio:2008zx",
    "ANASHIN.07" = "Anashin:2007zz",
    "AUBERT.07AP" = "Aubert:2007jh",
    "AUBERT.07BK" = "Aubert:2007pw",
    "AUBERT.07I" = "Aubert:2006cz",
    "BELOUS.07" = "Abe:2006vf",
    "EPIFANOV.07" = "Epifanov:2007rf",
    "MIYAZAKI.07" = "Miyazaki:2007jp",
    "ABDALLAH.06A" = "Abdallah:2003cw",
    "AUBERT.06C" = "Aubert:2005wa",
    "AUBERT,B.06" = "Aubert:2006hw",
    "INAMI.06" = "Inami:2006vd",
    "MIYAZAKI.06" = "Miyazaki:2005ng",
    "MIYAZAKI.06A" = "Miyazaki:2006sx",
    "YUSA.06" = "Yusa:2006qq",
    "ARMS.05" = "Arms:2005qg",
    "AUBERT,B.05A" = "Aubert:2005ye",
    "AUBERT,B.05F" = "Aubert:2005nx",
    "AUBERT,B.05W" = "Aubert:2005waa",
    "AUBERT,BE.05D" = "Aubert:2005tp",
    "ENARI.05" = "Enari:2005gc",
    "HAYASAKA.05" = "Hayasaka:2005xw",
    "SCHAEL.05C" = "Schael:2005am",
    "ABBIENDI.04J" = "Abbiendi:2004xa",
    "ABDALLAH.04K" = "Abdallah:2003xd",
    "ABDALLAH.04T" = "Abdallah:2003yq",
    "ABE.04B" = "Abe:2003sx",
    "ACHARD.04G" = "Achard:2004jj",
    "AUBERT.04J" = "Aubert:2003pc",
    "ENARI.04" = "Enari:2004ax",
    "YUSA.04" = "Yusa:2004gm",
    "ABBIENDI.03" = "Abbiendi:2002jw",
    "BRIERE.03" = "Briere:2003fr",
    "HEISTER.03F" = "Heister:2002ik",
    "INAMI.03" = "Inami:2002ah",
    "CHEN.02C" = "not found",
    "ABBIENDI.01J" = "Abbiendi:2000ee",
    "ABREU.01M" = "Abreu:2001wc",
    "ACCIARRI.01F" = "Acciarri:2001sg",
    "ACHARD.01D" = "Achard:2001pk",
    "ANASTASSOV.01" = "Anastassov:2000xu",
    "HEISTER.01E" = "Heister:2001me",
    "ABBIENDI.00A" = "Abbiendi:2000kz",
    "ABBIENDI.00C" = "Abbiendi:1999pm",
    "ABBIENDI.00D" = "Abbiendi:1999cq",
    "ABREU.00L" = "Abreu:2000sg",
    "ACCIARRI.00B" = "Acciarri:2000vq",
    "AHMED.00" = "not found",
    "ALBRECHT.00" = "Albrecht:2000yg",
    "ASNER.00" = "Asner:1999kj",
    "ASNER.00B" = "Asner:2000nx",
    "BERGFELD.00" = "Bergfeld:1999yh",
    "BROWDER.00" = "Browder:1999fr",
    "EDWARDS.00A" = "Edwards:1999fj",
    "GONZALEZ-SPRINBERG.00" = "GonzalezSprinberg:2000mk",
    "ABBIENDI.99H" = "Abbiendi:1998cx",
    "ABREU.99X" = "Abreu:1999rb",
    "ACKERSTAFF.99D" = "Ackerstaff:1998yk",
    "ACKERSTAFF.99E" = "Ackerstaff:1998ia",
    "BARATE.99K" = "Barate:1999hi",
    "BARATE.99R" = "Barate:1999hj",
    "BISHAI.99" = "Bishai:1998gf",
    "GODANG.99" = "Godang:1999ge",
    "RICHICHI.99" = "Richichi:1998bc",
    "ACCIARRI.98C" = "Acciarri:1998zc",
    "ACCIARRI.98E" = "Acciarri:1998iv",
    "ACCIARRI.98R" = "Acciarri:1998as",
    "ACKERSTAFF.98M" = "Ackerstaff:1997tx",
    "ACKERSTAFF.98N" = "Ackerstaff:1998mt",
    "ALBRECHT.98" = "Albrecht:1997gn",
    "BARATE.98" = "Barate:1997ma",
    "BARATE.98E" = "Barate:1997tt",
    "BLISS.98" = "Bliss:1997iq",
    "ABE.97O" = "Abe:1997dy",
    "ACKERSTAFF.97J" = "Ackerstaff:1997is",
    "ACKERSTAFF.97L" = "Ackerstaff:1996gy",
    "ACKERSTAFF.97R" = "Ackerstaff:1997dv",
    "ALEXANDER.97F" = "Alexander:1997bv",
    "AMMAR.97B" = "Ammar:1996xh",
    "ANASTASSOV.97" = "Anastassov:1996tc",
    "ANDERSON.97" = "Anderson:1997kb",
    "AVERY.97" = "not found",
    "BARATE.97I" = "Barate:1997hw",
    "BARATE.97R" = "Barate:1997be",
    "BERGFELD.97" = "Bergfeld:1997zt",
    "BONVICINI.97" = "Bonvicini:1997bw",
    "BUSKULIC.97C" = "Buskulic:1996qs",
    "COAN.97" = "Coan:1997am",
    "EDWARDS.97" = "not found",
    "EDWARDS.97B" = "not found",
    "ESCRIBANO.97" = "Escribano:1996wp",
    "ABREU.96B" = "Abreu:1995xt",
    "ACCIARRI.96H" = "Acciarri:1996ht",
    "ACCIARRI.96K" = "Acciarri:1996kc",
    "ALAM.96" = "Alam:1995mt",
    "ALBRECHT.96E" = "Albrecht:1996gr",
    "ALEXANDER.96D" = "Alexander:1995tw",
    "ALEXANDER.96E" = "Alexander:1996nc",
    "ALEXANDER.96S" = "Alexander:1996um",
    "BAI.96" = "Bai:1995hf",
    "BALEST.96" = "Balest:1996cr",
    "BARTELT.96" = "Bartelt:1996iv",
    "BUSKULIC.96" = "Buskulic:1995ty",
    "BUSKULIC.96C" = "Buskulic:1995rh",
    "COAN.96" = "Coan:1996iu",
    "ABE.95Y" = "Abe:1995yj",
    "ABREU.95T" = "Abreu:1995gr",
    "ABREU.95U" = "Abreu:1995gs",
    "ACCIARRI.95" = "Acciarri:1994vr",
    "ACCIARRI.95F" = "Acciarri:1995kx",
    "AKERS.95F" = "Akers:1994xf",
    "AKERS.95I" = "Akers:1995nh",
    "AKERS.95P" = "Akers:1995vy",
    "AKERS.95Y" = "Akers:1995ry",
    "ALBRECHT.95" = "Albrecht:1994nm",
    "ALBRECHT.95C" = "Albrecht:1995un",
    "ALBRECHT.95G" = "Albrecht:1995ht",
    "ALBRECHT.95H" = "Albrecht:1995wx",
    "BALEST.95C" = "Balest:1995kq",
    "BUSKULIC.95C" = "Buskulic:1994yh",
    "BUSKULIC.95D" = "Buskulic:1994hi",
    "ABREU.94K" = "Abreu:1994fi",
    "AKERS.94E" = "Akers:1994ex",
    "AKERS.94G" = "Akers:1994td",
    "ALBRECHT.94E" = "Albrecht:1994es",
    "ARTUSO.94" = "Artuso:1994ii",
    "BARTELT.94" = "Bartelt:1994hn",
    "BATTLE.94" = "Battle:1994by",
    "BAUER.94" = "Bauer:1993wn",
    "BUSKULIC.94D" = "Buskulic:1993mn",
    "BUSKULIC.94E" = "Buskulic:1994je",
    "BUSKULIC.94F" = "Buskulic:1994jf",
    "GIBAUT.94B" = "Gibaut:1994ik",
    "ADRIANI.93M" = "Adriani:1993gk",
    "ALBRECHT.93C" = "Albrecht:1992ka",
    "ALBRECHT.93G" = "Albrecht:1993fr",
    "BALEST.93" = "not found",
    "BEAN.93" = "Bean:1992hh",
    "BORTOLETTO.93" = "Bortoletto:1993px",
    "ESCRIBANO.93" = "Escribano:1993pq",
    "PROCARIO.93" = "Procario:1992hd",
    "ABREU.92N" = "Abreu:1992gn",
    "ACTON.92F" = "Acton:1992ff",
    "ACTON.92H" = "Acton:1992by",
    "AKERIB.92" = "Akerib:1992hb",
    "ALBRECHT.92D" = "Albrecht:1991rh",
    "ALBRECHT.92K" = "Albrecht:1992uba",
    "ALBRECHT.92M" = "Albrecht:1992td",
    "ALBRECHT.92Q" = "Albrecht:1992pe",
    "AMMAR.92" = "Ammar:1991jc",
    "ARTUSO.92" = "Artuso:1992qu",
    "BAI.92" = "Bai:1992bu",
    "BATTLE.92" = "Battle:1992qv",
    "BUSKULIC.92J" = "Buskulic:1992jj",
    "DECAMP.92C" = "Decamp:1991jp",
    "ADEVA.91F" = "Adeva:1991qq",
    "ALBRECHT.91D" = "Albrecht:1990nc",
    "ALEXANDER.91D" = "Alexander:1991am",
    "ANTREASYAN.91" = "Antreasyan:1990fq",
    "GRIFOLS.91" = "Grifols:1990ha",
    "ABACHI.90" = "Abachi:1989vr",
    "ALBRECHT.90E" = "Albrecht:1990zj",
    "BEHREND.90" = "Behrend:1989qx",
    "BOWCOCK.90" = "Bowcock:1989mq",
    "DELAGUILA.90" = "delAguila:1990jg",
    "GOLDBERG.90" = "Goldberg:1990ep",
    "WU.90" = "Wu:1989hx",
    "ABACHI.89B" = "Abachi:1988gx",
    "BEHREND.89B" = "Behrend:1989wc",
    "JANSSEN.89" = "Janssen:1989wg",
    "KLEINWORT.89" = "Kleinwort:1988sb",
    "ADEVA.88" = "Adeva:1988eq",
    "ALBRECHT.88B" = "Albrecht:1987zf",
    "ALBRECHT.88L" = "Albrecht:1988bt",
    "ALBRECHT.88M" = "Albrecht:1988ik",
    "AMIDEI.88" = "Amidei:1987zj",
    "BEHREND.88" = "Behrend:1987fa",
    "BRAUNSCHWEIG.88C" = "Braunschweig:1988yy",
    "KEH.88" = "Keh:1988gs",
    "TSCHIRHART.88" = "Tschirhart:1987uh",
    "ABACHI.87B" = "Abachi:1987hw",
    "ABACHI.87C" = "Abachi:1987nh",
    "ADLER.87B" = "Adler:1987bf",
    "AIHARA.87B" = "Aihara:1986mw",
    "AIHARA.87C" = "Aihara:1987zm",
    "ALBRECHT.87L" = "not found",
    "ALBRECHT.87P" = "Albrecht:1987fb",
    "BAND.87" = "Band:1987nv",
    "BAND.87B" = "Band:1987bm",
    "BARINGER.87" = "Baringer:1987tr",
    "BEBEK.87C" = "Bebek:1987nh",
    "BURCHAT.87" = "Burchat:1986yp",
    "BYLSMA.87" = "Bylsma:1986zy",
    "COFFMAN.87" = "Coffman:1987da",
    "DERRICK.87" = "Derrick:1987sp",
    "FORD.87" = "Ford:1986zk",
    "FORD.87B" = "Ford:1987ha",
    "GAN.87" = "Gan:1987mr",
    "GAN.87B" = "Gan:1987xs",
    "AIHARA.86E" = "Aihara:1986nj",
    "BARTEL.86D" = "Bartel:1986un",
    "RUCKSTUHL.86" = "Ruckstuhl:1986qg",
    "SCHMIDKE.86" = "Schmidke:1986gp",
    "YELTON.86" = "Yelton:1985sx",
    "ALTHOFF.85" = "Althoff:1984ja",
    "ASH.85B" = "Ash:1985hp",
    "BALTRUSAITIS.85" = "Baltrusaitis:1985fh",
    "BARTEL.85F" = "not found",
    "BEHRENDS.85" = "Behrends:1985pm",
    "BELTRAMI.85" = "Beltrami:1985vb",
    "BERGER.85" = "Berger:1985in",
    "BURCHAT.85" = "Burchat:1985mv",
    "FERNANDEZ.85" = "Fernandez:1984if",
    "MILLS.85" = "Mills:1985mh",
    "AIHARA.84C" = "Aihara:1984hr",
    "BEHREND.84" = "Behrend:1984up",
    "MILLS.84" = "Mills:1984dn",
    "BEHREND.83C" = "not found",
    "SILVERMAN.83" = "Silverman:1982ft",
    "BEHREND.82" = "not found",
    "BLOCKER.82B" = "Blocker:1982rw",
    "BLOCKER.82D" = "not found",
    "FELDMAN.82" = "Feldman:1981md",
    "HAYES.82" = "Hayes:1981bn",
    "BERGER.81B" = "not found",
    "DORFAN.81" = "Dorfan:1980mk",
    "BRANDELIK.80" = "not found",
    "BACINO.79B" = "Bacino:1979fz",
    "BACINO.78B" = "Bacino:1978gb",
    "BRANDELIK.78" = "not found",
    "JAROS.78" = "Jaros:1978um",
    "DAVIER.06" = "Davier:2005xq",
    "RAHAL-CALLOT.98" = "RahalCallot:1998ki",
    "GENTILE.96" = "Gentile:1995ue",
    "WEINSTEIN.93" = "not found",
    "PERL.92" = "Perl:1991gd",
    "PICH.90" = "Pich:1990gh",
    "BARISH.88" = "Barish:1987nj",
    "GAN.88" = "Gan:1987zi",
    "HAYES.88" = "Hayes:1988ut",
    "PERL.80" = "not found"
    )

  bib.conf.table = list(
    "BaBar.ICHEP08" = "Aubert:2008an",
    "BaBar.DPF09" = "Paramesvaran:2009ec",
    "Belle.PHIPSI11" = "SooRyu:2011aa"
    )

  ref.conf = bib.conf.table[[paste(tags[1], tail(tags, -3), collapse=".", sep=".")]]
  if (!is.null(ref.conf)) return(ref.conf)

  ref.insp = bib.pdg.table[[paste(tail(tags, -3), collapse=".")]]
  if (is.null(ref.insp) || ref.insp=="not found") {
    ref.insp = paste("not found:", paste(tags, collapse="."))
  }
  return(ref.insp)
}

##
##
##
get.tex.quant.descr = function(quant) {
  repeat {
    texdescr = quant$texdescr
    if (!is.null(texdescr) && texdescr != "") {
      ##+++ check for ratio of BR
      texdescr = paste("\\BRF{\\tau^-}{", texdescr, "}", sep="")
      texdescr = gsub("mathrm", "text", texdescr, fixed=TRUE)
      break
    }
    descr = quant$descr
    if (is.null(descr) || descr == "") {
      texdescr = ""
      break
    }
    descr = sub("\\s+/\\s+G\\(total\\)\\s*", "", descr)
    descr = gsub("ex[.]\\s+K\\(S\\)0 --> pi- pi[+]", "ex. K0", descr, perl=TRUE)
    descr = gsub("\\s*\\(``\\d-prong''\\)", "", descr, perl=TRUE)
    descr = gsub("-->", "\\to", descr, fixed=TRUE)
    descr = gsub("([+-])", "^\\1", descr)
    descr = gsub("\\^-prong", "\\\\text{-prong}", descr)
    descr = gsub("G(\\(((?:[^()]++|(?1))+)*+\\))", "\\\\BRF{tau^-}{\\2}", descr, perl=TRUE)
    descr = gsub(">=", "\\\\ge{}", descr)
    descr = gsub("nubar\\((mu|tau|e)\\)", "\\\\bar{nu}_\\1", descr)
    descr = gsub("nu\\((mu|tau|e)\\)", "nu_\\1", descr)
    descr = gsub("pi0", "pi^0", descr)
    descr = gsub("Kstar", "K^*", descr)
    descr = gsub("Kbar0", "\\\\bar{K}^0", descr)
    descr = gsub("K0", "K^0", descr)
    descr = gsub("K\\(S\\)0", "K_S^0", descr)
    descr = gsub("K\\(L\\)0", "K_L^0", descr)
    descr = gsub("(nu|tau|mu|pi|omega|gamma)", "\\\\\\1", descr)
    descr = gsub("\\s+\\(ex[.]", "\\\\;(\\\\text{ex.}", descr)
    descr = gsub("(neutrals)", "\\\\text{\\1}", descr)
    descr = gsub("(particles|strange|total)", "\\\\text{\\1}", descr)
    descr = gsub("(.*\\S)\\s*/\\s*(\\S.*)$", "\\\\frac{\\1}{\\2}", descr, perl=TRUE)
    texdescr = descr
    break
  }

  return(texdescr)
}

##
## return numeric value formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
get.tex.val = function(quant.val, precision, order) {
  quant.val = quant.val/10^order
  quant.err = quant.err/10^order
  if (order == 0) {
    rc = sprintf(paste("%", precision, "f", sep=""), quant.val)
  } else {
    rc = sprintf(paste("\\ensuremath{%", precision, "f\\cdot 10^{%d}}", sep=""), quant.val, order)
  }
  return(rc)
}

get.tex.val.auto = function(quant.val) {
  rc = get.precision.order(quant.val)
  precision = rc[1]
  order = rc[2]
  qvred = quant.val / 10^order
  if (qvred>=1) {
    precision = precision - 3.3
  } else {
    precision = precision - 2.2
  }
  quant.val = quant.val/10^order
  quant.err = quant.err/10^order
  if (order == 0) {
    rc = sprintf(paste("%", precision, "f", sep=""), quant.val)
  } else if (order == -2) {
    rc = sprintf(paste("%", precision, "f\\%%", sep=""), quant.val)    
  } else {
    rc = sprintf(paste("\\ensuremath{%", precision, "f\\cdot 10^{%d}}", sep=""), quant.val, order)
  }
  return(rc)
}

##
## return quantity formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
get.tex.quant.val = function(quant.val, quant.err, precision, order) {
  quant.val = quant.val/10^order
  quant.err = quant.err/10^order
  if (order == 0) {
    rc = sprintf(paste("\\(%", precision, "f \\pm %", precision, "f\\)", sep=""),
      quant.val, quant.err)
  } else {
    rc = sprintf(paste("\\((%", precision, "f \\pm %", precision, "f) \\cdot 10^{%d}\\)", sep=""),
      quant.val, quant.err, order)
  }
  return(rc)
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
get.tex.meas.val = function(meas, precision, order) {
  norm = 10^order
  str.val = sprintf(paste("%", precision, "f", sep=""), meas$value/norm)
  if (meas$stat.p == -meas$stat.n) {
    str.stat = sprintf(paste("\\pm %", precision, "f", sep=""), meas$stat/norm)
  } else {
    str.stat = sprintf(paste("{}^{%+", precision, "f}_{%+", precision, "f}", sep=""), meas$stat.p/norm, meas$stat.n/norm)
  }
  if (meas$syst.p == -meas$syst.n) {
    str.syst = sprintf(paste("\\pm %", precision, "f", sep=""), meas$syst/norm)
  } else {
    str.syst = sprintf(paste("{}^{%+", precision, "f}_{%+", precision, "f}", sep=""), meas$syst.p/norm, meas$syst.n/norm)
  }
  str.meas = paste(str.val, str.stat, str.syst)
  if (order == 0) {
    rc = paste("\\(", str.meas, "\\)", sep="")
  } else {
    rc = sprintf(paste("\\((", str.meas, ") \\cdot 10^{%d}\\) ", sep=""), order)
  }

  return(rc)
}

##
## for the spec. measurement return list with:
## - experiment
## - latex \cite{} reference
## - type: prelim or pub
## - formatter val +- stat +- syst
##
get.tex.meas = function(meas, precision, order) {
  meas.item = list()
  meas.item$exp = meas$tags[1]
  meas.item$ref = paste("\\cite{", get.reference(meas$tags), "}", sep="")
  meas.item$type = meas$tags[3]
  meas.item$value = get.tex.meas.val(meas, precision, order)
  return(meas.item)
}

##
## for the spec. quantity return the associated measurements
##
get.meas.quant = function(quant.name, delta) {
  delta.names = delta[, quant.name]
  meas.names = delta.names[delta.names!=0]
  return(names(meas.names))
}

##
## return tex label such as \Gamma_1 or \frac{\Gamma_1}{\Gamma_2}
##
alucomb2.gamma.texlabel.nv = function(gamma.name) {
  if (gamma.name == "GammaAll") return("1")
  str = str_match(gamma.name, "Gamma(\\d+)(by(\\d+))?")[1,]
  if (is.na(str[1])) return(gamma.name)
  if (str[4] == "") {
    return(paste("\\Gamma_{", str[2], "}", sep=""))
  }
  return(paste("\\frac{\\Gamma_{", str[2], "}}{\\Gamma_{", str[4], "}}", sep=""))
}
alucomb2.gamma.texlabel = Vectorize(alucomb2.gamma.texlabel.nv)

##
## for the spec vector of measurement, compute the
## appropriate precision and power-of-ten order
##
get.precision.order = function(vals) {
  max.abs.val = max(abs(vals))

  if (max.abs.val == 0) {
    not.higher.order = 0
  } else {
    not.higher.order = floor(log(max.abs.val*1.01)/log(10))
  }
  if (not.higher.order == -1) {
    order = -2
    precision = 5.2
  } else if (not.higher.order == -2) {
    order = -2
    precision = 5.3
  } else if (not.higher.order == -3) {
    order = -2
    precision = 5.3
  } else if (not.higher.order == 0) {
    order = 0
    precision = 5.3
  } else if (not.higher.order == 1) {
    order = 0
    precision = 5.2
  } else if (not.higher.order == 2) {
    order = 0
    precision = 5.1
  } else if (not.higher.order == 3) {
    order = 0
    precision = 5.0
  } else {
    order = not.higher.order
    precision = 5.3
  }
  precision = precision + 1.1

  return(c(precision=precision, order=order))
}

get.precision.order.quant = function(quant.name) {
  meas.names = get.meas.quant(quant.name, delta)
  
  vals = unlist(lapply(measurements[meas.names], function(m)
    c(m$value, m$stat.p, m$stat.n, m$syst.p, m$syst.n)))
  vals = c(vals, quant.val[quant.name], quant.err[quant.name])
  return(get.precision.order(vals))
}

##
## return body of latex tabular environment for the requested quantities
## including the quantities values followed by the related measurement values
##
get.tex.table = function(quant.names, with.meas=TRUE) {
  quant.order = order(alucomb2.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]
  rc = mapply(function(quant.name, quant) {
    rc = get.precision.order.quant(quant.name)
    precision = rc[1]
    order = rc[2]

    quant.descr = get.tex.quant.descr(quant)
    quant.descr = paste(alucomb2.gamma.texlabel(quant.name), "=", quant.descr)
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}\n", sep="")

    rc = paste(
      quant.descr, "&",
      get.tex.quant.val(quant.val[quant.name], quant.err[quant.name], precision, order),
      " & HFAG & Winter 2012 fit"
      )

    if (with.meas) {
      meas.names = get.meas.quant(quant.name, delta)
      meas.lines = sapply(meas.names, function(meas.name) {
        meas.item = get.tex.meas(measurements[[meas.name]], precision, order)
        return(paste(c("", meas.item$value, meas.item$exp, meas.item$ref), collapse=" & "))
      })
      if (length(meas.lines) > 0) {
        meas.lines = paste(meas.lines, collapse=" \\\\\n")
        rc = paste(rc, meas.lines, sep=" \\\\\n")
      }
    }
    return(rc)
  },
    quant.names, combination$quantities[quant.names])
  if (with.meas) {
    return(invisible(paste(rc, collapse=" \\\\\n\\hline\n")))
  } else {
    return(invisible(paste(rc, collapse=" \\\\\n")))
  }
}

get.tex.table.simple = function(quant.names, precision, order) {
  quant.order = order(alucomb2.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]
  rc = mapply(function(quant.name, quant) {
    quant.descr = get.tex.quant.descr(quant)
    quant.descr = paste(alucomb2.gamma.texlabel(quant.name), "=", quant.descr)
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}\n", sep="")    
    rc = paste(
      quant.descr, "&",
      get.tex.quant.val(quant.val[quant.name], quant.err[quant.name], precision, order)
      )
  }, quant.names, combination$quantities[quant.names])
  return(paste(rc, collapse=" \\\\\n"))
}

mkreport.tex.cmd = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{%\n", body, "%\n}\n", sep="")
}

mkreport.tex.cmd.short = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{", body, "\\xspace}\n", sep="")
}

##
## create files
##
mkreport = function(fname = "average2-aleph-hcorr.rdata") {
  load(fname, .GlobalEnv)

  fname = sub("[.][^.]*$", ".tex", fname, perl=TRUE)
  fname = file.path("../report", fname)
  cat("", file=fname)
  cat("file '", fname, "' created\n", sep="")

  tex.defs = paste(
    mkreport.tex.cmd.short("HfagTauMeasNum", as.character(meas.num)),
    ##--- subtract GammaAll and Gamma998
    mkreport.tex.cmd.short("HfagTauQuantNum", as.character(quant.num-2)),
    ##--- subtract GammaAll def and Unitarity constraint with Gamma998
    mkreport.tex.cmd.short("HfagTauConstrNum", as.character(constr.num-2)),
    mkreport.tex.cmd.short("HfagTauChisq", sprintf("%.1f", chisq)),
    mkreport.tex.cmd.short("HfagDof", as.character(dof)),
    mkreport.tex.cmd.short("HfagTauChisqProb", get.tex.val.auto(chisq.prob)),
    sep="")
  cat(tex.defs, file=fname, append=TRUE)
  cat("file '", fname, "', initial defs\n", sep="")

  quant.names = combination$combine
  quant.names = setdiff(quant.names, "GammaAll")
  tex.all.tau.br.val = mkreport.tex.cmd("HfagTauBrVal", get.tex.table(quant.names))
  cat(tex.all.tau.br.val, file=fname, append=TRUE)
  cat("file '", fname, "', BR val table content\n", sep="")

  gamma110.names = names(combination$constr.all.comb$Gamma110.c)
  gamma110.names = setdiff(gamma110.names, "Gamma110")
  tex.tau.br.strange.val = mkreport.tex.cmd("HfagTauBrStrangeVal", get.tex.table.simple(gamma110.names, 6.4, -2))
  cat(tex.tau.br.strange.val, file=fname, append=TRUE)
  cat("file '", fname, "', BR strange table content\n", sep="")

  gamma110.names = "Gamma110"
  tex.tau.br.strange.tot.val = mkreport.tex.cmd("HfagTauBrStrangeTotVal", get.tex.table.simple(gamma110.names, 6.4, -2))
  cat(tex.tau.br.strange.tot.val, file=fname, append=TRUE)
  cat("file '", fname, "', BR strange tot table content\n", sep="")

  ##
  ## dump correlation of base nodes
  ##
  tex.all.tau.br.corr = NULL  
  quant.names = combination$combine
  ##--- quantity names defined via constraint
  quant.constr.names = sub(".c", "", names(combination$constr.all.expr), fixed=TRUE)
  quant.names = setdiff(quant.names, quant.constr.names)
  quant.names = setdiff(quant.names, "Gamma998")
  quant.names = quant.names[order(alucomb2.gamma.num.id(quant.names))]

  corr.pre = c(
    "%%",
    "%% base nodes correlation, @@num@@",
    "%%",
    "\\ifhevea\\begin{table}\\fi%% otherwise cannot have normalsize caption",
    "\\begin{center}",
    "\\ifhevea",
    "\\caption{Base nodes correlation coefficients in percent, section @@num@@\\label{tab:br-fit-corr@@num@@}}%",
    "\\else",
    "\\begin{minipage}{\\linewidth}",
    "\\begin{center}",
    "\\captionof{table}{Base nodes correlation coefficients in percent, section @@num@@}\\label{tab:br-fit-corr@@num@@}%",
    "\\fi",
    "{\\ifhevea\\footnotesize\\else\\small\\fi",
    "\\renewcommand*{\\arraystretch}{1.1}%",
    "\\begin{tabular}{@@tabcols@@}",
    "\\hline")
  corr.post = c(
    "\\\\\\hline",
    "\\end{tabular}}",
    "\\ifhevea\\else",
    "\\end{center}",
    "\\end{minipage}",
    "\\fi",
    "\\end{center}",
    "\\ifhevea\\end{table}\\fi")
  
  quant.corr.base = quant.corr[quant.names, quant.names] * 100
  inum = 1
  for (j in seq(1, length(quant.names), by=10)) {
    for (i in seq(1, length(quant.names), by=10)) {
      submat.txt = NULL
      for(ii in i:min(i+10-1, length(quant.names))) {
        if (ii<=j) next
        maxcol = min(ii-1,j+10-1)
        row.txt = paste("\\(", alucomb2.gamma.texlabel(quant.names[ii]), "\\)")
        row.txt = c(row.txt, sprintf("%4.0f", quant.corr.base[ii, j:maxcol]))
        row.txt = c(row.txt, rep("", length.out=min(j+10-1, length(quant.names))-j+1-(maxcol-j+1)))
        submat.txt = c(submat.txt, paste(row.txt, collapse=" & "))
      }
      if (is.null(submat.txt)) next
      label.txt = alucomb2.gamma.texlabel(quant.names[j:min(j+10-1, length(quant.names))])
      label.num = min(j+10-1, length(quant.names)) - j + 1
      label.txt = paste("\\(", label.txt, "\\)", sep=" ", collapse=" & ")
      label.txt = paste("", label.txt, sep=" & ")
      submat.txt = paste(submat.txt, collapse=" \\\\\n")
      submat.txt = paste(submat.txt, label.txt, sep=" \\\\\n")
      corr.pre.mod = gsub("@@num@@", as.character(inum), corr.pre)
      corr.pre.mod = gsub("@@tabcols@@", paste(rep("r", length.out=label.num+1), collapse=""), corr.pre.mod)
      ## corr.pre.mod = paste(corr.pre.mod, collapse="\n")
      tex.all.tau.br.corr = c(tex.all.tau.br.corr, corr.pre.mod, submat.txt, corr.post)
      inum = inum + 1
    }
  }
  tex.all.tau.br.corr = mkreport.tex.cmd("HfagTauBrCorr", paste(tex.all.tau.br.corr, collapse="\n"))
  cat(tex.all.tau.br.corr, file=fname, append=TRUE)
  cat("file '", fname, "', BR correlations table content\n", sep="")
  
  if (FALSE) {
    rc = apply(quant.corr[quant.names, quant.names], 1, function(corr.val) {
      rc = paste(sprintf("%4.0f", 100*corr.val), collapse = " & ")
    })
    quant.texlabel = alucomb2.gamma.texlabel(quant.names)
    rc = paste(quant.texlabel, "&", rc)

    fname = "../report/tau-aleph-hcorr-corr.tex"
    corr.header = paste(quant.texlabel, collapse=" & ")
    cat(file=fname, c(corr.header, rc), sep=" \\\\")
    cat("produced file '", fname, "'\n", sep="")
  }

  ##
  ## dump constraint equations
  ##
  constr.order = order(alucomb2.gamma.num.id(names(combination$constr.all.str.expr)))
  comb.str = combination$constr.all.str.expr[constr.order]
  comb.str = comb.str[grep("Gamma\\d+by\\d+.*|Unitarity", names(comb.str), perl=TRUE, invert=TRUE)]
  comb.str.nl = intersect(names(comb.str), names(combination$constr.nl.str.expr))

  comb.str[comb.str.nl] = combination$constr.nl.str.expr[comb.str.nl]
  rc = str_match(unlist(comb.str), "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")
  constr.left = rc[,2]
  constr.right = rc[,3]
  constr.left = gsub("Gamma(\\d+)", "\\\\Gamma_{\\1}", constr.left, perl=TRUE)
  constr.left = gsub("GammaAll", "\\Gamma_{\\text{All}}", constr.left, fixed=TRUE)
  constr.right = sub("\\s*[+]\\s*", " + ", constr.right, perl=TRUE)
  constr.right = gsub("Gamma(\\d+)", "\\\\Gamma_{\\1}", constr.right, perl=TRUE)
  constr.right = gsub("GammaAll", "\\Gamma_{\\text{All}}", constr.right, fixed=TRUE)

  ##--- remove outer braces
  constr.right = gsub("(\\(((?:[^()]++|(?1))+)*+\\))", "\\2", constr.right, perl=TRUE)
  ##--- split long eqs
  constr.right = gsub("(([^+]+[+*]){4}[^+]+)[+]", "\\1 \\\\\\\\ \n  {}& +", constr.right, perl=TRUE)

  constr.right = gsub("*", "\\cdot{}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_eta_2gam", "\\Gamma_{\\eta\\to\\gamma\\gamma}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_eta_neutral", "\\Gamma_{\\eta\\to\\text{neutral}}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_eta_3piz", "\\Gamma_{\\eta\\to3\\pi^0}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_eta_pimpippiz", "\\Gamma_{\\eta\\to\\pi^+\\pi^-\\pi^0}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_eta_charged", "\\Gamma_{\\eta\\to\\text{charged}}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_KS_2piz", "\\Gamma_{K_S\\to\\pi^0\\pi^0}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_KS_pimpip", "\\Gamma_{K_S\\to\\pi^+\\pi^-}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_om_pimpippiz", "\\Gamma_{\\omega\\to\\pi^+\\pi^-\\pi^0}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_om_pimpip", "\\Gamma_{\\omega\\to\\pi^+\\pi^-}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_om_pizgamma", "\\Gamma_{\\omega\\to\\pi^0\\gamma}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_phi_KmKp", "\\Gamma_{\\phi\\to K^+K^-}", constr.right, fixed=TRUE)
  constr.right = gsub("BR_phi_KSKL", "\\Gamma_{\\phi\\to K_S K_L}", constr.right, fixed=TRUE)
  constr.right = gsub("BRA_Kz_KS_KET", "\\Gamma_{<K^0|K_S>}", constr.right, fixed=TRUE)
  constr.right = gsub("BRA_Kz_KL_KET", "\\Gamma_{<K^0|K_L>}", constr.right, fixed=TRUE)
  constr.right = gsub("BRA_Kzbar_KS_KET", "\\Gamma_{<\\bar{K}^0|K_S>}", constr.right, fixed=TRUE)
  constr.right = gsub("BRA_Kzbar_KL_KET", "\\Gamma_{<\\bar{K}^0|K_L>}", constr.right, fixed=TRUE)
  constr.right = gsub("BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET", "R_{0\\bar{0}SS/LL}", constr.right, fixed=TRUE)

  tex.constr.val = mkreport.tex.cmd("HfagConstrVal",
    paste("\\begin{align*}\n", constr.left, "={}&", constr.right, "\n\\end{align*}", collapse="\n"))
  cat(tex.constr.val, file=fname, append=TRUE)
  cat("file '", fname, "', constraint table content\n", sep="")
}

## ////////////////////////////////////////////////////////////////////////////
## code

args <- commandArgs(TRUE)
if (length(args) == 1) {
  rc = mkreport(fname = args[1])
} else {
  mkreport()
}
