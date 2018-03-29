"""
"""
import os
print "importing from ", os.path.abspath(__path__[0])
#from GaudiPython.Bindings import gbl
import cppyy
csm_template = cppyy.gbl.csm_template
csm_model = cppyy.gbl.csm_model
mclimit_csm = cppyy.gbl.mclimit_csm
csm = cppyy.gbl.csm
