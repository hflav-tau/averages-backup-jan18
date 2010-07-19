#
# source config.csh to get access to:
# - R on iris and flora at SLAC
# - root V5.26 on flora and yakut at SLAC
#

alias aluPathAppIfDU 'if ( -d \!:2 && $\!:1 !~ \!:2\:* && $\!:1 !~ *\:\!:2\:* && $\!:1 !~ *\:\!:2 ) setenv \!:1 ${\!:1}\:\!:2'
# --- append element to path if it does not exist, and if it is a directory
alias aluPathAppIfD 'if ( $?\!:1 ) eval aluPathAppIfDU \!:1 \!:2 ; if ( -d \!:2 && ! $?\!:1 ) setenv \!:1 \!:2'

alias aluPathPrepIfDU 'if ( -d \!:2 && $\!:1 !~ \!:2\:* && $\!:1 !~ *\:\!:2\:* && $\!:1 !~ *\:\!:2 ) setenv \!:1 \!:2\:${\!:1}; if ( -d \!:2 && $\!:1 !~ \!:2\:* ) setenv \!:1 \!:2\:`echo ${\!:1} | sed -e s%^\!:2\:%% -e s%:\!:2\:%:%g -e s%:\!:2\$%%`'
alias aluPathPrepIfD 'if ( $?\!:1 ) eval aluPathPrepIfDU \!:1 \!:2 ; if ( -d \!:2 && ! $?\!:1 ) setenv \!:1 \!:2'

aluPathAppIfD PATH /opt/TWWfsw/bin
if (`uname -s` == "SunOS") then
    aluPathPrepIfD PATH /afs/slac.stanford.edu/g/babar/package/root/5.26-00/SunOS5v10_sparc_Studio10/bin
endif

expr `uname -r` : '2\.4\.21\.*' >/dev/null
set sl3=$status
if ($sl3 == 0 && `uname -s` == "Linux") then
    aluPathPrepIfD PATH /afs/slac.stanford.edu/g/babar/package/root/5.26-00/Linux24SL3_i386_gcc323/bin
endif

alias make gmake

setenv ROOTSYS /afs/slac.stanford.edu/g/babar/package/root/5.26-00/Linux26SL5_i386_gcc412
if ( ${?LD_LIBRARY_PATH} ) then
   setenv LD_LIBRARY_PATH $ROOTSYS/lib:/usr/local/lib:/cern/2002/lib:${PWD}/tmp
else
   setenv LD_LIBRARY_PATH "/usr/lib:/usr/local/lib:${ROOTSYS}/lib:/cern/2002/lib"
endif
aluPathAppIfD PATH $ROOTSYS/bin
