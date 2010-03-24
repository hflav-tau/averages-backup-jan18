alias aluPathAppIfDU 'if ( -d \!:2 && $\!:1 !~ \!:2\:* && $\!:1 !~ *\:\!:2\:* && $\!:1 !~ *\:\!:2 ) setenv \!:1 ${\!:1}\:\!:2'
# --- append element to path if it does not exist, and if it is a directory
alias aluPathAppIfD 'if ( $?\!:1 ) eval aluPathAppIfDU \!:1 \!:2 ; if ( -d \!:2 && ! $?\!:1 ) setenv \!:1 \!:2'

aluPathAppIfD PATH /opt/TWWfsw/bin
