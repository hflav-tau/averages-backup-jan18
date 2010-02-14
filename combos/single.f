************************************************************************
*
*     Dummy combination routine of the COMBOS program
*     (can only combine  a single measurement)
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.30, May 12, 1997
*     -- original release
*
************************************************************************
*
      SUBROUTINE SINGLE(CVAL,ERR2P,ERR2N,CL,IERR)
*     ===========================================
*
*     Input data: available from a common block in file master.inc
*     ----------
*
      INCLUDE 'master.inc'
*      
*     Output data: returned in the arguments 
*     -----------
*
      DOUBLE PRECISION CVAL  ! combined value (or average value)
      DOUBLE PRECISION ERR2P ! positive total uncertainty**2 on CVAL
      DOUBLE PRECISION ERR2N ! negative total uncertainty**2 on CVAL
      DOUBLE PRECISION CL    ! confidence level of combination (e.g. chi2 prob.)
*                              or -1 if confidence level not available
      INTEGER IERR           ! error flag (non-zero if combination failed)
*
*     Note: this routine should only write output (including error messages)
*     ----  to logical unit LUNIT (variable LUNIT is available in a common block
*           described in master.inc); if LUNIT is zero or less, then this
*           routine should not write any output
*
********
*      
      INTEGER J
*
      IF(NMEAS.EQ.1) THEN
        IERR=0
        CVAL = MEAS(1)
        ERR2P = STATP(1)**2+USYSP(1)**2
        ERR2N = STATN(1)**2+USYSN(1)**2
        CL   = -1.D0
        DO J=1,NCSYS
          ERR2P=ERR2P+DMAX1(DMAX1(CSYSP(1,J),CSYSN(1,J)),0.D0)**2
          ERR2N=ERR2N+DMIN1(DMIN1(CSYSP(1,J),CSYSN(1,J)),0.D0)**2
        ENDDO
      ELSE ! combination failed, return dummy values
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'SINGLE: cannot combine more than 1 measurement'
        IERR=-1
        CVAL = 0.D0
        ERR2P = 0.D0
        ERR2N = 0.D0
        CL   =-1.D0
      ENDIF
      END
