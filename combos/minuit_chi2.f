************************************************************************
*
*     COMBOS combination routine to perform a NCSYS+1 parameter
*     chi2 fit using MINUIT, the (NCSYS+1)th parameter being the
*     combined value of the NMEAS measurements
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.20, December 13, 1996:
*     -- original release
*     Version 1.30, January 27, 1997:
*     -- extra argument SPECIAL_FLAG in FCN_INIT
*     -- allow for iterations on calls to FCN with a special argument
*     Version 2.00, February 20, 1997:
*     -- first argument of MINUIT_CHI2 removed (CHROUT now in master.inc)
*     -- parameter MINUIT_PRTLUN not used anymore (LUNIT now in master.inc)
*     -- default MINUIT printout level changed from -1 to +1 
*     -- improved printout (printout turned off if LUNIT is zero or less)
*     -- two first arguments of PREPARE_CHI2 and FCN_INIT removed
*         (CHROUT and LUNIT now in master.inc)
*     Version 2.30, May 12, 1997:
*     -- use NMEFF instead of NMEAS in computation of chi2 probability
*
************************************************************************
*
      SUBROUTINE MINUIT_CHI2(FCN_INIT,FCN,CVAL,ERR2P,ERR2N,CL,IERR)
*     =============================================================
*
*     Olivier Schneider, CERN/PPE-ALE
*     December 11, 1996
*
*     Input:   FCN_INIT = routine to be called before fit for initialization
*     ------   FCN      = "FCN" routine for MINUIT
*      
*     Output:  CVAL    = combined value (result of fit)
*     -------  ERR2P   = positive error**2 on CVAL
*              ERR2N   = negative error**2 on CVAL
*              CL      = confidence level of fit
*              IERR    = error flag (0 means OK)
*
      IMPLICIT NONE
*
*     Arguments
*      
      EXTERNAL FCN_INIT,FCN
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR
*
*     Externals
*
      REAL PROB
      INTEGER LENOCC
*
*     Local variables
*
      CHARACTER*10 CHNAM
      INTEGER I,NPARI,NPARX,ISTAT,IVARBL
      DOUBLE PRECISION FMIN,FEDM,ERRDEF,DUMMY,PRTLEV,START,CVALO
      INTEGER MITER,NITER,SPECIAL_FLAG
      PARAMETER(MITER=50) ! maximum number of iterations on special FCN calls
*
      INCLUDE 'master.inc'
*
*     Default output arguments
*
      CVAL = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      CL   = -1.D0
      IERR = 9999
*
*     Determine MINUIT printout level
*
      PRTLEV=1.D0 ! default printout level
      DO I=NCSYS+1,NPARA
        IF(CHPARA(I).EQ.'MINUIT_PRTLEV') PRTLEV=PARA(I)
      ENDDO
      IF(LUNIT.LE.0) THEN
        PRTLEV=-1.D0 ! switch MINUIT output off
      ELSE
        WRITE(LUNIT,2000) CHROUT(:LENOCC(CHROUT)),
     &                    CHMEAS(:LENOCC(CHMEAS)),
     &                    CHCOMB(:LENOCC(CHCOMB))
        IF(CHSTEP.NE.' '.OR.ISTEP.NE.1) WRITE(LUNIT,2001) ISTEP,
     &                    CHSTEP(:LENOCC(CHSTEP)),VSTEP
        WRITE(LUNIT,2002)
      ENDIF
 2000 FORMAT(/,1X,78('='),
     &       /,1X,A,':',T20,'MINUIT fit to find average ',A,
     &       /,1X,      T20,'for ',A)
 2001 FORMAT(  1X,      T20,'Step',I4,' at ',A,' = ',F10.4)
 2002 FORMAT(  1X,78('='),/)
*
*     Initialize MINUIT
*
      CALL MNINIT(0,IABS(LUNIT),0)
      CALL MNSETI(CHROUT(:LENOCC(CHROUT))//
     &            ': MINUIT fit to find average '//CHMEAS)
      CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
      IF(IERR.NE.0) RETURN
      CALL MNEXCM(0,'SET ERROR',1.D0,1,IERR,0) ! set ERRORDEF for chi2 fit
      IF(IERR.NE.0) RETURN
      CALL EXECUTE(FCN,'SET PRINT',PRTLEV,1,IERR,CHROUT,LUNIT)
      IF(IERR.NE.0) RETURN
*
*     Parameter definitions
*
      START=0.D0
      DO I=1,NMEAS
        START=START+MEAS(I)
      ENDDO
      START=START/DBLE(NMEAS)

 
*DR allow NQUAN quantities to be fitted     DO I=1,NCSYS+1
      DO I=1,NCSYS+NQUAN
        IF(I.LE.NCSYS) THEN
          CALL MNPARM(I,CHPARA(I),PARA(I),0.1D0*EXCU(I),
     &                0.D0,0.D0,IERR)
        ELSE
*DR          CALL MNPARM(I,CHMEAS,START,0.02D0*START,0.D0,0.D0,IERR)
          CALL MNPARM(I,CHQUAN(I-NCSYS),START,0.02D0*START,
     &           0.D0,0.D0,IERR)
        ENDIF
        IF(IERR.NE.0) THEN
          IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                   'MNPARM','IERR',IERR
 1000     FORMAT(1X,A,': error returned by ',A,' ! ',A,'=',I6)
          RETURN
        ENDIF
      ENDDO
*
*     Initialize FCN routine
*
      CALL FCN_INIT(SPECIAL_FLAG,IERR)
      IF(IERR.NE.0) RETURN 
      CVALO=0.D0
      IF(SPECIAL_FLAG.GT.0) THEN 
        NITER=MITER
      ELSE
        NITER=0
      ENDIF
*
*     Perform minimization
*
    1 CALL EXECUTE(FCN,'MINIMIZE',0.D0,0,IERR,CHROUT,LUNIT)

      IF(IERR.NE.0) RETURN
      CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
      IF(ISTAT.NE.3) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,1001) CHROUT(:LENOCC(CHROUT)),ISTAT
 1001   FORMAT(1X,A,': no proper convergence ! ISTAT=',I6)
        IERR=ISTAT
        IF(IERR.EQ.0) IERR=-9999
        RETURN
      ENDIF
*
*     Get result of fit
*
      CALL MNPOUT(NCSYS+1,CHNAM,CVAL,ERR2P,DUMMY,DUMMY,IVARBL)
      ERR2P=ERR2P**2
      IF(IVARBL.LE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'MNPOUT','IVARBL',IVARBL
        IERR=IVARBL 
        IF(IERR.EQ.0) IERR=9999
        RETURN
      ENDIF
*
*     Iterate
*
      IF(SPECIAL_FLAG.GT.0) THEN
        IF(NITER.EQ.MITER) THEN
          CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
          IF(IERR.NE.0) RETURN
          IF(LUNIT.GT.0) WRITE(LUNIT,*) ' '
        ENDIF
        IF(LUNIT.GT.0) WRITE(LUNIT,1003) CHNAM,CVAL,DSQRT(ERR2P)
 1003   FORMAT(1X,A,' = ',F14.10,' +-',F14.10)
        IF(DABS(CVAL-CVALO).GT.1.D-8) THEN
          CVALO=CVAL
          CALL EXECUTE(FCN,'CALL FCN',
     &                 DBLE(SPECIAL_FLAG),1,IERR,CHROUT,LUNIT)
          IF(IERR.NE.0) RETURN
          NITER=NITER-1
          IF(NITER.GT.0) GOTO 1
        ENDIF
        IF(NITER.GT.0) THEN 
          IF(LUNIT.GT.0) WRITE(LUNIT,3000) 
     &                   'Converged',MITER-NITER,SPECIAL_FLAG
 3000     FORMAT(/,1X,A,' after',I3,
     &           ' special calls to FCN with IFLAG =',I3,/)
        ELSE IF(LUNIT.GT.0) THEN
          WRITE(LUNIT,3000) 'Stop without convergence',
     &                      MITER-NITER,SPECIAL_FLAG
        ENDIF
        CALL EXECUTE(FCN,'SET PRINT',PRTLEV,1,IERR,CHROUT,LUNIT)
        IF(IERR.NE.0) RETURN
        SPECIAL_FLAG=0
        GOTO 1
      ENDIF
*
*     Finalize
*
      ERR2N=ERR2P ! return symmetric error
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(FMIN),NMEFF-1)) ! is this correct ?
      CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
      IF(IERR.NE.0) RETURN
      CALL EXECUTE(FCN,'CALL FCN',3.D0,1,IERR,CHROUT,LUNIT)
      END
*
************************************************************************
*
      SUBROUTINE EXECUTE(FCN,CHCOM,ARGLIS,NARG,IERR,CHROUT,LUNIT)
*     ===========================================================
*
*     Execute MINUIT command
*
      IMPLICIT NONE
*
*     Arguments
*
      EXTERNAL FCN
      CHARACTER*(*) CHCOM,CHROUT
      DOUBLE PRECISION ARGLIS(*)
      INTEGER NARG,IERR,LUNIT
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I
*
      CALL MNEXCM(FCN,CHCOM,ARGLIS,NARG,IERR,0)
      IF(IERR.NE.0.AND.LUNIT.GT.0) WRITE(LUNIT,1000)
     & CHROUT(:LENOCC(CHROUT)),IERR,CHCOM,(ARGLIS(I),I=1,NARG)
 1000 FORMAT(1X,A,': error',I6,' returned by ',A,10G10.4)
      END
*
************************************************************************
*
      SUBROUTINE PREPARE_CHI2(LCSYS,W,IERR)
*     =====================================
*
*     Prepare for chi2 fit.
*     Return number of unknown parameters and inverse of error matrix.
* DR Version 3.0 1 June 1999 modif: allow more than one param
      IMPLICIT NONE
      INCLUDE 'master.inc'
*
*     Arguments
*
      INTEGER LCSYS,IERR
      DOUBLE PRECISION W(MMEAS,MMEAS)
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I,J
*
*     Check that the number of correlated systematic uncertainties is 
*     less than the maximum
*
*DR      LCSYS=NCSYS+1
      LCSYS=NCSYS+NQUAN
      IF(LCSYS.GT.MCSYS) THEN
        IERR=MCSYS
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'parameter MCSYS too small',IERR
 1000   FORMAT(1X,A,': ',A,' ! IERR=',I6)
        RETURN
      ENDIF
*
*     Fill an additional row (row LCSYS) of matrix CSYS with -1's; 
*     this is against the philosophy that a combination routine should
*     not write into the common blocks listed in master.inc, but here it's OK
*     because:
*        - we write in a location which is not supposed to be read by other
*          combination routines
*        - routine MASTER knows about it, and will restore the overwritten 
*          stuff if necessary
*
*DR fill -1 and 0 at the end of the matrix
*      DO I=1,NMEAS
*        CSYS(I,LCSYS)=-1.D0
*      ENDDO
*DR
      IF (NQUAN.EQ.0) THEN
        CALL COMBOS_ERROR(-1,
     & 'MINUIT_CHI2:No quantities to be fit','cannot do anything')
      ENDIF
      DO I=1,NMEAS
        DO J=1,NQUAN
          IF (KQUAN(I).EQ.J) THEN
            CSYS(I,NCSYS+J)=-1.D0
          ELSE
            CSYS(I,NCSYS+J)=0.D0
          ENDIF
        ENDDO
      ENDDO
*DR end 
*
*     Build total error matrix
*
      DO I=1,NMEAS
        DO J=1,NMEAS
          W(I,J)=STACOR(I,J)*STAT(I)*STAT(J) ! statistical error matrix
        ENDDO
*DR XXXXX add correlated systematics here
        W(I,I)=W(I,I)+USYS(I)**2 ! add uncorrelated systematic error
      ENDDO


*
*     Invert total error matrix
*
      CALL DSINV(NMEAS,W,MMEAS,IERR)
      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'cannot invert error matrix',IERR
        RETURN 
      ENDIF 
*
*     At this point, the correspondance between the matrices and vectors used
*     in this routine and in the documentation (COMBOS note) is the following:
*
*                               In this routine     In the COMBOS note
*                               ---------------     ------------------
*     Number of measurements:   NMEAS               n
*     Vector of measurements:   MEAS                X
*     Inverse of error matrix:  W                   M**(-1)
*     Number of unknowns:       LCSYS               l
*     Vector of unknowns:       V                   V
*     Matrix of changes:        CSYS                Delta
*
      END
