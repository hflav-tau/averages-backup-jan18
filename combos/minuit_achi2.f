************************************************************************
*
*     COMBOS combination routine to perform a NCSYS+1 parameter
*     chi2 fit using MINUIT, the (NCSYS+1)th parameter being the
*     combined value of the NMEAS measurements
*
*     This routine prepares for fits using asymmetric errors
*     and is based on MINUIT_CHI2 of O. Schneider.
*     It also includes the option to calculate minos errors on the
*     fitted parameters.
*
*     Version 2.40, October 14, 1997:
*     -- original release (HCJS)
*
************************************************************************
*
       SUBROUTINE MINUIT_ACHI2(FCN_INIT,FCN,CVAL,ERR2P,ERR2N,CL,IERR,
     &   ICASE)
*     ===============================================================
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
*              ICASE   = 1 -- combination with systematics
*                        2 -- no systematics in combination
*
      IMPLICIT NONE
*
*     Arguments
*      
      EXTERNAL FCN_INIT,FCN
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR
      INTEGER ICASE
*
*     Externals
*
      REAL PROB
      INTEGER LENOCC
*
*     Control Commons
*
      INCLUDE 'master.inc'
      INCLUDE 'combos.inc'


*     Local variables
*
      CHARACTER*10 CHMNAM
      INTEGER I,NPARI,NPARX,ISTAT,IVARBL
      DOUBLE PRECISION FMIN,FEDM,ERRDEF,DUMMY,PRTLEV,START,CVALO
      INTEGER MITER,NITER,SPECIAL_FLAG
      PARAMETER(MITER=50)     ! maximum iterations on special FCN calls
      DOUBLE PRECISION ERR2PB ! Parabolic error
      LOGICAL XMINOS_DONE     ! True if minos errors allready calculated
      SAVE XMINOS_DONE
      DOUBLE PRECISION MINARG(MNAM) ! Arguments for MINOS
      INTEGER NMINARG               ! Num of Arguments
      LOGICAL XMINOSALL             ! Minos for all params
      INTEGER ICINQ,JMINOS          ! Character search function
*
*
*     Default output arguments
*
      CVAL  = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      ERR2PB= 0.0D0
      CL    = -1.D0
      IERR  = 9999

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
* Prepare MINOS argument list
      MINARG(1) = 5000.0D0      ! Number of calls.
* Check if MINOS for all params
      IF (ICINQ('ALL',CHMINOS,NMINOS) .GT. 0) THEN
        NMINARG = 1
        XMINOSALL = .TRUE.
      ELSE
        XMINOSALL = .FALSE.
        NMINARG = 1
      ENDIF

      DO I=1,NCSYS+1
        IF(I.LE.NCSYS) THEN
          CALL MNPARM(I,CHPARA(I),PARA(I),0.1D0*EXCU(I),
     &                0.D0,0.D0,IERR)
          IF (.NOT. XMINOSALL) THEN
*            Get the parameter numbers of those paramters for 
*            which MINOS errors are to be caclulcated
            JMINOS = ICINQ(CHPARA(I),CHMINOS,NMINOS)
            IF (JMINOS .GT. 0) THEN
              NMINARG = NMINARG + 1
              MINARG(NMINARG) = I
            ENDIF
          ENDIF
        ELSE
          CALL MNPARM(I,CHMEAS,START,DABS(0.02D0*START),0.D0,0.D0,IERR)
          IF (.NOT. XMINOSALL) THEN
*            Get the parameter numbers of those paramters for 
*            which MINOS errors are to be caclulcated
            JMINOS = ICINQ(CHMEAS,CHMINOS,NMINOS)
            IF (JMINOS .GT. 0) THEN
              NMINARG = NMINARG + 1
              MINARG(NMINARG) = I
            ENDIF
          ENDIF
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
* MINOS, perform only for full combination with systematics
*
      IF (XMINOS.AND.NCSYS.GT.0) THEN
        CALL EXECUTE(FCN,'SET PRINT',0.D0,1,IERR,CHROUT,LUNIT)
        CALL EXECUTE(FCN,'MINOS',MINARG,NMINARG,IERR,CHROUT,LUNIT)
        IF(IERR.NE.0) RETURN
        CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
        IF(ISTAT.NE.3) THEN 
          IF(LUNIT.GT.0) WRITE(LUNIT,1002) CHROUT(:LENOCC(CHROUT)),ISTAT
 1002     FORMAT(1X,A,': no proper MINOS convergence ! ISTAT=',I6)
          IERR=ISTAT
          IF(IERR.EQ.0) IERR=-9999
          RETURN
        ENDIF 
        XMINOS_DONE = .TRUE.
        CALL EXECUTE(FCN,'SET PRINT',4.D0,1,IERR,CHROUT,LUNIT)
      ELSE
        XMINOS_DONE = .FALSE.
      ENDIF

*
*     Get result of fit, and minos errors.
*    
      CALL MNPOUT(NCSYS+1,CHMNAM,CVAL,ERR2PB,DUMMY,DUMMY,IVARBL)
      ERR2PB=ERR2PB**2
      CALL MNERRS(NCSYS+1,ERR2P,ERR2N,ERR2PB,DUMMY)
      ERR2P=ERR2P**2
      ERR2N=ERR2P**2
      IF (ERR2P.EQ.0.0D0) ERR2P = ERR2PB
      IF (ERR2N.EQ.0.0D0) ERR2N = ERR2PB

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
        IF (.NOT. xminos_done .AND. XMINOS) THEN
C For the special case do minos even if ncsys == 0
          CALL EXECUTE(FCN,'SET PRINT',0.D0,1,IERR,CHROUT,LUNIT)
          CALL EXECUTE(FCN,'MINOS',MINARG,NMINARG,IERR,CHROUT,LUNIT)
          IF(IERR.NE.0) RETURN
          CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
          IF(ISTAT.NE.3) THEN 
            IF(LUNIT.GT.0) WRITE(LUNIT,1002)
     $         CHROUT(:LENOCC(CHROUT)),ISTAT
            IERR=ISTAT
            IF(IERR.EQ.0) IERR=-9999
            RETURN
          ENDIF 
          CALL MNERRS(NCSYS+1,ERR2P,ERR2N,ERR2PB,DUMMY)
          ERR2P=ERR2P**2
          ERR2N=ERR2N**2
          IF (ERR2P.EQ.0.0D0) ERR2P = ERR2PB
          IF (ERR2N.EQ.0.0D0) ERR2N = ERR2PB
          xminos_done = .TRUE.
        ENDIF
        IF(NITER.EQ.MITER) THEN
          CALL EXECUTE(FCN,'SET PRINT',1.D0,1,IERR,CHROUT,LUNIT)
          IF(IERR.NE.0) RETURN
          IF(LUNIT.GT.0) WRITE(LUNIT,*) ' '
        ENDIF
        IF(LUNIT.GT.0) WRITE(LUNIT,1003) CHMNAM,CVAL,DSQRT(ERR2PB),
     &                                   DSQRT(ERR2P),DSQRT(ERR2N)
 1003   FORMAT(1X,A,' = ',F14.10,' +-',F14.10,' +',F14.10,F14.10)
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
      IF (ICASE.EQ.1) THEN
* Now MINOS, only for full combination with systematics
        CALL MNERRS(NCSYS+1,ERR2P,ERR2N,ERR2PB,DUMMY)
        ERR2P=ERR2P**2
        ERR2N=ERR2N**2
        IF (ERR2P.EQ.0.0D0) ERR2P = ERR2PB
        IF (ERR2N.EQ.0.0D0) ERR2N = ERR2PB
        IF(NMEAS.GT.1) CL=DBLE(PROB(SNGL(FMIN),NMEAS-1)) ! is this correct ?
        CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
        IF(IERR.NE.0) RETURN
        CALL EXECUTE(FCN,'CALL FCN',3.D0,1,IERR,CHROUT,LUNIT)
      ELSE
        ERR2P=ERR2PB   ! return symmetric error
        ERR2N=ERR2PB   ! return symmetric error
        IF(NMEAS.GT.1) CL=DBLE(PROB(SNGL(FMIN),NMEAS-1)) ! is this correct ?
        CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
        IF(IERR.NE.0) RETURN
        CALL EXECUTE(FCN,'CALL FCN',3.D0,1,IERR,CHROUT,LUNIT)
      ENDIF

      END
*

************************************************************************
*
      SUBROUTINE PREPARE_ACHI2(LCSYS,IERR)
*     =====================================
*
*     Prepare for assymmetric error chi2 fit.
*     Return number of unknown parameters
*     Safety checks on array bounds etc.
*
      IMPLICIT NONE
      INCLUDE 'master.inc'
*
*     Arguments
*
      INTEGER LCSYS,IERR
*
*     Externals
*
      INTEGER LENOCC
*
*     Check that the number of correlated systematic uncertainties is 
*     less than the maximum
*
      LCSYS=NCSYS+1
      IF(LCSYS.GT.MCSYS) THEN
        IERR=MCSYS
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'parameter MCSYS too small',IERR
 1000   FORMAT(1X,A,': ',A,' ! IERR=',I6)
        RETURN
      ENDIF
      END

