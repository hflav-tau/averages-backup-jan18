************************************************************************
*
*     "CHI2_ASYM_MIN" combination routine for COMBOS program
*
*     Henry Seywerd -- Based on CHI2_SYM_MIN from
*     Olivier Schneider, CERN/PPE-ALE
*
*     First attempt to deal with asymmetric errors for combinations
*
*     Version 2.40, October 14, 1997:
*     -- original release (HCJS)
*
************************************************************************
*
      SUBROUTINE CHI2_ASYM_MIN(CVAL,ICASE,ERR2P,ERR2N,CL,IERR)
*     ========================================================
*
*     (HCJS)
*     May 1997
*
*     Combination routine performing a chi2 fit in MINUIT:
*       - can handle correlated statistical uncertainties with any
*          correlation coefficient
*       - can handle 100% or 0% systematic correlations
*       - asymmertric errors
*
      IMPLICIT NONE
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR, ICASE
      EXTERNAL FCN_CHI2_ASYM_MIN_INIT, FCN_CHI2_ASYM_MIN
*
      CALL MINUIT_ACHI2(FCN_CHI2_ASYM_MIN_INIT,FCN_CHI2_ASYM_MIN,
     &                 CVAL,ERR2P,ERR2N,CL,IERR,ICASE)
      END
*
************************************************************************
*
      SUBROUTINE FCN_CHI2_ASYM_MIN(NPAR,D,F,X,IFLAG,FUNC)
*     ==================================================
*
*     Routine to return value of function to be minimized by MINUIT.
*
*     Input:  NPAR  = number of parameters
*     -----   X     = values of parameters
*             IFLAG = input flag
*                     IFLAG=1 means initialization has to be performed
*                     IFLAG=2 means derivatives have to be computed
*                     IFLAG=3 means finalization has to be performed
*
*     Output: D     = derivatives (not computed in this routine)
*     ------  F     = value of function to be minimized
*
*
      IMPLICIT NONE
*
*     Arguments
*
      EXTERNAL FUNC ! not used
      DOUBLE PRECISION F,D(*),X(*)
      INTEGER NPAR, IFLAG
*
*     Argument of entry FCN_CHI2_ASYM_MIN_INIT
*
      INTEGER SPECIAL_FLAG,IERR
*
      INCLUDE 'master.inc'
*
*     External
* 
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I,J,IVARBL
      CHARACTER*10 CHNAM
      CHARACTER*16 CHNAME
      CHARACTER*1 COMM(2)
      DOUBLE PRECISION VAL,ERR2PB
      DOUBLE PRECISION ERR2P,ERR2N
      DOUBLE PRECISION CHI2,DUMMY
      INTEGER LCSYS
      SAVE LCSYS
      DATA LCSYS/0/
      DOUBLE PRECISION W(MMEAS,MMEAS), Z(MMEAS),V(MCSYS)
      SAVE             W
      LOGICAL NOT_READY
      SAVE    NOT_READY
      DATA    NOT_READY/.TRUE./
C Assymmetric error array
      DOUBLE PRECISION ASYS(MMEAS,MCSYS)
      INTEGER JERR
*
*     Check that routine was initialized
*
      IF(NOT_READY) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ' not initialized !'
        STOP 566
      ENDIF
*
*     Determine negative/positive errors based upon
*     relative value of parameters and current fitted values
*
      CALL ASYM_CALC_ERRS(LCSYS, X, ASYS)
*
*
*     Compute derivatives
*
      IF(IFLAG.EQ.2) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ': cannot compute derivatives'
        STOP 669
      ENDIF
*
*     Compute vector V = rescaled parameters
*     Rescale by positive/negative excursions.
*
      DO I=1,NCSYS
        V(I)=(X(I)-PARA(I))/EXCU(I)
*        IF (X(I) .GT. PARA(I)) THEN
*          V(I)=(X(I)-PARA(I))/EXCUP(I)
*        ELSE
*          V(I)=(X(I)-PARA(I))/EXCUN(I)
*        ENDIF
      ENDDO
      V(LCSYS)=X(LCSYS)
*
*     Compute Z = "Delta^T * V + X"
*
      CALL DVCPY(NMEAS,MEAS(1),MEAS(2),Z(1),Z(2))
      CALL DMMLA(NMEAS,LCSYS,1,ASYS(1,1),ASYS(1,2),ASYS(2,1), ! transposed
     &                         V   (  1),DUMMY    ,V   (  2),
     &                         Z   (  1),DUMMY    ,Z   (  2),DUMMY)
*
*     Compute CHI2
*
*
* Prepare the error matrix reflecting possible statistical corrleations.
* Do this on each iteration to get correct asymetric errors
      CALL PREPARE_ASYM_W(LCSYS,W,X(LCSYS),JERR)
*
      CHI2=0.D0
      DO I=1,NCSYS ! not LCSYS !
        CHI2=CHI2+V(I)**2
      ENDDO
      DO I=1,NMEAS
        DO J=1,NMEAS
          CHI2=CHI2+Z(I)*W(I,J)*Z(J)
        ENDDO
      ENDDO
*
*     Function to be minimized
*
      F=CHI2
*
      IF(IFLAG.NE.3.OR.LUNIT.LE.0) RETURN
*
*     Print PARAMETERS logical lines for possible use in subsequent COMBOS jobs
*
      WRITE(LUNIT,6001)
 6001 FORMAT(/,' PARAMETER logical lines',
     &         ' for possible use in subsequent COMBOS jobs:',/)
 6002 FORMAT(A1,'PARAMETER ',A,3(2X,F10.5,SP),1X,
     &       A1,' CHI2_ASYM_MIN output')
*OS 6004 FORMAT('*          EPLUS,EMINUS,EPARAB,GLOBCC=',SP,4(2X,F10.5))
      DO I=1,LCSYS
        IF(I.LT.LCSYS) THEN 
          CHNAME=CHPARA(I)
        ELSE
          CHNAME=CHMEAS
        ENDIF
 1      CALL MNPOUT(I,CHNAM,VAL,ERR2PB,DUMMY,DUMMY,IVARBL)
        ERR2PB=ERR2PB**2
        CALL MNERRS(I,ERR2P,ERR2N,ERR2PB,DUMMY)
        ERR2P=ERR2P**2
        ERR2N=ERR2N**2
        ERR2PB=ERR2PB**2
        IF (ERR2P.EQ.0.0D0) ERR2P = ERR2PB
        IF (ERR2N.EQ.0.0D0) ERR2N = ERR2PB
        IF(IVARBL.GE.0.AND.CHNAM(:5).NE.'lump ') THEN 
          COMM(1)=' '
        ELSE
          COMM(1)='*'
        ENDIF
        IF(IVARBL.GE.0.AND.CHNAM.EQ.CHNAME(:10)) THEN
          COMM(2)='!'
        ELSE
          COMM(2)='?' ! undefined or illdefined parameter (should not happen)
        ENDIF
        WRITE(LUNIT,6002) COMM(1),CHNAME,VAL,DSQRT(ERR2P),-DSQRT(ERR2N),
     &                    COMM(2)
      ENDDO
      WRITE(LUNIT,'(1X)')

      RETURN
*
      ENTRY FCN_CHI2_ASYM_MIN_INIT(SPECIAL_FLAG,IERR)
*     ==============================================
*
      SPECIAL_FLAG=0
      CALL PREPARE_ACHI2(LCSYS,IERR)
      NOT_READY=IERR.NE.0
      END

************************************************************************
      SUBROUTINE ASYM_CALC_ERRS(LCSYS, X, ASYS)
*     ===============================================
*
*     H. Seywerd
*     April 22, 1997
*
*     Fill the matrix asys with the positive/negative errors
*     as determined by the relative value of x(i) and para(i)
*
      IMPLICIT NONE

      INCLUDE 'master.inc'

      INTEGER LCSYS
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION ASYS(MMEAS,MCSYS)
      INTEGER IPAR, IMEAS
      DOUBLE PRECISION AERROR

      DO IMEAS = 1, NMEAS
        DO IPAR = 1, LCSYS
          IF (IPAR.ne.LCSYS) THEN
C For systematic contributions decide wether to use positive 
C or negative error depending on relative value of fitted parameter 
C and its data value.
            ASYS(IMEAS,IPAR) = AERROR(X(IPAR),PARA(IPAR), 
     $         (CSYSP(IMEAS,IPAR)), -(CSYSN(IMEAS,IPAR)))
C            print *, 'call cs errro', csysp(imeas,ipar),
C     $         -csysn(imeas,ipar),asys(imeas,ipar)
C Paolo.
C            ASYS(IMEAS,IPAR) = 0.5*DERF(X(IPAR)-PARA(IPAR))*(
C     $      (CSYSP(IMEAS,IPAR) + CSYSN(IMEAS,IPAR)) +
C     $         CSYS(IMEAS,IPAR))**2
          ELSE
C For row corresponding to fitted quantity (lcsys) just copy
            ASYS(IMEAS,IPAR) = CSYS(IMEAS,IPAR)
          ENDIF
        ENDDO
      ENDDO

      END


      SUBROUTINE PREPARE_ASYM_W(LCSYS,W,X,IERR)
*
*     Build the stat + uncorr systematic inverse error matrix
*     Use relative value of X, i.e. fitted value to measured values
*     to determine if positive or negative errors should be used
*     If fitted value > measured value use
*     positive error, otherwise negative.
*
      IMPLICIT NONE
      INCLUDE 'master.inc'
*
*     Arguments
*
      INTEGER LCSYS,IERR
      DOUBLE PRECISION W(MMEAS,MMEAS)
      DOUBLE PRECISION X

*
*     Externals
*
      INTEGER LENOCC
*
*     Locals
      INTEGER I,J
      DOUBLE PRECISION STATI, STATJ, USYSI
*
* Determine approriate error
      DOUBLE PRECISION AERROR

      DO I=1,NMEAS
        CSYS(I,LCSYS)=-1.D0
      ENDDO
*
*     Build total error matrices
*
      DO I=1,NMEAS
        STATI = AERROR(MEAS(I),X,DABS(STATP(I)),DABS(STATN(I)))
C        print *, 'call stat errro', STATP(I),STATN(I), stati
        STATI = STAT(I)
        USYSI = AERROR(MEAS(I),X,DABS(USYSP(I)),DABS(USYSN(I)))
C        print *, 'call usys errro', USYSP(I),USYSN(I),usysi
        DO J=1,NMEAS
          STATJ = AERROR(MEAS(J),X,DABS(STATP(J)),DABS(STATN(J)))
          STATJ = STAT(J)
          W(I,J)=STACOR(I,J)*STATI*STATJ ! statistical error matrix
        ENDDO
        W(I,I)=W(I,I)+USYSI**2 ! add uncorrelated systematic error
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

 1000 FORMAT(1X,A,': ',A,' ! IERR=',I6)
*
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


      DOUBLE PRECISION FUNCTION AERROR(X,X0,SIGMA_P,SIGMA_N)
*
* Returns the appropriate value of SIGMA_P or SIGMA_N 
* depending on relation of x to central value x0
*
      IMPLICIT NONE
      DOUBLE PRECISION X, X0, SIGMA_P, SIGMA_N
      DOUBLE PRECISION A, B
      DOUBLE PRECISION SIGMA


      CALL SYMMETRIZE(SIGMA_P,-SIGMA_N,SIGMA)

* Based on that suggested in PDG make smooth function between + and - sigma
      IF (X.GT.X0+SIGMA_P) THEN
        AERROR = SIGMA_P
      ELSEIF (X.LT.X0-SIGMA_N) THEN
        AERROR = SIGMA_N
      ELSE
        A = SIGMA_P + SIGMA_N
C        A = 2.0*SIGMA
        IF (A.EQ.0.D0) THEN
          AERROR = 0.D0
          RETURN
        ENDIF
        A = (SIGMA_P - SIGMA_N)/A
        B = SIGMA_N - A*(X0-SIGMA_N)
        AERROR = A*X + B
      ENDIF
C      PRINT 10, SIGMA_P, SIGMA_N, AERROR, X, X0
 10   FORMAT(6(F10.5,3X))
     
C Paolos error model
C      AERROR = 0.5*DERF(X-X0)*(SIGMA_P-SIGMA_N) + SIGMA
C      print 10, SIGMA_P, SIGMA_N, X, X0, SIGMA, AERROR

      END

