************************************************************************
*
*     Termination routine of the COMBOS program
*     to compute the final b-hadron fractions
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.52, July 12, 1999:
*     -- original release
*     Version 2.53, November 27, 1999:
*     -- print also a subset of the full error matrix, 
*        print also the Bs fraction from mixing only (including its
*        correlations with the Bs fraction and the b-baryon fraction 
*        obtained without mixing inofrmation) in a special format to 
*        be read by an external program to compute the final set of 
*        fractions
*
************************************************************************
*
      SUBROUTINE FINAL_FRACTIONS(NS,NC,MS,
     &                           STEP,CVAL,ERR2P,ERR2N,ERR2,CL,IERR)
*     ==============================================================
*
*     Input:  NS     = number of steps
*     -----   NC     = number of cases
*                      NC=1: means only total error available
*                      NC=2: means also statistical error available
*             MS     = first dimension of arrays STEP, CVAL, ... 
*             STEP   = array of step values
*             CVAL   = array of combined values of the measured quantity
*             ERR2P  = array of positive errors**2 on CVAL
*             ERR2N  = array of negative errors**2 on CVAL
*             ERR2   = array of symmetric errors**2 on CVAL
*             CL     = array of confidence levels of CVAL
*             IERR   = array of error codes (0 means OK)
*
*     For step IS=1,NS:
*             CVAL (IS,1) = combined value at step IS (including syst.)
*             ERR2P(IS,1) = positive error**2 on CVAL(IS,1)
*             ERR2N(IS,1) = negative error**2 on CVAL(IS,1)
*             ERR2 (IS,1) = symmetric error**2 on CVAL(IS,1)
*             CL   (IS,1) = confidence level of fit performed to get CVAL(IS,1)
*                           (or -1 if not available)
*
*     Other input (available from common block described in master.inc):
*     -----------
*              
*             CHROUT = name of combination routine 
*                       used to obtain CVAL,ERR2P,ERR2N,ERR2,CL and IERR
*             CHCOMB = name of combined analysis
*             CHMEAS = name of measured quantity
*             CHSTEP = name of step variable
*             LUNIT  = logical unit for printout
*                       produced by this termination routine
*
*     Output:  none
*     -------
*
*     Note: this routine only write output to logical unit LUNIT; 
*     ----- if LUNIT is zero or less, then this routine should not produce
*           any output
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER MC,NC,NS,MS
      PARAMETER(MC=2)
      DOUBLE PRECISION STEP(MS),CVAL(MS,MC),CL(MS,MC),
     &                 ERR2P(MS,MC),ERR2N(MS,MC),ERR2(MS,MC)
      INTEGER IERR(0:MS,MC)
*
      INCLUDE 'master.inc'
*
      INTEGER KERR
*
      KERR=0
      IF(NS.NE.1) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'FINAL_FRACTIONS: NS = ',NS,': Nothing done.'
        KERR=KERR+1
      ENDIF
      IF(NC.LE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'FINAL_FRACTIONS: NC = ',NC
        KERR=KERR+1
      ENDIF
      IF(IERR(1,1).NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'FINAL_FRACTIONS: IERR = ',IERR(1,1)
        KERR=KERR+1
      ENDIF
      IF(KERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'FINAL_FRACTIONS: Nothing done (1): KERR = ',KERR
        RETURN
      ENDIF
      CALL STEP3(CVAL(1,1),DSQRT(ERR2(1,1)),KERR)
      IF(KERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'FINAL_FRACTIONS: Nothing done (2): KERR = ',KERR
        RETURN
      ENDIF
      END
*
************************************************************************
*
      SUBROUTINE STEP3(DMD,EDMD)
*     ==========================
*
*     Olivier Schneider, CERN/PPE-ALE
*
*
*     Input:  DMD  = world average of dmd
*     -----   EDMD = error on dmd
*
*     Output: KERR = error flag (0 means OK)
*     ------ 
*
*     Arguments
*
      DOUBLE PRECISION DMD,EDMD
*
*     Common blocks
*
      INCLUDE 'master.inc'
      INCLUDE 'circularity_g.inc'
* 
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER NNEWPAR2
      PARAMETER(NNEWPAR2=NNEWPAR*NNEWPAR)
      INTEGER INDX(NNEWPAR)
      DOUBLE PRECISION XIN (NNEWPAR),EXIN (NNEWPAR,NNEWPAR)
      DOUBLE PRECISION ERR (NNEWPAR)
      DOUBLE PRECISION XOUT(NOLDPAR),EXOUT(NOLDPAR,NOLDPAR)
      DOUBLE PRECISION RHO,CORR(NOLDPAR,NNEWPAR)
*
      INTEGER KERR
      INTEGER I,J,K
      LOGICAL SHOW
      INTEGER K1,K2,J1,J2
      INTEGER NTOTPAR,MTOTPAR
      PARAMETER(MTOTPAR=NNEWPAR+NOLDPAR)
      REAL ETOT(MTOTpAR,MTOTPAR)
      INTEGER INDTOT(MTOTPAR)
      CHARACTER*10 CHERR(0:MTOTPAR),CHERN(0:MTOTPAR)
*
*     Input
*
      KERR=0
      DO I=1,NNEWPAR
        IF(CHNEWPAR(I).EQ.'DMD') THEN
          XIN(I)=DMD
          ERR(I)=EDMD
          INDX(I)=9999
        ELSE
          INDX(I)=0
          DO J=1,NPARA
            IF(CHPARA(J).EQ.CHNEWPAR(I)) THEN 
              XIN(I)=PARA(J)
              ERR(I)=EXCU(J)
              INDX(I)=J
            ENDIF
          ENDDO
          IF(INDX(I).EQ.0) THEN ! missing parameter
            XIN(I)=0.D0
            ERR(I)=0.D0
            IF(NEEDED(I)) THEN
              KERR=KERR+1
              IF(LUNIT.GT.0) WRITE(LUNIT,*) 'FINAL_FRACTIONS: ',
     &         'parameter '//CHNEWPAR(I)//' not found'
            ELSE IF(I.EQ.IR.OR.I.EQ.IRS.OR.I.EQ.IRL) THEN
              XIN(I)=1.D0
            ELSE IF(I.EQ.ICHIS) THEN
              XIN(I)=0.5D0
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF(LUNIT.LE.0) KERR=KERR+1
      IF(KERR.NE.0) RETURN
      RHO=0.D0 ! correlation between IFLAMB and IFSBR
      DO J=1,NPARA
        IF(CHPARA(J).EQ.'CORR_FLAMB_FBS_') RHO=PARA(J)
      ENDDO
* 
      DO I=1,NNEWPAR
        DO J=1,NNEWPAR
          IF(I.EQ.J) THEN
            IF(ERR(I).EQ.0.D0) THEN
              INDX(I)=0
            ELSE
              INDX(I)=1
            ENDIF
            EXIN(I,I)=ERR(I)**2
          ELSE
            EXIN(I,J)=0.D0
          ENDIF
        ENDDO
      ENDDO
      EXIN(IFLAMB,IFSBR)=RHO*ERR(IFLAMB)*ERR(IFSBR)
      EXIN(IFSBR,IFLAMB)=EXIN(IFLAMB,IFSBR)
*
*     Compute XOUT from XIN
*
      CALL CIRCULARITY_G(INDX,XIN,EXIN,XOUT,EXOUT,CORR)
*
*     Print final results
*
 3000 FORMAT(/,1X,A,
     &       /,1X,'------',/,
     &       /,1X,'       Parameter_name',T26,'  Value+-error')
 3001 FORMAT(1X,6X,1X,A,T26,F10.5,' +- ',F10.5)
 3011 FORMAT(1X,6X,1X,A,T26,F10.5)
*
      WRITE(LUNIT,3000) 'Input'
      NTOTPAR=0
      CHERN(NTOTPAR)=' '
      DO I=1,NNEWPAR
        IF(INDX(I).GT.0) THEN
          NTOTPAR=NTOTPAR+1
          INDTOT(NTOTPAR)=I
          CHERN(NTOTPAR)=CHNEWPAR(I)
          WRITE(LUNIT,3001) CHNEWPAR(I),XIN(I),SQRT(EXIN(I,I))
        ENDIF
      ENDDO
      WRITE(LUNIT,3011) 'CORR_FLAMB_FBS_',RHO
*
      WRITE(LUNIT,3000) 'Output'
      DO I=1,NOLDPAR
        SHOW=.TRUE.
*zzz        IF(INDX(ICHID4S).LE.0.AND.
*zzz     &     (I.EQ.ICHID1.OR.I.EQ.IXDW.OR.I.EQ.IDMDW)) SHOW=.FALSE.
*zzz        IF(INDX(IFSBR).LE.0.AND.I.EQ.IFS1) SHOW=.FALSE.
        IF(SHOW) THEN 
          NTOTPAR=NTOTPAR+1
          INDTOT(NTOTPAR)=-I
          CHERN(NTOTPAR)=CHOLDPAR(I)
          WRITE(LUNIT,3001) CHOLDPAR(I),XOUT(I),SQRT(EXOUT(I,I))
        ENDIF
      ENDDO
*
 3002 FORMAT(/,1X,'Full error/correlation matrix on the above ',I2,
     &           ' parameters:',
     &   /,1X,'(elements on or below diagonal are from error matrix,',
     &   /,1X,' elements above diagonal are from correlation matrix)')
 3003 FORMAT(/,1X,'Error/correlation matrix on the ',I2,
     &           ' (final and intermediate) b-hadron fractions:',
     &   /,1X,'(elements on or below diagonal are from error matrix,',
     &   /,1X,' elements above diagonal are from correlation matrix)')
 3004 FORMAT(1X,16A)
*
*
      DO K=1,2
        IF(K.EQ.1) THEN
          WRITE(LUNIT,3002) NTOTPAR
        ELSE
          CALL UZERO(INDTOT,1,MTOTPAR)
          NTOTPAR=5
          INDTOT(1)=+IFLAMB
          INDTOT(2)=+IFSBR
          INDTOT(3)=-IFS1
          INDTOT(4)=-IFS
          INDTOT(5)=-IFD
          DO K1=1,NTOTPAR
            IF(INDTOT(K1).GT.0) THEN
              CHERN(K1)=CHNEWPAR(INDTOT(K1))
            ELSE
              CHERN(K1)=CHOLDPAR(-INDTOT(K1))
            ENDIF
          ENDDO
          WRITE(LUNIT,3003) NTOTPAR
        ENDIF
*
*     Build error matrix
*       K=1: full error matrix
*       K=2: just error matrix for final fractions
*
      DO K1=1,NTOTPAR
        J1=INDTOT(K1)
        DO K2=1,NTOTPAR
          J2=INDTOT(K2)
          IF(J1.GT.0.AND.J2.GT.0) THEN
            ETOT(K1,K2)=SNGL(EXIN(J1,J2))
          ELSE IF(J1.LT.0.AND.J2.LT.0) THEN
            ETOT(K1,K2)=SNGL(EXOUT(-J1,-J2))
          ELSE
            ETOT(K1,K2)=SNGL(CORR(-MIN0(J1,J2),MAX0(J1,J2)))
          ENDIF
        ENDDO
      ENDDO
*
*     Print error and correlation matrices
*
      DO K2=0,NTOTPAR
        I=LENOCC(CHERN(K2))
        J=LEN   (CHERN(K2))
        IF(I.EQ.J) THEN 
          CHERN(K2)(J-3:)='... '
          I=1
        ELSE
          I=(J-I+1)/2
        ENDIF
        CHERR(K2)=CHERN(0)(:I)//CHERN(K2)
      ENDDO
      WRITE(LUNIT,3004) (CHERR(K2),K2=0,NTOTPAR)
      DO K1=1,NTOTPAR
        CHERR(0)=CHERN(K1)
        DO K2=1,NTOTPAR
          IF(K2.LE.K1) THEN ! write element of error matrix
            WRITE(CHERR(K2),'(E10.3)') ETOT(K1,K2)
          ELSE              ! write element of correlation matrix
            WRITE(CHERR(K2),'(SP,F9.4)') ETOT(K1,K2)/(SQRT(ETOT(K1,K1))*
     &                                                SQRT(ETOT(K2,K2)))
          ENDIF
        ENDDO
        WRITE(LUNIT,3004) (CHERR(K2),K2=0,NTOTPAR)
      ENDDO
*
      ENDDO
*
*     Final fractions
*
 4000 FORMAT(/
     & '* COMBOS cards for (final and intermediate) b-hadron fractions')
 4001 FORMAT(' PARAMETER ',A,T30,F8.4,SP,2F8.4,'  ! FINAL_FRACTIONS')
      WRITE(LUNIT,4000)
      DO I=1,NTOTPAR
        J=INDTOT(I)
        IF(J.GT.0) THEN
          WRITE(LUNIT,4001) CHERN(I),XIN(J),
     &     SQRT(EXIN(J,J)),-SQRT(EXIN(J,J))
        ELSE IF(J.LT.0) THEN
          WRITE(LUNIT,4001) CHERN(I),XOUT(-J),
     &     SQRT(EXOUT(-J,-J)),-SQRT(EXOUT(-J,-J))
        ENDIF
      ENDDO
*
*     Write info in special format for the program to calculate fractions
*     (note that all fractions are here printed in percents) 
*
      WRITE(LUNIT,*) ' '
      WRITE(LUNIT,9000) 'FBS_MIX',100.*XOUT(IFS1)
      WRITE(LUNIT,9000) 'ERR_FBS_MIX',100.*SQRT(EXOUT(IFS1,IFS1))
      WRITE(LUNIT,9001) 'COV_FBS_ST1_FBS_MIX',10000.*CORR(IFS1,IFSBR)
*    & ,CORR(IFS1,IFSBR)/
*    &  (SQRT(EXOUT(IFS1,IFS1))*SQRT(EXIN(IFSBR,IFSBR)))
      WRITE(LUNIT,9001) 'COV_FBB_ST1_FBS_MIX',10000.*CORR(IFS1,IFLAMB)
*    & ,CORR(IFS1,IFLAMB)/
*    &  (SQRT(EXOUT(IFS1,IFS1))*SQRT(EXIN(IFLAMB,IFLAMB)))
 9000 FORMAT(A,':',/,F12.8)
 9001 FORMAT(A,':',/,SP,2F12.8)
      END
