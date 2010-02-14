************************************************************************
*
*     "AVG_AMPLITUDE" termination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.00, February 20, 1997:
*     -- original release
*     Version 2.10, February 27, 1997:
*     -- print WARNING because this code is not tested yet
*
************************************************************************
*
      SUBROUTINE AVG_AMPLITUDE(NSTEP,NC,MSTEP,
     &                         DMS,AMP,ERR2P,ERR2N,ERR2,CL,IERR)
*     ==========================================================
*
*     Termination routine to compute average amplitude over the whole DMS range
*
*     Input:  NSTEP  = number of steps
*     -----   NC     = number of cases
*                      NC=1: means only total error available
*                      NC=2: means also statistical error available
*             MSTEP  = first dimension of arrays DMS, AMP, ... 
*             DMS    = array of DMS values (i.e. step values)
*             AMP    = array of amplitude values
*             ERR2P  = array of positive errors**2 on AMP
*             ERR2N  = array of negative errors**2 on AMP
*             ERR2   = array of symmetric errors**2 on AMP
*             CL     = array of confidence levels of AMP
*             IERR   = array of error codes (0 means OK)
*       
*     Output:  none
*     -------
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER MC,NC,NSTEP,MSTEP
      PARAMETER(MC=2)
      DOUBLE PRECISION DMS(MSTEP),AMP(MSTEP,MC),CL(MSTEP,MC),
     &                 ERR2P(MSTEP,MC),ERR2N(MSTEP,MC),ERR2(MSTEP,MC)
      INTEGER IERR(0:MSTEP,MC)
*
      INCLUDE 'master.inc' ! input data of last combination performed
*
*     Externals
*
      INTEGER LENOCC
      DOUBLE PRECISION DFREQ
      REAL PROB
*
*     Local variables
*
      INTEGER IC,LTOT,LSTA
      PARAMETER(LTOT=1,LSTA=2) ! cases
      CHARACTER*16 TEXT(MC)
      DATA TEXT(LTOT)/'(syst. included)'/
      DATA TEXT(LSTA)/'(stat. only)'/
      INTEGER I,J,IFAIL
      INTEGER MSTEPBIS
      PARAMETER(MSTEPBIS=100)
      DOUBLE PRECISION TAUBS,GAMMA,BW,MSYS,MATRIX(MSTEPBIS,MSTEPBIS)
      LOGICAL XCOSTA,XCOSYS,DEBUG
      DATA XCOSTA/.TRUE./  ! take statistical correlations into account ?
      DATA XCOSYS/.FALSE./  ! take systematic  correlations into account ?
      DATA DEBUG/.FALSE./  ! debug printout for average amplitude ?
      DOUBLE PRECISION SUMW,AVG,CHI2,EAVG,APROB,CHI2NDF,PCHI2
*os   DOUBLE PRECISION          CHI3,                   PCHI3
*
      IF(LUNIT.GT.0) WRITE(LUNIT,2000) CHROUT(:LENOCC(CHROUT))
 2000 FORMAT(/,' AVG_AMPLITUDE called for results obtained with ',A,':')
      IF(LUNIT.GT.0) WRITE(LUNIT,2222)
 2222 FORMAT(/,/,' WARNING: AVG_AMPLITUDE has not been tested',/,
     &           '          and is likely to produce wrong results',/)
*
*     Perform checks
*
      IFAIL=0
      IF(NC.NE.MC) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: invalid number of cases = ',NC
        IFAIL=IFAIL+1
      ENDIF
      IF(NSTEP.LE.1) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: less than 2 steps provided'
        IFAIL=IFAIL+1
      ENDIF
      IF(NSTEP.GT.MSTEP) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: wrong arguments: MSTEP =',MSTEP,
     &   ' less than NSTEP = ',NSTEP
        IFAIL=IFAIL+1
      ENDIF
      IF(NSTEP.GT.MSTEPBIS) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: parameter MSTEPBIS =',MSTEPBIS,
     &   ' less than NSTEP = ',NSTEP
        IFAIL=IFAIL+1
      ENDIF
      IF(CHSTEP.NE.'DMS') THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: step name is not DMS'
        IFAIL=IFAIL+1
      ENDIF
      IF(CHMEAS.NE.'A'.AND.CHMEAS.NE.'AMP'.AND.CHMEAS.NE.'AMPLITUDE')
     &THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &   'AVG_AMPLITUDE: amplitude name is not recognized'
        IFAIL=IFAIL+1
      ENDIF
      IF(IFAIL.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*)
     &                 'AVG_AMPLITUDE: cannot do the job !'
        RETURN
      ENDIF
*
*     Get Bs lifetime
*
      TAUBS=-1.D0
      DO I=1,NPARA
        IF(CHPARA(I).EQ.'TAUBS') TAUBS=PARA(I)
      ENDDO
      IF(TAUBS.LT.1.D0.OR.TAUBS.GT.2.D0) THEN 
        TAUBS=1.5D0
        IF(LUNIT.GT.0) WRITE(LUNIT,3001) TAUBS
 3001   FORMAT(/,' Bs lifetime not found or out of range;',
     &           ' assume TAUBS =',F7.4)
      ENDIF
*      IF(LUNIT.GT.0) WRITE(LUNIT,*) ' TAUBS = ',TAUBS
      GAMMA=1./TAUBS
*
*     Perform a "flat line" fit across the amplitude plot, i.e. determine 
*     the average amplitude (averaged in the whole DMS range)
*
      DO IC=1,NC
        IF(LUNIT.GT.0) WRITE(LUNIT,3000) CHMEAS(:LENOCC(CHMEAS)),
     &    CHSTEP(:LENOCC(CHSTEP)),DMS(1),DMS(NSTEP),
     &    TEXT(IC)(:LENOCC(TEXT(IC)))
 3000   FORMAT(/,' Determine average ',A,' in ',A,
     &           ' range from',F10.4,' to',F10.4,' ',A,':')
*
*       Build statistical error matrix
*
        DO I=1,NSTEP
          DO J=1,NSTEP
            IF(XCOSTA) THEN
              BW=GAMMA**2/(GAMMA**2+(DMS(I)-DMS(J))**2)
            ELSE IF(I.EQ.J) THEN
              BW=1.D0
            ELSE
              BW=0.D0
            ENDIF
            MATRIX(I,J)=DSQRT(ERR2(I,LSTA)*ERR2(J,LSTA))*BW
          ENDDO
        ENDDO
*
*       Build total error matrix (by addind the systematic uncertainties)
*
        DO I=1,NSTEP
          DO J=1,NSTEP
            IF(I.EQ.J.OR.XCOSYS) THEN
              MSYS=DSQRT(ERR2(I,IC)-ERR2(I,LSTA))*
     &             DSQRT(ERR2(J,IC)-ERR2(J,LSTA))
            ELSE
              MSYS=0.D0
            ENDIF
            MATRIX(I,J)=MATRIX(I,J)+MSYS
          ENDDO
        ENDDO
*
        IF(XCOSTA.AND.LUNIT.GT.0.AND.DEBUG) THEN 
          WRITE(LUNIT,4000) 
 4000     FORMAT(/,' Error matrix:')
          DO I=1,NSTEP
            WRITE(LUNIT,4001) (MATRIX(I,J),J=1,NSTEP)
 4001       FORMAT(30F6.2)
          ENDDO
        ENDIF
*
*       Invert error matrix
*
        CALL DSINV(NSTEP,MATRIX,MSTEPBIS,IFAIL)
        IF(IFAIL.NE.0) THEN  
          IF(LUNIT.GT.0) WRITE(LUNIT,4002) IFAIL
 4002     FORMAT(' DSINV cannot invert error matrix, IFAIL=',I8)
          RETURN
        ENDIF
        IF(XCOSTA.AND.LUNIT.GT.0.AND.DEBUG) THEN 
          WRITE(LUNIT,4003) 
 4003     FORMAT(/,' Inverse of error matrix (=weights):')
          DO I=1,NSTEP
            WRITE(*,4001) (MATRIX(I,J),J=1,NSTEP)
          ENDDO
        ENDIF
*
*       Compute weighted average and chi**2
*
        SUMW=0.D0
        AVG=0.D0
        CHI2=0.D0
        DO I=1,NSTEP
          DO J=1,NSTEP
            SUMW=SUMW+MATRIX(I,J)
            AVG=AVG+MATRIX(I,J)*0.5D0*(AMP(I,IC)+AMP(J,IC))
            CHI2=CHI2+AMP(I,IC)*MATRIX(I,J)*AMP(J,IC)
          ENDDO
        ENDDO
        IF(SUMW.EQ.0.D0) THEN 
          IF(LUNIT.GT.0) WRITE(LUNIT,*) 'ERROR: SUMW = 0'
          AVG=9999.D0
          EAVG=9999.D0
        ELSE
          AVG=AVG/SUMW
          EAVG=1.D0/DSIGN(DSQRT(DABS(SUMW)),SUMW)
        ENDIF
*
*OS     CHI3=0.D0
*OS     DO I=1,NSTEP
*OS       DO J=1,NSTEP
*OS         CHI3=CHI3+(AMP(I,IC)-AVG)*MATRIX(I,J)*(AMP(J,IC)-AVG)
*OS       ENDDO
*OS     ENDDO
*OS     PCHI3=DBLE(PROB(SNGL(CHI3),NSTEP-1))
*
        APROB=2.D0*(1.D0-DFREQ(DABS(AVG/EAVG)))
        CHI2NDF=CHI2/DBLE(NSTEP)
        PCHI2=DBLE(PROB(SNGL(CHI2),NSTEP))
        IF(LUNIT.GT.0) THEN 
          IF(DEBUG) WRITE(LUNIT,5000)
     &    SUMW,AVG,EAVG,AVG/EAVG,APROB,CHI2,NSTEP,CHI2NDF,PCHI2
 5000     FORMAT(/,'Sum of weights = ',G12.6,/,
     &         /,'Weighted average of amplitudes = ',G12.6,' +- ',G12.6,
     &         /,'  ---> Number of sigmas from zero = ',G12.6,
     &         /,'       corresponding to a probability of ',G12.6,/,
     &         /,'Chi**2 / NDF = ',G12.6,' / ',I3,' = ',G12.6,
     &         /,'  ---> Chi**2 probability = ',G12.6,/)
          WRITE(LUNIT,5001) CHMEAS(:LENOCC(CHMEAS)),
     &                      AVG,EAVG,TEXT(IC)(:LENOCC(TEXT(IC))),
     &                      CHI2,NSTEP,CHI2NDF,PCHI2
 5001     FORMAT(3X,'average ',A,' =',F10.4,' +- ',F10.4,1X,A,/,
     &           3X,'chi2/NDF = chi2/step =',F10.4,' /',I3,' =',F10.4,
     &           3X,'prob(chi2) =',F6.4)
*OS       PRINT * ,'Chi3 = ',CHI3,' NDF = ',NSTEP-1,
*OS  &             ' prob(chi3) = ',PCHI3
        ENDIF
      ENDDO
      END
