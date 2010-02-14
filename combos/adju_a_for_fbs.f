************************************************************************
*
*     Preparation routine of the COMBOS program
*
*     Adjust central values and statistical uncertainties
*     (due to the Bs fraction) of the measured amplitude to new value of FBS;
*     adjustment performed depends on the value of the parameter ADJU_A_FOR_FBS 
*
*     Olivier Schneider, CERN/PPE-ALE
*
*
*     Version 2.52, July 12, 1999:
*     -- original release as routine ADJUST_SIGASTAT
*     -- statistical amplitude rescaled assuming it is proportional to 1/FBS
*     Version 2.53, November 27, 1999:
*     -- routine renamed ADJU_A_FOR_FBS
*     -- both the central value and the statistical uncertainty can be 
*        adjusted for the common value of FBS according to an algorithm given 
*        by the value of the parameter ADJU_A_FOR_FBS.
*     Version 2.95, March 31, 2000
*     -- bug fixed (variable I was used at initialisation stage without 
*        being assigned a value)
*     Version 2.96, April 3, 2000
*     -- bug fixed in scaling of negative statistical uncertainty
*
************************************************************************
*
      SUBROUTINE ADJU_A_FOR_FBS(N,ICOMB,AMEAS)
*     ========================================
*
*     Input:  N      = number of analyses to combine
*     -----   ICOMB  = array of indexes to analyses
*                      ICOMB(I) = analysis number of Ith analysis to combine
*             AMEAS  = ajustments already performed on central values
*                      by routine MASTER
*                      AMEAS(I,ICONT) = adjustment due to systematic
*                                       contribution ICONT performed on
*                                       Ith analysis to combine
*
*     Output: none
*     ------
*
*     Input data:  available from combos.inc and master.inc
*     -----------
*
*     Output data: written in master.inc
*     ------------
*
*     Note: this routine only writes output to logical unit LUNIT; 
*     ----- if LUNIT is zero or less, then this routine should not produce
*           any output
*
*****************************
*
*     For each analysis, the adjustements are performed depending on the
*     value of the parameter ADJU_A_FOR_FBS.
*
*     If 'xy' is the value of the parameter ADJU_A_FOR_FBS, then
*     x = INT(ABS(xy)/10) is a code for the adjustement of the central value and
*     y = MOD(ABS(xy),10) is a code for the adjustement of the stat. uncertainty
*     sign(xy) determines wether the adjustments for FBS are performed before
*     (negative sign) or after (positive sign) all the other adjustments 
*
*     Meaning for codes x and/or y:
*     code=0 means default adjustement (defined by value of parameter 
*            ADJU_A_FOR_FBS given for the combined analysis, or code=2 if
*            no valid value given in the  combined analysis)
*     code=1 means no adjustment
*     code=2 original adjustment performed by COMBOS
*     code=3 adjustment assuming 1/FBS scaling
*     code=4 adjustment assuming A/sigma_A_strat = constant (indep. of FBS)
*
*     Note that the value 44 is not permitted for parameter ADJU_A_FOR_FBS
*     (in that case the default will be used).
*
*****************************
*
      IMPLICIT NONE
      INCLUDE 'combos.inc'
      INCLUDE 'master.inc'
*
*     Argument
*
      INTEGER N
      INTEGER ICOMB(N)
      DOUBLE PRECISION AMEAS(MMEAS,MCONT)
*
*     Local variable
*
      CHARACTER*6 CH6
      CHARACTER*8 APPLIED(MMEAS)
      INTEGER I,J,JJ,K,IANAL,ICONT,IS
      CHARACTER*(*) CHFBS
      PARAMETER(CHFBS='FBS')
      DOUBLE PRECISION FACT(MMEAS),ATOT,AFBS,OMEAS,OSTAT,FF
      CHARACTER*10 CHFACT
      LOGICAL BEFORE(0:MMEAS)
*
      INTEGER CODE(2,0:MMEAS) ! index 0 for combined analysis
      INTEGER MCODE,DEFCODE
      PARAMETER(MCODE=4) ! maximum value of code
      PARAMETER(DEFCODE=2) ! default code
      CHARACTER*60 CHCODE(MCODE) ! meaning of code
      DATA CHCODE/
     1 'no adjustment (not even regular combos adjustment)',
     2 'regular combos adjustment according to syst. uncertainty',
     3 'adjustment assuming 1/FBS scaling',
     4 'adjustment assuming A/sigma_A_stat independent of FBS'/
*
      SAVE FACT,BEFORE,CODE ! important
*
      IF(ISTEP.EQ.1) THEN ! initialize
        IF(N.NE.NMEAS) CALL COMBOS_ERROR(-1,
     &                 'Inconsistency in ADJU_A_FOR_FBS',
     &                 'please check input arguments')
        IF(LUNIT.GT.0) THEN ! write what the routine is doing  
          WRITE(LUNIT,1000) (CHFBS,JJ=1,3),(JJ,CHCODE(JJ),JJ=1,MCODE)
 1000     FORMAT(/,'Preparation routine ADJU_A_FOR_FBS will adjust',
     &           /,'central values and statistical uncertainties on',
     &           /,'the measured amplitudes for the new value of ',A,
     &     '.',/,/,'Meaning of the codes:',
     &           /,'  combined code = sign*(x*10+y)',
     &           /,'  sign=-1: adjustments for ',A,' applied before',
     &                       ' all other adjustements',
     &           /,'  sign=+1: adjustments for ',A,' applied after ',
     &                       ' all other adjustements',
     &           /,'  1st digit = x = code for central value',
     &                              ' (see meaning below)',
     &           /,'  2nd digit = y = code for statistical uncertainty',
     &                              ' (see meaning below)',
     &         /,/,('  code=',I1,': ',A))
          WRITE(LUNIT,1200) CHFBS
 1200     FORMAT(/,'Analysis',T50,A,T60,'Adj. code',T70,'    Factor')
        ENDIF
        K=0
        CH6='N/A'
        CODE(1,0)=DEFCODE
        CODE(2,0)=DEFCODE
        BEFORE(0)=.TRUE.
        DO J=1,NPARA
          IF(CHPARA(J).EQ.CHFBS) THEN
            K=J
            WRITE(CH6,'(F6.4)') PARA(J)
          ELSE IF(CHPARA(J).EQ.CHROUT) THEN
            BEFORE(0)=PARA(J).LE.0.D0
            CODE(2,0)=MIN0(NINT(ABS(PARA(J))),99)
            CODE(1,0)=CODE(2,0)/10      ! default code for stat. uncertainty
            CODE(2,0)=MOD(CODE(2,0),10) ! default code for central value
            DO JJ=1,2
              IF(CODE(JJ,0).LT.1.OR.CODE(JJ,0).GT.MCODE.OR.
     &           (CODE(1,0).EQ.4.AND.CODE(2,0).EQ.4))
     &         CODE(JJ,0)=DEFCODE ! ultimate default
            ENDDO
          ENDIF
        ENDDO
        IF(LUNIT.GT.0) WRITE(LUNIT,1001) CHCOMB,CH6
 1001   FORMAT(A50,T50,A,T60,A,T70,A)
        DO I=1,N
          IANAL=ICOMB(I)
          CH6='N/A'
          APPLIED(I)='NO'
          CODE(1,I)=DEFCODE
          CODE(2,I)=DEFCODE
          BEFORE(I)=.TRUE.
          FACT(I)=1.D0
          DO J=1,NPAR(IANAL)
            IF(CHNAM(KPAR(J,IANAL)).EQ.CHFBS) THEN
              WRITE(CH6,'(F6.4)') PAR(0,J,IANAL)
              IF(K.NE.0) FACT(I)=PAR(0,J,IANAL)/PARA(K)
            ELSE IF(CHNAM(KPAR(J,IANAL)).EQ.CHROUT) THEN
              BEFORE(I)=PAR(0,J,IANAL).LE.0.D0
              CODE(2,I)=MIN0(NINT(ABS(PAR(0,J,IANAL))),99)
              CODE(1,I)=CODE(2,I)/10      ! code for statistical uncertainty
              CODE(2,I)=MOD(CODE(2,I),10) ! code for central value
              DO JJ=1,2
                IF(CODE(JJ,I).LT.1.OR.CODE(JJ,I).GT.MCODE.OR.
     &             (CODE(1,I).EQ.4.AND.CODE(2,I).EQ.4)) THEN
                  CODE(JJ,I)=CODE(JJ,0) ! use default code
                  IF(JJ.EQ.1) BEFORE(I)=BEFORE(0)
                ENDIF
              ENDDO
              IF(BEFORE(I)) THEN
                IS=-1
              ELSE
                IS=+1
              ENDIF
              WRITE(APPLIED(I),'(''YES  '',SP,I3)')
     &         IS*(CODE(1,I)*10+CODE(2,I))
            ENDIF
          ENDDO
          IF(CH6.EQ.'N/A') APPLIED(I)='NO'
          IF(APPLIED(I).EQ.'NO') THEN 
            CODE(1,I)=DEFCODE
            CODE(2,I)=DEFCODE
            BEFORE(I)=.TRUE.
            FACT(I)=1.D0
          ENDIF
          IF(LUNIT.GT.0) THEN
            IF(CODE(1,I).EQ.3.OR.CODE(2,I).EQ.3) THEN
              WRITE(CHFACT,'(F10.4)') FACT(I)
            ELSE
              CHFACT=' '
            ENDIF
            WRITE(LUNIT,1001) CHANAL(I),CH6,APPLIED(I),CHFACT
          ENDIF
        ENDDO
      ENDIF
*
*     Now do the job ...
*
      DO I=1,NMEAS
        IANAL=ICOMB(I)
*
*       Find out adjustments already done by COMBOS on central value
*
        AFBS=0.D0
        ATOT=0.D0
        DO J=1,NPAR(IANAL)
          ICONT=PCONT(J,IANAL)
          IF(CHNAM(KPAR(J,IANAL)).EQ.CHFBS) THEN
            AFBS=AMEAS(I,ICONT)
          ELSE
            ATOT=ATOT+AMEAS(I,ICONT)
          ENDIF
        ENDDO
*
*       Undo regular COMBOS adjustments 
*
        MEAS(I)=MEAS(I)-AFBS                ! undo FBS adjustment
        IF(BEFORE(I)) MEAS(I)=MEAS(I)-ATOT  ! undo other adjustements
*       ...                                 ! nothing to undo on stat. error
*
*       Save original measurement
*
        OMEAS=MEAS(I)
        OSTAT=STAT(I)
*
*       Adjust central value (except if A/sigmas_A_stat = constant)
*
        IF(CODE(1,I).EQ.DEFCODE) THEN ! apply regular COMBOS adjustment
          MEAS(I)=MEAS(I)+AFBS
        ELSE IF(CODE(1,I).EQ.3) THEN ! scale like 1/FBS
          MEAS(I)=MEAS(I)*FACT(I)
        ENDIF ! do nothing if CODE(1,I)=1 or 4
*
*       Adjust statistical uncertainties
*
        IF(CODE(2,I).EQ.3) THEN ! scale like 1/FBS
          FF=FACT(I)
        ELSE IF(CODE(2,I).EQ.4) THEN ! assume A/sigma_A_stat = constant
          FF=MEAS(I)/OMEAS ! factor applied on original measurement
        ELSE ! CODE(2,I)=1 or 2
          FF=1.D0
        ENDIF 
        STATP(I)=STATP(I)*FF
        STATN(I)=STATN(I)*FF
        STAT (I)=STAT (I)*FF
*
*       Adjust central value (if A/sigma_A_stat = constant)
*
        IF(CODE(1,I).EQ.4) THEN
          FF=STAT(I)/OSTAT
          MEAS(I)=MEAS(I)*FF
        ENDIF
*
*       Restore all other adjustments if needed
*
        IF(BEFORE(I)) MEAS(I)=MEAS(I)+ATOT
*
      ENDDO
      END
