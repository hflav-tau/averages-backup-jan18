************************************************************************
*
*     "CHI2_SYM_RVK" combination routine for COMBOS program
*
*     Bob Kowalewski, U. of Victoria
*
************************************************************************
*
      SUBROUTINE CHI2_SYM_RVK(CVAL,ERR2P,ERR2N,CL,ES,EU,EC,IERR)
*     =============================================
*
*     Olivier Schneider, CERN/PPE-ALE
*     December 10, 1996
*
*     Combination routine using chi2:
*       - only symmetric uncertainties are used
*       - can handle correlated statistical uncertainties with any 
*          correlation coefficient
*       - can handle 100% or 0% systematic correlations
*
*     Input data: available from a common block in file master.inc
*     ----------
*
      INCLUDE 'master.inc'
      INCLUDE 'combos.inc'

*      
*     Output data: returned in the arguments 
*     -----------
* 
      DOUBLE PRECISION CVAL  ! combined value (or average value)
      DOUBLE PRECISION ERR2P ! positive total uncertainty**2 on CVAL
      DOUBLE PRECISION ERR2N ! negative total uncertainty**2 on CVAL
      DOUBLE PRECISION CL    ! confidence level of combination (e.g. chi2 prob.)
      INTEGER IERR           ! error flag (non-zero if combination failed)
*
********
*
*     Externals
*
      REAL PROB
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I,J,K,ITER,MXITER
      PARAMETER(MXITER=50)
      DOUBLE PRECISION W(MMEAS,MMEAS),WU(MMEAS,MMEAS),WS(MMEAS,MMEAS),
     +   WC(MMEAS,MMEAS,MCSYS),W0,WY,W1,TEMP
      DOUBLE PRECISION V(MCSYS),
     &                 CHI2,DUMMY(MMEAS),
     &                 T1(MMEAS,MMEAS),
     &                 T2(MMEAS,MMEAS)
      DOUBLE PRECISION ES,EU,EC(MCSYS)
      DOUBLE PRECISION SCALF(MMEAS,-1:NCSYS),CHI2L

*
*     Default output arguments
*
      ITER = 0
      CVAL = 0.D0
      TEMP = 0.D0
      DO I=1,NMEAS
         CVAL = CVAL + MEAS(I)/STAT(I)**2
         TEMP = TEMP + 1.0D0/STAT(I)**2
      ENDDO
      CVAL = CVAL/TEMP
      CL   = -1.D0
      IERR = 9999
      CHI2 = 0.D0
*
* Master loop - 1st combination is using absolute errors, subsequent
* combinations using specified dependence.
*
      DO ITER=0,MXITER
         ERR2P = 0.D0
         ERR2N = 0.D0
         CHI2L=CHI2
         DO I=1,NMEAS
            IF(ITER.eq.0) THEN
               DO J=-1,NCSYS
                  SCALF(I,J) = 1.D0
               ENDDO
            ELSE
               DO J=1,NCSYS
                  IF (XON(LSYRE)) THEN
                     SCALF(I,J) = CVAL/MEAS(I)
                  ELSEIF (XON(LSYAB)) THEN
                     SCALF(I,J) = 1.D0
                  ELSE
                     SCALF(I,J) = (CVAL/MEAS(I))**XSCSYS(I,J)
                  ENDIF
               ENDDO
               IF (XON(LSYRE)) THEN
                  SCALF(I,0) = CVAL/MEAS(I)
               ELSEIF (XON(LSYAB)) THEN
                  SCALF(I,0) = 1.D0
               ELSE
                  SCALF(I,0) = (CVAL/MEAS(I))**XSUSYS(I)
               ENDIF
               IF (XON(LSTSQ)) THEN
                  SCALF(I,-1) = CVAL/MEAS(I)**0.5
               ELSEIF (XON(LSTAB)) THEN
                  SCALF(I,-1) = 1.D0
               ELSE
                  SCALF(I,-1) = (CVAL/MEAS(I))**XSSTAT(I)
               ENDIF
            ENDIF
         ENDDO
CRVK         PRINT *, 'in chi2_sym_rvk, scalf',(scalf(1,j),j=-1,ncsys)
CRVK         PRINT *, 'in chi2_sym_rvk, scalf',(scalf(2,j),j=-1,ncsys)
         CALL PREPARE_CHI2_RVK(CVAL,W,WU,WS,WC,SCALF,IERR)
         IF(IERR.NE.0) RETURN
         W0=0
         WY=0
         DO I=1,NMEAS
            DO J=1,NMEAS
               W0 = W0 + W(I,J)
               WY = WY + MEAS(I)*W(I,J)
            ENDDO
         ENDDO
*     average
         CVAL = WY/W0
*     errors broken down per source - start with statistical
         CALL DMMLT(NMEAS,NMEAS,NMEAS,W(1,1),W(2,1),W(1,2)
     &        ,WS(1,1),WS(2,1),WS(1,2)
     &        ,T1(1,1),T1(2,1),T1(1,2),DUMMY)
         CALL DMMLT(NMEAS,NMEAS,NMEAS,T1(1,1),T1(2,1),T1(1,2)
     &        ,W(1,1),W(2,1),W(1,2)
     &        ,T2(1,1),T2(2,1),T2(1,2),W1)
         W1 = 0
         DO I=1,NMEAS
            DO J=1,NMEAS
               W1 = W1 + T2(I,J)
            ENDDO
         ENDDO
         ES = SQRT(W1)/W0
         ERR2P = ERR2P + W1/W0**2
*     uncorrelated systematic
         CALL DMMLT(NMEAS,NMEAS,NMEAS,W(1,1),W(2,1),W(1,2)
     &        ,WU(1,1),WU(2,1),WU(1,2)
     &        ,T1(1,1),T1(2,1),T1(1,2),DUMMY)
         CALL DMMLT(NMEAS,NMEAS,NMEAS,T1(1,1),T1(2,1),T1(1,2)
     &        ,W(1,1),W(2,1),W(1,2)
     &        ,T2(1,1),T2(2,1),T2(1,2),W1)
         W1=0
         DO I=1,NMEAS
            DO J=1,NMEAS
               W1 = W1 + T2(I,J)
            ENDDO
         ENDDO
         EU = SQRT(W1)/W0
         ERR2P = ERR2P + W1/W0**2
*     correlated systematics
         DO K=1,NCSYS
            CALL DMMLT(NMEAS,NMEAS,NMEAS,W(1,1),W(2,1),W(1,2)
     &           ,WC(1,1,K),WC(2,1,K),WC(1,2,K)
     &           ,T1(1,1),T1(2,1),T1(1,2),DUMMY)
            CALL DMMLT(NMEAS,NMEAS,NMEAS,T1(1,1),T1(2,1),T1(1,2)
     &           ,W(1,1),W(2,1),W(1,2)
     &           ,T2(1,1),T2(2,1),T2(1,2),W1)
            W1=0
            DO I=1,NMEAS
               DO J=1,NMEAS
                  W1 = W1 + T2(I,J)
               ENDDO
            ENDDO
            EC(K) = SQRT(W1)/W0
            ERR2P = ERR2P + W1/W0**2
         ENDDO
*     overall error
***         ERR2P = SQRT(ERR2P)
*     
*     Compute CHI2
*     
         CHI2=0.D0
         DO I=1,NMEAS
            DO J=1,NMEAS
               CHI2=CHI2+(MEAS(I)-CVAL)*W(I,J)*(MEAS(J)-CVAL)
            ENDDO
         ENDDO
CRVK         WRITE(6,*) 'ITER ',ITER, '  CHI2= ',CHI2, ' CVAL= ',CVAL
         IF(ABS(CHI2L-CHI2).LT.0.001) GOTO 100
      ENDDO
*
*     No convergence if we come here
      CHI2=-1.D0
      CL=-1.D0
      IERR = -1
*     
*     Return result
*
 100  ERR2N=ERR2P               ! return symmetric error
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(CHI2),NMEFF-1)) ! is this correct ?
* print out final scale factors
        DO I=1,NMEAS
           WRITE(6,1001) 'chi2_sym_rvk final scalefactors for ',
     +       CHANAL(I)(:LENOCC(CHANAL(I)))
           WRITE(6,1002) (SCALF(I,J),J=-1,NCSYS)
 1001      FORMAT(A,A)
 1002      FORMAT(12F6.3)
        ENDDO
      END

************************************************************************
*
      SUBROUTINE PREPARE_CHI2_RVK(VAL,W,WU,WS,WC,SCALF,IERR)
*     =====================================
*
*     Prepare for chi2 fit.
*     Return number of unknown parameters and inverse of error matrix.
      IMPLICIT NONE
      INCLUDE 'master.inc'
      INCLUDE 'combos.inc'
*
*     Arguments
*
      INTEGER IERR
      DOUBLE PRECISION VAL,W(MMEAS,MMEAS),WU(MMEAS,MMEAS),
     +   WS(MMEAS,MMEAS),WC(MMEAS,MMEAS,MCSYS)
      DOUBLE PRECISION SCALF(MMEAS,-1:NCSYS),STATI,USYSI,CSYSIK,
     &   STATJ,CSYSJK,S1,S2
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I,J,K
*
*     Build total error matrix and each submatrix.  Smoothly translate
*     from negative to positive errors in the case of asymmetric errors.
*
      DO I=1,NMEAS
* here we assume STATN and STATP have the same sign
        IF(STATN(I)*STATP(I).LT.0) WRITE(6,*) 'statn*statp<0',I
        STATI = STATN(I)+(1+TANH(2*(VAL-MEAS(I))/(STATP(I)+STATN(I))))
     +         /2*(STATP(I)-STATN(I))
* here we assume STATN and STATP have the same sign
        IF(USYSN(I)*USYSP(I).LT.0) WRITE(6,*) 'usysn*usysp<0',I
        USYSI = USYSN(I)+(1+TANH(2*(VAL-MEAS(I))/(USYSP(I)+USYSN(I))))
     +         /2*(USYSP(I)-USYSN(I))
        DO J=1,NMEAS
          STATJ = STATN(J)+(1+TANH(2*(VAL-MEAS(J))/(STATP(J)+STATN(J))))
     +         /2*(STATP(J)-STATN(J))
          WU(I,J)=0
          WS(I,J)=STACOR(I,J)*STATI*STATJ*SCALF(I,-1)*SCALF(J,-1)
          W(I,J)=WS(I,J)
          DO K=1,NCSYS
* here we assume CSYSN and CSYSP have the OPPOSITE sign
            IF(CSYSN(I,K)*CSYSP(I,K).GT.0)
     +    WRITE(6,*) K,CSYSN(I,K),CSYSP(I,K)
            CSYSIK = -CSYSN(I,K)+(1+TANH(2*(VAL-MEAS(I))/
     +         ABS(CSYSP(I,K)-CSYSN(I,K)))) /2 *(CSYSP(I,K)+CSYSN(I,K))
            CSYSJK = -CSYSN(J,K)+(1+TANH(2*(VAL-MEAS(J))/
     +         ABS(CSYSP(J,K)-CSYSN(J,K)))) /2 *(CSYSP(J,K)+CSYSN(J,K))
c          WRITE(6,*) I,K,VAL-MEAS(I),
c     +         2*(VAL-MEAS(I))/ABS(CSYSP(I,K)-CSYSN(I,K)),
c     +         CSYSN(I,K),CSYSP(I,K),CSYSIK

* The next bit is a bit tricky....  There is a problem with assuming
* two errors to be completely correlated when they are not really; you
* get SMALLER uncertainties on the average than you would at smaller
* values of correlation. In particular, the largest error comes when
* the correlation between two measurements is given by
* 
* rho = ( (s1^2+s2^2)/(s1*s2) - sqrt( (s1^2+s2^2)^2/(s1*s2)^2 - 4 ) )/2
* or, for the covariance element
* rho*s1*s2 = ( (s1^2+s2^2) - sqrt( (s1^2+s2^2)^2 - 4*(s1*s2)^2 ) )/2

            IF (I.NE.J .AND. XON(LMXCR) ) THEN
               S1=CSYSIK*SCALF(I,K)
               S2=CSYSJK*SCALF(J,K)
               WC(I,J,K)=0.5*( (S1**2+S2**2) -
     +              SQRT( (S1**2+S2**2)**2 - 4*(S1*S2)**2 ) )
C               IF(S1*S2.GT.0) WRITE(6,*) I,J,K,WC(I,J,K)/ABS(S1*S2)
            ELSE
               WC(I,J,K)=CSYSIK*SCALF(I,K)*CSYSJK*SCALF(J,K)
            ENDIF
            W(I,J)=W(I,J)+WC(I,J,K)
          ENDDO
        ENDDO
        WU(I,I)=USYSI**2*SCALF(I,0)*SCALF(I,0)
        W(I,I)=W(I,I)+WU(I,I)
      ENDDO
*
*     Invert total error matrix
*
      CALL DSINV(NMEAS,W,MMEAS,IERR)
      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'cannot invert error matrix',IERR
 1000   FORMAT(1X,A,': error returned by ',A,' ! ',A,'=',I6)
        RETURN 
      ENDIF 
*
      END
