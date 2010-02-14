************************************************************************
*
*     "COMBY_LIFETIME" combination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.00, February 20, 1997:
*     -- original release (from Hans-Guenther Moser)
*         based on COMBY program used by the LEP B lifetime working group
*
************************************************************************
*
      SUBROUTINE COMBY_LIFETIME(CVAL,ERR2P,ERR2N,CL,IERR)
*     ===================================================
*
*     combination routine for b-lifetime averages, derived from
*     'comby', see B-lifetime group:
* 
*              http://wwwcn.cern.ch/~claires/lepblife.html
*
*     H.-G. Moser 13.2.97
*
*
*
*     some conventiones:
*
*     statistical and uncorrelated systematic errors may be asymmetric and
*     are treated properly then
*
*     any correlated systematic error and the parameter connected to it
*     should have symmetric errors (at least the program treats them as
*     symmetric)
*
*
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
      DOUBLE PRECISION ERR2P ! positive total uncertainty on CVAL
      DOUBLE PRECISION ERR2N ! negative total uncertainty on CVAL
      DOUBLE PRECISION CL    ! confidence level of combination (e.g. chi2 prob.)
*                              not calculated in this routine
      INTEGER IERR           ! error flag (non-zero if combination failed)
*
*     Externals
*
      INTEGER LENOCC
      EXTERNAL FCN_comby
*
*     Local variables
*
      INTEGER NFIT,I,IER,IVAR
      PARAMETER (NFIT = MCSYS+1)
C
      DOUBLE PRECISION PRTLEV
      DOUBLE PRECISION ERROR2,BND1,BND2,GLOBCC
      DOUBLE PRECISION START(NFIT),STEP(NFIT),BOUNDH(NFIT),BOUNDL(NFIT)
      DOUBLE PRECISION S_MEAN,S_WH
      CHARACTER*16 CHNAM
C
********
C
      IERR = -1
      PRTLEV = 1.D0        ! default MINUIT printout level
      DO I = I,NPARA
       IF(CHPARA(I).EQ.'MINUIT_PRTLEV') PRTLEV=PARA(I)
      ENDDO
      IF(LUNIT.LE.0) THEN
       PRTLEV = -1.D0 ! switch MINUIT output off
      ELSE
       WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),' called'
      ENDIF
C
C     get reasonable starting value
C
      S_MEAN = 0.
      S_WH   = 0.
      DO I = 1,NMEAS
       S_MEAN = S_MEAN + MEAS(I)/STAT(I)**2
       S_WH   = S_WH   + 1.0D0/STAT(I)**2
      ENDDO
      S_MEAN = S_MEAN/S_WH
C
C     prepare start, step and boundary values
C
      START(1) = S_MEAN
      STEP(1)  = 0.01D0
      BOUNDL(1) = 0.0D0
      BOUNDH(1) = 0.0D0
      DO I = 1,NCSYS
       START(I+1) = PARA(I)
       STEP(I+1)  = EXCU(I)/5.0D0
       BOUNDL(I+1)= 0.0D0
       BOUNDH(I+1)= 0.0D0
      ENDDO
C
C     call MINUIT
C
      CALL MNINIT(0,IABS(LUNIT),0)
      CALL MNEXCM(FCN_comby,'SET PRI',PRTLEV,1,IER,0)
C
      CALL MNSETI('combine lifetime results')
C
      CALL MNPARM(1,'lifetime',START(1),STEP(1),BOUNDL(1),
     &            BOUNDH(1),IER)
      IF(IER.NE.0) THEN
       IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                ': MINUIT failed'
       RETURN
      ENDIF
C
      DO I = 1,NCSYS
       CALL MNPARM(I+1,CHPARA(I),START(I+1),STEP(I+1),
     &             BOUNDL(I+1),BOUNDH(I+1),IER)
       IF(IER.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ': MINUIT failed'
        RETURN
       ENDIF
      ENDDO
C
      CALL MNEXCM(FCN_comby,'SET ERR',0.5D0,1,IER,0)
      CALL MNEXCM(FCN_comby,'CALL FCN',1.D0,1,IER,0)
C
      CALL MNEXCM(FCN_comby,'MIGRAD',0.0D0,1,IER,0)
      IF(IER.NE.0) THEN
       IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                ': MINUIT failed'
       RETURN
      ENDIF
      CALL MNEXCM(FCN_comby,'MINOS',0.0D0,1,IER,0)
      IF(IER.NE.0) THEN
       IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                ': MINUIT failed'
       RETURN
      ENDIF
C
      IERR = IER     
C
      CALL MNPOUT(1,CHNAM,CVAL,ERROR2,BND1,BND2,IVAR)
      ERROR2=ERROR2**2
      CALL MNERRS(1,ERR2P,ERR2N,ERROR2,GLOBCC)
      ERROR2=ERROR2**2
      ERR2P=ERR2P**2
      ERR2N=ERR2N**2
C
      CL =-1.0D0 ! confidence level not available
      IF(IERR.NE.0) THEN ! combination failed, return dummy values
           CVAL = 0.D0
           ERR2P = 0.D0
           ERR2N = 0.D0
      ENDIF
      RETURN
      END
C
C================================================================
C
      SUBROUTINE FCN_comby(NPAR,G,F,P,IFLAG)
C
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      DOUBLE PRECISION P(*),G,F
C
      INCLUDE 'master.inc'
C
      INTEGER I,J
      DOUBLE PRECISION T(MMEAS),A(MMEAS),B(MMEAS),SP(MMEAS),
     &                 SM(MMEAS),TC(MMEAS),BM,BP
C
C     initialization
C
      IF(IFLAG.LT.2) THEN
C
       DO I = 1,NMEAS
C
C       combine statistical and uncorrelated systematical error
C
        SP(I) = DSQRT(STATP(I)**2 + USYSP(I)**2)
        SM(I) = DSQRT(STATN(I)**2 + USYSN(I)**2)
C
C       if negative error larger then positive: average
C
        IF(SP(I).LT.SM(I)) THEN
          SP(I) = 0.5D0 * (SP(I) + SM(I))
          SM(I) = SP(I)
        ENDIF
        T(I) = MEAS(I)
C
C       calculate asymmetric and symmetric coefficients
C
        BM = SM(I) / (T(I)-SM(I)) + DLOG((T(I)-SM(I))/T(I))
        BP =-SP(I) / (T(I)+SP(I)) + DLOG((T(I)+SP(I))/T(I))
        A(I) = 0.5D0*(SP(I)**2-SM(I)**2)/(BM*SP(I)**2-BP*SM(I)**2)
        B(I) = (0.5D0 - A(I)*BM) / SM(I)**2
       ENDDO
      ENDIF  ! IFLAG .LT. 2
C
C
      F = 0.0D0
C
      DO I = 1,NMEAS
C
       TC(I) = T(I)
C
       DO J =1, NCSYS
C
          TC(I) = TC(I) + CSYS(I,J)*(P(J+1)-PARA(J))/EXCU(J)
C
C          IF((P(J+1)-PARA(J)).GT.0.D0) 
C     &    TC(I) = TC(I) + CSYSP(I,J)*(P(J+1)-PARA(J))/EXCUP(J)
C          IF((P(J+1)-PARA(J)).LE.0.D0) 
C     &    TC(I) = TC(I) + CSYSN(I,J)*(P(J+1)-PARA(J))/EXCUN(J)
C
       ENDDO
C
        F = F + A(I)*( (TC(I)-P(1))/P(1) + DLOG(P(1)/TC(I)) )
     &        + B(I)*( P(1)-TC(I) )**2
C
      ENDDO
C
C     add constraints due to parameter
C
      DO I = 1,NCSYS
        F = F + 0.5D0 * ( (P(I+1)-PARA(I))/EXCU(I) )**2
C
C         IF((P(I+1)-PARA(I)).GT.0.D0)
C     &   F = F + 0.5D0 * ( (P(I+1)-PARA(I))/EXCUP(I) )**2
C         IF((P(I+1)-PARA(I)).LE.0.D0)
C     &   F = F + 0.5D0 * ( (P(I+1)-PARA(I))/EXCUN(I) )**2
      ENDDO
C
      RETURN
      END
