************************************************************************
*
*     "CHI2_SYM_CIRC_G" combination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.20, April 3, 1997:
*     -- original release (cloned from CHI2_SYM_CIRC)
*     Version 2.30, May 12, 1997:
*     -- use NMEFF instead of NMEAS in computation of chi2 probability
*     Version 2.53, November 27, 1999:
*     -- add correlation between parameters FLAMB and FBS_FROM_BR; 
*        this correlation is taken from "parameter" CORR_FLAMB_FBS_
*        (which is assumed to be zero if not defined)
*     -- add one more digit to DMD printout
*     Version 3.00, January 29, 2002:
*     -- parameter MAXPAR increased from 50 to 100
*
************************************************************************
*
      SUBROUTINE CHI2_SYM_CIRC_G(CVAL,ERR2P,ERR2N,CL,IERR)
*     ====================================================
*
*     Olivier Schneider, CERN/PPE-ALE
*     April 3, 1997
*
*     Combination routine performing a chi2 fit in MINUIT to 
*     determine the average value of DMD:
*       - only symmetric uncertainties are used
*       - can handle correlated statistical uncertainties with any
*          correlation coefficient
*       - can handle 100% or 0% systematic correlations
*       - CHID --> FS --> DMD circularity handled correctly
*
*     Note: this routine is similar to CHI2_SYM_MIN, except that the
*           fs circulartiy is taken into account
*
      IMPLICIT NONE
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR
      EXTERNAL FCN_CHI2_SYM_CIRC_G_INIT,FCN_CHI2_SYM_CIRC_G
*
      CALL MINUIT_CHI2(FCN_CHI2_SYM_CIRC_G_INIT,FCN_CHI2_SYM_CIRC_G,
     &                 CVAL,ERR2P,ERR2N,CL,IERR)
      END
*
************************************************************************
*
      SUBROUTINE FCN_CHI2_SYM_CIRC_G(NPAR,D,F,X,IFLAG,FUNC)
*     =====================================================
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
      INTEGER NPAR,IFLAG
*
*     Argument of entry FCN_CHI2_SYM_CIRC_G_INIT
*
      INTEGER SPECIAL,IERR
*
*     Externals
*
      INTEGER LENOCC
      REAL PROB
*
*     Common blocks
*
      INCLUDE 'master.inc'
      INCLUDE 'circularity_g.inc'
*
*     Other static (saved) local variables
*
      INTEGER MAXPAR
      PARAMETER(MAXPAR=100)
      DOUBLE PRECISION CONCEN(MAXPAR),CONSIG(MAXPAR),CORRF
      SAVE             CONCEN        ,CONSIG        ,CORRF
      DOUBLE PRECISION W(MMEAS,MMEAS)
      SAVE             W
      LOGICAL NOT_READY,XCORRF
      SAVE    NOT_READY,XCORRF
      DATA    NOT_READY,XCORRF/.TRUE.,.FALSE./
      INTEGER LCSYS,LPAR,MPAR
      SAVE    LCSYS,LPAR,MPAR
      DATA    LCSYS,LPAR,MPAR/0,0,0/
      INTEGER IODX(NOLDPAR) ! MINUIT indexes of replaced parameters 
      SAVE    IODX
      DATA    IODX/NOLDPAR*0/
      INTEGER INDX(NNEWPAR) ! External MINUIT indexes of parameters 
*                           !  on which replaced parameters depend
      INTEGER JNDX(NNEWPAR) ! Internal MINUIT indexes of parameters
*                           !  on which replaced parameters depend
      SAVE    INDX,JNDX
      DATA    INDX,JNDX/NNEWPAR*0,NNEWPAR*0/
      DOUBLE PRECISION REF,EREF2
      SAVE             REF,EREF2
      DATA             REF,EREF2/0.D0,0.D0/
      INTEGER NNEWPAR2
      PARAMETER(NNEWPAR2=NNEWPAR*NNEWPAR)
      DOUBLE PRECISION EXIN(NNEWPAR,NNEWPAR)
      SAVE             EXIN
      DATA             EXIN/NNEWPAR2*0.D0/
*
*     Other non-static local variables
*
      INTEGER SPECIAL_FLAG
      PARAMETER(SPECIAL_FLAG=9)
      INTEGER I,J,IPAR
      CHARACTER*10 CHNAM
      LOGICAL ADD,SHOW
      DOUBLE PRECISION XIN(NNEWPAR),XOUT(NOLDPAR)
      DOUBLE PRECISION EXOUT(NOLDPAR,NOLDPAR),CORR(NOLDPAR,NNEWPAR)
      DOUBLE PRECISION CHI2,DUMMY,CONTRIB
      DOUBLE PRECISION Z(MMEAS)
      DOUBLE PRECISION V(MAXPAR),EMAT(MAXPAR,MAXPAR)
      DOUBLE PRECISION DMDLEP(2),TAU(2),CHID4S(2),XD4S(2),DMD4S(2),
     &                 WLEP,W4S,SUMW,DMDW(2)
      INTEGER K1,K2,J1,J2
      INTEGER NTOTPAR,MTOTPAR
      PARAMETER(MTOTPAR=NNEWPAR+NOLDPAR)
      REAL ETOT(MTOTPAR,MTOTPAR)
      INTEGER INDTOT(MTOTPAR)
      CHARACTER*10 CHERR(0:MTOTPAR),CHERN(0:MTOTPAR)
*
*     Check that routine was initialized
*
      IF(NOT_READY) THEN
        PRINT * ,'FCN_CHI2_SYM_CIRC_G: not initialized !!!'
        STOP 566
      ENDIF
*
      IF(IFLAG.EQ.1) THEN ! dummy initialization of FCN
      ENDIF
*
*     Compute derivatives for MINUIT
*
      IF(IFLAG.EQ.2) THEN
        PRINT * ,'FCN_CHI2_SYM_CIRC_G: cannot compute derivatives'
        STOP 669
      ENDIF
*
*     Check number of active parameters
*
      IF(NPAR.NE.LPAR) THEN 
        PRINT * ,'FCN_CHI2_SYM_CIRC_G: problem NPAR,LPAR = ',NPAR,LPAR
        STOP 543
      ENDIF
*
*     Copy MINUIT parameters in array V
*
      CALL DVCPY(MPAR,X(1),X(2),V(1),V(2))
*
*     Compute the replaced parameters from the parameters on which they depend
*
      DO I=1,NNEWPAR
        IF(INDX(I).GT.0) THEN 
          XIN(I)=X(INDX(I))
        ELSE 
          XIN(I)=0.D0
        ENDIF
      ENDDO
*
      IF(IFLAG.EQ.SPECIAL_FLAG.OR.IFLAG.EQ.3) THEN
*
*       Get full covariance matrix of fit and build covariance matrix of XIN
*
        CALL MNEMAT(EMAT,MAXPAR)
        DO I=1,NNEWPAR
          DO J=1,NNEWPAR
            IF(JNDX(I).GT.0.AND.JNDX(J).GT.0) THEN
              EXIN(I,J)=EMAT(JNDX(I),JNDX(J))
            ELSE 
              EXIN(I,J)=0.D0
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*
*     Compute XOUT from XIN
*
      CALL CIRCULARITY_G(INDX,XIN,EXIN,XOUT,EXOUT,CORR)
*
*     Replace constant zero values of replaced parameters in vector V
*     with the results of the above computation
*
      DO I=1,NOLDPAR
        IF(IODX(I).GT.0) THEN
          IF(V(IODX(I)).NE.0.D0) PRINT *,'Problem @@@ ',I,V(IODX(I))
          V(IODX(I))=XOUT(I)
        ENDIF
      ENDDO
*
*     Normalize array V
*
      DO I=1,MPAR
        V(I)=(V(I)-CONCEN(I))/CONSIG(I)
      ENDDO
*
*     Compute Z = "Delta^T * V + X", i.e. adjust original measurements
*
      CALL DVCPY(NMEAS,MEAS(1),MEAS(2),Z(1),Z(2))
      CALL DMMLA(NMEAS,LCSYS,1,CSYS(1,1),CSYS(1,2),CSYS(2,1), ! transposed
     &                         V   (  1),DUMMY    ,V   (  2),
     &                         Z   (  1),DUMMY    ,Z   (  2),DUMMY)
*
*     Compute CHI2
*
      CHI2=0.D0
      DO I=1,NMEAS
        DO J=1,NMEAS
          CHI2=CHI2+Z(I)*W(I,J)*Z(J)
        ENDDO
      ENDDO
*
*     Add constraint for each parameter, except for ...
*
      DO I=1,MPAR
        ADD=I.NE.LCSYS              !               ... measured quantity (dmd)
        DO J=1,NOLDPAR
          ADD=ADD.AND.I.NE.IODX(J)  !               ... replaced parameters
        ENDDO
        IF(ADD) THEN
          IF(XCORRF.AND.(I.EQ.INDX(IFLAMB).OR.I.EQ.INDX(IFSBR))) THEN
            CHI2=CHI2+(V(I)**2-V(INDX(IFLAMB))*V(INDX(IFSBR))*CORRF)/
     &                                                (1.D0-CORRF**2)
          ELSE
            CHI2=CHI2+V(I)**2
          ENDIF
        ENDIF
      ENDDO
*
*     Function to be minimized
*
      F=CHI2
*
      IF(IFLAG.NE.3.OR.LUNIT.LE.0) RETURN
*
*     Finalization: print all special parameters and their uncertainties
*     =============
*
*     Print final results
*
      WRITE(LUNIT,3000) 
 3000 FORMAT(/,1X,'Circularity solution:',
     &       /,1X,'---------------------',/,
     &       /,1X,'Number Parameter_name',T26,'  Fitted_value',
     &                                    T55,'  Constraint_used')
 3001 FORMAT(1X,I6,1X,A,T26,F10.5,' +- ',F10.5,
     &                  T55,F10.5,' +- ',F10.5)
 3011 FORMAT(1X,6X,1X,'    ',A,T55,F10.5)
 3002 FORMAT(1X,I6,1X,A,T26,F10.5,' +- ',F10.5,
     &                  T53,'function of other parameters')
      NTOTPAR=0
      CHERN(NTOTPAR)=' '
      DO I=1,NNEWPAR1
        IF(I.EQ.ICORRF) THEN
          IF(XCORRF) WRITE(LUNIT,3011) CHNEWPAR(I),CORRF
        ELSE IF(INDX(I).GT.0) THEN
          NTOTPAR=NTOTPAR+1
          INDTOT(NTOTPAR)=I
          CHERN(NTOTPAR)=CHNEWPAR(I)
          IF(INDX(I).EQ.LCSYS) THEN 
            WRITE(LUNIT,3001) INDX(I),CHNEWPAR(I),XIN(I),SQRT(EXIN(I,I))
          ELSE
            WRITE(LUNIT,3001) INDX(I),CHNEWPAR(I),XIN(I),SQRT(EXIN(I,I))
     &                               ,CONCEN(INDX(I)),CONSIG(INDX(I))
          ENDIF
        ENDIF
      ENDDO
      DO I=1,NOLDPAR
        SHOW=.TRUE.
        IF(INDX(ICHID4S).LE.0.AND.
     &     (I.EQ.ICHID1.OR.I.EQ.IXDW.OR.I.EQ.IDMDW)) SHOW=.FALSE.
        IF(INDX(IFSBR).LE.0.AND.I.EQ.IFS1) SHOW=.FALSE.
        IF(SHOW) THEN 
          NTOTPAR=NTOTPAR+1
          INDTOT(NTOTPAR)=-I
          CHERN(NTOTPAR)=CHOLDPAR(I)
          WRITE(LUNIT,3002) IODX(I),CHOLDPAR(I),XOUT(I),SQRT(EXOUT(I,I))
        ENDIF
      ENDDO
*
      WRITE(LUNIT,3003) NTOTPAR
 3003 FORMAT(/,1X,'Full error/correlation matrix on the above ',I2,
     &           ' parameters:',
     &   /,1X,'(elements on or below diagonal are from error matrix,',
     &   /,1X,' elements above diagonal are from correlation matrix)')
 3004 FORMAT(1X,16A)
*
*     Build error matrix
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
      WRITE(LUNIT,4000)
 4000 FORMAT(/,1X,'       Measurement   ',T26,
     & '  Chi2_contrib  Adjusted',T55,'  Input_measurement (stat)')
 4001 FORMAT(8X,A,I3,T26,F10.5,'    ',F10.5,
     &               T55,F10.5,' +- ',F10.5)
 4002 FORMAT(T20,'Total',T26,F10.5)
      DO I=1,NMEAS
        CONTRIB=0.D0
        DO J=1,NMEAS
          CONTRIB=CONTRIB+Z(I)*W(I,J)*Z(J)
        ENDDO
        WRITE(LUNIT,4001) CHMEAS(:LENOCC(CHMEAS)),I,
     &                    CONTRIB,Z(I)+XIN(IDMD),MEAS(I),STAT(I)
      ENDDO
      WRITE(LUNIT,4002) CHI2
*
      WRITE(LUNIT,5000) 
 5000 FORMAT(/,1X,'Final CHI2_SYM_CIRC_G result:')
      WRITE(LUNIT,5001) 
     & 'Average '//CHNEWPAR(IDMD)(:LENOCC(CHNEWPAR(IDMD))),
     & XIN(IDMD),SQRT(EXIN(IDMD,IDMD))
 5001 FORMAT(8X,41('-'),/,8X,A,T26,F10.6,' +- ',F10.6,/,
     &       8X,41('-'))
      IF(NMEFF.GT.1) THEN
      WRITE(LUNIT,5002) CHI2,NMEFF-1,CHI2/FLOAT(NMEFF-1),
     &                  PROB(SNGL(CHI2),NMEFF-1)
 5002 FORMAT(/,8X,'chi2 / NDF = ',T26,F10.5,' /',I3,' =',F10.5,
     &       /,8X,'prob(chi2) = ',T26,F10.5)
      ENDIF
      WRITE(LUNIT,5003)
 5003 FORMAT(/,1X,'... to be compared with CHI2_SYM result:',/)
      WRITE(LUNIT,5004)
     & 'Average '//CHNEWPAR(IDMD)(:LENOCC(CHNEWPAR(IDMD))),
     & REF,DSQRT(EREF2)
 5004 FORMAT(8X,A,T26,F10.5,' +- ',F10.5)
*
*     Compare world average of dmd with naive computation
*
      IF(INDX(ICHID4S).GT.0) THEN
      WRITE(LUNIT,6000) 
 6000 FORMAT(/,1X,'Extended average result including Upsilon(4S) data:')
      WRITE(LUNIT,5001)
     & 'Extended average '//CHNEWPAR(IDMD)(:LENOCC(CHNEWPAR(IDMD))),
     & XOUT(IDMDW),SQRT(EXOUT(IDMDW,IDMDW))
      WRITE(LUNIT,6003) CHNEWPAR(IDMD)(:LENOCC(CHNEWPAR(IDMD)))
 6003 FORMAT(/,1X,'... to be compared with naive computation ignoring',
     &       /,1X,'    correlation between average ',A,
     &            ' and Upsilon(4S) data:')
      DMDLEP(1)=XIN(IDMD)
      DMDLEP(2)=DSQRT(EXIN(IDMD,IDMD))
      TAU(1)=CONCEN(INDX(ITAUBD))
      TAU(2)=CONSIG(INDX(ITAUBD))
      CHID4S(1)=CONCEN(INDX(ICHID4S))
      CHID4S(2)=CONSIG(INDX(ICHID4S))
      XD4S(1)=1.D0/DSQRT(0.5D0/CHID4S(1)-1.D0)
      XD4S(2)=CHID4S(2)/(XD4S(1)*(1.D0-2.D0*CHID4S(1))**2)
      DMD4S(1)=XD4S(1)/TAU(1)
      DMD4S(2)=(XD4S(2)/XD4S(1))**2+(TAU(2)/TAU(1))**2
      DMD4S(2)=DMD4S(1)*DSQRT(DMD4S(2))
      WLEP=1.D0/DMDLEP(2)**2
      W4S =1.D0/DMD4S(2)**2
      SUMW=WLEP+W4S
      DMDW(1)=(WLEP*DMDLEP(1)+W4S*DMD4S(1))/SUMW
      DMDW(2)=1.D0/DSQRT(SUMW)
      WRITE(LUNIT,7001) 'Using','TAUBD'         ,TAU,
     &                  'and'  ,'CHID_UPSILON4S',CHID4S,
     &                  'yields','DMD_UPSILON4S',DMD4S
      WRITE(LUNIT,7001) 'Averaging','DMD_UPSILON4S',DMD4S,
     &                  'with'  ,'DMD_AVERAGE',DMDLEP,
     &                  'yields','DMD_EXTENDED',DMDW
 7001 FORMAT(/,(8X,A,T20,A,T34,' =',F10.5,' +- ',F10.5))
      ENDIF
*
*     Print PARAMETERS logical lines for possible use in subsequent COMBOS jobs
*
      WRITE(LUNIT,6001)
 6001 FORMAT(/,' PARAMETER logical lines',
     &         ' for possible use in subsequent COMBOS jobs:',/)
 6002 FORMAT(' PARAMETER ',A,3(2X,F10.5,SP),' ! CHI2_SYM_CIRC_G output')
      DO I=1,NNEWPAR
        IF(INDX(I).GT.0) WRITE(LUNIT,6002)
     &   CHNEWPAR(I),XIN(I),SQRT(EXIN(I,I)),-SQRT(EXIN(I,I))
      ENDDO
      DO I=1,NOLDPAR
        SHOW=.TRUE.
        IF(INDX(ICHID4S).LE.0.AND.
     &     (I.EQ.ICHID1.OR.I.EQ.IXDW.OR.I.EQ.IDMDW)) SHOW=.FALSE.
        IF(INDX(IFSBR).LE.0.AND.I.EQ.IFS1) SHOW=.FALSE.
        J=INDEX(CHOLDPAR(I),' ')
        IF(J.GT.0.AND.J.LT.LENOCC(CHOLDPAR(I))) SHOW=.FALSE.
        IF(SHOW) WRITE(LUNIT,6002)
     &   CHOLDPAR(I),XOUT(I),SQRT(EXOUT(I,I)),-SQRT(EXOUT(I,I))
      ENDDO
      WRITE(LUNIT,'(1X)')
      RETURN
*
      ENTRY FCN_CHI2_SYM_CIRC_G_INIT(SPECIAL,IERR)
*     ============================================
*
*     Do first simple combination without solving circularity
*
      CALL CHI2_SYM(REF,EREF2,DUMMY,DUMMY,IERR)
      IF(IERR.NE.0) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'error returned from CHI2_SYM',IERR
 1000   FORMAT(1X,A,': ',A,' ! IERR=',I6)
        RETURN
      ENDIF
*
*     Prepare (W matrix, ...)
*
      CALL PREPARE_CHI2(LCSYS,W,IERR)
      IF(LCSYS.NE.NCSYS+1) THEN 
        IERR=1000
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'inconsistency',IERR
      ENDIF
      IF(IERR.NE.0) RETURN
*
      LPAR=NCSYS+1 ! number of active MINUIT parameters
      MPAR=LPAR    ! total number of MINUIT parameters
      IF(MPAR.GT.MAXPAR) THEN 
        IERR=IERR+1
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'parameter MAXPAR too small',IERR
        RETURN
      ENDIF
*
*     Build list of central values and sigmas 
*     (for the constraints and systematic corrections)
*
      DO I=1,NCSYS
        CONCEN(I)=PARA(I)
        CONSIG(I)=EXCU(I)
      ENDDO
      CONCEN(LCSYS)=0.D0
      CONSIG(LCSYS)=1.D0
*
*     Loop on parameters to be replaced;
*     if these parameters are MINUIT parameters make them constant
*
      DO I=1,NOLDPAR
        IODX(I)=0
        DO J=1,NCSYS
          IF(CHPARA(J).EQ.CHOLDPAR(I)) IODX(I)=J
        ENDDO
        IF(CHMEAS.EQ.CHOLDPAR(I)) IODX(I)=LCSYS
        IF(IODX(I).GT.0) THEN 
          CALL MNPOUT(IODX(I),CHNAM,DUMMY,DUMMY,DUMMY,DUMMY,IPAR)
          IF(IPAR.LE.0) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'error returned by MNPOUT',IPAR
          ELSE IF(CHNAM.NE.CHOLDPAR(I)(:LEN(CHNAM))) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'inconsistency in old parameter names',IERR
          ENDIF
*
*         Turn MINUIT parameter into a constant with zero value
*  
          LPAR=LPAR-1
          CALL MNPARM(IODX(I),'(dummy)',0.D0,0.D0,0.D0,0.D0,J)
          IF(J.NE.0) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'error from MNPARM (constant)',J
          ENDIF
        ENDIF
      ENDDO
*
*     Loop on parameters on which the replaced parameters depend;
*     if these parameters are not already MINUIT parameters, add them
*
      CALL UZERO(EXIN,1,2*NNEWPAR2) ! double precision
      XCORRF=.FALSE.
      DO I=1,NNEWPAR1
        INDX(I)=0
        DO J=1,NCSYS
          IF(CHPARA(J).EQ.CHNEWPAR(I)) INDX(I)=J
        ENDDO
        IF(CHMEAS.EQ.CHNEWPAR(I)) INDX(I)=LCSYS
        DO J=NCSYS+1,NPARA
          IF(CHPARA(J).EQ.CHNEWPAR(I)) INDX(I)=-J
        ENDDO
        IF(INDX(I).EQ.0) THEN ! missing parameter
          IF(NEEDED(I)) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'parameter '//CHNEWPAR(I)//' not found',IERR
          ENDIF
        ELSE IF(INDX(I).GT.0) THEN ! already a MINUIT parameter
          CALL MNPOUT(INDX(I),CHNAM,DUMMY,DUMMY,DUMMY,DUMMY,IPAR)
          IF(IPAR.LE.0) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'error returned by MNPOUT',IPAR
          ELSE IF(CHNAM.NE.CHNEWPAR(I)(:LEN(CHNAM))) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'inconsistency in new parameter names',IERR
          ENDIF
          IF(INDX(I).EQ.LCSYS) THEN
            EXIN(I,I)=EREF2 ! guess for error on dmd
          ELSE
            EXIN(I,I)=EXCU(INDX(I))**2 ! parameter error
          ENDIF
        ELSE IF(I.EQ.ICORRF) THEN ! correlation between
*                                 ! parameters IFLAMB and IFSBR
          XCORRF=.TRUE.
          CORRF=PARA(IABS(INDX(I)))
          INDX(I)=0
          IF(DABS(CORRF).GE.1.D0) THEN ! invalid correlation coefficient
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'invalid value for '//CHNEWPAR(ICORRF),IERR
          ENDIF
        ELSE ! not a MINUIT parameter yet
*
*         Add MINUIT parameter
*
          IPAR=IABS(INDX(I))
          LPAR=LPAR+1
          MPAR=MPAR+1
          INDX(I)=MPAR
          CALL MNPARM(MPAR,CHPARA(IPAR),PARA(IPAR),
     &                0.1D0*EXCU(IPAR),0.D0,0.D0,J)
          EXIN(I,I)=EXCU(IPAR)**2
          IF(J.NE.0) THEN
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'error from MNPARM',J
          ENDIF
          IF(MPAR.GT.MAXPAR) THEN 
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'parameter MAXPAR too small',IERR
            RETURN
          ENDIF
          CONCEN(MPAR)=PARA(IPAR)
          CONSIG(MPAR)=EXCU(IPAR)
        ENDIF
      ENDDO
*
*     Determine whether a correlation between FLAMB and FBS_FROM_BR 
*     has to be taken into account
*
      XCORRF=XCORRF.AND.INDX(IFLAMB).GT.0.AND.INDX(IFSBR).GT.0
*
*     Determine internal indexes and set special flag
*
*     Note: the special flag is the value of IFLAG for calls to FCN
*           when the full error matrix is available
*
      SPECIAL=0
      DO I=1,NNEWPAR
        IF(INDX(I).GT.0) THEN
          IF(I.EQ.ICHID4S.OR.I.EQ.IFSBR) SPECIAL=SPECIAL_FLAG
          CALL MNPOUT(INDX(I),CHNAM,DUMMY,DUMMY,DUMMY,DUMMY,JNDX(I))
          IF(JNDX(I).LE.0.OR.CHNAM.NE.CHNEWPAR(I)(:LEN(CHNAM))) THEN
*            PRINT *, 'Fatal problem ...',INDX(I),JNDX(I),CHNAM
            IERR=IERR+1
            IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                     'fatal problem ...',IERR
            RETURN
          ENDIF
        ELSE
          JNDX(I)=0
        ENDIF
      ENDDO
*
*     Is everything OK ?
*
      NOT_READY=IERR.NE.0
      IF(NOT_READY.OR.LUNIT.LE.0) RETURN
*
*     Printout
*
      WRITE(LUNIT,8000) 
 8000 FORMAT(/,'The following parameters:')
      J=0
      DO I=1,NOLDPAR
        IF(IODX(I).GT.0) THEN
          WRITE(LUNIT,'(20X,A)') CHOLDPAR(I)
          J=J+1
        ENDIF
      ENDDO
      WRITE(LUNIT,8001) 
 8001 FORMAT(/,'will be replaced by functions of:')
      DO I=1,NNEWPAR
        IF(INDX(I).GT.0) THEN
          WRITE(LUNIT,'(20X,A)') CHNEWPAR(I)
        ELSE IF(I.EQ.ICHIS) THEN
          WRITE(LUNIT,'(20X,2A)') CHNEWPAR(I),' = 1/2'
        ELSE IF(I.EQ.IR.OR.I.EQ.IRS.OR.I.EQ.IRL) THEN
          WRITE(LUNIT,'(20X,2A)') CHNEWPAR(I),' = 1'
        ENDIF
      ENDDO
      WRITE(LUNIT,8002) 
 8002 FORMAT(/,'in order to solve the circularity',/)
      IF(XCORRF) WRITE(LUNIT,8012)
     & CHNEWPAR(IFLAMB)(:LENOCC(CHNEWPAR(IFLAMB))),
     & CHNEWPAR(IFSBR )(:LENOCC(CHNEWPAR(IFSBR ))),
     & CHNEWPAR(ICORRF)(:LENOCC(CHNEWPAR(ICORRF))),CORRF
 8012 FORMAT('Correlation coefficient between ',A,' and ',A,': ',
     &       A,' = ',F7.4,/)
      IF(J.LE.0) WRITE(LUNIT,8003)
 8003 FORMAT('Since no parameters will be replaced,',/,
     &       'the combined result is expected to be',/,
     &       'the same as the one from routine CHI2_SYM',/)
      END
