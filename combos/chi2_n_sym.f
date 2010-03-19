************************************************************************
*
*     "CHI2_N_SYM" combination routine for COMBOS program
*       copied by David ROUSSEAU from 
*     CHI2_SYM Olivier Schneider, CERN/PPE-ALE
*       to be able to average measurement of different quantities

*     Version 2.40, June 1, 1999
************************************************************************
*
      SUBROUTINE CHI2_N_SYM(CVAL,ERR2P,ERR2N,CL,IERR)
*     ===============================================
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
*      David Rousseau CPPM, 
*        - Several experiments measure different quantities
*          June 1st 1999
*
*      Olivier Schneider: replaced 2nd and 3rd arguments with 
*                         real error matrices
*                         (with diagonal elements equal to error**2)
*
*     Input data: available from a common block in file master.inc
*     ----------
*
*DR  
      IMPLICIT NONE
      INCLUDE 'master.inc'
*      
*     Output data: returned in the arguments 
*     -----------
*
      DOUBLE PRECISION CVAL(MQUAN) ! combined value (or average value)
      DOUBLE PRECISION ERR2P(MQUAN,MQUAN) ! positive total uncertainty**2 on CVAL
      DOUBLE PRECISION ERR2N(MQUAN,MQUAN) ! negative total uncertainty**2 on CVAL
      DOUBLE PRECISION CL   ! confidence level of combination (e.g. chi2 prob.)
      INTEGER IERR          ! error flag (non-zero if combination failed)
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
      INTEGER LCSYS,I,J
      DOUBLE PRECISION V(MCSYS),W(MMEAS,MMEAS),
     &                 CHI2,DUMMY(MCSYS),Z(MMEAS),
     &                 TEMP(MMEAS,MCSYS),S(MCSYS,MCSYS)

*
*     Default output arguments
*
*DR      CVAL = 0.D0
*      ERR2P = 0.D0
*      ERR2N = 0.D0
*DR
      CALL VZERO(CVAL,MQUAN*2)
      CALL VZERO(ERR2P,MQUAN*MQUAN*2)
      CALL VZERO(ERR2N,MQUAN*MQUAN*2)
*DR end
      CL   = -1.D0
      IERR = 9999
      CALL PREPARE_CHI2(LCSYS,W,IERR)
      IF(IERR.NE.0) RETURN

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
*     Compute TEMP = CSYS*W = "Delta * M**(-1)"
*

      CALL DMMLT(LCSYS,NMEAS,NMEAS,CSYS(1,1),CSYS(2,1),CSYS(1,2),
     &                             W   (1,1),W   (2,1),W   (1,2),
     &                             TEMP(1,1),TEMP(2,1),TEMP(1,2),DUMMY)
*
*     Compute S = error matrix on V = "Delta*M**(-1)*Delta^T+Q)**(-1)"
*
      CALL UZERO(S,1,2*MCSYS*MCSYS) ! double precision
      DO I=1,NCSYS ! not LCSYS !
        S(I,I)=1.D0
      ENDDO
*should add correlation between systematics there ?

*DR      CALL DMMLA(LCSYS,NMEAS,LCSYS,TEMP(1,1),TEMP(2,1),TEMP(1,2),
      CALL DMMLT(LCSYS,NMEAS,LCSYS,TEMP(1,1),TEMP(2,1),TEMP(1,2),
     &                             CSYS(1,1),CSYS(1,2),CSYS(2,1), ! transposed
     &                             S   (1,1),S   (2,1),S   (1,2),DUMMY)

      DO I=1,NCSYS ! not LCSYS !
        S(I,I)=S(I,I)+1.D0
      ENDDO

      CALL DSINV(LCSYS,S,MCSYS,IERR)
      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1001) CHROUT(:LENOCC(CHROUT)),IERR
 1001   FORMAT(1X,A,': cannot invert matrix ! IERR=',I6)
        RETURN 
      ENDIF 
*
*     Compute TEMP = S*TEMP = "S * Delta * M**(-1)"
*
      CALL DMMLT(LCSYS,LCSYS,NMEAS,S   (1,1),S   (2,1),S   (1,2),
     &                             TEMP(1,1),TEMP(2,1),TEMP(1,2),
     &                             TEMP(1,1),TEMP(2,1),TEMP(1,2),DUMMY)
*
*     Compute V = TEMP*MEAS = "S * Delta * M**(-1) * X"
*     Note: this V here is corresponds to -V of the documentation
*
      CALL DMMLT(LCSYS,NMEAS,1,TEMP(1,1),TEMP(2,1),TEMP(1,2),
     &                         MEAS(  1),DUMMY    ,MEAS(  2),
     &                         V   (  1),DUMMY    ,V   (  2),DUMMY)
*
*     Compute Z = "Delta^T * V - X"
*     Note: this Z here corresponds to -(Delta^T * V + X) of the documentation
*
      CALL DVCPY(NMEAS,MEAS(1),MEAS(2),Z(1),Z(2))
      CALL DMMLS(NMEAS,LCSYS,1,CSYS(1,1),CSYS(1,2),CSYS(2,1), ! transposed
     &                         V   (  1),DUMMY    ,V   (  2),
     &                         Z   (  1),DUMMY    ,Z   (  2),DUMMY)
*
*     Compute CHI2
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
*     Return result
*

*DR      CVAL=-V(LCSYS) ! note the sign !
*      ERR2P=S(LCSYS,LCSYS)
*      ERR2N=ERR2P ! return symmetric error
*      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(CHI2),NMEFF-1)) ! is this correct ?
*DR
      DO I=1,NQUAN
        CVAL(I)=-V(NCSYS+I) ! note the sign !
        DO J=1,NQUAN
          ERR2P(I,J)=S(NCSYS+I,NCSYS+J)
          ERR2N(I,J)=ERR2P(I,J) ! return symmetric error
        ENDDO
      ENDDO
C      Print *, 'my chi2', chi2, NMEFF,NQUAn 
      IF (NMEFF.GT.1) THEN
         PRINT *, 'CHI2_N_SYM: CHI2, NMEFF, NQUAN, CHI2/NDOF  = ',
     &                         CHI2, NMEFF, NQUAN, CHI2/(NMEFF-NQUAN)
      ELSE
         PRINT *, 'CHI2_N_SYM: CHI2, NMEFF, NQUAN, CHI2/NDOF  = ',
     &             CHI2, NMEFF, NQUAN, ' set to 0.0 '
      ENDIF
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(CHI2),NMEFF-NQUAN)) ! is this correct ?
      CALL DUMP_FIT(V,S)
*DR end
      END

*DR routine to dump full results from covariance matrix
      SUBROUTINE DUMP_FIT(V,S)
* Input : V vector of fitted parameter (1:NCSYS) systematics
*                                     (NCSYS+1:NCSYS+NQUAN) quantities
*         S  error matrix
*         SINV inverse of error matrix 
      IMPLICIT NONE
      INCLUDE 'master.inc'
      INCLUDE 'combos.inc' 
 
      INTEGER IERR
      DOUBLE PRECISION V(MCSYS),S(MCSYS,MCSYS),a
      DOUBLE PRECISION SINV(MCSYS,MCSYS)
      DOUBLE PRECISION SAUX(MCSYS,MCSYS),SSYS(MCSYS)
      DOUBLE PRECISION SQINV(MQUAN,MQUAN),SQ(MQUAN,MQUAN) 
      DOUBLE PRECISION SYSSUMTOT(MQUAN),SYSTOT(MQUAN),ERRTOT(MQUAN)
      DOUBLE PRECISION STATTOT(MQUAN)
      DOUBLE PRECISION PRTLEV
      INTEGER LENOCC
      EXTERNAL LENOCC            

      INTEGER LCSYS,I,J,IQUAN,JQUAN,IVAR
      CHARACTER*16 CHAUX(MCSYS)

      LCSYS=NCSYS+NQUAN
*
*     Determine debug printout level
*
      PRTLEV=-1.D0 ! default printout level: none
      DO I=NCSYS+1,NPARA
         IF(CHPARA(I).EQ.'CHI2_N_SYM_PRT') PRTLEV=PARA(I)
      ENDDO
      IF (PRTLEV.GT.0) THEN 
         PRINT*,'SwB: NCSYS,NQUAN,LCSYS,LUNIT= ',NCSYS,NQUAN,LCSYS,LUNIT
         DO I=1,LCSYS 
            DO J=1,LCSYS 
               PRINT *, 'SwB: I, J, SIJ = ',I,J,S(I,J)
            ENDDO
         ENDDO
      ENDIF

      DO I=1,LCSYS
        IF (I.LE.NCSYS) THEN
          CHAUX(I)=CHPARA(I)
        ELSE
          CHAUX(I)=CHQUAN(I-NCSYS)
        ENDIF
        DO J=1,LCSYS
          IF (I.EQ.J) THEN
*diagonal element:error
            SAUX(I,J)=SQRT(S(I,J))
          ELSE 
*off-diagonal: correlation coeff.
            SAUX(I,J)=S(I,J)/
     &          SQRT(S(I,I)*S(J,J))
          ENDIF
        ENDDO
      ENDDO 

      IF (PRTLEV.GT.0) THEN
         DO I=1,LCSYS 
            DO J=1,LCSYS 
               PRINT *, 'SwB: I,J,SAUX(IJ) = ',I,J,SAUX(I,J)
            ENDDO
         ENDDO
      ENDIF

      WRITE(LUNLOG,*)' '
      WRITE(LUNLOG,*)
     &  ' DUMP FIT RESULTS: diag=toterror',
     &  ' off-diag=corr. coeff=corr_sys/toterror'

      WRITE(LUNLOG,'(T25,20A14)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)

      DO I=1,LCSYS
        WRITE(LUNLOG,1000)CHAUX(I)(:LENOCC(CHAUX(I))),
     &   -V(I),(SAUX(I,J),J=1,LCSYS)
 1000  FORMAT(1X,A,T10,F14.7,' | ',T25,85F14.7)
      ENDDO
      WRITE(LUNLOG,'(T25,20A14)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)
      WRITE(LUNLOG,*)' '


*DR compute S inverse matrix
      CALL UCOPY(S,SINV,MCSYS*MCSYS*2)
      CALL DSINV(LCSYS,SINV,MCSYS,IERR)

      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1001) CHROUT(:LENOCC(CHROUT)),IERR
 1001   FORMAT(1X,A,': cannot reinvert error  matrix ! IERR=',I6)
        RETURN 
      ENDIF 

*compute intrinsic quantity error
*extract quantitiy matrix from inverse error matrix
      a=0.
      DO I=1,NQUAN
        DO J=1,NQUAN
          SQINV(I,J)=SINV(NCSYS+I,NCSYS+J)
        ENDDO
      ENDDO        
     
*inverse
      CALL UCOPY(SQINV,SQ,MQUAN*MQUAN*2)
      CALL DSINV(NQUAN,SQ,MQUAN,IERR)
      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1002) CHROUT(:LENOCC(CHROUT)),IERR
 1002   FORMAT(1X,A,': cannot reinvert Q error  matrix ! IERR=',I6)
        RETURN 
      ENDIF 
*now SQ is the statistical error matrix


*DR compute all error breakdown
      IF (NCSYS.EQ.0) THEN
        WRITE(LUNLOG,*)'--Averaged parameters, using only stat---'
      ELSE
        WRITE(LUNLOG,*)'--Averaged parameters, using stat+syst---'
      ENDIF
      WRITE(LUNLOG,*)'--breakdown of errors: syst then stat-----'
      WRITE(LUNLOG,'(T25,20A14)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)

      DO IQUAN=1,NQUAN
         IF (PRTLEV.GT.0) PRINT *, ' '
*statistical part
        IVAR=NCSYS+IQUAN
        DO JQUAN=1,NQUAN
          IF (JQUAN.EQ.IQUAN) THEN
             SSYS(NCSYS+JQUAN)=SQRT(SQ(IQUAN,JQUAN))
             IF (PRTLEV.GT.0)
     &        PRINT*, 'SwB: SQ(IQUAN=',IQUAN,',JQUAN=',JQUAN,')=',
     &                     SQ(IQUAN,JQUAN),' ',
     &                'SSYS(',NCSYS+JQUAN,')=',SSYS(NCSYS+JQUAN)
          ELSE
             SSYS(NCSYS+JQUAN)=SQ(IQUAN,JQUAN)/
     &                 SQRT(SQ(IQUAN,IQUAN)*SQ(JQUAN,JQUAN))
             IF (PRTLEV.GT.0) 
     &       PRINT*, 'SwB: SQ(IQUAN=',IQUAN,',JQUAN=',JQUAN,')=',
     &                     SQ(IQUAN,JQUAN),' ',
     &               'SSYS(',NCSYS+JQUAN,')=',SSYS(NCSYS+JQUAN)
          ENDIF
        ENDDO

*total error
        ERRTOT(IQUAN)=SAUX(IVAR,IVAR)

*statistical error
        STATTOT(IQUAN)=SSYS(IVAR)

*systematic part
*compute systematics sum
        SYSSUMTOT(IQUAN)=0.D0
        DO I=1,NCSYS
          SSYS(I)=SAUX(IVAR,I)*SAUX(IVAR,IVAR)
          IF (PRTLEV.GT.0) 
     &       PRINT*, 'SwB: ',
     &         'SSYS(I=',I,')=',
     &         'SAUX(',IVAR,',',I,')*',
     &         'SAUX(',IVAR,',',IVAR,')=',SSYS(I)
          SYSSUMTOT(IQUAN)=SYSSUMTOT(IQUAN)+SSYS(I)**2 
          IF (PRTLEV.GT.0) 
     &     PRINT*,'IQUAN,SYSSUMTOT(IQUAN),SSYS(I)**2 = ',
     &            IQUAN,SYSSUMTOT(IQUAN),SSYS(I)**2
        ENDDO
        SYSSUMTOT(IQUAN)=SQRT(SYSSUMTOT(IQUAN))
        IF (PRTLEV.GT.0) THEN
        PRINT *, 'SYSSUMTOT(IQUAN) = sqrt()=',SYSSUMTOT(IQUAN)
        PRINT *, ' '
        ENDIF

*compute sytematics from total-statistical (should be more or less 
* equal to SYSUMTOT but not exactly, since systematics are now correlated
        SYSTOT(IQUAN)=SQRT(MAX(0.,ERRTOT(IQUAN)**2-STATTOT(IQUAN)**2))

        WRITE(LUNLOG,1000)CHAUX(IVAR)(:LENOCC(CHAUX(IVAR))),
     &     -V(IVAR),(SSYS(I),I=1,LCSYS)

      ENDDO

      WRITE(LUNLOG,'(T25,20A14)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)
      WRITE(LUNLOG,*)' '

*dump errors
c      IF (NCSYS.NE.0) THEN ! why not print results anyway?
        DO IQUAN=1,NQUAN
          IVAR=NCSYS+IQUAN
          WRITE(LUNLOG,1100)CHAUX(IVAR)(:LENOCC(CHAUX(IVAR))),
     &       -V(IVAR),STATTOT(IQUAN),SYSTOT(IQUAN),
     &       ERRTOT(IQUAN),SYSSUMTOT(IQUAN)
 1100     FORMAT(1X,A20,1X,F14.7,' +- ',F14.7,' +- ',F14.7, 
     &  ' Tot Err:',F14.7,' Check Sys:',F14.7)
        ENDDO
c      ENDIF


*compute correlation
      DO IQUAN=1,NQUAN
        DO JQUAN=1,NQUAN
          IF (IQUAN.LT.JQUAN) 
     &    WRITE(LUNLOG,
     & '('' Correlation between '',A,'' and '',A,'' = '',F14.7)')
     &          CHAUX(NCSYS+IQUAN)(:LENOCC(CHAUX(NCSYS+IQUAN))),
     &          CHAUX(NCSYS+JQUAN)(:LENOCC(CHAUX(NCSYS+JQUAN))),
     &          SAUX(NCSYS+IQUAN,NCSYS+JQUAN)
        ENDDO
      ENDDO
      
      WRITE(LUNLOG,*)' '

      RETURN
      END

 
*check error matrix is sensible
      SUBROUTINE CHECK_ERRMAT(NDEC,N,MAT)
      IMPLICIT NONE
      INTEGER NDEC,N,I,J,IER
      DOUBLE PRECISION MAT(NDEC,NDEC),AUX
      LOGICAL PROBLEM

      PROBLEM=.FALSE.
      
*check symmetry
      IER=0
      DO I=1,N
        DO J=1,I-1
          IF (MAT(I,J).NE.0.OR.MAT(J,I).NE.0) THEN
            AUX=(MAT(I,J)-MAT(J,I))/(MAT(I,J)+MAT(J,I))
            IF (AUX.GT.1D-6) THEN
         WRITE(*,*)'*DR WARNING NOT SYM MAT',I,J,AUX,MAT(I,J),MAT(J,I)
         PROBLEM=.TRUE.
            ENDIF
          ENDIF
        ENDDO
      ENDDO




      IF (PROBLEM) THEN
        CALL DUMMAT(MAT,N,N,NDEC,NDEC)
      ENDIF

      RETURN
      END

      SUBROUTINE DUMMAT(A,N,M,NT,MT)
C dump matrix A(N,M) with a nice format (NT,MT) TAILLE MATRICE
      IMPLICIT NONE
      DOUBLE PRECISION A(*)
      INTEGER NT,MT,N,M,I,J,IF,IL,NN,MM

      NN=MIN(20,N)
      MM=MIN(20,M)

      DO I=1,NN
       IF=1
       IL=MM
       WRITE(*,1000) (A((J-1)*NT+I),J=1,IL)
 1000  FORMAT (1X,20E10.3)
*1000  FORMAT (1X,20F6.3)
      ENDDO

      RETURN
      END
