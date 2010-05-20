*
***********************************************************************
*
*     "CHI2_N_SYM" combination routine for COMBOS program
*       copied by David ROUSSEAU from 
*     CHI2_SYM Olivier Schneider, CERN/PPE-ALE
*       to be able to average measurement of different quantities
*
*     Version 2.40, June 1, 1999
*
*     Version 3.33, March 20, 2010  
*     -- improved printout in DUMP_MASTER_INC and in CHI2_N_SYM
*        (provided by Swagato)
*
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
      INCLUDE 'combos.inc' 
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
     &                ,SCOPY(MCSYS,MCSYS)
      INTEGER          ErrorFlag, INVOPT, NMEFFNEW, NQUANNEW, NDOFNEW
      DOUBLE PRECISION PRTLEV
*
*     Determine debug printout level
*
      PRTLEV=-1.D0 ! default printout level: none
      DO I=1,NPARA
        IF(CHPARA(I).EQ.'CHI2_N_SYM_PRT') PRTLEV=PARA(I)
      ENDDO
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

      INVOPT=0 ! default option for matrix inversion using DSINV
      DO I=1,NPARA
        IF(CHPARA(I).EQ.'CHI2_N_SYM_INV') INVOPT=INT(PARA(I))
      ENDDO
      IF (INVOPT.EQ.0) THEN
        CALL DSINV(LCSYS,S,MCSYS,IERR)
        PRINT *, 'CHI2_N_SYM: DSINV: S->S: IERR = ',IERR 
      ELSE
        ErrorFlag=1
        DO I=1,LCSYS
          DO J=1,LCSYS
            if (ErrorFlag.eq.1.and.i.eq.j) print *, 
     &      'CHI2: i,j,s(i,j) = ',i,j,s(i,j),sqrt(max(0,s(i,j)))
            SCOPY(I,J)=S(I,J)
          ENDDO
        ENDDO
        CALL FindInv(SCOPY,S,LCSYS,MCSYS,ErrorFlag)
        PRINT *, 'CHI2_N_SYM: FindInv: S->S: ErrorFlag = ',ErrorFlag
        IERR=ErrorFlag
      ENDIF

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
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(CHI2),NMEFF-NQUAN)) ! is this correct ?
C     IF(PRTLEV.GT.0) THEN 
      IF (PRTLEV.GT.-999) THEN  ! SwB [default changed]
        WRITE(LUNLOG,
     &       '(''CHI2_N_SYM: CHI2, NMEFF, NQUAN, NDOF    = '','//
     &       'G10.5,3(1X,I3))') CHI2, NMEFF, NQUAN, NMEFF-NQUAN
        IF (NMEFF.GT.NQUAN) THEN
          WRITE (LUNLOG,
     &         '(''CHI2_N_SYM: CHI2/NDOF, SCALE FAC, CL    = '','//
     &         '3(1X,G10.5))')
     &         CHI2/(NMEFF-NQUAN),SQRT(MAX(0.,CHI2/(NMEFF-NQUAN))),
     &         CL
        ELSE
          WRITE (LUNLOG, *)
     &      'CHI2_N_SYM: CHI2/NDOF, SCALE FAC, CL set to 0.0'
        ENDIF
C     
        NMEFFNEW=NMEFF ! initialize
        NQUANNEW=NQUAN ! initialize
        DO I=1,NPARA
          IF(CHPARA(I).EQ.'CHI2_N_SYM_NDOF') THEN
            NMEFFNEW=INT(PARA(I))
            NQUANNEW=INT(EXCUP(I))
          ENDIF
        ENDDO
        NDOFNEW = NMEFFNEW - NQUANNEW
C
        WRITE (LUNLOG,
     &       '(''CHI2_N_SYM: CHI2, NMEFFNEW, NQUANNEW, NDOFNEW = '','//
     &       'G10.5,3(1X,I3))') CHI2, NMEFFNEW, NQUANNEW, NDOFNEW
        IF (NMEFFNEW.GT.NQUANNEW) THEN
          WRITE (LUNLOG,
     &         '(''CHI2_N_SYM: CHI2/NDOFNEW, SCALE FAC, CL = '','//
     &         '3(1X,G10.5))')
     &         CHI2/NDOFNEW,SQRT(MAX(0.,CHI2/NDOFNEW)),
     &         PROB(SNGL(CHI2),NDOFNEW)
        ELSE
          WRITE (LUNLOG,*)
     &      'CHI2_N_SYM: CHI2/NDOFNEW, SCALE FAC, CL set to 0.0'
        ENDIF
C     
      ENDIF
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
 
      INTEGER IERR, ErrorFlag, INVOPT
      DOUBLE PRECISION V(MCSYS),S(MCSYS,MCSYS)
      DOUBLE PRECISION SINV(MCSYS,MCSYS)
      DOUBLE PRECISION SAUX(MCSYS,MCSYS),SSYS(MCSYS)
      DOUBLE PRECISION SQINV(MQUAN,MQUAN),SQ(MQUAN,MQUAN) 
      DOUBLE PRECISION SYSSUMTOT(MQUAN),SYSTOT(MQUAN),ERRTOT(MQUAN)
      DOUBLE PRECISION STATTOT(MQUAN)
      DOUBLE PRECISION PRTLEV
      DOUBLE PRECISION AVESUM, AVEERR
      INTEGER LENOCC
      EXTERNAL LENOCC            

      INTEGER LCSYS,I,J,IQUAN,JQUAN,IVAR
      CHARACTER*16 CHAUX(MCSYS)

      INTEGER PSUM, ISUM, NSUM, II
      LOGICAL USEQUAN(MQUAN)
      DOUBLE PRECISION COEFF(MQUAN)
      CHARACTER*13 CHSUM
      CHARACTER*16 CHSUM2

      LCSYS=NCSYS+NQUAN
*
*     Determine debug printout level
*
      PRTLEV=-1.D0 ! default printout level: none
      DO I=1,NPARA
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
      WRITE(LUNLOG,'(T30,10A12)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)
      DO I=1,LCSYS
        WRITE(LUNLOG,1000)CHAUX(I)(:LENOCC(CHAUX(I))),
     &   -V(I),(SAUX(I,J),J=1,LCSYS)
 1000  FORMAT(1X,A12,T14,1X,F11.6,' | ',20(T30,10(1X,F11.6),/))
      ENDDO
      WRITE(LUNLOG,'(T30,10A12)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)
      WRITE(LUNLOG,*)' '

      INVOPT=0 ! default option for matrix inversion using DSINV
      DO I=1,NPARA
        IF(CHPARA(I).EQ.'CHI2_N_SYM_INV') INVOPT=INT(PARA(I))
      ENDDO
*DR compute S inverse matrix
      IF (INVOPT.EQ.0) THEN
        CALL UCOPY(S,SINV,MCSYS*MCSYS*2)
        CALL DSINV(LCSYS,SINV,MCSYS,IERR)
        PRINT *, 'DUMP_FIT: DSINV: S->SINV: IERR = ',IERR
      ELSE
        CALL FindInv(S,SINV,LCSYS,MCSYS,ErrorFlag)
        PRINT *, 'DUMP_FIT: FindInv: S->SINV: ErrorFlag = ',ErrorFlag
        IERR=ErrorFlag
      ENDIF

      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1001) CHROUT(:LENOCC(CHROUT)),IERR
 1001   FORMAT(1X,A,': cannot reinvert error  matrix ! IERR=',I6)
        RETURN 
      ENDIF 

*compute intrinsic quantity error
*extract quantitiy matrix from inverse error matrix
      DO I=1,NQUAN
        DO J=1,NQUAN
          SQINV(I,J)=SINV(NCSYS+I,NCSYS+J)
        ENDDO
      ENDDO        
     
*inverse
      IF (INVOPT.EQ.0) THEN
        CALL UCOPY(SQINV,SQ,MQUAN*MQUAN*2)
        CALL DSINV(NQUAN,SQ,MQUAN,IERR)
        PRINT *, 'DUMP_FIT: DSINV: SQINV->SQ: IERR = ',IERR
      ELSE
        CALL FindInv(SQINV,SQ,NQUAN,MQUAN,ErrorFlag)
        PRINT *, 'DUMP_FIT: FindInv: SQINV->SQ: ErrorFlag = ',ErrorFlag
        IERR=ErrorFlag
      ENDIF

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
      WRITE(LUNLOG,'(T30,10A12)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)

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

      WRITE(LUNLOG,'(T30,10A12)')(CHAUX(I)(:LENOCC(CHAUX(I))),I=1,LCSYS)
      WRITE(LUNLOG,*)' '

*dump errors
c      IF (NCSYS.NE.0) THEN ! why not print results anyway?
        DO IQUAN=1,NQUAN
          IVAR=NCSYS+IQUAN
          WRITE(LUNLOG,1100)CHAUX(IVAR)(:LENOCC(CHAUX(IVAR))),
     &       -V(IVAR),STATTOT(IQUAN),SYSTOT(IQUAN),
     &       ERRTOT(IQUAN),SYSSUMTOT(IQUAN)
 1100     FORMAT(1X,A20,1X,F14.7,' +/- ',F14.7,' +/- ',F14.7, 
     &  ' Tot Err:',F14.7,' Check Sys:',F14.7)
        ENDDO
c      ENDIF


*compute correlation
      DO IQUAN=1,NQUAN
        DO JQUAN=1,NQUAN
          IF (IQUAN.LT.JQUAN) 
     &    WRITE(LUNLOG,
     & '('' Correlation between '',A20,'' and '',A20,'' = '',F14.7)')
     &      CHAUX(NCSYS+IQUAN)(:AMIN0(20,LENOCC(CHAUX(NCSYS+IQUAN)))),
     &      CHAUX(NCSYS+JQUAN)(:AMIN0(20,LENOCC(CHAUX(NCSYS+JQUAN)))),
     &      SAUX(NCSYS+IQUAN,NCSYS+JQUAN)
        ENDDO
      ENDDO

*print average sum of multiple quantities

      AVESUM=0.0
      AVEERR=0.0
      DO IQUAN=1,NQUAN
        IVAR=NCSYS+IQUAN
        AVESUM = AVESUM - V (IVAR)
        AVEERR = AVEERR + ERRTOT(IQUAN)**2
        DO JQUAN=1,NQUAN
          IF (IQUAN.LT.JQUAN) THEN
            AVEERR = AVEERR + 2 * ERRTOT(IQUAN) * ERRTOT(JQUAN) 
     &                          * SAUX(NCSYS+IQUAN,NCSYS+JQUAN) 
          ENDIF
        ENDDO
      ENDDO
      AVEERR = SQRT(MAX(0,AVEERR))
      WRITE(LUNLOG, '(''All: Sum of Multi-Quantity Average = '','//
     &              'F14.7,'' +- '',F14.7)')
     &               AVESUM, AVEERR
      WRITE(LUNLOG,*)' '

* SwB begin
      PSUM = 0
      DO I=NCSYS+1,NPARA
        IF(CHPARA(I).EQ.'CHI2_N_SYM_PSUM') PSUM=INT(PARA(I))
      ENDDO
      DO ISUM=1,PSUM            ! next line assumes number of sums to compute <= 9
        WRITE(CHSUM,'("CHI2_N_SYM_P",I1.1)') ISUM
        NSUM=0
        DO I=NCSYS+1,NPARA
          IF (CHPARA(I).EQ.CHSUM) THEN
            NSUM=INT(PARA(I))
          ENDIF
        ENDDO
        DO IQUAN=1,NQUAN
          USEQUAN(IQUAN)=.FALSE.
          COEFF(IQUAN)=0.0D0
        ENDDO
        DO II=1,NSUM             ! next line assumes number of quantities <=99
          WRITE(CHSUM2,'("CHI2_N_SYM_P",I1.1,"_",I2.2)') ISUM, II
          DO I=NCSYS+1,NPARA
            IF (CHPARA(I).EQ.CHSUM2) THEN
              IQUAN=INT(PARA(I))
              USEQUAN(IQUAN)=.TRUE.
              COEFF(IQUAN)=EXCUP(I)
            ENDIF
          ENDDO
        ENDDO
*     
        AVESUM=0.0
        AVEERR=0.0
        DO IQUAN=1,NQUAN
          IF (USEQUAN(IQUAN)) THEN
            IVAR=NCSYS+IQUAN
            AVESUM = AVESUM - V (IVAR) * COEFF(IQUAN)
            AVEERR = AVEERR + ERRTOT(IQUAN)**2 * COEFF(IQUAN)**2
            DO JQUAN=1,NQUAN
              IF (USEQUAN(JQUAN)) THEN
                IF (IQUAN.LT.JQUAN) THEN
                  AVEERR = AVEERR + 2 * ERRTOT(IQUAN) * ERRTOT(JQUAN) *
     &    SAUX(NCSYS+IQUAN,NCSYS+JQUAN) * COEFF(IQUAN) * COEFF(JQUAN)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        AVEERR = SQRT(MAX(0,AVEERR))
        WRITE(LUNLOG, '('' P'','//
     &                'I1.1,'': Sum of Multi-Quantity Average = '','//
     &                'F14.7,'' +- '',F14.7)')
     &                 ISUM, AVESUM, AVEERR
        WRITE(LUNLOG,*)' '
*
      ENDDO
* SwB end
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
      INTEGER NT,MT,N,M,I,J,IJ,IL,NN,MM

      NN=MIN(30,N)
      MM=MIN(30,M)
      DO I=1,NN
       IJ=1
       IL=MM
       WRITE(*,1000) (A((J-1)*NT+I),J=IJ,IL)
      ENDDO

      IF (N.GT.30.OR.M.GT.30) THEN
      WRITE (*,*) ' '
      NN=MIN(30,N)
      MM=MIN(60,M)
      DO I=1,NN
       IJ=31
       IL=MM
       WRITE(*,1000) (A((J-1)*NT+I),J=IJ,IL)
      ENDDO

      WRITE (*,*) ' '
      NN=MIN(60,N)
      MM=MIN(30,M)
      DO I=31,NN
       IJ=1
       IL=MM
       WRITE(*,1000) (A((J-1)*NT+I),J=IJ,IL)
      ENDDO

      WRITE (*,*) ' '
      NN=MIN(60,N)
      MM=MIN(60,M)
      DO I=31,NN
       IJ=31
       IL=MM
       WRITE(*,1000) (A((J-1)*NT+I),J=IJ,IL)
      ENDDO
      ENDIF

      
1000  FORMAT (30(1X,G7.1))
*1000 FORMAT (1X,20E10.3)
*1000 FORMAT (1X,20F6.3)
      RETURN
      END
