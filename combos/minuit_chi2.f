************************************************************************
*
*     COMBOS combination routine to perform a NCSYS+1 parameter
*     chi2 fit using MINUIT, the (NCSYS+1)th parameter being the
*     combined value of the NMEAS measurements
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.20, December 13, 1996:
*     -- original release
*     Version 1.30, January 27, 1997:
*     -- extra argument SPECIAL_FLAG in FCN_INIT
*     -- allow for iterations on calls to FCN with a special argument
*     Version 2.00, February 20, 1997:
*     -- first argument of MINUIT_CHI2 removed (CHROUT now in master.inc)
*     -- parameter MINUIT_PRTLUN not used anymore (LUNIT now in master.inc)
*     -- default MINUIT printout level changed from -1 to +1 
*     -- improved printout (printout turned off if LUNIT is zero or less)
*     -- two first arguments of PREPARE_CHI2 and FCN_INIT removed
*         (CHROUT and LUNIT now in master.inc)
*     Version 2.30, May 12, 1997:
*     -- use NMEFF instead of NMEAS in computation of chi2 probability
*     Version 3.33, March 21, 2010:
*     -- hack (provided by Swagato) to handle the sum of 
*        two quantities in an average (see parameter CHI2_N_SYM_NSUM)
*
************************************************************************
*
      SUBROUTINE MINUIT_CHI2(FCN_INIT,FCN,CVAL,ERR2P,ERR2N,CL,IERR)
*     =============================================================
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
*
      IMPLICIT NONE
*
*     Arguments
*      
      EXTERNAL FCN_INIT,FCN
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR
*
*     Externals
*
      REAL PROB
      INTEGER LENOCC
*
*     Local variables
*
      CHARACTER*10 CHNAM
      INTEGER I,NPARI,NPARX,ISTAT,IVARBL
      DOUBLE PRECISION FMIN,FEDM,ERRDEF,DUMMY,PRTLEV,START,CVALO
      INTEGER MITER,NITER,SPECIAL_FLAG
      PARAMETER(MITER=50) ! maximum number of iterations on special FCN calls
*
      INCLUDE 'master.inc'
*
*     Default output arguments
*
      CVAL = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      CL   = -1.D0
      IERR = 9999
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

 
*DR allow NQUAN quantities to be fitted     DO I=1,NCSYS+1
      DO I=1,NCSYS+NQUAN
        IF(I.LE.NCSYS) THEN
          CALL MNPARM(I,CHPARA(I),PARA(I),0.1D0*EXCU(I),
     &                0.D0,0.D0,IERR)
        ELSE
*DR          CALL MNPARM(I,CHMEAS,START,0.02D0*START,0.D0,0.D0,IERR)
          CALL MNPARM(I,CHQUAN(I-NCSYS),START,0.02D0*START,
     &           0.D0,0.D0,IERR)
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
*     Get result of fit
*
      CALL MNPOUT(NCSYS+1,CHNAM,CVAL,ERR2P,DUMMY,DUMMY,IVARBL)
      ERR2P=ERR2P**2
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
        IF(NITER.EQ.MITER) THEN
          CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
          IF(IERR.NE.0) RETURN
          IF(LUNIT.GT.0) WRITE(LUNIT,*) ' '
        ENDIF
        IF(LUNIT.GT.0) WRITE(LUNIT,1003) CHNAM,CVAL,DSQRT(ERR2P)
 1003   FORMAT(1X,A,' = ',F14.10,' +-',F14.10)
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
      ERR2N=ERR2P ! return symmetric error
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(FMIN),NMEFF-1)) ! is this correct ?
      CALL EXECUTE(FCN,'SET PRINT',-1.D0,1,IERR,CHROUT,LUNIT)
      IF(IERR.NE.0) RETURN
      CALL EXECUTE(FCN,'CALL FCN',3.D0,1,IERR,CHROUT,LUNIT)
      END
*
************************************************************************
*
      SUBROUTINE EXECUTE(FCN,CHCOM,ARGLIS,NARG,IERR,CHROUT,LUNIT)
*     ===========================================================
*
*     Execute MINUIT command
*
      IMPLICIT NONE
*
*     Arguments
*
      EXTERNAL FCN
      CHARACTER*(*) CHCOM,CHROUT
      DOUBLE PRECISION ARGLIS(*)
      INTEGER NARG,IERR,LUNIT
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I
*
      CALL MNEXCM(FCN,CHCOM,ARGLIS,NARG,IERR,0)
      IF(IERR.NE.0.AND.LUNIT.GT.0) WRITE(LUNIT,1000)
     & CHROUT(:LENOCC(CHROUT)),IERR,CHCOM,(ARGLIS(I),I=1,NARG)
 1000 FORMAT(1X,A,': error',I6,' returned by ',A,10G10.4)
      END
*
************************************************************************
*
      SUBROUTINE PREPARE_CHI2(LCSYS,W,IERR)
*     =====================================
*
*     Prepare for chi2 fit.
*     Return number of unknown parameters and inverse of error matrix.
* DR Version 3.0 1 June 1999 modif: allow more than one param
      IMPLICIT NONE
      INCLUDE 'master.inc'
      INCLUDE 'combos.inc'
*
*     Arguments
*
      INTEGER LCSYS,IERR
      DOUBLE PRECISION W(MMEAS,MMEAS)
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER IADJ
      DOUBLE PRECISION ADJ(2)
      INTEGER INVOPT, ErrorFlag
      INTEGER I, J, II
      INTEGER NSUM, ISUM, IMEAS, IQUAN, ISUMOVER
      DOUBLE PRECISION    COEFF
      CHARACTER*13 CHSUM
      CHARACTER*16 CHSUM1
      CHARACTER*16 CHSUM2
      DOUBLE PRECISION WCOPY(MMEAS,MMEAS)
      DOUBLE PRECISION MEAS_SV(1000)
      LOGICAL IFIRST
      DATA IFIRST /.TRUE./
      SAVE MEAS_SV, IFIRST
*
*     Check that the number of correlated systematic uncertainties is 
*     less than the maximum
*
*DR      LCSYS=NCSYS+1
      LCSYS=NCSYS+NQUAN
      IF(LCSYS.GT.MCSYS) THEN
        IERR=MCSYS
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'parameter MCSYS too small',IERR
 1000   FORMAT(1X,A,': ',A,' ! IERR=',I6)
        RETURN
      ENDIF
*
*     Fill an additional row (row LCSYS) of matrix CSYS with -1's; 
*     this is against the philosophy that a combination routine should
*     not write into the common blocks listed in master.inc, but here it's OK
*     because:
*        - we write in a location which is not supposed to be read by other
*          combination routines
*        - routine MASTER knows about it, and will restore the overwritten 
*          stuff if necessary
*
*DR fill -1 and 0 at the end of the matrix
*      DO I=1,NMEAS
*        CSYS(I,LCSYS)=-1.D0
*      ENDDO
*DR
      IF (NQUAN.EQ.0) THEN
        CALL COMBOS_ERROR(-1,
     & 'MINUIT_CHI2:No quantities to be fit','cannot do anything')
      ENDIF
      DO I=1,NMEAS
        DO J=1,NQUAN
          IF (KQUAN(I).EQ.J) THEN
            CSYS(I,NCSYS+J)=-1.D0
          ELSE
            CSYS(I,NCSYS+J)=0.D0
          ENDIF
        ENDDO
      ENDDO
*DR end 

* SwB begin
      IF (IFIRST.EQV..TRUE.) THEN
        IFIRST = .FALSE.
        DO IMEAS = 1, NMEAS
          MEAS_SV(IMEAS)=MEAS(IMEAS)
        ENDDO
      ENDIF
*
      NSUM = 0
      DO I=1,NSPAR
        IF (CHSPAR(I).EQ.'CHI2_N_SYM_NSUM') NSUM=INT(SPAR1(I))
      ENDDO
*
      DO ISUM=1,NSUM            ! next line assumes number of linearized measurements <= 99
*
        WRITE (CHSUM,'("CHI2_N_SYM_",I2.2)') ISUM
        IMEAS=0
        ISUMOVER=0
        DO I=1,NSPAR
          IF (CHSPAR(I).EQ.CHSUM) THEN
            IMEAS=INT(SPAR1(I))
            ISUMOVER=INT(SPAR2(I))
          ENDIF
        ENDDO
*
        WRITE (CHSUM1,'("CHI2_N_SYM_",I2.2,"_AD")') ISUM
        IADJ=0
        ADJ(1)=0.0
        ADJ(2)=0.0
        DO I=1,NSPAR
          IF (CHSPAR(I).EQ.CHSUM1) THEN
            IADJ=1
            ADJ(1)=SPAR1(I)
            ADJ(2)=SPAR2(I)
          ENDIF
        ENDDO
* 
        IF (IADJ.NE.0) THEN
          WRITE (LUNLOG,
     &         '("ADJUST IMEAS = ",I3," MEAS = ",G14.7," to ",G14.7)')
     &         IMEAS, MEAS_SV(IMEAS), ADJ(1) + ADJ(2)*MEAS_SV(IMEAS)
          MEAS(IMEAS) = ADJ(1) + ADJ(2)*MEAS_SV(IMEAS)
        ENDIF
*        
        IF (IMEAS.GT.0.AND.ISUMOVER.EQ.0) THEN ! special flag for sum over all quantities
          WRITE (LUNLOG, '("OBSOLETE FEATURE")')
          STOP
*         DO IQUAN=1,NQUAN 
*            CSYS(IMEAS,NCSYS+IQUAN) = -1.D0
*         ENDDO
        ENDIF
*
        DO II=1,ISUMOVER        ! next line assumes number of quantities <=99
          WRITE(CHSUM2,'("CHI2_N_SYM_",I2.2,"_",I2.2)') ISUM, II
          DO I=1,NSPAR
            IF (IMEAS.GT.0.AND.CHSPAR(I).EQ.CHSUM2) THEN
              IQUAN=INT(SPAR1(I))
              COEFF=SPAR2(I)
              CSYS(IMEAS,NCSYS+IQUAN) = -1.D0*COEFF
            ENDIF
          ENDDO
        ENDDO 
*
      ENDDO     
*
      DO IMEAS=1,NMEAS
        WRITE (LUNLOG,'(''IMEAS,CSYS = '',I4,5(T20,20(1X,F6.2),/))')
     &          IMEAS,(CSYS(IMEAS,NCSYS+IQUAN),IQUAN=1,NQUAN)
      ENDDO
* SwB end

*
*     Build total error matrix
*
      DO I=1,NMEAS
        DO J=1,NMEAS
          W(I,J)=STACOR(I,J)*STAT(I)*STAT(J) ! statistical error matrix
        ENDDO
*DR XXXXX add correlated systematics here
        W(I,I)=W(I,I)+USYS(I)**2 ! add uncorrelated systematic error
      ENDDO


*
*     Invert total error matrix
*
      INVOPT=0 ! default option for matrix inversion using DSINV
      DO I=1,NSPAR
        IF(CHSPAR(I).EQ.'CHI2_N_SYM_INV') INVOPT=INT(SPAR1(I))
      ENDDO
      IF (INVOPT.EQ.0) THEN
        CALL MY_DSINV(NMEAS,W,MMEAS,IERR)
        PRINT *, 'MINUIT_CHI2: DSINV: W->W: IERR = ', IERR
      ELSE
        ErrorFlag=1
        DO I=1,NMEAS
          DO J=1,NMEAS
            if (ErrorFlag.eq.1.and.i.eq.j) print *, 
     &      'CHI2: i,j,w(i,j) = ',i,j,w(i,j),sqrt(max(0,w(i,j)))
            WCOPY(i,j)=W(i,j)
          ENDDO
        ENDDO
        CALL FindInv(WCOPY,W,NMEAS,MMEAS,ErrorFlag)
        PRINT *, 'MINUIT_CHI2: FindInv: W->W: ErrorFlag = ',ErrorFlag
        IERR=ErrorFlag
        ErrorFlag=1 ! <-
        DO I=1,NMEAS
          DO J=1,NMEAS
            if (ErrorFlag.eq.1.and.i.eq.j) print *, 
     &      'DONE: i,j,w(i,j) = ',i,j,w(i,j),sqrt(max(0,w(i,j)))
          ENDDO
        ENDDO
      ENDIF
*     
      IF(IERR.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,1000) CHROUT(:LENOCC(CHROUT)),
     &                 'cannot invert error matrix',IERR
        RETURN 
      ENDIF 
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
      
!     Subroutine to find the inverse of a square matrix
!     Author : Louisda16th a.k.a Ashwith J. Rego
!     Reference : Algorithm has been well explained in:
!     http://math.uww.edu/~mcfarlat/inverse.htm           
!     http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
      SUBROUTINE FindInv(matrix, inverse, n, d, errorflag)
      IMPLICIT NONE
!     Declarations
      INTEGER N,D
      INTEGER ErrorFlag !Return error status. -1 for error, 0 for normal
      DOUBLE PRECISION matrix(d,d) !Input matrix
      DOUBLE PRECISION inverse(d,d) !Inverted matrix
      INTEGER I,J,K,L
      DOUBLE PRECISION m
      DOUBLE PRECISION augmatrix(d,2*d) !augmented matrix
      LOGICAL FLAG
      DATA    FLAG /.TRUE./
      SAVE    FLAG

!     Augment input matrix with an identity matrix
      DO i = 1, n
        DO j = 1, 2*n
          IF (j .le. n ) THEN
            augmatrix(i,j) = matrix(i,j)
            if (ErrorFlag.eq.1.and.i.eq.j) 
     &      print *, 'Inv: i,j,matrix(i,j) = ',i,j,matrix(i,j)
          ELSE IF ((i+n) .eq. j) THEN
            augmatrix(i,j) = 1.d0
          Else
            augmatrix(i,j) = 0.d0
          ENDIF
        ENDDO
      ENDDO
      
!Reduce augmented matrix to upper traingular form
      DO k =1, n-1
        IF (augmatrix(k,k) .eq. 0.d0) THEN
          FLAG = .FALSE.
          DO i = k+1, n
            IF (augmatrix(i,k) .ne. 0) THEN
              DO j = 1,2*n
                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
              ENDDO
              FLAG = .TRUE.
              EXIT
            ENDIF
            IF (FLAG .EQV. .FALSE.) THEN
              PRINT*, "Matrix is non - invertible : 1 : check k = ",k,
     &                " augmatrix(k,k) = ",augmatrix(k,k),
     &                " maxtrix(k,k) = ",matrix(k,k)
              ErrorFlag = -1
              return
            ENDIF
          ENDDO
        ENDIF
        DO j = k+1, n			
          m = augmatrix(j,k)/augmatrix(k,k)
          DO i = k, 2*n
            augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
          ENDDO
        ENDDO
      ENDDO
      
!Test for invertibility
      DO i = 1, n
        IF (augmatrix(i,i) .eq. 0.d0) THEN
          PRINT*, "Matrix is non - invertible : 2 : check i = ",i
          ErrorFlag = -1
          return
        ENDIF
      ENDDO
      
!Make diagonal elements as 1
      DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)				
          augmatrix(i,j) = (augmatrix(i,j) / m)
        ENDDO
      ENDDO
      
!Reduced right side half of augmented matrix to identity matrix
      DO k = n-1, 1, -1
        DO i =1, k
          m = augmatrix(i,k+1)
          DO j = k, (2*n)
            augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
          ENDDO
        ENDDO
      ENDDO				
      
!store answer
      DO i =1, n
        DO j = 1, n
          inverse(i,j) = augmatrix(i,j+n)
        ENDDO
      ENDDO
      ErrorFlag = 0
      END ! SUBROUTINE FindInv
      
*
* $Id: dsinv.F,v 1.2 1999/09/08 08:05:11 mclareni Exp $
*
* $Log: dsinv.F,v $
* Revision 1.2  1999/09/08 08:05:11  mclareni
* A problem was reported in DSINV which failed on very small numbers, probably
* due to converting to single before a test. The conversion has been removed here
* and also in DSFACT. This resulted in mods to sfact.inc and sfactd.inc which
* meant that some other routines had to be tidied also.
*
* Revision 1.1.1.1  1996/02/15 17:49:05  mclareni
* Kernlib
*
*
*#include "kernnum/pilot.h"
          SUBROUTINE          MY_DSINV(N,A,IDIM,IFAIL)
          DOUBLE PRECISION    A(IDIM,*),  ZERO,  ONE,  X, Y
          CHARACTER*6         HNAME
          DOUBLE PRECISION    S1, S31, S32, S33,  DOTF
          DOTF(X,Y,S1)  =  X * Y + S1
          DATA      HNAME               /  'DSINV '  /
          DATA      ZERO, ONE           /  0.D0, 1.D0 /
          IF(IDIM .LT. N  .OR.  N .LE. 0)  GOTO 900
*#include "sfact.inc"
*
* $Id: sfact.inc,v 1.2 1999/09/08 08:05:21 mclareni Exp $
*
* $Log: sfact.inc,v $
* Revision 1.2  1999/09/08 08:05:21  mclareni
* A problem was reported in DSINV which failed on very small numbers, probably
* due to converting to single before a test. The conversion has been removed here
* and also in DSFACT. This resulted in mods to sfact.inc and sfactd.inc which
* meant that some other routines had to be tidied also.
*
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
*
* sfact.inc
*
          IFAIL  =  0
          DO 144    J  =  1, N
C             PRINT *, 'DSINV: J, A(J,J) = ',J,A(J,J)
             IF((A(J,J)) .LE. ZERO) GOTO 150
             A(J,J)  =  ONE / A(J,J)
             IF(J .EQ. N)  GOTO 199
 140         JP1  =  J+1
             DO 143   L  =  JP1, N
                A(J,L)  =  A(J,J)*A(L,J)
                S1      =  -A(L,J+1)
                DO 141  I  =  1, J
                   S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141               CONTINUE
                A(L,J+1)  =  -S1
 143            CONTINUE
 144         CONTINUE
 150      IFAIL  =  -1
          PRINT *, 'dsinv: J, A(J,J), IFAIL = ',J,A(J,J),IFAIL
          RETURN
 199      CONTINUE

*#include "sfinv.inc"
* $Id: sfinv.inc,v 1.1.1.1 1996/02/15 17:49:04 mclareni Exp $
*
* $Log: sfinv.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
*
* sfinv.inc
*
          IF(N .EQ. 1)  GOTO 399
          A(1,2)  =  -A(1,2)
          A(2,1)  =   A(1,2)*A(2,2)
          IF(N .EQ. 2)  GOTO 320
          DO 314    J  =  3, N
             JM2  =  J - 2
             DO 312 K  =  1, JM2
                S31  =  A(K,J)
                DO 311  I  =  K, JM2
                   S31  =  DOTF(A(K,I+1),A(I+1,J),S31)
 311               CONTINUE
                A(K,J)  =  -S31
                A(J,K)  =  -S31*A(J,J)
 312            CONTINUE
             A(J-1,J)  =  -A(J-1,J)
             A(J,J-1)  =   A(J-1,J)*A(J,J)
 314         CONTINUE
 320      J  =  1
 323         S33  =  A(J,J)
             IF(J .EQ. N)  GOTO 325
             JP1  =  J + 1
             DO 324 I  =  JP1, N
                S33  =  DOTF(A(J,I),A(I,J),S33)
 324            CONTINUE
 325         A(J,J)  =  S33
          JM1  =  J
          J    =  JP1
             DO 328 K  =  1, JM1
                S32  =  ZERO
                DO 327  I  =  J, N
                   S32  =  DOTF(A(K,I),A(I,J),S32)
 327               CONTINUE
                A(K,J)  =  S32
                A(J,K)  =  S32
 328            CONTINUE
          IF(J .LT. N)  GOTO 323
 399      CONTINUE

          RETURN
 900      CALL TMPRNT(HNAME,N,IDIM,0)
          RETURN
          END
