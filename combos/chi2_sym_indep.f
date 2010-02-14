************************************************************************
*
*     "CHI2_SYM_INDEP" combination routine for COMBOS program
*
*     Henry Seywerd and Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.00, December  6, 1996:
*     -- original release (was called CHI2_SYM)
*     Version 1.10, December 10, 1996:
*     -- renamed to CHI2_SYM_INDEP
*     -- local variables and subroutines renamed to avoid conflicts with 
*        new include file master.inc
*     -- bug fixed in chi2 calculation
*     Version 1.20, December 13, 1996:
*     -- protect PROB in case of zero degrees of freedom
*     Version 2.00, February 20, 1997:
*     -- write output on logical unit LUNIT (no output if LUNIT is zero or less)
*     Version 2.30, May 12, 1997:
*     -- use NMEFF instead of NMEAS in computation of chi2 probability
*
************************************************************************
*
*     ===================================================
      SUBROUTINE chi2_sym_indep(cval,err2p,err2n,cl,ierr)
*     ===================================================
*
*     Combination routine using simple chi2 fit
*
*     Olivier Schneider and Henry Seywerd, CERN/PPE-ALE
*     December 4, 1996
*
*     Input data: available from a common block in file master.inc
*     ----------
*
      IMPLICIT NONE
      INCLUDE 'master.inc'
*      
*     Output data: returned in the arguments 
*     -----------
*
      DOUBLE PRECISION cval  ! combined value (or average value)
      DOUBLE PRECISION err2p ! positive total uncertainty**2 on CVAL
      DOUBLE PRECISION err2n ! negative total uncertainty**2 on CVAL
      DOUBLE PRECISION cl    ! confidence level of combination (e.g. chi2 prob.)
      INTEGER ierr           ! error flag (non-zero if combination failed)
*
********
*      
C Wrapper routine for calculating the mean value of the parameters in
C the array values.
C Sets up the covariance matrix from the input parameters
C Calls the calculation routine.
C
C -- This version for using the chi**2 method of D Brown etc. 
C

      DOUBLE PRECISION V(0:mcsys,0:mcsys), c(0:mcsys)
      DOUBLE PRECISION state(0:mcsys), estate(0:mcsys)
C$$$  INTEGER imeas, icsys
      DOUBLE PRECISION chi2
C Symmetrized errors.
      DOUBLE PRECISION my_stat(mmeas),my_usys(mmeas),
     &                 my_csys(mmeas,mcsys)

      ierr = 0

C$$$        print *, 'nmeas, ncsys', nmeas, ncsys
C$$$        DO imeas=1,nmeas

C$$$        print *, meas(imeas)   
C$$$        print *, statp(imeas), statn(imeas)   
C$$$        print *, usysp(imeas), usysn(imeas)   

C$$$          DO icsys=1,ncsys
C$$$            print *, 'anal, icsys, sys', 
C$$$       &         imeas, icsys, csysp(imeas, icsys),csysn(imeas, icsys)
C$$$        ENDDO
C$$$        ENDDO
C Symmetrize the errors
      CALL my_symmetrize(my_usys,my_stat,my_csys)

C Build the covariance matrix from the tables of uncorrelated and 
C correlated errors.
      CALL corrcv(my_stat,my_usys,my_csys,v, c)
C$$$      PRINT *,'c = ',c
C$$$      PRINT *,'v = ',v

C Perform the fit and get the return values.
      CALL corrfit(ierr, v, c, state, estate)
C Errors at symmetric in this version.
      cval = state(0)
      err2p = estate(0)**2
      err2n = estate(0)**2

      CALL corrchisq(v, c, state, estate, my_stat, my_usys, my_csys,
     $   chi2, cl)
C$$$      print *,'debug: ',cval,dsqrt(err2p),dsqrt(err2n),chi2
      END

C==================================================================
      SUBROUTINE my_symmetrize(my_stat,my_usys,my_csys)

C Symmetrize the errors for later use.
      IMPLICIT NONE
      INCLUDE 'master.inc'
      DOUBLE PRECISION my_stat(mmeas),my_usys(mmeas),
     &                 my_csys(mmeas,mcsys)
      INTEGER imeas, icsys

      DO imeas=1,nmeas
C Statistics
        my_stat(imeas)  = SQRT(0.5*(statp(imeas)**2 + statn(imeas)**2))
**        my_stat(imeas)  = 0.5*(ABS(statp(imeas)) + ABS(statn(imeas)))
C Uncorr systematics
        my_usys(imeas) = SQRT(0.5*(usysp(imeas)**2 + usysn(imeas)**2))
**        my_usys(imeas) = 0.5*(ABS(usysp(imeas)) + ABS(usysn(imeas)))
C Correlated Systematics
        DO icsys = 1,ncsys
          my_csys(imeas,icsys) =  SQRT(0.5*(csysp(imeas,icsys)**2 +
     $       csysn(imeas,icsys)**2))
**          my_csys(imeas,icsys) =  0.5*(ABS(csysp(imeas,icsys)) +
**     $       ABS(csysn(imeas,icsys)))
        ENDDO
      ENDDO
      END
C==================================================================
      SUBROUTINE corrcv(my_stat,my_usys,my_csys,v, c)
C
C Constuct the matrix V of the errors on the parameters
C V is a k+1 x k+1 matrix with k the number of correlated systematic errors.
C The elements of V are:
C     V_lm:
C     V_00 = Sum_i 1/sig_i^2
C     V_ll = [Sum_i (sys_li/sig_i)^2] + 1
C     V_0l = V_l0 = Sum_i sys_li/sig_i^2
C     V_lm = V_ml = Sum_i (sys_li sys_mi)/sig_i^2
C
C Where: sig_i is the sum of the statistical and uncorrelated systematics
C for the i^th measurement
C      : sys_ij is the j^th correlated systematic for the i^it measurement
C
C Similarly construct the vector C_l of size k+1
C Where:
C     C_0 = Sum_i A_i/sig_i^2
C     C_l = Sum_i (A_i sys_li)/sig_i^2

C
      IMPLICIT NONE
      INCLUDE 'master.inc'

C Covariance Matrix and vector of values.
      DOUBLE PRECISION v(0:mcsys,0:mcsys)
      DOUBLE PRECISION c(0:mcsys)

      INTEGER imeas, icsys, jcsys
      DOUBLE PRECISION uncor, cori, corj
      DOUBLE PRECISION my_stat(mmeas),my_usys(mmeas),
     &                 my_csys(mmeas,mcsys)

C Mean value from weighted uncorrelated errors.
C$$$      REAL wmean

      CALL vzero(v(0,0),2*(mcsys+1)*(mcsys+1))
      CALL vzero(c(0),2*(mcsys+1))

C Build the matrix and vector
C For each analysis do:
      DO imeas=1,nmeas
C Uncorrelated statistical error + uncorr systematic
C Symmetrize the errors
        uncor = my_stat(imeas)**2 + my_usys(imeas)**2
        V(0,0) = V(0,0) + 1.0D0/uncor

        C(0) = C(0) + meas(imeas)/uncor

C Build covar in two loops over the systemtics 
        DO icsys = 1,ncsys
          cori = my_csys(imeas,icsys)
          cori=DSIGN(cori,csysp(imeas,icsys))
          V(icsys,icsys) = V(icsys,icsys) + cori**2/uncor
          DO jcsys = icsys+1,ncsys
            corj = my_csys(imeas,jcsys)
            corj=DSIGN(corj,csysp(imeas,jcsys))
            V(icsys,jcsys) = V(icsys,jcsys) + cori*corj/uncor
          ENDDO
          V(0,icsys) = V(0,icsys) + cori/uncor
          C(icsys) = C(icsys) + meas(imeas)*cori/uncor
        ENDDO
      ENDDO


c$$$C For interests sake calculate mean weighted by uncorrelated errors
c$$$      wmean = 0.
c$$$      DO imeas=1,nmeas
c$$$        uncor = (0.5D0*(statp(imeas) + statn(imeas)))**2
c$$$     &         +(0.5D0*(usysp(imeas) +  usysn(imeas)))**2
c$$$        print 10, imeas, SQRT(uncor), (1.0/uncor)/V(0,0)
c$$$ 10     FORMAT(' Uncorr error for anal ', I2, ' is ', F6.4,
c$$$     $     ' and wt. ', F6.4) 
c$$$        wmean = wmean + (1.0/uncor)/V(0,0)*meas(imeas)
c$$$      ENDDO
c$$$      print 20, wmean, sqrt(1.0/SNGL(V(0,0)))
c$$$ 20   FORMAT(' Mean from uncorrelated errors ', 2F6.4)


C Loop again to finalize the matrices
      DO icsys = 1,ncsys
        V(0,icsys) = V(0,icsys)
        V(icsys,0) = V(0,icsys)
        V(icsys,icsys) = V(icsys,icsys) + 1.0D0
        DO jcsys = 1,icsys-1
          V(icsys,jcsys) = V(jcsys,icsys)
        ENDDO
        C(icsys) = C(icsys)
      ENDDO

      END

C==================================================================
      SUBROUTINE corrfit(ierr, v, c, state, estate)
      IMPLICIT NONE
C Given the input covariance matrix V and the normalized input
C values in C, calculate a state vector: state and its
C error: estate

      INCLUDE 'master.inc'
      INTEGER ierr,LENOCC
      DOUBLE PRECISION V(0:mcsys,0:mcsys)
      DOUBLE PRECISION C(0:mcsys)
      DOUBLE PRECISION state(0:mcsys), estate(0:mcsys)

      INTEGER ifail
      INTEGER icsys, jcsys

C Invert V and multiply by C to get mean value of input.
C Error from SQRT of V^-1(i,i)

      ierr = 0   
      CALL DSINV(ncsys+1,v(0,0),mcsys+1,ifail)
      IF(IFAIL.NE.0) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &   ': FATAL ERROR in corrfit matrix inversion failed', ifail
        ierr = -1
        RETURN 
      ENDIF 

C Calculate the state vector, and its error. 
      CALL vzero(state(0), 2*(mcsys+1))
      CALL vzero(estate(0), 2*(mcsys+1))
      DO icsys=0,ncsys
        DO jcsys=0,ncsys
          state(icsys) = state(icsys) + V(icsys,jcsys)*c(jcsys)
        ENDDO
        estate(icsys) = SQRT(v(icsys,icsys))
      ENDDO

      END
C==================================================================
      SUBROUTINE corrchisq(v, c, state, estate, 
     &                     my_stat, my_usys, my_csys,
     &                     chi2, probchi2)

      IMPLICIT NONE
      INCLUDE 'master.inc'
C Covariance Matrix and vector of values.
      DOUBLE PRECISION v(0:mcsys,0:mcsys)
      DOUBLE PRECISION c(0:mcsys)
      DOUBLE PRECISION state(0:mcsys), estate(0:mcsys)
      DOUBLE PRECISION my_stat(mmeas),my_usys(mmeas),
     $                 my_csys(mmeas,mcsys)
      DOUBLE PRECISION chi2, probchi2

      INTEGER icsys, imeas
      REAL chisq, chidof, probchi
      REAL prob
      REAL xxx

      chisq = 0.

C Compute chi**2 from state vector + measured quantities and error.
      DO imeas=1,nmeas
        xxx = state(0) - meas(imeas)
        DO icsys=1,ncsys
*OS Bug fixed Dec 10, 1996
*OS       xxx = xxx + state(icsys)*my_csys(imeas,icsys)
          xxx = xxx + state(icsys)*
     &                DSIGN(csys(imeas,icsys),csysp(imeas,icsys))
        ENDDO
        chisq = chisq + xxx**2/(my_stat(imeas)**2+my_usys(imeas)**2)
      ENDDO
      DO icsys=1,ncsys
        chisq = chisq + state(icsys)**2
      ENDDO

      if(nmeff.gt.1) then
        chidof = chisq/REAL(nmeff-1)
        probchi = PROB(chisq,nmeff-1)
      else
        chidof = -1.
        probchi = -1.
      endif
C      PRINT *, 'Chisq, from state per dof, and prob: ',
C     $   chisq, chidof, probchi
C Convert to DP
      chi2 = DBLE(chisq) 
      probchi2 = DBLE(probchi)
      END
