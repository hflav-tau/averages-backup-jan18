************************************************************************
*
*     "CHI2_SYM" combination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.10, December 10, 1996:
*     -- original release
*     Version 1.20, December 13, 1996:
*     -- error matrix construction and inversion, and other preparation
*        for chi2 fit are now performed in routine PREPARE_CHI2
*     -- sign convention for quantities "Y" has been changed (corrected)
*     -- protect PROB in case of zero degrees of freedom
*     Version 2.00, February 20, 1997:
*     -- two first arguments removed in call to PREPARE_CHI2
*     -- write output on logical unit LUNIT (no output if LUNIT is zero or less)
*     Version 2.30, May 12, 1997:
*     -- use NMEFF instead of NMEAS in computation of chi2 probability
*
************************************************************************
*
      SUBROUTINE CHI2_SYM(CVAL,ERR2P,ERR2N,CL,IERR)
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
      INTEGER LCSYS,I,J
      DOUBLE PRECISION V(MCSYS),W(MMEAS,MMEAS),
     &                 CHI2,DUMMY(MCSYS),Z(MMEAS),
     &                 TEMP(MMEAS,MCSYS),S(MCSYS,MCSYS)
*
*     Default output arguments
*
      CVAL = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      CL   = -1.D0
      IERR = 9999
      CALL PREPARE_CHI2(LCSYS,W,IERR)
      IF(IERR.NE.0) RETURN
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
      CALL DMMLA(LCSYS,NMEAS,LCSYS,TEMP(1,1),TEMP(2,1),TEMP(1,2),
     &                             CSYS(1,1),CSYS(1,2),CSYS(2,1), ! transposed
     &                             S   (1,1),S   (2,1),S   (1,2),DUMMY)
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
      CVAL=-V(LCSYS) ! note the sign !
      ERR2P=S(LCSYS,LCSYS)
      ERR2N=ERR2P ! return symmetric error
      IF(NMEFF.GT.1) CL=DBLE(PROB(SNGL(CHI2),NMEFF-1)) ! is this correct ?
      IF (NMEFF.GT.1) THEN
         PRINT *, 'CHI2_SYM: CHI2, NMEFF, CHI2/NDOF = ',
     &                       CHI2, NMEFF, CHI2/(NMEFF-1)
      ELSE
         PRINT *, 'CHI2_SYM: CHI2, NMEFF, CHI2/NDOF = ',
     &                       CHI2, NMEFF, ' set to 0.0 '
      ENDIF
      END
