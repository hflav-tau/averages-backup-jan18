************************************************************************
*
*     Template for a combination routine of the COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 2.00, February 20, 1997:
*     -- routine renamed from NEW_ROUTINE to NEW_COMBINATION
*     -- write output on logical unit LUNIT (no output if LUNIT is zero or less)
*
************************************************************************
*
      SUBROUTINE NEW_COMBINATION(CVAL,ERR2P,ERR2N,CL,IERR)
*     ====================================================
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
*                              or -1 if confidence level not available
      INTEGER IERR           ! error flag (non-zero if combination failed)
*
*     Note: this routine should only write output (including error messages)
*     ----  to logical unit LUNIT (variable LUNIT is available in a common block
*           described in master.inc); if LUNIT is zero or less, then this
*           routine should not write any output
*
********
*      
*     Add below code to perform the combination.
*     Set IERR to zero if algorithm is successful.
*     For the time being this is a dummy routine, so IERR=-1 is returned.
*
      IERR = -1
      IF(LUNIT.GT.0) WRITE(LUNIT,*)
     & 'Dummy routine NEW_COMBINATION called'
*
*     Return result 
*
      IF(IERR.EQ.0) THEN
*       CVAL  = ...
*       ERR2P = ...
*       ERR2N = ...
*       CL    = ...
      ELSE ! combination failed, return dummy values
        CVAL  = 0.D0
        ERR2P = 0.D0
        ERR2N = 0.D0
        CL    =-1.D0
      ENDIF
      END
