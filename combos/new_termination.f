************************************************************************
*
*     Template for a termination routine of the COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.00, February 20, 1997:
*     -- original release
*
************************************************************************
*
      SUBROUTINE NEW_TERMINATION
     &           (NC,NS,STEP,CVAL,ERR2P,ERR2N,ERR2,CL,IERR)
*     =====================================================
*
*     Input:  NC     = number of cases
*     -----            NC=1: means only total error available
*                      NC=2: means also statistical error available
*             NS     = number of steps
*             STEP   = array of step values
*             CVAL   = array of combined values of the measured quantity
*             ERR2P  = array of positive errors**2 on CVAL
*             ERR2N  = array of negative errors**2 on CVAL
*             ERR2   = array of symmetric errors**2 on CVAL
*             CL     = array of confidence levels of CVAL
*             IERR   = array of error codes (0 means OK)
*
*     For step IS=1,NS:
*             CVAL (1,IS) = combined value at step IS (including syst.)
*             ERR2P(1,IS) = positive error**2 on CVAL(1,IS)
*             ERR2N(1,IS) = negative error**2 on CVAL(1,IS)
*             ERR2 (1,IS) = symmetric error**2 on CVAL(1,IS)
*             CL   (1,IS) = confidence level of fit performed to get CVAL(1,IS)
*                           (or -1 if not available)
*
*     Other input (available from common block described in master.inc):
*     -----------
*              
*             CHROUT = name of combination routine 
*                       used to obtain CVAL,ERR2P,ERR2N,ERR2,CL and IERR
*             CHCOMB = name of combined analysis
*             CHMEAS = name of measured quantity
*             CHSTEP = name of step variable
*             LUNIT  = logical unit for printout
*                       produced by this termination routine
*
*     Output:  none
*     -------
*
*     Note: this routine only write output to logical unit LUNIT; 
*     ----- if LUNIT is zero or less, then this routine should not produce
*           any output
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER MC,NC,NS
      PARAMETER(MC=2)
      DOUBLE PRECISION STEP(NS),CVAL(MC,NS),CL(MC,NS),
     &                 ERR2P(MC,NS),ERR2N(MC,NS),ERR2(MC,NS)
      INTEGER IERR(MC,NS)
*
      INCLUDE 'master.inc'
*
      IF(LUNIT.GT.0) WRITE(LUNIT,*)
     & 'Dummy routine NEW_TERMINATION called'
      END
