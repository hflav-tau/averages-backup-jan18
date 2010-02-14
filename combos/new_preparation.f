************************************************************************
*
*     Template for a preparation routine of the COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.52, July 12, 1999:
*     -- original release
*
************************************************************************
*
      SUBROUTINE NEW_PREPARATION(N,ICOMB)
*     ===================================
*
*     Input:  N      = number of analyses to combine
*     -----   ICOMB  = array of indexes to analyses
*                      ICOMB(I) = analysis number of Ith analysis to combine
*
*     Output: none
*     ------
*
*     Input data:  available from combos.inc and master.inc
*     -----------
*
*     Output data: written in master.inc
*     ------------
*
*     Note: this routine only writes output to logical unit LUNIT; 
*     ----- if LUNIT is zero or less, then this routine should not produce
*           any output
*
      IMPLICIT NONE
      INCLUDE 'master.inc'
*
*     Argument
*
      INTEGER N
      INTEGER ICOMB(N)
*
*     Local variable

*
      IF(LUNIT.GT.0) WRITE(LUNIT,*)
     & 'Dummy routine NEW_PREPARATION called'
      END
