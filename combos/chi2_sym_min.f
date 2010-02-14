************************************************************************
*
*     "CHI2_SYM_MIN" combination routine for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.20, December 13, 1996:
*     -- original release
*     Version 1.30, January 27, 1997:
*     -- extra argument SPECIAL_FLAG in FCN_CHI2_SYM_CIRC_INIT
*     Version 2.00, February 20, 1997:
*     -- first argument removed in call to MINUIT_CHI2
*     -- two first arguments removed in call to PREPARE_CHI2
*     -- two first arguments removed in entry FCN_CHI2_SYM_MIN
*     -- write output on logical unit LUNIT (no output if LUNIT is zero or less)
*     Version 2.30, May 12, 1997:
*     -- write PARAMETER logical lines with fit results
*
************************************************************************
*
      SUBROUTINE CHI2_SYM_MIN(CVAL,ERR2P,ERR2N,CL,IERR)
*     =================================================
*
*     Olivier Schneider, CERN/PPE-ALE
*     December 11, 1996
*
*     Combination routine performing a chi2 fit in MINUIT:
*       - only symmetric uncertainties are used
*       - can handle correlated statistical uncertainties with any
*          correlation coefficient
*       - can handle 100% or 0% systematic correlations
*
*     Note: this routine is equivalent to CHI2_SYM, except that it uses
*           MINUIT rather than linear algebra to solve the problem !
*
      IMPLICIT NONE
      DOUBLE PRECISION CVAL,ERR2P,ERR2N,CL
      INTEGER IERR
      EXTERNAL FCN_CHI2_SYM_MIN_INIT,FCN_CHI2_SYM_MIN
*
      CALL MINUIT_CHI2(FCN_CHI2_SYM_MIN_INIT,FCN_CHI2_SYM_MIN,
     &                 CVAL,ERR2P,ERR2N,CL,IERR)
      END
*
************************************************************************
*
      SUBROUTINE FCN_CHI2_SYM_MIN(NPAR,D,F,X,IFLAG,FUNC)
*     ==================================================
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
      INTEGER NPAR, IFLAG
*
*     Argument of entry FCN_CHI2_SYM_MIN_INIT
*
      INTEGER SPECIAL_FLAG,IERR
*
      INCLUDE 'master.inc'
*
*     External
* 
      INTEGER LENOCC
*
*     Local variables
*
      INTEGER I,J,IVARBL
      CHARACTER*10 CHNAM
      CHARACTER*16 CHNAME
      CHARACTER*1 COMM(2)
      DOUBLE PRECISION VAL,ERROR
*OS      DOUBLE PRECISION EPLUS,EMINUS,EPARAB,GLOBCC
      DOUBLE PRECISION CHI2,DUMMY
      INTEGER LCSYS
      SAVE LCSYS
      DATA LCSYS/0/
      DOUBLE PRECISION W(MMEAS,MMEAS),Z(MMEAS),V(MCSYS)
      SAVE             W
      LOGICAL NOT_READY
      SAVE    NOT_READY
      DATA    NOT_READY/.TRUE./
*
*     Check that routine was initialized
*
      IF(NOT_READY) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ' not initialized !'
        STOP 566
      ENDIF
*
*OS      IF(IFLAG.EQ.1) THEN ! initialization
*OS      ENDIF
*
*     Compute derivatives
*
      IF(IFLAG.EQ.2) THEN
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ': cannot compute derivatives'
        STOP 669
      ENDIF
*
*     Compute vector V = rescaled parameters
*
      IF(NPAR.NE.LCSYS) THEN 
        IF(LUNIT.GT.0) WRITE(LUNIT,*) CHROUT(:LENOCC(CHROUT)),
     &                 ': problem NPAR,LCSYS = ',NPAR,LCSYS
        STOP 543
      ENDIF
      DO I=1,NCSYS
        V(I)=(X(I)-PARA(I))/EXCU(I)
      ENDDO
*DR
*      V(LCSYS)=X(LCSYS)
      DO I=NCSYS+1,LCSYS
        V(I)=X(I)
      ENDDO
*DR end     
*
*     Compute Z = "Delta^T * V + X"
*
      CALL DVCPY(NMEAS,MEAS(1),MEAS(2),Z(1),Z(2))
      CALL DMMLA(NMEAS,LCSYS,1,CSYS(1,1),CSYS(1,2),CSYS(2,1), ! transposed
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
*     Function to be minimized
*
      F=CHI2
*
      IF(IFLAG.NE.3.OR.LUNIT.LE.0) RETURN
*
*     Print PARAMETERS logical lines for possible use in subsequent COMBOS jobs
*
      WRITE(LUNIT,6001)
 6001 FORMAT(/,' PARAMETER logical lines',
     &         ' for possible use in subsequent COMBOS jobs:',/)
 6002 FORMAT(A1,'PARAMETER ',A,3(2X,F10.5,SP),1X,
     &       A1,' CHI2_SYM_MIN output')
*OS 6004 FORMAT('*          EPLUS,EMINUS,EPARAB,GLOBCC=',SP,4(2X,F10.5))
      DO I=1,LCSYS
* DR        IF(I.LT.LCSYS) THEN 
*          CHNAME=CHPARA(I)
*        ELSE
*          CHNAME=CHMEAS
*        ENDIF
        IF(I.LE.NCSYS) THEN 
          CHNAME=CHPARA(I)
        ELSE
          CHNAME=CHQUAN(I-NCSYS)
        ENDIF
*DR end
        CALL MNPOUT(I,CHNAM,VAL,ERROR,DUMMY,DUMMY,IVARBL)
        IF(IVARBL.GE.0.AND.CHNAM(:5).NE.'lump ') THEN 
          COMM(1)=' '
        ELSE
          COMM(1)='*'
        ENDIF
        IF(IVARBL.GE.0.AND.CHNAM.EQ.CHNAME(:10)) THEN
          COMM(2)='!'
        ELSE
          COMM(2)='?' ! undefined or illdefined parameter (should not happen)
        ENDIF
        WRITE(LUNIT,6002) COMM(1),CHNAME,VAL,ERROR,-ERROR,COMM(2)
*OS        CALL MNERRS(I,EPLUS,EMINUS,EPARAB,GLOBCC)
*OS        WRITE(LUNIT,6004) EPLUS,EMINUS,EPARAB,GLOBCC
      ENDDO
      WRITE(LUNIT,'(1X)')
      RETURN
*
      ENTRY FCN_CHI2_SYM_MIN_INIT(SPECIAL_FLAG,IERR)
*     ==============================================
*
      SPECIAL_FLAG=0
      CALL PREPARE_CHI2(LCSYS,W,IERR)
      NOT_READY=IERR.NE.0
      END
