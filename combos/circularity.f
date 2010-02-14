************************************************************************
*
*     "Circularity" routine called by CHI2_SYM_CIRC combination routine
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.20, December 13, 1996:
*     -- original release
*     Version 1.21, January 17, 1997:
*     -- new argument INDX
*     -- handle new parameters XD, CHID and CHIS
*     Version 1.30, January 27, 1997:
*     -- add old parameters CHID4S and FSBR to allow external constraints on 
*        chid and fs; compute weighted averages for chid and fs
*     -- compute world averages for xd and dmd
*     Version 1.50, February 12, 1997:
*     -- add output argument CORR
*
************************************************************************
*
      SUBROUTINE CIRCULARITY(INDX,XIN,EXIN,XOUT,EXOUT,CORR)
*     =====================================================
*
*     Input:
*     ------
*     INDX(I) > 0 if input value available for parameter I
*     XIN(I) = value of input parameter I
*     EXIN = error matrix on XIN with 
*            EXIN(I,J)=0 if INDX(I)=0 or INDX(J)=0
*
*     Output:
*     -------
*     XOUT(I) = value of output parameter I
*     EXOUT = error matrix on XOUT
*     CORR  = correlation matrix between XIN and XOUT
*
      IMPLICIT NONE
*
      INCLUDE 'circularity.inc'
*
*     Arguments
*
      INTEGER INDX(NNEWPAR)
      DOUBLE PRECISION XIN(NNEWPAR),XOUT(NOLDPAR)
      DOUBLE PRECISION EXIN(NNEWPAR,NNEWPAR),EXOUT(NOLDPAR,NOLDPAR),
     &                 CORR(NOLDPAR,NNEWPAR)
*
*     Local variables
*
      INTEGER I
      DOUBLE PRECISION CHIS,FMESON,CHID1_XD,CHID_CHID4S,CHID_CHID1,
     &                 FMESON_FLAMB,FS1_FMESON,FS1_CHID,FS_FSBR,FS_FS1,
     &                 FD_FMESON,FD_FS,XDW_CHID,DMDW_XDW
      DOUBLE PRECISION DERIV(NOLDPAR,NNEWPAR)
*
*     Convention for partial derivatives:
*        YYY_XXX          = partial derivative of YYY with respect to XXX
*        DERIV(IYYY,IXXX) = partial derivative of YYY with respect to XXX
*
      CALL UZERO(DERIV,1,2*NOLDPAR*NNEWPAR) ! double precision
*
*     Compute xd from dmd and taubd
*     
      XOUT(IXD)=XIN(IDMD)*XIN(ITAUBD)
*
      DERIV(IXD,IDMD)=XIN(ITAUBD)
      DERIV(IXD,ITAUBD)=XIN(IDMD)
*
*     Compute chid from xd
*     
      XOUT(ICHID1)=0.5D0*XOUT(IXD)**2/(1.0D0+XOUT(IXD)**2)
      CHID1_XD=XOUT(IXD)/(1.D0+XOUT(IXD)**2)**2
*
      DERIV(ICHID1,IDMD)=CHID1_XD*DERIV(IXD,IDMD)
      DERIV(ICHID1,ITAUBD)=CHID1_XD*DERIV(IXD,ITAUBD)
*
*     Compute average value of chid from xd and Upsilon(4S) result
*
      IF(INDX(ICHID4S).GT.0) THEN
        CALL COMPUTE_ERROR(NNEWPAR,NOLDPAR,DERIV,EXIN,EXOUT,CORR)
        CALL WEIGHTED_AVERAGE(XIN(ICHID4S),EXIN(ICHID4S,ICHID4S),
     &                        XOUT(ICHID1),EXOUT(ICHID1,ICHID1),
     &                        CORR(ICHID1,ICHID4S),
     &                        XOUT(ICHID),CHID_CHID4S,CHID_CHID1)
      ELSE ! case where Upsilon(4S) measurement is not available
        XOUT(ICHID)=XOUT(ICHID1)
        CHID_CHID4S=0.D0
        CHID_CHID1=1.D0
      ENDIF
*
      DERIV(ICHID,IDMD)=CHID_CHID1*DERIV(ICHID1,IDMD)
      DERIV(ICHID,ITAUBD)=CHID_CHID1*DERIV(ICHID1,ITAUBD)
      DERIV(ICHID,ICHID4S)=CHID_CHID4S
*
*     Compute fs from average chid, chis, chibar and flamb
*
      IF(INDX(ICHIS).GT.0) THEN 
        CHIS=XIN(ICHIS)
      ELSE
        CHIS=0.5D0 ! default value of chi_s
      ENDIF
      FMESON=1.D0-XIN(IFLAMB)
      XOUT(IFS1)=(2.D0*XIN(ICHIB)-XOUT(ICHID)*FMESON)/
     &           (2.D0*     CHIS -XOUT(ICHID))
*
      FS1_FMESON=-XOUT(ICHID)/(2.D0*CHIS-XOUT(ICHID))
      FMESON_FLAMB=-1.D0
      FS1_CHID=(2.D0*XIN(ICHIB)-FMESON*2.D0*CHIS)/
     &         (2.D0*     CHIS -XOUT(ICHID))**2
      DERIV(IFS1,IDMD)=FS1_CHID*DERIV(ICHID,IDMD)
      DERIV(IFS1,ITAUBD)=FS1_CHID*DERIV(ICHID,ITAUBD)
      DERIV(IFS1,ICHID4S)=FS1_CHID*DERIV(ICHID,ICHID4S)
      DERIV(IFS1,ICHIS)=-2.D0*XOUT(IFS1)/(2.D0*CHIS-XOUT(ICHID))
      DERIV(IFS1,ICHIB)=2.D0/(2.D0*CHIS-XOUT(ICHID))
      DERIV(IFS1,IFLAMB)=FS1_FMESON*FMESON_FLAMB
*
*     Compute average value of fs from integrated mixing mixing measurements
*     and from branching ratio measurements
*
      IF(INDX(IFSBR).GT.0) THEN
        CALL COMPUTE_ERROR(NNEWPAR,NOLDPAR,DERIV,EXIN,EXOUT,CORR)
        CALL WEIGHTED_AVERAGE(XIN(IFSBR),EXIN(IFSBR,IFSBR),
     &                        XOUT(IFS1),EXOUT(IFS1,IFS1),
     &                        CORR(IFS1,IFSBR),
     &                        XOUT(IFS),FS_FSBR,FS_FS1)
      ELSE ! case where branching ratio measurement measurement is not available
        XOUT(IFS)=XOUT(IFS1)
        FS_FSBR=0.D0
        FS_FS1=1.D0
      ENDIF
*
      DERIV(IFS,IDMD)=FS_FS1*DERIV(IFS1,IDMD)
      DERIV(IFS,ITAUBD)=FS_FS1*DERIV(IFS1,ITAUBD)
      DERIV(IFS,ICHID4S)=FS_FS1*DERIV(IFS1,ICHID4S)
      DERIV(IFS,ICHIS)=FS_FS1*DERIV(IFS1,ICHIS)
      DERIV(IFS,ICHIB)=FS_FS1*DERIV(IFS1,ICHIB)
      DERIV(IFS,IFLAMB)=FS_FS1*DERIV(IFS1,IFLAMB)
      DERIV(IFS,IFSBR)=FS_FSBR
*
*     Compute fd from average fs and flamb
*
      XOUT(IFD)=0.5D0*(FMESON-XOUT(IFS))
*
      FD_FS=-0.5D0
      FD_FMESON=0.5D0+FD_FS*FS_FS1*FS1_FMESON
      DERIV(IFD,IDMD)=FD_FS*DERIV(IFS,IDMD)
      DERIV(IFD,ITAUBD)=FD_FS*DERIV(IFS,ITAUBD)
      DERIV(IFD,ICHID4S)=FD_FS*DERIV(IFS,ICHID4S)
      DERIV(IFD,ICHIS)=FD_FS*DERIV(IFS,ICHIS)
      DERIV(IFD,ICHIB)=FD_FS*DERIV(IFS,ICHIB)
      DERIV(IFD,IFLAMB)=FD_FMESON*FMESON_FLAMB
      DERIV(IFD,IFSBR)=FD_FS*DERIV(IFS,IFSBR)
*
*     Compute world average xd from chid
*
      XOUT(IXDW)=1.D0/DSQRT(0.5D0/XOUT(ICHID)-1.D0)
*
      XDW_CHID=1.D0/(XOUT(IXDW)*(1.D0-2.D0*XOUT(ICHID))**2)
      DO I=1,NNEWPAR
        DERIV(IXDW,I)=XDW_CHID*DERIV(ICHID,I)
      ENDDO
*
*     Compute world average dmd
*
      XOUT(IDMDW)=XOUT(IXDW)/XIN(ITAUBD)
      DMDW_XDW=1.D0/XIN(ITAUBD)
*
      DO I=1,NNEWPAR
        DERIV(IDMDW,I)=DMDW_XDW*DERIV(IXDW,I)
      ENDDO
      DERIV(IDMDW,ITAUBD)=DERIV(IDMDW,ITAUBD)-DMDW_XDW*XOUT(IDMDW)
*
*     Compute full error matrix
*
      CALL COMPUTE_ERROR(NNEWPAR,NOLDPAR,DERIV,EXIN,EXOUT,CORR)
      END
*
************************************************************************
*
      SUBROUTINE COMPUTE_ERROR(NIN,NOUT,DERIV,EXIN,EXOUT,CORR)
*     ========================================================
*
*     Input:  NIN   = number of input parameters
*     ------  NOUT  = number of output parameters
*             DERIV = matrix of partial derivatives:
*                       DERIV(I,J) = partial derivative of output parameter I
*                                    with respect to input parameter J
*             EXIN  = error matrix on input parameters
*
*     Output: EXOUT = error matrix on output parameters
*     ------- CORR  = correlations between output and input parameters:
*                       CORR(I,J) = correlation between output parameter I
*                                   and input parameter J
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER NIN,NOUT
      DOUBLE PRECISION EXIN(NIN,NIN),EXOUT(NOUT,NOUT)
      DOUBLE PRECISION DERIV(NOUT,NIN),CORR(NOUT,NIN)
*
*     Local variables
*
      DOUBLE PRECISION DUMMY
*
*     Compute CORR = EXIN * DERIV
*
      CALL DMMLT(NIN,NIN,NOUT,EXIN (1,1),EXIN (2,1),EXIN (1,2),
     &                        DERIV(1,1),DERIV(2,1),DERIV(1,2),
     &                        CORR (1,1),CORR (2,1),CORR (1,2),DUMMY)
*
*     Compute EXOUT = DERIV^T * EXIN * DERIV = DERIV^T * CORR
*
      CALL DMMLT(NOUT,NIN,NOUT,DERIV(1,1),DERIV(1,2),DERIV(2,1), ! transposed
     &                         CORR (1,1),CORR (2,1),CORR (1,2),
     &                         EXOUT(1,1),EXOUT(2,1),EXOUT(1,2),DUMMY)
      END
*
************************************************************************
*
      SUBROUTINE WEIGHTED_AVERAGE(X1,S2X1,X2,S2X2,CORR,AVG,W1,W2)
*     ================================================================
*
*     Compute the weighted average of two measurements
*
*     Input:  X1       = first measurement
*     ------  S2X1     = square of the error on X1
*             X2       = second measurement
*             S2X2     = square of the error on X2
*             CORR     = correlation term = sqrt(S2X1*S2X2)*rho
*
*     Output: AVG      = weighted average of X1 and X2
*     ------- W1       = weight of X1
*             W2       = weight of X2 = 1-W1
*
      IMPLICIT NONE
*
*     Arguments
*
      DOUBLE PRECISION X1,S2X1,X2,S2X2,CORR,AVG,W1,W2
*
*     Local variables
*
      INTEGER IER
      DOUBLE PRECISION W(2,2),SUMW
*
*     Build and invert error matrix
*
      W(1,1)=S2X1
      W(2,2)=S2X2
      W(1,2)=CORR
      W(2,1)=CORR
      CALL DSINV(2,W,2,IER) ! weight matrix
      IF(IER.NE.0) THEN
        PRINT *,'WEIGHTED_AVERAGE: cannot invert 2x2 matrix ! IER=',IER
        STOP 42
      ENDIF 
*
*     Compute weighted average
*
      SUMW=W(1,1)+W(1,2)+W(2,1)+W(2,2)
      W1  =(W(1,1)+0.5*(W(1,2)+W(2,1)))/SUMW
      W2  =(W(2,2)+0.5*(W(1,2)+W(2,1)))/SUMW
      AVG=W1*X1+W2*X2
      END
