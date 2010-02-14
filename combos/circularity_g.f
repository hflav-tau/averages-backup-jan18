************************************************************************
*
*     "Circularity" routine called by CHI2_SYM_CIRC_G combination routine,
*     using the "g" fractions rather than the "f" fractions in the 
*     definition of chibar.
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.20, April 3, 1997:
*     -- original release (cloned from circularity.f)
*
************************************************************************
*
      SUBROUTINE CIRCULARITY_G(INDX,XIN,EXIN,XOUT,EXOUT,CORR)
*     =======================================================
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
      INCLUDE 'circularity_g.inc'
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
*start new
      DOUBLE PRECISION R,TWO,RS,RL,GMESON
*end new
      DOUBLE PRECISION CHIS,CHID1_XD,CHID_CHID4S,CHID_CHID1,
*old     &                 FMESON,FMESON_FLAMB,
*old     &                 FS1_FMESON,FS1_CHID,FS_FSBR,FS_FS1,
*old     &                 FD_FMESON,FD_FS,XDW_CHID,DMDW_XDW
*start new
     &                 GMESON_RL,GMESON_FLAMB,FS1_GMESON,FS1_CHID,
     &                 FS_FSBR,FS_FS1,
     &                 FD_FS,FD_FLAMB,XDW_CHID,DMDW_XDW
*end new
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
*start new
      IF(INDX(IR).GT.0) THEN
        R=XIN(IR)
      ELSE
        R=1.D0 ! default value of r=TAUBU/TAUBD
      ENDIF
      TWO=1.D0+R
      IF(INDX(IRS).GT.0) THEN
        RS=XIN(IRS)
      ELSE
        RS=1.D0 ! default value of Rs=TAUBS/TAUB
      ENDIF
      IF(INDX(IRL).GT.0) THEN
        RL=XIN(IRL)
      ELSE
        RL=1.D0 ! default value of Rlambda=TAULAMB/TAUB
      ENDIF
      GMESON=1.D0-XIN(IFLAMB)*RL
      XOUT(IFS1)=(1.D0/RS)*(TWO*XIN(ICHIB)-XOUT(ICHID)*GMESON)/
     &                     (TWO*     CHIS -XOUT(ICHID))
*end new
*old      FMESON=1.D0-XIN(IFLAMB)
*old      XOUT(IFS1)=(2.D0*XIN(ICHIB)-XOUT(ICHID)*FMESON)/
*old     &           (2.D0*     CHIS -XOUT(ICHID))
*
*old      FS1_FMESON=-XOUT(ICHID)/(2.D0*CHIS-XOUT(ICHID))
*old      FMESON_FLAMB=-1.D0
*old      FS1_CHID=(2.D0*XIN(ICHIB)-FMESON*2.D0*CHIS)/
*old     &         (2.D0*     CHIS -XOUT(ICHID))**2
*start new
      FS1_GMESON=(1.D0/RS)*(-XOUT(ICHID))/(TWO*CHIS-XOUT(ICHID))
      GMESON_RL=-XIN(IFLAMB)
      GMESON_FLAMB=-RL
      FS1_CHID=(1.D0/RS)*(TWO*XIN(ICHIB)-GMESON*TWO*CHIS)/
     &                   (TWO*     CHIS -XOUT(ICHID))**2
*end new
      DERIV(IFS1,IDMD)=FS1_CHID*DERIV(ICHID,IDMD)
      DERIV(IFS1,ITAUBD)=FS1_CHID*DERIV(ICHID,ITAUBD)
      DERIV(IFS1,ICHID4S)=FS1_CHID*DERIV(ICHID,ICHID4S)
*old      DERIV(IFS1,ICHIS)=-2.D0*XOUT(IFS1)/(2.D0*CHIS-XOUT(ICHID))
*old      DERIV(IFS1,ICHIB)=2.D0/(2.D0*CHIS-XOUT(ICHID))
*old      DERIV(IFS1,IFLAMB)=FS1_FMESON*FMESON_FLAMB
*start new
      DERIV(IFS1,ICHIS)=-TWO*XOUT(IFS1)/(TWO*CHIS-XOUT(ICHID))
      DERIV(IFS1,ICHIB)=(1.D0/RS)*TWO/(TWO*CHIS-XOUT(ICHID))
      DERIV(IFS1,IR)=(XOUT(ICHID)/RS)*
     &               (GMESON*CHIS-XIN(ICHIB))/(TWO*CHIS-XOUT(ICHID))**2
      DERIV(IFS1,IRS)=-XOUT(IFS1)/RS
      DERIV(IFS1,IRL)=FS1_GMESON*GMESON_RL
      DERIV(IFS1,IFLAMB)=FS1_GMESON*GMESON_FLAMB
*end new
*
*     Compute average value of fs from integrated mixing measurements
*     and from branching ratio measurements
*
      IF(INDX(IFSBR).GT.0) THEN
        CALL COMPUTE_ERROR(NNEWPAR,NOLDPAR,DERIV,EXIN,EXOUT,CORR)
        CALL WEIGHTED_AVERAGE(XIN(IFSBR),EXIN(IFSBR,IFSBR),
     &                        XOUT(IFS1),EXOUT(IFS1,IFS1),
     &                        CORR(IFS1,IFSBR),
     &                        XOUT(IFS),FS_FSBR,FS_FS1)
      ELSE ! case where branching ratio measurement is not available
        XOUT(IFS)=XOUT(IFS1)
        FS_FSBR=0.D0
        FS_FS1=1.D0
      ENDIF
*
      DO I=1,NNEWPAR
        DERIV(IFS,I)=FS_FS1*DERIV(IFS1,I)
      ENDDO
      DERIV(IFS,IFSBR)=FS_FSBR
*
*     Compute fd from average fs and flamb
*
      XOUT(IFD)=0.5D0*(1.D0-XIN(IFLAMB)-XOUT(IFS))
*
*old      FD_FMESON=0.5D0+FD_FS*FS_FS1*FS1_FMESON
      FD_FS=-0.5D0
      FD_FLAMB=-0.5D0
      DO I=1,NNEWPAR
        DERIV(IFD,I)=FD_FS*DERIV(IFS,I)
      ENDDO
      DERIV(IFD,IFLAMB)=DERIV(IFD,IFLAMB)+FD_FLAMB
*
*     Compute world average xd from chid
*
      XOUT(IXDW)=1.D0/DSQRT(0.5D0/XOUT(ICHID)-1.D0)
*
      XDW_CHID=1.D0/(XOUT(IXDW)*(1.D0-2.D0*XOUT(ICHID))**2)
*     XDW_CHID=0.25D0*XOUT(IXDW)**3/XOUT(ICHID)**2
*      PRINT * ,'XDW_CHID=',XDW_CHID-0.25D0*XOUT(IXDW)**3/XOUT(ICHID)**2
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
