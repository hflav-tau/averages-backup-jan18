************************************************************************
*
*     "Circularity" include file used by CHI2_SYM_CIRC_G combination algorithm
*     (routines CHI2_SYM_CIRC_G and CIRCULARITY_G).
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 2.20, April 3, 1997:
*     -- original release (cloned from circularity.inc)
*     Version 2.53, November 27, 1999:
*     -- add correlation between parameters FLAMB and FBS_FROM_BR; 
*        this correlation is taken from "parameter" CORR_FLAMB_FBS_
*        (which is assumed to be zero if not defined)
*
************************************************************************
*
*     Parameters to be replaced and expressed as a function of other parameters)
*     e.g. fs and fd=fu
*
      INTEGER NOLDPAR ! number of parameters to be replaced
      INTEGER IXD,ICHID1,ICHID,IFS1,IFS,IFD,IXDW,IDMDW!parameters to be replaced
      PARAMETER(IXD=1,ICHID1=2,ICHID=3,IFS1=4,IFS=5,IFD=6,
     &          IXDW=7,IDMDW=8,NOLDPAR=8)
      CHARACTER*16 CHOLDPAR(NOLDPAR)
      DATA CHOLDPAR(IXD   )/'XD = DMD * TAUBD'/
      DATA CHOLDPAR(ICHID1)/'1 / (2+2*XD**-2)'/
      DATA CHOLDPAR(ICHID )/'CHID            '/
      DATA CHOLDPAR(IFS1  )/'FBS (int mixing)'/
      DATA CHOLDPAR(IFS   )/'FBS             '/
      DATA CHOLDPAR(IFD   )/'FBD             '/
      DATA CHOLDPAR(IXDW  )/'XD  (extended)  '/
      DATA CHOLDPAR(IDMDW )/'DMD (extended)  '/
*
*     Parameters on which the replaced parameters depend
*
      INTEGER NNEWPAR ! number of param. on which the replaced param. depend
      INTEGER NNEWPAR1 ! NNEWPAR+1 to allow for CORR_FLAMB_FBS_
      INTEGER IDMD,ITAUBD,ICHID4S,ICHIS,ICHIB,
     &        IR,IRS,IRL,IFLAMB,IFSBR ! these parameters
      INTEGER ICORRF ! correlation coefficient between FLAMB and FBS_FROM_BR
      PARAMETER(IDMD=1,ITAUBD=2,ICHID4S=3,ICHIS=4,
     &          ICHIB=5,IR=6,IRS=7,IRL=8,IFLAMB=9,IFSBR=10,NNEWPAR=10,
     &          NNEWPAR1=NNEWPAR+1,ICORRF=NNEWPAR1)
      CHARACTER*16 CHNEWPAR(NNEWPAR1)
      LOGICAL NEEDED(NNEWPAR1)
      DATA CHNEWPAR(IDMD   ),NEEDED(IDMD   )/'DMD           ',.TRUE./
      DATA CHNEWPAR(ITAUBD ),NEEDED(ITAUBD )/'TAUBD         ',.TRUE./
      DATA CHNEWPAR(ICHID4S),NEEDED(ICHID4S)/'CHID_UPSILON4S',.FALSE./
      DATA CHNEWPAR(ICHIS  ),NEEDED(ICHIS  )/'CHIS          ',.FALSE./
      DATA CHNEWPAR(ICHIB  ),NEEDED(ICHIB  )/'CHIBAR        ',.TRUE./
      DATA CHNEWPAR(IR     ),NEEDED(IR     )/'TAUBU/TAUBD   ',.FALSE./
      DATA CHNEWPAR(IRS    ),NEEDED(IRS    )/'TAUBS/TAUB    ',.FALSE./
      DATA CHNEWPAR(IRL    ),NEEDED(IRL    )/'TAULAMB/TAUB  ',.FALSE./
      DATA CHNEWPAR(IFLAMB ),NEEDED(IFLAMB )/'FLAMB         ',.TRUE./
      DATA CHNEWPAR(IFSBR  ),NEEDED(IFSBR  )/'FBS_FROM_BR   ',.FALSE./
      DATA CHNEWPAR(ICORRF ),NEEDED(ICORRF) /'CORR_FLAMB_FBS_',.FALSE./
