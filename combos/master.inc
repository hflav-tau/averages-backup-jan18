************************************************************************
*
*     "Master" include file for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 1.10, December 13, 1996:
*     -- symmetric uncertainties, symmetric statistical correlations
*        and the name of the measured quantity now made available
*     Version 1.20, December 13, 1996:
*     -- symmetric parameter excursions, and the name of the combined
*        analysis now made available 
*     -- the excursions of the parameter related to a systematic contribution
*        not associated to a parameter at the input file level are now set
*        to +-1 instead of 0
*     Version 1.21, January 17, 1996:
*     -- parameter MPARA increased from 20 to 30
*     Version 1.40, February 6, 1997:
*     -- parameter MCSYS increased from 20 to 30
*     Version 1.50, February 12, 1997:
*     -- parameter MCSYS increased from 30 to 50
*     Version 2.00, February 20, 1997:
*     -- added variables CHROUT, CHSTEP, ISTEP, VSTEP, and LUNIT
*     -- arrays in common COMINP reordered: order now matters for routine MASTER
*     Version 2.10, February 27, 1997:
*     -- new array of variables CHANAL containing the names of the individual
*        analyses in the combination
*     Version 2.30, May 12, 1997:
*    -- new variable NMEFF added (NMEFF = number of "effective" measurements = 
*        number of measurements which have not been de-weighted by COMBOS)
*     Version 2.31, Jun 19, 1997:
*     -- parameter MMEAS increased from 20 to 30
*     Version 2.97, July 21, 2000:
*     -- parameter MQUAN increased from 10 to 12
*     Version 3.00, January 29, 2002:
*     -- parameter MMEAS increased from 30 to 40
*     -- parameter MPARA increased from MCSYS+10 to MCSYS+20 (60 to 70)
*     -- parameter MCSYS increased from 50 to 80
*     Version 3.10, January 21, 2003 (E. Barberio):
*     -- parameter MMEAS increased from 40 to 85
*     -- parameter MCSYS increased from 80 to 100
*     Version 3.??, May 28, 2005:
*     -- add flags for error scaling
*
************************************************************************
*
*  Input for a combination routine:  
*  --------------------------------
*
*  CHROUT = name of combination routine; CHCOMB = name of combined analysis
*  CHMEAS = name of measured quantity  ; CHSTEP = name of step variable (or ' ')
*  ISTEP  = step number                ; VSTEP  = step value
*  LUNIT  = logical unit for printout
*
*  NMEAS = number of measurements to combine (I=1,NMEAS)
*  NMEFF = number of measurements which have not been de-weighted (< or = NMEAS)
*  NCSYS = number of correlated systematic uncertainties (J=1,NCSYS)
*  NPARA = total number of parameters (K=1,NPARA)
*DR
*  NQUAN = number of different measured quantities
*
*  For each measurement I=1,NMEAS:
*  CHANAL(I)  = name of analysis producing measurement I
*  MEAS(I)    = central value of measurement I
*  STATP(I)   = absolute value of the positive  statistical uncertainty 
*  STATN(I)   = absolute value of the negative  statistical uncertainty
*  STAT (I)   = absolute value of the symmetric statistical uncertainty
*  USYSP(I)   = absolute value of the positive  uncorrelated syst. uncertainty
*  USYSN(I)   = absolute value of the negative  uncorrelated syst. uncertainty
*  USYS (I)   = absolute value of the symmetric uncorrelated syst. uncertainty
*DR 
*  KQUAN(I) = index of the quantity which is measured
*
*  For each pair of measurements I1=1,NMEAS and I2=1,NMEAS:
*  STACOR(I1,I2) = statistical correlation between measurements I1 and I2
*
*  For each measurement I=1,NMEAS and each correlated systematic uncertainty J:
*  CSYSP(I,J) = change in central value of measurement I for positive excursion
*                 of the parameter associated to systematic uncertainty J
*  CSYSN(I,J) = change in central value of measurement I for negative excursion
*                 of the parameter associated to systematic uncertainty J
*  CSYS (I,J) = symmetric change in central value of measurement I for excursion
*                 of the parameter associated to systematic uncertainty J
*
*    Note: the signs of CSYSP(I,J), CSYSN(I,J) and CSYS(I,J) are significant;
*          CSYSP(I,J)=CSYSN(I,J)=CSYS(I,J)=0 if parameter J does not affect
*          measurement I.
*
*  For each parameter K=1,NPARA:
*  PARA(K)    = value of parameter K
*  EXCUP(K)   = absolute value of the positive excursion of parameter K
*  EXCUN(K)   = absolute value of the negative excursion of parameter K
*  EXCU (K)   = absolute value of the symmetric excursion of parameter K
*  CHPARA(K)  = name of parameter K
*
*    Notes: NPARA is always at least equal to NCSYS, and if K is smaller than
*           or equal to NCSYS, K refers to the parameter associated with the
*           correlated systematic uncertainty J=K.
*
*DR add MQUAN (warning MQUAN should equal MQUANT in combos.inc
      INTEGER MMEAS,MCSYS,MPARA,MQUAN
      PARAMETER(MMEAS=85,MCSYS=100,MPARA=MCSYS+20,MQUAN=12)
*DR add NQUAN,KQUAN
      INTEGER ISTEP,LUNIT,NMEAS,NMEFF,NCSYS,NPARA,NQUAN,KQUAN(MMEAS)
      DOUBLE PRECISION VSTEP,MEAS(MMEAS),
     & STAT(MMEAS),STATP(MMEAS),STATN(MMEAS),STACOR(MMEAS,MMEAS),
     & USYS(MMEAS),USYSP(MMEAS),USYSN(MMEAS),
     & CSYS(MMEAS,MCSYS),CSYSP(MMEAS,MCSYS),CSYSN(MMEAS,MCSYS),
     & PARA(MPARA),EXCU(MPARA),EXCUP(MPARA),EXCUN(MPARA)
      CHARACTER*80 CHCOMB,CHANAL(MMEAS)
*DR add chquan
      CHARACTER*16 CHPARA(MPARA),CHROUT,CHMEAS,CHSTEP,CHQUAN(MQUAN)
*
      COMMON/COMINP/ISTEP,LUNIT,MEAS,STAT,STATP,STATN,STACOR,
     &       USYS,USYSP,USYSN,CSYS,CSYSP,CSYSN,PARA,EXCU,EXCUP,EXCUN,
     &       CHPARA,CHROUT,CHCOMB,CHANAL,CHMEAS,CHSTEP,CHQUAN,
     &       VSTEP,NMEAS,NMEFF,NCSYS,NPARA,NQUAN,KQUAN
*
*     Important note: routine MASTER relies on the fact that arrays
*     --------------  USYS, USYSP, USYSN, and CSYS are stored next to each 
*                     other in the common block
*RVK add flags for error scaling
      DOUBLE PRECISION XSSTAT(MMEAS),XSUSYS(MMEAS),XSCSYS(MMEAS,MCSYS)
      COMMON/RVK/XSSTAT,XSUSYS,XSCSYS
