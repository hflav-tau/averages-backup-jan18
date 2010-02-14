************************************************************************
*
*     "Writer" routines for COMBOS program
*
*     Olivier Schneider, CERN/PPE-ALE
*
*     Version 1.00, December  6, 1996:
*     -- original release
*     Version 1.10, December 10, 1996:
*     -- added 4 keywords: SUPERSEDE, STAT_CORR_WITH, NOSYSTEMATICS, SYMMETRIZE
*     -- changed characteristics of NOCORRELATIONS keyword
*     Version 1.30, January 27, 1997:
*     -- added new keyword: SYNONYMS
*     Version 1.40, February 6, 1997:
*     -- added new keyword: LUMP
*     Version 1.50, February 12, 1997:
*     -- new format of STEP logical line for combined analyses
*     Version 2.00, February 20, 1997:
*     -- new format for MEASUREMENTS and CALL logical lines
*         (termination routines can now be specified on CALL logical lines)
*     -- added new keyword: SWITCH
*     -- removed 4 obsolete keywords: NOADJUSTMENTS, NOCORRELATIONS,
*                                     NOSYSTEMATICS, SYMMETRIZE
*     Version 2.30, May 12, 1997:
*     -- new routine WRITE_KUMAC to write kumac file
*     Version 2.50, March 9, 1998:
*     -- added new keyword: DISPLAY
*     -- print display names and values in KUMAC file (only when no STEP name)
*     Version 2.52, July 12, 1999:
*     -- introduction of the concept of preparation routines
*     Versions 3.25, May 20, 2005:
*     -- correct small bug in the computation of SYTOT in WRITE_KUMAC
*
************************************************************************
*
      SUBROUTINE WRITE_CARDS(IA)
*     ==========================
*
*     Routine called for the "END" keyword after the input analyses 
*     have been checked and after the combination has been performed:
*     write a "cards" file, i.e. a file similar to the input of the job
*
      IMPLICIT NONE
*
*     Argument
*
      INTEGER IA ! analysis number
*
*     Externals
*  
      INTEGER LENOCC
*
*     Local variable
*
      INTEGER MCOLS,NITEM,NCOLS
      PARAMETER(MCOLS=132) ! maximum number of columns in output file
      PARAMETER(NITEM=(MCOLS-1)/(1+16))! number of items per line in output file
*         WARNING: NITEM should be at least 2
      PARAMETER(NCOLS=1+NITEM*(1+16))  ! number of columns in output file
      CHARACTER*1 COMMENT
      CHARACTER*16 CH16(NITEM)
      INTEGER IO(3),K,I,IC,ICSTART,J,IS,NCOMB,ITEM,FCALL
      INTEGER IANAL,NL,ICONT,JCONT,ILUMP,NS
      DATA IO/0,+1,-1/
      INTEGER KEXPER,KMETHO,KSTATU
      LOGICAL ALREADY
*
      INCLUDE 'combos.inc'
      INTEGER IL(MCONT),ICOMB(MANAL)
*
*     BEGIN logical lines
*
      WRITE(LUNOUT,'(1X)') 
      WRITE(LUNOUT,1001) 
     & 'BEGIN',CHEXP(KEXP(IA)),CHMETH(KMETH(IA)),CHSTAT(KSTAT(IA)),
     & CHREF(IA)(:MAX0(1,LENOCC(CHREF(IA))))
 1001 FORMAT(A,T17,4(1X,A))
*
      IF(XCOMB(0)) THEN ! combination
        NS=NDEST
*
*     COMBINE logical lines
*
        WRITE(LUNOUT,'(1X)') 
        NCOMB=0
        DO IANAL=1,NANAL-1
          IF(XCOMB(IANAL)) THEN
            NCOMB=NCOMB+1
            ICOMB(NCOMB)=IANAL
            DO 77 ICONT=FCONT,NCONT(IANAL)
              IF(KLUMP(ICONT,IANAL).NE.KCONT(ICONT,IANAL)) THEN
                IF(NCOMB.GT.0) WRITE(LUNOUT,1012) 'COMBINE',
     &           (CHEXP(KEXP(ICOMB(I))),CHMETH(KMETH(ICOMB(I))),
     &            CHSTAT(KSTAT(ICOMB(I))),I=1,NCOMB)
 1012           FORMAT(A,(T17,3(1X,A)))
                NCOMB=0
                ILUMP=KLUMP(ICONT,IANAL)
                NL=0
                DO JCONT=FCONT,NCONT(IANAL)
                  IF(KLUMP(JCONT,IANAL).EQ.ILUMP) THEN
                    IF(JCONT.LT.ICONT) GOTO 77
                    NL=NL+1
                    IL(NL)=JCONT
                  ENDIF
                ENDDO
                WRITE(LUNOUT,1112) CHNAM(ILUMP),
     &           (CHNAM(KCONT(IL(I),IANAL)),I=1,NL)
 1112           FORMAT(T18,A,(T34,5(1X,A)))
              ENDIF
   77       ENDDO
          ENDIF
        ENDDO
        IF(NCOMB.GT.0) WRITE(LUNOUT,1012) 'COMBINE',
     &   (CHEXP(KEXP(ICOMB(I))),CHMETH(KMETH(ICOMB(I))),
     &    CHSTAT(KSTAT(ICOMB(I))),I=1,NCOMB)
      ELSE
*
*     SUPERCEDE and STAT_CORR_WITH logical lines
*
        NS=0
        DO I=1,NLINK(IA)
          IF(RHOSTA(I,IA).GT.1.) THEN
            CH16(1)='SUPERSEDE'
            CH16(2)=' '
          ELSE
            CH16(1)='STAT_CORR_WITH'
            WRITE(CH16(2),'(F16.6)',IOSTAT=J) RHOSTA(I,IA)
          ENDIF
          KEXPER=    LINK(I,IA)         /1000000
          KMETHO=MOD(LINK(I,IA),1000000)/1000
          KSTATU=MOD(LINK(I,IA),1000)
          WRITE(LUNOUT,1001) CH16(1),
     &     CHEXP(KEXPER),CHMETH(KMETHO),CHSTAT(KSTATU),CH16(2)
        ENDDO
      ENDIF
*
*     MEASUREMENTS logical lines
*
      WRITE(LUNOUT,'(1X)') 
      NL=0
      DO K=1,MIN0(FCONT-1,NCONT(NANAL))
        IF(KCONT(K,IA).NE.0) NL=K
      ENDDO
      WRITE(LUNOUT,1001) 'MEASUREMENT',(CHNAM(IABS(KCONT(K,IA))),K=0,NL)
*OSOS start
      DO I=2,NQUANT
        WRITE(LUNOUT,1001) 'MEASUREMENT',CHNAM(IABS(KCONT(1-I,IA)))
      ENDDO
*OSOS end
*
*     DISPLAY logical lines
*
      IF(NDISP(NANAL).GE.0) THEN
        IF(KDISP(0,NANAL).EQ.-1) THEN
          WRITE(LUNOUT,1001) 'DISPLAY',SUPERSEDED
        ELSE
          NL=0
          DO K=0,NDISP(NANAL)
            IF(KDISP(K,IA).GT.0) NL=K
          ENDDO
          WRITE(LUNOUT,1001) 'DISPLAY',(CHNAM(IABS(KDISP(K,IA))),K=0,NL)
        ENDIF
      ENDIF
*
*     STEPS logical lines
*
      IF(KSTEP(IA).GT.0) WRITE(LUNOUT,1011) 'STEP',CHNAM(KSTEP(IA)),
     &                                      (DEST(I),I=1,NS)
 1011 FORMAT(A,T18,A,(T34,8F10.4))
*
*     PARAMETERS logical lines
*
      WRITE(LUNOUT,1002) 
     & (CHNAM(KPAR(K,IA)),(PAR(IO(I),K,IA),I=1,3),K=1,NPAR(IA))
 1002 FORMAT(:,/,'PARAMETERS',(T18,A,1X,SS,F16.6,SP,2(1X,F16.6)))
*
*     SWITCH logical lines
*
      IF(XCOMB(0)) THEN ! combination
        CALL PRINT_SWIT
        DO I=1,NNAM
          ITEM=2
          CH16(1)='SYNONYMS'
          CH16(2)=CHNAM(I)
          DO J=1,NNAM
            IF(J.NE.I.AND.KSYN(J).EQ.I) THEN
              IF(ITEM.EQ.NITEM) THEN
                WRITE(LUNOUT,1004) CH16
                ITEM=2
                CH16(1)=' '
                CH16(2)=' '
              ENDIF
              ITEM=ITEM+1
              CH16(ITEM)=CHNAM(J)
            ENDIF
          ENDDO
          IF(ITEM.GT.2) WRITE(LUNOUT,1004) (CH16(J),J=1,ITEM)
        ENDDO
*
*       CALL logical lines
*
        WRITE(LUNOUT,'(1X)')
        IF(NCALLP.GT.0) THEN
          CH16(1)='CALL'
          FCALL=1
          DO I=1,NCALLP
            IF(LCALLP(I).NE.LCALLP(FCALL)) THEN
              WRITE(CH16(2),'(I2)') LCALLP(FCALL)
              WRITE(LUNOUT,1013)
     &         CH16(1),(CHCALP(KCALLP(J)),J=FCALL,I-1),CH16(2)
 1013         FORMAT(A,(T17,6(1X,A)))
              CH16(1)=' '
              FCALL=I
            ENDIF
          ENDDO
          WRITE(CH16(2),'(I2)') LCALLP(FCALL)
          WRITE(LUNOUT,1013)
     &     CH16(1),(CHCALP(KCALLP(J)),J=FCALL,NCALLP),CH16(2)
        ENDIF
        IF(NCALLC.GT.0) THEN
          CH16(1)='CALL'
          FCALL=1
          DO I=1,NCALLC
            IF(LCALLC(I).NE.LCALLC(FCALL)) THEN
              WRITE(CH16(2),'(I2)') LCALLC(FCALL)
              WRITE(LUNOUT,1013)
     &         CH16(1),(CHCALC(KCALLC(J)),J=FCALL,I-1),CH16(2)
              CH16(1)=' '
              FCALL=I
            ENDIF
          ENDDO
          WRITE(CH16(2),'(I2)') LCALLC(FCALL)
          WRITE(LUNOUT,1013)
     &     CH16(1),(CHCALC(KCALLC(J)),J=FCALL,NCALLC),CH16(2)
        ENDIF
        IF(NCALLT.GT.0) THEN
          CH16(1)='CALL'
          FCALL=1
          DO I=1,NCALLT
            IF(LCALLT(I).NE.LCALLT(FCALL)) THEN
              WRITE(CH16(2),'(I2)') LCALLT(FCALL)
              WRITE(LUNOUT,1013) 
     &         CH16(1),(CHCALT(KCALLT(J)),J=FCALL,I-1),CH16(2)
              CH16(1)=' '
              FCALL=I
            ENDIF
          ENDDO
          WRITE(CH16(2),'(I2)') LCALLT(FCALL)
          WRITE(LUNOUT,1013)
     &     CH16(1),(CHCALT(KCALLT(J)),J=FCALL,NCALLT),CH16(2)
        ENDIF
*
*     Comments or DATA logical lines
*
        COMMENT='*'
        WRITE(LUNOUT,1003) '* Combined results from routine '//
     &        CHCALC(KCALLC(1))(:LENOCC(CHCALC(KCALLC(1))))//':'
      ELSE
        COMMENT=' '
        WRITE(LUNOUT,1003) 'DATA'
*     &   'DATA ! (measurements to be used in combinations)'
      ENDIF
 1003 FORMAT(/,A)
      IC=-1
    1 IF(NCONT(IA).GT.IC) THEN
        IF(IC.GT.0) WRITE(LUNOUT,'(1X)')
        ITEM=0
        IF(KSTEP(IA).GT.0) THEN
          ITEM=ITEM+1
          CH16(ITEM)=CHNAM(KSTEP(IA))
        ENDIF
        ICSTART=IC
    2   IC=IC+1
        IF(KCONT(IC,IA).LE.0) THEN
        ELSE IF(XSYMC(IC,IA)) THEN
          ITEM=ITEM+1
          CH16(ITEM)=CHNAM(KCONT(IC,IA))(:LENOCC(
     &               CHNAM(KCONT(IC,IA))))
        ELSE IF(ITEM+1.EQ.NITEM) THEN 
          ITEM=ITEM+1
          CH16(ITEM)=' '
          IC=IC-1
        ELSE
          ITEM=ITEM+1
          CH16(ITEM)=CHNAM(KCONT(IC,IA))(:LENOCC(
     &               CHNAM(KCONT(IC,IA))))//'+'
          ITEM=ITEM+1
          CH16(ITEM)=CHNAM(KCONT(IC,IA))(:LENOCC(
     &               CHNAM(KCONT(IC,IA))))//'-'
        ENDIF
        IF(ITEM.LT.NITEM.AND.IC.LT.NCONT(IA)) GOTO 2
        DO I=1,ITEM
          CALL CRIGHT(CH16(I),1,LEN(CH16(I)))
        ENDDO
        WRITE(LUNOUT,1004) COMMENT,(CH16(I),I=1,ITEM)
        DO IS=1,NSTEP(IA)
          ITEM=0
          IF(KSTEP(IA).GT.0) THEN
            ITEM=ITEM+1
            WRITE(CH16(ITEM),'(F16.6)',IOSTAT=J) STEP(IS,IA)
          ENDIF
          IC=ICSTART
    3     IC=IC+1
          IF(KCONT(IC,IA).LE.0) THEN
          ELSE IF(XSYMC(IC,IA)) THEN
            ITEM=ITEM+1
            IF(IC.LE.FCONT-1) THEN
              WRITE(CH16(ITEM),'(   F16.6)',IOSTAT=J) CONT(IS,IC,IA)
            ELSE
              WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) CONT(IS,IC,IA)
            ENDIF
          ELSE IF(ITEM+1.EQ.NITEM) THEN 
            ITEM=ITEM+1
            CH16(ITEM)=' '
            IC=IC-1
          ELSE
            ITEM=ITEM+1
            WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) CONT(IS,+IC,IA)
            ITEM=ITEM+1
            WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) CONT(IS,-IC,IA)
          ENDIF
          IF(ITEM.LT.NITEM.AND.IC.LT.NCONT(IA)) GOTO 3
          WRITE(LUNOUT,1004) COMMENT,(CH16(I),I=1,ITEM)
        ENDDO
        GOTO 1
      ENDIF
 1004 FORMAT(A,7(1X,A16  ))
*
*     Extra DATA logical lines (for display data)
*
      IF(NDISP(IA).GE.0.AND.KDISP(0,IA).NE.-1) THEN
        WRITE(LUNOUT,1003)
     &   'DATA ! (for display only, not to be used in combinations)'
        IC=-1
   11   IF(NDISP(IA).GT.IC) THEN
          IF(IC.GT.0) WRITE(LUNOUT,'(1X)')
          ITEM=0
          IF(KSTEP(IA).GT.0) THEN
            ITEM=ITEM+1
            CH16(ITEM)=CHNAM(KSTEP(IA))
          ENDIF
          ICSTART=IC
   12     IC=IC+1
*         Check that this data was not already output
          ALREADY=.FALSE.
          DO ICONT=0,NCONT(IA)
            ALREADY=ALREADY.OR.KCONT(ICONT,IA).EQ.KDISP(IC,IA)
          ENDDO
          IF(KDISP(IC,IA).LE.0.OR.ALREADY) THEN
          ELSE IF(XSYMD(IC,IA)) THEN
            ITEM=ITEM+1
            CH16(ITEM)=CHNAM(KDISP(IC,IA))(:LENOCC(
     &                 CHNAM(KDISP(IC,IA))))
          ELSE IF(ITEM+1.EQ.NITEM) THEN 
            ITEM=ITEM+1
            CH16(ITEM)=' '
            IC=IC-1
          ELSE
            ITEM=ITEM+1
            CH16(ITEM)=CHNAM(KDISP(IC,IA))(:LENOCC(
     &                 CHNAM(KDISP(IC,IA))))//'+'
            ITEM=ITEM+1
            CH16(ITEM)=CHNAM(KDISP(IC,IA))(:LENOCC(
     &                 CHNAM(KDISP(IC,IA))))//'-'
          ENDIF
          IF(ITEM.LT.NITEM.AND.IC.LT.NDISP(IA)) GOTO 12
          DO I=1,ITEM
            CALL CRIGHT(CH16(I),1,LEN(CH16(I)))
          ENDDO
          WRITE(LUNOUT,1004) COMMENT,(CH16(I),I=1,ITEM)
          DO IS=1,NSTEP(IA)
            ITEM=0
            IF(KSTEP(IA).GT.0) THEN
              ITEM=ITEM+1
              WRITE(CH16(ITEM),'(F16.6)',IOSTAT=J) STEP(IS,IA)
            ENDIF
            IC=ICSTART
   13       IC=IC+1
*           Check that this data was not already output
            ALREADY=.FALSE.
            DO ICONT=0,NCONT(IA)
              ALREADY=ALREADY.OR.KCONT(ICONT,IA).EQ.KDISP(IC,IA)
            ENDDO
            IF(KDISP(IC,IA).LE.0.OR.ALREADY) THEN
            ELSE IF(XSYMD(IC,IA)) THEN
              ITEM=ITEM+1
              IF(IC.LE.MDISP) THEN
                WRITE(CH16(ITEM),'(   F16.6)',IOSTAT=J) DISP(IS,IC,IA)
              ELSE
                WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) DISP(IS,IC,IA)
              ENDIF
            ELSE IF(ITEM+1.EQ.NITEM) THEN 
              ITEM=ITEM+1
              CH16(ITEM)=' '
              IC=IC-1
            ELSE
              ITEM=ITEM+1
              WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) DISP(IS,+IC,IA)
              ITEM=ITEM+1
              WRITE(CH16(ITEM),'(SP,F16.6)',IOSTAT=J) DISP(IS,-IC,IA)
            ENDIF
            IF(ITEM.LT.NITEM.AND.IC.LT.NDISP(IA)) GOTO 13
            WRITE(LUNOUT,1004) COMMENT,(CH16(I),I=1,ITEM)
          ENDDO
          GOTO 11
        ENDIF
      ENDIF
*
*     END logical lines
*
      WRITE(LUNOUT,1005) ('*',I=1,NCOLS)
 1005 FORMAT('END',/,/,512A)
      END
*
************************************************************************
*
      SUBROUTINE WRITE_KUMAC(N,ICOMB,
     &                       MS,NC,CVAL,ERR2P,ERR2N,ERR2,CL,IERR,NS)
*     ==============================================================
*
*     This routine write a kumac file on logical unit LUNKUM
*
*     Input:  N       = number of analysis in combination
*     -----   ICOMB   = array of indexes to analyses
*                      ICOMB(I) = analysis number of Ith analysis to combine
*             MS      = 1st dimension of arrays CVAL,ERR2P,ERR2N,ERR2,CL and IERR
*             NC      = 2nd dimension of arrays CVAL,ERR2P,ERR2N,ERR2,CL and IERR
*                       NC=1: means only fit with syst. included is available
*                       NC=2: means fit with stat. only is also available
*             CVAL (1,I) = combined value (with syst. if I=1, stat. only if I=2)
*             ERR2P(1,I) = positive error**2 on CVAL(1,I)
*             ERR2N(1,I) = negative error**2 on CVAL(1,I)
*             ERR2 (1,I) = symmetric error**2 on CVAL(1,I)
*             CL   (1,I) = confidence level of fit performed to get CVAL(1,I)
*                       (or -1 if not available)
*             IERR(1,I)= error code from fit performed to get CVAL(1,I);
*                        0 means OK
*             NS     = number of steps
*
*     Output:  none
*     -------
*
      IMPLICIT NONE
*
*     Arguments
*
      INTEGER N,ICOMB(N),MS,NC,IERR(0:MS,NC),NS
      DOUBLE PRECISION CVAL(MS,NC),CL(MS,NC),
     &                 ERR2P(MS,NC),ERR2N(MS,NC),ERR2(MS,NC)
*
      INCLUDE 'combos.inc'
      INCLUDE 'master.inc'
*
      INTEGER MDISPLAY
      PARAMETER(MDISPLAY=MMEAS*2)
      INTEGER NUSED,JU(MMEAS),NDISPLAY,ID(MDISPLAY),ND(MMEAS)
*
*     Externals
*
      INTEGER LENOCC
*
*     Local variables
*
      DOUBLE PRECISION SYTOTP(MMEAS),SYTOTN(MMEAS),SYTOT(MMEAS)
      INTEGER IC,LTOT,LSTA,MC
      PARAMETER(LTOT=1,LSTA=2,MC=2) ! cases
      INTEGER I,J,K,ISUP,NDIS,IANAL
      CHARACTER*20 CH20
      CHARACTER*80 CH80
      DOUBLE PRECISION INFINITY
      PARAMETER(INFINITY=1.73205D+19/MMEAS) ! = sqrt(3.e+38)/MMEAS
*
*     Determine which measurements are actually used (not deweighted)
*
      J=0
      DO I=1,NMEAS
        IF(STAT(I).LT.INFINITY) THEN
          J=J+1
          JU(J)=I
        ENDIF
      ENDDO
      NUSED=J
      DO I=1,NMEAS
        IF(STAT(I).GE.INFINITY) THEN
          J=J+1
          JU(J)=I
        ENDIF
      ENDDO
*
*     Determine list of analyses for display
*
      NDISPLAY=0
      DO I=1,NMEAS
        ISUP=-1
        NDIS=NDISPLAY
        IF(KDISP(0,ICOMB(I)).EQ.-1) THEN ! this analysis supersedes
          DO K=1,NLINK(ICOMB(I)) ! loop on links
            IF(RHOSTA(I,ICOMB(I)).GT.1.) THEN ! supercede
              ISUP=0
              DO IANAL=1,NANAL ! search for superseded analysis
                IF(LINK(K,ICOMB(I)).EQ.LINK(0,IANAL)) ISUP=IANAL
              ENDDO
              IF(ISUP.EQ.0) THEN 
                CALL ANALYSIS_NAME(.TRUE.,ICOMB(I),CH80)
                CALL COMBOS_ERROR(1,
     &           'At least one analysis superseded by '//
     &           CH80(:LENOCC(CH80))//' not found','measurement used '//
     &           'in the combination will be displayed instead')
                NDISPLAY=NDIS
                GOTO 1
              ENDIF
              IF(NDISPLAY.EQ.MDISPLAY) CALL COMBOS_ERROR(-1,
     &         'Too many analyses to display',
     &         'please increase parameter MDISPLAY in WRITE_KUMAC')
              IF(NDISP(ISUP).GE.0.AND.KDISP(0,ISUP).EQ.-1) THEN
                CALL ANALYSIS_NAME(.TRUE.,ISUP,CH80)
                CH20=SUPERSEDED
                CALL COMBOS_ERROR(1,
     &           'Nesting of ''DISPLAY '//CH20(:LENOCC(CH20))//
     &           ''' logical lines not allowed','measurement of '//
     &           CH80(:LENOCC(CH80))//' will be displayed')
              ENDIF
              NDISPLAY=NDISPLAY+1
              ID(NDISPLAY)=ISUP
            ENDIF
          ENDDO
        ENDIF
    1   IF(KDISP(0,ICOMB(I)).NE.-1.OR.ISUP.EQ.0) THEN ! anal. doesn't supersede
          IF(NDISPLAY.EQ.MDISPLAY) CALL COMBOS_ERROR(-1,
     &     'Too many analyses to display',
     &     'please increase parameter MDISPLAY in WRITE_KUMAC')
          NDISPLAY=NDISPLAY+1
          ID(NDISPLAY)=ICOMB(I)
        ENDIF
        ND(I)=NDISPLAY-NDIS
      ENDDO
*
*     Perform checks
*
      IF(N.NE.NMEAS.OR.NUSED.NE.NMEFF.OR.
     &   J.NE.N.OR.(NC.NE.1.AND.NC.NE.2)) THEN 
        PRINT * ,'Fatal error in WRITE_KUMAC:'
        PRINT * ,'N    ,NMEAS=',N    ,NMEAS
        PRINT * ,'NUSED,NMEFF=',NUSED,NMEFF
        PRINT * ,'J    ,NC   =',J    ,NC
        STOP 7779
      ENDIF
*
*     Write KUMAC file
*
      WRITE(LUNKUM,6001) CHCOMB(:LENOCC(CHCOMB))
      IF(CHSTEP.EQ.' ') THEN
        WRITE(LUNKUM,6002) CHROUT(:LENOCC(CHROUT))
      ELSE
        WRITE(LUNKUM,6003) CHROUT(:LENOCC(CHROUT)),ISTEP
      ENDIF
      IF(CHSTEP.NE.' ')
     &WRITE(LUNKUM,6004) 'combos_step',CHSTEP(:LENOCC(CHSTEP))
      WRITE(LUNKUM,6004) 'combos_meas',CHMEAS(:LENOCC(CHMEAS))
      WRITE(LUNKUM,6004) 'combos_rout',CHROUT(:LENOCC(CHROUT))
      WRITE(LUNKUM,6004) 'combos_comb',CHCOMB(:LENOCC(CHCOMB))
*
*     Write list of analyses included in combination
*
      DO J=1,NMEAS
        WRITE(CH20,'(I20)') J
        CALL CLEFT(CH20,1,LEN(CH20))
        CH20='combos_'//CH20
        WRITE(LUNKUM,6004) CH20(:LENOCC(CH20)),
     &                     CHANAL(JU(J))(:LENOCC(CHANAL(JU(J))))
*OS        IF(NDISP(ICOMB(JU(J))).GE.0.AND.
*OS     &     KDISP(0,ICOMB(JU(J))).EQ.-1) THEN
*OS          WRITE(LUNKUM,6555) JU(J),ND(J)
*OS 6555     FORMAT('* analysis ',I5,' supersedes',I2,' analyses')
*OS        ENDIF
      ENDDO
*
*     Write list of analyses to be displayed
*
      DO J=1,NDISPLAY
        WRITE(CH20,'(I20)') J
        CALL CLEFT(CH20,1,LEN(CH20))
        CH20='combos_display_'//CH20
        CALL ANALYSIS_NAME(.FALSE.,ID(J),CH80)
        WRITE(LUNKUM,6004) CH20(:LENOCC(CH20)),CH80(:LENOCC(CH80))
      ENDDO
*
*     Write number of analyses in combinations  (c_no),
*           number of analyses effectively used (c_na), 
*           number of analyses to be displayed  (c_nd),
*           number of analyses to be displayed for each analysis (c_ndisp),
*           value of step variable at this step (c_step)
*
      WRITE(LUNKUM,6006) 'c_no',NMEAS,'c_na',NUSED,'c_nd',NDISPLAY
      WRITE(LUNKUM,6005) 'c_ndisp',NMEAS,(FLOAT(ND(J)),J=1,NMEAS)
      IF(CHSTEP.NE.' ') WRITE(LUNKUM,6007) 'c_step',VSTEP
*
*     Original measurements and measurements to be displayed
*
      IF(CHSTEP.EQ.' ') THEN
        WRITE(LUNKUM,6005) 'c_omeas',N,( CONT(1, 0,ICOMB(JU(J))),J=1,N)
        WRITE(LUNKUM,6005) 'c_ostap',N,(+CONT(1,+1,ICOMB(JU(J))),J=1,N)
        WRITE(LUNKUM,6005) 'c_ostan',N,(-CONT(1,-1,ICOMB(JU(J))),J=1,N)
        WRITE(LUNKUM,6005) 'c_osysp',N,(+CONT(1,+2,ICOMB(JU(J))),J=1,N)
        WRITE(LUNKUM,6005) 'c_osysn',N,(-CONT(1,-2,ICOMB(JU(J))),J=1,N)
        WRITE(LUNKUM,6005)
     &              'c_dmeas',NDISPLAY,( DISP(1, 0,ID(J)),J=1,NDISPLAY)
        WRITE(LUNKUM,6005)
     &              'c_dstap',NDISPLAY,(+DISP(1,+1,ID(J)),J=1,NDISPLAY)
        WRITE(LUNKUM,6005)
     &              'c_dstan',NDISPLAY,(-DISP(1,-1,ID(J)),J=1,NDISPLAY)
        WRITE(LUNKUM,6005)
     &              'c_dsysp',NDISPLAY,(+DISP(1,+2,ID(J)),J=1,NDISPLAY)
        WRITE(LUNKUM,6005)
     &              'c_dsysn',NDISPLAY,(-DISP(1,-2,ID(J)),J=1,NDISPLAY)
      ENDIF
*
*     Adjusted measurements
*
      DO I=1,NMEAS
        SYTOTP(I)=USYSP(I)**2
        SYTOTN(I)=USYSN(I)**2
        SYTOT(I)=USYS(I)**2
        DO J=1,NCSYS
          SYTOTP(I)=SYTOTP(I)+DMAX1(CSYSP(I,J),CSYSN(I,J),0.D0)**2
          SYTOTN(I)=SYTOTN(I)+DMIN1(CSYSP(I,J),CSYSN(I,J),0.D0)**2
          SYTOT(I)=SYTOT(I)+CSYS(I,J)**2
        ENDDO
        SYTOTP(I)=DSQRT(SYTOTP(I))
        SYTOTN(I)=DSQRT(SYTOTN(I))
        SYTOT(I)=DSQRT(SYTOT(I))
      ENDDO
      WRITE(LUNKUM,6005) 'c_ameas',NUSED,(MEAS(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_astap',NUSED,(STATP(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_astan',NUSED,(STATN(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_astas',NUSED,(STAT(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_asysp',NUSED,(SYTOTP(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_asysn',NUSED,(SYTOTN(JU(J)),J=1,NUSED)
      WRITE(LUNKUM,6005) 'c_asyss',NUSED,(SYTOT(JU(J)),J=1,NUSED)
*
*     Results of combination
*     Loop on cases (total error, then stat. only)
*    
      WRITE(LUNKUM,6005) 'c_resu',NC,(       CVAL(1,IC) ,IC=1,NC)
      WRITE(LUNKUM,6005) 'c_errp',NC,(DSQRT(ERR2P(1,IC)),IC=1,NC)
      WRITE(LUNKUM,6005) 'c_errn',NC,(DSQRT(ERR2N(1,IC)),IC=1,NC)
      WRITE(LUNKUM,6005) 'c_errs',NC,(DSQRT(ERR2 (1,IC)),IC=1,NC)
      WRITE(LUNKUM,6005) 'c_prob',NC,(         CL(1,IC) ,IC=1,NC)
      WRITE(LUNKUM,6005) 'c_ierr',NC,(FLOAT(IERR(1,IC)) ,IC=1,NC)
      WRITE(LUNKUM,6008) NS
      RETURN
*
      ENTRY WRITE_KUMAC_INIT
*     ======================
*
      WRITE(LUNKUM,6000) 
 6000 FORMAT(
     & /,'macro combos_output comb='' '' rout='' '' step=1',/,
     & /,'* Input arguments:',
     & /,'*',
     & /,'*   comb = name of combined result',
     & /,'*           (default is '' '' to indicate first result)',
     & /,'*   rout = name of routine used to perform combination',
     & /,'*           (default is '' '' to indicate first routine)',
     & /,'*   step = step number',
     & /,'*           (default is 1 to indicate first step)'
     & /,'*',
     & /,'* After a call to this macro, variable [@] contains',
     &    ' the maximum step number',
     & /,'* valid for [comb] and [rout],',
     &    ' or 0 if the input arguments don''t refer to data',
     & /,'* available in this file.',/,
     & /,'g/delete combos_* ; sigma c_=0 ; v/delete c_*',/)
 6001 FORMAT('if (([comb].eq.'' '').OR.([comb].eq.''',A,''')).AND. _')
 6002 FORMAT('   (([rout].eq.'' '').OR.([rout].eq.''',A,''')) then')
 6003 FORMAT('   (([rout].eq.'' '').OR.([rout].eq.''',A,''')).AND.',
     &           '([step].eq.',I4.4,') then')
 6004 FORMAT(2X,'g/create ',A,' ''',A,'''')
 6005 FORMAT(2X,'v/create ',A,'(',I3.3,') r ',(T27,4(G12.5,1X),:,'_'))
 6006 FORMAT(2X,3('sigma ',A,' = ',I3,:,' ; '))
 6007 FORMAT(2X,'sigma ',A,' = ',G12.5)
 6008 FORMAT(2X,'exitm ',I4,' ; endif')
      END
