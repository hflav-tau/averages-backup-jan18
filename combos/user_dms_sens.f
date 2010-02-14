************************************************************************
*
*     Preparation routine of the COMBOS program
*
*     Get user value for dms sensitivities
*
*     Olivier Schneider, IPHE/Lausanne
*
*     Version 2.98, September 7, 2000:
*     -- original release as routine USER_DMS_SENS
*
************************************************************************
*
      SUBROUTINE USER_DMS_SENS(N,ICOMB,AMEAS)
*     =======================================
*
*     Input:  N      = number of analyses to combine
*     -----   ICOMB  = array of indexes to analyses
*                      ICOMB(I) = analysis number of Ith analysis to combine
*             AMEAS  = ajustments already performed on central values
*                      by routine MASTER
*                      AMEAS(I,ICONT) = adjustment due to systematic
*                                       contribution ICONT performed on
*                                       Ith analysis to combine
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
*****************************
*
      IMPLICIT NONE
      INCLUDE 'combos.inc'
      INCLUDE 'master.inc'
*
*     Argument
*
      INTEGER N
      INTEGER ICOMB(N)
      DOUBLE PRECISION AMEAS(MMEAS,MCONT)
*
*     External 
*  
      INTEGER LENOCC
*
*     Local variable
*
      LOGICAL XDEF
      REAL VDEF(3)
      INTEGER J,I,COUNT,IANAL
*
      IF(ISTEP.EQ.1) THEN ! initialize
        IF(N.NE.NMEAS) CALL COMBOS_ERROR(-1,
     &                 'Inconsistency in USER_DMS_SENS',
     &                 'please check input arguments')
        XDEF=.FALSE.
        DO J=1,NPARA
          IF(CHPARA(J).EQ.CHROUT) THEN
            VDEF(1)=PARA(J)
            VDEF(2)=+EXCUP(J)
            VDEF(3)=-EXCUN(J)
            XDEF=.TRUE.
          ENDIF
        ENDDO
        IF(LUNIT.GT.0) THEN ! write what the routine is doing  
          IF(XDEF) THEN 
            WRITE(LUNIT,1000) 
            WRITE(LUNIT,2000) 'default',VDEF
          ELSE
            WRITE(LUNIT,1001) CHROUT(:LENOCC(CHROUT))
            RETURN
          ENDIF
 1000     FORMAT(/,'Preparation routine USER_DMS_SENS will look for',
     &           /,'for user values of the dms sensitivities in the',
     &           /,'individual analyses and use these values for the',
     &           /,'combined analysis:')
 2000     FORMAT(3X,A,' values = ',3(2X,G10.4))
 1001     FORMAT(/,'Preparation routine USER_DMS_SENS will do nothing',
     &           /,'because no parameter ',A,
     &           /,'is present in the combined analysis.',/)
        ENDIF
        COUNT=0
        DO I=1,N
          IANAL=ICOMB(I)
          DO J=1,NPAR(IANAL)
            IF(CHNAM(KPAR(J,IANAL)).EQ.CHROUT) THEN
              COUNT=COUNT+1
              VDEF(1)=PAR(0,J,IANAL)
              VDEF(2)=PAR(+1,J,IANAL)
              VDEF(3)=PAR(-1,J,IANAL)
              WRITE(LUNIT,2000) 'user   ',VDEF
            ENDIF
          ENDDO
        ENDDO
        IF(COUNT.GT.1) CALL COMBOS_ERROR(1,
     &   'More than one set of user values found in USER_DMS_SENS',
     &   'arbitrarily use latest values; check results')
        CALL SET_DMS_LIMIT_USER(.TRUE.,VDEF)
      ENDIF
      END
