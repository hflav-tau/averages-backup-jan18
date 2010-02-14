************************************************************************
*
*     Error-Matrix inversion combination routine of the COMBOS program
*
*     Paolo Checchia, INFN/PD
*
*     Version 2.30, May 12, 1997:
*     -- original release
*
************************************************************************
*
      SUBROUTINE BLUE(CVAL,ERR2P,ERR2N,CL,IERR)
*     =========================================
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
*     External
*
      REAL prob
*
*     Local variables
*
      DOUBLE PRECISION ave,w,chi2,sigma2
      real*8 csyss(mmeas,mmeas)
      real*8 ww(mmeas,mmeas),cov(mmeas,mmeas),err2(mmeas,mmeas)
      real*4 dum(mmeas*2)
      INTEGER n,i,j,ii
      REAL chisq
*      
      IERR = -1
      CVAL = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      CL   = 0.D0
*OS      IF(LUNIT.GT.0) WRITE(LUNIT,*) 'BLUE called'
*
*     Return result 
*

      n=nmeas

      call vzero(cov,mmeas*mmeas*2)
      call vzero(csyss,mmeas*mmeas*2)
      call vzero(err2,mmeas*mmeas*2)

*
****  define systematics error matrix
*
      do i=1,n
       do ii=1,n
        if(i.eq.ii) then
* diagonal terms
         do j=1,ncsys
          csyss(i,ii)=csyss(i,ii)+csys(i,j)**2
         enddo
         csyss(i,i)=csyss(i,i)+usys(i)**2
        else
* non-diagonal terms
         do j=1,ncsys
          csyss(i,ii)=csyss(i,ii) +
     +    csys(i,j)*csys(ii,j)
         enddo
        endif
       enddo
      enddo

*
****  define statistics error matrix
*
      do i=1,n
       do j=1,n
        cov(i,j)=stat(i)*stat(j)*stacor(i,j)
       enddo
      enddo


*
****  define total error matrix
*
      do i=1,n
       do j=1,n
        err2(i,j)=cov(i,j)+csyss(i,j)
       enddo
      enddo 
*
****  compute the weight matrix
*
      call ucopy(err2,ww,mmeas*mmeas*2)
      call dinv(n,ww,mmeas,dum,ierr)      

*
*     Return result 
*
      IF(IERR.EQ.0) THEN
*DR if several parameters plug here (not done yet)
       w=0.D0
       ave=0.
       do i=1,n
        do j=1,n
         ave=ave+meas(i)*ww(i,j)
         w=w+ww(i,j)
        enddo
       enddo
       ave=ave/w
       sigma2=1.D0/w

       chi2=0.d0
       do i=1,n
        do j=1,n
         chi2=chi2+(ave-meas(i))*(ave-meas(j))*ww(i,j)
        enddo
       enddo

       CVAL = ave
       ERR2P = sigma2
       ERR2N = sigma2
       chisq= sngl(chi2)
       if(nmeff.gt.1) CL = dble(prob(chisq,nmeff-1))
      ELSE ! combination failed, return dummy values
        CVAL = 0.D0
        ERR2P = 0.D0
        ERR2N = 0.D0
        CL   = 0.D0
      ENDIF
      END
