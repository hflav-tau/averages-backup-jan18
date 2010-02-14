************************************************************************
*
*     Error-matrix inversion combination routine of the COMBOS program
*     with handling of asymmetric uncertainties.
*
*     Paolo Checchia, INFN/PD
*
*     Version 2.32, July 16, 1997:
*     -- original release
*
************************************************************************
*
      SUBROUTINE BLUE_ASYM(CVAL,ERR2P,ERR2N,CL,IERR)
*     ==============================================
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
      DOUBLE PRECISION ERR2P ! positive total uncertainty on CVAL
      DOUBLE PRECISION ERR2N ! negative total uncertainty on CVAL
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
      DOUBLE PRECISION derf
*
*     Local variables
*
      DOUBLE PRECISION ave,w,chi2,sigma2,disc
      real*8 csyss(mmeas,mmeas),rhotot(mmeas,mmeas)
      real*8 ww(mmeas,mmeas),cov(mmeas,mmeas),err2(mmeas,mmeas)
      real*4 dum(mmeas*2)
      INTEGER n,i,j,ii,iter,maxiter
      REAL chisq
*      
      IERR = -1
      CVAL = 0.D0
      ERR2P = 0.D0
      ERR2N = 0.D0
      CL   = 0.D0
      IF(LUNIT.GT.0) WRITE(LUNIT,*) 'BLUE_ASYM called'
*
*     Return result 
*

* values to be given as input (for iteration)
      maxiter=10
      disc=0.00001
*
      n=nmeas

* 3 iteration for asymmetric errors
      do iter=1,maxiter
    
      call vzero(cov,mmeas*mmeas*2)
      call vzero(csyss,mmeas*mmeas*2)
      call vzero(err2,mmeas*mmeas*2)
      call vzero(rhotot,mmeas*mmeas*2)

*
****  define systematics error matrix
*
      do i=1,n
       do ii=1,n
        if(i.eq.ii) then
* diagonal terms
         if(iter.eq.1) ave=meas(i)
         do j=1,ncsys
          csyss(i,ii)=csyss(i,ii)
     +     +(.5*derf(ave-meas(i))*(csysp(i,j)+csysn(i,j))+csys(i,j))**2
         enddo
         csyss(i,i)=csyss(i,i)
     +  +(.5*derf(ave-meas(i))
     *  *(dabs(usysp(i))-dabs(usysn(i)))+usys(i))**2
     
        else
* non-diagonal terms
         do j=1,ncsys
          csyss(i,ii)=csyss(i,ii)
     +     +(.5*derf(ave-meas(i))*(csysp(i,j)+csysn(i,j))+csys(i,j))
     +     *(.5*derf(ave-meas(ii))*(csysp(ii,j)+csysn(ii,j))+csys(ii,j))
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


* printout for correlation coefficient
c      do i=1,n
c       do j=1,n
c        rhotot(i,j)=err2(i,j)/dsqrt(err2(i,i)*err2(j,j))
c       enddo
c       write(lunit,*), (rhotot(i,j),j=1,n)
c      enddo
      
*
****  compute the weight matrix
*
      call ucopy(err2,ww,mmeas*mmeas*2)
      call dinv(n,ww,mmeas,dum,ierr)      

*
*     Return result 
*
      IF(IERR.EQ.0) THEN

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

* check the variation at the present iteration

       if(cval.ne.0..and.abs(cval-ave)/dsqrt(sigma2).lt.disc) then
        CVAL = ave
        ERR2P = sigma2
        ERR2N = sigma2
        chisq= sngl(chi2)
        if(nmeff.gt.1) CL = dble(prob(chisq,nmeff-1))
        go to 99
       else
        CVAL = ave
        ERR2P = sigma2
        ERR2N = sigma2
        chisq= sngl(chi2)
        if(nmeff.gt.1) CL = dble(prob(chisq,nmeff-1))
       endif
      ELSE ! combination failed, return dummy values
        CVAL = 0.D0
        ERR2P = 0.D0
        ERR2N = 0.D0
        CL   = 0.D0
      ENDIF

      IF(LUNIT.GT.0)
     & WRITE(LUNIT,*) iter,'result', ave,'+-',dsqrt(sigma2)

      enddo
 99   continue

      END
