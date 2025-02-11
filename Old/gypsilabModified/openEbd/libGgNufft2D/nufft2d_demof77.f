cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory). 
cc
      program testfft
      implicit none
c
      integer i,ier,iflag,j,k1,k2,lused,lw,mx,ms,mt,n1,n2,nj,nk
      parameter (mx=256*256)
      parameter (lw=1000 000)
      real*8 xj(mx),yj(mx)
      real *8 sk(mx),tk(mx)
      real*8 err,pi,eps,salg,ealg
      real*8 fw(0:lw-1)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(mx),cj0(mx),cj1(mx)
      complex*16 fk0(mx),fk1(mx)
c
c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
c
      n1 = 36
      n2 = 40
      ms = 32
      mt = 30
      nj = n1*n2
      do k1 = -n1/2, (n1-1)/2
         do k2 = -n2/2, (n2-1)/2
            j = (k2+n2/2+1) + (k1+n1/2)*n2
            xj(j) = pi*dcos(-pi*k1/n1)
            yj(j) = pi*dcos(-pi*k2/n2)
            cj(j) = dcmplx(dsin(pi*j/n1),dcos(pi*j/n2))
         enddo
      enddo
c
c     -----------------------
c     start tests
c     -----------------------
c
      iflag = 1
      print*,'Starting 2D testing: ', ' nj =',nj, ' ms,mt =',ms,mt
      do i = 1,4
         if (i.eq.1) eps=1d-4
         if (i.eq.2) eps=1d-8
         if (i.eq.3) eps=1d-12
         if (i.eq.4) eps=1d-16
c extended/quad precision tests
         if (i.eq.5) eps=1d-20
         if (i.eq.6) eps=1d-24
         if (i.eq.7) eps=1d-28
         if (i.eq.8) eps=1d-32
	 print*,' '
	 print*,' Requested precision eps =',eps
	 print*,' '
c
c     -----------------------
c     call 2D Type 1 method
c     -----------------------
c
         call dirft2d1(nj,xj,yj,cj,iflag,ms,mt,fk0)
         call nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt,fk1,fw,lw,lused,ier)
         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' iflag = ',iflag
         print *, ' ier = ',ier
         print *, ' meomory used = ',lused
         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' type 1 err = ',err
c
c     -----------------------
c      call 2D Type 2 method
c     -----------------------
         call dirft2d2(nj,xj,yj,cj0,iflag,ms,mt,fk0)
         call nufft2d2(nj,xj,yj,cj1,iflag,eps,ms,mt,fk1,
     1                 fw,lw,lused,ier)
         print *, ' ier = ',ier
         print *, ' meomory used = ',lused
         call errcomp(cj0,cj1,nj,err)
         print *, ' type 1 err = ',err
c
c     -----------------------
c      call 2D Type3 method
c     -----------------------
         nk = ms*mt
         do k1 = 1, nk
            sk(k1) = 64*(dcos(k1*pi/nk))
            tk(k1) = 48*(dsin(-pi/2+k1*pi/nk))
         enddo

         call dirft2d3(nj,xj,yj,cj,iflag, nk,sk,tk,fk0)
         call nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk,fk1,
     1                 fw,lw,lused,ier)
c
         print *, ' ier = ',ier
         print *, ' meomory used = ',lused
         call errcomp(fk0,fk1,nk,err)
         print *, ' type 1 err = ',err
      enddo 
      stop
      end
c
c
c
c
c
      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err =sqrt(ealg/salg)
      return
      end
