c-----------------------------------------------------------------------
      program mpfft2csxtest
      implicit none
      integer indx, indy, indnvp, mshare
      integer nx, ny, nx1, ny1, nxh, nyh, nvp, kyp, kxp2, kblok, j2blok
      integer nmx, nxy, nmxh, nxhy, nxv, nyv, nxyh, nxvh, nyvh, nxyh2
      integer kyp1, kxp21, ndim, nmt
c indnvp = exponent which determines number of virtual processors
c mshare = (0,1) = (no,yes) architecture is shared memory
      parameter( indx =   7, indy =   8, indnvp =   2, mshare =   0)
      parameter(nx=2**indx,ny=2**indy)
      parameter(nxh=nx/2,nyh=ny/2,nx1=nx+1,ny1=ny+1)
      parameter(nvp=2**indnvp,kyp=(ny-1)/nvp+1,kxp2=(nx-1)/nvp+1)
      parameter(kblok=1+mshare*(ny/kyp-1),j2blok=1+mshare*(nx/kxp2-1))
      parameter(nmx=nx*(ny/nx)+ny*(nx/ny),nxy=nmx/(2-nx/nmx-ny/nmx))
      parameter(nmxh=nxh*(ny/nxh)+ny*(nxh/ny))
      parameter(nxhy=nmxh/(2-nxh/nmxh-ny/nmxh))
      parameter(nxv=nx+2,nyv=ny+2,nxvh=nxv/2,nyvh=nyv/2)
      parameter(nxyh=nxy/2,nxyh2=2*nxyh)
      parameter(kyp1=kyp+1,kxp21=kxp2+1)
      parameter(ndim=2)
      parameter(nmt=1)
c
      integer ntpose, idproc, kstrt, ks, isign, nproc, ntasks, ierr
      integer i, j, k, l, kk, k1, j1, joff, koff
      real epsmax, eps, time, ttp
      double precision ranorm
      real f, t, g, h, bs, br
      real ff, tt, gg, hh
      integer mixup
      complex sct2, ss, z1, z2
      dimension f(nxv,kyp1,kblok), t(nyv,kxp21,j2blok)
      dimension g(nxv,ny1), h(nyv,nx1)
      dimension mixup(nxhy), sct2(nxyh2)
      dimension bs(ndim,kxp21,kyp1,kblok), br(ndim,kxp21,kyp1,j2blok)
      dimension ff(ndim,nxv,kyp1,kblok), tt(ndim,nyv,kxp21,j2blok)
      dimension gg(ndim,nxv,ny1), hh(ndim,nyv,nx1)
      dimension ss(ndim,nxv)
      complex gc, ggc, ttc
      dimension gc(nxv,nyh), ggc(ndim,nxv,nyh)
      dimension ttc(ndim,nyvh,kxp21,j2blok)
      equivalence (g,gc), (gg,ggc), (tt,ttc)
c
      integer iftask, nxyip
      dimension iftask(nmt), nxyip(nmt)
c ntpose = (0,1) = (no,yes) input, output data are transposed in pfft2r
      data ntpose /0/
c initialize for multiprocessing
      call MP_INIT(nproc)
      if (nproc.eq.0) then
         write (*,*) 'MPLibrary not installed'
         ntasks = 0
      else
c return the number of processors on host computer
         ntasks = min(nmt,nproc-1)
         write (*,*) nproc, ' processors found, ', ntasks+1, ' used'
      endif
c debug
c     ntasks = 0
c
c initialize for parallel processing
      call PPINIT(idproc,nvp)
      kstrt = idproc + 1
      ks = kstrt - 2
      write (70+kstrt,*) nproc, ' processors found, ', ntasks+1, ' used'
c prepare fft tables
      call WFST2RINIT(mixup,sct2,indx,indy,nxhy,nxy)
c create test function
      do 30 k = 1, ny
      kk = (k - 1)/kyp
      k1 = k - kyp*kk
      do 20 j = 1, nx1
      g(j,k) = ranorm()
      h(k,j) = g(j,k)
      do 10 l = 1, kblok
      if (kk.eq.(l+ks)) f(j,k1,l) = g(j,k)
   10 continue
   20 continue
   30 continue
c last point is special
      do 50 j = 1, nx1
      g(j,ny1) = ranorm()
      h(ny1,j) = g(j,ny1)
      do 40 l = 1, kblok
      if ((nvp-1).eq.(l+ks)) f(j,kyp1,l) = g(j,ny1)
   40 continue
   50 continue
c create real vector test function
      do 45 k = 1, ny
      kk = (k - 1)/kyp
      k1 = k - kyp*kk
      do 35 j = 1, nx1
      do 25 i = 1, ndim
      gg(i,j,k) = ranorm()
      hh(i,k,j) = gg(i,j,k)
      do 15 l = 1, kblok
      if (kk.eq.(l+ks)) ff(i,j,k1,l) = gg(i,j,k)
   15 continue
   25 continue
   35 continue
   45 continue
c last point is special
      do 75 j = 1, nx1
      do 65 i = 1, ndim
      gg(i,j,ny1) = ranorm()
      hh(i,ny1,j) = gg(i,j,ny1)
      do 55 l = 1, kblok
      if ((nvp-1).eq.(l+ks)) ff(i,j,kyp1,l) = gg(i,j,ny1)
   55 continue
   65 continue
   75 continue
      call TIMERA(-1,'total   ',time)
c start special test case
      isign = -1
      write (71,*) 'kblok=',kblok
c-----------------------------------------------------------------------
      call MPFSFT2R(f,t,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,kstr
     1t,nxvh,nyvh,kxp2,kyp,kyp1,kxp21,j2blok,kblok,nxhy,nxy,nxyip,iftask
     2,ntasks,ierr)
c     call MPFCFT2R(f,t,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,kstr
c    1t,nxvh,nyvh,kxp2,kyp,kyp1,kxp21,j2blok,kblok,nxhy,nxy,nxyip,iftask
c    2,ntasks,ierr)
      if (ierr.ne.0) then
         call MP_END()
         stop
      endif
      isign = -1
c-----------------------------------------------------------------------
      call WFSFT2RX(g,ss,isign,mixup,sct2,indx,indy,nxvh,nyvh,nxhy,nxy)
c     call WFCFT2RX(g,ss,isign,mixup,sct2,indx,indy,nxvh,nyvh,nxhy,nxy)
      do 95 k = 1, ny1
      do 85 j = 1, nx1
      h(k,j) = g(j,k)
   85 continue
   95 continue
c end special test case
      call TIMERA(1,'total   ',time)
      epsmax = 0.0
      if (kstrt.gt.ny) go to 100
      do 90 l = 1, kblok
      koff = kyp*(l + ks)
      do 80 k = 1, kyp
      k1 = k + koff
      do 70 j = 1, nx1
      eps = abs(f(j,k,l) - g(j,k1))
      write (60+kstrt,*) j,k,k1,f(j,k,l),g(j,k1),eps
      if (eps.gt.epsmax) then
         write (70+kstrt,*) j,k,k1,f(j,k,l),g(j,k1),eps
         epsmax = eps
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      write (70+kstrt,*) 'local epsmax=',epsmax
      call PSUM(epsmax,eps,1,1)
      write (70+kstrt,*) 'global epsmax=',epsmax
      epsmax = 0.0
      if (kstrt.gt.nx) go to 140
      do 130 l = 1, j2blok
      joff = kxp2*(l + ks)
      do 120 j = 1, kxp21
      if (((l+ks).lt.(nvp-1)).and.(j.eq.kxp21)) go to 120
      j1 = j + joff
      do 110 k = 1, nyh
      z1 = cmplx(t(2*k-1,j,l),t(2*k,j,l))
      z2 = gc(j1,k)
      eps = abs(z1 - z2)
      if (eps.gt.epsmax) then
         write (70+kstrt,*) k,j,j1,z1,z2,eps
         epsmax = eps
      endif
  110 continue
  120 continue
  130 continue
  140 continue
      write (70+kstrt,*) 'local transpose epsmax=',epsmax
      call PSUM(epsmax,eps,1,1)
      write (70+kstrt,*) 'global transpose epsmax=',epsmax
c verify real vector data
c     isign = -1
c     call MPFCST2R2(ff,tt,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,k
c    1strt,nxvh,nyv,kxp2,kyp,kyp1,kxp21,jblok,kblok,nxhy,nxy,nxyip,iftas
c    2k,ntasks,ierr)
c     call MPFSCT2R2(ff,tt,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,k
c    1strt,nxvh,nyv,kxp2,kyp,kyp1,kxp21,jblok,kblok,nxhy,nxy,nxyip,iftas
c    2k,ntasks,ierr)
c     call WFCST2R2(gg,isign,mixup,sct2,indx,indy,nxvh,nyv,nxhy,nxy)
c     call WFSCT2R2(gg,isign,mixup,sct2,indx,indy,nxvh,nyv,nxhy,nxy)
c
c     call MPFCST2R3(ff,tt,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,k
c    1strt,nxvh,nyv,kxp2,kyp,kyp1,kxp21,jblok,kblok,nxhy,nxy,nxyip,iftas
c    2k,ntasks,ierr)
c     call MPFSCT2R3(ff,tt,bs,br,isign,ntpose,mixup,sct2,ttp,indx,indy,k
c    1strt,nxvh,nyv,kxp2,kyp,kyp1,kxp21,jblok,kblok,nxhy,nxy,nxyip,iftas
c    2k,ntasks,ierr)
c     call WFCST2R3(gg,isign,mixup,sct2,indx,indy,nxvh,nyv,nxhy,nxy)
c     call WFSCT2R3(gg,isign,mixup,sct2,indx,indy,nxvh,nyv,nxhy,nxy)
c     epsmax = 0.0
c     if (kstrt.gt.ny) go to 190
c     do 180 l = 1, kblok
c     koff = kyp*(l + ks)
c     do 170 k = 1, kyp
c     k1 = k + koff
c     do 160 j = 1, nx1
c     do 150 i = 1, ndim
c     eps = abs(ff(i,j,k,l) - gg(i,j,k1))
c     if (eps.gt.epsmax) then
c        write (70+kstrt,*) i,j,k,k1,ff(i,j,k,l),gg(i,j,k1),eps
c        epsmax = eps
c     endif
c 150 continue
c 160 continue
c 170 continue
c 180 continue
c 190 continue
c     write (70+kstrt,*) 'local vector epsmax=',epsmax
c     call PSUM(epsmax,eps,1,1)
c     write (70+kstrt,*) 'global vector epsmax=',epsmax
c     epsmax = 0.0
c     if (kstrt.gt.nx) go to 240
c     do 230 l = 1, j2blok
c     joff = kxp2*(l + ks)
c     do 220 j = 1, kxp21
c     if (((l+ks).lt.(nvp-1)).and.(j.eq.kxp21)) go to 220
c     j1 = j + joff
c     do 210 k = 1, nyh
c     do 200 i = 1, ndim
c     z1 = ttc(i,k,j,l)
c     z2 = ggc(i,j1,k)
c     eps = abs(z1 - z2)
c     if (eps.gt.epsmax) then
c        write (70+kstrt,*) i,k,j,j1,z1,z2,eps
c        epsmax = eps
c     endif
c 200 continue
c 210 continue
c 220 continue
c 230 continue
c 240 continue
c     write (70+kstrt,*) 'local vector transpose epsmax=',epsmax
c     call PSUM(epsmax,eps,1,1)
c     write (70+kstrt,*) 'global vector transpose epsmax=',epsmax
      call PPEXIT
      call MP_END()
      stop
      end
c-----------------------------------------------------------------------
      subroutine PPINIT(idproc,nvp)
c this subroutine initializes parallel processing
c input: nvp, output: idproc
c idproc = processor id
c nvp = number of real or virtual processors requested
      implicit none
      integer idproc, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror, ndprec
      save /PPARMS/
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      data ndprec /1/
c this segment is used for shared memory computers
c     nproc = nvp
c     idproc = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c initialize the MPI execution environment
      call MPI_INIT(ierror)
      if (ierror.ne.0) stop
      lgrp = MPI_COMM_WORLD
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
c set default datatypes
         mint = MPI_INTEGER
c single precision
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c requested number of processors not obtained
      if (nproc.ne.nvp) then
         write (2,*) ' processor number error: nvp, nproc=', nvp, nproc
         call PPEXIT
         stop
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPEXIT
c this subroutine terminates parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      integer ierror
c synchronize processes
      call MPI_BARRIER(lgrp,ierror)
c terminate MPI execution environment
      call MPI_FINALIZE(ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
c this subroutine calculates tables needed by a two dimensional
c fast real sine and cosine transforms and their inverses.
c input: indx, indy, nxhyd, nxyd
c output: mixup, sctd
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyd
      integer mixup
      complex sctd
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nxy
      dnxy = 0.5*6.28318530717959/float(nxy)
      do 30 j = 1, nxy
      arg = dnxy*float(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kxyi
     2p,iftask,nmt,ierr)
c multi-tasking real sine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFST2RXX, PFDFT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x sine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         write (*,*) 'task kblok=',kblok
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish x sine transform
         write (*,*) 'sub kblok=',kblok
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c debug:
      if (indx.ge.0) return
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y ffts
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y ffts
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x sine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kxyi
     2p,iftask,nmt,ierr)
c multi-tasking real cosine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFCT2RXX, PFDFT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y ffts
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y ffts
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x cosine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,ky
     1pp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c f(n,k,i) = (1/nx*ny)*sum(f(j,k,i)*sin(pi*n*j/nx))
c if isign = 1, a forward sine transform is performed
c f(j,k,i) = sum(f(n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = .5*at2
      f(j,k,l) = at1 + at2
      f(nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,k,l) = 0.0
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 60 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 150 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = ani*(t2 + t4)
         f(2*j,k,l) = ani*(t3 + t5)
         f(nx3-2*j,k,l) = ani*(t2 - t4)
         f(nx3-2*j+1,k,l) = ani*(t5 - t3)
  120    continue
  130    continue
         ani = 2.*ani
         do 140 k = kypi, kypt
         f(nxh+1,k,l) = ani*f(nxh+1,k,l)
         f(nxh+2,k,l) = -ani*f(nxh+2,k,l)
         t2 = ani*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = ani*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = ani*f(nx+1,k,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 170 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 160 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = t2 + t4
         f(2*j,k,l) = t3 + t5
         f(nx3-2*j,k,l) = t2 - t4
         f(nx3-2*j+1,k,l) = t5 - t3
  160    continue
  170    continue
         do 180 k = kypi, kypt
         f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
         f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
         t2 = 2.0*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  180    continue
  190    continue
      endif
c perform recursion for sine transform
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      sum1 = .5*f(1,k,l)
      f(1,k,l) = 0.0
      f(2,k,l) = sum1
      do 200 j = 2, nxh
      sum1 = sum1 + f(2*j-1,k,l)
      f(2*j-1,k,l) = -f(2*j,k,l)
      f(2*j,k,l) = sum1
  200 continue
      f(nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,ky
     1pp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c f(n,k,i) = (1/nx*ny)*(.5*f(1,k,i) + ((-1)**n)*f(nx+1,k,i)
c            + sum(f(j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, a forward cosine transform is performed
c f(j,k,i) = 2*(.5*f(1,k,i) + ((-1)**j)*f(n+1,k,i) + sum(f(n,k,i)*
c            cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp+1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,k,l) - f(nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = .5*at1
      f(j,k,l) = at1 - at2
      f(nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,k,l) = .5*(f(1,k,l) + f(nx+1,k,l))
      f(nx+1,k,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 60 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 150 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = ani*(t2 + t4)
         f(2*j,k,l) = ani*(t3 + t5)
         f(nx3-2*j,k,l) = ani*(t2 - t4)
         f(nx3-2*j+1,k,l) = ani*(t5 - t3)
  120    continue
  130    continue
         ani = 2.*ani
         do 140 k = kypi, kypt
         f(nxh+1,k,l) = ani*f(nxh+1,k,l)
         f(nxh+2,k,l) = -ani*f(nxh+2,k,l)
         t2 = ani*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = ani*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = ani*f(nx+1,k,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 170 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 160 k = kypi, kypt
         t4 = f(nx3-2*j,k,l)
         t5 = -f(nx3-2*j+1,k,l)
         t2 = f(2*j-1,k,l) + t4
         t3 = f(2*j,k,l) + t5
         t6 = f(2*j-1,k,l) - t4
         t5 = f(2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k,l) = t2 + t4
         f(2*j,k,l) = t3 + t5
         f(nx3-2*j,k,l) = t2 - t4
         f(nx3-2*j+1,k,l) = t5 - t3
  160    continue
  170    continue
         do 180 k = kypi, kypt
         f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
         f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
         t2 = 2.0*(f(1,k,l) + f(2,k,l))
         f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
         f(1,k,l) = t2
         f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  180    continue
  190    continue
      endif
c perform recursion for cosine transform
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      sum1 = f(nx+1,k,l)
      f(nx+1,k,l) = f(2,k,l)
      f(2,k,l) = sum1
      do 200 j = 2, nxh
      sum1 = sum1 - f(2*j,k,l)
      f(2*j,k,l) = sum1
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,
     1kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c when the x transform is a real sine or cosine transform
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c g(k,j,i) = sum(g(m,n,i)*exp(sqrt(-1)*2pi*m*k/ny)
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyvh = first dimension of field arrays, must be >= nyh
c kxp2d = second dimension of field arrays, must be >= kxp2+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1), 0 < k < ny/2
c g(1,j,i) = real part of mode jj-1,ky=0
c g(2,j,i) = real part of mode jj-1,ky=ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp2, kxpi, kxpp
      integer nyvh, kxp2d, j2blok, nxhyd, nxyd
      complex g, sctd
      dimension g(nyvh,kxp2d,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, nyh2
      integer nxy, nxhy, ks, kxps, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.gt.0) go to 130
c inverse fourier transform
      nry = nxhy/nyh
c bit-reverse array elements in y
      kxpt = kxps
      do 30 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
   10 continue
   20 continue
   30 continue
c first transform in y
      nry = nxy/nyh
      do 80 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 70 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sctd(1+kmr*(j-1))
      do 40 i = kxpi, kxpt
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes ky = 0, ny/2 and normalize
      kmr = nxy/nyh
      kxpt = kxps
      ani = 0.5
      do 120 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 100 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 90 j = kxpi, kxpt
      t = conjg(g(nyh2-k,j,l))
      s = g(k,j,l) + t
      t = (g(k,j,l) - t)*t1
      g(k,j,l) = ani*(s + t)
      g(nyh2-k,j,l) = ani*conjg(s - t)
   90 continue
  100 continue
      do 110 j = kxpi, kxpt
      g(nyhh+1,j,l) = conjg(g(nyhh+1,j,l))
      g(1,j,l) = cmplx(real(g(1,j,l)) + aimag(g(1,j,l)),real(g(1,j,l)) -
     1 aimag(g(1,j,l)))
  110 continue
  120 continue
      return
c forward fourier transform
  130 kmr = nxy/nyh
      kxpt = kxps
      do 190 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
c scramble coefficients
      do 150 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 140 j = kxpi, kxpt
      t = conjg(g(nyh2-k,j,l))
      s = g(k,j,l) + t
      t = (g(k,j,l) - t)*t1
      g(k,j,l) = s + t
      g(nyh2-k,j,l) = conjg(s - t)
  140 continue
  150 continue
      do 160 j = kxpi, kxpt
      g(nyhh+1,j,l) = 2.*conjg(g(nyhh+1,j,l))
      g(1,j,l) = cmplx(real(g(1,j,l)) + aimag(g(1,j,l)),real(g(1,j,l)) -
     1 aimag(g(1,j,l)))
  160 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      do 180 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 180
      do 170 j = kxpi, kxpt
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  170 continue
  180 continue
  190 continue
c then transform in y
      nry = nxy/nyh
      do 240 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 230 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sctd(1+kmr*(j-1)))
      do 200 i = kxpi, kxpt
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  200 continue
  210 continue
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer i, j, l, nx, ny, nx1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nx1
	 do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y ffts
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y ffts
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 60 l = 1, j2blok
         do 50 j = 1, nx1
	 do 40 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer i, j, l, nx, ny, nx1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nx1
	 do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y ffts
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y ffts
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 60 l = 1, j2blok
         do 50 j = 1, nx1
	 do 40 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,k
     1ypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + ((-1)**n)*f(1,nx+1,k,i)
c              + sum(f(1,j,k,i)*cos(pi*n*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*sum(f(2,j,k,i)*sin(pi*n*j/nx))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = 2*(.5*f(1,1,k,i) + ((-1)**j)*f(1,n+1,k,i)
c              + sum(f(1,n,k,i)*cos(pi*n*j/nx))
c f(2,j,k,i) = sum(f(2,n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,1,k,l) - f(1,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+2-j,k,l) = at1 + at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,1,k,l) = .5*(f(1,1,k,l) + f(1,nx+1,k,l))
      f(1,nx+1,k,l) = sum1
      f(2,1,k,l) = 0.0
      f(2,nxh+1,k,l) = 2.0*f(2,nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 2
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 2
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = f(1,nx+1,k,l)
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = sum1
      sum2 = .5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      f(2,2,k,l) = sum2
      do 260 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 + f(2,2*j-1,k,l)
      f(2,2*j-1,k,l) = -f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
  260 continue
      f(2,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,k
     1ypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*sum(f(1,j,k,i)*sin(pi*n*j/nx))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,1,k,i) + ((-1)**n)*f(2,nx+1,k,i)
c              + sum(f(2,j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = sum(f(1,n,k,i)*sin(pi*n*j/nx))
c f(2,j,k,i) = 2*(.5*f(2,1,k,i) + ((-1)**j)*f(2,n+1,k,i)
c              + sum(f(2,n,k,i)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(2,1,k,l) - f(2,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+2-j,k,l) = at1 - at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,1,k,l) = 0.0
      f(1,nxh+1,k,l) = 2.0*f(1,nxh+1,k,l)
      f(2,1,k,l) = .5*(f(2,1,k,l) + f(2,nx+1,k,l))
      f(2,nx+1,k,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 2
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 2
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 2
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = .5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      f(1,2,k,l) = sum1
      sum2 = f(2,nx+1,k,l)
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = sum2
      do 260 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k,l)
      f(1,2*j-1,k,l) = -f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 - f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
  260 continue
      f(1,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,
     1kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms and its inverses, for a subset of x,
c when the x transform is a real sine or cosine transform
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:2,n,m,i) = sum(f(1:2,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c g(1:2,k,j,i) = sum(g(1:2,m,n,i)*exp(sqrt(-1)*2pi*m*k/ny)
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyvh = first dimension of field arrays, must be >= nyh
c kxp2d = second dimension of field arrays, must be >= kxp2+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1), 0 < k < ny/2
c g(1,j,i) = real part of mode jj-1,ky=0
c g(2,j,i) = real part of mode jj-1,ky=ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp2, kxpi, kxpp
      integer nyvh, kxp2d, j2blok, nxhyd, nxyd
      complex g, sctd
      dimension g(2,nyvh,kxp2d,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, nyh2
      integer nxy, nxhy, ks, kxps, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c swap complex components
      kxpt = kxps
      do 30 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      do 10 k = 1, nyh
      at1 = aimag(g(1,k,j,l))
      g(1,k,j,l) = cmplx(real(g(1,k,j,l)),real(g(2,k,j,l)))
      g(2,k,j,l) = cmplx(at1,aimag(g(2,k,j,l)))
   10 continue
   20 continue
   30 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      kxpt = kxps
      do 60 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
   40 continue
   50 continue
   60 continue
c first transform in y
      nry = nxy/nyh
      do 110 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 100 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble modes ky = 0, ny/2 and normalize
      kmr = nxy/nyh
      kxpt = kxps
      ani = 0.5
      do 170 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 140 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 130 j = kxpi, kxpt
      do 120 jj = 1, 2
      t = conjg(g(jj,nyh2-k,j,l))
      s = g(jj,k,j,l) + t
      t = (g(jj,k,j,l) - t)*t1
      g(jj,k,j,l) = ani*(s + t)
      g(jj,nyh2-k,j,l) = ani*conjg(s - t)
  120 continue
  130 continue
  140 continue
      do 160 j = kxpi, kxpt
      do 150 jj = 1, 2
      g(jj,nyhh+1,j,l) = conjg(g(jj,nyhh+1,j,l))
      g(jj,1,j,l) = cmplx(real(g(jj,1,j,l)) + aimag(g(jj,1,j,l)),real(g(
     1jj,1,j,l)) - aimag(g(jj,1,j,l)))
  150 continue
  160 continue
  170 continue
      return
c forward fourier transform
  180 kmr = nxy/nyh
      kxpt = kxps
      do 260 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
c scramble coefficients
      do 210 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 200 j = kxpi, kxpt
      do 190 jj = 1, 2
      t = conjg(g(jj,nyh2-k,j,l))
      s = g(jj,k,j,l) + t
      t = (g(jj,k,j,l) - t)*t1
      g(jj,k,j,l) = s + t
      g(jj,nyh2-k,j,l) = conjg(s - t)
  190 continue
  200 continue
  210 continue
      do 230 j = kxpi, kxpt
      do 220 jj = 1, 2
      g(jj,nyhh+1,j,l) = 2.*conjg(g(jj,nyhh+1,j,l))
      g(jj,1,j,l) = cmplx(real(g(jj,1,j,l)) + aimag(g(jj,1,j,l)),real(g(
     1jj,1,j,l)) - aimag(g(jj,1,j,l)))
  220 continue
  230 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      do 250 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
  240 continue
  250 continue
  260 continue
c then transform in y
      nry = nxy/nyh
      do 310 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 300 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sctd(1+kmr*(j-1)))
      do 270 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c swap complex components
      kxpt = kxps
      do 340 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 330 j = kxpi, kxpt
      do 320 k = 1, nyh
      at1 = aimag(g(1,k,j,l))
      g(1,k,j,l) = cmplx(real(g(1,k,j,l)),real(g(2,k,j,l)))
      g(2,k,j,l) = cmplx(at1,aimag(g(2,k,j,l)))
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer i, j, l, nx, ny, nx1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nx1
	 do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y ffts
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y ffts
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 60 l = 1, j2blok
         do 50 j = 1, nx1
	 do 40 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer i, j, l, nx, ny, nx1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nx1
	 do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y ffts
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y ffts
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 60 l = 1, j2blok
         do 50 j = 1, nx1
	 do 40 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,
     1kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + ((-1)**n)*f(1,nx+1,k,i)
c              + sum(f(1,j,k,i)*cos(pi*n*j/nx)))
c f(2:3,n,k,i) = (1/nx*ny)*sum(f(2:3,j,k,i)*sin(pi*n*j/nx))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = 2*(.5*f(1,1,k,i) + ((-1)**j)*f(1,n+1,k,i)
c              + sum(f(1,n,k,i)*cos(pi*n*j/nx))
c f(2:3,j,k,i) = sum(f(2:3,n,k,i)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(1,1,k,l) - f(1,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+2-j,k,l) = at1 + at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+2-j,k,l) = at1 - at2
      at2 = f(3,nx+2-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(3,j,k,l) = at1 + at2
      f(3,nx+2-j,k,l) = at1 - at2
   10 continue
      f(1,1,k,l) = .5*(f(1,1,k,l) + f(1,nx+1,k,l))
      f(1,nx+1,k,l) = sum1
      f(2,1,k,l) = 0.0
      f(2,nxh+1,k,l) = 2.0*f(2,nxh+1,k,l)
      f(3,1,k,l) = 0.0
      f(3,nxh+1,k,l) = 2.0*f(3,nxh+1,k,l)
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 3
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 3
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = f(1,nx+1,k,l)
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = sum1
      sum2 = .5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      f(2,2,k,l) = sum2
      sum3 = .5*f(3,1,k,l)
      f(3,1,k,l) = 0.0
      f(3,2,k,l) = sum3
      do 260 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 + f(2,2*j-1,k,l)
      f(2,2*j-1,k,l) = -f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
      sum3 = sum3 + f(3,2*j-1,k,l)
      f(3,2*j-1,k,l) = -f(3,2*j,k,l)
      f(3,2*j,k,l) = sum3
  260 continue
      f(2,nx+1,k,l) = 0.0
      f(3,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,
     1kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k,i) = (1/nx*ny)*sum(f(1,j,k,i)*sin(pi*n*j/nx))
c f(2:3,n,k,i) = (1/nx*ny)*(.5*f(2:3,1,k,i) + ((-1)**n)*f(2:3,nx+1,k,i)
c              + sum(f(2:3,j,k,i)*cos(pi*n*j/nx)))
c if isign = 1, forward sine transforms are performed
c f(1,j,k,i) = sum(f(1,n,k,i)*sin(pi*n*j/nx))
c f(2:3,j,k,i) = 2*(.5*f(2:3,1,k,i) + ((-1)**j)*f(2:3,n+1,k,i)
c              + sum(f(2:3,n,k,i)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kyp = number of data values per block in y
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f >= nx/2 + 1
c kypd = second dimension of f >= kyp + 1
c kblok = number of data blocks in y
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kyp, kypi, kypp
      integer nxvh, kypd, kblok, nxhyd, nxyd
      real f
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kyps = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = .5*(f(2,1,k,l) - f(2,nx+1,k,l))
      sum2 = .5*(f(3,1,k,l) - f(3,nx+1,k,l))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+2-j,k,l) = at1 - at2
      at2 = f(2,nx+2-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+2-j,k,l) = at1 + at2
      at2 = f(3,nx+2-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      sum2 = sum2 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k,l) = at1 - at2
      f(3,nx+2-j,k,l) = at1 + at2
   10 continue
      f(1,1,k,l) = 0.0
      f(1,nxh+1,k,l) = 2.0*f(1,nxh+1,k,l)
      f(2,1,k,l) = .5*(f(2,1,k,l) + f(2,nx+1,k,l))
      f(2,nx+1,k,l) = sum1
      f(3,1,k,l) = .5*(f(3,1,k,l) + f(3,nx+1,k,l))
      f(3,nx+1,k,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 120 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kypi, kypt
      do 80 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         kypt = kyps
         do 190 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 160 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 150 k = kypi, kypt
         do 140 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = ani*(t2 + t4)
         f(jj,2*j,k,l) = ani*(t3 + t5)
         f(jj,nx3-2*j,k,l) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         ani = 2.*ani
         do 180 k = kypi, kypt
         do 170 jj = 1, 3
         f(jj,nxh+1,k,l) = ani*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -ani*f(jj,nxh+2,k,l)
         t2 = ani*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = ani*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = ani*f(jj,nx+1,k,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         kypt = kyps
         do 250 l = 1, kblok
         if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
         do 220 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 210 k = kypi, kypt
         do 200 jj = 1, 3
         t4 = f(jj,nx3-2*j,k,l)
         t5 = -f(jj,nx3-2*j+1,k,l)
         t2 = f(jj,2*j-1,k,l) + t4
         t3 = f(jj,2*j,k,l) + t5
         t6 = f(jj,2*j-1,k,l) - t4
         t5 = f(jj,2*j,k,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k,l) = t2 + t4
         f(jj,2*j,k,l) = t3 + t5
         f(jj,nx3-2*j,k,l) = t2 - t4
         f(jj,nx3-2*j+1,k,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 k = kypi, kypt
         do 230 jj = 1, 3
         f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
         f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
         t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
         f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
         f(jj,1,k,l) = t2
         f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for cosine-sine transform
      kypt = kyps
      do 280 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 270 k = kypi, kypt
      sum1 = .5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      f(1,2,k,l) = sum1
      sum2 = f(2,nx+1,k,l)
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = sum2
      sum3 = f(3,nx+1,k,l)
      f(3,nx+1,k,l) = f(3,2,k,l)
      f(3,2,k,l) = sum3
      do 260 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k,l)
      f(1,2*j-1,k,l) = -f(1,2*j,k,l)
      f(1,2*j,k,l) = sum1
      sum2 = sum2 - f(2,2*j,k,l)
      f(2,2*j,k,l) = sum2
      sum3 = sum3 - f(3,2*j,k,l)
      f(3,2*j,k,l) = sum3
  260 continue
      f(1,nx+1,k,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,
     1kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and its inverses, for a subset of x,
c when the x transform is a real sine or cosine transform
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n,m,i) = sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c g(1:3,k,j,i) = sum(g(1:3,m,n,i)*exp(sqrt(-1)*2pi*m*k/ny)
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of x indices used
c nyvh = first dimension of field arrays, must be >= nyh
c kxp2d = second dimension of field arrays, must be >= kxp2+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1), 0 < k < ny/2
c g(1,j,i) = real part of mode jj-1,ky=0
c g(2,j,i) = real part of mode jj-1,ky=ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp2, kxpi, kxpp
      integer nyvh, kxp2d, j2blok, nxhyd, nxyd
      complex g, sctd
      dimension g(3,nyvh,kxp2d,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, nyh2
      integer nxy, nxhy, ks, kxps, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c swap complex components
      kxpt = kxps
      do 30 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 i = kxpi, kxpt
      do 10 k = 1, nyh
      at1 = real(g(3,k,i,l))
      g(3,k,i,l) = cmplx(real(g(2,k,i,l)),aimag(g(3,k,i,l)))
      at2 = aimag(g(2,k,i,l))
      g(2,k,i,l) = cmplx(aimag(g(1,k,i,l)),at1)
      g(1,k,i,l) = cmplx(real(g(1,k,i,l)),at2)
   10 continue
   20 continue
   30 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      kxpt = kxps
      do 60 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in y
      nry = nxy/nyh
      do 110 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 100 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble modes ky = 0, ny/2 and normalize
      kmr = nxy/nyh
      kxpt = kxps
      ani = 0.5
      do 170 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 140 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 130 j = kxpi, kxpt
      do 120 jj = 1, 3
      t = conjg(g(jj,nyh2-k,j,l))
      s = g(jj,k,j,l) + t
      t = (g(jj,k,j,l) - t)*t1
      g(jj,k,j,l) = ani*(s + t)
      g(jj,nyh2-k,j,l) = ani*conjg(s - t)
  120 continue
  130 continue
  140 continue
      do 160 j = kxpi, kxpt
      do 150 jj = 1, 3
      g(jj,nyhh+1,j,l) = conjg(g(jj,nyhh+1,j,l))
      g(jj,1,j,l) = cmplx(real(g(jj,1,j,l)) + aimag(g(jj,1,j,l)),real(g(
     1jj,1,j,l)) - aimag(g(jj,1,j,l)))
  150 continue
  160 continue
  170 continue
      return
c forward fourier transform
  180 kmr = nxy/nyh
      kxpt = kxps
      do 260 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
c scramble coefficients
      do 210 k = 2, nyhh
      t1 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 200 j = kxpi, kxpt
      do 190 jj = 1, 3
      t = conjg(g(jj,nyh2-k,j,l))
      s = g(jj,k,j,l) + t
      t = (g(jj,k,j,l) - t)*t1
      g(jj,k,j,l) = s + t
      g(jj,nyh2-k,j,l) = conjg(s - t)
  190 continue
  200 continue
  210 continue
      do 230 j = kxpi, kxpt
      do 220 jj = 1, 3
      g(jj,nyhh+1,j,l) = 2.*conjg(g(jj,nyhh+1,j,l))
      g(jj,1,j,l) = cmplx(real(g(jj,1,j,l)) + aimag(g(jj,1,j,l)),real(g(
     1jj,1,j,l)) - aimag(g(jj,1,j,l)))
  220 continue
  230 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      do 250 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
  240 continue
  250 continue
  260 continue
c then transform in y
      nry = nxy/nyh
      do 310 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 300 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sctd(1+kmr*(j-1)))
      do 270 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c swap complex components
      kxpt = kxps
      do 340 l = 1, j2blok
      if ((kxps+kxp2*(l+ks)).eq.nx) kxpt = kxps + 1
      do 330 i = kxpi, kxpt
      do 320 k = 1, nyh
      at1 = real(g(3,k,i,l))
      g(3,k,i,l) = cmplx(aimag(g(2,k,i,l)),aimag(g(3,k,i,l)))
      at2 = real(g(2,k,i,l))
      g(2,k,i,l) = cmplx(at1,aimag(g(1,k,i,l)))
      g(1,k,i,l) = cmplx(real(g(1,k,i,l)),at2)
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
c this subroutine calculates tables needed by a two dimensional
c fast real sine and cosine transforms and their inverses.
c input: indx, indy, nxhyd, nxyd
c output: mixup, sctd
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyd
      integer mixup
      complex sctd
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nxy
      dnxy = 0.5*6.28318530717959/float(nxy)
      do 30 j = 1, nxy
      arg = dnxy*float(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WFSFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd)
c wrapper function for real sine/periodic transform
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c debug:
         if (indx.ge.0) return
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFCFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd)
c wrapper function for real cosine/periodic transform
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c f(n,k) = (1/nx*ny)*sum(f(j,k)*sin(pi*n*j/nx))
c if isign = 1, a forward sine transform is performed
c f(j,k) = sum(f(n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f => nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = .5*at2
      f(j,k) = at1 + at2
      f(nx+2-j,k) = at1 - at2
   10 continue
      f(1,k) = 0.0
      f(nxh+1,k) = 2.0*f(nxh+1,k)
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 100 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 90 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = ani*(t2 + t4)
         f(2*j,k) = ani*(t3 + t5)
         f(nx3-2*j,k) = ani*(t2 - t4)
         f(nx3-2*j+1,k) = ani*(t5 - t3)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 k = nyi, nyt
         f(nxh+1,k) = ani*f(nxh+1,k)
         f(nxh+2,k) = -ani*f(nxh+2,k)
         t2 = ani*(f(1,k) + f(2,k))
         f(2,k) = ani*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = ani*f(nx+1,k)
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = t2 + t4
         f(2*j,k) = t3 + t5
         f(nx3-2*j,k) = t2 - t4
         f(nx3-2*j+1,k) = t5 - t3
  120    continue
  130    continue
         do 140 k = nyi, nyt
         f(nxh+1,k) = 2.0*f(nxh+1,k)
         f(nxh+2,k) = -2.0*f(nxh+2,k)
         t2 = 2.0*(f(1,k) + f(2,k))
         f(2,k) = 2.0*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = 2.0*f(nx+1,k)
  140    continue
      endif
c perform recursion for sine transform
      do 160 k = nyi, nyt
      sum1 = .5*f(1,k)
      f(1,k) = 0.0
      f(2,k) = sum1
      do 150 j = 2, nxh
      sum1 = sum1 + f(2*j-1,k)
      f(2*j-1,k) = -f(2*j,k)
      f(2*j,k) = sum1
  150 continue
      f(nx+1,k) = 0.0
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c f(n,k) = (1/nx*ny)*(.5*f(1,k) + ((-1)**n)*f(nx+1,k) + sum(f(j,k)*
c       cos(pi*n*j/nx)))
c if isign = 1, a forward cosine transform is performed
c f(j,k) = 2*(.5*f(1,k) + ((-1)**j)*f(n+1,k) + sum(f(n,k)*
c       cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,k) - f(nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = .5*at1
      f(j,k) = at1 - at2
      f(nx+2-j,k) = at1 + at2
   10 continue
      f(1,k) = .5*(f(1,k) + f(nx+1,k))
      f(nx+1,k) = sum1
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 100 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 90 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = ani*(t2 + t4)
         f(2*j,k) = ani*(t3 + t5)
         f(nx3-2*j,k) = ani*(t2 - t4)
         f(nx3-2*j+1,k) = ani*(t5 - t3)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 k = nyi, nyt
         f(nxh+1,k) = ani*f(nxh+1,k)
         f(nxh+2,k) = -ani*f(nxh+2,k)
         t2 = ani*(f(1,k) + f(2,k))
         f(2,k) = ani*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = ani*f(nx+1,k)
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = t2 + t4
         f(2*j,k) = t3 + t5
         f(nx3-2*j,k) = t2 - t4
         f(nx3-2*j+1,k) = t5 - t3
  120    continue
  130    continue
         do 140 k = nyi, nyt
         f(nxh+1,k) = 2.0*f(nxh+1,k)
         f(nxh+2,k) = -2.0*f(nxh+2,k)
         t2 = 2.0*(f(1,k) + f(2,k))
         f(2,k) = 2.0*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = 2.0*f(nx+1,k)
  140    continue
      endif
c perform recursion for cosine transform
      do 160 k = nyi, nyt
      sum1 = f(nx+1,k)
      f(nx+1,k) = f(2,k)
      f(2,k) = sum1
      do 150 j = 2, nxh
      sum1 = sum1 - f(2*j,k)
      f(2*j,k) = sum1
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(j,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3
      dimension f(nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t2 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 80 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 70 j = nxi, nxt
      t2 = conjg(f(j,nyh2-k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      do 90 j = nxi, nxt
      f(j,nyhh+1) = conjg(f(j,nyhh+1))
      f(j,1) = cmplx(real(f(j,1)) + aimag(f(j,1)),real(f(j,1)) - aimag(f
     1(j,1)))
   90 continue
      return
c forward fourier transform
c scramble coefficients
  100 kmr = nxy/nyh
      do 120 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 110 j = nxi, nxt
      t2 = conjg(f(j,nyh2-k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(j,nyh2-k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 j = nxi, nxt
      f(j,nyhh+1) = 2.0*conjg(f(j,nyhh+1))
      f(j,1) = cmplx(real(f(j,1)) + aimag(f(j,1)),real(f(j,1)) - aimag(f
     1(j,1)))
  130 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 150 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 150
      do 140 j = nxi, nxt
      t2 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t2
  140 continue
  150 continue
c then transform in y
      nry = nxy/nyh
      do 190 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 160 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RLTOCX2(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c nxd = half of the first dimension of f
c nyd = second dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(2*nxd*nyd), s(2*nxd)
c local data
      integer j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 30 k = 1, nyh
         joff = nx2d*(k - 1)
         do 10 j = nxi, nxt
         s(2*j-1) = f(j+joff)
         s(2*j) = f(j+joff+nxd)
   10    continue
         do 20 j = nxi, nxt
         f(2*j+joff-1) = s(2*j-1)
         f(2*j+joff) = s(2*j)
   20    continue
   30    continue
c complex to real
      else if (isign.gt.0) then
         do 60 k = 1, nyh
         joff = nx2d*(k - 1)
         do 40 j = nxi, nxt
         s(2*j-1) = f(2*j+joff-1)
         s(2*j) = f(2*j+joff)
   40    continue
         do 50 j = nxi, nxt
         f(j+joff) = s(2*j-1)
         f(j+joff+nxd) = s(2*j)
   50    continue
   60    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFCSFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 2 real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFSCFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 2 real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
c             + sum(f(1,j,k)*cos(pi*n*j/nx)))
c f(2,n,k) = (1/nx*ny)*sum(f(2,j,k)*sin(pi*n*j/nx))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k) + sum(f(1,n,k)*
c       cos(pi*n*j/nx))
c f(2,j,k) = sum(f(2,n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,1,k) - f(1,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+2-j,k) = at1 + at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+2-j,k) = at1 - at2
   10 continue
      f(1,1,k) = .5*(f(1,1,k) + f(1,nx+1,k))
      f(1,nx+1,k) = sum1
      f(2,1,k) = 0.0
      f(2,nxh+1,k) = 2.0*f(2,nxh+1,k)   
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 2
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 2
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = f(1,nx+1,k)
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = sum1
      sum2 = .5*f(2,1,k)
      f(2,1,k) = 0.0
      f(2,2,k) = sum2
      do 210 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 + f(2,2*j-1,k)
      f(2,2*j-1,k) = -f(2,2*j,k)
      f(2,2*j,k) = sum2
  210 continue
      f(2,nx+1,k) = 0.0
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
c             + sum(f(2,j,k)*cos(pi*n*j/nx)))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
c f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
c             + sum(f(2,n,k)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(2,1,k) - f(2,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+2-j,k) = at1 - at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+2-j,k) = at1 + at2
   10 continue
      f(1,1,k) = 0.0
      f(1,nxh+1,k) = 2.0*f(1,nxh+1,k)
      f(2,1,k) = .5*(f(2,1,k) + f(2,nx+1,k))
      f(2,nx+1,k) = sum1
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 2
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 2
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = .5*f(1,1,k)
      f(1,1,k) = 0.0
      f(1,2,k) = sum1
      sum2 = f(2,nx+1,k)
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = sum2
      do 210 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k)
      f(1,2*j-1,k) = -f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 - f(2,2*j,k)
      f(2,2*j,k) = sum2
  210 continue
      f(1,nx+1,k) = 0.0
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:2,j,m) = sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(1:2,j,k) = sum(f(1:2,j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3
      dimension f(2,nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 90 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 80 j = nxi, nxt
      do 70 jj = 1, 2
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
   90 continue
      do 110 j = nxi, nxt
      do 100 jj = 1, 2
      f(jj,j,nyhh+1) = conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  100 continue
  110 continue
      return
c forward fourier transform
c scramble coefficients
  120 kmr = nxy/nyh
      do 150 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 140 j = nxi, nxt
      do 130 jj = 1, 2
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,j,nyh2-k) = conjg(t1 - t2)
  130 continue
  140 continue
  150 continue
      do 170 j = nxi, nxt
      do 160 jj = 1, 2
      f(jj,j,nyhh+1) = 2.0*conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  160 continue
  170 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 190 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  180 continue
  190 continue
c then transform in y
      nry = nxy/nyh
      do 230 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 200 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RLTOCX22(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms 2 real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c ndim = first dimension of f
c nxd = half of the second dimension of f
c nyd = third dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(2,2*nxd*nyd), s(2,2*nxd)
c local data
      integer i, j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 40 k = 1, nyh
         joff = nx2d*(k - 1)
         do 20 j = nxi, nxt
         do 10 i = 1, 2
         s(i,2*j-1) = f(i,j+joff)
         s(i,2*j) = f(i,j+joff+nxd)
   10    continue
   20    continue
         do 30 j = nxi, nxt
         f(1,2*j+joff-1) = s(1,2*j-1)
         f(2,2*j+joff-1) = s(1,2*j)
         f(1,2*j+joff) = s(2,2*j-1)
         f(2,2*j+joff) = s(2,2*j)
   30    continue
   40    continue
c complex to real
      else if (isign.gt.0) then
         do 80 k = 1, nyh
         joff = nx2d*(k - 1)
         do 50 j = nxi, nxt
         s(1,2*j-1) = f(1,2*j+joff-1)
         s(1,2*j) = f(2,2*j+joff-1)
         s(2,2*j-1) = f(1,2*j+joff)
         s(2,2*j)= f(2,2*j+joff)
   50    continue
         do 70 j = nxi, nxt
         do 60 i = 1, 2
         f(i,j+joff) = s(i,2*j-1)
         f(i,j+joff+nxd) = s(i,2*j)
   60    continue
   70    continue
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFCSFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 3 real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFSCFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 3 real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
c             + sum(f(1,j,k)*cos(pi*n*j/nx)))
c f(2:3,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k) + sum(f(1,n,k)*
c       cos(pi*n*j/nx))
c f(2:3,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,1,k) - f(1,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+2-j,k) = at1 + at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+2-j,k) = at1 - at2
      at2 = f(3,nx+2-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(3,j,k) = at1 + at2
      f(3,nx+2-j,k) = at1 - at2
   10 continue
      f(1,1,k) = .5*(f(1,1,k) + f(1,nx+1,k))
      f(1,nx+1,k) = sum1
      f(2,1,k) = 0.0
      f(2,nxh+1,k) = 2.0*f(2,nxh+1,k)   
      f(3,1,k) = 0.0
      f(3,nxh+1,k) = 2.0*f(3,nxh+1,k)
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 3
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 3
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = f(1,nx+1,k)
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = sum1
      sum2 = .5*f(2,1,k)
      f(2,1,k) = 0.0
      f(2,2,k) = sum2
      sum3 = .5*f(3,1,k)
      f(3,1,k) = 0.0
      f(3,2,k) = sum3
      do 210 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 + f(2,2*j-1,k)
      f(2,2*j-1,k) = -f(2,2*j,k)
      f(2,2*j,k) = sum2
      sum3 = sum3 + f(3,2*j-1,k)
      f(3,2*j-1,k) = -f(3,2*j,k)
      f(3,2*j,k) = sum3
  210 continue
      f(2,nx+1,k) = 0.0
      f(3,nx+1,k) = 0.0
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
c f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
c             + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
c f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
c             + sum(f(2:3,n,k)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(2,1,k) - f(2,nx+1,k))
      sum2 = .5*(f(3,1,k) - f(3,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at4 = real(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+2-j,k) = at1 - at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+2-j,k) = at1 + at2
      at2 = f(3,nx+2-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      sum2 = sum2 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k) = at1 - at2
      f(3,nx+2-j,k) = at1 + at2
   10 continue
      f(1,1,k) = 0.0
      f(1,nxh+1,k) = 2.0*f(1,nxh+1,k)
      f(2,1,k) = .5*(f(2,1,k) + f(2,nx+1,k))
      f(2,nx+1,k) = sum1
      f(3,1,k) = .5*(f(3,1,k) + f(3,nx+1,k))
      f(3,nx+1,k) = sum2
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 3
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 3
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = .5*f(1,1,k)
      f(1,1,k) = 0.0
      f(1,2,k) = sum1
      sum2 = f(2,nx+1,k)
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = sum2
      sum3 = f(3,nx+1,k)
      f(3,nx+1,k) = f(3,2,k)
      f(3,2,k) = sum3
      do 210 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k)
      f(1,2*j-1,k) = -f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 - f(2,2*j,k)
      f(2,2*j,k) = sum2
      sum3 = sum3 - f(3,2*j,k)
      f(3,2*j,k) = sum3
  210 continue
      f(1,nx+1,k) = 0.0
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:3,j,m) = sum(f(1:3,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(1:3,j,k) = sum(f(1:3,j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3, t4
      dimension f(3,nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 90 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 80 j = nxi, nxt
      do 70 jj = 1, 3
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
   90 continue
      do 110 j = nxi, nxt
      do 100 jj = 1, 3
      f(jj,j,nyhh+1) = conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  100 continue
  110 continue
      return
c forward fourier transform
c scramble coefficients
  120 kmr = nxy/nyh
      do 150 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 140 j = nxi, nxt
      do 130 jj = 1, 3
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,j,nyh2-k) = conjg(t1 - t2)
  130 continue
  140 continue
  150 continue
      do 170 j = nxi, nxt
      do 160 jj = 1, 3
      f(jj,j,nyhh+1) = 2.0*conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  160 continue
  170 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 190 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  180 continue
  190 continue
c then transform in y
      nry = nxy/nyh
      do 230 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 200 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine RLTOCX23(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms 2 real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c ndim = first dimension of f
c nxd = half of the second dimension of f
c nyd = third dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(3,2*nxd*nyd), s(3,2*nxd)
c local data
      integer i, j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 40 k = 1, nyh
         joff = nx2d*(k - 1)
         do 20 j = nxi, nxt
         do 10 i = 1, 3
         s(i,2*j-1) = f(i,j+joff)
         s(i,2*j) = f(i,j+joff+nxd)
   10    continue
   20    continue
         do 30 j = nxi, nxt
         f(1,2*j+joff-1) = s(1,2*j-1)
         f(2,2*j+joff-1) = s(1,2*j)
         f(3,2*j+joff-1) = s(2,2*j-1)
         f(1,2*j+joff) = s(2,2*j)
         f(2,2*j+joff) = s(3,2*j-1)
         f(3,2*j+joff) = s(3,2*j)
   30    continue
   40    continue
c complex to real
      else if (isign.gt.0) then
         do 80 k = 1, nyh
         joff = nx2d*(k - 1)
         do 50 j = nxi, nxt
         s(1,2*j-1) = f(1,2*j+joff-1)
         s(1,2*j) = f(2,2*j+joff-1)
         s(2,2*j-1) = f(3,2*j+joff-1)
         s(2,2*j) = f(1,2*j+joff)
         s(3,2*j-1) = f(2,2*j+joff)
         s(3,2*j) = f(3,2*j+joff)
   50    continue
         do 70 j = nxi, nxt
         do 60 i = 1, 3
         f(i,j+joff) = s(i,2*j-1)
         f(i,j+joff+nxd) = s(i,2*j)
   60    continue
   70    continue
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      integer r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      subroutine TIMERA(icntrl,chr,time)
c this subroutine performs timing
c input: icntrl, chr
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c chr = character variable for labeling timings
c time = elapsed time in seconds
c written for mpi
      implicit none
      integer icntrl
      character*8 chr
      real time
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx
c lgrp = current communicator
c mreal = default datatype for reals
      common /pparms/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer idproc, ierr
      real nclock, mclock
      double precision jclock
      save jclock
   91 format (1x,a8,1x,'max/min real time = ',e14.7,1x,e14.7,1x,'sec')
      data jclock /0.0d0/
      if (icntrl.eq.0) return
      if (icntrl.eq.1) go to 10
c initialize clock
      call MPI_BARRIER(lgrp,ierr)
      jclock = MPI_WTIME()
      return
c read clock and write time difference from last clock initialization
   10 nclock = real(MPI_WTIME() - jclock)
      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MIN,lgrp,ierr)
      mclock = time
      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MAX,lgrp,ierr)
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      if (idproc.eq.0) write (6,91) chr, time, mclock
      return
      end
c-----------------------------------------------------------------------
      subroutine PSUM(f,g,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output data
c g = scratch array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx
c lstat = length of status array
      parameter(lstat=8)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /pparms/ nproc, lgrp, mreal, mint, mcplx
c local data
      integer istatus
      integer idproc, ierr, kstrt, ks, l, kxs, k, kb, lb, msid, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine PWTIMERA(icntrl,time,dtime)
c this subroutine performs local wall clock timing
c input: icntrl, dtime
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c time = elapsed time in seconds
c dtime = current time
c written for mpi
      implicit none
      integer icntrl
      real time
      double precision dtime
c local data
      double precision jclock
      double precision MPI_WTIME
      external MPI_WTIME
c initialize clock
      if (icntrl.eq.(-1)) then
         dtime = MPI_WTIME()
c read clock and write time difference from last clock initialization
      else if (icntrl.eq.1) then
         jclock = dtime
         dtime = MPI_WTIME()
         time = real(dtime - jclock)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PRTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for each row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(nxv,kypd,kblok), g(nyv,kxpd,jblok)
      dimension s(kxp+1,kyp+1,kblok), t(kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 20 k = 1, kyp1
c     do 10 j = 1, kxp1
c     g(k+koff,j,l) = f(j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,l),kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp,msi
     1d,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(j,k,l) = f(j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,l),kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgrp,ie
     1rr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(k+koff,j,l) = t(j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:2,k+kyp*(m-1),j,l) = f(1:2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(2,nxv,kypd,kblok), g(2,nyv,kxpd,jblok)
      dimension s(2,kxp+1,kyp+1,kblok), t(2,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),2*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),2*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:3,k+kyp*(m-1),j,l) = f(1:3,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(3,nxv,kypd,kblok), g(3,nyv,kxpd,jblok)
      dimension s(3,kxp+1,kyp+1,kblok), t(3,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c     g(3,k+koff,j,l) = f(3,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),3*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
         s(3,j,k,l) = f(3,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),3*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
         g(3,k+koff,j,l) = t(3,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
