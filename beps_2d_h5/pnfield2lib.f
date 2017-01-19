c 2d parallel PIC library for solving field equations with neumann
c boundary conditions
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: february 6, 2005
c-----------------------------------------------------------------------
      subroutine PPOISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny
     12d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with neumann
c boundary conditions (zero normal efield), for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,fy,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fx, fy, ffd, zero
      dimension q(ny2d,kxp2,j2blok)
      dimension fx(ny2d,kxp2,j2blok), fy(ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fx(k,j,l) = cmplx(0.,at2)
         fx(k1,j,l) = cmplx(0.,at2)
         fy(k,j,l) = cmplx(0.,at3)
         fy(k1,j,l) = cmplx(0.,-at3)
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
c mode numbers ky = 0, ny
         at1 = real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at3 = -at1*real(q(1,j,l))
         at2 = dkx*at3
         fx(1,j,l) = cmplx(0.,at2)
         fx(ny+1,j,l) = zero
         fy(1,j,l) = zero
         fy(ny+1,j,l) = zero
         wp = wp + .5*at1*real(q(1,j,l))**2
      endif
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at3 = -at1*real(q(k,1,l))
         at2 = dny*float(k - 1)*at3
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
         fy(k,1,l) = cmplx(0.,at2)
         fy(k1,1,l) = zero
         wp = wp + .5*at1*real(q(k,1,l))**2
   70    continue
         fx(1,1,l) = zero
         fx(ny+1,1,l) = zero
         fy(1,1,l) = zero
         fy(ny+1,1,l) = zero
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         at3 = at2*real(q(k,j,l))
         fx(k,j,l) = cmplx(at3,0.)
         fx(k1,j,l) = cmplx(at3,0.)
         wp = wp + at1*real(q(k,j,l))**2
  110    continue
c mode numbers ky = 0, ny
         at2 = real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         at3 = at2*real(q(1,j,l))
         fx(1,j,l) = cmplx(at3,0.)
         fx(ny+1,j,l) = zero
         wp = wp + .5*at1*real(q(1,j,l))**2
      endif
  120 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         at2 = real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         at3 = at2*real(q(k,1,l))
         fx(k,1,l) = cmplx(at3,0.)
         fx(k1,1,l) = zero
         wp = wp + .5*at1*real(q(k,1,l))**2
  130    continue
         fx(1,1,l) = zero
         fx(ny+1,1,l) = zero
      endif
  140 continue
  150 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         at2 = at1*real(q(k,j,l))
         fy(k,j,l) = cmplx(at2,0.)
         fy(k1,j,l) = cmplx(at2,0.)
  170    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         at2 = at1*real(q(1,j,l))
         fy(1,j,l) = cmplx(at2,0.)
         fy(ny+1,j,l) = zero
      endif
  180 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,1,l))
         at2 = at1*real(q(k,1,l))
         fy(k,1,l) = cmplx(at2,0.)
         fy(k1,1,l) = zero
  190    continue
         at1 = aimag(ffd(1,1,l))
         fy(1,1,l) = cmplx(at1*real(q(1,1,l)),0.)
         fy(ny+1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2
     1d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(ny2d,kxp2,j2blok), fxy(2,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = cmplx(0.,at2)
         fxy(2,k,j,l) = cmplx(0.,at3)
         fxy(1,k1,j,l) = cmplx(0.,at2)
         fxy(2,k1,j,l) = cmplx(0.,-at3)
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
c mode numbers ky = 0, ny
         at1 = real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at3 = -at1*real(q(1,j,l))
         at2 = dkx*at3
         fxy(1,1,j,l) = cmplx(0.,at2)
         fxy(2,1,j,l) = zero
         fxy(1,ny+1,j,l) = zero
         fxy(2,ny+1,j,l) = zero
         wp = wp + .5*at1*real(q(1,j,l))**2
      endif
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at3 = -at1*real(q(k,1,l))
         at2 = dny*float(k - 1)*at3
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = cmplx(0.,at2)
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         wp = wp + .5*at1*real(q(k,1,l))**2
   70    continue
         fxy(1,1,1,l) = zero
         fxy(2,1,1,l) = zero
         fxy(1,ny+1,1,l) = zero
         fxy(2,ny+1,1,l) = zero
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2
     1d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with neumann boundary conditions (zero normal efield),
c for distributed data.  Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(ny2d,kxp2,j2blok), fxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = cmplx(0.,at2)
         fxy(2,k,j,l) = cmplx(0.,at3)
         fxy(3,k,j,l) = zero
         fxy(1,k1,j,l) = cmplx(0.,at2)
         fxy(2,k1,j,l) = cmplx(0.,-at3)
         fxy(3,k1,j,l) = zero
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
c mode numbers ky = 0, ny
         at1 = real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at3 = -at1*real(q(1,j,l))
         at2 = dkx*at3
         fxy(1,1,j,l) = cmplx(0.,at2)
         fxy(2,1,j,l) = zero
         fxy(3,1,j,l) = zero
         fxy(1,ny+1,j,l) = zero
         fxy(2,ny+1,j,l) = zero
         fxy(3,ny+1,j,l) = zero
         wp = wp + .5*at1*real(q(1,j,l))**2
      endif
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at3 = -at1*real(q(k,1,l))
         at2 = dny*float(k - 1)*at3
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = cmplx(0.,at2)
         fxy(3,k,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
         wp = wp + .5*at1*real(q(k,1,l))**2
   70    continue
         fxy(1,1,1,l) = zero
         fxy(2,1,1,l) = zero
         fxy(3,1,1,l) = zero
         fxy(1,ny+1,1,l) = zero
         fxy(2,ny+1,1,l) = zero
         fxy(3,ny+1,1,l) = zero
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPNX2(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c j2blok = number of data blocks
c kxp2 = number of data values per block
      complex cu, zero
      dimension cu(3,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*aimag(cu(1,k,j,l)) + dky*aimag(cu(2,k,j,l)))
         at3 = aimag(cu(1,k,j,l)) - dkx*at2
         at4 = aimag(cu(2,k,j,l)) - dky*at2
         cu(1,k,j,l) = cmplx(0.,at3)
         cu(2,k,j,l) = cmplx(0.,at4)
         cu(1,k1,j,l) = cmplx(0.,at3)
         cu(2,k1,j,l) = cmplx(0.,-at4)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPN2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c j2blok = number of data blocks
c kxp2 = number of data values per block
      dimension cu(3,nyv,kxp2,j2blok)
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         at3 = cu(1,k,j,l) - dkx*at2
         at4 = cu(2,k,j,l) - dky*at2
         cu(1,k,j,l) = at3
         cu(2,k,j,l) = at4
   10    continue
c mode numbers ky = 0, ny
         cu(1,1,j,l) = 0.
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         cu(2,k,1,l) = 0.
   30    continue
         cu(1,1,1,l) = 0.
         cu(2,1,1,l) = 0.
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPNX22(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c j2blok = number of data blocks
c kxp2 = number of data values per block
      complex cu, zero
      dimension cu(2,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*aimag(cu(1,k,j,l)) + dky*aimag(cu(2,k,j,l)))
         at3 = aimag(cu(1,k,j,l)) - dkx*at2
         at4 = aimag(cu(2,k,j,l)) - dky*at2
         cu(1,k,j,l) = cmplx(0.,at3)
         cu(2,k,j,l) = cmplx(0.,at4)
         cu(1,k1,j,l) = cmplx(0.,at3)
         cu(2,k1,j,l) = cmplx(0.,-at4)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISNX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstr
     1t,ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 15*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 8*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j,l) = cmplx(0.,at2*real(cu(3,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,-at3*real(cu(3,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,
     1l)),0.)
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
         bxy(3,k1,j,l) = -bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
   50    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,-at2*real(cu(3,1,j,l)))
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bxy(1,k,1,l) = cmplx(0.,at2*real(cu(3,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
   70    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at2*aimag(cu(2,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*real(cu(3,k,j,l)),0.)
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
         bxy(3,k1,j,l) = bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
  110    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(at2*real(cu(3,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at2*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at2*real(cu(3,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
  130    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at1*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at1*aimag(cu(2,k,j,l)))
         bxy(3,k,j,l) = cmplx(at1*real(cu(3,k,j,l)),0.)
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
         bxy(3,k1,j,l) = bxy(3,k,j,l)
  170    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(at1*real(cu(3,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at1*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at1*real(cu(3,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
  190    continue
         k1 = ny + 1
         at1 = aimag(ffd(1,1,l))
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = cmplx(at1*real(cu(3,1,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISN23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstrt
     1,nyv,kxp2,j2blok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 15*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 8*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(3,nyv,kxp2,j2blok), bxy(3,nyv,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j,l) = at2*cu(3,k,j,l)
         bxy(2,k,j,l) = -at3*cu(3,k,j,l)
         bxy(3,k,j,l) = at2*cu(1,k,j,l) - at3*cu(2,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = -at2*cu(3,1,j,l)
         bxy(3,1,j,l) = -at2*cu(2,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bxy(1,k,1,l) = at2*cu(3,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at2*cu(1,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
   70    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(3,k,j,l) = at2*cu(3,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
  110    continue
c mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(3,1,j,l) = at2*cu(3,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at2*cu(3,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
  130    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
  170    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(3,1,j,l) = at1*cu(3,1,j,l)
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at1*cu(3,k,1,l)
  190    continue
         at1 = aimag(ffd(1,1,l))
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = at1*cu(3,1,1,l)
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISNX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,k
     1strt,ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bz,wm
c approximate flop count is: 15*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 11*nxc*nyc + 6*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 6*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, bz, ffd, zero
      dimension cu(2,ny2d,kxp2,j2blok), bxy(2,ny2d,kxp2,j2blok)
      dimension bz(ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bz(k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,l))
     1,0.)
         bz(k1,j,l) = -bz(k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12)
   50    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bz(k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bz(k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bz(k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2)
   70    continue
         k1 = ny + 1
         bz(1,1,l) = zero
         bz(k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at2*aimag(cu(2,k,j,l)))
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12)
  110    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at2*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2)
  130    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at1*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at1*aimag(cu(2,k,j,l)))
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
  170    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at1*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
  190    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISNX23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,ny2d,kxp2,j2bl
     1ok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd, output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*,
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*real(cu(3,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,-at3*real(cu(3,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,
     1l)),0.)
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
         bxy(3,k1,j,l) = -bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,-at2*real(cu(3,1,j,l)))
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at2*real(cu(3,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
   30    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wm = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blok
     1,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with neumann boundary conditions (zero normal efield),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd, output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(3,nyv,kxp2,j2blok), bxy(3,nyv,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at2*cu(3,k,j,l)
         bxy(2,k,j,l) = -at3*cu(3,k,j,l)
         bxy(3,k,j,l) = at2*cu(1,k,j,l) - at3*cu(2,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
   10    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = -at2*cu(3,1,j,l)
         bxy(3,1,j,l) = -at2*cu(2,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at2*cu(3,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at2*cu(1,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
   30    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
   40 continue
   50 continue
      wm = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELNX2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,
     1ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with neumann boundary
c conditions (zero normal efield).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 62*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp, ws
      complex exy, bxy, cu, ffd
      complex zero
      dimension exy(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension cu(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      if (ci.le.0.) return
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
c calculate the electromagnetic fields
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,j,l))
         at7 = aimag(exy(1,k,j,l))
         at8 = aimag(exy(2,k,j,l))
         at9 = real(exy(3,k,j,l))
c update magnetic field half time step, ky > 0
         at4 = aimag(bxy(1,k,j,l)) - dth*(dky*at9)
         at5 = aimag(bxy(2,k,j,l)) + dth*(dkx*at9)
         at6 = real(bxy(3,k,j,l)) + dth*(dkx*at8 - dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,k,j,l))
         at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,k,j,l))
         at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*real(cu(3,k,j,l))
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8 - dky*at7)
         ws = ws + 2.0*anorm*(at7*at7 + at8*at8 + at9*at9)
         wp = wp + 2.0*anorm*(at4*at4 + at5*at5 + at6*at6)
         exy(1,k,j,l) = cmplx(0.,at7)
         exy(2,k,j,l) = cmplx(0.,at8)
         exy(3,k,j,l) = cmplx(at9,0.)
         bxy(1,k,j,l) = cmplx(0.,at4)
         bxy(2,k,j,l) = cmplx(0.,at5)
         bxy(3,k,j,l) = cmplx(at6,0.)
c update electric and magnetic fields, ky < 0
         exy(1,k1,j,l) = cmplx(0.,at7)
         exy(2,k1,j,l) = cmplx(0.,-at8)
         exy(3,k1,j,l) = cmplx(at9,0.)
         bxy(1,k1,j,l) = cmplx(0.,-at4)
         bxy(2,k1,j,l) = cmplx(0.,at5)
         bxy(3,k1,j,l) = cmplx(-at6,0.)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         afdt = adt*aimag(ffd(1,j,l))
         at8 = aimag(exy(2,1,j,l))
         at9 = real(exy(3,1,j,l))
c update magnetic field half time step, ky > 0
         at5 = aimag(bxy(2,1,j,l)) + dth*(dkx*at9)
         at6 = real(bxy(3,1,j,l)) + dth*(dkx*at8)
c update electric field whole time step
         at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,1,j,l))
         at9 = at9 - cdt*(dkx*at5) - afdt*real(cu(3,1,j,l))
c update magnetic field half time step and store electric field
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = zero
         exy(2,1,j,l) = cmplx(0.,at8)
         exy(3,1,j,l) = cmplx(at9,0.)
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at5)
         bxy(3,1,j,l) = cmplx(at6,0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,1,l))
         at7 = aimag(exy(1,k,1,l))
         at9 = real(exy(3,k,1,l))
c update magnetic field half time step, ky > 0
         at4 = aimag(bxy(1,k,1,l)) - dth*(dky*at9)
         at6 = real(bxy(3,k,1,l)) - dth*(dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,k,1,l))
         at9 = at9 + cdt*(dky*at4) - afdt*real(cu(3,k,1,l))
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at6 = at6 - dth*(dky*at7)
         ws = ws + anorm*(at7*at7 + at9*at9)
         wp = wp + anorm*(at4*at4 + at6*at6)
         exy(1,k,1,l) = cmplx(0.,at7)
         exy(2,k,1,l) = zero
         exy(3,k,1,l) = cmplx(at9,0.)
         bxy(1,k,1,l) = cmplx(0.,at4)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at6,0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELN2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,n
     1yv,kxp2,j2blok,nyd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with neumann boundary
c conditions (zero normal efield).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 62*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp, ws
      complex ffd
      dimension exy(3,nyv,kxp2,j2blok), bxy(3,nyv,kxp2,j2blok)
      dimension cu(3,nyv,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      if (ci.le.0.) return
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
c calculate the electromagnetic fields
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,j,l))
         at7 = exy(1,k,j,l)
         at8 = exy(2,k,j,l)
         at9 = exy(3,k,j,l)
c update magnetic field half time step, ky > 0
         at4 = bxy(1,k,j,l) - dth*(dky*at9)
         at5 = bxy(2,k,j,l) + dth*(dkx*at9)
         at6 = bxy(3,k,j,l) + dth*(dkx*at8 - dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*cu(1,k,j,l)
         at8 = at8 - cdt*(dkx*at6) - afdt*cu(2,k,j,l)
         at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,k,j,l)
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8 - dky*at7)
         ws = ws + 2.0*anorm*(at7*at7 + at8*at8 + at9*at9)
         wp = wp + 2.0*anorm*(at4*at4 + at5*at5 + at6*at6)
         exy(1,k,j,l) = at7
         exy(2,k,j,l) = at8
         exy(3,k,j,l) = at9
         bxy(1,k,j,l) = at4
         bxy(2,k,j,l) = at5
         bxy(3,k,j,l) = at6
   10    continue
c mode numbers ky = 0, ny
         afdt = adt*aimag(ffd(1,j,l))
         at8 = exy(2,1,j,l)
         at9 = exy(3,1,j,l)
c update magnetic field half time step, ky > 0
         at5 = bxy(2,1,j,l) + dth*(dkx*at9)
         at6 = bxy(3,1,j,l) + dth*(dkx*at8)
c update electric field whole time step
         at8 = at8 - cdt*(dkx*at6) - afdt*cu(2,1,j,l)
         at9 = at9 - cdt*(dkx*at5) - afdt*cu(3,1,j,l)
c update magnetic field half time step and store electric field
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = 0.
         exy(2,1,j,l) = at8
         exy(3,1,j,l) = at9
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at5
         bxy(3,1,j,l) = at6
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,1,l))
         at7 = exy(1,k,1,l)
         at9 = exy(3,k,1,l)
c update magnetic field half time step, ky > 0
         at4 = bxy(1,k,1,l) - dth*(dky*at9)
         at6 = bxy(3,k,1,l) - dth*(dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*cu(1,k,1,l)
         at9 = at9 + cdt*(dky*at4) - afdt*cu(3,k,1,l)
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at6 = at6 - dth*(dky*at7)
         ws = ws + anorm*(at7*at7 + at9*at9)
         wp = wp + anorm*(at4*at4 + at6*at6)
         exy(1,k,1,l) = at7
         exy(2,k,1,l) = 0.
         exy(3,k,1,l) = at9
         bxy(1,k,1,l) = at4
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at6
   30    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
         exy(1,1,1,l) = 0.
         exy(2,1,1,l) = 0.
         exy(3,1,1,l) = 0.
      endif
   40 continue
   50 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDMFIELDN2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies the charge density into a smaller array
      implicit none
      integer nx, ny, kstrt, nyv, ny2d, kxp2, j2blok
      complex q2
      real q
      dimension q2(ny2d,kxp2,j2blok), q(nyv,kxp2,j2blok)
      integer j, k, l
      if (kstrt.gt.nx) return
      do 30 l = 1, j2blok
      do 20 j = 1, kxp2
      do 10 k = 1, ny
      q(k,j,l) = real(q2(k,j,l))
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCMFIELDN2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies the current into a smaller array
      implicit none
      integer nx, ny, kstrt, nyv, ny2d, kxp2, j2blok
      complex cu2
      real cu
      dimension cu2(3,ny2d,kxp2,j2blok), cu(3,nyv,kxp2,j2blok)
      integer j, k, l
      if (kstrt.gt.nx) return
      do 30 l = 1, j2blok
      do 20 j = 1, kxp2
      do 10 k = 1, ny
      cu(1,k,j,l) = aimag(cu2(1,k,j,l))
      cu(2,k,j,l) = aimag(cu2(2,k,j,l))
      cu(3,k,j,l) = real(cu2(3,k,j,l))
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELDN2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,
     1j2blok,nyd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
      complex fxy, ffd
      real exy
      dimension fxy(3,ny2d,kxp2,j2blok), exy(3,nyv,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      integer j, k, l, ny2, ks, joff, k1
      real at1
      complex zero
      ny2 = 2*ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
c add the fields
      if (isign.gt.0) then
         do 40 l = 1, j2blok
         do 20 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 10 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,j,l))
            fxy(1,k,j,l) = fxy(1,k,j,l) + cmplx(0.,exy(1,k,j,l)*at1)
            fxy(2,k,j,l) = fxy(2,k,j,l) + cmplx(0.,exy(2,k,j,l)*at1)
            fxy(3,k,j,l) = fxy(3,k,j,l) + cmplx(exy(3,k,j,l)*at1,0.)
            fxy(1,k1,j,l) = fxy(1,k1,j,l) + cmplx(0.,exy(1,k,j,l)*at1)
            fxy(2,k1,j,l) = fxy(2,k1,j,l) + cmplx(0.,-exy(2,k,j,l)*at1)
            fxy(3,k1,j,l) = fxy(3,k1,j,l) + cmplx(exy(3,k,j,l)*at1,0.)
   10       continue
            at1 = aimag(ffd(1,j,l))
            fxy(1,1,j,l) = fxy(1,1,j,l) + cmplx(0.,exy(1,1,j,l)*at1)
            fxy(2,1,j,l) = fxy(2,1,j,l) + cmplx(0.,exy(2,1,j,l)*at1)
            fxy(3,1,j,l) = fxy(3,1,j,l) + cmplx(exy(3,1,j,l)*at1,0.)
         endif
   20    continue
         if ((l+ks).eq.0) then
            do 30 k = 1, ny
            at1 = aimag(ffd(k,1,l))
            fxy(1,k,1,l) = fxy(1,k,1,l) + cmplx(0.,exy(1,k,1,l)*at1)
            fxy(2,k,1,l) = fxy(2,k,1,l) + cmplx(0.,exy(2,k,1,l)*at1)
            fxy(3,k,1,l) = fxy(3,k,1,l) + cmplx(exy(3,k,1,l)*at1,0.)
   30       continue
         endif
   40    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 80 l = 1, j2blok
         do 60 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 50 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,j,l))
            fxy(1,k,j,l) = cmplx(0.,exy(1,k,j,l)*at1)
            fxy(2,k,j,l) = cmplx(0.,exy(2,k,j,l)*at1)
            fxy(3,k,j,l) = cmplx(exy(3,k,j,l)*at1,0.)
            fxy(1,k1,j,l) = cmplx(0.,-exy(1,k,j,l)*at1)
            fxy(2,k1,j,l) = cmplx(0.,exy(2,k,j,l)*at1)
            fxy(3,k1,j,l) = cmplx(-exy(3,k,j,l)*at1,0.)
   50       continue
            k1 = ny + 1
            at1 = aimag(ffd(1,j,l))
            fxy(1,1,j,l) = cmplx(0.,exy(1,1,j,l)*at1)
            fxy(2,1,j,l) = cmplx(0.,exy(2,1,j,l)*at1)
            fxy(3,1,j,l) = cmplx(exy(3,1,j,l)*at1,0.)
            fxy(1,k1,j,l) = zero
            fxy(2,k1,j,l) = zero
            fxy(3,k1,j,l) = zero
         endif
   60    continue
         if ((l+ks).eq.0) then
            do 70 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,1,l))
            fxy(1,k,1,l) = cmplx(0.,exy(1,k,1,l)*at1)
            fxy(2,k,1,l) = cmplx(0.,exy(2,k,1,l)*at1)
            fxy(3,k,1,l) = cmplx(exy(3,k,1,l)*at1,0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
   70       continue
            k1 = ny + 1
            at1 = aimag(ffd(1,1,l))
            fxy(1,1,1,l) = cmplx(0.,exy(1,1,1,l)*at1)
            fxy(2,1,1,l) = cmplx(0.,exy(2,1,1,l)*at1)
            fxy(3,1,1,l) = cmplx(exy(3,1,1,l)*at1,0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
         endif
   80    continue
c copy the electric fields
      else
         do 120 l = 1, j2blok
         do 100 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 90 k = 2, ny
            k1 = ny2 - k
            fxy(1,k,j,l) = cmplx(0.,exy(1,k,j,l))
            fxy(2,k,j,l) = cmplx(0.,exy(2,k,j,l))
            fxy(3,k,j,l) = cmplx(exy(3,k,j,l),0.)
            fxy(1,k1,j,l) = cmplx(0.,exy(1,k,j,l))
            fxy(2,k1,j,l) = cmplx(0.,-exy(2,k,j,l))
            fxy(3,k1,j,l) = cmplx(exy(3,k,j,l),0.)
   90       continue
            k1 = ny + 1
            fxy(1,1,j,l) = cmplx(0.,exy(1,1,j,l))
            fxy(2,1,j,l) = cmplx(0.,exy(2,1,j,l))
            fxy(3,1,j,l) = cmplx(exy(3,1,j,l),0.)
            fxy(1,k1,j,l) = zero
            fxy(2,k1,j,l) = zero
            fxy(3,k1,j,l) = zero
         endif
  100    continue
         if ((l+ks).eq.0) then
            do 110 k = 2, ny
            k1 = ny2 - k
            fxy(1,k,1,l) = cmplx(0.,exy(1,k,1,l))
            fxy(2,k,1,l) = cmplx(0.,exy(2,k,1,l))
            fxy(3,k,1,l) = cmplx(exy(3,k,1,l),0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
  110       continue
            k1 = ny + 1
            fxy(1,1,1,l) = cmplx(0.,exy(1,1,1,l))
            fxy(2,1,1,l) = cmplx(0.,exy(2,1,1,l))
            fxy(3,1,1,l) = cmplx(exy(3,1,1,l),0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
         endif
  120    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMFIELDN2(pot2,pot,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
      complex pot2
      real pot
      dimension pot2(ny2d,kxp2,j2blok), pot(nyv,kxp2,j2blok)
      integer j, k, l, ny2, ks, joff, k1
      complex zero
      ny2 = 2*ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
      do 20 j = 1, kxp2
      joff = kxp2*(l + ks) - 1
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         pot2(k,j,l) = cmplx(pot(k,j,l),0.)
         pot2(k1,j,l) = cmplx(-pot(k,j,l),0.)
   10    continue
         k1 = ny + 1
         pot2(1,j,l) = cmplx(pot(1,j,l),0.)
         pot2(k1,j,l) = zero
      endif
   20 continue
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         pot2(k,1,l) = cmplx(pot(k,1,l),0.)
         pot2(k1,1,l) = zero
   30    continue
         k1 = ny + 1
         pot2(1,1,l) = cmplx(pot(1,1,l),0.)
         pot2(k1,1,l) = zero
         endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCPFIELDN2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
      complex fxy
      real exy
      dimension fxy(3,ny2d,kxp2,j2blok), exy(3,nyv,kxp2,j2blok)
c local data
      integer isign
      complex ffd
      dimension ffd(1,1,1)
      isign = 0
      call PEMFIELDD2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok
     1,nyd)
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTNX23(bxy,axy,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with neumann boundary conditions (zero normal efield)
c for distributed data.
c input: bxy,nx,ny,kstrt,ny2d,kxp2,j2blok, output: axy
c approximate flop count is: 14*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
      complex bxy, axy, zero
      dimension bxy(3,ny2d,kxp2,j2blok), axy(3,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         at4 = real(bxy(3,k,j,l))
         at5 = aimag(bxy(2,k,j,l))
         at6 = aimag(bxy(1,k,j,l))
         at6 = at2*at6 - at3*at5
         at5 = -at3*at4
         at4 = at2*at4
         axy(1,k,j,l) = cmplx(0.,at4)
         axy(2,k,j,l) = cmplx(0.,at5)
         axy(3,k,j,l) = cmplx(at6,0.)
         axy(1,k1,j,l) = cmplx(0.,at4)
         axy(2,k1,j,l) = cmplx(0.,-at5)
         axy(3,k1,j,l) = cmplx(at6,0.)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = ny + 1
         at2 = 1.0/dkx
         at4 = real(bxy(3,1,j,l))
         at5 = aimag(bxy(2,1,j,l))
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = cmplx(0.,-at2*at4)
         axy(3,1,j,l) = cmplx(-at2*at5,0.)
         axy(1,k1,j,l) = zero
         axy(2,k1,j,l) = zero
         axy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         at4 = real(bxy(3,k,1,l))
         at6 = aimag(bxy(1,k,1,l))
         axy(1,k,1,l) = cmplx(0.,at2*at4)
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = cmplx(at2*at6,0.)
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         axy(1,1,1,l) = zero
         axy(2,1,1,l) = zero
         axy(3,1,1,l) = zero
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTN23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with neumann boundary conditions (zero normal efield)
c for distributed data.
c input: bxy,nx,ny,kstrt,ny2d,kxp2,j2blok, output: axy
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
      dimension bxy(3,nyv,kxp2,j2blok), axy(3,nyv,kxp2,j2blok)
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate vector potential
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         at4 = bxy(3,k,j,l)
         at5 = bxy(2,k,j,l)
         at6 = bxy(1,k,j,l)
         axy(1,k,j,l) = at2*at4
         axy(2,k,j,l) = -at3*at4
         axy(3,k,j,l) = at2*at6 - at3*at5
   10    continue
c mode numbers ky = 0, ny/2
         at2 = 1.0/dkx
         at4 = bxy(3,1,j,l)
         at5 = bxy(2,1,j,l)
         axy(1,1,j,l) = 0.
         axy(2,1,j,l) = -at2*at4
         axy(3,1,j,l) = -at2*at5
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         at4 = bxy(3,k,1,l)
         at6 = bxy(1,k,1,l)
         axy(1,k,1,l) = at2*at4
         axy(2,k,1,l) = 0.
         axy(3,k,1,l) = at2*at6
   30    continue
         axy(1,1,1,l) = 0.
         axy(2,1,1,l) = 0.
         axy(3,1,1,l) = 0.
      endif
   40 continue
   50 continue
      return
      end
