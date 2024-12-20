c 2d parallel PIC library for fast fourier transforms
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: june 16, 2008
c-----------------------------------------------------------------------
      subroutine PFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstrt
     1,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs a two dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(n,m,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(m,n,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, a forward fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(j,k,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(j,k,i) = sum(g(m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c bs, br = scratch arrays
c kstrt = starting data block number
c nxvh/nyv = first dimension of f/g
c kypd = second dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension bs(kxp,kyp,kblok), br(kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2
      real dnxy, arg, ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 300
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 180
      nrx = nxhy/nxh
      do 80 l = 1, kblok
c bit-reverse array elements in x
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
   60 continue
   70 continue
   80 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 120 l = 1, kblok
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 90 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 170 l = 1, kblok
      do 150 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 140 k = 1, kyp
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = ani*(s + t)
      f(nxh2-j,k,l) = ani*conjg(s - t)
  140 continue
  150 continue
      do 160 k = 1, kyp
      f(nxhh+1,k,l) = 2.*ani*conjg(f(nxhh+1,k,l))
      f(1,k,l) = 2.*ani*cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,
     1k,l)) - aimag(f(1,k,l)))
  160 continue
  170 continue
c transpose f array to g
  180 call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok
     1,kblok)
      if (kstrt.gt.nxh) go to 290
      nry = nxhy/ny
      do 210 l = 1, jblok
c bit-reverse array elements in y
      do 200 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 200
      do 190 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  190 continue
  200 continue
  210 continue
c then transform in y
      nry = nxy/ny
      do 260 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 250 l = 1, jblok
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 220 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
c unscramble modes kx = 0, nx/2
      do 280 l = 1, jblok
      if ((l+ks).gt.0) go to 280
      do 270 k = 2, nyh
      s = g(ny2-k,1,l)
      g(ny2-k,1,l) = .5*cmplx(aimag(g(k,1,l) + s),real(g(k,1,l) - s))
      g(k,1,l) = .5*cmplx(real(g(k,1,l) + s),aimag(g(k,1,l) - s))
  270 continue
  280 continue
c transpose g array to f
  290 if (ntpose.eq.0) call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,k
     1xp,kypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  300 if (ntpose.eq.0) call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,k
     1yp,kxp,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 410
      nry = nxhy/ny
      do 350 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 320
      do 310 k = 2, nyh
      s = cmplx(aimag(g(ny2-k,1,l)),real(g(ny2-k,1,l)))
      g(ny2-k,1,l) = conjg(g(k,1,l) - s)
      g(k,1,l) = g(k,1,l) + s
  310 continue
c bit-reverse array elements in y
  320 do 340 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 340
      do 330 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  330 continue
  340 continue
  350 continue
c first transform in y
      nry = nxy/ny
      do 400 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 390 l = 1, jblok
      do 380 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 370 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 360 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  360 continue
  370 continue
  380 continue
  390 continue
  400 continue
c transpose g array to f
  410 call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok
     1,jblok)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 470 l = 1, kblok
c scramble coefficients
      do 430 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 420 k = 1, kyp
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = s + t
      f(nxh2-j,k,l) = conjg(s - t)
  420 continue
  430 continue
      do 440 k = 1, kyp
      f(nxhh+1,k,l) = 2.*conjg(f(nxhh+1,k,l))
      f(1,k,l) = cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,k,l)) -
     1 aimag(f(1,k,l)))
  440 continue
c bit-reverse array elements in x
      do 460 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 460
      do 450 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
  450 continue
  460 continue
  470 continue
c then transform in x
      nrx = nxy/nxh
      do 520 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 510 l = 1, kblok
      do 500 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 490 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 480 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs 2 two dimensional complex to real fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, two inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:2,n,m,l) = (1/nx*ny)*sum(f(1:2,j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(1:2,m,n,l) = (1/nx*ny)*sum(f(1:2,j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, two forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:2,j,k,i) = sum(f(1:2,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(1:2,j,k,i) = sum(g(1:2,m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c bs, br = scratch arrays
c kstrt = starting data block number
c nxvh/nyv = second dimension of f/g
c kypd = third dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      complex f, g, bs, br, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension bs(2,kxp,kyp,kblok), br(2,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real dnxy, arg, ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 360
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 230
c swap complex components
      do 80 l = 1, kblok
      do 70 k = 1, kyp
      do 60 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
   60 continue
   70 continue
   80 continue
      nrx = nxhy/nxh
      do 110 l = 1, kblok
c bit-reverse array elements in x
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
   90 continue
  100 continue
  110 continue
c first transform in x
      nrx = nxy/nxh
      do 160 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 150 l = 1, kblok
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 120 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 220 l = 1, kblok
      do 190 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 180 k = 1, kyp
      do 170 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = 1, kyp
      do 200 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  200 continue
  210 continue
  220 continue
c transpose f array to g
  230 call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblo
     1k,kblok)
      if (kstrt.gt.nxh) go to 350
      nry = nxhy/ny
      do 260 l = 1, jblok
c bit-reverse array elements in y
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = 1, kxp
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
      nry = nxy/ny
      do 310 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 l = 1, jblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 270 i = 1, kxp
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
c unscramble modes kx = 0, nx/2
      do 340 l = 1, jblok
      if ((l+ks).gt.0) go to 340
      do 330 k = 2, nyh
      do 320 jj = 1, 2
      s = g(jj,ny2-k,1,l)
      g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1,l)
     1- s))
      g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) - s
     1))
  320 continue
  330 continue
  340 continue
c transpose g array to f
  350 if (ntpose.eq.0) call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,
     1kxp,kypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  360 if (ntpose.eq.0) call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,
     1kyp,kxp,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 480
      nry = nxhy/ny
      do 420 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 390
      do 380 k = 2, nyh
      do 370 jj = 1, 2
      s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
      g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
      g(jj,k,1,l) = g(jj,k,1,l) + s
  370 continue
  380 continue
c bit-reverse array elements in y
  390 do 410 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 410
      do 400 j = 1, kxp
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
  400 continue
  410 continue
  420 continue
c first transform in y
      nry = nxy/ny
      do 470 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 460 l = 1, jblok
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, kxp
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
c transpose g array to f
  480 call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblo
     1k,jblok)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 560 l = 1, kblok
c scramble coefficients
      do 510 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 500 k = 1, kyp
      do 490 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  490 continue
  500 continue
  510 continue
      do 530 k = 1, kyp
      do 520 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  520 continue
  530 continue
c bit-reverse array elements in x
      do 550 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 550
      do 540 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
  540 continue
  550 continue
  560 continue
c then transform in x
      nrx = nxy/nxh
      do 610 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 600 l = 1, kblok
      do 590 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 580 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 570 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
  570 continue
  580 continue
  590 continue
  600 continue
  610 continue
c swap complex components
      do 640 l = 1, kblok
      do 630 k = 1, kyp
      do 620 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
  620 continue
  630 continue
  640 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs 3 two dimensional complex to real fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, three inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:3,n,m,l) = (1/nx*ny)*sum(f(1:3,j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(1:3,m,n,l) = (1/nx*ny)*sum(f(1:3,j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, three forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:3,j,k,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(1:3,j,k,i) = sum(g(1:3,m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c bs, br = scratch arrays
c kstrt = starting data block number
c nxvh/nyv = second dimension of f/g
c kypd = third dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension bs(3,kxp,kyp,kblok), br(3,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real dnxy, arg, ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 360
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 230
c swap complex components
      do 80 l = 1, kblok
      do 70 i = 1, kyp
      do 60 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(real(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = aimag(f(2,j,i,l))
      f(2,j,i,l) = cmplx(aimag(f(1,j,i,l)),at1)
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
   60 continue
   70 continue
   80 continue
      nrx = nxhy/nxh
      do 110 l = 1, kblok
c bit-reverse array elements in x
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
   90 continue
  100 continue
  110 continue
c first transform in x
      nrx = nxy/nxh
      do 160 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 150 l = 1, kblok
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 120 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 220 l = 1, kblok
      do 190 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 180 k = 1, kyp
      do 170 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = 1, kyp
      do 200 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  200 continue
  210 continue
  220 continue
c transpose f array to g
  230 call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblo
     1k,kblok)
      if (kstrt.gt.nxh) go to 350
      nry = nxhy/ny
      do 260 l = 1, jblok
c bit-reverse array elements in y
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = 1, kxp
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
      nry = nxy/ny
      do 310 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 l = 1, jblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 270 i = 1, kxp
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
c unscramble modes kx = 0, nx/2
      do 340 l = 1, jblok
      if ((l+ks).gt.0) go to 340
      do 330 k = 2, nyh
      do 320 jj = 1, 3
      s = g(jj,ny2-k,1,l)
      g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1,l)
     1- s))
      g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) - s
     1))
  320 continue
  330 continue
  340 continue
c transpose g array to f
  350 if (ntpose.eq.0) call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,
     1kxp,kypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  360 if (ntpose.eq.0) call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,
     1kyp,kxp,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 480
      nry = nxhy/ny
      do 420 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 390
      do 380 k = 2, nyh
      do 370 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
      g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
      g(jj,k,1,l) = g(jj,k,1,l) + s
  370 continue
  380 continue
c bit-reverse array elements in y
  390 do 410 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 410
      do 400 j = 1, kxp
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
  400 continue
  410 continue
  420 continue
c first transform in y
      nry = nxy/ny
      do 470 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 460 l = 1, jblok
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, kxp
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
c transpose g array to f
  480 call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblo
     1k,jblok)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 560 l = 1, kblok
c scramble coefficients
      do 510 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 500 k = 1, kyp
      do 490 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  490 continue
  500 continue
  510 continue
      do 530 k = 1, kyp
      do 520 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  520 continue
  530 continue
c bit-reverse array elements in x
      do 550 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 550
      do 540 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
  540 continue
  550 continue
  560 continue
c then transform in x
      nrx = nxy/nxh
      do 610 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 600 l = 1, kblok
      do 590 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 580 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 570 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
  570 continue
  580 continue
  590 continue
  600 continue
  610 continue
c swap complex components
      do 640 l = 1, kblok
      do 630 i = 1, kyp
      do 620 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(aimag(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = real(f(2,j,i,l))
      f(2,j,i,l) = cmplx(at1,aimag(f(1,j,i,l)))
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
  620 continue
  630 continue
  640 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RX(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,nxvh
     1,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs a two dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(n,m,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(m,n,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, a forward fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(j,k,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(j,k,i) = sum(g(m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c kstrt = starting data block number
c nxvh/nyv = first dimension of f/g
c kypd = second dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      complex f, g, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2
      real dnxy, arg, ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 300
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 180
      nrx = nxhy/nxh
      do 80 l = 1, kblok
c bit-reverse array elements in x
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
   60 continue
   70 continue
   80 continue
c first transform in x
      nrx = nxy/nxh
      do 130 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 120 l = 1, kblok
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 90 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 170 l = 1, kblok
      do 150 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 140 k = 1, kyp
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = ani*(s + t)
      f(nxh2-j,k,l) = ani*conjg(s - t)
  140 continue
  150 continue
      do 160 k = 1, kyp
      f(nxhh+1,k,l) = 2.*ani*conjg(f(nxhh+1,k,l))
      f(1,k,l) = 2.*ani*cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,
     1k,l)) - aimag(f(1,k,l)))
  160 continue
  170 continue
c transpose f array to g
  180 call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,kblo
     1k)
      if (kstrt.gt.nxh) go to 290
      nry = nxhy/ny
      do 210 l = 1, jblok
c bit-reverse array elements in y
      do 200 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 200
      do 190 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  190 continue
  200 continue
  210 continue
c then transform in y
      nry = nxy/ny
      do 260 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 250 l = 1, jblok
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 220 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
c unscramble modes kx = 0, nx/2
      do 280 l = 1, jblok
      if ((l+ks).gt.0) go to 280
      do 270 k = 2, nyh
      s = g(ny2-k,1,l)
      g(ny2-k,1,l) = .5*cmplx(aimag(g(k,1,l) + s),real(g(k,1,l) - s))
      g(k,1,l) = .5*cmplx(real(g(k,1,l) + s),aimag(g(k,1,l) - s))
  270 continue
  280 continue
c transpose g array to f
  290 if (ntpose.eq.0) call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,ky
     1pd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  300 if (ntpose.eq.0) call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kx
     1p,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 410
      nry = nxhy/ny
      do 350 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 320
      do 310 k = 2, nyh
      s = cmplx(aimag(g(ny2-k,1,l)),real(g(ny2-k,1,l)))
      g(ny2-k,1,l) = conjg(g(k,1,l) - s)
      g(k,1,l) = g(k,1,l) + s
  310 continue
c bit-reverse array elements in y
  320 do 340 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 340
      do 330 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  330 continue
  340 continue
  350 continue
c first transform in y
      nry = nxy/ny
      do 400 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 390 l = 1, jblok
      do 380 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 370 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 360 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  360 continue
  370 continue
  380 continue
  390 continue
  400 continue
c transpose g array to f
  410 call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,jblo
     1k)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 470 l = 1, kblok
c scramble coefficients
      do 430 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 420 k = 1, kyp
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = s + t
      f(nxh2-j,k,l) = conjg(s - t)
  420 continue
  430 continue
      do 440 k = 1, kyp
      f(nxhh+1,k,l) = 2.*conjg(f(nxhh+1,k,l))
      f(1,k,l) = cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,k,l)) -
     1 aimag(f(1,k,l)))
  440 continue
c bit-reverse array elements in x
      do 460 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 460
      do 450 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
  450 continue
  460 continue
  470 continue
c then transform in x
      nrx = nxy/nxh
      do 520 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 510 l = 1, kblok
      do 500 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 490 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 480 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RX2(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,nxv
     1h,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs 2 two dimensional complex to real fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, two inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:2,n,m,l) = (1/nx*ny)*sum(f(1:2,j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(1:2,m,n,l) = (1/nx*ny)*sum(f(1:2,j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, two forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:2,j,k,i) = sum(f(1:2,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(1:2,j,k,i) = sum(g(1:2,m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c kstrt = starting data block number
c nxvh/nyv = second dimension of f/g
c kypd = third dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, optimized version
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      complex f, g, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real dnxy, arg, ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 360
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 230
c swap complex components
      do 80 l = 1, kblok
      do 70 k = 1, kyp
      do 60 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
   60 continue
   70 continue
   80 continue
      nrx = nxhy/nxh
      do 110 l = 1, kblok
c bit-reverse array elements in x
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
   90 continue
  100 continue
  110 continue
c first transform in x
      nrx = nxy/nxh
      do 160 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 150 l = 1, kblok
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 120 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 220 l = 1, kblok
      do 190 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 180 k = 1, kyp
      do 170 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = 1, kyp
      do 200 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  200 continue
  210 continue
  220 continue
c transpose f array to g
  230 call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,kbl
     1ok)
      if (kstrt.gt.nxh) go to 350
      nry = nxhy/ny
      do 260 l = 1, jblok
c bit-reverse array elements in y
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = 1, kxp
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
      nry = nxy/ny
      do 310 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 l = 1, jblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 270 i = 1, kxp
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
c unscramble modes kx = 0, nx/2
      do 340 l = 1, jblok
      if ((l+ks).gt.0) go to 340
      do 330 k = 2, nyh
      do 320 jj = 1, 2
      s = g(jj,ny2-k,1,l)
      g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1,l)
     1- s))
      g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) - s
     1))
  320 continue
  330 continue
  340 continue
c transpose g array to f
  350 if (ntpose.eq.0) call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,k
     1ypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  360 if (ntpose.eq.0) call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,k
     1xp,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 480
      nry = nxhy/ny
      do 420 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 390
      do 380 k = 2, nyh
      do 370 jj = 1, 2
      s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
      g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
      g(jj,k,1,l) = g(jj,k,1,l) + s
  370 continue
  380 continue
c bit-reverse array elements in y
  390 do 410 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 410
      do 400 j = 1, kxp
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
  400 continue
  410 continue
  420 continue
c first transform in y
      nry = nxy/ny
      do 470 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 460 l = 1, jblok
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, kxp
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
c transpose g array to f
  480 call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,jbl
     1ok)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 560 l = 1, kblok
c scramble coefficients
      do 510 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 500 k = 1, kyp
      do 490 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  490 continue
  500 continue
  510 continue
      do 530 k = 1, kyp
      do 520 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  520 continue
  530 continue
c bit-reverse array elements in x
      do 550 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 550
      do 540 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
  540 continue
  550 continue
  560 continue
c then transform in x
      nrx = nxy/nxh
      do 610 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 600 l = 1, kblok
      do 590 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 580 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 570 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
  570 continue
  580 continue
  590 continue
  600 continue
  610 continue
c swap complex components
      do 640 l = 1, kblok
      do 630 k = 1, kyp
      do 620 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
  620 continue
  630 continue
  640 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RX3(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,nxv
     1h,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c this subroutine performs 3 two dimensional complex to real fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxhyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, three inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:3,n,m,l) = (1/nx*ny)*sum(f(1:3,j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(1:3,m,n,l) = (1/nx*ny)*sum(f(1:3,j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, three forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(1:3,j,k,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(1:3,j,k,i) = sum(g(1:3,m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c kstrt = starting data block number
c nxvh/nyv = second dimension of f/g
c kypd = third dimension of f
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c if ntpose = 1,
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, optimized version
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      complex f, g, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nxy
      integer nxhy, ks, j, k, lb, ll, jb, it, nxyh, nrx, nry, l, i, m
      integer ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real dnxy, arg, ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      if (isign) 50, 10, 360
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 230
c swap complex components
      do 80 l = 1, kblok
      do 70 i = 1, kyp
      do 60 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(real(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = aimag(f(2,j,i,l))
      f(2,j,i,l) = cmplx(aimag(f(1,j,i,l)),at1)
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
   60 continue
   70 continue
   80 continue
      nrx = nxhy/nxh
      do 110 l = 1, kblok
c bit-reverse array elements in x
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
   90 continue
  100 continue
  110 continue
c first transform in x
      nrx = nxy/nxh
      do 160 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 150 l = 1, kblok
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 120 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 220 l = 1, kblok
      do 190 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 180 k = 1, kyp
      do 170 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = 1, kyp
      do 200 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  200 continue
  210 continue
  220 continue
c transpose f array to g
  230 call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,kbl
     1ok)
      if (kstrt.gt.nxh) go to 350
      nry = nxhy/ny
      do 260 l = 1, jblok
c bit-reverse array elements in y
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = 1, kxp
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
      nry = nxy/ny
      do 310 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 l = 1, jblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 270 i = 1, kxp
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
c unscramble modes kx = 0, nx/2
      do 340 l = 1, jblok
      if ((l+ks).gt.0) go to 340
      do 330 k = 2, nyh
      do 320 jj = 1, 3
      s = g(jj,ny2-k,1,l)
      g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1,l)
     1- s))
      g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) - s
     1))
  320 continue
  330 continue
  340 continue
c transpose g array to f
  350 if (ntpose.eq.0) call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,k
     1ypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  360 if (ntpose.eq.0) call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,k
     1xp,kypd,jblok,kblok)
      if (kstrt.gt.nxh) go to 480
      nry = nxhy/ny
      do 420 l = 1, jblok
c scramble modes kx = 0, nx/2
      if ((l+ks).gt.0) go to 390
      do 380 k = 2, nyh
      do 370 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
      g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
      g(jj,k,1,l) = g(jj,k,1,l) + s
  370 continue
  380 continue
c bit-reverse array elements in y
  390 do 410 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 410
      do 400 j = 1, kxp
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
  400 continue
  410 continue
  420 continue
c first transform in y
      nry = nxy/ny
      do 470 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 460 l = 1, jblok
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, kxp
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
c transpose g array to f
  480 call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,jbl
     1ok)
      if (kstrt.gt.ny) return
      nrx = nxhy/nxh
      kmr = nxy/nx
      do 560 l = 1, kblok
c scramble coefficients
      do 510 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 500 k = 1, kyp
      do 490 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  490 continue
  500 continue
  510 continue
      do 530 k = 1, kyp
      do 520 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  520 continue
  530 continue
c bit-reverse array elements in x
      do 550 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 550
      do 540 k = 1, kyp
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
  540 continue
  550 continue
  560 continue
c then transform in x
      nrx = nxy/nxh
      do 610 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 600 l = 1, kblok
      do 590 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 580 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 570 i = 1, kyp
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
  570 continue
  580 continue
  590 continue
  600 continue
  610 continue
c swap complex components
      do 640 l = 1, kblok
      do 630 i = 1, kyp
      do 620 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(aimag(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = real(f(2,j,i,l))
      f(2,j,i,l) = cmplx(at1,aimag(f(1,j,i,l)))
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
  620 continue
  630 continue
  640 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2C(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstrt
     1,nxv,nyv,kxp,kyp,kypd,jblok,kblok,nxyd,nxyhd)
c this subroutine performs a two dimensional complex to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks
c for isign = 0, input: isign, indx, indy, kstrt, nxyd, nxyhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c approximate flop count: 5*N*log2(N)/nvp
c where N = nx*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(n,m,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, f is the input and g is the output
c g(m,n,l) = (1/nx*ny)*sum(f(j,k,i)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c if isign = 1, a forward fourier transform is performed
c if ntpose = 0, f is the input and output array, g is a scratch array
c f(j,k,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny))
c where mm = m + kyp*(l - 1) and kk = k + kyp*(i - 1)
c if ntpose = 1, g is the input and f is the output
c f(j,k,i) = sum(g(m,n,l)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*m*kk/ny))
c where nn = n + kxp*(l - 1) and kk = k + kyp*(i - 1)
c bs, br = scratch arrays
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd = second dimension of f 
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxy = maximum of (nx,ny)
c nxyh = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxv, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxyd, nxyhd, mixup
      complex f, g, bs, br, sct
      dimension f(nxv,kypd,kblok), g(nyv,kxp,jblok)
      dimension bs(kxp,kyp,kblok), br(kxp,kyp,jblok)
      dimension mixup(nxyd), sct(nxyhd)
c local data
      integer indxy, nx, nxh, ny, nyh, nxy, j, k, lb, ll, jb, it, nxyh
      integer nrx, nry, l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real dnxy, arg, ani
      complex s, t
      indxy = max0(indx,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nxy = 2**indxy
      if (isign) 50, 10, 270
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxy
      lb = j - 1
      ll = 0
      do 20 k = 1, indxy
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 if (kstrt.gt.ny) go to 140
      nrx = nxy/nx
      do 80 l = 1, kblok
c bit-reverse array elements in x
      do 70 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
   60 continue
   70 continue
   80 continue
c first transform in x
      do 130 m = 1, indx
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 120 l = 1, kblok
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 90 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c transpose f array to g
  140 call PTPOSE(f,g,bs,br,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxp,kypd,jblok,k
     1blok)
      if (kstrt.gt.nx) go to 260
      nry = nxy/ny
      do 170 l = 1, jblok
c bit-reverse array elements in y
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 160
      do 150 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  150 continue
  160 continue
  170 continue
c then transform in y
      do 220 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 210 l = 1, jblok
      do 200 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 190 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 180 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  180 continue
  190 continue
  200 continue
  210 continue
  220 continue
c normalize result
      ani = 1./float(nx*ny)
      do 250 l = 1, jblok
      do 240 j = 1, kxp
      do 230 k = 1, ny
      g(k,j,l) = g(k,j,l)*ani
  230 continue
  240 continue
  250 continue
c transpose g array to f
  260 if (ntpose.eq.0) call PTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,nxv,kyp,kxp
     1,kypd,kxp,kblok,jblok)
      return
c forward fourier transform
c transpose f array to g
  270 if (ntpose.eq.0) call PTPOSE(f,g,bs,br,nx,ny,kstrt,nxv,nyv,kxp,kyp
     1,kxp,kypd,jblok,kblok)
      if (kstrt.gt.nx) go to 360
      nry = nxy/ny
      do 300 l = 1, jblok
c bit-reverse array elements in y
      do 290 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 290
      do 280 j = 1, kxp
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  280 continue
  290 continue
  300 continue
c first transform in y
      do 350 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 340 l = 1, jblok
      do 330 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 320 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 310 i = 1, kxp
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
c transpose g array to f
  360 call PTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,nxv,kyp,kxp,kypd,kxp,kblok,j
     1blok)
      if (kstrt.gt.ny) return
      nrx = nxy/nx
      do 390 l = 1, kblok
c bit-reverse array elements in x
      do 380 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 380
      do 370 k = 1, kyp
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
  370 continue
  380 continue
  390 continue
c then transform in x
      do 440 m = 1, indx
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 430 l = 1, kblok
      do 420 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 410 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 400 i = 1, kyp
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
  400 continue
  410 continue
  420 continue
  430 continue
  440 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
c this subroutine calculates tables needed by a two dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, nxhyd, nxyhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy, nxyh
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
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
     1kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension bs(kxp,kyp,kblok), br(kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,k
     1ypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jb
     1lok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,kx
     1p,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp
     1,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd
     1,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,kx
     1p,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kb
     1lok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,k
     1ypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RX(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt
     1,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      real ttp
      complex f, g, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,k
     1ypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,k
     1blok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,kx
     1p,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblo
     1k,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblo
     1k,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,kx
     1p,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,j
     1blok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,k
     1ypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension bs(2,kxp,kyp,kblok), br(2,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension bs(3,kxp,kyp,kblok), br(3,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RX2(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      real ttp
      complex f, g, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RX3(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      real ttp
      complex f, g, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,i
     1ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, ndim, nxhyd, nxyhd, mixup
      real ttp
      complex f, g, bs, br, ss, sct
      dimension f(ndim,nxvh,kypd,kblok), g(ndim,nyv,kxp,jblok)
      dimension bs(ndim,kxp,kyp,kblok), br(ndim,kxp,kyp,jblok)
      dimension ss(ndim,nxvh)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2RNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nx
     1vh,kypd,kblok,ndim,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,ndim,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok,ndim)
            call PWTIMERA(1,tf,dtime)
          endif
c perform y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,ndim,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2RNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nx
     1vh,kypd,kblok,ndim,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT2RXN(f,g,ss,isign,ntpose,mixup,sct,ttp,indx,indy,k
     1strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
c wrapper function for parallel real to complex fft
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, ndim, nxhyd, nxyhd
      real ttp
      complex f, g, ss, sct
      dimension f(ndim,nxvh,kypd,kblok), g(ndim,nyv,kxp,jblok)
      dimension ss(ndim,nxvh)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2RNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nx
     1vh,kypd,kblok,ndim,nxhyd,nxyhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,ndim,nxhyd,nxyhd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,k
     1xp,jblok,ndim,nxhyd,nxyhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2RNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,nx
     1vh,kypd,kblok,ndim,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WP2FFT2RN(f1,f2,g1,g2,bs,br,ss,isign,ntpose,mixup,sct,t
     1tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim1,ndim2,n
     2xhyd,nxyhd)
c wrapper function for two parallel real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, ndim1, ndim2, nxhyd, nxyhd, mixup
      real ttp
      complex f1, f2, g1, g2, bs, br, ss, sct
      dimension f1(ndim1,nxvh,kypd,kblok), f2(ndim2,nxvh,kypd,kblok)
      dimension g1(ndim1,nyv,kxp,jblok), g2(ndim2,nyv,kxp,jblok)
      dimension bs(ndim1+ndim2,kxp,kyp,kblok)
      dimension br(ndim1+ndim2,kxp,kyp,jblok)
      dimension ss(ndim1+ndim2,nxvh)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT2RNXX(f1,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,n
     1xvh,kypd,kblok,ndim1,nxhyd,nxyhd)
         call PFFT2RNXX(f2,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,n
     1xvh,kypd,kblok,ndim2,nxhyd,nxyhd)
c transpose f arrays to g
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOSE(f1,f2,g1,g2,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,k
     1xp,kypd,jblok,kblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT2RNXY(g1,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,
     1kxp,jblok,ndim1,nxhyd,nxyhd)
         call PFFT2RNXY(g2,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,
     1kxp,jblok,ndim2,nxhyd,nxyhd)
c transpose g arrays to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOSE(g1,g2,f1,f2,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kx
     1p,kypd,kxp,kblok,jblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f arrays to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOSE(f1,f2,g1,g2,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,ky
     1p,kxp,kypd,jblok,kblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
          endif
c perform y fft
         call PFFT2RNXY(g1,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,
     1kxp,jblok,ndim1,nxhyd,nxyhd)
         call PFFT2RNXY(g2,isign,mixup,sct,indx,indy,kstrt,kxpi,kxp,nyv,
     1kxp,jblok,ndim2,nxhyd,nxyhd)
c transpose g arrays to f
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOSE(g1,g2,f1,f2,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,k
     1ypd,kxp,kblok,jblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT2RNXX(f1,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,n
     1xvh,kypd,kblok,ndim1,nxhyd,nxyhd)
         call PFFT2RNXX(f2,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp,n
     1xvh,kypd,kblok,ndim2,nxhyd,nxyhd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,nx
     1vh,kypd,kblok,nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = (1/nx*ny)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k,i) = sum(f(n,m,i)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c kblok = number of data blocks in y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, kblok, nxhyd, nxyhd
      complex f, sct
      dimension f(nxvh,kypd,kblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, ks, kypt, j, k, nrx
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      real ani
      complex s, t, t1
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 130
c inverse fourier transform
      ani = 1./float(2*nx*ny)
      nrx = nxhy/nxh
      do 30 l = 1, kblok
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = kypi, kypt
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
   10 continue
   20 continue
   30 continue
c first transform in x
      nrx = nxy/nxh
      do 80 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 l = 1, kblok
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 i = kypi, kypt
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 120 l = 1, kblok
      do 100 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 90 k = kypi, kypt
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = ani*(s + t)
      f(nxh2-j,k,l) = ani*conjg(s - t)
   90 continue
  100 continue
      do 110 k = kypi, kypt
      f(nxhh+1,k,l) = 2.*ani*conjg(f(nxhh+1,k,l))
      f(1,k,l) = 2.*ani*cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,
     1k,l)) - aimag(f(1,k,l)))
  110 continue
  120 continue
      return
c forward fourier transform
  130 kmr = nxy/nx
      do 190 l = 1, kblok
c scramble coefficients
      do 150 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 140 k = kypi, kypt
      t = conjg(f(nxh2-j,k,l))
      s = f(j,k,l) + t
      t = (f(j,k,l) - t)*t1
      f(j,k,l) = s + t
      f(nxh2-j,k,l) = conjg(s - t)
  140 continue
  150 continue
      do 160 k = kypi, kypt
      f(nxhh+1,k,l) = 2.*conjg(f(nxhh+1,k,l))
      f(1,k,l) = cmplx(real(f(1,k,l)) + aimag(f(1,k,l)),real(f(1,k,l)) -
     1 aimag(f(1,k,l)))
  160 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 180 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 180
      do 170 k = kypi, kypt
      t = f(j1,k,l)
      f(j1,k,l) = f(j,k,l)
      f(j,k,l) = t
  170 continue
  180 continue
  190 continue
c then transform in x
      nrx = nxy/nxh
      do 240 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 230 l = 1, kblok
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 200 i = kypi, kypt
      t = s*f(j2,i,l)
      f(j2,i,l) = f(j1,i,l) - t
      f(j1,i,l) = f(j1,i,l) + t
  200 continue
  210 continue
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,ny
     1v,kxp,jblok,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n,i) = sum(g(k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j,i) = sum(g(m,n,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of   indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c jblok = number of data blocks in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, jblok, nxhyd, nxyhd
      complex g, sct
      dimension g(nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2
      complex s, t
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 110
c inverse fourier transform
      nry = nxhy/ny
      do 30 l = 1, jblok
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxy/ny
      do 80 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 l = 1, jblok
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 i = kxpi, kxpt
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 100 l = 1, jblok
      if ((l+ks).gt.0) go to 100
      do 90 k = 2, nyh
      if (kxpi.eq.1) then
         s = g(ny2-k,1,l)
         g(ny2-k,1,l) = .5*cmplx(aimag(g(k,1,l) + s),real(g(k,1,l) - s))
         g(k,1,l) = .5*cmplx(real(g(k,1,l) + s),aimag(g(k,1,l) - s))
      endif
   90 continue
  100 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  110 nry = nxhy/ny
      do 160 l = 1, jblok
      if ((l+ks).gt.0) go to 130
      do 120 k = 2, nyh
      if (kxpi.eq.1) then
         s = cmplx(aimag(g(ny2-k,1,l)),real(g(ny2-k,1,l)))
         g(ny2-k,1,l) = conjg(g(k,1,l) - s)
         g(k,1,l) = g(k,1,l) + s
      endif
  120 continue
c bit-reverse array elements in y
  130 do 150 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 150
      do 140 j = kxpi, kxpt
      t = g(k1,j,l)
      g(k1,j,l) = g(k,j,l)
      g(k,j,l) = t
  140 continue
  150 continue
  160 continue
c first transform in y
      nry = nxy/ny
      do 210 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 200 l = 1, jblok
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 170 i = kxpi, kxpt
      t = s*g(j2,i,l)
      g(j2,i,l) = g(j1,i,l) - t
      g(j1,i,l) = g(j1,i,l) + t
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,n
     1xvh,kypd,kblok,nxhyd,nxyhd)
c this subroutine performs the x part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = (1/nx*ny)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k,i) = sum(f(n,m,i)*exp(sqrt(-1)*2pi*n*j/nx)*
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c kblok = number of data blocks in y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, kblok, nxhyd, nxyhd
      complex f, sct
      dimension f(2,nxvh,kypd,kblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, ks, kypt, j, k, nrx
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real ani, at1
      complex s, t, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 180
c inverse fourier transform
      ani = 1./float(2*nx*ny)
c swap complex components
      do 30 l = 1, kblok
      do 20 k = kypi, kypt
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
   10 continue
   20 continue
   30 continue
      nrx = nxhy/nxh
      do 60 l = 1, kblok
c bit-reverse array elements in x
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 100 l = 1, kblok
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 170 l = 1, kblok
      do 140 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 130 k = kypi, kypt
      do 120 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  120 continue
  130 continue
  140 continue
      do 160 k = kypi, kypt
      do 150 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  150 continue
  160 continue
  170 continue
      return
c forward fourier transform
  180 kmr = nxy/nx
      do 260 l = 1, kblok
c scramble coefficients
      do 210 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 200 k = kypi, kypt
      do 190 jj = 1, 2
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  190 continue
  200 continue
  210 continue
      do 230 k = kypi, kypt
      do 220 jj = 1, 2
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  220 continue
  230 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 250 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 250
      do 240 k = kypi, kypt
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
  240 continue
  250 continue
  260 continue
c then transform in x
      nrx = nxy/nxh
      do 310 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 300 l = 1, kblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 270 i = kypi, kypt
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c swap complex components
      do 340 l = 1, kblok
      do 330 k = kypi, kypt
      do 320 j = 1, nxh
      at1 = aimag(f(1,j,k,l))
      f(1,j,k,l) = cmplx(real(f(1,j,k,l)),real(f(2,j,k,l)))
      f(2,j,k,l) = cmplx(at1,aimag(f(2,j,k,l)))
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,n
     1yv,kxp,jblok,nxhyd,nxyhd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n,i) = sum(g(k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j,i) = sum(g(m,n,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of   indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c jblok = number of data blocks in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, jblok, nxhyd, nxyhd
      complex g, sct
      dimension g(2,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      complex s, t1, t2
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 120
c inverse fourier transform
      nry = nxhy/ny
      do 30 l = 1, jblok
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxy/ny
      do 80 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 l = 1, jblok
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 110 l = 1, jblok
      if ((l+ks).gt.0) go to 110
      do 100 k = 2, nyh
      if (kxpi.eq.1) then
         do 90 jj = 1, 2
         s = g(jj,ny2-k,1,l)
         g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1
     1,l)- s))
         g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) 
     1- s))
   90    continue
      endif
  100 continue
  110 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  120 nry = nxhy/ny
      do 180 l = 1, jblok
      if ((l+ks).gt.0) go to 150
      do 140 k = 2, nyh
      if (kxpi.eq.1) then
         do 130 jj = 1, 2
         s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
         g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
         g(jj,k,1,l) = g(jj,k,1,l) + s
  130    continue
      endif
  140 continue
c bit-reverse array elements in y
  150 do 170 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 170
      do 160 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
  160 continue
  170 continue
  180 continue
c first transform in y
      nry = nxy/ny
      do 230 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 220 l = 1, jblok
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 190 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypp,n
     1xvh,kypd,kblok,nxhyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = (1/nx*ny)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k,i) = sum(f(n,m,i)*exp(sqrt(-1)*2pi*n*j/nx)*
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c kblok = number of data blocks in y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, kblok, nxhyd, nxyhd
      complex f, sct
      dimension f(3,nxvh,kypd,kblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, ks, kypt, j, k, nrx
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 180
c inverse fourier transform
      ani = 1./float(2*nx*ny)
c swap complex components
      do 30 l = 1, kblok
      do 20 i = kypi, kypt
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(real(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = aimag(f(2,j,i,l))
      f(2,j,i,l) = cmplx(aimag(f(1,j,i,l)),at1)
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
   10 continue
   20 continue
   30 continue
      nrx = nxhy/nxh
      do 60 l = 1, kblok
c bit-reverse array elements in x
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = kypi, kypt
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
   40 continue
   50 continue
   60 continue
c first transform in x
      nrx = nxy/nxh
      do 110 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 100 l = 1, kblok
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 70 i = kypi, kypt
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 170 l = 1, kblok
      do 140 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 130 k = kypi, kypt
      do 120 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  120 continue
  130 continue
  140 continue
      do 160 k = kypi, kypt
      do 150 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  150 continue
  160 continue
  170 continue
      return
c forward fourier transform
  180 kmr = nxy/nx
      do 260 l = 1, kblok
c scramble coefficients
      do 210 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 200 k = kypi, kypt
      do 190 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  190 continue
  200 continue
  210 continue
      do 230 k = kypi, kypt
      do 220 jj = 1, 3
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  220 continue
  230 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 250 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 250
      do 240 k = kypi, kypt
      t1 = f(1,j1,k,l)
      t2 = f(2,j1,k,l)
      t3 = f(3,j1,k,l)
      f(1,j1,k,l) = f(1,j,k,l)
      f(2,j1,k,l) = f(2,j,k,l)
      f(3,j1,k,l) = f(3,j,k,l)
      f(1,j,k,l) = t1
      f(2,j,k,l) = t2
      f(3,j,k,l) = t3
  240 continue
  250 continue
  260 continue
c then transform in x
      nrx = nxy/nxh
      do 310 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 300 l = 1, kblok
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 270 i = kypi, kypt
      t1 = s*f(1,j2,i,l)
      t2 = s*f(2,j2,i,l)
      t3 = s*f(3,j2,i,l)
      f(1,j2,i,l) = f(1,j1,i,l) - t1
      f(2,j2,i,l) = f(2,j1,i,l) - t2
      f(3,j2,i,l) = f(3,j1,i,l) - t3
      f(1,j1,i,l) = f(1,j1,i,l) + t1
      f(2,j1,i,l) = f(2,j1,i,l) + t2
      f(3,j1,i,l) = f(3,j1,i,l) + t3
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c swap complex components
      do 340 l = 1, kblok
      do 330 i = kypi, kypt
      do 320 j = 1, nxh
      at1 = real(f(3,j,i,l))
      f(3,j,i,l) = cmplx(aimag(f(2,j,i,l)),aimag(f(3,j,i,l)))
      at2 = real(f(2,j,i,l))
      f(2,j,i,l) = cmplx(at1,aimag(f(1,j,i,l)))
      f(1,j,i,l) = cmplx(real(f(1,j,i,l)),at2)
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,n
     1yv,kxp,jblok,nxhyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 10)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 8)/nvp
c where N = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n,i) = sum(g(k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j,i) = sum(g(m,n,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of   indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c jblok = number of data blocks in x
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, jblok, nxhyd, nxyhd
      complex g, sct
      dimension g(3,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      complex s, t1, t2, t3
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 120
c inverse fourier transform
      nry = nxhy/ny
      do 30 l = 1, jblok
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxy/ny
      do 80 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 l = 1, jblok
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 110 l = 1, jblok
      if ((l+ks).gt.0) go to 110
      do 100 k = 2, nyh
      if (kxpi.eq.1) then
         do 90 jj = 1, 3
         s = g(jj,ny2-k,1,l)
         g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1
     1,l)- s))
         g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) 
     1- s))
   90    continue
      endif
  100 continue
  110 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  120 nry = nxhy/ny
      do 180 l = 1, jblok
      if ((l+ks).gt.0) go to 150
      do 140 k = 2, nyh
      if (kxpi.eq.1) then
         do 130 jj = 1, 3
         s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
         g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
         g(jj,k,1,l) = g(jj,k,1,l) + s
  130    continue
      endif
  140 continue
c bit-reverse array elements in y
  150 do 170 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 170
      do 160 j = kxpi, kxpt
      t1 = g(1,k1,j,l)
      t2 = g(2,k1,j,l)
      t3 = g(3,k1,j,l)
      g(1,k1,j,l) = g(1,k,j,l)
      g(2,k1,j,l) = g(2,k,j,l)
      g(3,k1,j,l) = g(3,k,j,l)
      g(1,k,j,l) = t1
      g(2,k,j,l) = t2
      g(3,k,j,l) = t3
  160 continue
  170 continue
  180 continue
c first transform in y
      nry = nxy/ny
      do 230 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 220 l = 1, jblok
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 190 i = kxpi, kxpt
      t1 = s*g(1,j2,i,l)
      t2 = s*g(2,j2,i,l)
      t3 = s*g(3,j2,i,l)
      g(1,j2,i,l) = g(1,j1,i,l) - t1
      g(2,j2,i,l) = g(2,j1,i,l) - t2
      g(3,j2,i,l) = g(3,j1,i,l) - t3
      g(1,j1,i,l) = g(1,j1,i,l) + t1
      g(2,j1,i,l) = g(2,j1,i,l) + t2
      g(3,j1,i,l) = g(3,j1,i,l) + t3
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RNXX(f,ss,isign,mixup,sct,indx,indy,kstrt,kypi,kyp
     1p,nxvh,kypd,kblok,ndim,nxhyd,nxyhd)
c this subroutine performs the x part of N two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 10)/nvp
c for isign = 1,  approximate flop count: M*(5*log2(M) + 8)/nvp
c where M = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:N,n,m,i) = (1/nx*ny)*sum(f(1:N,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:N,j,k,i) = sum(f(1:N,n,m,i)*exp(sqrt(-1)*2pi*n*j/nx))
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kypd = second dimension of f
c kblok = number of data blocks in y
c ndim = leading dimension of arrays f and g
c ss = scratch array
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,i) = mode j-1,kk-1, where kk = k + kyp*(i - 1)
c 1 <= j <= nx/2 and 1 <= kk <= ny, except for
c f(1,k,i) = mode nx/2,kk-1, where ny/2+2 <= kk <= ny, and
c imaginary part of f(1,1,1) = real part of mode nx/2,0 and
c imaginary part of f(1,1,(ny/2)/kyp+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, nxvh, kypi, kypp
      integer kypd, kblok, ndim, nxhyd, nxyhd
      complex f, ss, sct
      dimension f(ndim,nxvh,kypd,kblok), ss(ndim,nxvh)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny
      integer nxy, nxhy, ks, kypt, j, k, nrx
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kypt = kypi + kypp - 1
      if (kstrt.gt.ny) return
      if (isign.gt.0) go to 170
c inverse fourier transform
      ani = 1./float(2*nx*ny)
c swap complex components
      call PSWAPC2N(f,ss,isign,nxh,kypi,kypt,nxvh,kypd,kblok,ndim)
      nrx = nxhy/nxh
      do 40 l = 1, kblok
c bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 30
      do 20 k = kypi, kypt
      do 10 jj = 1, ndim
      t1 = f(jj,j1,k,l)
      f(jj,j1,k,l) = f(jj,j,k,l)
      f(jj,j,k,l) = t1
   10 continue
   20 continue
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 100 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 90 l = 1, kblok
      do 80 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 70 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 60 i = kypi, kypt
      do 50 jj = 1, ndim
      t1 = s*f(jj,j2,i,l)
      f(jj,j2,i,l) = f(jj,j1,i,l) - t1
      f(jj,j1,i,l) = f(jj,j1,i,l) + t1
   50 continue
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      do 160 l = 1, kblok
      do 130 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 120 k = kypi, kypt
      do 110 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = ani*(s + t)
      f(jj,nxh2-j,k,l) = ani*conjg(s - t)
  110 continue
  120 continue
  130 continue
      do 150 k = kypi, kypt
      do 140 jj = 1, ndim
      f(jj,nxhh+1,k,l) = 2.*ani*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = 2.*ani*cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),
     1real(f(jj,1,k,l)) - aimag(f(jj,1,k,l)))
  140 continue
  150 continue
  160 continue
      return
c forward fourier transform
  170 kmr = nxy/nx
      do 260 l = 1, kblok
c scramble coefficients
      do 200 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 190 k = kypi, kypt
      do 180 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,k,l))
      s = f(jj,j,k,l) + t
      t = (f(jj,j,k,l) - t)*t1
      f(jj,j,k,l) = s + t
      f(jj,nxh2-j,k,l) = conjg(s - t)
  180 continue
  190 continue
  200 continue
      do 220 k = kypi, kypt
      do 210 jj = 1, ndim
      f(jj,nxhh+1,k,l) = 2.*conjg(f(jj,nxhh+1,k,l))
      f(jj,1,k,l) = cmplx(real(f(jj,1,k,l)) + aimag(f(jj,1,k,l)),real(f(
     1jj,1,k,l)) - aimag(f(jj,1,k,l)))
  210 continue
  220 continue
      nrx = nxhy/nxh
c bit-reverse array elements in x
      do 250 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 250
      do 240 k = kypi, kypt
      do 230 jj = 1, ndim
      t1 = f(jj,j1,k,l)
      f(jj,j1,k,l) = f(jj,j,k,l)
      f(jj,j,k,l) = t1
  230 continue
  240 continue
  250 continue
  260 continue
c then transform in x
      nrx = nxy/nxh
      do 320 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 310 l = 1, kblok
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 280 i = kypi, kypt
      do 270 jj = 1, ndim
      t1 = s*f(jj,j2,i,l)
      f(jj,j2,i,l) = f(jj,j1,i,l) - t1
      f(jj,j1,i,l) = f(jj,j1,i,l) + t1
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
  320 continue
c swap complex components
      call PSWAPC2N(f,ss,isign,nxh,kypi,kypt,nxvh,kypd,kblok,ndim)
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpp,n
     1yv,kxp,jblok,ndim,nxhyd,nxyhd)
c this subroutine performs the y part of N two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: M*(5*log2(M) + 10)/nvp
c for isign = 1,  approximate flop count: M*(5*log2(M) + 8)/nvp
c where M = (nx/2)*ny, and nvp = number of procs
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(1:N,m,n,i) = sum(g(1:N,k,j,i)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(1:N,k,j,i) = sum(g(1:N,m,n,i)*exp(sqrt(-1)*2pi*m*k/ny))
c kstrt = starting data block number
c kxpi = initial x index used
c kxpp = number of   indices used
c nyv = first dimension of g
c kxp = number of data values per block in x
c jblok = number of data blocks in x
c ndim = leading dimension of arrays f and g
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c the real data is stored in a complex array of length nx/2, ny
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c g(k,j,i) = mode jj-1,k-1, where jj = j + kxp*(i - 1)
c 1 <= jj <= nx/2 and 1 <= k <= ny, except for
c g(k,1,1) = mode nx/2,k-1, where ny/2+2 <= k <= ny, and
c imaginary part of g(1,1,1) = real part of mode nx/2,0 and
c imaginary part of g(ny/2+1,1,1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxpi, kxpp, nyv
      integer kxp, jblok, ndim, nxhyd, nxyhd
      complex g, sct
      dimension g(ndim,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, ny, nyh, ny2
      integer nxy, nxhy, ks, kxpt, j, k, nry
      integer l, i, m, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      complex s, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxpt = kxpi + kxpp - 1
      if (kstrt.gt.nxh) return
      if (isign.gt.0) go to 140
c inverse fourier transform
      nry = nxhy/ny
      do 40 l = 1, jblok
c bit-reverse array elements in y
      do 30 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 30
      do 20 j = kxpi, kxpt
      do 10 jj = 1, ndim
      t1 = g(jj,k1,j,l)
      g(jj,k1,j,l) = g(jj,k,j,l)
      g(jj,k,j,l) = t1
   10 continue
   20 continue
   30 continue
   40 continue
c then transform in y
      nry = nxy/ny
      do 100 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 90 l = 1, jblok
      do 80 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 70 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 60 i = kxpi, kxpt
      do 50 jj = 1, ndim
      t1 = s*g(jj,j2,i,l)
      g(jj,j2,i,l) = g(jj,j1,i,l) - t1
      g(jj,j1,i,l) = g(jj,j1,i,l) + t1
   50 continue
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble modes kx = 0, nx/2
      do 130 l = 1, jblok
      if ((l+ks).gt.0) go to 130
      do 120 k = 2, nyh
      if (kxpi.eq.1) then
         do 110 jj = 1, ndim
         s = g(jj,ny2-k,1,l)
         g(jj,ny2-k,1,l) = .5*cmplx(aimag(g(jj,k,1,l) + s),real(g(jj,k,1
     1,l)- s))
         g(jj,k,1,l) = .5*cmplx(real(g(jj,k,1,l) + s),aimag(g(jj,k,1,l) 
     1- s))
  110    continue
      endif
  120 continue
  130 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  140 nry = nxhy/ny
      do 210 l = 1, jblok
      if ((l+ks).gt.0) go to 170
      do 160 k = 2, nyh
      if (kxpi.eq.1) then
         do 150 jj = 1, ndim
         s = cmplx(aimag(g(jj,ny2-k,1,l)),real(g(jj,ny2-k,1,l)))
         g(jj,ny2-k,1,l) = conjg(g(jj,k,1,l) - s)
         g(jj,k,1,l) = g(jj,k,1,l) + s
  150    continue
      endif
  160 continue
c bit-reverse array elements in y
  170 do 200 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 200
      do 190 j = kxpi, kxpt
      do 180 jj = 1, ndim
      t1 = g(jj,k1,j,l)
      g(jj,k1,j,l) = g(jj,k,j,l)
      g(jj,k,j,l) = t1
  180 continue
  190 continue
  200 continue
  210 continue
c first transform in y
      nry = nxy/ny
      do 270 m = 1, indy
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 260 l = 1, jblok
      do 250 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 240 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 230 i = kxpi, kxpt
      do 220 jj = 1, ndim
      t1 = s*g(jj,j2,i,l)
      g(jj,j2,i,l) = g(jj,j1,i,l) - t1
      g(jj,j1,i,l) = g(jj,j1,i,l) + t1
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSWAPC2N(f,s,isign,nxh,kypi,kypt,nxvh,kypd,kblok,ndim)
c this subroutine swaps components for multiple ffts
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nxh = complex dimension in x direction
c kypi/kypt = initial/final y index used
c nxvh = half of the second dimension of f
c kypd = third dimension of f
c kblok = number of data blocks in y
c ndim = leading dimension of array f
      implicit none
      integer isign, nxh, kypi, kypt, nxvh, kypd, kblok, ndim
      real f, s
      dimension f(ndim,2*nxvh,kypd,kblok), s(2*ndim*nxvh)
c local data
      integer i, j, k, l, ioff
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 70 l = 1, kblok
         do 60 k = kypi, kypt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1) = f(i,2*j-1,k,l)
         s(2*i+ioff) = f(i,2*j,k,l)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k,l) = s(i+ioff)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k,l) = s(i+ioff)
   40    continue
   50    continue
   60    continue
   70    continue
      else if (isign.gt.0) then
c swap complex components
         do 140 l = 1, kblok
         do 130 k = kypi, kypt
         do 100 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 80 i = 1, ndim
         s(i+ioff) = f(i,2*j-1,k,l)
   80    continue
         ioff = ioff + ndim
         do 90 i = 1, ndim
         s(i+ioff) = f(i,2*j,k,l)
   90    continue
  100    continue
         do 120 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 110 i = 1, ndim
         f(i,2*j-1,k,l) = s(2*i+ioff-1)
         f(i,2*j,k,l) = s(2*i+ioff)
  110    continue
  120    continue
  130    continue
  140    continue
      endif
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
      subroutine WPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real sine/sine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real sine/cosine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real cosine/sine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for parallel real cosine/cosine transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp2
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
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
      subroutine PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,kx
     1pp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c sine transform and its inverse, for a subset of x,
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
c g(m,n,i) = sum(g(k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine transform is performed
c g(k,n,i) = sum(g(m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,j,l)
      at1 = g(k,j,l) + at2
      at2 = g(k,j,l) - at2
      at1 = -aimag(sctd(k1))*at1
      at2 = .5*at2
      g(k,j,l) = at1 + at2
      g(ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,j,l) = 0.0
      g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 60 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t2 = g(2*k1-1,j,l)
      t3 = g(2*k1,j,l)
      g(2*k1-1,j,l) = g(2*k-1,j,l)
      g(2*k1,j,l) = g(2*k,j,l)
      g(2*k-1,j,l) = t2
      g(2*k,j,l) = t3
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
      do 100 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t2 = real(t1)*g(2*j2-1,i,l) - aimag(t1)*g(2*j2,i,l)
      t3 = aimag(t1)*g(2*j2-1,i,l) + real(t1)*g(2*j2,i,l)
      g(2*j2-1,i,l) = g(2*j1-1,i,l) - t2
      g(2*j2,i,l) = g(2*j1,i,l) - t3
      g(2*j1-1,i,l) = g(2*j1-1,i,l) + t2
      g(2*j1,i,l) = g(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 150 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = ani*(t2 + t4)
         g(2*k,j,l) = ani*(t3 + t5)
         g(ny3-2*k,j,l) = ani*(t2 - t4)
         g(ny3-2*k+1,j,l) = ani*(t5 - t3)
  120    continue
  130    continue
         do 140 j = kxpi, kxpt
         g(nyh+1,j,l) = g(nyh+1,j,l)
         g(nyh+2,j,l) = -g(nyh+2,j,l)
         t2 = g(1,j,l) + g(2,j,l)
         g(2,j,l) = g(1,j,l) - g(2,j,l)
         g(1,j,l) = t2
         g(ny+1,j,l) = g(ny+1,j,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 170 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 160 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = t2 + t4
         g(2*k,j,l) = t3 + t5
         g(ny3-2*k,j,l) = t2 - t4
         g(ny3-2*k+1,j,l) = t5 - t3
  160    continue
  170    continue
         do 180 j = kxpi, kxpt
         g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
         g(nyh+2,j,l) = -2.0*g(nyh+2,j,l)
         t2 = 2.0*(g(1,j,l) + g(2,j,l))
         g(2,j,l) = 2.0*(g(1,j,l) - g(2,j,l))
         g(1,j,l) = t2
         g(ny+1,j,l) = 2.0*g(ny+1,j,l)
  180    continue
  190    continue
      endif
c perform recursion for sine transform
      kxpt = kxps
      do 220 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 210 j = kxpi, kxpt
      sum1 = .5*g(1,j,l)
      g(1,j,l) = 0.0
      g(2,j,l) = sum1
      do 200 k = 2, nyh
      sum1 = sum1 + g(2*k-1,j,l)
      g(2*k-1,j,l) = -g(2*k,j,l)
      g(2*k,j,l) = sum1
  200 continue
      g(ny+1,j,l) = 0.0
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,kx
     1pp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c cosine transform and its inverse, for a subset of x,
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
c g(m,n,i) = (.5*g(1,n,i) + ((-1)**m)*g(ny+1,n,i)
c            + sum(g(k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward cosine transform is performed
c g(k,n,i) = 2*(.5*g(1,n,i) + ((-1)**m)*g(ny+1,n,i) + sum(g(m,n,i)*
c            cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt
      real at1, at2, t2, t3, t4, t5, t6, ani, sum1
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,j,l) - g(ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = g(ny+2-k,j,l)
      at1 = g(k,j,l) + at2
      at2 = g(k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = -aimag(sctd(k1))*at2
      at1 = .5*at1
      g(k,j,l) = at1 - at2
      g(ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,j,l) = .5*(g(1,j,l) + g(ny+1,j,l))
      g(ny+1,j,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 60 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 50 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = kxpi, kxpt
      t2 = g(2*k1-1,j,l)
      t3 = g(2*k1,j,l)
      g(2*k1-1,j,l) = g(2*k-1,j,l)
      g(2*k1,j,l) = g(2*k,j,l)
      g(2*k-1,j,l) = t2
      g(2*k,j,l) = t3
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
      do 100 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = kxpi, kxpt
      t2 = real(t1)*g(2*j2-1,i,l) - aimag(t1)*g(2*j2,i,l)
      t3 = aimag(t1)*g(2*j2-1,i,l) + real(t1)*g(2*j2,i,l)
      g(2*j2-1,i,l) = g(2*j1-1,i,l) - t2
      g(2*j2,i,l) = g(2*j1,i,l) - t3
      g(2*j1-1,i,l) = g(2*j1-1,i,l) + t2
      g(2*j1,i,l) = g(2*j1,i,l) + t3
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 150 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = ani*(t2 + t4)
         g(2*k,j,l) = ani*(t3 + t5)
         g(ny3-2*k,j,l) = ani*(t2 - t4)
         g(ny3-2*k+1,j,l) = ani*(t5 - t3)
  120    continue
  130    continue
         do 140 j = kxpi, kxpt
         g(nyh+1,j,l) = g(nyh+1,j,l)
         g(nyh+2,j,l) = -g(nyh+2,j,l)
         t2 = g(1,j,l) + g(2,j,l)
         g(2,j,l) = g(1,j,l) - g(2,j,l)
         g(1,j,l) = t2
         g(ny+1,j,l) = g(ny+1,j,l)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 170 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 160 j = kxpi, kxpt
         t4 = g(ny3-2*k,j,l)
         t5 = -g(ny3-2*k+1,j,l)
         t2 = g(2*k-1,j,l) + t4
         t3 = g(2*k,j,l) + t5
         t6 = g(2*k-1,j,l) - t4
         t5 = g(2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(2*k-1,j,l) = t2 + t4
         g(2*k,j,l) = t3 + t5
         g(ny3-2*k,j,l) = t2 - t4
         g(ny3-2*k+1,j,l) = t5 - t3
  160    continue
  170    continue
         do 180 j = kxpi, kxpt
         g(nyh+1,j,l) = 2.0*g(nyh+1,j,l)
         g(nyh+2,j,l) = -2.0*g(nyh+2,j,l)
         t2 = 2.0*(g(1,j,l) + g(2,j,l))
         g(2,j,l) = 2.0*(g(1,j,l) - g(2,j,l))
         g(1,j,l) = t2
         g(ny+1,j,l) = 2.0*g(ny+1,j,l)
  180    continue
  190    continue
      endif
c perform recursion for cosine transform
      kxpt = kxps
      do 220 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 210 j = kxpi, kxpt
      sum1 = g(ny+1,j,l)
      g(ny+1,j,l) = g(2,j,l)
      g(2,j,l) = sum1
      do 200 k = 2, nyh
      sum1 = sum1 - g(2*k,j,l)
      g(2*k,j,l) = sum1
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real sine/sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine-cosine transform
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine-cosine transform
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 2 parallel real sine/sine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine-sine transform
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine-sine transform
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
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
      subroutine PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,k
     1xpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = sum(g(1,k,n,i)*sin(pi*m*k/ny))
c g(2,m,n,i) = (.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = sum(g(1,m,n,i)*sin(pi*m*k/ny))
c g(2,k,n,i) = 2*(.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,m,n,i)*cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(2,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(2,1,j,l) - g(2,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(1,k,j,l) = at1 + at2
      g(1,ny+2-k,j,l) = at1 - at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(2,k,j,l) = at1 - at2
      g(2,ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,1,j,l) = 0.0
      g(1,nyh+1,j,l) = 2.0*g(1,nyh+1,j,l)
      g(2,1,j,l) = .5*(g(2,1,j,l) + g(2,ny+1,j,l))
      g(2,ny+1,j,l) = sum1
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 2
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 2
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 2
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = .5*g(1,1,j,l)
      g(1,1,j,l) = 0.0
      g(1,2,j,l) = sum1
      sum2 = g(2,ny+1,j,l)
      g(2,ny+1,j,l) = g(2,2,j,l)
      g(2,2,j,l) = sum2
      do 260 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,j,l)
      g(1,2*k-1,j,l) = -g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 - g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
  260 continue
      g(1,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,k
     1xpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = (.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,k,n,i)*cos(pi*m*k/ny))
c g(2,m,n,i) = sum(g(2,k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = 2*(.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,m,n,i)*cos(pi*m*k/ny))
c g(2,k,n,i) = sum(g(2,m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(2,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,1,j,l) - g(1,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(1,k,j,l) = at1 - at2
      g(1,ny+2-k,j,l) = at1 + at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(2,k,j,l) = at1 + at2
      g(2,ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,1,j,l) = .5*(g(1,1,j,l) + g(1,ny+1,j,l))
      g(1,ny+1,j,l) = sum1
      g(2,1,j,l) = 0.0
      g(2,nyh+1,j,l) = 2.0*g(2,nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 2
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 2
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 2
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 2
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 2
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = g(1,ny+1,j,l)
      g(1,ny+1,j,l) = g(1,2,j,l)
      g(1,2,j,l) = sum1
      sum2 = .5*g(2,1,j,l)
      g(2,1,j,l) = 0.0
      g(2,2,j,l) = sum2
      do 260 k = 2, nyh
      sum1 = sum1 - g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 + g(2,2*k-1,j,l)
      g(2,2*k-1,j,l) = -g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
  260 continue
      g(2,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine/sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y sine-cosine transform
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y sine-cosine transform
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c perform x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
c wrapper function for 3 parallel real sine/sine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y cosine-sine transform
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform y cosine-sine transform
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1p2,nyv,kxp2d,jblok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
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
      subroutine PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,
     1kxpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = sum(g(1,k,n,i)*sin(pi*m*k/ny))
c g(2,m,n,i) = (.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,k,n,i)*cos(pi*m*k/ny))
c g(3,m,n,i) = sum(g(3,k,n,i)*sin(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = sum(g(1,m,n,i)*sin(pi*m*k/ny))
c g(2,k,n,i) = 2*(.5*g(2,1,n,i) + ((-1)**m)*g(2,ny+1,n,i)
c              + sum(g(2,m,n,i)*cos(pi*m*k/ny))
c g(3,k,n,i) = sum(g(3,m,n,i)*sin(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(3,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(2,1,j,l) - g(2,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(1,k,j,l) = at1 + at2
      g(1,ny+2-k,j,l) = at1 - at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(2,k,j,l) = at1 - at2
      g(2,ny+2-k,j,l) = at1 + at2
      at2 = g(3,ny+2-k,j,l)
      at1 = g(3,k,j,l) + at2
      at2 = g(3,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(3,k,j,l) = at1 + at2
      g(3,ny+2-k,j,l) = at1 - at2
   10 continue
      g(1,1,j,l) = 0.0
      g(1,nyh+1,j,l) = 2.0*g(1,nyh+1,j,l)
      g(2,1,j,l) = .5*(g(2,1,j,l) + g(2,ny+1,j,l))
      g(2,ny+1,j,l) = sum1
      g(3,1,j,l) = 0.0
      g(3,nyh+1,j,l) = 2.0*g(3,nyh+1,j,l)
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 3
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 3
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 3
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = .5*g(1,1,j,l)
      g(1,1,j,l) = 0.0
      g(1,2,j,l) = sum1
      sum2 = g(2,ny+1,j,l)
      g(2,ny+1,j,l) = g(2,2,j,l)
      g(2,2,j,l) = sum2
      sum3 = .5*g(3,1,j,l)
      g(3,1,j,l) = 0.0
      g(3,2,j,l) = sum3
      do 260 k = 2, nyh
      sum1 = sum1 + g(1,2*k-1,j,l)
      g(1,2*k-1,j,l) = -g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 - g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
      sum3 = sum3 + g(3,2*k-1,j,l)
      g(3,2*k-1,j,l) = -g(3,2*k,j,l)
      g(3,2*k,j,l) = sum3
  260 continue
      g(1,ny+1,j,l) = 0.0
      g(3,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp,kxpi,
     1kxpp,nyv,kxpd,jblok,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transform are performed
c g(1,m,n,i) = (.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,k,n,i)*cos(pi*m*k/ny))
c g(2,m,n,i) = sum(g(2,k,n,i)*sin(pi*m*k/ny))
c g(3,m,n,i) = (.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,k,n,i)*cos(pi*m*k/ny))
c if isign = 1, a forward sine-cosine transforms are performed
c g(1,k,n,i) = 2*(.5*g(1,1,n,i) + ((-1)**m)*g(1,ny+1,n,i)
c              + sum(g(1,m,n,i)*cos(pi*m*k/ny))
c g(2,k,n,i) = sum(g(2,m,n,i)*sin(pi*m*k/ny))
c g(3,k,n,i) = 2*(.5*g(3,1,n,i) + ((-1)**m)*g(3,ny+1,n,i)
c              + sum(g(3,m,n,i)*cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c kstrt = starting data block number
c kxp = number of data values per block in x
c kxpi = initial x index used
c kxpp = number of x indices used
c nyv = first dimension of g >= ny + 1
c kxpd = second dimension of g >= kxp + 1
c jblok = number of data blocks in x
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, kstrt, kxp, kxpi, kxpp
      integer nyv, kxpd, jblok, nxhyd, nxyd
      real g
      complex sctd
      dimension g(3,nyv,kxpd,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indy1, indx1y, nx, ny, nyh, nyhh, ny3, nxy, nxhy
      integer i, j, k, l, m, ks, km, kmr, nry, j1, j2, ns, ns2, k1, k2
      integer kxps, kxpt, jj
      real at1, at2, at3, t2, t3, t4, t5, t6, ani, sum1, sum2, sum3
      complex t1
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny3 = ny + 3
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      ks = kstrt - 2
      kxps = kxpi + kxpp - 1
      if (kstrt.gt.nx) return
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      kxpt = kxps
      do 30 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 20 j = kxpi, kxpt
      sum1 = .5*(g(1,1,j,l) - g(1,ny+1,j,l))
      sum2 = .5*(g(3,1,j,l) - g(3,ny+1,j,l))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = g(1,ny+2-k,j,l)
      at1 = g(1,k,j,l) + at2
      at2 = g(1,k,j,l) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(1,k,j,l) = at1 - at2
      g(1,ny+2-k,j,l) = at1 + at2
      at2 = g(2,ny+2-k,j,l)
      at1 = g(2,k,j,l) + at2
      at2 = g(2,k,j,l) - at2
      at1 = at3*at1
      at2 = .5*at2
      g(2,k,j,l) = at1 + at2
      g(2,ny+2-k,j,l) = at1 - at2
      at2 = g(3,ny+2-k,j,l)
      at1 = g(3,k,j,l) + at2
      at2 = g(3,k,j,l) - at2
      sum2 = sum2 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      g(3,k,j,l) = at1 - at2
      g(3,ny+2-k,j,l) = at1 + at2
   10 continue
      g(1,1,j,l) = .5*(g(1,1,j,l) + g(1,ny+1,j,l))
      g(1,ny+1,j,l) = sum1
      g(2,1,j,l) = 0.0
      g(2,nyh+1,j,l) = 2.0*g(2,nyh+1,j,l)
      g(3,1,j,l) = .5*(g(3,1,j,l) + g(3,ny+1,j,l))
      g(3,ny+1,j,l) = sum2
   20 continue
   30 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      kxpt = kxps
      do 70 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 60 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 60
      do 50 j = kxpi, kxpt
      do 40 jj = 1, 3
      t2 = g(jj,2*k1-1,j,l)
      t3 = g(jj,2*k1,j,l)
      g(jj,2*k1-1,j,l) = g(jj,2*k-1,j,l)
      g(jj,2*k1,j,l) = g(jj,2*k,j,l)
      g(jj,2*k-1,j,l) = t2
      g(jj,2*k,j,l) = t3
   40 continue
   50 continue
   60 continue
   70 continue
c first transform in y
      nry = nxy/nyh
      do 130 m = 1, indy1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      kxpt = kxps
      do 120 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 90 i = kxpi, kxpt
      do 80 jj = 1, 3
      t2 = real(t1)*g(jj,2*j2-1,i,l) - aimag(t1)*g(jj,2*j2,i,l)
      t3 = aimag(t1)*g(jj,2*j2-1,i,l) + real(t1)*g(jj,2*j2,i,l)
      g(jj,2*j2-1,i,l) = g(jj,2*j1-1,i,l) - t2
      g(jj,2*j2,i,l) = g(jj,2*j1,i,l) - t3
      g(jj,2*j1-1,i,l) = g(jj,2*j1-1,i,l) + t2
      g(jj,2*j1,i,l) = g(jj,2*j1,i,l) + t3
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         kxpt = kxps
         do 190 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 160 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 150 j = kxpi, kxpt
         do 140 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = ani*(t2 + t4)
         g(jj,2*k,j,l) = ani*(t3 + t5)
         g(jj,ny3-2*k,j,l) = ani*(t2 - t4)
         g(jj,ny3-2*k+1,j,l) = ani*(t5 - t3)
  140    continue
  150    continue
  160    continue
         do 180 j = kxpi, kxpt
         do 170 jj = 1, 3
         g(jj,nyh+1,j,l) = g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -g(jj,nyh+2,j,l)
         t2 = g(jj,1,j,l) + g(jj,2,j,l)
         g(jj,2,j,l) = g(jj,1,j,l) - g(jj,2,j,l)
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = g(jj,ny+1,j,l)
  170    continue
  180    continue
  190    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         kxpt = kxps
         do 250 l = 1, jblok
         if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
         do 220 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 210 j = kxpi, kxpt
         do 200 jj = 1, 3
         t4 = g(jj,ny3-2*k,j,l)
         t5 = -g(jj,ny3-2*k+1,j,l)
         t2 = g(jj,2*k-1,j,l) + t4
         t3 = g(jj,2*k,j,l) + t5
         t6 = g(jj,2*k-1,j,l) - t4
         t5 = g(jj,2*k,j,l) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         g(jj,2*k-1,j,l) = t2 + t4
         g(jj,2*k,j,l) = t3 + t5
         g(jj,ny3-2*k,j,l) = t2 - t4
         g(jj,ny3-2*k+1,j,l) = t5 - t3
  200    continue
  210    continue
  220    continue
         do 240 j = kxpi, kxpt
         do 230 jj = 1, 3
         g(jj,nyh+1,j,l) = 2.0*g(jj,nyh+1,j,l)
         g(jj,nyh+2,j,l) = -2.0*g(jj,nyh+2,j,l)
         t2 = 2.0*(g(jj,1,j,l) + g(jj,2,j,l))
         g(jj,2,j,l) = 2.0*(g(jj,1,j,l) - g(jj,2,j,l))
         g(jj,1,j,l) = t2
         g(jj,ny+1,j,l) = 2.0*g(jj,ny+1,j,l)
  230    continue
  240    continue
  250    continue
      endif
c perform recursion for sine-cosine transform
      kxpt = kxps
      do 280 l = 1, jblok
      if ((kxps+kxp*(l+ks)).eq.nx) kxpt = kxps + 1
      do 270 j = kxpi, kxpt
      sum1 = g(1,ny+1,j,l)
      g(1,ny+1,j,l) = g(1,2,j,l)
      g(1,2,j,l) = sum1
      sum2 = .5*g(2,1,j,l)
      g(2,1,j,l) = 0.0
      g(2,2,j,l) = sum2
      sum3 = g(3,ny+1,j,l)
      g(3,ny+1,j,l) = g(3,2,j,l)
      g(3,2,j,l) = sum3
      do 260 k = 2, nyh
      sum1 = sum1 - g(1,2*k,j,l)
      g(1,2*k,j,l) = sum1
      sum2 = sum2 + g(2,2*k-1,j,l)
      g(2,2*k-1,j,l) = -g(2,2*k,j,l)
      g(2,2*k,j,l) = sum2
      sum3 = sum3 - g(3,2*k,j,l)
      g(3,2*k,j,l) = sum3
  260 continue
      g(2,ny+1,j,l) = 0.0
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFSFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for parallel real sine/periodic transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
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
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 40 l = 1, j2blok
         do 30 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   30    continue
   40    continue
c perform x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFCFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
c wrapper function for parallel real cosine/periodic transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
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
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 40 l = 1, j2blok
         do 30 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   30    continue
   40    continue
c perform x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp,n
     1xvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
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
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
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
         do 50 j = 1, nxh1
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
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
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
         do 50 j = 1, nxh1
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
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
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
         do 50 j = 1, nxh1
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
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
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
         do 50 j = 1, nxh1
         do 40 i = 1, 3
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
      subroutine WPFDT2RINIT(sctdx,indx,nxd)
c this subroutine calculates extra table needed by a one dimensional
c fast real sine and cosine transforms for mixed boundary conditions
c and their inverses.
c input: indx, nxd
c output: sctdx
c sctdx = sine/cosine table
c indx = exponent which determines length in x direction,
c where nx=2**indx
c nxd = must be >= nx
c written by viktor k. decyk, ucla
      implicit none
      integer indx, nxd
      complex sctdx
      dimension sctdx(nxd)
c local data
      integer j, nx
      real dnx, arg, at1, at2, at3, at4
      nx = 2**indx
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/float(nx + nx)
      arg = 0.5*dnx
      at3 = cos(arg)
      at4 = sin(arg)
      do 10 j = 1, nx
      arg = dnx*float(j - 1)
      at1 = cos(arg)
      at2 = sin(arg)
      at1 = at1*at4 + at2*at3
      sctdx(j) = cmplx(at1,1.0/at1)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDSFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp,
     1indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,n
     2xyd)
c wrapper function for parallel real sine/periodic transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c perform x sine transform
         call PFDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
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
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 40 l = 1, j2blok
         do 30 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   30    continue
   40    continue
c perform x sine transform
         call PFDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDCFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp,
     1indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,n
     2xyd)
c wrapper function for parallel real cosine/periodic transform
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c perform x sine transform
         call PFDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
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
c perform y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     12,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 40 l = 1, j2blok
         do 30 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   30    continue
   40    continue
c perform x sine transform
         call PFDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,
     1kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is similar to the one described in Numerical Recipies in
c Fortran, Second Ed.,by W. H. Press, B. P. Flannery, S. A. Teukolsky,
c and W. T. Vetterling, [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform DST-III is performed
c f(n,k,i) = (1/nx*ny)*(.5*f(nx+1,k,i)*(-1)**n + sum(f(j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward sine transform DST-II is performed
c f(j+1,k,i) = 2.0*sum(f(n,k,i)*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1
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
      if (isign.gt.0) go to 190
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      f(1,k,l) = 2.0*f(2,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(2*j1,k,l) - f(2*j1-2,k,l)
      at3 = f(2*j1-1,k,l)*at2 + at4*at1
      f(2*j1,k,l) = at4*at2 - f(2*j1-1,k,l)*at1
      f(2*j1-1,k,l) = at3
   10 continue
      f(2,k,l) = f(nx+1,k,l)
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = kypi, kypt
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
   40 continue
   50 continue
      do 60 k = kypi, kypt
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
      f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
      t2 = f(1,k,l) + f(2,k,l)
      f(2,k,l) = f(1,k,l) - f(2,k,l)
      f(1,k,l) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 90
      do 80 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 140 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 110 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) + aimag(t1)*f(2*j2,i,l)
      t3 = real(t1)*f(2*j2,i,l) - aimag(t1)*f(2*j2-1,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for sine transform and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = kypi, kypt
      do 160 j = 1, nxh
      at2 = f(nx+1-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(j,k,l) = ani*(at1 + at2)
      f(nx+1-j,k,l) = ani*(at1 - at2)
  160 continue
      f(nx+1,k,l) = 0.0
  170 continue
  180 continue
      return
c forward fourier transform
c create auxiliary array in x
  190 kmr = nxy/nx
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(nx+1-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(j,k,l) = at1 + at2
      f(nx+1-j,k,l) = at1 - at2
  200 continue
  210 continue
  220 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 250 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 240 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 240
      do 230 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
  230 continue
  240 continue
  250 continue
c first transform in x
      nrx = nxy/nxh
      do 300 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 290 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 280 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 270 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 260 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 340 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 320 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 310 k = kypi, kypt
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
  310 continue
  320 continue
      do 330 k = kypi, kypt
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
      f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
      t2 = 2.0*(f(1,k,l) + f(2,k,l))
      f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
      f(1,k,l) = t2
      f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  330 continue
  340 continue
c perform recursion for sine transform
      kmr = nxy/nx
      kypt = kyps
      do 370 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 360 k = kypi, kypt
      f(nx+1,k,l) = f(2,k,l)
      f(2,k,l) = 0.5*f(1,k,l)
      f(1,k,l) = 0.0
      sum1 = f(2,k,l)
      do 350 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2*j-1,k,l)*at2 - f(2*j,k,l)*at1
      at2 = f(2*j-1,k,l)*at1 + f(2*j,k,l)*at2
      f(2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(2*j,k,l) = sum1
  350 continue
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,
     1kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform DCT-III is performed
c f(n,k,i) = (1/nx*ny)*(.5*f(1,k,i) + sum(f(j,k,i)*cos(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward cosine transform DCT-II is performed
c f(j,k,i) = 2.0*sum(f(n,k,i)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1
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
      if (isign.gt.0) go to 190
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = f(nx,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(2*j1,k,l) - f(2*j1-2,k,l)
      at3 = f(2*j1-1,k,l)*at1 + at4*at2
      f(2*j1,k,l) = f(2*j1-1,k,l)*at2 - at4*at1
      f(2*j1-1,k,l) = at3
   10 continue
      f(2,k,l) = -2.0*sum1
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 70 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = kypi, kypt
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
   40 continue
   50 continue
      do 60 k = kypi, kypt
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
      f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
      t2 = f(1,k,l) + f(2,k,l)
      f(2,k,l) = f(1,k,l) - f(2,k,l)
      f(1,k,l) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 100 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 90
      do 80 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 140 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 110 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) + aimag(t1)*f(2*j2,i,l)
      t3 = real(t1)*f(2*j2,i,l) - aimag(t1)*f(2*j2-1,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for cosine transform and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = kypi, kypt
      do 160 j = 1, nxh
      at2 = f(nx+1-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(j,k,l) = ani*(at1 - at2)
      f(nx+1-j,k,l) = ani*(at1 + at2)
  160 continue
      f(nx+1,k,l) = 0.0
  170 continue
  180 continue
      return
c forward fourier transform
c create auxiliary array in x
  190 kmr = nxy/nx
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(nx+1-j,k,l)
      at1 = f(j,k,l) + at2
      at2 = f(j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(j,k,l) = at1 - at2
      f(nx+1-j,k,l) = at1 + at2
  200 continue
  210 continue
  220 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 250 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 240 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 240
      do 230 k = kypi, kypt
      t2 = f(2*j1-1,k,l)
      t3 = f(2*j1,k,l)
      f(2*j1-1,k,l) = f(2*j-1,k,l)
      f(2*j1,k,l) = f(2*j,k,l)
      f(2*j-1,k,l) = t2
      f(2*j,k,l) = t3
  230 continue
  240 continue
  250 continue
c first transform in x
      nrx = nxy/nxh
      do 300 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 290 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 280 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 270 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 260 i = kypi, kypt
      t2 = real(t1)*f(2*j2-1,i,l) - aimag(t1)*f(2*j2,i,l)
      t3 = aimag(t1)*f(2*j2-1,i,l) + real(t1)*f(2*j2,i,l)
      f(2*j2-1,i,l) = f(2*j1-1,i,l) - t2
      f(2*j2,i,l) = f(2*j1,i,l) - t3
      f(2*j1-1,i,l) = f(2*j1-1,i,l) + t2
      f(2*j1,i,l) = f(2*j1,i,l) + t3
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 340 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 320 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 310 k = kypi, kypt
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
  310 continue
  320 continue
      do 330 k = kypi, kypt
      f(nxh+1,k,l) = 2.0*f(nxh+1,k,l)
      f(nxh+2,k,l) = -2.0*f(nxh+2,k,l)
      t2 = 2.0*(f(1,k,l) + f(2,k,l))
      f(2,k,l) = 2.0*(f(1,k,l) - f(2,k,l))
      f(1,k,l) = t2
      f(nx+1,k,l) = 2.0*f(nx+1,k,l)
  330 continue
  340 continue
c perform recursion for cosine transform
      kmr = nxy/nx
      kypt = kyps
      do 370 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 360 k = kypi, kypt
      sum1 = -0.5*f(2,k,l)
      do 350 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2*j1-1,k,l)*at1 + f(2*j1,k,l)*at2
      at2 = -f(2*j1-1,k,l)*at2 + f(2*j1,k,l)*at1
      f(2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2*j1,k,l) = at1
  350 continue
      f(2,k,l) = sum1
      f(nx+1,k,l) = 0.0
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd)
c wrapper function for 2 parallel real cosine-sine/periodic transforms
c for the electric field with mixed dirichlet-neumann or magnetic field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x cosine-sine transform
         call PFDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
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
         do 50 j = 1, nxh1
         do 40 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x cosine-sine transform
         call PFDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd)
c wrapper function for 2 parallel real sine-cosine/periodic transforms
c for the magnetic field with mixed dirichlet-neumann or electric field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x sine-cosine transform
         call PFDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
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
         do 50 j = 1, nxh1
         do 40 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x sine-cosine transform
         call PFDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp
     1,kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + sum(f(1,j,k,i)*
c cos(pi*(n+1/2)*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,nx+1,k,i)*(-1)**n + sum(f(2,j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j,k,i) = 2.0*sum(f(1,n,k,i)*cos(pi*n*(j+1/2)/nx))
c f(2,j+1,k,i) = 2.0*sum(f(2,n,k),i*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1
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
      if (isign.gt.0) go to 230
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = f(1,nx,k,l)
      f(2,1,k,l) = 2.0*f(2,2,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k,l) - f(1,2*j1-2,k,l)
      at3 = f(1,2*j1-1,k,l)*at1 + at4*at2
      f(1,2*j1,k,l) = f(1,2*j1-1,k,l)*at2 - at4*at1
      f(1,2*j1-1,k,l) = at3
      at4 = f(2,2*j1,k,l) - f(2,2*j1-2,k,l)
      at3 = f(2,2*j1-1,k,l)*at2 + at4*at1
      f(2,2*j1,k,l) = at4*at2 - f(2,2*j1-1,k,l)*at1
      f(2,2*j1-1,k,l) = at3
   10 continue
      f(1,2,k,l) = -2.0*sum1
      f(2,2,k,l) = f(2,nx+1,k,l)
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 90 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
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
   40 continue
   50 continue
   60 continue
      do 80 k = kypi, kypt
      do 70 jj = 1, 2
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = f(jj,1,k,l) + f(jj,2,k,l)
      f(jj,2,k,l) = f(jj,1,k,l) - f(jj,2,k,l)
      f(jj,1,k,l) = t2
   70 continue
   80 continue
   90 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 130 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 120 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 120
      do 110 k = kypi, kypt
      do 100 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  100 continue
  110 continue
  120 continue
  130 continue
c first transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 150 i = kypi, kypt
      do 140 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) + aimag(t1)*f(jj,2*j2,i,l)
      t3 = real(t1)*f(jj,2*j2,i,l) - aimag(t1)*f(jj,2*j2-1,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(1,j,k,l) = ani*(at1 - at2)
      f(1,nx+1-j,k,l) = ani*(at1 + at2)
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(2,j,k,l) = ani*(at1 + at2)
      f(2,nx+1-j,k,l) = ani*(at1 - at2)
  200 continue
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
c forward fourier transform
c create auxiliary array in x
  230 kmr = nxy/nx
      kypt = kyps
      do 260 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 250 k = kypi, kypt
      do 240 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+1-j,k,l) = at1 + at2
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+1-j,k,l) = at1 - at2
  240 continue
  250 continue
  260 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 300 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 290 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 290
      do 280 k = kypi, kypt
      do 270 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  270 continue
  280 continue
  290 continue
  300 continue
c first transform in x
      nrx = nxy/nxh
      do 360 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 350 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 320 i = kypi, kypt
      do 310 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 420 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 390 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 380 k = kypi, kypt
      do 370 jj = 1, 2
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
  370 continue
  380 continue
  390 continue
      do 410 k = kypi, kypt
      do 400 jj = 1, 2
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
      f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
      f(jj,1,k,l) = t2
      f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  400 continue
  410 continue
  420 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      kypt = kyps
      do 460 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 450 k = kypi, kypt
      sum1 = -0.5*f(1,2,k,l)
      do 430 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(1,2*j1-1,k,l)*at1 + f(1,2*j1,k,l)*at2
      at2 = -f(1,2*j1-1,k,l)*at2 + f(1,2*j1,k,l)*at1
      f(1,2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(1,2*j1,k,l) = at1
  430 continue
      f(1,2,k,l) = sum1
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = 0.5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      sum1 = f(2,2,k,l)
      do 440 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2,2*j-1,k,l)*at2 - f(2,2*j,k,l)*at1
      at2 = f(2,2*j-1,k,l)*at1 + f(2,2*j,k,l)*at2
      f(2,2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(2,2*j,k,l) = sum1
  440 continue
  450 continue
  460 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp
     1,kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,nx+1,k,i)*(-1)**n + sum(f(1,j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,1,k,i) + sum(f(2,j,k,i)*
c cos(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j+1,k,i) = 2.0*sum(f(1,n,k,i)*sin(pi*n*(j+1/2)/nx)))
c f(2,j,k,i) = 2.0*sum(f(2,n,k,i)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(2,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1
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
      if (isign.gt.0) go to 230
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      f(1,1,k,l) = 2.0*f(1,2,k,l)
      sum1 = f(2,nx,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k,l) - f(1,2*j1-2,k,l)
      at3 = f(1,2*j1-1,k,l)*at2 + at4*at1
      f(1,2*j1,k,l) = at4*at2 - f(1,2*j1-1,k,l)*at1
      f(1,2*j1-1,k,l) = at3
      at4 = f(2,2*j1,k,l) - f(2,2*j1-2,k,l)
      at3 = f(2,2*j1-1,k,l)*at1 + at4*at2
      f(2,2*j1,k,l) = f(2,2*j1-1,k,l)*at2 - at4*at1
      f(2,2*j1-1,k,l) = at3
   10 continue
      f(1,2,k,l) = f(1,nx+1,k,l)
      f(2,2,k,l) = -2.0*sum1
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 90 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 50 k = kypi, kypt
      do 40 jj = 1, 2
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
   40 continue
   50 continue
   60 continue
      do 80 k = kypi, kypt
      do 70 jj = 1, 2
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = f(jj,1,k,l) + f(jj,2,k,l)
      f(jj,2,k,l) = f(jj,1,k,l) - f(jj,2,k,l)
      f(jj,1,k,l) = t2
   70 continue
   80 continue
   90 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 130 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 120 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 120
      do 110 k = kypi, kypt
      do 100 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  100 continue
  110 continue
  120 continue
  130 continue
c first transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 150 i = kypi, kypt
      do 140 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) + aimag(t1)*f(jj,2*j2,i,l)
      t3 = real(t1)*f(jj,2*j2,i,l) - aimag(t1)*f(jj,2*j2-1,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(1,j,k,l) = ani*(at1 + at2)
      f(1,nx+1-j,k,l) = ani*(at1 - at2)
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(2,j,k,l) = ani*(at1 - at2)
      f(2,nx+1-j,k,l) = ani*(at1 + at2)
  200 continue
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
c forward fourier transform
c create auxiliary array in x
  230 kmr = nxy/nx
      kypt = kyps
      do 260 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 250 k = kypi, kypt
      do 240 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+1-j,k,l) = at1 - at2
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+1-j,k,l) = at1 + at2
  240 continue
  250 continue
  260 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 300 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 290 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 290
      do 280 k = kypi, kypt
      do 270 jj = 1, 2
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  270 continue
  280 continue
  290 continue
  300 continue
c first transform in x
      nrx = nxy/nxh
      do 360 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 350 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 320 i = kypi, kypt
      do 310 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 420 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 390 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 380 k = kypi, kypt
      do 370 jj = 1, 2
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
  370 continue
  380 continue
  390 continue
      do 410 k = kypi, kypt
      do 400 jj = 1, 2
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
      f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
      f(jj,1,k,l) = t2
      f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  400 continue
  410 continue
  420 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      kypt = kyps
      do 460 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 450 k = kypi, kypt
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = 0.5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      sum1 = f(1,2,k,l)
      do 430 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(1,2*j-1,k,l)*at2 - f(1,2*j,k,l)*at1
      at2 = f(1,2*j-1,k,l)*at1 + f(1,2*j,k,l)*at2
      f(1,2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(1,2*j,k,l) = sum1
  430 continue
      sum1 = -0.5*f(2,2,k,l)
      do 440 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2,2*j1-1,k,l)*at1 + f(2,2*j1,k,l)*at2
      at2 = -f(2,2*j1-1,k,l)*at2 + f(2,2*j1,k,l)*at1
      f(2,2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2,2*j1,k,l) = at1
  440 continue
      f(2,2,k,l) = sum1
      f(2,nx+1,k,l) = 0.0
  450 continue
  460 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd)
c wrapper function for 3 parallel real cosine-sine/periodic transforms
c for the electric field with mixed dirichlet-neumann or magnetic field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x cosine-sine transform
         call PFDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
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
         do 50 j = 1, nxh1
         do 40 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x cosine-sine transform
         call PFDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFDSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd)
c wrapper function for 3 parallel real cosine-sine/periodic transforms
c for the magnetic field with mixed dirichlet-neumann or electric field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer i, j, l, nx, ny, nxh1, kxpi, kypi
      real tf
      double precision dtime
      data kxpi, kypi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c perform x cosine-sine transform
         call PFDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
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
         do 50 j = 1, nxh1
         do 40 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   40    continue
   50    continue
   60    continue
c perform x cosine-sine transform
         call PFDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kyp,nxvh,kypd,kblok,nxhyd,nxyd)
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,ky
     1p,kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,1,k,i) + sum(f(1,j,k,i)*
c cos(pi*(n+1/2)*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,nx+1,k,i)*(-1)**n + sum(f(2,j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c f(3,n,k,i) = (1/nx*ny)*(.5*f(3,nx+1,k,i)*(-1)**n + sum(f(3,j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j,k,i) = 2.0*sum(f(1,n,k,i)*cos(pi*n*(j+1/2)/nx))
c f(2,j+1,k,i) = 2.0*sum(f(2,n,k,i)*sin(pi*n*(j+1/2)/nx)))
c f(3,j+1,k,i) = 2.0*sum(f(3,n,k),i*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1
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
      if (isign.gt.0) go to 230
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      sum1 = f(1,nx,k,l)
      f(2,1,k,l) = 2.0*f(2,2,k,l)
      f(3,1,k,l) = 2.0*f(3,2,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k,l) - f(1,2*j1-2,k,l)
      at3 = f(1,2*j1-1,k,l)*at1 + at4*at2
      f(1,2*j1,k,l) = f(1,2*j1-1,k,l)*at2 - at4*at1
      f(1,2*j1-1,k,l) = at3
      at4 = f(2,2*j1,k,l) - f(2,2*j1-2,k,l)
      at3 = f(2,2*j1-1,k,l)*at2 + at4*at1
      f(2,2*j1,k,l) = at4*at2 - f(2,2*j1-1,k,l)*at1
      f(2,2*j1-1,k,l) = at3
      at4 = f(3,2*j1,k,l) - f(3,2*j1-2,k,l)
      at3 = f(3,2*j1-1,k,l)*at2 + at4*at1
      f(3,2*j1,k,l) = at4*at2 - f(3,2*j1-1,k,l)*at1
      f(3,2*j1-1,k,l) = at3
   10 continue
      f(1,2,k,l) = -2.0*sum1
      f(2,2,k,l) = f(2,nx+1,k,l)
      f(3,2,k,l) = f(3,nx+1,k,l)
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 90 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
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
   40 continue
   50 continue
   60 continue
      do 80 k = kypi, kypt
      do 70 jj = 1, 3
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = f(jj,1,k,l) + f(jj,2,k,l)
      f(jj,2,k,l) = f(jj,1,k,l) - f(jj,2,k,l)
      f(jj,1,k,l) = t2
   70 continue
   80 continue
   90 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 130 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 120 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 120
      do 110 k = kypi, kypt
      do 100 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  100 continue
  110 continue
  120 continue
  130 continue
c first transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 150 i = kypi, kypt
      do 140 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) + aimag(t1)*f(jj,2*j2,i,l)
      t3 = real(t1)*f(jj,2*j2,i,l) - aimag(t1)*f(jj,2*j2-1,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(1,j,k,l) = ani*(at1 - at2)
      f(1,nx+1-j,k,l) = ani*(at1 + at2)
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(2,j,k,l) = ani*(at1 + at2)
      f(2,nx+1-j,k,l) = ani*(at1 - at2)
      at2 = f(3,nx+1-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(3,j,k,l) = ani*(at1 + at2)
      f(3,nx+1-j,k,l) = ani*(at1 - at2)
  200 continue
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = 0.0
      f(3,nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
c forward fourier transform
c create auxiliary array in x
  230 kmr = nxy/nx
      kypt = kyps
      do 260 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 250 k = kypi, kypt
      do 240 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(1,j,k,l) = at1 - at2
      f(1,nx+1-j,k,l) = at1 + at2
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(2,j,k,l) = at1 + at2
      f(2,nx+1-j,k,l) = at1 - at2
      at2 = f(3,nx+1-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(3,j,k,l) = at1 + at2
      f(3,nx+1-j,k,l) = at1 - at2
  240 continue
  250 continue
  260 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 300 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 290 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 290
      do 280 k = kypi, kypt
      do 270 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  270 continue
  280 continue
  290 continue
  300 continue
c first transform in x
      nrx = nxy/nxh
      do 360 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 350 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 320 i = kypi, kypt
      do 310 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 420 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 390 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 380 k = kypi, kypt
      do 370 jj = 1, 3
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
  370 continue
  380 continue
  390 continue
      do 410 k = kypi, kypt
      do 400 jj = 1, 3
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
      f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
      f(jj,1,k,l) = t2
      f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  400 continue
  410 continue
  420 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      kypt = kyps
      do 470 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 460 k = kypi, kypt
      sum1 = -0.5*f(1,2,k,l)
      do 430 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(1,2*j1-1,k,l)*at1 + f(1,2*j1,k,l)*at2
      at2 = -f(1,2*j1-1,k,l)*at2 + f(1,2*j1,k,l)*at1
      f(1,2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(1,2*j1,k,l) = at1
  430 continue
      f(1,2,k,l) = sum1
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = f(2,2,k,l)
      f(2,2,k,l) = 0.5*f(2,1,k,l)
      f(2,1,k,l) = 0.0
      sum1 = f(2,2,k,l)
      do 440 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2,2*j-1,k,l)*at2 - f(2,2*j,k,l)*at1
      at2 = f(2,2*j-1,k,l)*at1 + f(2,2*j,k,l)*at2
      f(2,2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(2,2*j,k,l) = sum1
  440 continue
      f(3,nx+1,k,l) = f(3,2,k,l)
      f(3,2,k,l) = 0.5*f(3,1,k,l)
      f(3,1,k,l) = 0.0
      sum1 = f(3,2,k,l)
      do 450 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(3,2*j-1,k,l)*at2 - f(3,2*j,k,l)*at1
      at2 = f(3,2*j-1,k,l)*at1 + f(3,2*j,k,l)*at2
      f(3,2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(3,2*j,k,l) = sum1
  450 continue
  460 continue
  470 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,ky
     1p,kypi,kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic, for data which is distributed in blocks
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)/nvp
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k,i) = (1/nx*ny)*(.5*f(1,nx+1,k,i)*(-1)**n + sum(f(1,j+1,k,i)*
c sin(pi*(n+1/2)*j/nx)))
c f(2,n,k,i) = (1/nx*ny)*(.5*f(2,1,k,i) + sum(f(2,j,k,i)*
c cos(pi*(n+1/2)*j/nx)))
c f(3,n,k,i) = (1/nx*ny)*(.5*f(3,1,k,i) + sum(f(3,j,k,i)*
c cos(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j+1,k,i) = 2.0*sum(f(1,n,k,i)*sin(pi*n*(j+1/2)/nx)))
c f(2,j,k,i) = 2.0*sum(f(2,n,k,i)*cos(pi*n*(j+1/2)/nx))
c f(3,j,k,i) = 2.0*sum(f(3,n,k,i)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
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
      complex sctd, sctdx
      dimension f(3,2*nxvh,kypd,kblok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nx3, ny, nxy, nxhy, ks, kypt
      integer i, j, k, l, m, km, kmr, nrx, j1, j2, ns, ns2, k1, k2, kyps
      integer jj
      real at1, at2, at3, at4, t2, t3, t4, t5, t6, ani, sum1, sum2
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
      if (isign.gt.0) go to 230
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      kypt = kyps
      do 30 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 20 k = kypi, kypt
      f(1,1,k,l) = 2.0*f(1,2,k,l)
      sum1 = f(2,nx,k,l)
      sum2 = f(3,nx,k,l)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k,l) - f(1,2*j1-2,k,l)
      at3 = f(1,2*j1-1,k,l)*at2 + at4*at1
      f(1,2*j1,k,l) = at4*at2 - f(1,2*j1-1,k,l)*at1
      f(1,2*j1-1,k,l) = at3
      at4 = f(2,2*j1,k,l) - f(2,2*j1-2,k,l)
      at3 = f(2,2*j1-1,k,l)*at1 + at4*at2
      f(2,2*j1,k,l) = f(2,2*j1-1,k,l)*at2 - at4*at1
      f(2,2*j1-1,k,l) = at3
      at4 = f(3,2*j1,k,l) - f(3,2*j1-2,k,l)
      at3 = f(3,2*j1-1,k,l)*at1 + at4*at2
      f(3,2*j1,k,l) = f(3,2*j1-1,k,l)*at2 - at4*at1
      f(3,2*j1-1,k,l) = at3
   10 continue
      f(1,2,k,l) = f(1,nx+1,k,l)
      f(2,2,k,l) = -2.0*sum1
      f(3,2,k,l) = -2.0*sum2
   20 continue
   30 continue
c scramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 90 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 60 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 50 k = kypi, kypt
      do 40 jj = 1, 3
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
   40 continue
   50 continue
   60 continue
      do 80 k = kypi, kypt
      do 70 jj = 1, 3
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = f(jj,1,k,l) + f(jj,2,k,l)
      f(jj,2,k,l) = f(jj,1,k,l) - f(jj,2,k,l)
      f(jj,1,k,l) = t2
   70 continue
   80 continue
   90 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 130 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 120 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 120
      do 110 k = kypi, kypt
      do 100 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  100 continue
  110 continue
  120 continue
  130 continue
c first transform in x
      nrx = nxy/nxh
      do 190 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 180 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 170 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 160 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 150 i = kypi, kypt
      do 140 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) + aimag(t1)*f(jj,2*j2,i,l)
      t3 = real(t1)*f(jj,2*j2,i,l) - aimag(t1)*f(jj,2*j2-1,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      kypt = kyps
      do 220 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 210 k = kypi, kypt
      do 200 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(1,j,k,l) = ani*(at1 + at2)
      f(1,nx+1-j,k,l) = ani*(at1 - at2)
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(2,j,k,l) = ani*(at1 - at2)
      f(2,nx+1-j,k,l) = ani*(at1 + at2)
      at2 = f(3,nx+1-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(3,j,k,l) = ani*(at1 - at2)
      f(3,nx+1-j,k,l) = ani*(at1 + at2)
  200 continue
      f(1,nx+1,k,l) = 0.0
      f(2,nx+1,k,l) = 0.0
      f(3,nx+1,k,l) = 0.0
  210 continue
  220 continue
      return
c forward fourier transform
c create auxiliary array in x
  230 kmr = nxy/nx
      kypt = kyps
      do 260 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 250 k = kypi, kypt
      do 240 j = 1, nxh
      at2 = f(1,nx+1-j,k,l)
      at1 = f(1,j,k,l) + at2
      at2 = f(1,j,k,l) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(1,j,k,l) = at1 + at2
      f(1,nx+1-j,k,l) = at1 - at2
      at2 = f(2,nx+1-j,k,l)
      at1 = f(2,j,k,l) + at2
      at2 = f(2,j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(2,j,k,l) = at1 - at2
      f(2,nx+1-j,k,l) = at1 + at2
      at2 = f(3,nx+1-j,k,l)
      at1 = f(3,j,k,l) + at2
      at2 = f(3,j,k,l) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(3,j,k,l) = at1 - at2
      f(3,nx+1-j,k,l) = at1 + at2
  240 continue
  250 continue
  260 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      kypt = kyps
      do 300 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 290 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 290
      do 280 k = kypi, kypt
      do 270 jj = 1, 3
      t2 = f(jj,2*j1-1,k,l)
      t3 = f(jj,2*j1,k,l)
      f(jj,2*j1-1,k,l) = f(jj,2*j-1,k,l)
      f(jj,2*j1,k,l) = f(jj,2*j,k,l)
      f(jj,2*j-1,k,l) = t2
      f(jj,2*j,k,l) = t3
  270 continue
  280 continue
  290 continue
  300 continue
c first transform in x
      nrx = nxy/nxh
      do 360 m = 1, indx1
      ns = 2**(m - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      kypt = kyps
      do 350 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 320 i = kypi, kypt
      do 310 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i,l) - aimag(t1)*f(jj,2*j2,i,l)
      t3 = aimag(t1)*f(jj,2*j2-1,i,l) + real(t1)*f(jj,2*j2,i,l)
      f(jj,2*j2-1,i,l) = f(jj,2*j1-1,i,l) - t2
      f(jj,2*j2,i,l) = f(jj,2*j1,i,l) - t3
      f(jj,2*j1-1,i,l) = f(jj,2*j1-1,i,l) + t2
      f(jj,2*j1,i,l) = f(jj,2*j1,i,l) + t3
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
c unscramble coefficients
      kmr = nxy/nxh
      kypt = kyps
      do 420 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 390 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 380 k = kypi, kypt
      do 370 jj = 1, 3
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
  370 continue
  380 continue
  390 continue
      do 410 k = kypi, kypt
      do 400 jj = 1, 3
      f(jj,nxh+1,k,l) = 2.0*f(jj,nxh+1,k,l)
      f(jj,nxh+2,k,l) = -2.0*f(jj,nxh+2,k,l)
      t2 = 2.0*(f(jj,1,k,l) + f(jj,2,k,l))
      f(jj,2,k,l) = 2.0*(f(jj,1,k,l) - f(jj,2,k,l))
      f(jj,1,k,l) = t2
      f(jj,nx+1,k,l) = 2.0*f(jj,nx+1,k,l)
  400 continue
  410 continue
  420 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      kypt = kyps
      do 470 l = 1, kblok
      if ((kyps+kyp*(l+ks)).eq.ny) kypt = kyps + 1
      do 460 k = kypi, kypt
      f(1,nx+1,k,l) = f(1,2,k,l)
      f(1,2,k,l) = 0.5*f(1,1,k,l)
      f(1,1,k,l) = 0.0
      sum1 = f(1,2,k,l)
      do 430 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(1,2*j-1,k,l)*at2 - f(1,2*j,k,l)*at1
      at2 = f(1,2*j-1,k,l)*at1 + f(1,2*j,k,l)*at2
      f(1,2*j-1,k,l) = at3
      sum1 = sum1 + at2
      f(1,2*j,k,l) = sum1
  430 continue
      sum1 = -0.5*f(2,2,k,l)
      do 440 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2,2*j1-1,k,l)*at1 + f(2,2*j1,k,l)*at2
      at2 = -f(2,2*j1-1,k,l)*at2 + f(2,2*j1,k,l)*at1
      f(2,2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2,2*j1,k,l) = at1
  440 continue
      f(2,2,k,l) = sum1
      f(2,nx+1,k,l) = 0.0
      sum1 = -0.5*f(3,2,k,l)
      do 450 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(3,2*j1-1,k,l)*at1 + f(3,2*j1,k,l)*at2
      at2 = -f(3,2*j1-1,k,l)*at2 + f(3,2*j1,k,l)*at1
      f(3,2*j1-1,k,l) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(3,2*j1,k,l) = at1
  450 continue
      f(3,2,k,l) = sum1
      f(3,nx+1,k,l) = 0.0
  460 continue
  470 continue
      return
      end
