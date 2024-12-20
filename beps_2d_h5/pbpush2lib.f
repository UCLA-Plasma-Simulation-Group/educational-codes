c 2d parallel PIC library for pushing particles with magnetic field
c and depositing current
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: june 9, 2008
c-----------------------------------------------------------------------
      subroutine PJDOST2(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npmax,
     1nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries
c and distributed data.
c in addition, particle positions are advanced a half time-step
c baseline scalar distributed version
c 76 flops/particle, 32 loads, 29 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cui(n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cui(n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cui(n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cui(n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cui(n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cui(n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cui(n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cui(n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cui(j,k,l) = ith component of current density at grid point (j,kk),
c where n = j + nxv*(kk - 1) and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), cux(nxv,nypmx,nblok)
      dimension cuy(nxv,nypmx,nblok), cuz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qmh = .5*qm
      do 20 l = 1, nblok
      mnoff = noff(l) - 2
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      amx = qm*(.75 - dxp*dxp)
      mm = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = qmh*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = qmh*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      cux(nl,mm,l) = cux(nl,mm,l) + vx*dx
      cux(nn,mm,l) = cux(nn,mm,l) + vx*dy
      cux(np,mm,l) = cux(np,mm,l) + vx*dz
      cuy(nl,mm,l) = cuy(nl,mm,l) + vy*dx
      cuy(nn,mm,l) = cuy(nn,mm,l) + vy*dy
      cuy(np,mm,l) = cuy(np,mm,l) + vy*dz
      cuz(nl,mm,l) = cuz(nl,mm,l) + vz*dx
      cuz(nn,mm,l) = cuz(nn,mm,l) + vz*dy
      cuz(np,mm,l) = cuz(np,mm,l) + vz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cux(nl,ml,l) = cux(nl,ml,l) + vx*dx
      cux(nn,ml,l) = cux(nn,ml,l) + vx*dy
      cux(np,ml,l) = cux(np,ml,l) + vx*dz
      cuy(nl,ml,l) = cuy(nl,ml,l) + vy*dx
      cuy(nn,ml,l) = cuy(nn,ml,l) + vy*dy
      cuy(np,ml,l) = cuy(np,ml,l) + vy*dz
      cuz(nl,ml,l) = cuz(nl,ml,l) + vz*dx
      cuz(nn,ml,l) = cuz(nn,ml,l) + vz*dy
      cuz(np,ml,l) = cuz(np,ml,l) + vz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cux(nl,mp,l) = cux(nl,mp,l) + vx*dx
      cux(nn,mp,l) = cux(nn,mp,l) + vx*dy
      cux(np,mp,l) = cux(np,mp,l) + vx*dz
      cuy(nl,mp,l) = cuy(nl,mp,l) + vy*dx
      cuy(nn,mp,l) = cuy(nn,mp,l) + vy*dy
      cuy(np,mp,l) = cuy(np,mp,l) + vy*dz
      cuz(nl,mp,l) = cuz(nl,mp,l) + vz*dx
      cuz(nn,mp,l) = cuz(nn,mp,l) + vz*dy
      cuz(np,mp,l) = cuz(np,mp,l) + vz*dz
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok
     1,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, for distributed data
c 76 flops/particle, 32 loads, 29 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      ml = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + vx*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + vy*dx
      cu(3,nl,mm,l) = cu(3,nl,mm,l) + vz*dx
      dx = dxl*dyl
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
      cu(3,nn,mm,l) = cu(3,nn,mm,l) + vz*dy
      dy = amx*dyl
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dz
      cu(3,np,mm,l) = cu(3,np,mm,l) + vz*dz
      dz = dxp*dyl
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + vx*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + vy*dx
      cu(3,nl,ml,l) = cu(3,nl,ml,l) + vz*dx
      dx = dxl*dyp
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + vx*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + vy*dy
      cu(3,nn,ml,l) = cu(3,nn,ml,l) + vz*dy
      dy = amx*dyp
      cu(1,np,ml,l) = cu(1,np,ml,l) + vx*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + vy*dz
      cu(3,np,ml,l) = cu(3,np,ml,l) + vz*dz
      dz = dxp*dyp
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + vx*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + vy*dx
      cu(3,nl,mp,l) = cu(3,nl,mp,l) + vz*dx
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      cu(3,nn,mp,l) = cu(3,nn,mp,l) + vz*dy
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dz
      cu(3,np,mp,l) = cu(3,np,mp,l) + vz*dz
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblo
     1k,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 76 flops/particle, 32 loads, 29 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyp = first actual dimension of current array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l) + .5
      mmn = part(2,j,l) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      ml = mm + nn
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      vz = part(5,j-1,l)
      dx1 = cu(1,mn,l) + vx*dx
      dy1 = cu(2,mn,l) + vy*dx
      amy = cu(3,mn,l) + vz*dx
      dx2 = cu(1,mn+1,l) + vx*dy
      dy2 = cu(2,mn+1,l) + vy*dy
      dx3 = cu(3,mn+1,l) + vz*dy
      dx = cu(1,mn+2,l) + vx*dz
      dy = cu(2,mn+2,l) + vy*dz
      dz = cu(3,mn+2,l) + vz*dz
      cu(1,mn,l) = dx1
      cu(2,mn,l) = dy1
      cu(3,mn,l) = amy
      cu(1,mn+1,l) = dx2
      cu(2,mn+1,l) = dy2
      cu(3,mn+1,l) = dx3
      cu(1,mn+2,l) = dx
      cu(2,mn+2,l) = dy
      cu(3,mn+2,l) = dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = cu(1,ml,l) + vx*dx
      dy1 = cu(2,ml,l) + vy*dx
      amy = cu(3,ml,l) + vz*dx
      dx2 = cu(1,ml+1,l) + vx*dy
      dy2 = cu(2,ml+1,l) + vy*dy
      dyl = cu(3,ml+1,l) + vz*dy
      dx = cu(1,ml+2,l) + vx*dz
      dy = cu(2,ml+2,l) + vy*dz
      dz = cu(3,ml+2,l) + vz*dz
      cu(1,ml,l) = dx1
      cu(2,ml,l) = dy1
      cu(3,ml,l) = amy
      cu(1,ml+1,l) = dx2
      cu(2,ml+1,l) = dy2
      cu(3,ml+1,l) = dyl
      cu(1,ml+2,l) = dx
      cu(2,ml+2,l) = dy
      cu(3,ml+2,l) = dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = cu(1,mp,l) + vx*dx
      dy1 = cu(2,mp,l) + vy*dx
      amy = cu(3,mp,l) + vz*dx
      dxl = cu(1,mp+1,l) + vx*dy
      amx = cu(2,mp+1,l) + vy*dy
      dxp = cu(3,mp+1,l) + vz*dy
      dx = cu(1,mp+2,l) + vx*dz
      dy = cu(2,mp+2,l) + vy*dz
      dz = cu(3,mp+2,l) + vz*dz
      cu(1,mp,l) = dx1
      cu(2,mp,l) = dy1
      cu(3,mp,l) = amy
      cu(1,mp+1,l) = dxl
      cu(2,mp+1,l) = amx
      cu(3,mp+1,l) = dxp
      cu(1,mp+2,l) = dx
      cu(2,mp+2,l) = dy
      cu(3,mp+2,l) = dz
c advance position half a time-step
      dx = part(1,j-1,l) + vx*dt
      dy = part(2,j-1,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,l)
            part(4,j-1,l) = -part(4,j-1,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
      endif
c set new position
      part(1,j-1,l) = dx
      part(2,j-1,l) = dy
   10 continue
      nop = npp(l)
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      cu(1,mn,l) = cu(1,mn,l) + vx*dx
      cu(2,mn,l) = cu(2,mn,l) + vy*dx
      cu(3,mn,l) = cu(3,mn,l) + vz*dx
      cu(1,mn+1,l) = cu(1,mn+1,l) + vx*dy
      cu(2,mn+1,l) = cu(2,mn+1,l) + vy*dy
      cu(3,mn+1,l) = cu(3,mn+1,l) + vz*dy
      cu(1,mn+2,l) = cu(1,mn+2,l) + vx*dz
      cu(2,mn+2,l) = cu(2,mn+2,l) + vy*dz
      cu(3,mn+2,l) = cu(3,mn+2,l) + vz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml,l) = cu(1,ml,l) + vx*dx
      cu(2,ml,l) = cu(2,ml,l) + vy*dx
      cu(3,ml,l) = cu(3,ml,l) + vz*dx
      cu(1,ml+1,l) = cu(1,ml+1,l) + vx*dy
      cu(2,ml+1,l) = cu(2,ml+1,l) + vy*dy
      cu(3,ml+1,l) = cu(3,ml+1,l) + vz*dy
      cu(1,ml+2,l) = cu(1,ml+2,l) + vx*dz
      cu(2,ml+2,l) = cu(2,ml+2,l) + vy*dz
      cu(3,ml+2,l) = cu(3,ml+2,l) + vz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp,l) = cu(1,mp,l) + vx*dx
      cu(2,mp,l) = cu(2,mp,l) + vy*dx
      cu(3,mp,l) = cu(3,mp,l) + vz*dx
      cu(1,mp+1,l) = cu(1,mp+1,l) + vx*dy
      cu(2,mp+1,l) = cu(2,mp+1,l) + vy*dy
      cu(3,mp+1,l) = cu(3,mp+1,l) + vz*dy
      cu(1,mp+2,l) = cu(1,mp+2,l) + vx*dz
      cu(2,mp+2,l) = cu(2,mp+2,l) + vy*dz
      cu(3,mp+2,l) = cu(3,mp+2,l) + vz*dz
c advance position half a time-step
      dx = part(1,nop,l) + vx*dt
      dy = part(2,nop,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,idimp,npmax,
     1nblok,nxv,nxvyp,npd,n27)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with 1d addressing
c 76 flops/particle, 87 loads, 83 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,i)=qci*(.75-dx**2)*(.75-dy**2)
c cu(n+1,m,i)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(n-1,m,i)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(n,m+1,i)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(n+1,m+1,i)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(n-1,m+1,i)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(n,m-1,i)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(n+1,m-1,i)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(n-1,m-1,i)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(n,i,l) = ith component of current density at grid point (j,kk),
c where n = j + nxv*(kk - 1) and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of current array
c npd = size of scratch buffers for vectorized current deposition
c n27 = number of independent weights
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n27,npd,nblok), amxy(n27,npd,nblok)
      nxvyp2 = nxvyp + nxvyp
      zero = 0.
      anx = float(nx)
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l) - 1
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l) + .5
      dxl = part(1,i+jb,l) - float(n)
      n = n + 1
      if (n.gt.nx) n = n - nx
      amx = qm*(.75 - dxl*dxl)
      np = n + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*qm*(.5 + dxl)**2
      nl = n - 1
      if (nl.lt.1) nl = nl + nx
      dxl = .5*qm*(.5 - dxl)**2
      m = part(2,i+jb,l) + .5
      dyl = part(2,i+jb,l) - float(m)
      m = nxv*(m - mnoff)
      amy = .75 - dyl*dyl
      mp = m + nxv
      dyp = .5*(.5 + dyl)**2
      ml = m - nxv
      dyl = .5*(.5 - dyl)**2
      in = n + m
      nn(1,i,l) = in
      nn(1+9,i,l) = in + nxvyp
      nn(1+18,i,l) = in + nxvyp2
      in = np + m
      nn(2,i,l) = in
      nn(2+9,i,l) = in + nxvyp
      nn(2+18,i,l) = in + nxvyp2
      in = nl + m
      nn(3,i,l) = in
      nn(3+9,i,l) = in + nxvyp
      nn(3+18,i,l) = in + nxvyp2
      in = n + mp
      nn(4,i,l) = in
      nn(4+9,i,l) = in + nxvyp
      nn(4+18,i,l) = in + nxvyp2
      in = np + mp
      nn(5,i,l) = in
      nn(5+9,i,l) = in + nxvyp
      nn(5+18,i,l) = in + nxvyp2
      in = nl + mp
      nn(6,i,l) = in
      nn(6+9,i,l) = in + nxvyp
      nn(6+18,i,l) = in + nxvyp2
      in = n + ml
      nn(7,i,l) = in
      nn(7+9,i,l) = in + nxvyp
      nn(7+18,i,l) = in + nxvyp2
      in = np + ml
      nn(8,i,l) = in
      nn(8+9,i,l) = in + nxvyp
      nn(8+18,i,l) = in + nxvyp2
      in = nl + ml
      nn(9,i,l) = in
      nn(9+9,i,l) = in + nxvyp
      nn(9+18,i,l) = in + nxvyp2
      dx = amx*amy
      dy = dxp*amy
      dz = dxl*amy
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      vz = part(5,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vx*dy
      amxy(3,i,l) = vx*dz
      amxy(1+9,i,l) = vy*dx
      amxy(2+9,i,l) = vy*dy
      amxy(3+9,i,l) = vy*dz
      amxy(1+18,i,l) = vz*dx
      amxy(2+18,i,l) = vz*dy
      amxy(3+18,i,l) = vz*dz
      dx = amx*dyp
      dy = dxp*dyp
      dz = dxl*dyp
      amxy(4,i,l) = vx*dx
      amxy(5,i,l) = vx*dy
      amxy(6,i,l) = vx*dz
      amxy(4+9,i,l) = vy*dx
      amxy(5+9,i,l) = vy*dy
      amxy(6+9,i,l) = vy*dz
      amxy(4+18,i,l) = vz*dx
      amxy(5+18,i,l) = vz*dy
      amxy(6+18,i,l) = vz*dz
      dx = amx*dyl
      dy = dxp*dyl
      dz = dxl*dyl
      amxy(7,i,l) = vx*dx
      amxy(8,i,l) = vx*dy
      amxy(9,i,l) = vx*dz
      amxy(7+9,i,l) = vy*dx
      amxy(8+9,i,l) = vy*dy
      amxy(9+9,i,l) = vy*dz
      amxy(7+18,i,l) = vz*dx
      amxy(8+18,i,l) = vz*dy
      amxy(9+18,i,l) = vz*dz
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 27
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,np
     1max,nblok,nxv,nxvyp,npd,n27,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 76 flops/particle, 87 loads, 83 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,n,l) = ith component of current density at grid point (j,kk),
c where n = j + nxv*kk + 1 and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of current array
c npd = size of scratch buffers for vectorized current deposition
c n27 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n27,npd,nblok), amxy(n27,npd,nblok)
      nxv3 = 3*nxv
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l) + .5
      m = part(2,i+jb,l) + .5
      dxp = part(1,i+jb,l) - float(n)
      dyp = part(2,i+jb,l) - float(m)
      n3 = 3*n + 1
      m = nxv3*(m - mnoff)
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n3
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mn = ml + nxv3
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv3
      nn(1,i,l) = mn
      nn(2,i,l) = mn + 1
      nn(3,i,l) = mn + 2
      nn(4,i,l) = mn + 3
      nn(5,i,l) = mn + 4
      nn(6,i,l) = mn + 5
      nn(7,i,l) = mn + 6
      nn(8,i,l) = mn + 7
      nn(9,i,l) = mn + 8
      nn(10,i,l) = ml
      nn(11,i,l) = ml + 1
      nn(12,i,l) = ml + 2
      nn(13,i,l) = ml + 3
      nn(14,i,l) = ml + 4
      nn(15,i,l) = ml + 5
      nn(16,i,l) = ml + 6
      nn(17,i,l) = ml + 7
      nn(18,i,l) = ml + 8
      nn(19,i,l) = mp
      nn(20,i,l) = mp + 1
      nn(21,i,l) = mp + 2
      nn(22,i,l) = mp + 3
      nn(23,i,l) = mp + 4
      nn(24,i,l) = mp + 5
      nn(25,i,l) = mp + 6
      nn(26,i,l) = mp + 7
      nn(27,i,l) = mp + 8
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      vz = part(5,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vy*dx
      amxy(3,i,l) = vz*dx
      dx = dxl*dyl
      amxy(4,i,l) = vx*dy
      amxy(5,i,l) = vy*dy
      amxy(6,i,l) = vz*dy
      dy = amx*dyl
      amxy(7,i,l) = vx*dz
      amxy(8,i,l) = vy*dz
      amxy(9,i,l) = vz*dz
      dz = dxp*dyl
      amxy(10,i,l) = vx*dx
      amxy(11,i,l) = vy*dx
      amxy(12,i,l) = vz*dx
      dx = dxl*dyp
      amxy(13,i,l) = vx*dy
      amxy(14,i,l) = vy*dy
      amxy(15,i,l) = vz*dy
      dy = amx*dyp
      amxy(16,i,l) = vx*dz
      amxy(17,i,l) = vy*dz
      amxy(18,i,l) = vz*dz
      dz = dxp*dyp
      amxy(19,i,l) = vx*dx
      amxy(20,i,l) = vy*dx
      amxy(21,i,l) = vz*dx
      amxy(22,i,l) = vx*dy
      amxy(23,i,l) = vy*dy
      amxy(24,i,l) = vz*dy
      amxy(25,i,l) = vx*dz
      amxy(26,i,l) = vy*dz
      amxy(27,i,l) = vz*dz
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,l)
            part(4,i+jb,l) = -part(4,i+jb,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
      endif
c set new position
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 27
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PJDOST2L(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npmax
     1,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries
c and distributed data.
c in addition, particle positions are advanced a half time-step
c baseline scalar distributed version
c 35 flops/particle, 17 loads, 14 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m)=qci*(1.-dx)*(1.-dy)
c cui(n+1,m)=qci*dx*(1.-dy)
c cui(n,m+1)=qci*(1.-dx)*dy
c cui(n+1,m+1)=qci*dx*dy
c and qci = qm*vi, where i = x,y,z
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cui(j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), cux(nxv,nypmx,nblok)
      dimension cuy(nxv,nypmx,nblok), cuz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = qm*(part(1,j,l) - float(nn))
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      cux(nn,mm,l) = cux(nn,mm,l) + vx*dx
      cux(np,mm,l) = cux(np,mm,l) + vx*dy
      cuy(nn,mm,l) = cuy(nn,mm,l) + vy*dx
      cuy(np,mm,l) = cuy(np,mm,l) + vy*dy
      cuz(nn,mm,l) = cuz(nn,mm,l) + vz*dx
      cuz(np,mm,l) = cuz(np,mm,l) + vz*dy
      dx = amx*dyp
      dy = dxp*dyp
      cux(nn,mp,l) = cux(nn,mp,l) + vx*dx
      cux(np,mp,l) = cux(np,mp,l) + vx*dy
      cuy(nn,mp,l) = cuy(nn,mp,l) + vy*dx
      cuy(np,mp,l) = cuy(np,mp,l) + vy*dy
      cuz(nn,mp,l) = cuz(nn,mp,l) + vz*dx
      cuz(np,mp,l) = cuz(np,mp,l) + vz*dy
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblo
     1k,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, for distributed data
c 35 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = qm*(part(1,j,l) - float(nn))
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dx
      cu(3,np,mp,l) = cu(3,np,mp,l) + vz*dx
      dx = dxp*amy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      cu(3,nn,mp,l) = cu(3,nn,mp,l) + vz*dy
      dy = amx*amy
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dx
      cu(3,np,mm,l) = cu(3,np,mm,l) + vz*dx
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
      cu(3,nn,mm,l) = cu(3,nn,mm,l) + vz*dy
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nbl
     1ok,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 35 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l)
      mmn = part(2,j,l)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      vz = part(5,j-1,l)
      dx1 = cu(1,mp+1,l) + vx*dx
      dy1 = cu(2,mp+1,l) + vy*dx
      dyp = cu(3,mp+1,l) + vz*dx
      dx = cu(1,mp,l) + vx*dz
      dy = cu(2,mp,l) + vy*dz
      dz = cu(3,mp,l) + vz*dz
      cu(1,mp+1,l) = dx1
      cu(2,mp+1,l) = dy1
      cu(3,mp+1,l) = dyp
      cu(1,mp,l) = dx
      cu(2,mp,l) = dy
      cu(3,mp,l) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1,l) + vx*dx
      amx = cu(2,mm+1,l) + vy*dx
      dyp = cu(3,mm+1,l) + vz*dx
      dx = cu(1,mm,l) + vx*dz
      dy = cu(2,mm,l) + vy*dz
      dz = cu(3,mm,l) + vz*dz
      cu(1,mm+1,l) = dxp
      cu(2,mm+1,l) = amx
      cu(3,mm+1,l) = dyp
      cu(1,mm,l) = dx
      cu(2,mm,l) = dy
      cu(3,mm,l) = dz
c advance position half a time-step
      dx = part(1,j-1,l) + vx*dt
      dy = part(2,j-1,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,l)
            part(4,j-1,l) = -part(4,j-1,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
      endif
c set new position
      part(1,j-1,l) = dx
      part(2,j-1,l) = dy
   10 continue
      nop = npp(l)
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      cu(1,mp+1,l) = cu(1,mp+1,l) + vx*dx
      cu(2,mp+1,l) = cu(2,mp+1,l) + vy*dx
      cu(3,mp+1,l) = cu(3,mp+1,l) + vz*dx
      cu(1,mp,l) = cu(1,mp,l) + vx*dy
      cu(2,mp,l) = cu(2,mp,l) + vy*dy
      cu(3,mp,l) = cu(3,mp,l) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1,l) = cu(1,mm+1,l) + vx*dx
      cu(2,mm+1,l) = cu(2,mm+1,l) + vy*dx
      cu(3,mm+1,l) = cu(3,mm+1,l) + vz*dx
      cu(1,mm,l) = cu(1,mm,l) + vx*dy
      cu(2,mm,l) = cu(2,mm,l) + vy*dy
      cu(3,mm,l) = cu(3,mm,l) + vz*dy
c advance position half a time-step
      dx = part(1,nop,l) + vx*dt
      dy = part(2,nop,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,idimp,npmax
     1,nblok,nxv,nxvyp,npd,n12)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with 1d addressing
c 35 flops/particle, 41 loads, 38 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,i)=qci*(1.-dx)*(1.-dy)
c cu(n+1,m,i)=qci*dx*(1.-dy)
c cu(n,m+1,i)=qci*(1.-dx)*dy
c cu(n+1,m+1,i)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(n,i,l) = ith component of current density at grid point (j,kk),
c where n = j + nxv*(kk - 1) and  kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized current deposition
c amxy = scratch weight array for vectorized current deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of current array
c npd = size of scratch buffers for vectorized current deposition
c n12 = number of independent weights
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n12,npd,nblok), amxy(n12,npd,nblok)
      nxvyp2 = nxvyp + nxvyp
      zero = 0.
      anx = float(nx)
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l)
      dxp = qm*(part(1,i+jb,l) - float(n))
      n = n + 1
      amx = qm - dxp
      np = n + 1
      if (np.gt.nx) np = np - nx
      m = part(2,i+jb,l)
      dyp = part(2,i+jb,l) - float(m)
      m = nxv*(m - mnoff)
      amy = 1. - dyp
      mp = m + nxv
      in = n + m
      nn(1,i,l) = in
      nn(1+4,i,l) = in + nxvyp
      nn(1+8,i,l) = in + nxvyp2
      in = np + m
      nn(2,i,l) = in
      nn(2+4,i,l) = in + nxvyp
      nn(2+8,i,l) = in + nxvyp2
      in = n + mp
      nn(3,i,l) = in
      nn(3+4,i,l) = in + nxvyp
      nn(3+8,i,l) = in + nxvyp2
      in = np + mp
      nn(4,i,l) = in
      nn(4+4,i,l) = in + nxvyp
      nn(4+8,i,l) = in + nxvyp2
      dx = amx*amy
      dy = dxp*amy
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      vz = part(5,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vx*dy
      amxy(1+4,i,l) = vy*dx
      amxy(2+4,i,l) = vy*dy
      amxy(1+8,i,l) = vz*dx
      amxy(2+8,i,l) = vz*dy
      dx = amx*dyp
      dy = dxp*dyp
      amxy(3,i,l) = vx*dx
      amxy(4,i,l) = vx*dy
      amxy(3+4,i,l) = vy*dx
      amxy(4+4,i,l) = vy*dy
      amxy(3+8,i,l) = vz*dx
      amxy(4+8,i,l) = vz*dy
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 12
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,n
     1pmax,nblok,nxv,nxvyp,npd,n12,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with guard cells and 1d addressing
c 35 flops/particle, 41 loads, 38 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c part(5,n,l) = z velocity of particle n in partition l
c cu(i,n,l) = charge density at grid point (j,kk),
c where n = j + nxv*(kk - 1) and  kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c n12 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n12,npd,nblok), amxy(n12,npd,nblok)
      nxv3 = 3*nxv
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l)
      m = part(2,i+jb,l)
      dxp = qm*(part(1,i+jb,l) - float(n))
      dyp = part(2,i+jb,l) - float(m)
      n3 = 3*n + 1
      mm = nxv3*(m - mnoff)
      amx = qm - dxp
      mm = mm + n3
      amy = 1. - dyp
      mp = mm + nxv3
      nn(10,i,l) = mm
      nn(11,i,l) = mm + 1
      nn(12,i,l) = mm + 2
      nn(7,i,l) = mm + 3
      nn(8,i,l) = mm + 4
      nn(9,i,l) = mm + 5
      nn(4,i,l) = mp
      nn(5,i,l) = mp + 1
      nn(6,i,l) = mp + 2
      nn(1,i,l) = mp + 3
      nn(2,i,l) = mp + 4
      nn(3,i,l) = mp + 5
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      vz = part(5,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vy*dx
      amxy(3,i,l) = vz*dx
      amxy(4,i,l) = vx*dy
      amxy(5,i,l) = vy*dy
      amxy(6,i,l) = vz*dy
      dx = dxp*amy
      dy = amx*amy
      amxy(7,i,l) = vx*dx
      amxy(8,i,l) = vy*dx
      amxy(9,i,l) = vz*dx
      amxy(10,i,l) = vx*dy
      amxy(11,i,l) = vy*dy
      amxy(12,i,l) = vz*dy
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,l)
            part(4,i+jb,l) = -part(4,i+jb,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
      endif
c set new position
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 12
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblo
     1k,nxv,nypmx,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, for distributed data
c 58 flops/particle, 22 loads, 20 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      ml = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + vx*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + vy*dx
      dx = dxl*dyl
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
      dy = amx*dyl
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dz
      dz = dxp*dyl
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + vx*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + vy*dx
      dx = dxl*dyp
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + vx*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + vy*dy
      dy = amx*dyp
      cu(1,np,ml,l) = cu(1,np,ml,l) + vx*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + vy*dz
      dz = dxp*dyp
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + vx*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + vy*dx
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dz
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nbl
     1ok,nxv,nxyp,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 58 flops/particle, 22 loads, 20 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyp = first actual dimension of current array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l) + .5
      mmn = part(2,j,l) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      ml = mm + nn
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      dx1 = cu(1,mn,l) + vx*dx
      dy1 = cu(2,mn,l) + vy*dx
      dx2 = cu(1,mn+1,l) + vx*dy
      dy2 = cu(2,mn+1,l) + vy*dy
      dx = cu(1,mn+2,l) + vx*dz
      dy = cu(2,mn+2,l) + vy*dz
      cu(1,mn,l) = dx1
      cu(2,mn,l) = dy1
      cu(1,mn+1,l) = dx2
      cu(2,mn+1,l) = dy2
      cu(1,mn+2,l) = dx
      cu(2,mn+2,l) = dy
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = cu(1,ml,l) + vx*dx
      dy1 = cu(2,ml,l) + vy*dx
      dx2 = cu(1,ml+1,l) + vx*dy
      dy2 = cu(2,ml+1,l) + vy*dy
      dx = cu(1,ml+2,l) + vx*dz
      dy = cu(2,ml+2,l) + vy*dz
      cu(1,ml,l) = dx1
      cu(2,ml,l) = dy1
      cu(1,ml+1,l) = dx2
      cu(2,ml+1,l) = dy2
      cu(1,ml+2,l) = dx
      cu(2,ml+2,l) = dy
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = cu(1,mp,l) + vx*dx
      dy1 = cu(2,mp,l) + vy*dx
      dxl = cu(1,mp+1,l) + vx*dy
      amx = cu(2,mp+1,l) + vy*dy
      dx = cu(1,mp+2,l) + vx*dz
      dy = cu(2,mp+2,l) + vy*dz
      cu(1,mp,l) = dx1
      cu(2,mp,l) = dy1
      cu(1,mp+1,l) = dxl
      cu(2,mp+1,l) = amx
      cu(1,mp+2,l) = dx
      cu(2,mp+2,l) = dy
c advance position half a time-step
      dx = part(1,j-1,l) + vx*dt
      dy = part(2,j-1,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,l)
            part(4,j-1,l) = -part(4,j-1,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
      endif
c set new position
      part(1,j-1,l) = dx
      part(2,j-1,l) = dy
   10 continue
      nop = npp(l)
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      cu(1,mn,l) = cu(1,mn,l) + vx*dx
      cu(2,mn,l) = cu(2,mn,l) + vy*dx
      cu(1,mn+1,l) = cu(1,mn+1,l) + vx*dy
      cu(2,mn+1,l) = cu(2,mn+1,l) + vy*dy
      cu(1,mn+2,l) = cu(1,mn+2,l) + vx*dz
      cu(2,mn+2,l) = cu(2,mn+2,l) + vy*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml,l) = cu(1,ml,l) + vx*dx
      cu(2,ml,l) = cu(2,ml,l) + vy*dx
      cu(1,ml+1,l) = cu(1,ml+1,l) + vx*dy
      cu(2,ml+1,l) = cu(2,ml+1,l) + vy*dy
      cu(1,ml+2,l) = cu(1,ml+2,l) + vx*dz
      cu(2,ml+2,l) = cu(2,ml+2,l) + vy*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp,l) = cu(1,mp,l) + vx*dx
      cu(2,mp,l) = cu(2,mp,l) + vy*dx
      cu(1,mp+1,l) = cu(1,mp+1,l) + vx*dy
      cu(2,mp+1,l) = cu(2,mp+1,l) + vy*dy
      cu(1,mp+2,l) = cu(1,mp+2,l) + vx*dz
      cu(2,mp+2,l) = cu(2,mp+2,l) + vy*dz
c advance position half a time-step
      dx = part(1,nop,l) + vx*dt
      dy = part(2,nop,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,n
     1pmax,nblok,nxv,nxvyp,npd,n18,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 58 flops/particle, 58 loads, 56 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,n,l) = ith component of current density at grid point (j,kk),
c where n = j + nxv*kk + 1 and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of current array
c npd = size of scratch buffers for vectorized current deposition
c n18 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n18,npd,nblok), amxy(n18,npd,nblok)
      nxv2 = 2*nxv
      qmh = .5*qm
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l) + .5
      m = part(2,i+jb,l) + .5
      dxp = part(1,i+jb,l) - float(n)
      dyp = part(2,i+jb,l) - float(m)
      n2 = 2*n + 1
      m = nxv2*(m - mnoff)
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n2
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mn = ml + nxv2
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv2
      nn(1,i,l) = mn
      nn(2,i,l) = mn + 1
      nn(3,i,l) = mn + 2
      nn(4,i,l) = mn + 3
      nn(5,i,l) = mn + 4
      nn(6,i,l) = mn + 5
      nn(7,i,l) = ml
      nn(8,i,l) = ml + 1
      nn(9,i,l) = ml + 2
      nn(10,i,l) = ml + 3
      nn(11,i,l) = ml + 4
      nn(12,i,l) = ml + 5
      nn(13,i,l) = mp
      nn(14,i,l) = mp + 1
      nn(15,i,l) = mp + 2
      nn(16,i,l) = mp + 3
      nn(17,i,l) = mp + 4
      nn(18,i,l) = mp + 5
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vy*dx
      dx = dxl*dyl
      amxy(3,i,l) = vx*dy
      amxy(4,i,l) = vy*dy
      dy = amx*dyl
      amxy(5,i,l) = vx*dz
      amxy(6,i,l) = vy*dz
      dz = dxp*dyl
      amxy(7,i,l) = vx*dx
      amxy(8,i,l) = vy*dx
      dx = dxl*dyp
      amxy(9,i,l) = vx*dy
      amxy(10,i,l) = vy*dy
      dy = amx*dyp
      amxy(11,i,l) = vx*dz
      amxy(12,i,l) = vy*dz
      dz = dxp*dyp
      amxy(13,i,l) = vx*dx
      amxy(14,i,l) = vy*dx
      amxy(15,i,l) = vx*dy
      amxy(16,i,l) = vy*dy
      amxy(17,i,l) = vx*dz
      amxy(18,i,l) = vy*dz
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,l)
            part(4,i+jb,l) = -part(4,i+jb,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
      endif
c set new position
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 18
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nbl
     1ok,nxv,nypmx,ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, for distributed data
c 27 flops/particle, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = qm*(part(1,j,l) - float(nn))
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j,l)
      vy = part(4,j,l)
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dx
      dx = dxp*amy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      dy = amx*amy
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dx
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
c advance position half a time-step
      dx = part(1,j,l) + vx*dt
      dy = part(2,j,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nb
     1lok,nxv,nxyp,ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation, and distributed data.
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 27 flops/particle, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qm*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qm*dx*(1.-dy)
c cu(i,n,m+1)=qm*(1.-dx)*dy
c cu(i,n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l)
      mmn = part(2,j,l)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      dx1 = cu(1,mp+1,l) + vx*dx
      dy1 = cu(2,mp+1,l) + vy*dx
      dx = cu(1,mp,l) + vx*dz
      dy = cu(2,mp,l) + vy*dz
      cu(1,mp+1,l) = dx1
      cu(2,mp+1,l) = dy1
      cu(1,mp,l) = dx
      cu(2,mp,l) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1,l) + vx*dx
      amx = cu(2,mm+1,l) + vy*dx
      dx = cu(1,mm,l) + vx*dz
      dy = cu(2,mm,l) + vy*dz
      cu(1,mm+1,l) = dxp
      cu(2,mm+1,l) = amx
      cu(1,mm,l) = dx
      cu(2,mm,l) = dy
c advance position half a time-step
      dx = part(1,j-1,l) + vx*dt
      dy = part(2,j-1,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,l)
            part(4,j-1,l) = -part(4,j-1,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,l)
            part(3,j-1,l) = -part(3,j-1,l)
         endif
      endif
c set new position
      part(1,j-1,l) = dx
      part(2,j-1,l) = dy
   10 continue
      nop = npp(l)
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      cu(1,mp+1,l) = cu(1,mp+1,l) + vx*dx
      cu(2,mp+1,l) = cu(2,mp+1,l) + vy*dx
      cu(1,mp,l) = cu(1,mp,l) + vx*dy
      cu(2,mp,l) = cu(2,mp,l) + vy*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1,l) = cu(1,mm+1,l) + vx*dx
      cu(2,mm+1,l) = cu(2,mm+1,l) + vy*dx
      cu(1,mm,l) = cu(1,mm,l) + vx*dy
      cu(2,mm,l) = cu(2,mm,l) + vy*dy
c advance position half a time-step
      dx = part(1,nop,l) + vx*dt
      dy = part(2,nop,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,
     1npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with guard cells and 1d addressing
c 27 flops/particle, 29 loads, 26 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = x velocity of particle n in partition l
c part(4,n,l) = y velocity of particle n in partition l
c cu(i,n,l) = charge density at grid point (j,kk),
c where n = j + nxv*(kk - 1) and  kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c n8 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(n8,npd,nblok), amxy(n8,npd,nblok)
      nxv2 = 2*nxv
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l)
      m = part(2,i+jb,l)
      dxp = qm*(part(1,i+jb,l) - float(n))
      dyp = part(2,i+jb,l) - float(m)
      n2 = 2*n + 1
      mm = nxv2*(m - mnoff)
      amx = qm - dxp
      mm = mm + n2
      amy = 1. - dyp
      mp = mm + nxv2
      nn(7,i,l) = mm
      nn(8,i,l) = mm + 1
      nn(5,i,l) = mm + 2
      nn(6,i,l) = mm + 3
      nn(3,i,l) = mp
      nn(4,i,l) = mp + 1
      nn(1,i,l) = mp + 2
      nn(2,i,l) = mp + 3
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,i+jb,l)
      vy = part(4,i+jb,l)
      amxy(1,i,l) = vx*dx
      amxy(2,i,l) = vy*dx
      amxy(3,i,l) = vx*dy
      amxy(4,i,l) = vy*dy
      dx = dxp*amy
      dy = amx*amy
      amxy(5,i,l) = vx*dx
      amxy(6,i,l) = vy*dx
      amxy(7,i,l) = vx*dy
      amxy(8,i,l) = vy*dy
c advance position half a time-step
      dx = part(1,i+jb,l) + vx*dt
      dy = part(2,i+jb,l) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,l)
            part(4,i+jb,l) = -part(4,i+jb,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,l)
            part(3,i+jb,l) = -part(3,i+jb,l)
         endif
      endif
c set new position
      part(1,i+jb,l) = dx
      part(2,i+jb,l) = dy
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 8
      cu(nn(k,i,l),l) = cu(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH2(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,idimp
     1,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover.
c baseline scalar distributed version.
c 186 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bx(j,k,l) = x component of magnetic field at grid (j,kk)
c by(j,k,l) = y component of magnetic field at grid (j,kk)
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      do 20 l = 1, nblok
      mnoff = noff(l) - 2
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      mm = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = .5*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fx(nl,mm,l) + amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + d
     1yl*(dxl*fx(nl,ml,l) + amx*fx(nn,ml,l) + dxp*fx(np,ml,l)) + dyp*(dx
     2l*fx(nl,mp,l) + amx*fx(nn,mp,l) + dxp*fx(np,mp,l))
      dy = amy*(dxl*fy(nl,mm,l) + amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + d
     1yl*(dxl*fy(nl,ml,l) + amx*fy(nn,ml,l) + dxp*fy(np,ml,l)) + dyp*(dx
     2l*fy(nl,mp,l) + amx*fy(nn,mp,l) + dxp*fy(np,mp,l))
c find magnetic field
      ox = amy*(dxl*bx(nl,mm,l) + amx*bx(nn,mm,l) + dxp*bx(np,mm,l)) + d
     1yl*(dxl*bx(nl,ml,l) + amx*bx(nn,ml,l) + dxp*bx(np,ml,l)) + dyp*(dx
     2l*bx(nl,mp,l) + amx*bx(nn,mp,l) + dxp*bx(np,mp,l))
      oy = amy*(dxl*by(nl,mm,l) + amx*by(nn,mm,l) + dxp*by(np,mm,l)) + d
     1yl*(dxl*by(nl,ml,l) + amx*by(nn,ml,l) + dxp*by(np,ml,l)) + dyp*(dx
     2l*by(nl,mp,l) + amx*by(nn,mp,l) + dxp*by(np,mp,l))
      oz = amy*(dxl*bz(nl,mm,l) + amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + d
     1yl*(dxl*bz(nl,ml,l) + amx*bz(nn,ml,l) + dxp*bz(np,ml,l)) + dyp*(dx
     2l*bz(nl,mp,l) + amx*bz(nn,mp,l) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,idim
     1p,npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field.  Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 186 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
c find magnetic field
      ox = amy*(dxl*bxy(1,nl,mm,l) + amx*bxy(1,nn,mm,l) + dxp*bxy(1,np,m
     1m,l)) + dyl*(dxl*bxy(1,nl,ml,l) + amx*bxy(1,nn,ml,l) + dxp*bxy(1,n
     2p,ml,l)) + dyp*(dxl*bxy(1,nl,mp,l) + amx*bxy(1,nn,mp,l) + dxp*bxy(
     31,np,mp,l))
      oy = amy*(dxl*bxy(2,nl,mm,l) + amx*bxy(2,nn,mm,l) + dxp*bxy(2,np,m
     1m,l)) + dyl*(dxl*bxy(2,nl,ml,l) + amx*bxy(2,nn,ml,l) + dxp*bxy(2,n
     2p,ml,l)) + dyp*(dxl*bxy(2,nl,mp,l) + amx*bxy(2,nn,mp,l) + dxp*bxy(
     32,np,mp,l))
      oz = amy*(dxl*bxy(3,nl,mm,l) + amx*bxy(3,nn,mm,l) + dxp*bxy(3,np,m
     1m,l)) + dyl*(dxl*bxy(3,nl,ml,l) + amx*bxy(3,nn,ml,l) + dxp*bxy(3,n
     2p,ml,l)) + dyp*(dxl*bxy(3,nl,mp,l) + amx*bxy(3,nn,mp,l) + dxp*bxy(
     33,np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 186 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l) + .5
      mmn = part(2,j+1,l) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      mn = ml + nxv
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mp = mn + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mmn = mmn - mnoff
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
      acz = part(5,nop,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,nop,l) = dx
      part(4,nop,l) = dy
      part(5,nop,l) = dz
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH2L(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,idim
     1p,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover.
c baseline scalar distributed version
c 105 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n),l = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bx(j,k,l) = x component of magnetic field at grid (j,kk)
c by(j,k,l) = y component of magnetic field at grid (j,kk)
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
c find electric field
      dx = amy*(amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + dyp*(amx*fx(nn,mp,l
     1) + dxp*fx(np,mp,l))
      dy = amy*(amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + dyp*(amx*fy(nn,mp,l
     1) + dxp*fy(np,mp,l))
c find magnetic field
      ox = amy*(amx*bx(nn,mm,l) + dxp*bx(np,mm,l)) + dyp*(amx*bx(nn,mp,l
     1) + dxp*bx(np,mp,l))
      oy = amy*(amx*by(nn,mm,l) + dxp*by(np,mm,l)) + dyp*(amx*by(nn,mp,l
     1) + dxp*by(np,mp,l))
      oz = amy*(amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + dyp*(amx*bz(nn,mp,l
     1) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 105 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp,l) + amx*bxy(1,nn,mp,l)) + amy*(dxp*bxy(
     11,np,mm,l) + amx*bxy(1,nn,mm,l))
      oy = dyp*(dxp*bxy(2,np,mp,l) + amx*bxy(2,nn,mp,l)) + amy*(dxp*bxy(
     12,np,mm,l) + amx*bxy(2,nn,mm,l))
      oz = dyp*(dxp*bxy(3,np,mp,l) + amx*bxy(3,nn,mp,l)) + amy*(dxp*bxy(
     13,np,mm,l) + amx*bxy(3,nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,id
     1imp,npmax,nblok,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 105 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l)
      mmn = part(2,j+1,l)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxp*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyp*(dxp*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxp*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyp*(dxp*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxp*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxn*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyn*(dxn*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxn*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyn*(dxn*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxn*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
      acz = part(5,nop,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,nop,l) = dx
      part(4,nop,l) = dy
      part(5,nop,l) = dz
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH2CQ(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,idi
     1mp,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions.  Using the Boris Mover with correction.
c baseline scalar distributed version
c 202 flops/particle, 2 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bx(j,k,l) = x component of magnetic field at grid (j,kk)
c by(j,k,l) = y component of magnetic field at grid (j,kk)
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      do 20 l = 1, nblok
      mnoff = noff(l) - 2
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      mm = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = .5*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fx(nl,mm,l) + amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + d
     1yl*(dxl*fx(nl,ml,l) + amx*fx(nn,ml,l) + dxp*fx(np,ml,l)) + dyp*(dx
     2l*fx(nl,mp,l) + amx*fx(nn,mp,l) + dxp*fx(np,mp,l))
      dy = amy*(dxl*fy(nl,mm,l) + amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + d
     1yl*(dxl*fy(nl,ml,l) + amx*fy(nn,ml,l) + dxp*fy(np,ml,l)) + dyp*(dx
     2l*fy(nl,mp,l) + amx*fy(nn,mp,l) + dxp*fy(np,mp,l))
c find magnetic field
      ox = amy*(dxl*bx(nl,mm,l) + amx*bx(nn,mm,l) + dxp*bx(np,mm,l)) + d
     1yl*(dxl*bx(nl,ml,l) + amx*bx(nn,ml,l) + dxp*bx(np,ml,l)) + dyp*(dx
     2l*bx(nl,mp,l) + amx*bx(nn,mp,l) + dxp*bx(np,mp,l))
      oy = amy*(dxl*by(nl,mm,l) + amx*by(nn,mm,l) + dxp*by(np,mm,l)) + d
     1yl*(dxl*by(nl,ml,l) + amx*by(nn,ml,l) + dxp*by(np,ml,l)) + dyp*(dx
     2l*by(nl,mp,l) + amx*by(nn,mp,l) + dxp*by(np,mp,l))
      oz = amy*(dxl*bz(nl,mm,l) + amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + d
     1yl*(dxl*bz(nl,ml,l) + amx*bz(nn,ml,l) + dxp*bz(np,ml,l)) + dyp*(dx
     2l*bz(nl,mp,l) + amx*bz(nn,mp,l) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,l) + (dx*dtt + omxt*omt)
      dy = part(2,j,l) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH2CQ(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idimp,
     1npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field.  Using the Boris Mover
c with correction.
c scalar version using guard cells, for distributed data
c 202 flops/particle, 2 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c nxv = first dimension of field arrays, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
c find magnetic field
      ox = amy*(dxl*bxy(1,nl,mm,l) + amx*bxy(1,nn,mm,l) + dxp*bxy(1,np,m
     1m,l)) + dyl*(dxl*bxy(1,nl,ml,l) + amx*bxy(1,nn,ml,l) + dxp*bxy(1,n
     2p,ml,l)) + dyp*(dxl*bxy(1,nl,mp,l) + amx*bxy(1,nn,mp,l) + dxp*bxy(
     31,np,mp,l))
      oy = amy*(dxl*bxy(2,nl,mm,l) + amx*bxy(2,nn,mm,l) + dxp*bxy(2,np,m
     1m,l)) + dyl*(dxl*bxy(2,nl,ml,l) + amx*bxy(2,nn,ml,l) + dxp*bxy(2,n
     2p,ml,l)) + dyp*(dxl*bxy(2,nl,mp,l) + amx*bxy(2,nn,mp,l) + dxp*bxy(
     32,np,mp,l))
      oz = amy*(dxl*bxy(3,nl,mm,l) + amx*bxy(3,nn,mm,l) + dxp*bxy(3,np,m
     1m,l)) + dyl*(dxl*bxy(3,nl,ml,l) + amx*bxy(3,nn,ml,l) + dxp*bxy(3,n
     2p,ml,l)) + dyp*(dxl*bxy(3,nl,mp,l) + amx*bxy(3,nn,mp,l) + dxp*bxy(
     33,np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,l) + (dx*dtt + omxt*omt)
      dy = part(2,j,l) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            omt = .5*omzt*dt
            dx = part(1,j,l)
            dy = dy + omt*part(3,j,l)
            if ((dy.lt.edgely).or.(dy.ge.edgery)) then
               dx = dx - omt*part(4,j,l)
               dy = part(2,j,l) + omt*part(3,j,l)
               part(4,j,l) = -part(4,j,l)
            endif
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            omt = .5*omzt*dt
            dy = part(2,j,l)
            dx = dx - omt*part(4,j,l)
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j,l) - omt*part(4,j,l)
               dy = dy + omt*part(3,j,l)
               part(3,j,l) = -part(3,j,l)
            endif
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            dy = dy + .5*omzt*dt*part(3,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH2CL(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,idi
     1mp,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover with correction.
c baseline scalar distributed version
c 121 flops/particle, 2 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bx(j,k,l) = x component of magnetic field at grid (j,kk)
c by(j,k,l,l) = y component of magnetic field at grid (j,kk)
c bz(j,k) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
c find electric field
      dx = amy*(amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + dyp*(amx*fx(nn,mp,l
     1) + dxp*fx(np,mp,l))
      dy = amy*(amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + dyp*(amx*fy(nn,mp,l
     1) + dxp*fy(np,mp,l))
c find magnetic field
      ox = amy*(amx*bx(nn,mm,l) + dxp*bx(np,mm,l)) + dyp*(amx*bx(nn,mp,l
     1) + dxp*bx(np,mp,l))
      oy = amy*(amx*by(nn,mm,l) + dxp*by(np,mm,l)) + dyp*(amx*by(nn,mp,l
     1) + dxp*by(np,mp,l))
      oz = amy*(amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + dyp*(amx*bz(nn,mp,l
     1) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,l) + (dx*dtt + omxt*omt)
      dy = part(2,j,l) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH2CL(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idimp,
     1npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover
c with correction.
c scalar version using guard cells, for distributed data
c 121 flops/particle, 2 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n
c part(2,n,l) = position y of particle n
c part(3,n,l) = velocity vx of particle n
c part(4,n,l) = velocity vy of particle n
c part(5,n,l) = velocity vz of particle n
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp,l) + amx*bxy(1,nn,mp,l)) + amy*(dxp*bxy(
     11,np,mm,l) + amx*bxy(1,nn,mm,l))
      oy = dyp*(dxp*bxy(2,np,mp,l) + amx*bxy(2,nn,mp,l)) + amy*(dxp*bxy(
     12,np,mm,l) + amx*bxy(2,nn,mm,l))
      oz = dyp*(dxp*bxy(3,np,mp,l) + amx*bxy(3,nn,mp,l)) + amy*(dxp*bxy(
     13,np,mm,l) + amx*bxy(3,nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,l) + (dx*dtt + omxt*omt)
      dy = part(2,j,l) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            omt = .5*omzt*dt
            dx = part(1,j,l)
            dy = dy + omt*part(3,j,l)
            if ((dy.lt.edgely).or.(dy.ge.edgery)) then
               dx = dx - omt*part(4,j,l)
               dy = part(2,j,l) + omt*part(3,j,l)
               part(4,j,l) = -part(4,j,l)
            endif
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            omt = .5*omzt*dt
            dy = part(2,j,l)
            dx = dx - omt*part(4,j,l)
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j,l) - omt*part(4,j,l)
               dy = dy + omt*part(3,j,l)
               part(3,j,l) = -part(3,j,l)
            endif
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            dy = dy + .5*omzt*dt*part(3,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field.  Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 209 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j+1,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
      dz = amy*(dxl*fxy(3,nl,mm,l) + amx*fxy(3,nn,mm,l) + dxp*fxy(3,np,m
     1m,l)) + dyl*(dxl*fxy(3,nl,ml,l) + amx*fxy(3,nn,ml,l) + dxp*fxy(3,n
     2p,ml,l)) + dyp*(dxl*fxy(3,nl,mp,l) + amx*fxy(3,nn,mp,l) + dxp*fxy(
     33,np,mp,l))
c find magnetic field
      ox = amy*(dxl*bxy(1,nl,mm,l) + amx*bxy(1,nn,mm,l) + dxp*bxy(1,np,m
     1m,l)) + dyl*(dxl*bxy(1,nl,ml,l) + amx*bxy(1,nn,ml,l) + dxp*bxy(1,n
     2p,ml,l)) + dyp*(dxl*bxy(1,nl,mp,l) + amx*bxy(1,nn,mp,l) + dxp*bxy(
     31,np,mp,l))
      oy = amy*(dxl*bxy(2,nl,mm,l) + amx*bxy(2,nn,mm,l) + dxp*bxy(2,np,m
     1m,l)) + dyl*(dxl*bxy(2,nl,ml,l) + amx*bxy(2,nn,ml,l) + dxp*bxy(2,n
     2p,ml,l)) + dyp*(dxl*bxy(2,nl,mp,l) + amx*bxy(2,nn,mp,l) + dxp*bxy(
     32,np,mp,l))
      oz = amy*(dxl*bxy(3,nl,mm,l) + amx*bxy(3,nn,mm,l) + dxp*bxy(3,np,m
     1m,l)) + dyl*(dxl*bxy(3,nl,ml,l) + amx*bxy(3,nn,ml,l) + dxp*bxy(3,n
     2p,ml,l)) + dyp*(dxl*bxy(3,nl,mp,l) + amx*bxy(3,nn,mp,l) + dxp*bxy(
     33,np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,id
     1imp,npmax,nblok,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 209 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j+1,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l) + .5
      mmn = part(2,j+1,l) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      mn = ml + nxv
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mp = mn + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mmn = mmn - mnoff
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
      dz = amy*(dxl*fxy(3,mn,l) + amx*fxy(3,mn+1,l) + dxp*fxy(3,mn+2,l))
     1 + dyl*(dxl*fxy(3,ml,l) + amx*fxy(3,ml+1,l) + dxp*fxy(3,ml+2,l)) +
     2 dyp*(dxl*fxy(3,mp,l) + amx*fxy(3,mp+1,l) + dxp*fxy(3,mp+2,l)) 
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
      dz = amy*(dxl*fxy(3,mn,l) + amx*fxy(3,mn+1,l) + dxp*fxy(3,mn+2,l))
     1 + dyl*(dxl*fxy(3,ml,l) + amx*fxy(3,ml+1,l) + dxp*fxy(3,ml+2,l)) +
     2 dyp*(dxl*fxy(3,mp,l) + amx*fxy(3,mp+1,l) + dxp*fxy(3,mp+2,l)) 
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
      acz = part(5,nop,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop,l) = dx
      part(4,nop,l) = dy
      part(5,nop,l) = dz
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,id
     1imp,npmax,nblok,nxv,nypmx,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 117 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
      dz = dyp*(dxp*fxy(3,np,mp,l) + amx*fxy(3,nn,mp,l)) + amy*(dxp*fxy(
     13,np,mm,l) + amx*fxy(3,nn,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp,l) + amx*bxy(1,nn,mp,l)) + amy*(dxp*bxy(
     11,np,mm,l) + amx*bxy(1,nn,mm,l))
      oy = dyp*(dxp*bxy(2,np,mp,l) + amx*bxy(2,nn,mp,l)) + amy*(dxp*bxy(
     12,np,mm,l) + amx*bxy(2,nn,mm,l))
      oz = dyp*(dxp*bxy(3,np,mp,l) + amx*bxy(3,nn,mp,l)) + amy*(dxp*bxy(
     13,np,mm,l) + amx*bxy(3,nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,i
     1dimp,npmax,nblok,nxv,nxyp,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 117 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l)
      mmn = part(2,j+1,l)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
      dz = dyp*(dxp*fxy(3,mp+1,l) + amx*fxy(3,mp,l)) + amy*(dxp*fxy(3,mm
     1+1,l) + amx*fxy(3,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxp*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyp*(dxp*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxp*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyp*(dxp*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxp*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
      acz = part(5,j,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j,l) = dx
      part(4,j,l) = dy
      part(5,j,l) = dz
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
      dz = dyn*(dxn*fxy(3,mp+1,l) + amx*fxy(3,mp,l)) + amy*(dxn*fxy(3,mm
     1+1,l) + amx*fxy(3,mm,l))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxn*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyn*(dxn*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxn*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyn*(dxn*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxn*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
      acz = part(5,nop,l) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop,l) = dx
      part(4,nop,l) = dy
      part(5,nop,l) = dz
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,idim
     1p,npmax,nblok,nxv,nypmx,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field in z direction.
c Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 141 flops/particle, 1 divide, 31 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
c find magnetic field
      oz = amy*(dxl*bz(nl,mm,l) + amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + d
     1yl*(dxl*bz(nl,ml,l) + amx*bz(nn,ml,l) + dxp*bz(np,ml,l)) + dyp*(dx
     2l*bz(nl,mp,l) + amx*bz(nn,mp,l) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nxyp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space with magnetic field in z direction.
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 141 flops/particle, 1 divide, 31 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l) + .5
      mmn = part(2,j+1,l) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      mn = ml + nxv
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mp = mn + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mmn = mmn - mnoff
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field 
      oz = amy*(dxl*bz(mn,l) + amx*bz(mn+1,l) + dxp*bz(mn+2,l)) + dyl*(d
     1xl*bz(ml,l) + amx*bz(ml+1,l) + dxp*bz(ml+2,l)) + dyp*(dxl*bz(mp,l)
     2 + amx*bz(mp+1,l) + dxp*bz(mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field 
      oz = amy*(dxl*bz(mn,l) + amx*bz(mn+1,l) + dxp*bz(mn+2,l)) + dyl*(d
     1xl*bz(ml,l) + amx*bz(ml+1,l) + dxp*bz(ml+2,l)) + dyp*(dxl*bz(mp,l)
     2 + amx*bz(mp+1,l) + dxp*bz(mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,nop,l) = dx
      part(4,nop,l) = dy
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field in z direction.
c Using the Boris Mover.
c scalar version using guard cells, for distributed data
c 62 flops/particle, 1 divide, 16 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c find magnetic field
      oz = dyp*(dxp*bz(np,mp,l) + amx*bz(nn,mp,l)) + amy*(dxp*bz(np,mm,l
     1) + amx*bz(nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,id
     1imp,npmax,nblok,nxv,nxyp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space with magnetic field in z direction.
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 62 flops/particle, 1 divide, 16 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l)
      mmn = part(2,j+1,l)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      oz = dyp*(dxp*bz(mp+1,l) + amx*bz(mp,l)) + amy*(dxp*bz(mm+1,l) + a
     1mx*bz(mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j,l) + dx
      acy = part(4,j,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dtc
      dy = part(2,j,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      oz = dyn*(dxn*bz(mp+1,l) + amx*bz(mp,l)) + amy*(dxn*bz(mm+1,l) + a
     1mx*bz(mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop,l) + dx
      acy = part(4,nop,l) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,nop,l) = dx
      part(4,nop,l) = dy
c new position
      dx = part(1,nop,l) + dx*dtc
      dy = part(2,nop,l) + dy*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PRETARD2(part,npp,dtc,nx,ny,idimp,npmax,nblok,ipbc)
c for 2-1/2d code, particle positions are retarded a half time-step
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - vx(t+dt/2)*dtc, y(t+dt) = y(t) - vy(t+dt/2)*dtc,
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c npp(l) = number of particles in partition l
c dtc = time interval between successive co-ordinate calculations
c nx/ny = system length in x/y direction
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
      do 20 l = 1, nblok
c retard position half a time-step for current deposit
      do 10 j = 1, npp(l)
      dx = part(1,j,l) - part(3,j,l)*dtc
      dy = part(2,j,l) - part(4,j,l)*dtc
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
      return
      end
