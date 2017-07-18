      subroutine PGSPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,np
     1max,nblok,nxv,nxyp,ipbc)
c modified to store vdotE
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = second virtual dimension of field array, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
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
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = .5*(.5 - dxp)**2
      dxp = .5*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c find acceleration
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1+ dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) + 
     2dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l))
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1+ dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) + 
     2dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
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
c find acceleration
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1+ dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) + 
     2dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l))
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1+ dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) + 
     2dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
c new velocity
      dx = part(3,nop,l) + qtm*dx
      dy = part(4,nop,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop,l))**2 + (dy + part(4,nop,l))**2
      part(3,nop,l) = dx
      part(4,nop,l) = dy
c new position
      dx = part(1,nop,l) + dx*dt
      dy = part(2,nop,l) + dy*dt
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
      ek = ek + .125*sum1
      return
      end subroutine PGSPUSH2_jf
      
c-----------------------------------------------------------------------
      subroutine PGPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm
     1ax,nblok,nxv,nypmx,ipbc,nvexloc,nveyloc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, for distributed data
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
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
c find acceleration
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))

c store old velocities
      old_vx = part(3,j,l)
      old_vy = part(4,j,l)
c store field
      old_dx = dx
      old_dy = dy
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
c store vdote
      part(nvexloc,j,l) = 0.5*(part(3,j,l)+dx)*old_dx
      part(nveyloc,j,l) = 0.5*(part(4,j,l)+dy)*old_dy
      part(3,j,l) = dx
      part(4,j,l) = dy

c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
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
      ek = ek + .125*sum1
      return
      end subroutine PGPUSH2_jf
      
c-----------------------------------------------------------------------
      subroutine GETVPHI(part,fxy,phi,npp,noff,qbm,dt,ek,nx,ny,idimp,npm
     1ax,nblok,nxv,nypmx,ipbc,nvexloc,nveyloc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, for distributed data
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension phi(nxv,nypmx,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
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
c find acceleration
      dp = amy*(dxl*phi(nl,mm,l) + amx*phi(nn,mm,l) + dxp*phi(np,mm,l)) 
     1+ dyl*(dxl*phi(nl,ml,l) + amx*phi(nn,ml,l) + dxp*phi(np,ml,l)) + d
     2yp*(dxl*phi(nl,mp,l) + amx*phi(nn,mp,l) + dxp*phi(np,mp,l))
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))

c store old velocities
c      old_vx = part(3,j,l)
c      old_vy = part(4,j,l)
c store field
c      old_dx = dx
c      old_dy = dy
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
c store vdote
      part(nvexloc,j,l) = 0.5*(part(3,j,l)+dx)
      part(nveyloc,j,l) = 0.5*(part(4,j,l)+dy)
      part(3,j,l) = dx
      part(4,j,l) = dy

c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
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
      ek = ek + .125*sum1
      return
      end subroutine GETVPHI
      
      
c-----------------------------------------------------------------------
      subroutine PGPUSH2_jf_sum(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,&
     1npmax,nblok,nxv,nypmx,ipbc,nvexloc,nveyloc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, for distributed data
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
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
c find acceleration
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))

c store old velocities
      old_vx = part(3,j,l)
      old_vy = part(4,j,l)
c store field
      old_dx = dx
      old_dy = dy
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
      
c store vdotE
      part(nvexloc,j,l) = part(nvexloc,j,l) + 0.5*(part(3,j,l)+old_vx)*o
     1ld_dx
      part(nveyloc,j,l) = part(nveyloc,j,l) + 0.5*(part(4,j,l)+old_vy)*o
     1ld_dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
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
      ek = ek + .125*sum1
      return
      end subroutine PGPUSH2_jf_sum
      
c modified deposit for depositing vdotE_part data
c-----------------------------------------------------------------------
      subroutine PGSJPOST22_jf(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,
     1nblok,nxv,nxyp,ipbc,nvexloc,nveyloc)
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
      vx = part(nvexloc,j-1,l)
      vy = part(nveyloc,j-1,l)
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
      vx = part(nvexloc,nop,l)
      vy = part(nveyloc,nop,l)
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
      end subroutine PGSJPOST22_jf
c-----------------------------------------------------------------------
      subroutine PGCJEJPOST2_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp,
     1npmax,nblok,nxv,nypmx)
c Modification to PGCJPOST2 to deposit j.E too, into jE array
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation.
c scalar version using guard cells, for distributed data
c Flop count is now incorrect
c 96 flops/particle, 40 loads, 18 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, cu, jE, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension jE(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy
      real t1, t2
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      t1 = vx*dx
      t2 = vy*dx
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + t1!vx*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + t2!vy*dx
      jE(1,nl,mm,l) = jE(1,nl,mm,l) + t1*fxy(1,nl,mm,l)
      jE(2,nl,mm,l) = jE(2,nl,mm,l) + t2*fxy(2,nl,mm,l)
c      jE(1,nl,mm,l) = jE(1,nl,mm,l) + vx*dx*fxy(1,nl,mm,l)
c      jE(2,nl,mm,l) = jE(2,nl,mm,l) + vy*dx*fxy(2,nl,mm,l)
      dx = dxl*dyl
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + t1!vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + t2!vy*dy
      jE(1,nn,mm,l) = jE(1,nn,mm,l) + t1*fxy(1,nn,mm,l)
      jE(2,nn,mm,l) = jE(2,nn,mm,l) + t2*fxy(2,nn,mm,l)
c      jE(1,nn,mm,l) = jE(1,nn,mm,l) + vx*dy*fxy(1,nn,mm,l)
c      jE(2,nn,mm,l) = jE(2,nn,mm,l) + vy*dy*fxy(2,nn,mm,l)
      dy = amx*dyl
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,mm,l) = cu(1,np,mm,l) + t1!vx*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + t2!vy*dz
      jE(1,np,mm,l) = jE(1,np,mm,l) + t1*fxy(1,np,mm,l)
      jE(2,np,mm,l) = jE(2,np,mm,l) + t2*fxy(2,np,mm,l)
c      jE(1,np,mm,l) = jE(1,np,mm,l) + vx*dz*fxy(1,np,mm,l)
c      jE(2,np,mm,l) = jE(2,np,mm,l) + vy*dz*fxy(2,np,mm,l)
      dz = dxp*dyl
      t1 = vx*dx
      t2 = dy*dx
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + t1!vx*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + t2!vy*dx
      jE(1,nl,ml,l) = jE(1,nl,ml,l) + t1*fxy(1,nl,ml,l)
      jE(2,nl,ml,l) = jE(2,nl,ml,l) + t2*fxy(2,nl,ml,l)
c      jE(1,nl,ml,l) = jE(1,nl,ml,l) + vx*dx*fxy(1,nl,ml,l)
c      jE(2,nl,ml,l) = jE(2,nl,ml,l) + vy*dx*fxy(2,nl,ml,l)
      dx = dxl*dyp
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + t1!vx*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + t2!vy*dy
      jE(1,nn,ml,l) = jE(1,nn,ml,l) + t1*fxy(1,nn,ml,l)
      jE(2,nn,ml,l) = jE(2,nn,ml,l) + t2*fxy(2,nn,ml,l)
c      jE(1,nn,ml,l) = jE(1,nn,ml,l) + vx*dy*fxy(1,nn,ml,l)
c      jE(2,nn,ml,l) = jE(2,nn,ml,l) + vy*dy*fxy(2,nn,ml,l)
      dy = amx*dyp
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,ml,l) = cu(1,np,ml,l) + t1!vx*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + t2!vy*dz
      jE(1,np,ml,l) = jE(1,np,ml,l) + t1*fxy(1,np,ml,l)
      jE(2,np,ml,l) = jE(2,np,ml,l) + t2*fxy(2,np,ml,l)
c      jE(1,np,ml,l) = jE(1,np,ml,l) + vx*dz*fxy(1,np,ml,l)
c      jE(2,np,ml,l) = jE(2,np,ml,l) + vy*dz*fxy(2,np,ml,l)
      dz = dxp*dyp
      t1 = vx*dx
      t2 = vy*dx
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + t1!vx*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + t2!vy*dx
      jE(1,nl,mp,l) = jE(1,nl,mp,l) + t1*fxy(1,nl,mp,l)
      jE(2,nl,mp,l) = jE(2,nl,mp,l) + t2*fxy(2,nl,mp,l)
c      jE(1,nl,mp,l) = jE(1,nl,mp,l) + vx*dx*fxy(1,nl,mp,l)
c      jE(2,nl,mp,l) = jE(2,nl,mp,l) + vy*dx*fxy(2,nl,mp,l)
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + t1!vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + t2!vy*dy
      jE(1,nn,mp,l) = jE(1,nn,mp,l) + t1*fxy(1,nn,mp,l)
      jE(2,nn,mp,l) = jE(2,nn,mp,l) + t2*fxy(2,nn,mp,l)
c      jE(1,nn,mp,l) = jE(1,nn,mp,l) + vx*dy*fxy(1,nn,mp,l)
c      jE(2,nn,mp,l) = jE(2,nn,mp,l) + vy*dy*fxy(2,nn,mp,l)
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,mp,l) = cu(1,np,mp,l) + t1!vx*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + t2!vy*dz
      jE(1,np,mp,l) = jE(1,np,mp,l) + t1*fxy(1,np,mp,l)
      jE(2,np,mp,l) = jE(2,np,mp,l) + t2*fxy(2,np,mp,l)
c      jE(1,np,mp,l) = jE(1,np,mp,l) + vx*dz*fxy(1,np,mp,l)
c      jE(2,np,mp,l) = jE(2,np,mp,l) + vy*dz*fxy(2,np,mp,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGCJEJPOST2L_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp
     1,npmax,nblok,nxv,nypmx)
c for 2d code, this subroutine calculates particle current density
c using first-order spline interpolation.
c scalar version using guard cells, for distributed data
c 52 flops/particle, 20 loads, 0 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1.m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j,k,l) = ith component of current density
c at grid point j,kk for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp,noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, cu, jE, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension jE(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, vx, vy, t1, t2
      qtmh = 0.5*qbm*dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
      dx = dxp*dyp
      dy = amx*dyp
      t1 = vx*dx
      t2 = vy*dx
      cu(1,np,mp,l) = cu(1,np,mp,l) + t1!vx*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + t2!vy*dx
      jE(1,np,mp,l) = jE(1,np,mp,l) + t1*fxy(1,np,mp,l)! vx*dx
      jE(2,np,mp,l) = jE(2,np,mp,l) + t2*fxy(2,np,mp,l)! vy*dx
      dx = dxp*amy
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + t1!vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + t2!vy*dy
      jE(1,nn,mp,l) = jE(1,nn,mp,l) + t1*fxy(1,nn,mp,l)!  vx*dy
      jE(2,nn,mp,l) = jE(2,nn,mp,l) + t2*fxy(2,nn,mp,l)!  vy*dy
      dy = amx*amy
      t1 = vx*dx
      t2 = vy*dx
      cu(1,np,mm,l) = cu(1,np,mm,l) + t1!vx*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + t2!vy*dx
      jE(1,np,mm,l) = jE(1,np,mm,l) + t1*fxy(1,np,mm,l)!  vx*dx
      jE(2,np,mm,l) = jE(2,np,mm,l) + t2*fxy(2,np,mm,l)!  vy*dx
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + t1!vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + t2!vy*dy
      jE(1,nn,mm,l) = jE(1,nn,mm,l) + t1*fxy(1,nn,mm,l)!  vx*dy
      jE(2,nn,mm,l) = jE(2,nn,mm,l) + t2*fxy(2,nn,mm,l)!  vy*dy
   10 continue
   20 continue
      return
      end
      
c-----------------------------------------------------------------------
      subroutine PGCJEKEJPOST2_jf(part,fxy,npp,noff,cu,jE,kE,qm,qbm,dt,i
     1dimp,npmax,nblok,nxv,nypmx)
c Modification to PGCJPOST2 to deposit j.E too, into jE array
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation.
c scalar version using guard cells, for distributed data
c Flop count is now incorrect
c 96 flops/particle, 40 loads, 18 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, cu, jE, qm, qbm, dt, kE
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension jE(2,nxv,nypmx,nblok)
      dimension kE(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, ene
      real t1, t2
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
      ene = 0.5*(vx*vx + vy*vy)
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      t1 = vx*dx
      t2 = vy*dx
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + t1!vx*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + t2!vy*dx
      jE(1,nl,mm,l) = jE(1,nl,mm,l) + t1*fxy(1,nl,mm,l)
      jE(2,nl,mm,l) = jE(2,nl,mm,l) + t2*fxy(2,nl,mm,l)
      kE(nl,mm,l) = kE(nl,mm,l) + ene*dx
c      jE(1,nl,mm,l) = jE(1,nl,mm,l) + vx*dx*fxy(1,nl,mm,l)
c      jE(2,nl,mm,l) = jE(2,nl,mm,l) + vy*dx*fxy(2,nl,mm,l)
      dx = dxl*dyl
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + t1!vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + t2!vy*dy
      jE(1,nn,mm,l) = jE(1,nn,mm,l) + t1*fxy(1,nn,mm,l)
      jE(2,nn,mm,l) = jE(2,nn,mm,l) + t2*fxy(2,nn,mm,l)
      kE(nn,mm,l) = kE(nn,mm,l) + ene*dy
c      jE(1,nn,mm,l) = jE(1,nn,mm,l) + vx*dy*fxy(1,nn,mm,l)
c      jE(2,nn,mm,l) = jE(2,nn,mm,l) + vy*dy*fxy(2,nn,mm,l)
      dy = amx*dyl
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,mm,l) = cu(1,np,mm,l) + t1!vx*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + t2!vy*dz
      jE(1,np,mm,l) = jE(1,np,mm,l) + t1*fxy(1,np,mm,l)
      jE(2,np,mm,l) = jE(2,np,mm,l) + t2*fxy(2,np,mm,l)
      kE(np,mm,l) = kE(np,mm,l) + ene*dz
c      jE(1,np,mm,l) = jE(1,np,mm,l) + vx*dz*fxy(1,np,mm,l)
c      jE(2,np,mm,l) = jE(2,np,mm,l) + vy*dz*fxy(2,np,mm,l)
      dz = dxp*dyl
      t1 = vx*dx
      t2 = vy*dx
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + t1!vx*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + t2!vy*dx
      jE(1,nl,ml,l) = jE(1,nl,ml,l) + t1*fxy(1,nl,ml,l)
      jE(2,nl,ml,l) = jE(2,nl,ml,l) + t2*fxy(2,nl,ml,l)
      kE(nl,ml,l) = kE(nl,ml,l) + ene*dx
c      jE(1,nl,ml,l) = jE(1,nl,ml,l) + vx*dx*fxy(1,nl,ml,l)
c      jE(2,nl,ml,l) = jE(2,nl,ml,l) + vy*dx*fxy(2,nl,ml,l)
      dx = dxl*dyp
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + t1!vx*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + t2!vy*dy
      jE(1,nn,ml,l) = jE(1,nn,ml,l) + t1*fxy(1,nn,ml,l)
      jE(2,nn,ml,l) = jE(2,nn,ml,l) + t2*fxy(2,nn,ml,l)
      kE(nn,ml,l) = kE(nn,ml,l) + ene*dy
c      jE(1,nn,ml,l) = jE(1,nn,ml,l) + vx*dy*fxy(1,nn,ml,l)
c      jE(2,nn,ml,l) = jE(2,nn,ml,l) + vy*dy*fxy(2,nn,ml,l)
      dy = amx*dyp
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,ml,l) = cu(1,np,ml,l) + t1!vx*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + t2!vy*dz
      jE(1,np,ml,l) = jE(1,np,ml,l) + t1*fxy(1,np,ml,l)
      jE(2,np,ml,l) = jE(2,np,ml,l) + t2*fxy(2,np,ml,l)
      kE(np,ml,l) = kE(np,ml,l) + ene*dz
c      jE(1,np,ml,l) = jE(1,np,ml,l) + vx*dz*fxy(1,np,ml,l)
c      jE(2,np,ml,l) = jE(2,np,ml,l) + vy*dz*fxy(2,np,ml,l)
      dz = dxp*dyp
      t1 = vx*dx
      t2 = vy*dx
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + t1!vx*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + t2!vy*dx
      jE(1,nl,mp,l) = jE(1,nl,mp,l) + t1*fxy(1,nl,mp,l)
      jE(2,nl,mp,l) = jE(2,nl,mp,l) + t2*fxy(2,nl,mp,l)
      kE(nl,mp,l) = kE(nl,mp,l) + ene*dx
c      jE(1,nl,mp,l) = jE(1,nl,mp,l) + vx*dx*fxy(1,nl,mp,l)
c      jE(2,nl,mp,l) = jE(2,nl,mp,l) + vy*dx*fxy(2,nl,mp,l)
      t1 = vx*dy
      t2 = vy*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + t1!vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + t2!vy*dy
      jE(1,nn,mp,l) = jE(1,nn,mp,l) + t1*fxy(1,nn,mp,l)
      jE(2,nn,mp,l) = jE(2,nn,mp,l) + t2*fxy(2,nn,mp,l)
      kE(nn,mp,l) = kE(nn,mp,l) + ene*dy
c      jE(1,nn,mp,l) = jE(1,nn,mp,l) + vx*dy*fxy(1,nn,mp,l)
c      jE(2,nn,mp,l) = jE(2,nn,mp,l) + vy*dy*fxy(2,nn,mp,l)
      t1 = vx*dz
      t2 = vy*dz
      cu(1,np,mp,l) = cu(1,np,mp,l) + t1!vx*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + t2!vy*dz
      jE(1,np,mp,l) = jE(1,np,mp,l) + t1*fxy(1,np,mp,l)
      jE(2,np,mp,l) = jE(2,np,mp,l) + t2*fxy(2,np,mp,l)
      kE(np,mp,l) = kE(np,mp,l) + ene*dz
c      jE(1,np,mp,l) = jE(1,np,mp,l) + vx*dz*fxy(1,np,mp,l)
c      jE(2,np,mp,l) = jE(2,np,mp,l) + vy*dz*fxy(2,np,mp,l)
   10 continue
   20 continue
      return
      end
      
      subroutine PGCJEPOST2_onlytrackjf(part,fxy,npp,noff,jE,qm,qbm,dt,
     1idimp,npmax,nblok,nxv,nypmx,kivx_loc)
c Modification to PGCJEPOST2_J to deposit j.E too, into jE array
c This one only deposits jE for particles whose vx initial is < 10000.
c all others are those not being tracked
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation.
c scalar version using guard cells, for distributed data
c Flop count is now incorrect
c 96 flops/particle, 40 loads, 18 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx,kivx_loc
      real part, fxy, cu, jE, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension jE(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy
      real t1, t2
      qtmh = 0.5*qbm*dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
      if (part(kivx_loc,j,1) < 10000.) then
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
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      t1 = vx*dx
      t2 = vy*dx
      jE(1,nl,mm,l) = jE(1,nl,mm,l) + t1*fxy(1,nl,mm,l)
      jE(2,nl,mm,l) = jE(2,nl,mm,l) + t2*fxy(2,nl,mm,l)
      dx = dxl*dyl
      t1 = vx*dy
      t2 = vy*dy
      jE(1,nn,mm,l) = jE(1,nn,mm,l) + t1*fxy(1,nn,mm,l)
      jE(2,nn,mm,l) = jE(2,nn,mm,l) + t2*fxy(2,nn,mm,l)
      dy = amx*dyl
      t1 = vx*dz
      t2 = vy*dz
      jE(1,np,mm,l) = jE(1,np,mm,l) + t1*fxy(1,np,mm,l)
      jE(2,np,mm,l) = jE(2,np,mm,l) + t2*fxy(2,np,mm,l)
      dz = dxp*dyl
      t1 = vx*dx
      t2 = dy*dx
      jE(1,nl,ml,l) = jE(1,nl,ml,l) + t1*fxy(1,nl,ml,l)
      jE(2,nl,ml,l) = jE(2,nl,ml,l) + t2*fxy(2,nl,ml,l)
      dx = dxl*dyp
      t1 = vx*dy
      t2 = vy*dy
      jE(1,nn,ml,l) = jE(1,nn,ml,l) + t1*fxy(1,nn,ml,l)
      jE(2,nn,ml,l) = jE(2,nn,ml,l) + t2*fxy(2,nn,ml,l)
      dy = amx*dyp
      t1 = vx*dz
      t2 = vy*dz
      jE(1,np,ml,l) = jE(1,np,ml,l) + t1*fxy(1,np,ml,l)
      jE(2,np,ml,l) = jE(2,np,ml,l) + t2*fxy(2,np,ml,l)
      dz = dxp*dyp
      t1 = vx*dx
      t2 = vy*dx
      jE(1,nl,mp,l) = jE(1,nl,mp,l) + t1*fxy(1,nl,mp,l)
      jE(2,nl,mp,l) = jE(2,nl,mp,l) + t2*fxy(2,nl,mp,l)
      t1 = vx*dy
      t2 = vy*dy
      jE(1,nn,mp,l) = jE(1,nn,mp,l) + t1*fxy(1,nn,mp,l)
      jE(2,nn,mp,l) = jE(2,nn,mp,l) + t2*fxy(2,nn,mp,l)
      t1 = vx*dz
      t2 = vy*dz
      jE(1,np,mp,l) = jE(1,np,mp,l) + t1*fxy(1,np,mp,l)
      jE(2,np,mp,l) = jE(2,np,mp,l) + t2*fxy(2,np,mp,l)
      endif
   10 continue
   20 continue
      return
      end