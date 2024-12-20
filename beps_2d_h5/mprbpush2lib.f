c 2d PIC parallel multi-tasking library for pushing relativistic
c particles with magnetic field and depositing current
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: june 9, 2008
c-----------------------------------------------------------------------
      subroutine MPGRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,np
     1max,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRJPOST2
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRJPOST2,nargs,part(1,npo,1),cup(1,1,
     11,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRJPOST2(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,npma
     1x,nblok,nxv,nypmx,ipbc)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 3
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,n
     1pmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRJPOST2
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 3
      cup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSRJPOST2,nargs,part(1,npo,1),cup(1,1
     1,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSRJPOST2(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm
     1ax,nblok,nxv,nxyp,ipbc)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 3
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,nx,ny
     1,idimp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n27,npd,nblok,nmt+1), amxy(n27,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSRJOST2X
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, 3*nxvyp
      cup(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSRJOST2X,nargs,part(1,npo,1),cup(1,1
     1,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,ci,nx,ny,idimp,np
     2max,nblok,nxv,nxvyp,npd,n27,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSRJOST2X(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,i
     1dimp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 l = 1, nblok
      do 70 j = 1, 3*nxvyp
      cu(j,l) = cu(j,l) + cup(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,n
     1pmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRJPOST2L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRJPOST2L,nargs,part(1,npo,1),cup(1,1
     1,1,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRJPOST2L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm
     1ax,nblok,nxv,nypmx,ipbc)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 3
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,
     1npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRJPOST2L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 3
      cup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSRJPOST2L,nargs,part(1,npo,1),cup(1,
     11,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSRJPOST2L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,np
     1max,nblok,nxv,nxyp,ipbc)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 3
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,nx,n
     1y,idimp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n12,npd,nblok,nmt+1), amxy(n12,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSRJOST2XL
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, 3*nxvyp
      cup(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSRJOST2XL,nargs,part(1,npo,1),cup(1,
     11,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,ci,nx,ny,idimp,n
     2pmax,nblok,nxv,nxvyp,npd,n12,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSRJOST2XL(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,
     1idimp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 l = 1, nblok
      do 70 j = 1, 3*nxvyp
      cu(j,l) = cu(j,l) + cup(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,n
     1pmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRJPOST22
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRJPOST22,nargs,part(1,npo,1),cup(1,1
     1,1,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRJPOST22(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm
     1ax,nblok,nxv,nypmx,ipbc)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,
     1npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRJPOST22
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      cup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSRJPOST22,nargs,part(1,npo,1),cup(1,
     11,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSRJPOST22(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,np
     1max,nblok,nxv,nxyp,ipbc)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,nx,n
     1y,idimp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n18,npd,nblok,nmt+1), amxy(n18,npd,nblok,nmt+1)
      dimension cup(2*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSRJOST22X
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, 2*nxvyp
      cup(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSRJOST22X,nargs,part(1,npo,1),cup(1,
     11,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,ci,nx,ny,idimp,n
     2pmax,nblok,nxv,nxvyp,npd,n18,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSRJOST22X(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,
     1idimp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 l = 1, nblok
      do 70 j = 1, 2*nxvyp
      cu(j,l) = cu(j,l) + cup(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp,
     1npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRJPOST22L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRJPOST22L,nargs,part(1,npo,1),cup(1,
     11,1,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRJPOST22L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,np
     1max,nblok,nxv,nypmx,ipbc)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp
     1,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRJPOST22L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      cup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSRJPOST22L,nargs,part(1,npo,1),cup(1
     1,1,1,i),nps,noff,qm,dt,ci,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSRJPOST22L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,idimp,n
     1pmax,nblok,nxv,nxyp,ipbc)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,nx,
     1ny,idimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n8,npd,nblok,nmt+1), amxy(n8,npd,nblok,nmt+1)
      dimension cup(2*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSRJOST22XL
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, 2*nxvyp
      cup(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSRJOST22XL,nargs,part(1,npo,1),cup(1
     1,1,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,ci,nx,ny,idimp,
     2npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSRJOST22XL(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny
     1,idimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 l = 1, nblok
      do 70 j = 1, 2*nxvyp
      cu(j,l) = cu(j,l) + cup(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,idim
     1p,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRPUSH2
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRPUSH2,nargs,part(1,npo,1),fxy,nps,n
     1off,qbm,dt,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRPUSH2(part(1,npo,1),fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp,
     1npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRPUSH2
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRPUSH2,nargs,part(1,npo,1),fxy,nps,
     1noff,qbm,dt,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRPUSH2(part(1,npo,1),fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp
     1,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRPUSH2L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRPUSH2L,nargs,part(1,npo,1),fxy,nps,
     1noff,qbm,dt,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRPUSH2L(part(1,npo,1),fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp
     1,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,id
     1imp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRPUSH2L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRPUSH2L,nargs,part(1,npo,1),fxy,nps
     1,noff,qbm,dt,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRPUSH2L(part(1,npo,1),fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idim
     1p,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH2
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH2,nargs,part(1,npo,1),fxy,bxy,
     1nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ip
     2bc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH2(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH2
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH2,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ip
     2bc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH2(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH2L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH2L,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH2L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH2L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH2L,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH2L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH23
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH23,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH23(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH23
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH23,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH23(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH23L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH23L,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,
     2ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH23L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,e
     1k,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH23L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH23L,nargs,part(1,npo,1),fxy,b
     1xy,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,
     2ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH23L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH22
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH22,nargs,part(1,npo,1),fxy,bz,
     1nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ip
     2bc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH22(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH22
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH22,nargs,part(1,npo,1),fxy,bz
     1,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ip
     2bc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH22(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGRBPUSH22L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGRBPUSH22L,nargs,part(1,npo,1),fxy,bz
     1,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGRBPUSH22L(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSRBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSRBPUSH22L
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSRBPUSH22L,nargs,part(1,npo,1),fxy,b
     1z,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,i
     2pbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSRBPUSH22L(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPRPUSH2ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,nx,n
     1y,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ci, ek, ekp
      integer npp, nps
      integer idimp, npmax, nblok, nx, ny, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok), nps(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PRPUSH2ZF
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PRPUSH2ZF,nargs,part(1,npo,1),nps,dt,c
     1i,ekp(i),idimp,npmax,nblok,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PRPUSH2ZF(part(1,npo,1),npp,dt,ci,ek,idimp,npmax,nblok,nx,ny,
     1ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPRPUSH23ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,nx,
     1ny,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ci, ek, ekp
      integer npp, nps
      integer idimp, npmax, nblok, nx, ny, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok), nps(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PRPUSH23ZF
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PRPUSH23ZF,nargs,part(1,npo,1),nps,dt,
     1ci,ekp(i),idimp,npmax,nblok,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PRPUSH23ZF(part(1,npo,1),npp,dt,ci,ek,idimp,npmax,nblok,nx,ny
     1,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,idimp
     1,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current density deposit
c for relativistic particles
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, ci, cup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRCJPOST2
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRCJPOST2,nargs,part(1,npo,1),fxy,nps
     1,noff,cup(1,1,1,1,i),qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRCJPOST2(part(1,npo,1),fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,n
     1pmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGRCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,idim
     1p,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current density deposit
c for relativistic particles
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, ci, cup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRCJPOST2L
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGRCJPOST2L,nargs,part(1,npo,1),fxy,np
     1s,noff,cup(1,1,1,1,i),qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGRCJPOST2L(part(1,npo,1),fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,
     1npmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
