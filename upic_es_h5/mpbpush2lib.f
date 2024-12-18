c 2d PIC parallel multi-tasking library for pushing particles
c with magnetic field and depositing current
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: june 10, 2008
c-----------------------------------------------------------------------
      subroutine MPJDOST2(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idimp,n
     1pmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cux(nxv,nypmx,nblok)
      dimension cuy(nxv,nypmx,nblok), cuz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(nxv,nypmx,3,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PJDOST2
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
      do 40 m = 1, 3
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      cup(j,k,m,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PJDOST2,nargs,part(1,npo,1),cup(1,1,1,
     11,i),cup(1,1,2,1,i),cup(1,1,3,1,i),nps,noff,qm,dt,nx,idimp,npmax,n
     2blok,nxv,nypmx)
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
      call PJDOST2(part(1,npo,1),cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npm
     1ax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 120 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 110 l = 1, nblok
      do 100 k = 1, nypmx
      do 90 j = 1, nx
      cux(j,k,l) = cux(j,k,l) + cup(j,k,1,l,i)
      cuy(j,k,l) = cuy(j,k,l) + cup(j,k,2,l,i)
      cuz(j,k,l) = cuz(j,k,l) + cup(j,k,3,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npmax,
     1nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGJPOST2
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
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGJPOST2,nargs,part(1,npo,1),cup(1,1,1
     1,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGJPOST2(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nb
     1lok,nxv,nypmx,ipbc)
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
      subroutine MPGSJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npmax
     1,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSJPOST2
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
      call MP_TASKSTART(idtask(i),PGSJPOST2,nargs,part(1,npo,1),cup(1,1,
     11,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSJPOST2(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n
     1blok,nxv,nxyp,ipbc)
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
      subroutine MPSJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,idimp,n
     1pmax,nblok,nxv,nxvyp,npd,n27,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, nblok, nxv, nxvyp, npd, n27
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n27,npd,nblok,nmt+1), amxy(n27,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PSJOST2X
      data nargs /16/
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
      call MP_TASKSTART(idtask(i),PSJOST2X,nargs,part(1,npo,1),cup(1,1,i
     1),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,idimp,npmax,nblo
     2k,nxv,nxvyp,npd,n27)
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
      call PSJOST2X(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,idimp,npm
     1ax,nblok,nxv,nxvyp,npd,n27)
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
      subroutine MPGSJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny,idi
     1mp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n27,npd,nblok,nmt+1), amxy(n27,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSJOST2X
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSJOST2X,nargs,part(1,npo,1),cup(1,1,
     1i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,ny,idimp,npmax,
     2nblok,nxv,nxvyp,npd,n27,ipbc)
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
      call PGSJOST2X(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp
     1,npmax,nblok,nxv,nxvyp,npd,n27,ipbc)
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
      subroutine MPJDOST2L(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idimp,
     1npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cux(nxv,nypmx,nblok)
      dimension cuy(nxv,nypmx,nblok), cuz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(nxv,nypmx,3,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PJDOST2L
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
      do 40 m = 1, 3
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      cup(j,k,m,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PJDOST2L,nargs,part(1,npo,1),cup(1,1,1
     1,1,i),cup(1,1,2,1,i),cup(1,1,3,1,i),nps,noff,qm,dt,nx,idimp,npmax,
     2nblok,nxv,nypmx)
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
      call PJDOST2L(part(1,npo,1),cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,np
     1max,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 120 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 110 l = 1, nblok
      do 100 k = 1, nypmx
      do 90 j = 1, nx
      cux(j,k,l) = cux(j,k,l) + cup(j,k,1,l,i)
      cuy(j,k,l) = cuy(j,k,l) + cup(j,k,2,l,i)
      cuz(j,k,l) = cuz(j,k,l) + cup(j,k,3,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npmax
     1,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGJPOST2L
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
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGJPOST2L,nargs,part(1,npo,1),cup(1,1,
     11,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGJPOST2L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n
     1blok,nxv,nypmx,ipbc)
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
      subroutine MPGSJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npma
     1x,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSJPOST2L
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
      call MP_TASKSTART(idtask(i),PGSJPOST2L,nargs,part(1,npo,1),cup(1,1
     1,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSJPOST2L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,
     1nblok,nxv,nxyp,ipbc)
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
      subroutine MPSJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,idimp,
     1npmax,nblok,nxv,nxvyp,npd,n12,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, nblok, nxv, nxvyp, npd, n12
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n12,npd,nblok,nmt+1), amxy(n12,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PSJOST2XL
      data nargs /16/
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
      call MP_TASKSTART(idtask(i),PSJOST2XL,nargs,part(1,npo,1),cup(1,1,
     1i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,idimp,npmax,nbl
     2ok,nxv,nxvyp,npd,n12)
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
      call PSJOST2XL(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,idimp,np
     1max,nblok,nxv,nxvyp,npd,n12)
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
      subroutine MPGSJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny,id
     1imp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(3*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n12,npd,nblok,nmt+1), amxy(n12,npd,nblok,nmt+1)
      dimension cup(3*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSJOST2XL
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSJOST2XL,nargs,part(1,npo,1),cup(1,1
     1,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,ny,idimp,npmax
     2,nblok,nxv,nxvyp,npd,n12,ipbc)
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
      call PGSJOST2XL(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,ny,idim
     1p,npmax,nblok,nxv,nxvyp,npd,n12,ipbc)
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
      subroutine MPGJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npmax
     1,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGJPOST22
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
      call MP_TASKSTART(idtask(i),PGJPOST22,nargs,part(1,npo,1),cup(1,1,
     11,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGJPOST22(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n
     1blok,nxv,nypmx,ipbc)
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
      subroutine MPGSJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npma
     1x,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSJPOST22
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
      call MP_TASKSTART(idtask(i),PGSJPOST22,nargs,part(1,npo,1),cup(1,1
     1,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSJPOST22(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,
     1nblok,nxv,nxyp,ipbc)
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
      subroutine MPGSJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny,id
     1imp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n18,npd,nblok,nmt+1), amxy(n18,npd,nblok,nmt+1)
      dimension cup(2*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSJOST22X
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSJOST22X,nargs,part(1,npo,1),cup(1,1
     1,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,ny,idimp,npmax
     2,nblok,nxv,nxvyp,npd,n18,ipbc)
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
      call PGSJOST22X(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,ny,idim
     1p,npmax,nblok,nxv,nxvyp,npd,n18,ipbc)
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
      subroutine MPGJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npma
     1x,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGJPOST22L
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
      call MP_TASKSTART(idtask(i),PGJPOST22L,nargs,part(1,npo,1),cup(1,1
     1,1,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGJPOST22L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax,
     1nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npm
     1ax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSJPOST22L
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
      call MP_TASKSTART(idtask(i),PGSJPOST22L,nargs,part(1,npo,1),cup(1,
     11,1,i),nps,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSJPOST22L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,idimp,npmax
     1,nblok,nxv,nxyp,ipbc)
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
      subroutine MPGSJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny,i
     1dimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxy, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), cu(2*nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(n8,npd,nblok,nmt+1), amxy(n8,npd,nblok,nmt+1)
      dimension cup(2*nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSJOST22XL
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSJOST22XL,nargs,part(1,npo,1),cup(1,
     11,i),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,dt,nx,ny,idimp,npma
     2x,nblok,nxv,nxvyp,npd,n8,ipbc)
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
      call PGSJOST22XL(part(1,npo,1),cu,npp,noff,nn,amxy,qm,dt,nx,ny,idi
     1mp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
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
      subroutine MPBPUSH2(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek,nx,
     1idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PBPUSH2
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PBPUSH2,nargs,part(1,npo,1),fx,fy,bx,b
     1y,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
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
      call PBPUSH2(part(1,npo,1),fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,id
     1imp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH2
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH2,nargs,part(1,npo,1),fxy,bxy,n
     1ps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH2(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,i
     1dimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH2
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH2,nargs,part(1,npo,1),fxy,bxy,
     1nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSBPUSH2(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      subroutine MPBPUSH2L(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek,nx
     1,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PBPUSH2L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PBPUSH2L,nargs,part(1,npo,1),fx,fy,bx,
     1by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
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
      call PBPUSH2L(part(1,npo,1),fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,i
     1dimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH2L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH2L,nargs,part(1,npo,1),fxy,bxy,
     1nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH2L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH2L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH2L,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSBPUSH2L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      subroutine MPBPUSH2CQ(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek,n
     1x,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PBPUSH2CQ
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PBPUSH2CQ,nargs,part(1,npo,1),fx,fy,bx
     1,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
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
      call PBPUSH2CQ(part(1,npo,1),fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,
     1idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGBPUSH2CQ(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,ny,i
     1dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH2CQ
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
      call MP_TASKSTART(idtask(i),PGBPUSH2CQ,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH2CQ(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPBPUSH2CL(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek,n
     1x,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension bx(nxv,nypmx,nblok), by(nxv,nypmx,nblok)
      dimension bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PBPUSH2CL
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PBPUSH2CL,nargs,part(1,npo,1),fx,fy,bx
     1,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
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
      call PBPUSH2CL(part(1,npo,1),fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,
     1idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGBPUSH2CL(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,ny,i
     1dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH2CL
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
      call MP_TASKSTART(idtask(i),PGBPUSH2CL,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH2CL(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idi
     1mp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH23
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH23,nargs,part(1,npo,1),fxy,bxy,
     1nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH23(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH23
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH23,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSBPUSH23(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      subroutine MPGBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH23L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH23L,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc
     2)
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
      call PGBPUSH23L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx
     1,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH23L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH23L,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc
     2)
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
      call PGSBPUSH23L(part(1,npo,1),fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      subroutine MPGBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH22
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH22,nargs,part(1,npo,1),fxy,bz,n
     1ps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH22(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,i
     1dimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH22
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH22,nargs,part(1,npo,1),fxy,bz,
     1nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSBPUSH22(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      subroutine MPGBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,nx,n
     1y,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGBPUSH22L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGBPUSH22L,nargs,part(1,npo,1),fxy,bz,
     1nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      call PGBPUSH22L(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1idimp,npmax,nblok,nxv,nypmx,ipbc)
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
      subroutine MPGSBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,nx,
     1ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSBPUSH22L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSBPUSH22L,nargs,part(1,npo,1),fxy,bz
     1,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
      call PGSBPUSH22L(part(1,npo,1),fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny
     1,idimp,npmax,nblok,nxv,nxyp,ipbc)
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
