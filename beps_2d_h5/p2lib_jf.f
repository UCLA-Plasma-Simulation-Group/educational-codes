c			use par_track_jf

! These functions have been modified by JF to do particle tracking
      subroutine PMOVEH2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idp
     1s,ntmax)
c this subroutine determines list of particles which are leaving this
c processor
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c ihole = location of holes left in particle arrays
c jss(1,l) = number of particles leaving, for particle partition l
c jss(2,l) = (0,1) = (no,yes) ihole overflowed, for particle partition l
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
      implicit none
      real part, edges
      integer npp, ihole, jss
      integer idimp, npmax, nblok, idps, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
c      type (t_track_set) :: tracks
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer j, l
      real yt
c buffer outgoing particles
      do 20 l = 1, nblok
      jss(1,l) = 0
      jss(2,l) = 0
c find particles out of bounds
      do 10 j = 1, npp(l)
      yt = part(iy,j,l)
      if ((yt.ge.edges(2,l)).or.(yt.lt.edges(1,l))) then
         if (jss(1,l).lt.ntmax) then
            jss(1,l) = jss(1,l) + 1
            ihole(jss(1,l),l) = j
         else
            jss(2,l) = 1
            go to 20
         endif
      endif
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------