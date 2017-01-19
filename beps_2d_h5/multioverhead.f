c-----------------------------------------------------------------------
      program main
c This program calculates the variance of the first 10**9 integers
c using MPI with multiple dual processor Macintoshes.
      implicit none
c get definition of MPI constants
      include 'mpif.h'
      integer NINT
      parameter(NINT=1000000000)
      integer n, n1, n2, start1, start2
      double precision sum11, sum12, sum21, sum22, mean, mesq
      external adder
      integer i, ierror, nproc, idproc, root, idtask, nargs
      double precision partial(2), result(2)
      real time
      double precision dtime
      root = 0
      nargs = 4
c initialize the MPI environment
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,idproc,ierror)
c set arguments for each node, last node might have more work
      n = NINT/nproc
      start1 = n*idproc
      n1 = n + start1
      if (idproc.eq.(nproc-1)) n1 = NINT

c initialize the MP environment
      call MP_INIT(nproc)
c set arguments for each task
      n = n/2
      n2 = n1
      n1 = n1 - n
      start2 = start1 + n

      call PWTIMERA(-1,time,dtime)
      do 10 i = 1, 1000000
      call MP_TASKSTART(idtask,adder,nargs,n2,start2,sum21,sum22)
c     call adder(n1,start1,sum11,sum12)
      call MP_TASKWAIT(idtask)
   10 continue
      call PWTIMERA(1,time,dtime)

c terminate the MP environment
      call MP_END()

c add up all the partial sums from each node
      partial(1) = sum11 + sum21
      partial(2) = sum12 + sum22
c     call MPI_REDUCE(partial,result,2,MPI_DOUBLE_PRECISION,MPI_SUM,root
c    &,MPI_COMM_WORLD,ierror)

c node 0 has the results
      if (idproc.eq.root) then
         n = NINT
         mean = result(1) / n
         mesq = result(2) / n
         write (*,'(i15,1x,f20.3,1x,f21.3)') n, mean, mesq-(mean * mean)
      endif

      write (*,*)
      write (*,'("time in sec: ",f7.2)') time

c terminate MPI environment
c     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FINALIZE(ierror)

      stop
      end


      subroutine adder(n,start,s1,s2)
      integer n, start
      double precision s1, s2
      integer i
      double precision ls1, ls2
c     ls1 = 0
c     ls2 = 0
c     do 10 i = start, n-1
c     ls1 = ls1 + dble(i)
c     ls2 = ls2 + dble(i)*i
c  10 continue
c     s1 = ls1
c     s2 = ls2
      return
      end

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