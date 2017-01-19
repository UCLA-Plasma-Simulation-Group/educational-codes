#Makefile for 2D parallel PIC codes in new_pbeps2.source

GOBJS = nullpgks2.o nullpgks1.o

# Makefile Absoft compiler with MacOS X

#MPIFC = f90
#FC90 = f90
#FC77 = f77
#CC = gcc

#OPTS90 = -O3 -N113 -cpu:g5
#OPTS77 = -O -N113 -cpu:g5
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s -N11
#LOPTS = -plainappl
LEGACY =

#MPIOBJS = MacMPIcf.o MacMPI_S.o
#MPOBJS = MacMPcf.o MacMP.o
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libmcX.o

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile Nag compiler with MacOS X

#MPIFC = f95
#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS = -framework carbon 
#LEGACY = -dusty

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS =

# Makefile IBM compiler with MacOS X

#MPIFC = xlf90
#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -C -qautodbl=dbl4 -qarch=g5
#OPTS77 = -C -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = MacMPIxlf.o MacMPI_S.o
#MPOBJS = MacMPxlf.o MacMP.o

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
#		-L/Users/uclapic/hdflib -lz -ljpeg -ldf -lmfhdf 

# Makefile mpif90 compiler with LAM MPI and MacOS X

MPIFC = /usr/local/apps/openmpi/current/bin/mpif90
FC90 = /usr/local/apps/openmpi/current/bin/mpif90
FC77 = /usr/local/apps/openmpi/current/bin/mpif77
CC = /usr/local/apps/openmpi/current/bin/mpicc

OPTS90 = 
OPTS77 = 
CCOPTS = -O
MOPTS = -s
MBOPTS = -s
LOPTS =
LEGACY =

MPIOBJS = nullLOG.o
MPOBJS = MacMPxlf.o LnxMP.o

H5_DIR = /usr/local/apps/hdf5/current
SZ_DIR = /usr/local/apps/szip/current
HDF_DIR = /Users/uclapic/hdflib

INCPATH = -I$(HDF_DIR)/include -I$(H5_DIR)/include -I$(H5_DIR)/lib

LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
	-L$(H5_DIR)/lib -lhdf5 -lhdf5_fortran -lsz \
	-L$(HDF_DIR) -lz -ljpeg -ldf -lmfhdf \
	-L$(SZ_DIR)/lib


# Makefile IBM compiler with LAM MPI and MacOS X

#MPIFC = mpif77
#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPxlf.o LnxMP.o

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile g95 compiler with MacOS X

#MPIFC = g95
#FC90 = g95
#FC77 = g95
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -fstatic
#MBOPTS = -fstatic
#LOPTS =
#LEGACY =

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon -lSystemStubs

# Makefile Intel compiler with Linux

#MPIFC = ifc
#FC90 = ifc
#FC77 = ifc
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPIOBJS = MacMPIf77.o LnxMPI_S.o
#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

# Makefile Intel compiler with MPI and Linux

#MPIFC = mpif90
#FC90 = ifc
#FC77 = ifc
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

#

ESOBJS = globals.o pinit2mod.o prbpush2mod.o ppush2mod.o \
pfft2mod.o pfield2mod.o pdiag2mod.o p2mod.o p0mod.o mp0mod.o \
pinit2lib.o prbpush2lib.o ppush2lib.o pfft2lib.o pfield2lib.o \
pdiag2lib.o p2lib.o p0lib.o diag_jf.o pinit2mod_jf.o ampere_jf.o \
ext_driver_jf.o p2mod_jf.o hdf_write_jf.o par_track_jf.o

MESOBJS = mprbpush2mod.o mppush2mod.o mpfft2mod.o mp2mod.o \
mprbpush2lib.o mppush2lib.o mpfft2lib.o mp2lib.o 

EMOBJS = pbpush2mod.o pbpush2lib.o

MEMOBJS = mpbpush2mod.o mpbpush2lib.o

DESOBJS = pdfield2mod.o pbfield2mod.o pcfield2mod.o pnfield2mod.o \
pnpfield2mod.o pdfield2lib.o pbfield2lib.o pcfield2lib.o pnfield2lib.o 

NPOBJS = nullMP.o

# Linkage rules

all : new_d0_mpbbeps2.out new_d0_pbbeps2.out new_d0_mpbeps2.out new_d0_pbeps2.out \
      new_mpbbeps2.out new_pbbeps2.out new_mpbeps2.out new_pbeps2_jf.out

new_d0_mpbbeps2.out : new_d0_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(EMOBJS) $(MEMOBJS) \
                     $(DESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbbeps2.out \
        new_d0_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(EMOBJS) $(MEMOBJS) $(DESOBJS) \
        $(MPIOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_d0_pbbeps2.out : new_d0_pbbeps2.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(MPIOBJS) \
                     $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbbeps2.out \
	new_d0_pbbeps2.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) \
        $(GOBJS) $(LIBS)

new_d0_mpbeps2.out : new_d0_mpbeps2.o $(ESOBJS) $(DESOBJS) $(MESOBJS) $(MPIOBJS) \
                     $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbeps2.out \
	new_d0_mpbeps2.o $(ESOBJS) $(DESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) \
        $(GOBJS) $(LIBS)

new_d0_pbeps2.out : new_d0_pbeps2.o $(ESOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) \
                    $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbeps2.out \
	new_d0_pbeps2.o $(ESOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_mpbbeps2.out : new_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(EMOBJS) $(MEMOBJS) \
                   $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbbeps2.out \
        new_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(EMOBJS) $(MEMOBJS) $(MPIOBJS) \
        $(MPOBJS) $(GOBJS) $(LIBS)

new_pbbeps2.out : new_pbbeps2.o $(ESOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(OPTS90) $(LOPTS) -o new_pbbeps2.out \
	new_pbbeps2.o $(ESOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_mpbeps2.out : new_mpbeps2.o $(ESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbeps2.out \
        new_mpbeps2.o $(ESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_pbeps2_jf.out : new_pbeps2_jf.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) $(EMOBJS) mkdir_f_fxns.o -o new_pbeps2_jf.out \
	new_pbeps2_jf.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS) $(INCPATH)

# Compilation rules

MacMPIcf.o : MacMPIcf.c
	$(CC) $(CCOPTS) -c MacMPIcf.c

MacMPIf77.o : MacMPIf77.c
	$(CC) $(CCOPTS) -c MacMPIf77.c

MacMPIxlf.o : MacMPIxlf.c
	$(CC) $(CCOPTS) -c MacMPIxlf.c

MacMPI_S.o : MacMPI_S.c
	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMPI_S.c

LnxMPI_S.o : LnxMPI_S.c
	$(CC) $(CCOPTS) -c LnxMPI_S.c

MacMPcf.o : MacMPcf.c
	$(CC) $(CCOPTS) -c MacMPcf.c

MacMPf77.o : MacMPf77.c
	$(CC) $(CCOPTS) -c MacMPf77.c

MacMPxlf.o : MacMPxlf.c
	$(CC) $(CCOPTS) -c MacMPxlf.c

MacMP.o : MacMP.c
	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMP.c

LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c
	
mkdir_f_fxns.o : mkdir_f_fxns.c
	$(CC) $(CCOPTS) -c mkdir_f_fxns.c

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

libgks2.o : libgks2.f
	$(FC77) $(OPTS77) -c libgks2.f

dlibgks2.o : dlibgks2.f
	$(FC77) $(OPTS77) -c dlibgks2.f

nullpgks1.o : nullpgks1.f
	$(FC77) $(OPTS77) -c nullpgks1.f

nullpgks2.o : nullpgks2.f
	$(FC77) $(OPTS77) -c nullpgks2.f

nullLOG.o : nullLOG.f
	$(FC77) $(OPTS77) -c nullLOG.f

nullMP.o : nullMP.f
	$(FC77) $(OPTS77) -c nullMP.f

p0lib.o : p0lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c p0lib.f

mp2lib.o : mp2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mp2lib.f

p2lib.o : p2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c p2lib.f

pinit2lib.o : pinit2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pinit2lib.f

mppush2lib.o : mppush2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c mppush2lib.f

ppush2lib.o : ppush2lib.f
	$(FC77) $(OPTS77) -c ppush2lib.f

mpbpush2lib.o : mpbpush2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c mpbpush2lib.f

pbpush2lib.o : pbpush2lib.f
	$(FC77) $(OPTS77) -c pbpush2lib.f

mprbpush2lib.o : mprbpush2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c mprbpush2lib.f

prbpush2lib.o : prbpush2lib.f
	$(FC77) $(OPTS77) -c prbpush2lib.f

mpfft2lib.o : mpfft2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c mpfft2lib.f

pfft2lib.o : pfft2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pfft2lib.f

pfield2lib.o : pfield2lib.f
	$(FC77) $(OPTS77) -c pfield2lib.f

pdfield2lib.o : pdfield2lib.f
	$(FC90) $(OPTS77) -c pdfield2lib.f

pbfield2lib.o : pbfield2lib.f
	$(FC77) $(OPTS77) -c pbfield2lib.f

pcfield2lib.o : pcfield2lib.f
	$(FC77) $(OPTS77) -c pcfield2lib.f

pnfield2lib.o : pnfield2lib.f
	$(FC77) $(OPTS77) -c pnfield2lib.f

pdiag2lib.o : pdiag2lib.f
	$(FC77) $(OPTS77) -c pdiag2lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

pinit2mod.o : pinit2mod.f globals.o
	$(FC90) $(OPTS90) -c pinit2mod.f

pinit2mod_jf.o : pinit2mod_jf.f
	$(FC90) $(OPTS90) -c pinit2mod_jf.f
	
ampere_jf.o : ampere_jf.f
	$(FC90) $(OPTS90) -c ampere_jf.f
	
ext_driver_jf.o : ext_driver_jf.f
	$(FC90) $(OPTS90) -c ext_driver_jf.f

mppush2mod.o : mppush2mod.f ppush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mppush2mod.f

ppush2mod.o : ppush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c ppush2mod.f

mpbpush2mod.o : mpbpush2mod.f pbpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpbpush2mod.f

pbpush2mod.o : pbpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pbpush2mod.f

mprbpush2mod.o : mprbpush2mod.f prbpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mprbpush2mod.f

prbpush2mod.o : prbpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c prbpush2mod.f

mpfft2mod.o : mpfft2mod.f pfft2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpfft2mod.f

pfft2mod.o : pfft2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pfft2mod.f

pfield2mod.o : pfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pfield2mod.f

pdfield2mod.o : pdfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pdfield2mod.f

pbfield2mod.o : pbfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pbfield2mod.f

pcfield2mod.o : pcfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pcfield2mod.f

pnfield2mod.o : pnfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pnfield2mod.f

pnpfield2mod.o : pnpfield2mod.f pfield2mod.o pdfield2mod.o pbfield2mod.o \
                 pcfield2mod.o pnfield2mod.o
	$(FC90) $(OPTS90) -c pnpfield2mod.f

pdiag2mod.o : pdiag2mod.f pinit2mod.o
	$(FC90) $(OPTS90) -c pdiag2mod.f

diag_jf.o : diag_jf.f pinit2mod.o pinit2mod_jf.o p2mod_jf.o hdf_write_jf.o 
	$(FC90) $(OPTS90) -c diag_jf.f

p0mod.o : p0mod.f globals.o
	$(FC90) $(OPTS90) -c p0mod.f

mp2mod.o : mp2mod.f p2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mp2mod.f

p2mod.o : p2mod.f p0mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c p2mod.f

p2mod_jf.o : p2mod_jf.f p0mod.o p2mod.o par_track_jf.o
	$(FC90) $(OPTS90) $(LEGACY) -c p2mod_jf.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f
	
par_track_jf.o : par_track_jf.f pinit2mod_jf.o
	$(FC90) $(OPTS90) $(INCPATH) -c par_track_jf.f

hdf_write_jf.o : hdf_write_jf.f
	$(FC90) $(OPTS90) -c hdf_write_jf.f

new_pbeps2_jf.o : new_pbeps2_jf.f ppush2mod.o prbpush2mod.o pfft2mod.o \
               pbpush2mod.o pbpush2lib.o mkdir_f_fxns.o \
               pfield2mod.o pdiag2mod.o p2mod.o mp0mod.o p2mod_jf.o \
               diag_jf.o ampere_jf.o ext_driver_jf.o hdf_write_jf.o \
               par_track_jf.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbeps2_jf.f

new_mpbeps2.o : new_mpbeps2.f mprbpush2mod.o mppush2mod.o mpfft2mod.o \
                pfield2mod.o pdiag2mod.o mp2mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_mpbeps2.f

new_pbbeps2.o : new_pbbeps2.f prbpush2mod.o pbpush2mod.o \
                ppush2mod.o pfft2mod.o pfield2mod.o \
                pdiag2mod.o p2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbbeps2.f

new_mpbbeps2.o : new_mpbbeps2.f mprbpush2mod.o mpbpush2mod.o \
                 mppush2mod.o mpfft2mod.o pfield2mod.o \
                 pdiag2mod.o mp2mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_mpbbeps2.f

new_d0_pbeps2.o : new_d0_pbeps2.f ppush2mod.o prbpush2mod.o pfft2mod.o \
                  pnpfield2mod.o pdiag2mod.o p2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_pbeps2.f

new_d0_mpbeps2.o : new_d0_mpbeps2.f mprbpush2mod.o mppush2mod.o mpfft2mod.o \
                   pnpfield2mod.o pdiag2mod.o mp2mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_mpbeps2.f

new_d0_pbbeps2.o : new_d0_pbbeps2.f prbpush2mod.o pbpush2mod.o ppush2mod.o \
                   pfft2mod.o pnpfield2mod.o pdiag2mod.o p2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_pbbeps2.f

new_d0_mpbbeps2.o : new_d0_mpbbeps2.f mprbpush2mod.o mpbpush2mod.o \
                    mppush2mod.o mpfft2mod.o pnpfield2mod.o pdiag2mod.o \
                    mp2mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_mpbbeps2.f

clean :
# Delete -i new_beps2.out {OBJS}
