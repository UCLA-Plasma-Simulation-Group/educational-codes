#Makefile for 2D parallel PIC codes in new_pbeps2.source

GOBJS = nullpgks2.o nullpgks1.o


MPIFC = mpif90
FC90 = mpif90
FC77 = mpif77
CC = mpicc

OPTS90 = -fopenmp -fdefault-real-8 -O4
OPTS77 = -fopenmp -fdefault-real-8 -O4
#OPTS90 = -r8 -DFORTRANSINGLEUNDERSCORE -O3
#OPTS77 = -r8 -DFORTRANSINGLEUNDERSCORE -O3

CCOPTS = -O4
MOPTS = -s
MBOPTS = -s
LOPTS =
LEGACY =

MPIOBJS = nullLOG.o
MPOBJS = MacMPxlf.o LnxMP.o

HDF_DIR = /usr/local/
H5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/openmpi
SZ_DIR = /usr/local/

INCPATH = -I$(H5_DIR)/include -L$(H5_DIR)/lib 
#INCPATH = -I$(HDF_DIR)/include -I$(H5_DIR)/lib 

# LIBS = -L$(HDF_DIR)/lib -lz -ljpeg -ldf -lmfhdf \
#	-L$(H5_DIR)/lib -lhdf5 -lhdf5_fortran -lsz \
#	-L$(SZ_DIR)/lib
LIBS = -ljpeg -lz \
	-L$(H5_DIR)/lib -lsz -lhdf5_fortran -lhdf5 \
#	-L$(SZ_DIR)/lib

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
#LEGACY =

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

# Dawson Makefile mpif90 compiler with LAM MPI and MacOS X

#MPIFC = mpif77
#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = MacMPIxlf.o MacMPI_S.o
#MPOBJS = MacMPxlf.o MacMP.o

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPxlf.o LnxMP.o

#H5_DIR = /Users/uclapic/libs/hdf5/current
#SZ_DIR = /Users/uclapic/libs/szip/current
#JP_DIR = /Users/uclapic/libs/jpeg/current
#HDF_DIR = /Users/uclapic/libs/hdf4/current

#H5_DIR = /u/home2/ist/epp/hdf5/current
#SZ_DIR = /u/home2/ist/epp/szip/current/lib
#JP_DIR = /u/home2/ist/epp/jpeg/current/lib
#HDF_DIR = /u/home2/ist/epp/hdf4/current

#SZ_DIR = /Users/uclapic/libs/szip/current
#JP_DIR = /Users/uclapic/libs/jpeg/current
#HDF_DIR = /u/local/apps/hdf
#H5_DIR = /u/local/apps/hdf5/hdf5_old
#SZ_DIR = /u/local/apps/szip/current

#INCPATH = -I$(HDF_DIR)/include -I$(H5_DIR)/include -I$(H5_DIR)/lib \
#	-I/u/local/mpi/mpilam/include

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
#	-L$(HDF_DIR)/lib -lz -ljpeg -ldf -lmfhdf \
#	-L$(H5_DIR)/lib -lhdf5 -lhdf5_fortran -lsz \
#	-L$(SZ_DIR)/lib -L$(JP_DIR)/lib


#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
#	-L$(H5_DIR)/lib -lhdf5 -lhdf5_fortran -lsz \
#	-L$(HDF_DIR) -lz -ljpeg -ldf -lmfhdf \
#	-L$(SZ_DIR)/lib

#LIBS = /System/Library/Frameworks/Carbon.framework/Carbon \
#	-L$(H5_DIR)/lib -lhdf5 -lhdf5_fortran -lsz \
#	-L$(HDF_DIR) -lz -ljpeg -ldf -lmfhdf \
#	-L$(SZ_DIR)/lib -L$(JP_DIR)/lib


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
ext_driver_jf.o p2mod_jf.o h5lib_beps.o hdf_write_jf.o par_track_new_jf.o \
ppush2mod_jf.o ppush2lib_jf.o pfield2mod_jf.o

MESOBJS = mprbpush2mod.o mppush2mod.o mpfft2mod.o mp2mod.o \
mprbpush2lib.o mppush2lib.o mpfft2lib.o mp2lib.o 

EMOBJS = pbpush2mod.o pbpush2lib.o

MEMOBJS = mpbpush2mod.o mpbpush2lib.o

DESOBJS = pdfield2mod.o pbfield2mod.o pcfield2mod.o pnfield2mod.o \
pnpfield2mod.o pdfield2lib.o pbfield2lib.o pcfield2lib.o pnfield2lib.o 

NPOBJS = nullMP.o

# Linkage rules

all : upic-es.out upic-es-arb.out upic-es-mag.out


upic-es.out : upic-es.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) $(EMOBJS) mkdir_f_fxns.o -o upic-es.out \
	upic-es.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS) $(INCPATH)

upic-es-arb.out : upic-es-arb.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) $(EMOBJS) mkdir_f_fxns.o -o upic-es-arb.out \
	upic-es-arb.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS) $(INCPATH)

upic-es-mag.out : upic-es-mag.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) $(EMOBJS) mkdir_f_fxns.o -o upic-es-mag.out \
	upic-es-mag.o $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS) $(INCPATH)

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
	$(FC77) $(OPTS77) $(LEGACY) $(INCPATH) -c p0lib.f

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

ppush2lib_jf.o : ppush2lib_jf.f
	$(FC77) $(OPTS77) -c ppush2lib_jf.f

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
	$(FC90) $(OPTS90) -c -ffree-form pinit2mod_jf.f
	
ampere_jf.o : ampere_jf.f
	$(FC90) $(OPTS90) -c -ffree-form ampere_jf.f
	
ext_driver_jf.o : ext_driver_jf.f
	$(FC90) $(OPTS90) -c -ffree-form ext_driver_jf.f

mppush2mod.o : mppush2mod.f ppush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mppush2mod.f

ppush2mod.o : ppush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c ppush2mod.f

ppush2mod_jf.o : ppush2mod_jf.f p0mod.o
	$(FC90) $(OPTS90) -c -ffree-form ppush2mod_jf.f

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

pfield2mod_jf.o : pfield2mod_jf.f globals.o
	$(FC90) $(OPTS90) -c -ffree-form pfield2mod_jf.f

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
	$(FC90) $(OPTS90) $(INCPATH) -c -ffree-form diag_jf.f

p0mod.o : p0mod.f globals.o
	$(FC90) $(OPTS90) -c p0mod.f

mp2mod.o : mp2mod.f p2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mp2mod.f

p2mod.o : p2mod.f p0mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c p2mod.f

p2mod_jf.o : p2mod_jf.f p0mod.o p2mod.o par_track_new_jf.o
	$(FC90) $(OPTS90) $(LEGACY) -c -ffree-form p2mod_jf.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f
	
par_track_new_jf.o : par_track_new_jf.f pinit2mod_jf.o
	$(FC90) $(OPTS90) $(INCPATH) -c -ffree-form par_track_new_jf.f

h5lib_beps.o : h5lib_beps.f
	$(FC90) $(OPTS90) $(INCPATH) -c -ffree-form h5lib_beps.f

hdf_write_jf.o : h5lib_beps.o hdf_write_jf.f
	$(FC90) $(OPTS90) $(INCPATH) -c -ffree-form hdf_write_jf.f

upic-es.o : upic-es.f ppush2mod.o prbpush2mod.o pfft2mod.o \
               pbpush2mod.o pbpush2lib.o mkdir_f_fxns.o \
               pfield2mod.o pdiag2mod.o p2mod.o mp0mod.o h5lib_beps.o p2mod_jf.o \
               diag_jf.o ampere_jf.o ext_driver_jf.o hdf_write_jf.o \
               par_track_new_jf.o ppush2mod_jf.o pfield2mod_jf.o
	$(FC90) $(OPTS90) $(MOPTS) $(INCPATH) -c -ffree-form upic-es.f

upic-es-arb.o : upic-es-arb.f ppush2mod.o prbpush2mod.o pfft2mod.o \
               pbpush2mod.o pbpush2lib.o mkdir_f_fxns.o \
               pfield2mod.o pdiag2mod.o p2mod.o mp0mod.o h5lib_beps.o p2mod_jf.o \
               diag_jf.o ampere_jf.o ext_driver_jf.o hdf_write_jf.o \
               par_track_new_jf.o ppush2mod_jf.o pfield2mod_jf.o
	$(FC90) $(OPTS90) $(MOPTS) $(INCPATH) -c -ffree-form upic-es-arb.f

upic-es-mag.o : upic-es-mag.f ppush2mod.o prbpush2mod.o pfft2mod.o \
               pbpush2mod.o pbpush2lib.o mkdir_f_fxns.o \
               pfield2mod.o pdiag2mod.o p2mod.o mp0mod.o h5lib_beps.o p2mod_jf.o \
               diag_jf.o ampere_jf.o ext_driver_jf.o hdf_write_jf.o \
               par_track_new_jf.o ppush2mod_jf.o pfield2mod_jf.o
	$(FC90) $(OPTS90) $(MOPTS) $(INCPATH) -c -ffree-form upic-es-mag.f

clean :
# Delete -i new_beps2.out {OBJS}
