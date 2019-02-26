# educational-codes
Repository for illustrating and exploring the physics of educational plasma physics problems.

Initially the intent is to include ~1-page pdf files, each treating a specific problem in plasma physics (e.g., landau damping, two-stream instability, etc.).  Eventually we are going to create codes that allow for the variability of certain parameters to explore the physics behind each of these problems.

## Compiling the code on a Mac

The codes have a makefile for the GNU compiler, so if you are using UNIX you should use the linux-gnu.make makefile.  However, on the Mac the default "gcc" compiler for c is the LLVM compiler and behaves slightly different than the GNU compier and it also does not come  with Fortran.  This is the source of many minor errors when you compile under OS X.  To resolve this, make sure you know what c compiler you are using.  You can do this by typing the UNIX command:

which gcc

if it says "/usr/bin/gcc" then you are using LLVM C compiler.  If your gcc is located in a different folder, then you are using the GNU c compiler.  The makefile for OS X assumes you are using LLVM for C and gfortran for Fortran.  The following instructions assume you have installed XCODE and gcc in /usr/bin/gcc and installed gfortran in anoter location (for this example, we are assuming that you have installed it @ /usr/local/bin/gfortran).  

Running these codes require some flavor of MPI and HDF5.  Because we use parallel I/O in HDF5, HDF5 must be installed after MPI so the parallel option can be turned on.  So MPI has to be installed before HDF5.

To install MPI, download the source file, and go to the MPI folder and type the following commands:

./configure CC=/usr/bin/gcc F77=/usr/local/bin/gfortran F90=/usr/local/bin/gfortran --prefix=/usr/local \
make all \
make test \
sudo make install 

(the SUDO is important as you cannot put sofware to /usr/local without root priviledges).

After you've installed MPI, you should have the parallel compilers MPICC and MPIF90, to check it, type

which mpicc \
which mpif90 

and make sure these commands return something.  Then you can install the HDF5 library.  To do that, download HDF5 source and go to the HDF5 folder and type:

./configure CC=mpicc F77=mpif90 F90=mpif90  --enable-parallel --enable-fortran --prefix=/usr/local \
make all \
make test \
sudo make install 

Then you should be able to compile the educational softwares on your Mac.


FST 2/23/2019
