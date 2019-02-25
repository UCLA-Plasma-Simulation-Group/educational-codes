# educational-codes
Repository for illustrating and exploring the physics of educational plasma physics problems.

Initially the intent is to include ~1-page pdf files, each treating a specific problem in plasma physics (e.g., landau damping, two-stream instability, etc.).  Eventually we are going to create codes that allow for the variability of certain parameters to explore the physics behind each of these problems.

## Compiling the code on a Mac

The codes have a makefile for the GNU compiler, so if you are using UNIX you should use the linux-gnu.make makefile.  However, on the Mac the gcc compiler for c is the LLVM compiler and behaves slightly different than the GNU compier and does not come with Fortran.  So in order to compile the code, you will need to install gfortran and by "accident" you have 2 different compilers on the same computer and that can be more trouble than it's worth.

