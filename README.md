# SPH-thixotropic
Code for simulating thixo-viscoplastic fluids using the SPH method.

 Code for simulating thixo-viscoplastic fluids using the SPH method.
 Reference:
   - Emanuele   Rossi;    Imanol   García   de    Beristáin,   Adolfo
     Vazquez-Quesada, José  Esteban López-Aguilar, Marco  Ellero. SPH
     simulations   of   thixo-viscoplastic    fluid   flow   past   a
     cylinder.  Journal   of  Non-Newtonian  Fluid   Mechanics,  308,
     104891. 2022

 Note  that  there  are  no inputs  for  the  non-Newtonian  viscosity
 equation.   The   viscosity   is   calculated   in   the   subroutine
 particles_compute_transport_coefficients,  and  any  changes  to  the
 viscosity model should be made directly in that subroutine.
 
 This code is  based on the original MCF code  developed by Xin Bian.
 The  current version  has  been developed  in collaboration  between
 - Marco Ellero,  leader of the  CFD Modelling and Simulation  group at
 BCAM (Basque Center  for Applied Mathematics) in  Bilbao, Spain
 - Emanuele Rossi, from Marco Ellero's group.
 - Adolfo Vazquez-Quesada from  the Department of Fundamental Physics
 at UNED, in Madrid, Spain.

 Developers:
     Xin Bian.
     Adolfo Vazquez-Quesada.
     Emanuele Rossi

 Contact: a.vazquez-quesada@fisfun.uned.es
          mellero@bcamath.org
--------------------------------------------------------------------

-------------------------- INSTALLATION -------------------------------
After download MCF/mcf/ folder,
I suppose you have already installed PPM library
in a proper path,
which needs to be given when you configure mcf client.
I also suppose you have installed makedepf90,
which is to check routines dependency of mcf/src/*.F90.

1) run script
 
  ./clean.sh

to clean all old files, which may be machine dependant.

2) run script

  ./bootstrap.sh

to make use of GNU autoconf + automake,
which generates configure script.

3) run configuration.

./configure --prefix=$HOME/MCF/mcf/mcf_install/ FC=ifort MPIFC=mpif90 LDFLAGS=-L$HOME/MCF/ppm/local/lib/ FCFLAGS=-I$HOME/MCF/ppm/local/include/ MAKEDEPF90=$HOME/MCF/ppm/local/bin/makedepf90

to generate Makefile, which is used to compile mcf code.

LDFLAGS: indicate ppm library object files.
FCFLAGS: indicate ppm library header files.


4) run compiler

  make -j 8

to compile mcf code,
-j 8 specify to use 8 processors to accelerate compiling.

5) run installation

  make install

to install mcf executable binary at ...../mcf/mcf_install/ folder.

-------------- USE -----------------------------------------------
Three input files are required to launch a simulation. Examples of
these input files can be found in the 'mcf_config' directory. The
details of the inputs are explained within those files.
-------------------------------------------------------------------
