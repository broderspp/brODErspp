# README for brODErs++

  For questions, write at: stefano.boccelli@polimi.it

##### ABOUT #####

  BrODErs++ is a collection of "brother" ODE solvers for chemically reacting
  hypersonic flows, developed at von Karman Institute for Fluid Dynamics (VKI).
  
  brODErs++ contains these solvers:
  
  > Shocking++   - thermochemical relaxation past a shockwave

  > LARSEN       - thermochemical refinement along a streamline

  > MultiLARSEN  - ***NOT INCLUDED*** thermochemical refinement along a multiple streamlines
                   ***ATTENTION:*** MultiLARSEN is not included in this release. Contact the authors
                   if you are interested to use it. We are always open to collaborations!

##### INSTALLING and RUNNING ####

  Please, refer to the .pdf file in the "doc" folder.
  Anyway, to compile:
  1) mkdir build; cd build
  2) ccmake ..   (and fill it)
  3) make install

##### CITING THIS SOFTWARE ######

-> If you use LARSEN/MultiLARSEN results in your publication, we kindly ask 
   you to cite:
  
   Boccelli, S., Bariselli, F., Dias, B., & Magin, T. E. (2019). Lagrangian 
   diffusive reactor for detailed thermochemical computations of plasma flows,
   Plasma Sources Science and Technology.
  
-> If you use Shocking results in your publication, a good paper to cite is:
  
   Magin, T. E., Caillault, L., Bourdon, A., & Laux, C. O. (2006). 
   Nonequilibrium radiative heat flux modeling for the Huygens entry probe. 
   Journal of Geophysical Research: Planets, 111(E7).

