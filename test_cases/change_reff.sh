#!/bin/bash

myhome=$PWD


for i in 2D_Poiseuille_gradP_xy 2D_Poiseuille_gradP_xz 2D_Poiseuille_gradP_yz 2D_Poiseuille_xy 2D_Poiseuille_xz 2D_Poiseuille_yz test 2D_Shear_xy test 2D_Shear_xz test 2D_Shear_yz 3D_Spinodal 3D_Particle_pbc 3D_Rotating_Particle 3D_Shear_Particle 3D_Particle_SC 3D_Particles_box ; do
  if [ -d "$myhome/$i" ]; then
    cd $myhome/$i
    echo $PWD
    rm *.reff.dat
    if [ -f "dumpGlobal.save.dat" ]; then
      mv dumpGlobal.save.dat dumpGlobal.reff.dat
    fi
    
    if [ -f "dumpisFluid.save.dat" ]; then
      mv dumpisFluid.save.dat dumpisFluid.reff.dat
    fi
    
    if [ -f "dumpPart.save.dat" ]; then
      mv dumpPart.save.dat dumpPart.reff.dat
    fi
    
    if [ -f "dumpPopsR.save.dat" ]; then
      rm dumpPopsR.save.dat 
    fi
    
    if [ -f "dumpPopsB.save.dat" ]; then
      rm dumpPopsB.save.dat 
    fi
    
    if [ -f "dumpPopsR.reff.dat" ]; then
      rm dumpPopsR.reff.dat
    fi
    
    if [ -f "dumpPopsB.reff.dat" ]; then
      rm dumpPopsB.reff.dat
    fi
    
    if [ -f "dumpRhoR.save.dat" ]; then
      mv dumpRhoR.save.dat dumpRhoR.reff.dat
    fi
    
    if [ -f "dumpRhoB.save.dat" ]; then
      mv dumpRhoB.save.dat dumpRhoB.reff.dat
    fi
    
    if [ -f "dumpU.save.dat" ]; then
      mv dumpU.save.dat dumpU.reff.dat
    fi
    
    if [ -f "dumpV.save.dat" ]; then
      mv dumpV.save.dat dumpV.reff.dat
    fi
    
    if [ -f "dumpW.save.dat" ]; then
      mv dumpW.save.dat dumpW.reff.dat
    fi
    
  fi
  cd $myhome
  
done


exit


