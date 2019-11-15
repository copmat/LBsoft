# LBsoft

LBsoft is a specific-purpose open-source software for soft glassy       
emulsion simulations. The code was originally written by                
                                                                        
Fabio Bonaccorso         IIT-CLNS,Rome                     Italy        
Marco Lauricella         IAC-CNR, Rome                     Italy        
Andrea Montessori        IAC-CNR, Rome                     Italy      
                                                                        
                                                                                                                                                
with contributions from:                                                
                                                                        
Giorgio Amati         CINECA-CED, Rome                     Italy        
Massimo Bernaschi        IAC-CNR, Rome                     Italy        
Sauro Succi              IAC-CNR, Rome                     Italy        
                                                                        
The software development process has received funding from the          
European Research Council under the Horizon 2020 Programme              
Grant Agreement n. 739964 (COPMAT).                                     
                                                                        
This is an experimental code. The authors accept no responsibility      
for the performance of the code or for the correctness of the results.  
                                                                        
The code is licensed under the 3-Clause BSD License (BSD-3-Clause).     
                                                                        
Structure                                                               
                                                                        
LBsoft is supplied as a main UNIX directory with 3 subdirectories.      
All the source code files are contained in the 'source' sub-directory.  
The 'test_cases' sub-directory contains different test cases that can   
help the user to edit new input files. For further details on the input  
directives the reader is referred to the article:                       
                                                                        
Compiling LBsoft                                                        
                                                                        
The 'source' sub-directory stores a UNIX makefile that assembles the    
executable versions of the code both in serial and parallel version     
with different compilers. Note the makefile could be copied eventually  
modified for several common workstations and parallel computers into the  
'source' sub-directory, where the code is compiled and linked. Finally, 
the binary executable file can be run in the 'execute' sub-directory,   
which is intended to be the working directory from which jobs are       
submitted for execution and the data files manipulated. A list of       
targets for several common workstations and parallel computers can be   
used by the command "make target", where target is one of the           
following options:                                                      
gfortran          ---> compile in serial mode using the GFortran          
                   compiler.                                              
gfortran-mpi      ---> compile in serial mode using the GFortran          
                   compiler and the Open Mpi library.                     
intel             ---> compile in serial mode using the Intel compiler.   
intel-mpi         ---> compile in parallel mode using the Intel compiler  
                   and the Intel Mpi library.                           
intel-openmpi     ---> compile in parallel mode using the Intel compiler  
                   and the Open Mpi library.                            
intel-mpi-skylake ---> compile in parallel mode using Intel compiler    
                   and the Intel Mpi library with flags for             
                   Skylake processor features (AVX-512 activated).      
intel-mpi-knl     ---> compile in parallel mode using the Intel compiler  
                   and the Intel Mpi library with flags for Xeon Phi    
                   processor features (AVX-512 activated).              
help              ---> return the list of possible target choices.      
On Windows system we advice the user to compile LBsoft under the        
Windows subsystem for Linux (WSL). Note that a Message Passing Interface  
Implementation is necessary to compile and run LBsoft in parallel       
mode.                                                                   
                                                                        
Executing LBsoft                                                        
                                                                        
To run LBsoft, it is necessary first to ensure that the program is      
compiled (from the source sub-directory) and that the file input.dat    
is present in the execute subdirectory.                                 
All output data files will be returned to the 'execute' sub-directory.  
Remember that the input file HAS TO BE NAMED 'input.dat' and put in     
the 'execute' Sub-directory. Few example input files are contained in   
the 'test_cases' sub-directory which can be used as test cases.         
                                                                        
Example command to run LBsoft in serial mode:                           
                                                                        
./main.x                                                                
                                                                        
Example command to run LBsoft in parallel mode on 4 CPUs:               
                                                                        
mpirun -np 8 ./main_mpi.x                                               
                                                                        
