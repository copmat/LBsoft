
#include <default_macro.h>
 program LBsoft
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! LBsoft is a specific-purpose open-source software for soft glassy       
! emulsion simulations. The code was originally written by  
!
! This is an experimental code. The authors accept no responsibility
! for the performance of the code or for the correctness of the results.
!
! The code is licensed under Open Software License v. 3.0 (OSL-3.0). 
! The full text of the licence can be found on the website:
! http://opensource.org/licenses/OSL-3.0 
!
! A brief explanation of this license is available on the website:  
! http://rosenlaw.com/OSL3.0-explained.htm  
!
! The software development process has received funding from the 
! European Research Council under the Horizon 2020 Programme              
! Grant Agreement n. 739964 (COPMAT).   
!
! If results obtained with this code are published, an
! appropriate citation would be:
! 
!The code was originally written by                
!                                                                        
!Fabio Bonaccorso         IIT-CLNS, Rome                    Italy        
!Marco Lauricella         IAC-CNR, Rome                     Italy        
!Andrea Montessori        IAC-CNR, Rome                     Italy      
!                                                                        
!                                                                                                                                                
!with contributions from:                                                
!                                                                               
!Massimo Bernaschi        IAC-CNR, Rome                     Italy        
!Sauro Succi              IAC-CNR, Rome                     Italy        
!                                                                        
!
!                        LBsoft VERSION 0.01
!
! (July 2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use version_mod,    only : init_world,get_rank_world,get_size_world,&
                       time_world,time_world,finalize_world,idrank
  use profiling_mod,  only : get_memory,timer_init,itime_start, &
                       startPreprocessingTime,print_timing_partial, &
                       reset_timing_partial,printSimulationTime, &
                       print_timing_final,itime_counter,idiagnostic, &
                       ldiagnostic,start_timing2,end_timing2
  use utility_mod,    only : init_random_seed
  use fluids_mod,     only : allocate_fluids,initialize_fluids
  use integrator_mod, only : initime,endtime,tstep,set_nstep, &
                       update_nstep,nstep,driver_integrator,nstepmax
  use statistic_mod,         only : statistic_driver
  use io_mod
  
  
  implicit none
  
  
  real(kind=PRC) :: mytime=0.d0
  real(kind=PRC) :: itime,ctime,ftime
  
  real(kind=PRC) :: mymemory
  
  logical :: ladd,lrem,lremdat,ldorefinment,lrecycle
  
  integer :: i,j,k,atype
  
  

! set up the communications 
  call init_world()
  
! determine processor identities
  call get_rank_world()
  call get_size_world()
  
  !call alloc_domain()
  
! start clock
  call time_world(itime)
  
! print logo on terminal
  call print_logo(6)
  
! read the input file
  call read_input(600,'input.dat')
  
! set the seed
  call init_random_seed(init_seed)

! allocate arrays of the nanofiber quantities
  call allocate_fluids
  
! allocate service arrays for printing modules
  call allocate_print()
  
! initialize the counter of the integration steps  
  call set_nstep(0)
  
! initialize the time of the integration steps  
  mytime=initime
  
! initialize and read the restart file if requested
  call initialize_fluids
  stop
! print memory
  call get_memory(mymemory)
  call print_memory_registration(6,'memory occupied after allocation', &
   mymemory)
  
! open the  output 'statdat.dat' file and print first record on terminal
! and output file
  !call outprint_driver(nstep,mytime)
  
  
! open the binary file (only for developers) 
  !call open_dat_file(lprintdat,130,'traj.dat')
  
! start diagnostic if requested
  if(ldiagnostic)then
    call timer_init()
    call startPreprocessingTime()  
    call print_timing_partial(1,1,itime_start,IOOUT)
    call reset_timing_partial()
  endif
  
! initialize lrecycle 
  lrecycle=.true.
  
!***********************************************************************
!     start the time integration
!***********************************************************************
  do while (lrecycle)
  
!   update the counter
    call update_nstep
    
    mytime=real(nstep,kind=PRC)*tstep
    
!   check recycle loop
    lrecycle=(nstep<nstepmax)
    
!   integrate the system
    call driver_integrator(mytime)
    
    
!   compute statistical quanities
    call statistic_driver(nstep,mytime)
    
!   print on the binary file the jet bead which have hit the collector 
!   (only for developers)
    !call write_datrem_frame(lprintdat,122,nstep,mytime,iprintdat, &
    ! inpjet,npjet,sprintdat,systype,lremdat,nremoved)
    
    
!   print data on terminal
    if(ldiagnostic)call start_timing2("IO","outprint_driver")
    call outprint_driver(nstep,mytime)
    if(ldiagnostic)call end_timing2("IO","outprint_driver")
     
!   print the jet geometry on the binary file (only for developers)
    !call write_dat_frame(lprintdat,130,nstep,mytime,iprintdat, &
    ! inpjet,npjet,sprintdat,systype,linserted)
     
!   print restart file
    !call write_restart_file(nrestartdump,135,'save.dat',nstep,mytime)
    
!   cycle time check
    call time_world(ctime)
    
!   check recycle loop
    lrecycle=(lrecycle .and. timjob-ctime>timcls)
    
  enddo
!***********************************************************************
!     end of the time integration
!***********************************************************************
  
! stop clock
  call time_world(ftime)
  
! print final memory
  call get_memory(mymemory)
  call print_memory_registration(6,'memory occupied at the end', &
   mymemory)
  
! print last record on terminal and close the output 'statdat.dat' file
  !call finish_print(nstep,mytime,itime,ftime)
  
! print restart file
  !call write_restart_file(1,135,'save.dat',nstep,mytime)
  
! close the XYZ formatted output file 
  !call close_xyz_file(lprintxyz,120)

! close the binary file for hit beads (only for developers) 
  !call close_datrem_file(lprintdat,122)
    
! close the binary file (only for developers) 
  !call close_dat_file(lprintdat,130)
  
! finalize and print the diagnostic data
  if(ldiagnostic)then
    call printSimulationTime()
    call print_timing_final(idiagnostic,itime_counter,itime_start,1,1,IOOUT)
  endif
  
  if(idrank==0)write(6,'(a)')'Programm exit correctly'
  
! close the communications
  call finalize_world()
  
  stop

 end program LBsoft
  

