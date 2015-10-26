
!Crown Copyright 2014 AWE.
!
! This file is part of Bookleaf.
!
! Bookleaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! Bookleaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! Bookleaf. If not, see http://www.gnu.org/licenses/.

PROGRAM main

! Internal
  USE kinds_mod,    ONLY: ink
  USE error_mod,    ONLY: halt
  USE logicals_mod,  ONLY: zparallel,zmprocw
  USE integers_mod, ONLY: Nthread
  USE timing_mod,   ONLY: bookleaf_times, get_time
  USE timers_mod,   ONLY: start_timers

  USE write_mod,    ONLY: write_sprint,write_iprint
  USE mesh_mod,     ONLY: mesh_gen,mesh_transfer,regions
  use op2_bookleaf
#ifdef SILO
  USE silo_mod,     ONLY: write_silo_dump
#endif
! External
!#ifndef NOOMP
  USE omp_lib
!#endif

  IMPLICIT NONE

  ! mesh data
  TYPE(regions),DIMENSION(:),ALLOCATABLE :: reg
  INTEGER(kind=ink)                      :: nk,nl

! ###################
! Parallelism
! ###################

! MPI
!  CALL init_parallel()
  call op_init(0)
  IF (op_is_root()) THEN
    zmprocw = .TRUE._lok
  ENDIF

! OpenMP
#ifdef NOOMP
  Nthread=1_ink
#else
  Nthread=OMP_Get_Max_Threads()
#endif

! ###################
! TIMERS
! ###################

! start timers
  CALL start_timers()

! ###################
! BANNER
! ###################

! welcome banner
  IF (zmprocw) THEN
    CALL banner()
  ENDIF

! ###################
! DEFAULTS
! ###################

! initialise input
  CALL init_defaults()

! ###################
! INPUT
! ###################

! read command line
  CALL read_command()

! read input files
  CALL read_files()

! generate mesh from input
  CALL mesh_gen(reg,nk,nl)

! check / correct input
  CALL init_check()

! print input
  CALL write_iprint(reg)

! ###################
! INITIALISATION
! ###################

! setup run parameters from input
  CALL init_parameters()

! setup memory
!  CALL init_mesh_memory()

! Transfer mesh onto solution arrays, populate connectivity arrays
!  CALL mesh_transfer(reg)

! setup memory
!  CALL init_memory()

! main initialisation
  CALL init()

! problem specific modifications
#ifdef MODY
  CALL modify()
#endif

  bookleaf_times%time_end_init=get_time()

! print initial totals
  CALL write_sprint()

! Dump initial graphics file
#ifdef SILO
  CALL write_silo_dump("initial_dump")
#endif
#ifdef TIO
  CALL write_tio_dump("initial_dump.h5")
#endif

! ###################
! SOLVER
! ###################

! hydrodynamics
  CALL hydro()

! ###################
! FINISH
! ###################

  CALL halt("End time reached, terminating cleanly",1,zend=.true.)

END PROGRAM main
