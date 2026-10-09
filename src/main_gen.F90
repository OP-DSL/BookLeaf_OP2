
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
  USE kinds_mod,    ONLY: ink,lok,rlk
  USE error_mod,    ONLY: halt
  USE logicals_mod,  ONLY: zparallel,zmprocw
  USE integers_mod, ONLY: nshape,nel,nnod,nel2,nnod2,Nthread
  USE timing_mod,   ONLY: bookleaf_times, get_time
  USE timers_mod,   ONLY: start_timers
  USE pointers_mod, ONLY: ielmat,rho,ein,elmass,elvol,qq,qx,qy,pre,     &
&                         csqrd,ndx,ndy,elx,ely,ielel,ielnd,ielsd,cnwt, &
&                         cnmass,spmass,indtype
  USE geometry_mod, ONLY: getgeom,getgeom2
  USE utilities_mod,ONLY: getconn,getsconn,corrconn
  USE write_mod,    ONLY: write_sprint,write_iprint
  USE mesh_mod,     ONLY: mesh_gen,mesh_transfer,regions
  use op2_bookleaf
! External
!#ifndef NOOMP
  USE omp_lib
!#endif

  IMPLICIT NONE

  ! mesh data
  TYPE(regions),DIMENSION(:),ALLOCATABLE :: reg
  INTEGER(kind=ink)                      :: nk,nl
  integer(4) :: ii,jj,j1,j2,iel
  INTEGER(KIND=ink),DIMENSION(0:3) :: nodes

! ###################
! Parallelism
! ###################

! MPI
!  CALL init_parallel()
  !call op_init_base(0,0)
  IF (op_is_root().EQ.1_ink) THEN
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
  CALL init_mesh_memory()

! Transfer mesh onto solution arrays, populate connectivity arrays
  CALL mesh_transfer(reg)

! setup memory
  CALL init_memory()

! initialise connectivity
  ielel(1:,1:nel2)=getconn(nel2,nshape,ielnd(1:,1:nel2))
  ielsd(1:,1:nel2)=getsconn(nel2,nshape,ielel(1:,1:nel2))
  CALL corrconn(nel2,nshape,ielel(1:,1:nel2),ielsd(1:,1:nel2))

   ! initialise node type
  DO iel=1,nel2
   nodes(0:nshape-1)=ielnd(1:nshape,iel)
   IF (COUNT(indtype(nodes).LT.0_ink).EQ.3_ink) THEN
     l1:DO ii=0,nshape-1
       IF (indtype(nodes(ii)).GT.0_ink) EXIT l1
     ENDDO l1
     ii=MOD(ii+2_ink,nshape)
     jj=nodes(ii)
     IF (jj.LE.nnod) THEN
       j1=nodes(MOD(ii+1_ink,nshape))
       j2=nodes(MOD(ii+3_ink,nshape))
       IF (((indtype(j1).EQ.-2_ink).AND.(indtype(j2).EQ.-1_ink)).OR.     &
  &           ((indtype(j2).EQ.-2_ink).AND.(indtype(j1).EQ.-1_ink))) THEN
         indtype(jj)=-3_ink
       ENDIF
     ENDIF
   ENDIF
  ENDDO

  !Now everything is declared hopefully, we can pass it on to OP2
  call op2_bookleaf_declare
  call op_dump_to_hdf5 ("mesh_hdf5"//CHAR(0))

  CALL op_exit()

END PROGRAM main
