
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

MODULE kinds_mod

  INTEGER,PARAMETER :: ink=4,rlk=8,lok=4

END MODULE kinds_mod

MODULE parameters_mod

  USE kinds_mod,ONLY: ink,rlk

  INTEGER(KIND=ink),PARAMETER :: LN=80_ink
  INTEGER(KIND=ink),PARAMETER :: LI=100_ink
  INTEGER(KIND=ink),PARAMETER :: MAX_NAMELIST_SIZE=100_ink
  REAL(KIND=rlk),   PARAMETER :: ONEBYNINE=1.0_rlk/9.0_rlk
  INTEGER(KIND=ink),PARAMETER :: N_SHAPE=4_ink
  REAL(KIND=rlk),   PARAMETER :: pi       =3.1415926535897932385_rlk
  REAL(KIND=rlk),   PARAMETER :: two_pi   =6.2831853071795864770_rlk

END MODULE parameters_mod

MODULE integers_mod

  USE kinds_mod,     ONLY: ink
  USE parameters_mod,ONLY: LI

  INTEGER(KIND=ink)               :: nel,nnod,nshape,nmat,nreg,nstep,   &
&                                    nel1,nnod1,idtel,max_seg,max_subseg
  INTEGER(KIND=ink),DIMENSION(LI) :: eos_type

END MODULE integers_mod

MODULE reals_mod

  USE kinds_mod,     ONLY: rlk
  USE parameters_mod,ONLY: LI

  ! time
  REAL(KIND=rlk)                 :: time,time_start,time_end,dt_min,    &
&                                   dt_initial,dt_max,cfl_sf,div_sf,dt_g
  ! cut-off
  REAL(KIND=rlk)                 :: ccut,zcut,zerocut,pcut,dencut,accut
  ! q
  REAL(KIND=rlk)                 :: cq1,cq2
  ! eos
  REAL(KIND=rlk),DIMENSION(LI)   :: mat_rho,mat_ein
  REAL(KIND=rlk),DIMENSION(6,LI) :: eos_param
  ! hourglass
  REAL(KIND=rlk)                 :: kappaall,pmeritall
  REAL(KIND=rlk),DIMENSION(LI)   :: kappareg,pmeritreg

END MODULE reals_mod

MODULE strings_mod

  USE parameters_mod,ONLY: LN

  CHARACTER(LEN=LN) :: sfile

END MODULE strings_mod

MODULE logicals_mod

  USE kinds_mod,     ONLY: lok
  USE parameters_mod,ONLY: LI

  LOGICAL(KIND=lok)               :: zhg,zsp
  LOGICAL(KIND=lok),DIMENSION(LI) :: zdtnotreg,zmidlength

END MODULE logicals_mod

MODULE paradef_mod

  USE kinds_mod,ONLY: ink,lok,rlk

  INTEGER(KIND=ink)                            :: NprocW,rankW,CommS,   &
&                                                 CommW,Nthread
  LOGICAL(KIND=lok)                            :: zparallel,MprocW
  INTEGER(KIND=ink),DIMENSION(:),  ALLOCATABLE :: e_loc_glob,n_loc_glob,&
&                                                 ielsort1
  INTEGER(KIND=ink),DIMENSION(:,:),ALLOCATABLE :: e_owner_proc,         &
&                                                 n_owner_proc

END MODULE paradef_mod

MODULE pointers_mod

  USE kinds_mod,ONLY: ink,rlk

  INTEGER(KIND=ink),DIMENSION(:),  ALLOCATABLE        :: ielreg,ielreg2,ielmat, &
&                                                        indtype,elidx, &
&                                                        zdtnotreg2,zmidlength2
  INTEGER(KIND=ink),DIMENSION(:,:),ALLOCATABLE        :: ielel,ielsd
  INTEGER(KIND=ink),DIMENSION(:,:),ALLOCATABLE,TARGET :: ielnod
  INTEGER(KIND=ink), DIMENSION(:,:), ALLOCATABLE :: ielnod2,ielel2
  REAL(KIND=rlk),   DIMENSION(:),  ALLOCATABLE        :: rho,qq,csqrd,  &
&                                                        pre,ein,elmass,&
&                                                        elvol,a1,a2,a3,&
&                                                        b1,b2,b3,ndx,  &
&                                                        ndy,ndu,ndv, &
&                                                        ndmass,ndarea
  REAL(KIND=rlk),   DIMENSION(:,:),ALLOCATABLE        :: elx,ely,cnwt,  &
&                                                        qx,qy,spmass,  &
&                                                        cnmass

END MODULE pointers_mod

MODULE scratch_mod

  USE kinds_mod,ONLY: rlk

  REAL(KIND=rlk),DIMENSION(:),  ALLOCATABLE,TARGET :: rscratch11,       &
&                                                     rscratch12,       &
&                                                     rscratch13,       &
&                                                     rscratch14,       &
&                                                     rscratch15
  REAL(KIND=rlk),DIMENSION(:,:),ALLOCATABLE,TARGET :: rscratch21,       &
&                                                     rscratch22,       &
&                                                     rscratch23,       &
&                                                     rscratch24,       &
&                                                     rscratch25,       &
&                                                     rscratch26,       &
&                                                     rscratch27

END MODULE scratch_mod

MODULE timing_mod

  USE kinds_mod, ONLY: rlk

  TYPE time_stats
     REAL(KIND=rlk) :: time_start
     REAL(KIND=rlk) :: time_end
     REAL(KIND=rlk) :: time_end_main
     REAL(KIND=rlk) :: time_total
     REAL(KIND=rlk) :: time_end_init
     REAL(KIND=rlk) :: time_hydro
     REAL(KIND=rlk) :: time_in_lag
     REAL(KIND=rlk) :: time_in_getdt
     REAL(KIND=rlk) :: time_in_io
     REAL(KIND=rlk) :: time_step_io
     REAL(KIND=rlk) :: time_in_getq
     REAL(KIND=rlk) :: time_in_gethg
     REAL(KIND=rlk) :: time_in_getsp
     REAL(KIND=rlk) :: time_in_getacc
     REAL(KIND=rlk) :: time_in_getfrc
     REAL(KIND=rlk) :: time_in_getein
     REAL(KIND=rlk) :: time_in_eos
     REAL(KIND=rlk) :: time_in_geom
     REAL(KIND=rlk) :: time_in_comreg
     REAL(KIND=rlk) :: time_in_comms
     REAL(KIND=rlk) :: time_in_colls
  END TYPE time_stats
  TYPE(time_stats) :: bookleaf_times, get_time

  CONTAINS

  real(kind=rlk) function get_time()

  call CPU_TIME(get_time)

  end function get_time

END MODULE timing_mod

MODULE op2_bookleaf

  USE OP2_Fortran_Reference
  use OP2_Fortran_RT_Support
  use, intrinsic :: ISO_C_BINDING
  type(op_set) :: s_nodes, s_elements, s_mat, s_reg
  type(op_map) :: m_el2node,m_el2el, m_el2reg
  type(op_dat) :: d_rho,d_qq,d_elmass,d_elvol,d_ielmat,d_ielreg,d_ein,d_pre,d_csqrd,      &
&                           d_ndx,d_ndy,d_elx,d_ely,d_ndu,d_ndv, &
&                           d_a1,d_a2,d_a3,d_b1,d_b2,d_b3,d_cnwt,d_cnmass,d_qx,d_qy,d_indtype, &
&                           d_spmass,d_ielsd,d_ielel, d_elidx,&
&                           d_ndmass,d_ndarea, d_zdtnotreg, d_zmidlength
  type(op_dat) :: d_rscratch11,d_rscratch12,d_rscratch13,d_rscratch14, &
&                         d_rscratch15,d_rscratch21,d_rscratch22,d_rscratch23, &
&                         d_rscratch24,d_rscratch25,d_rscratch26,d_rscratch27

  PUBLIC :: op2_bookleaf_declare

  CONTAINS

  SUBROUTINE op2_bookleaf_declare

  USE kinds_mod,    ONLY: rlk,ink,lok
  USE parameters_mod,ONLY: LI
  USE logicals_mod,ONLY: zsp,zdtnotreg,zmidlength
  USE reals_mod,    ONLY: mat_rho,mat_ein,eos_param
  USE integers_mod, ONLY: nel,nnod,nshape,nel1,nnod1,eos_type,nreg
  USE pointers_mod, ONLY: rho,elmass,elvol,ielmat,ein,pre,csqrd,      &
&                         ndx,ndy,elx,ely,ndu,ndv,ielnod, ielnod2, &
&                         ielel, ielel2,qq,a1,a2,a3,b1,b2,b3,cnwt,cnmass, &
&                         qx,qy,indtype,spmass,ielsd,ndmass,ndarea, &
&                         ielreg2,zdtnotreg2,zmidlength2,ielreg,elidx
  USE scratch_mod,  ONLY: rscratch11,rscratch12,rscratch13,rscratch14, & !11-13 on elem 14,15 on nodes
&                         rscratch15,rscratch21,rscratch22,rscratch23, & !21-27 on elem
&                         rscratch24,rscratch25,rscratch26,rscratch27
  USE OP2_Fortran_hdf5_Declarations

  INTEGER(kind=ink) :: ii,jj
  INTEGER(KIND=ink) :: ierr, status
  character(kind=c_char,len=10) :: fileName
  character(kind=c_char,len=9) :: libName = C_CHAR_'PTSCOTCH'//C_NULL_CHAR
  character(kind=c_char,len=5) :: routineName = C_CHAR_'KWAY'//C_NULL_CHAR

  fileName = "mesh_hdf5"// CHAR(0)
  ! Let's declare OP2 stuff
  call op_decl_set_hdf5(s_nodes,fileName,'nodes')
  call op_decl_set_hdf5(s_elements,fileName,'elements')
  call op_decl_set_hdf5(s_reg,fileName,'reg')

  call op_decl_map_hdf5(s_elements,s_nodes,nshape,m_el2node,fileName,'el2node',status)
  call op_decl_map_hdf5(s_elements,s_elements,nshape,m_el2el,fileName,'el2el',status)
!   call op_decl_map_hdf5(s_elements,s_mat,1,ielmat2,m_el2mat,'el2mat')
  call op_decl_map_hdf5(s_elements,s_reg,1,m_el2reg,fileName,'el2reg',status)


  call op_decl_dat_hdf5(s_elements,1,d_rho,'real(8)',fileName,'rho',status)
  call op_decl_dat_hdf5(s_elements,1,d_qq,'real(8)',fileName,'qq',status)
  call op_decl_dat_hdf5(s_elements,1,d_csqrd,'real(8)',fileName,'csqrd',status)
  call op_decl_dat_hdf5(s_elements,1,d_pre,'real(8)',fileName,'pre',status)
  call op_decl_dat_hdf5(s_elements,1,d_ein,'real(8)',fileName,'ein',status)
  call op_decl_dat_hdf5(s_elements,1,d_elmass,'real(8)',fileName,'elmass',status)
  call op_decl_dat_hdf5(s_elements,1,d_elvol,'real(8)',fileName,'elvol',status)
  call op_decl_dat_hdf5(s_elements,1,d_a1,'real(8)',fileName,'a1',status)
  call op_decl_dat_hdf5(s_elements,1,d_a2,'real(8)',fileName,'a2',status)
  call op_decl_dat_hdf5(s_elements,1,d_a3,'real(8)',fileName,'a3',status)
  call op_decl_dat_hdf5(s_elements,1,d_b1,'real(8)',fileName,'b1',status)
  call op_decl_dat_hdf5(s_elements,1,d_b2,'real(8)',fileName,'b2',status)
  call op_decl_dat_hdf5(s_elements,1,d_b3,'real(8)',fileName,'b3',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_cnwt,'real(8)',fileName,'cnwt',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_cnmass,'real(8)',fileName,'cnmass',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_elx,'real(8)',fileName,'elx',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_ely,'real(8)',fileName,'ely',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_qx,'real(8)',fileName,'qx',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_qy,'real(8)',fileName,'qy',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_ielsd,'integer(4)',fileName,'ielsd',status)
  call op_decl_dat_hdf5(s_elements,nshape,d_ielel,'integer(4)',fileName,'ielel',status)
  call op_decl_dat_hdf5(s_elements,1,d_elidx,'integer(4)',fileName,'elidx',status)
  call op_decl_dat_hdf5(s_elements,1,d_ielmat,'integer(4)',fileName,'ielmat',status)
  call op_decl_dat_hdf5(s_elements,1,d_ielreg,'integer(4)',fileName,'ielreg',status)
  IF (ZSP) THEN
    call op_decl_dat_hdf5(s_elements,nshape,d_spmass,'real(8)',fileName,'spmass',status)
  ENDIF

  call op_decl_dat_hdf5(s_nodes,1,d_ndu,'real(8)',fileName,'ndu',status)
  call op_decl_dat_hdf5(s_nodes,1,d_ndv,'real(8)',fileName,'ndv',status)
  call op_decl_dat_hdf5(s_nodes,1,d_ndx,'real(8)',fileName,'ndx',status)
  call op_decl_dat_hdf5(s_nodes,1,d_ndy,'real(8)',fileName,'ndy',status)
  call op_decl_dat_hdf5(s_nodes,1,d_ndmass,'real(8)',fileName,'ndmass',status)
  call op_decl_dat_hdf5(s_nodes,1,d_ndarea,'real(8)',fileName,'ndarea',status)
  call op_decl_dat_hdf5(s_nodes,1,d_indtype,'integer(4)',fileName,'indtype',status)

  call op_decl_dat_hdf5(s_elements,1,d_rscratch11,'real(8)',fileName,'rscratch11',status)
  call op_decl_dat_hdf5(s_elements,1,d_rscratch12,'real(8)',fileName,'rscratch12',status)
  call op_decl_dat_hdf5(s_elements,1,d_rscratch13,'real(8)',fileName,'rscratch13',status)
  call op_decl_dat_hdf5(s_nodes,1,d_rscratch14,'real(8)',fileName,'rscratch14',status)
  call op_decl_dat_hdf5(s_nodes,1,d_rscratch15,'real(8)',fileName,'rscratch15',status)

  call op_decl_dat_hdf5(s_reg,1,d_zdtnotreg,'integer(4)',fileName,'d_zdtnotreg',status)
  call op_decl_dat_hdf5(s_reg,1,d_zmidlength,'integer(4)',fileName,'d_zmidlength',status)

  call op_decl_dat_hdf5(s_elements,4,d_rscratch21,'real(8)',fileName,'rscratch21',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch22,'real(8)',fileName,'rscratch22',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch23,'real(8)',fileName,'rscratch23',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch24,'real(8)',fileName,'rscratch24',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch25,'real(8)',fileName,'rscratch25',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch26,'real(8)',fileName,'rscratch26',status)
  call op_decl_dat_hdf5(s_elements,4,d_rscratch27,'real(8)',fileName,'rscratch27',status)

  call op_partition(libName, routineName, s_elements, &
     & m_el2node,  d_elx)

  END SUBROUTINE op2_bookleaf_declare
END MODULE op2_bookleaf

MODULE op2_constants

  USE reals_mod
  USE integers_mod

#ifdef OP2_ENABLE_CUDA
  USE CUDAFOR

  REAL(KIND=rlk),constant        :: dt_min_OP2, dt_initial_OP2,   &
&                                   dt_max_OP2,cfl_sf_OP2,div_sf_OP2,dt_g_OP2
  ! cut-off
  REAL(KIND=rlk), constant       :: ccut_OP2,zcut_OP2,zerocut_OP2, &
 &                                  pcut_OP2,dencut_OP2,accut_OP2
  INTEGER(KIND=ink), constant :: elements_stride_OP2, nodes_stride_OP2, &
&                                 reg_stride_OP2


  CONTAINS

  SUBROUTINE bookleaf_op2_init_const

    USE OP2_Fortran_Declarations
    USE op2_bookleaf

    implicit none

    dt_min_OP2 = dt_min
    dt_initial_OP2 = dt_initial
    dt_max_OP2 = dt_max
    cfl_sf_OP2 = cfl_sf
    div_sf_OP2 = div_sf
    dt_g_OP2 = dt_g
    ccut_OP2 = ccut
    zcut_OP2 = zcut
    zerocut_OP2 = zerocut
    pcut_OP2 = pcut
    dencut_OP2 = dencut
    accut_OP2 = accut
    elements_stride_OP2 = s_elements%setPtr%size + &
&    s_elements%setPtr%exec_size + s_elements%setPtr%nonexec_size
    nodes_stride_OP2 = s_nodes%setPtr%size + &
&    s_nodes%setPtr%exec_size + s_nodes%setPtr%nonexec_size
    reg_stride_OP2 = s_reg%setPtr%size + &
&    s_reg%setPtr%exec_size + s_reg%setPtr%nonexec_size
    call op_decl_const(dt_min, 1, 'dt_min')
    call op_decl_const(dt_initial, 1, 'dt_initial')
    call op_decl_const(dt_max, 1, 'dt_max')
    call op_decl_const(cfl_sf, 1, 'cfl_sf')
    call op_decl_const(div_sf, 1, 'div_sf')
    call op_decl_const(dt_g, 1, 'dt_g')
    call op_decl_const(ccut, 1, 'ccut')
    call op_decl_const(zcut, 1, 'zcut')
    call op_decl_const(zerocut, 1, 'zerocut')
    call op_decl_const(pcut, 1, 'pcut')
    call op_decl_const(dencut, 1, 'dencut')
    call op_decl_const(accut, 1, 'accut')

  END SUBROUTINE bookleaf_op2_init_const
#else

  CONTAINS
  SUBROUTINE bookleaf_op2_init_const
  END SUBROUTINE bookleaf_op2_init_const
#endif

END MODULE op2_constants
