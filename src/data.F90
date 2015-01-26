
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
  TYPE(time_stats) :: bookleaf_times

END MODULE timing_mod

MODULE op2_bookleaf

  USE OP2_Fortran_Reference
  use, intrinsic :: ISO_C_BINDING
  type(op_set) :: s_nodes, s_elements, s_mat, s_reg
  type(op_map) :: m_el2node,m_el2el, m_el2reg
  type(op_dat) :: d_rho,d_qq,d_elmass,d_elvol,d_ielmat,d_ein,d_pre,d_csqrd,      &
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

  INTEGER(kind=ink) :: ii,jj
  INTEGER(KIND=ink) :: ierr

  ! Let's declare OP2 stuff
  call op_init(0)
  call op_decl_set(nnod,s_nodes,'nodes')
  call op_decl_set(nel,s_elements,'elements')
  call op_decl_set(nreg,s_reg,'reg')

  ! Make a copy of the original mapping and decrement by 1, due to OP2's C backend
  ALLOCATE(ielnod2(nshape,1:nel))
  ielnod2(:,:) = ielnod(:,1:) - 1
  call op_decl_map(s_elements,s_nodes,nshape,ielnod2,m_el2node,'el2node')
  ALLOCATE(ielel2(nshape,1:nel))
  ielel2(:,:) = ielel(:,1:) - 1
  DO ii=1,nel
    DO jj=1,nshape
      IF (ielel2(jj,ii).eq.-1_ink) ielel2(jj,ii)= ii-1
    END DO
  END DO
  call op_decl_map(s_elements,s_elements,nshape,ielel2,m_el2el,'el2el')
!   ALLOCATE(ielmat2(1:nel))
!   ielmat2(:) = ielmat(1:) - 1
!   call op_decl_map(s_elements,s_mat,1,ielmat2,m_el2mat,'el2mat')
  ALLOCATE(ielreg2(1:nel))
  ielreg2(:) = ielreg(1:) - 1
  call op_decl_map(s_elements,s_reg,1,ielreg2,m_el2reg,'el2reg')

  ! This was originally not defined on reg, but had fixed size
  ALLOCATE(zdtnotreg2(1:nreg))
  zdtnotreg2 = zdtnotreg(1:nreg)
  ALLOCATE(zmidlength2(1:nreg))
  zmidlength2 = zmidlength(1:nreg)

  ! Generate a dataset for element indices
  ALLOCATE(elidx(1:nel))
  DO ii=1,nel
    elidx(ii) = ii
  END DO



!&          ielsort1(nel1),              &

  call op_decl_dat(s_elements,1,'real(8)',rho,d_rho,'rho')
  call op_decl_dat(s_elements,1,'real(8)',qq,d_qq,'qq')
  call op_decl_dat(s_elements,1,'real(8)',csqrd,d_csqrd,'csqrd')
  call op_decl_dat(s_elements,1,'real(8)',pre,d_pre,'pre')
  call op_decl_dat(s_elements,1,'real(8)',ein,d_ein,'ein')
  call op_decl_dat(s_elements,1,'real(8)',elmass,d_elmass,'elmass')
  call op_decl_dat(s_elements,1,'real(8)',elvol,d_elvol,'elvol')
  call op_decl_dat(s_elements,1,'real(8)',a1,d_a1,'a1')
  call op_decl_dat(s_elements,1,'real(8)',a2,d_a2,'a2')
  call op_decl_dat(s_elements,1,'real(8)',a3,d_a3,'a3')
  call op_decl_dat(s_elements,1,'real(8)',b1,d_b1,'b1')
  call op_decl_dat(s_elements,1,'real(8)',b2,d_b2,'b2')
  call op_decl_dat(s_elements,1,'real(8)',b3,d_b3,'b3')
  call op_decl_dat(s_elements,nshape,'real(8)',cnwt,d_cnwt,'cnwt')
  call op_decl_dat(s_elements,nshape,'real(8)',cnmass,d_cnmass,'cnmass')
  call op_decl_dat(s_elements,nshape,'real(8)',elx,d_elx,'elx')
  call op_decl_dat(s_elements,nshape,'real(8)',ely,d_ely,'ely')
  call op_decl_dat(s_elements,nshape,'real(8)',qx,d_qx,'qx')
  call op_decl_dat(s_elements,nshape,'real(8)',qy,d_qy,'qy')
  call op_decl_dat(s_elements,nshape,'integer(4)',ielsd,d_ielsd,'ielsd')
  call op_decl_dat(s_elements,nshape,'integer(4)',ielel,d_ielel,'ielel')
  call op_decl_dat(s_elements,1,'integer(4)',elidx,d_elidx,'elidx')
  call op_decl_dat(s_elements,1,'integer(4)',ielmat,d_ielmat,'ielmat')
  IF (ZSP) THEN
    call op_decl_dat(s_elements,nshape,'real(8)',spmass,d_spmass,'spmass')
  ENDIF

  call op_decl_dat(s_nodes,1,'real(8)',ndu,d_ndu,'ndu')
  call op_decl_dat(s_nodes,1,'real(8)',ndv,d_ndv,'ndv')
  call op_decl_dat(s_nodes,1,'real(8)',ndx,d_ndx,'ndx')
  call op_decl_dat(s_nodes,1,'real(8)',ndy,d_ndy,'ndy')
  call op_decl_dat(s_nodes,1,'real(8)',ndmass,d_ndmass,'ndmass')
  call op_decl_dat(s_nodes,1,'real(8)',ndarea,d_ndarea,'ndarea')
  call op_decl_dat(s_nodes,1,'integer(4)',indtype,d_indtype,'indtype')

  call op_decl_dat(s_elements,1,'real(8)',rscratch11,  d_rscratch11,'rscratch11')
  call op_decl_dat(s_elements,1,'real(8)',rscratch12,  d_rscratch12,'rscratch12')
  call op_decl_dat(s_elements,1,'real(8)',rscratch13,  d_rscratch13,'rscratch13')
  call op_decl_dat(s_nodes   ,1,'real(8)',rscratch14,  d_rscratch14,'rscratch14')
  call op_decl_dat(s_nodes   ,1,'real(8)',rscratch15,  d_rscratch15,'rscratch15')

!   call op_decl_dat(s_mat   ,1,'real(8)',mat_rho,  d_mat_rho,'d_mat_rho')
!   call op_decl_dat(s_mat   ,1,'real(8)',mat_ein,  d_mat_ein,'d_mat_ein')
!   call op_decl_dat(s_mat   ,1,'integer(4)',eos_type,  d_eos_type,'d_eos_type')
!   call op_decl_dat(s_mat   ,6,'real(8)',eos_param,  d_eos_param,'d_eos_param')

  call op_decl_dat(s_reg   ,1,'integer(4)',zdtnotreg2,  d_zdtnotreg,'d_zdtnotreg')
  call op_decl_dat(s_reg   ,1,'integer(4)',zmidlength2, d_zmidlength,'d_zmidlength')

  call op_decl_dat(s_elements,4,'real(8)',rscratch21,d_rscratch21,'rscratch21')
  call op_decl_dat(s_elements,4,'real(8)',rscratch22,d_rscratch22,'rscratch22')
  call op_decl_dat(s_elements,4,'real(8)',rscratch23,d_rscratch23,'rscratch23')
  call op_decl_dat(s_elements,4,'real(8)',rscratch24,d_rscratch24,'rscratch24')
  call op_decl_dat(s_elements,4,'real(8)',rscratch25,d_rscratch25,'rscratch25')
  call op_decl_dat(s_elements,4,'real(8)',rscratch26,d_rscratch26,'rscratch26')
  call op_decl_dat(s_elements,4,'real(8)',rscratch27,d_rscratch27,'rscratch27')

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
