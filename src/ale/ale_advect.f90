
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

MODULE ale_advect_mod

  USE kinds_mod,    ONLY: ink,rlk,lok
  USE timing_mod,   ONLY: bookleaf_times,get_time


  IMPLICIT NONE

  PRIVATE :: update_el_basis,update_el_var,aleadvect_el,                &
&            update_nd_basis,update_nd_var,aleadvect_nd
  PUBLIC  :: aleadvect

CONTAINS

  SUBROUTINE aleadvect(id1,id2,nshape,nel,nel1,nel2,nnod,nnod1,nnod2,   &
&                      nsz,ielel,ielsd,ielsrt1,ielsrt2,ielnd,indstatus, &
&                      indtype,dencut,cut,cutv,cutm,elv0ndm1,elm0ndm0,  &
&                      elr0ndv0,ndv1,elv1,elm1,elr1,cnv0,cnm1,dfv,dfm,  &
&                      cnm0,eluv,elvv,flux,work1,work2,work11,work21,zactive, &
&                      d_ielel,d_ielsd,d_indstatus, &
&                      d_indtype,d_cutv,d_cutm,d_elv0ndm1,d_elm0ndm0,  &
&                      d_elr0ndv0,d_ndv1,d_elv1,d_elm1,d_elr1,d_cnv0,d_cnm1,d_dfv,d_dfm,  &
&                      d_cnm0,d_eluv,d_elvv,d_flux,d_work1,d_work2,d_work11,d_work21,d_zactive)

    USE logicals_mod, ONLY: zparallel
    USE pointers_mod, ONLY: ndu,ndv
    USE utilities_mod,ONLY: gather,gather2
    use op2_bookleaf
    USE ale_advect_kernels
    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: id1,id2,   &
&                                                            nshape,nsz,&
&                                                            nel,nel1,  &
&                                                            nel2,nnod, &
&                                                            nnod1,nnod2
    REAL(KIND=rlk),                          INTENT(IN)   :: dencut,cut
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielel,     &
&                                                            ielsd,     &
&                                                            ielnd
    type(op_dat) :: d_ielel,d_ielsd,d_indstatus, &
&                      d_indtype,d_cutv,d_cutm,d_elv0ndm1,d_elm0ndm0,  &
&                      d_elr0ndv0,d_ndv1,d_elv1,d_elm1,d_elr1,d_cnv0,d_cnm1,d_dfv,d_dfm,  &
&                      d_cnm0,d_eluv,d_elvv,d_flux,d_work1,d_work2,d_work11,d_work21,d_zactive
    INTEGER(KIND=ink),DIMENSION(nel1),       INTENT(IN)   :: ielsrt1
    INTEGER(KIND=ink),DIMENSION(nel2),       INTENT(IN)   :: ielsrt2
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(IN)   :: indstatus, &
&                                                            indtype
    REAL(KIND=rlk),   DIMENSION(nsz),        INTENT(OUT)  :: cutv,cutm, &
&                                                            elv0ndm1,  &
&                                                            elm0ndm0,  &
&                                                            elr0ndv0,  &
&                                                            ndv1
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(INOUT):: elv1,elm1, &
&                                                            elr1,work11,work21
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(INOUT):: cnv0,cnm1, &
&                                                            dfv,dfm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT)  :: cnm0,work1,&
&                                                            work2,flux,&
&                                                            eluv,elvv
    !LOGICAL(KIND=lok),DIMENSION(nnod2),      INTENT(OUT)  :: zactive
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(OUT)  :: zactive
    ! Local
    REAL(KIND=rlk)                                        :: t0,t1

    ! Timer
    t0=get_time()

    ! Advect element quantities
    CALL aleadvect_el(id1,id2,nshape,nel,nel1,nel2,elv0ndm1(1),         &
&                     elm0ndm0(1),elr0ndv0(1),elv1(1),elm1(1),elr1(1),  &
&                     cutv(1),cutm(1),cnv0(1,1),cnm1(1,1),dfv(1,1),     &
&                     dfm(1,1),flux(1,1),ielel(1,1),ielsd(1,1),         &
&                     work11(1),work21(1), &
&                     d_elv0ndm1,         &
&                     d_elm0ndm0,d_elr0ndv0,d_elv1,d_elm1,d_elr1,  &
&                     d_cutv,d_cutm,d_cnv0,d_cnm1,d_dfv,     &
&                     d_dfm,d_flux,d_ielel,d_ielsd,         &
&                     d_work11,d_work21)


    CALL gather2(s_elements,m_el2node,d_ndu,d_eluv)
    CALL gather2(s_elements,m_el2node,d_ndv,d_elvv)

    ! Advect nodal quantities
    CALL aleadvect_nd(id1,id2,nshape,nel,nel1,nel2,nnod,nnod1,nnod2,nsz,&
&                     ielel(1,1),ielsd(1,1),ielsrt1(1),ielsrt2(1),      &
&                     zparallel,ielnd(1,1),indstatus(1),indtype(1),     &
&                     dencut,cut,cutv(1),cutm(1),elr0ndv0(1),ndv1(1),   &
&                     elm0ndm0(1),elv0ndm1(1),elv1(1),cnv0(1,1),        &
&                     cnm0(1,1),cnm1(1,1),dfv(1,1),dfm(1,1),eluv(1,1),  &
&                     elvv(1,1),work1(1,1),work2(1,1),flux(1,1),        &
&                     zactive(1), &
&                     d_ielel,d_ielsd,d_indstatus,d_indtype,    &
&                     d_cutv,d_cutm,d_elr0ndv0,d_ndv1,  &
&                     d_elm0ndm0,d_elv0ndm1,d_elv1,d_cnv0,       &
&                     d_cnm0,d_cnm1,d_dfv,d_dfm,d_eluv, &
&                     d_elvv,d_work1,d_work2,d_flux,       &
&                     d_zactive)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_aleadvect=bookleaf_times%time_in_aleadvect+t1

  END SUBROUTINE aleadvect

  SUBROUTINE aleadvect_el(id1,id2,nshape,nel,nel1,nel2,elvpr,elmpr,     &
&                         elrpr,elv,elm,elr,cutv,cutm,cnv,cnm,delv,delm,&
&                         flux,ielel,ielsd,work1,work2, &
&                     d_elvpr,         &
&                     d_elmpr,d_elrpr,d_elv,d_elm,d_elr,  &
&                     d_cutv,d_cutm,d_cnv,d_cnm,d_delv,     &
&                     d_delm,d_flux,d_ielel,d_ielsd,         &
&                     d_work1,d_work2)

    use op2_bookleaf
    USE ale_advect_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: id1,id2,   &
&                                                            nel,nel1,  &
&                                                            nel2,nshape
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielel,ielsd
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(OUT)  :: elvpr,     &
&                                                            cutv,cutm, &
&                                                            elrpr,     &
&                                                            elmpr
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(INOUT):: elv,elm,elr,work1,work2
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(IN)   :: delv,cnv,  &
&                                                            cnm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT)  :: delm,flux
    type(op_dat) ::      d_elvpr, d_elmpr,d_elrpr,d_elv,d_elm,d_elr,  &
&                     d_cutv,d_cutm,d_cnv,d_cnm,d_delv,     &
&                     d_delm,d_flux,d_ielel,d_ielsd,         &
&                     d_work1,d_work2
    ! Local
    REAL(KIND=rlk)                                        :: t0,t1

    ! Timer
    t0=get_time()

    ! update element basis variables
    CALL update_el_basis(id1,id2,nshape,nel,nel1,nel2,elvpr(1),elmpr(1),&
&                        elrpr(1),elv(1),elm(1),elr(1),cutv(1),cutm(1), &
&                        cnv(1,1),cnm(1,1),delv(1,1),delm(1,1),         &
&                        ielel(1,1),ielsd(1,1),work1(1),work2(1), &
&                        d_elvpr,d_elmpr,&
&                        d_elrpr,d_elv,d_elm,d_elr,d_cutv,d_cutm, &
&                        d_cnv,d_cnm,d_delv,d_delm,         &
&                        d_ielel,d_ielsd,d_work1,d_work2)

    ! update element independent variables
    CALL update_el_var(id1,id2,nshape,nel,nel1,nel2,ielel(1,1),         &
&                      ielsd(1,1),elvpr(1),elmpr(1),elv(1),elm(1),      &
&                      cutv(1),cutm(1),cnv(1,1),cnm(1,1),delv(1,1),     &
&                      delm(1,1),flux(1,1),work1(1), &
&                      d_ielel,         &
&                      d_ielsd,d_elvpr,d_elmpr,d_elv,d_elm,      &
&                      d_cutv,d_cutm,d_cnv,d_cnm,d_delv,     &
&                      d_delm,d_flux,d_work1)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_aleadvect_el=                                &
&    bookleaf_times%time_in_aleadvect_el+t1

  END SUBROUTINE aleadvect_el

  SUBROUTINE aleadvect_nd(id1,id2,nshape,nel,nel1,nel2,nnod,nnod1,nnod2,&
&                         nsz,ielel,ielsd,ielsrt1,ielsrt2,zparallel,    &
&                         ielnd,indstatus,indtype,dencut,cut,cutv,cutm, &
&                         ndv0,ndv1,ndm0,elv0ndm1,elv1,cnv0,cnm0,cnm1,  &
&                         dfv,dfm,eluv,elvv,dcv,dcm,flux,zactive, &
&                     d_ielel,d_ielsd,d_indstatus,d_indtype,    &
&                     d_cutv,d_cutm,d_ndv0,d_ndv1,  &
&                     d_ndm0,d_elv0ndm1,d_elv1,d_cnv0,       &
&                     d_cnm0,d_cnm1,d_dfv,d_dfm,d_eluv, &
&                     d_elvv,d_dcv,d_dcm,d_flux,       &
&                     d_zactive)

    use op2_bookleaf
    USE ale_advect_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: nshape,nel,&
&                                                            nel1,nel2, &
&                                                            nnod,nnod1,&
&                                                            nnod2,nsz, &
&                                                            id1,id2
    REAL(KIND=rlk),                          INTENT(IN)   :: dencut,cut
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielel,     &
&                                                            ielsd,     &
&                                                            ielnd
    INTEGER(KIND=ink),DIMENSION(nel1),       INTENT(IN)   :: ielsrt1
    INTEGER(KIND=ink),DIMENSION(nel2),       INTENT(IN)   :: ielsrt2
    LOGICAL(KIND=lok),                       INTENT(IN)   :: zparallel
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(IN)   :: indstatus, &
&                                                            indtype
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(IN)   :: elv1
    REAL(KIND=rlk),   DIMENSION(nnod2),      INTENT(OUT)  :: ndv0,ndm0, &
&                                                            ndv1,cutv, &
&                                                            cutm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(IN)   :: cnv0,eluv, &
&                                                            elvv
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(INOUT):: dfv,dfm,   &
&                                                            cnm1
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT)  :: dcv,dcm,   &
&                                                            cnm0,flux
    REAL(KIND=rlk),   DIMENSION(nsz),        INTENT(INOUT):: elv0ndm1
    !LOGICAL(KIND=lok),DIMENSION(nnod2),      INTENT(OUT)  :: zactive
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(OUT)  :: zactive
    type(op_dat) ::      d_ielel,d_ielsd,d_indstatus,d_indtype,    &
&                     d_cutv,d_cutm,d_ndv0,d_ndv1,  &
&                     d_ndm0,d_elv0ndm1,d_elv1,d_cnv0,       &
&                     d_cnm0,d_cnm1,d_dfv,d_dfm,d_eluv, &
&                     d_elvv,d_dcv,d_dcm,d_flux,       &
&                     d_zactive
    ! Local
    REAL(KIND=rlk)                                        :: t0,t1

    ! Timer
    t0=get_time()

    ! update nodal basis variables
    CALL update_nd_basis(id1,id2,nshape,nel2,nnod2,nsz,dencut,cut,      &
&                        ielel(1,1),ielsd(1,1),ielsrt2(1),zparallel,    &
&                        ielnd(1,1),dfv(1,1),dfm(1,1),dcv(1,1),dcm(1,1),&
&                        cnm0(1,1),cnm1(1,1),cutv(1),cutm(1),ndv0(1),   &
&                        ndv1(1),ndm0(1),elv0ndm1(1),elv1(1),flux(1,1), &
&                        d_ielel,d_ielsd,  &
&                        d_dfv,d_dfm,d_dcv,d_dcm, &
&                        d_cnm0,d_cnm1,d_cutv,d_cutm,d_ndv0,  &
&                        d_ndv1,d_ndm0,d_elv0ndm1,d_elv1,d_flux)

    ! update nodal independent variables
    CALL update_nd_var(nshape,nel,nel1,nel2,nnod,nnod2,ielel(1,1),      &
&                      ielsd(1,1),ielnd(1,1),ielsrt1(1),zparallel,      &
&                      indstatus(1),indtype(1),ndv0(1),ndm0(1),ndv1(1), &
&                      elv0ndm1(1),cutv(1),cutm(1),cnv0(1,1),cnm0(1,1), &
&                      dcv(1,1),dcm(1,1),flux(1,1),eluv(1,1),elvv(1,1), &
&                      dfm(1,1),zactive(1), &
&                      d_ielel,d_ielsd,      &
&                      d_indstatus,d_indtype,d_ndv0,d_ndm0,d_ndv1, &
&                      d_elv0ndm1,d_cutv,d_cutm,d_cnv0,d_cnm0, &
&                      d_dcv,d_dcm,d_flux,d_eluv,d_elvv, &
&                      d_dfm,d_zactive)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_aleadvect_nd=                                &
&    bookleaf_times%time_in_aleadvect_nd+t1

  END SUBROUTINE aleadvect_nd

  SUBROUTINE update_el_basis(id1,id2,nshape,nel,nel1,nel2,elvpr,elmpr,  &
&                            elrpr,elv,elm,elr,cutv,cutm,cnv,cnm,delv,  &
&                            delm,ielel,ielsd,totv,totm, &
&                        d_elvpr,d_elmpr,&
&                        d_elrpr,d_elv,d_elm,d_elr,d_cutv,d_cutm, &
&                        d_cnv,d_cnm,d_delv,d_delm,         &
&                        d_ielel,d_ielsd,d_totv,d_totm)

    USE reals_mod,        ONLY: dencut,zerocut
    USE ale_advectors_mod,ONLY: flux_c1_VL,sum_flux
    USE op2_bookleaf
    USE ale_advect_kernels
    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: id1,id2,   &
&                                                            nshape,nel,&
&                                                            nel1,nel2
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielel,ielsd
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(OUT)  :: cutv,cutm, &
&                                                            elvpr,totv,&
&                                                            elmpr,totm,&
&                                                            elrpr
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(INOUT):: elv,elm,elr
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(IN)   :: delv,cnv,  &
&                                                            cnm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT)  :: delm
    type(op_dat) ::         d_elvpr,d_elmpr,&
&                        d_elrpr,d_elv,d_elm,d_elr,d_cutv,d_cutm, &
&                        d_cnv,d_cnm,d_delv,d_delm,         &
&                        d_ielel,d_ielsd,d_totv,d_totm
    ! Local
    INTEGER(KIND=ink) :: iel
    REAL(KIND=rlk)    :: t0,t1

    ! Timer
    t0=get_time()

    ! calculate total volume flux to nel
    CALL sum_flux(id1,id2,nshape,nel,nel1,ielel(1,1),ielsd(1,1),   &
&                 delv(1,1),totv(1), &
&                 d_ielel,d_ielsd,   &
&                 d_delv,d_totv)

    ! construct mass flux top nel1
    CALL flux_c1_VL(id1,id2,nshape,nel1,nel2,ielel(1,1),ielsd(1,1),     &
&                   cnv(1,1),delv(1,1),elr(1),delm(1,1), &
&                   d_ielel,d_ielsd,     &
&                   d_cnv,d_delv,d_elr,d_delm)

    ! calculate total mass flux to nel
    CALL sum_flux(id1,id2,nshape,nel,nel1,ielel(1,1),ielsd(1,1),   &
&                 delm(1,1),totm(1), &
&                 d_ielel,d_ielsd,d_delm,d_totm)

    ! update
    call op_par_loop_10(ale_advect_update,s_elements, &
&        op_arg_dat(d_elv,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_elm,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_elr,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_elvpr,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_elmpr,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_elrpr,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_cutv,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_cutm,-1,OP_ID,1,'real(8)',OP_RW), &
&        op_arg_dat(d_totv,-1,OP_ID,1,'real(8)',OP_READ), &
&        op_arg_dat(d_totm,-1,OP_ID,1,'real(8)',OP_READ))

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_update_el_basis=                             &
&    bookleaf_times%time_in_update_el_basis+t1

  END SUBROUTINE update_el_basis

  SUBROUTINE update_el_var(id1,id2,nshape,nel,nel1,nel2,ielel,ielsd,    &
&                          elvpr,elmpr,elv,elm,cutv,cutm,cnv,cnm,delv,  &
&                          delm,flux,tflux, &
&                      d_ielel,         &
&                      d_ielsd,d_elvpr,d_elmpr,d_elv,d_elm,      &
&                      d_cutv,d_cutm,d_cnv,d_cnm,d_delv,     &
&                      d_delm,d_flux,d_tflux)

    USE ale_advectors_mod,ONLY: flux_c1_VL,update_c1
    USE pointers_mod,     ONLY: ein
    USE op2_bookleaf
    USE ale_advect_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN) :: id1,id2,nel, &
&                                                          nel1,nel2,   &
&                                                          nshape
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN) :: ielel,ielsd
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(IN) :: elvpr,elmpr, &
&                                                          elv,elm,cutv,&
&                                                          cutm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(IN) :: cnv,cnm,delv,&
&                                                          delm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT):: flux
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(OUT):: tflux
    type(op_dat) ::                       d_ielel,         &
&                      d_ielsd,d_elvpr,d_elmpr,d_elv,d_elm,      &
&                      d_cutv,d_cutm,d_cnv,d_cnm,d_delv,     &
&                      d_delm,d_flux,d_tflux
    ! Local
    REAL(KIND=rlk)                                      :: t0,t1

    ! Timer
    t0=get_time()

    ! internal energy (mass weighted)
    CALL flux_c1_VL(id1,id2,nshape,nel1,nel2,ielel(1,1),ielsd(1,1),     &
&                   cnm(1,1),delm(1,1),ein(1),flux(1,1), &
&                   d_ielel,d_ielsd,     &
&                   d_cnm,d_delm,d_ein,d_flux)
    CALL update_c1(id1,id2,nshape,nel,nel2,ielel(1,1),ielsd(1,1),       &
&                  elmpr(1),elm(1),cutm(1),flux(1,1),tflux(1),ein(1), &
&                  d_ielel,d_ielsd,       &
&                  d_elmpr,d_elm,d_cutm,d_flux,d_tflux,d_ein)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_update_el_var=                               &
&    bookleaf_times%time_in_update_el_var+t1

  END SUBROUTINE update_el_var

  SUBROUTINE update_nd_basis(id1,id2,nshape,nel2,nnod2,nsz,dencut,cut,  &
&                            ielel,ielsd,ielsrt,zparallel,ielnd,delv,   &
&                            delm,dndv,dndm,cnm0,cnm1,cutv,cutm,ndv0,   &
&                            ndv1,ndm0,elv0ndm1,elv1,flux, &
&                        d_ielel,d_ielsd,  &
&                        d_delv,d_delm,d_dndv,d_dndm, &
&                        d_cnm0,d_cnm1,d_cutv,d_cutm,d_ndv0,  &
&                        d_ndv1,d_ndm0,d_elv0ndm1,d_elv1,d_flux)
    USE op2_bookleaf
    USE common_kernels
    USE ale_advect_kernels
    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: id1,id2,   &
&                                                            nel2,nnod2,&
&                                                            nsz,nshape
    REAL(KIND=rlk),                          INTENT(IN)   :: dencut,cut
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielel,     &
&                                                            ielsd,     &
&                                                            ielnd
    INTEGER(KIND=ink),DIMENSION(nel2),       INTENT(IN)   :: ielsrt
    LOGICAL(KIND=lok),                       INTENT(IN)   :: zparallel
    REAL(KIND=rlk),   DIMENSION(nel2),       INTENT(IN)   :: elv1
    REAL(KIND=rlk),   DIMENSION(nshape,nsz), INTENT(IN)   :: delv,delm
    REAL(KIND=rlk),   DIMENSION(nshape,nsz), INTENT(OUT)  :: dndv,dndm, &
&                                                            cnm0,flux
    REAL(KIND=rlk),   DIMENSION(nsz),        INTENT(OUT)  :: ndv0,ndv1, &
&                                                            ndm0,cutv, &
&                                                            cutm
    REAL(KIND=rlk),   DIMENSION(nsz),        INTENT(INOUT):: elv0ndm1
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(INOUT):: cnm1
    type(op_dat)         d_ielel,d_ielsd,  &
&                        d_delv,d_delm,d_dndv,d_dndm, &
&                        d_cnm0,d_cnm1,d_cutv,d_cutm,d_ndv0,  &
&                        d_ndv1,d_ndm0,d_elv0ndm1,d_elv1,d_flux

    ! Local
    INTEGER(KIND=ink) :: ind,iel,ii,i1,i2,ie1,ie2,is1,is2
    REAL(KIND=rlk)    :: w1,w2,w3,w4,t0,t1

    ! Timer
    t0=get_time()

    ! initialise
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndv0,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndv1,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndm0,-1,OP_ID,1,'real(8)',OP_WRITE))


    ! construct pre/post nodal volumes and pre nodal/corner mass
    call op_par_loop_16(ale_advect_prevolmass,s_elements, &
&           op_arg_dat(d_elv0ndm1,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_elv1,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_cnm0,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_dat(d_cnm1,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_ndv0, 1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv1, 1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndm0, 1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv0, 2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv1, 2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndm0, 2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv0, 3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv1, 3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndm0, 3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv0, 4,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndv1, 4,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndm0, 4,m_el2node,1,'real(8)',OP_INC))

    ! construct volume and mass flux
    call op_par_loop_1(set_zero4,s_elements, &
&           op_arg_dat(d_flux,-1,OP_ID,4,'real(8)',OP_WRITE))
    DO i1=id1,id2
    call op_par_loop_17(ale_advect_volmass, s_elements, &
&           op_arg_dat(d_delv,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delv, 1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delv, 2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delv, 3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delv, 4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delm,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delm, 1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delm, 2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delm, 3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_delm, 4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_dndv,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_dat(d_dndm,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_dat(d_flux,-1,OP_ID,4,'real(8)',OP_INC), &
&           op_arg_dat(d_ielsd,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_dat(d_ielel,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_dat(d_elidx,-1,OP_ID,1,'integer(4)',OP_READ), &
&           op_arg_gbl(i1,1,'integer(4)',OP_READ))
    ENDDO

    ! construct post nodal/corner mass

    call op_par_loop_2(a_eq_b,s_nodes, &
&           op_arg_dat(d_elv0ndm1,-1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_ndm0,    -1,OP_ID,1,'real(8)',OP_READ))
    call op_par_loop_6(ale_advect_postmass, s_elements, &
&           op_arg_dat(d_cnm1,   -1,OP_ID,4,'real(8)',OP_INC), &
&           op_arg_dat(d_flux,   -1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_elv0ndm1,1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_elv0ndm1,2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_elv0ndm1,3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_elv0ndm1,4,m_el2node,1,'real(8)',OP_INC))

    ! construct cut-offs
    call op_par_loop_4(ale_advect_cutoff,s_nodes, &
&           op_arg_dat(d_cutv,   -1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_cutm,   -1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_ndv0,   -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_gbl(cut,1,'real(8)',OP_READ))

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_update_nd_basis=                             &
&    bookleaf_times%time_in_update_nd_basis+t1

  END SUBROUTINE update_nd_basis

  SUBROUTINE update_nd_var(nshape,nel,nel1,nel2,nnod,nnod2,ielel,ielsd, &
&                          ielnd,ielsrt,zparallel,indstatus,indtype,    &
&                          ndv0,ndm0,ndv1,ndm1,cutv,cutm,cnv,cnm,delv,  &
&                          delm,flux,eluv,elvv,tflux,zactive, &
&                      d_ielel,d_ielsd,      &
&                      d_indstatus,d_indtype,d_ndv0,d_ndm0,d_ndv1, &
&                      d_ndm1,d_cutv,d_cutm,d_cnv,d_cnm, &
&                      d_delv,d_delm,d_flux,d_eluv,d_elvv, &
&                      d_tflux,d_zactive)

    USE ale_advectors_mod,ONLY: flux_n1_VL,update_n1
    USE pointers_mod,     ONLY: ndu,ndv
    USE op2_bookleaf
    USE ale_advect_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN) :: nshape,nel,  &
&                                                          nel1,nel2,   &
&                                                          nnod,nnod2
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN) :: ielel,ielsd, &
&                                                          ielnd
    INTEGER(KIND=ink),DIMENSION(nel1),       INTENT(IN) :: ielsrt
    LOGICAL(KIND=lok),                       INTENT(IN) :: zparallel
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(IN) :: indstatus,   &
&                                                          indtype
    REAL(KIND=rlk),   DIMENSION(nnod2),      INTENT(IN) :: ndv0,ndm0,   &
&                                                          ndv1,ndm1,   &
&                                                          cutv,cutm
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(IN) :: cnv,cnm,delv,&
&                                                          delm,eluv,   &
&                                                          elvv
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT):: flux
    REAL(KIND=rlk),   DIMENSION(nnod2),      INTENT(OUT):: tflux
    !LOGICAL(KIND=lok),DIMENSION(nnod2),      INTENT(OUT):: zactive
    INTEGER(KIND=ink),DIMENSION(nnod2),      INTENT(OUT)  :: zactive
    type(op_dat) ::    d_ielel,d_ielsd,      &
&                      d_indstatus,d_indtype,d_ndv0,d_ndm0,d_ndv1, &
&                      d_ndm1,d_cutv,d_cutm,d_cnv,d_cnm, &
&                      d_delv,d_delm,d_flux,d_eluv,d_elvv, &
&                      d_tflux,d_zactive
    ! Local
    INTEGER(KIND=ink) :: ind
    REAL(KIND=rlk)    :: t0,t1

    ! Timer
    t0=get_time()

    ! momentum (mass weighted)

!! gather here must happen before comms and out to nel (needs nnod1) and comm eluv
!    CALL gather(nshape,nel,nnod,ielnd(1,1),ndu(1),eluv(1,1))
    call op_par_loop_3(ale_advect_markactive,s_nodes, &
&           op_arg_dat(d_indstatus, -1, OP_ID, 1, 'integer(4)', OP_READ), &
&           op_arg_dat(d_indtype,   -1, OP_ID, 1, 'integer(4)', OP_READ), &
&           op_arg_dat(d_zactive,   -1, OP_ID, 1, 'integer(4)', OP_WRITE))

    CALL flux_n1_VL(nshape,nel1,nel2,ielel(1,1),ielsd(1,1),cnm(1,1),    &
&                   delm(1,1),eluv(1,1),flux(1,1), &
&                   d_ielel,d_ielsd,d_cnm,    &
&                   d_delm,d_eluv,d_flux)
    CALL update_n1(nshape,nnod,nel1,nel1,nnod2,ielnd(1,1),ielsrt(1),    &
&                  zparallel,ndm0(1),ndm1(1),cutm(1),zactive(1),        &
&                  flux(1,1),tflux(1),ndu(1), &
&                  d_ndm0,d_ndm1,d_cutm,d_zactive,        &
&                  d_flux,d_tflux,d_ndu)

! gather here must happen before comms and can't reuse eluv
!    CALL gather(nshape,nel,nnod,ielnd(1,1),ndv(1),eluv(1,1))
    call op_par_loop_3(ale_advect_markactive2,s_nodes, &
&           op_arg_dat(d_indstatus, -1, OP_ID, 1, 'integer(4)', OP_READ), &
&           op_arg_dat(d_indtype,   -1, OP_ID, 1, 'integer(4)', OP_READ), &
&           op_arg_dat(d_zactive,   -1, OP_ID, 1, 'integer(4)', OP_WRITE))

    CALL flux_n1_VL(nshape,nel1,nel2,ielel(1,1),ielsd(1,1),cnm(1,1),    &
&                   delm(1,1),elvv(1,1),flux(1,1), &
&                   d_ielel,d_ielsd,d_cnm,    &
&                   d_delm,d_elvv,d_flux)
    CALL update_n1(nshape,nnod,nel1,nel1,nnod2,ielnd(1,1),ielsrt(1),    &
&                  zparallel,ndm0(1),ndm1(1),cutm(1),zactive(1),        &
&                  flux(1,1),tflux(1),ndv(1), &
&                  d_ndm0,d_ndm1,d_cutm,d_zactive,        &
&                  d_flux,d_tflux,d_ndv)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_update_nd_var=                               &
&    bookleaf_times%time_in_update_nd_var+t1

  END SUBROUTINE update_nd_var

END MODULE ale_advect_mod
