
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

MODULE getq_mod

  IMPLICIT NONE

  PUBLIC :: getq

CONTAINS

  SUBROUTINE getq(d_elx,d_ely,d_elu,d_elv,d_rho,d_pre, &
&                 d_dx,d_dy,d_du,d_dv,d_scratch)

    USE kinds_mod,    ONLY: ink,rlk
    USE reals_mod,    ONLY: zerocut,cq1,cq2
    USE comms_mod,    ONLY: exchange,VISCOSITY
    USE paradef_mod,  ONLY: zparallel
    USE timing_mod,   ONLY: bookleaf_times
    USE typh_util_mod,ONLY: get_time
    USE integers_mod,ONLY: nshape
    USE OP2_Fortran_Reference
    use, intrinsic :: ISO_C_BINDING
    USE op2_bookleaf, ONLY: m_el2node,m_el2el,s_elements, &
&                           d_qq,d_qx,d_qy,d_csqrd,d_ielsd,d_indtype,d_ielel
    USE common_kernels, ONLY:set_zero1,set_zero4
    USE getq_kernels

    ! Argument list
    type(op_dat) d_elx,d_ely,d_elu,d_elv,d_rho,d_pre,d_dx,d_dy,d_du,d_dv, &
&                 d_scratch

    ! Local
    INTEGER(KIND=ink)                                  :: iel,iside
    REAL(KIND=rlk)                                     :: t0,t1
    ! Timer
    t0=get_time()

    !# Missing code here that can't be merged

!     ! initialisation
    call op_par_loop_1(set_zero1,s_elements, &
  &                    op_arg_dat(d_qq,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero4,s_elements, &
  &                    op_arg_dat(d_qx,-1,OP_ID,4,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero4,s_elements, &
  &                    op_arg_dat(d_qy,-1,OP_ID,4,'real(8)',OP_WRITE))


    ! gradient construction
    call op_par_loop_8(getq_gradcon,s_elements, &
  &                    op_arg_dat(d_du, -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_dv, -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_dx, -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_dy, -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_elu,-1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_elv,-1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_elx,-1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_ely,-1,OP_ID,4,'real(8)',OP_READ))


    ! MPI parallelism
    IF (zparallel) THEN
      CALL exchange(VISCOSITY)
    ENDIF

        DO iside=1,nshape/2_ink

    ! Christiensen monotonic limit
      call op_par_loop_23(getq_christiensen1,s_elements, &
  &                    op_arg_dat(d_du, -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dv, -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dx, -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dy, -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_du,  1,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_du,  2,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_du,  3,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_du,  4,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dv,  1,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dv,  2,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dv,  3,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dv,  4,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dx,  1,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dx,  2,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dx,  3,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dx,  4,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dy,  1,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dy,  2,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dy,  3,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_dy,  4,m_el2el,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_ielsd,-1,OP_ID,4,'integer(4)',OP_READ), &
  &                    op_arg_dat(d_scratch,-1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_gbl(iside,1,'integer(4)',OP_READ))

      ! BC
      call op_par_loop_7(getq_christiensen_bc,s_elements, &
  &                    op_arg_dat(d_ielel,   -1,OP_ID,4,'integer(4)',OP_READ), &
  &                    op_arg_dat(d_scratch, -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_indtype,  1,m_el2node,1,'integer(4)',OP_READ), &
  &                    op_arg_dat(d_indtype,  2,m_el2node,1,'integer(4)',OP_READ), &
  &                    op_arg_dat(d_indtype,  3,m_el2node,1,'integer(4)',OP_READ), &
  &                    op_arg_dat(d_indtype,  4,m_el2node,1,'integer(4)',OP_READ), &
  &                    op_arg_gbl(iside,1,'integer(4)',OP_READ))

      ! Apply limiter
      call op_par_loop_10(getq_christiensen_limiter,s_elements, &
  &                    op_arg_dat(d_du,   -1,OP_ID,4,'real(8)',OP_READ), & !du
  &                    op_arg_dat(d_dv,   -1,OP_ID,4,'real(8)',OP_READ), & !dv
  &                    op_arg_dat(d_qx,   -1,OP_ID,4,'real(8)',OP_WRITE), & !qx
  &                    op_arg_dat(d_qy,   -1,OP_ID,4,'real(8)',OP_WRITE), & !qy
  &                    op_arg_dat(d_csqrd,-1,OP_ID,1,'real(8)',OP_READ), &
  &                    op_arg_dat(d_rho,  -1,OP_ID,1,'real(8)',OP_READ), &
  &                    op_arg_dat(d_scratch, -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_gbl(cq1,1,'real(8)',OP_READ), &
  &                    op_arg_gbl(cq2,1,'real(8)',OP_READ), &
  &                    op_arg_gbl(iside,1,'integer(4)',OP_READ))
    ENDDO

    ! Final Q calculation
      call op_par_loop_7(getq_christiensen_q,s_elements, &
  &                    op_arg_dat(d_elx,   -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_ely,   -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_elu,   -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_elv,   -1,OP_ID,4,'real(8)',OP_READ), &
  &                    op_arg_dat(d_qx,    -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_qy,    -1,OP_ID,4,'real(8)',OP_WRITE), &
  &                    op_arg_dat(d_qq,    -1,OP_ID,1,'real(8)',OP_INC))

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_getq=bookleaf_times%time_in_getq+t1

  END SUBROUTINE getq

END MODULE getq_mod
