
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

MODULE getacc_mod

  IMPLICIT NONE

  PUBLIC :: getacc

CONTAINS

  SUBROUTINE getacc(dt05,dt,d_cnfx,d_cnfy,d_ndub,d_ndvb,d_elu,d_elv,d_rho)

    USE kinds_mod,    ONLY: ink,rlk
    USE reals_mod,    ONLY: zerocut,dencut,accut
    USE comms_mod,    ONLY: HALFSTEP,exchange
    USE paradef_mod,  ONLY: zparallel,ielsort1
    USE utilities_mod,ONLY: gather,gather2
    USE timing_mod,   ONLY: bookleaf_times
    USE TYPH_util_mod,ONLY: get_time
    USE OP2_Fortran_Reference
    use, intrinsic :: ISO_C_BINDING
    USE op2_bookleaf, ONLY: m_el2node,m_el2el,s_elements, s_nodes,d_ndmass,d_ndarea, &
&                           d_cnmass,d_cnwt,d_indtype,d_ndu,d_ndv,d_ndx,d_ndy
    USE getacc_kernels
    USE common_kernels, ONLY:set_zero1

    ! Argument list
    REAL(KIND=rlk),                      INTENT(IN)  :: dt05,dt
    type(op_dat), INTENT(INOUT) :: d_cnfx,d_cnfy,d_ndub,d_ndvb,d_elu,d_elv,d_rho
    ! Local
    INTEGER(KIND=ink)                                :: jj
    REAL(KIND=rlk)                                   :: t0,t1


    ! Timer
    t0=get_time()

    ! MPI parallelism
    IF (zparallel) THEN
      call exchange(HALFSTEP)
    ENDIF

    ! Construct nodal mass and scatter force to nodes
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndmass,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndarea,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndub,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_ndvb,-1,OP_ID,1,'real(8)',OP_WRITE))

    call op_par_loop_21(getacc_scatter,s_elements, &
&           op_arg_dat(d_cnmass,  -1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rho,     -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_cnwt,    -1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_cnfx,    -1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_cnfy,    -1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_ndmass,   1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndmass,   2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndmass,   3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndmass,   4,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndarea,   1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndarea,   2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndarea,   3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndarea,   4,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndub,     1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndub,     2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndub,     3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndub,     4,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndvb,     1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndvb,     2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndvb,     3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_ndvb,     4,m_el2node,1,'real(8)',OP_INC))

!This is silly, order of operations shouldn't matter, but it does, so in the meanwhile...
!For actual parallel execution it will be different anyway, will use the above
!    DO jj=1,4
!      call DISABLED(getacc_scatter2,s_elements, &
!&           op_arg_dat(d_cnmass,  -1,OP_ID,4,'real(8)',OP_READ), &
!&           op_arg_dat(d_rho,     -1,OP_ID,1,'real(8)',OP_READ), &
!&           op_arg_dat(d_cnwt,    -1,OP_ID,4,'real(8)',OP_READ), &
!&           op_arg_dat(d_cnfx,    -1,OP_ID,4,'real(8)',OP_READ), &
!&           op_arg_dat(d_cnfy,    -1,OP_ID,4,'real(8)',OP_READ), &
!&           op_arg_dat(d_ndmass,   jj,m_el2node,1,'real(8)',OP_INC), &
!&           op_arg_dat(d_ndarea,   jj,m_el2node,1,'real(8)',OP_INC), &
!&           op_arg_dat(d_ndub,   jj,m_el2node,1,'real(8)',OP_INC), &
!&           op_arg_dat(d_ndvb,   jj,m_el2node,1,'real(8)',OP_INC), &
!&           op_arg_gbl(jj,1,'integer(4)',OP_READ))
!    ENDDO


    ! Calculate acceleration

    call op_par_loop_4(getacc_accel,s_nodes, &
&           op_arg_dat(d_ndarea,  -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndmass,  -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndub,    -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndvb,    -1,OP_ID,1,'real(8)',OP_RW))


    !# Missing code here that can't be merged
    ! Boundary conditions
    call op_par_loop_3(getacc_bc,s_nodes, &
&           op_arg_dat(d_ndub,    -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndvb,    -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_indtype, -1,OP_ID,1,'integer(4)',OP_READ))

    !# Missing code here that can't be merged
    ! Calculate average velocity
    call op_par_loop_6(getacc_avgvel,s_nodes, &
&           op_arg_dat(d_ndub,    -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndvb,    -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndu,     -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndv,     -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_gbl(dt,  1,'real(8)',OP_READ), &
&           op_arg_gbl(dt05,1,'real(8)',OP_READ))

    CALL gather2(s_elements,m_el2node,d_ndub,d_elu)
    CALL gather2(s_elements,m_el2node,d_ndvb,d_elv)
    !# Missing code here that can't be merged
    ! Update position
    call op_par_loop_5(getacc_updpos,s_nodes, &
&           op_arg_dat(d_ndub,    -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndvb,    -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx,     -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_dat(d_ndy,     -1,OP_ID,1,'real(8)',OP_RW), &
&           op_arg_gbl(dt,  1,'real(8)',OP_READ))

    !# Missing code here that can't be merged

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_getacc=bookleaf_times%time_in_getacc+t1

  END SUBROUTINE getacc

END MODULE getacc_mod
