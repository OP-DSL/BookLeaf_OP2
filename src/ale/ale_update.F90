
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

MODULE ale_update_mod

#ifdef OP2_TRANSLATOR
  USE op2_kernels
#endif

  IMPLICIT NONE

  PUBLIC :: aleupdate

CONTAINS

  SUBROUTINE aleupdate(nshape,nel,nnod,d_ndx,d_ndy,d_elx,d_ely,d_elmass,d_rho,&
&                      d_pre,d_ein,d_csqrd)

    USE kinds_mod,    ONLY: ink,rlk
    USE pointers_mod, ONLY: elvol
    USE geometry_mod, ONLY: getgeom,getgeom2
    USE getpc_mod,    ONLY: getpc
    USE timing_mod,   ONLY: timer=>bookleaf_times,get_time
    USE op2_bookleaf, ONLY:s_elements,s_nodes,m_el2node,m_el2el,d_elidx,d_elvol
    USE op2_bookleaf_api
    USE ale_update_kernels

    ! Argument list
    INTEGER(KIND=ink),                      INTENT(IN)  :: nshape,nel,  &
&                                                           nnod
    type(op_dat) :: d_ndx,d_ndy,d_elx,d_ely,d_elmass,d_rho,&
&                      d_pre,d_ein,d_csqrd
    ! Local
    INTEGER(KIND=ink) :: iel
    REAL(KIND=rlk)    :: t0,t1

    ! Timer
    t0=get_time()

    ! update geometry
    CALL getgeom2(d_ndx,d_ndy,d_elx,d_ely,timer%time_in_getgeoma)

    ! update density to be consistent with geometry
    call op_par_loop_3(ale_update_rho, s_elements, &
&           op_arg_dat(d_rho,-1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_elmass,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_elvol,-1,OP_ID,1,'real(8)',OP_READ))
    ! update EoS
    CALL getpc(d_rho,d_ein,d_pre,d_csqrd,timer%time_in_getpca)

    ! Timing data
    t1=get_time()
    t1=t1-t0
    timer%time_in_aleupdate=timer%time_in_aleupdate+t1

  END SUBROUTINE aleupdate

END MODULE ale_update_mod
