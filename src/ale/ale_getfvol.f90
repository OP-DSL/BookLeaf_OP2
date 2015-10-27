
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

MODULE ale_getfvol_mod

  USE kinds_mod,ONLY: ink,lok,rlk
!  USE op2_bookleaf, ONLY: op_dat,s_elements,s_nodes,m_el2node
  IMPLICIT NONE

  PRIVATE :: fvol
  PUBLIC  :: alegetfvol

CONTAINS

  SUBROUTINE alegetfvol(nshape,nnod,nel,nel2,dt,cut,indstatus,ielnd,    &
&                       ndx,ndy,ndux,ndvy,rdelv, &
&                       d_indstatus,d_ndx,d_ndy,d_ndux,d_ndvy,d_rdelv)

    USE kinds_mod,    ONLY: ink,rlk
    USE logicals_mod, ONLY: zeul
    USE timers_mod,   ONLY: bookleaf_times,get_time
    USE op2_bookleaf
    USE common_kernels
    USE ale_getfvol_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN)   :: nshape,nel,&
&                                                            nnod,nel2
    REAL(KIND=rlk),                          INTENT(IN)   :: dt,cut
    INTEGER(KIND=ink),DIMENSION(nnod),       INTENT(IN)   :: indstatus
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN)   :: ielnd
    REAL(KIND=rlk),   DIMENSION(nnod),       INTENT(INOUT):: ndx,ndy,   &
&                                                            ndux,ndvy
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT)  :: rdelv
    type(op_dat) :: d_indstatus,d_ndx,d_ndy,d_ndux,d_ndvy,d_rdelv
    ! Local
    INTEGER(KIND=ink) :: iNd
    REAL(KIND=rlk)    :: t0,t1,m1


    ! Timer
    t0=get_time()


    ! calculate mesh velocity
    IF (zeul) THEN
      m1 = -1_rlk
      call op_par_loop_3(a_eq_b_times_c, s_nodes, &
&             op_arg_dat(d_ndux,-1,OP_ID,1,'real(8)',OP_WRITE), &
&             op_arg_dat(d_ndux,-1,OP_ID,1,'real(8)',OP_READ), &
&             op_arg_gbl(m1,1,'real(8)',OP_READ))
      call op_par_loop_3(a_eq_b_times_c, s_nodes, &
&             op_arg_dat(d_ndvy,-1,OP_ID,1,'real(8)',OP_WRITE), &
&             op_arg_dat(d_ndvy,-1,OP_ID,1,'real(8)',OP_READ), &
&             op_arg_gbl(m1,1,'real(8)',OP_READ))
!       ndux=-ndux
!       ndvy=-ndvy
    ELSE
      !# Exchange MESH_MOTION ndx,ndy, Ndux,ndvy to nnod1
      ! Other options
    ENDIF

    ! construct new position
    call op_par_loop_5(ale_getfvol_newpos, s_nodes, &
&             op_arg_dat(d_ndux,-1,OP_ID,1,'real(8)',OP_RW), &
&             op_arg_dat(d_ndvy,-1,OP_ID,1,'real(8)',OP_RW), &
&             op_arg_dat(d_ndx,-1,OP_ID,1,'real(8)',OP_READ), &
&             op_arg_dat(d_ndy,-1,OP_ID,1,'real(8)',OP_READ), &
&             op_arg_gbl(dt,1,'real(8)',OP_READ))
!     DO iNd=1,nNod
!       ndux(iNd)=ndx(iNd)+dt*ndux(iNd)
!       ndvy(iNd)=ndy(iNd)+dt*ndvy(iNd)
!     ENDDO

    ! construct flux volumes
    CALL fvol(nshape,nnod,nel,nel2,cut,ielnd(1,1),ndx(1),ndy(1),ndux(1),&
&             ndvy(1),rDelV(1,1), &
&             d_ndx,d_ndy,d_ndux, &
&             d_ndvy,d_rDelV )

    ! update position
    call op_par_loop_2(a_eq_b, s_nodes, &
&             op_arg_dat(d_ndx,-1,OP_ID,1,'real(8)',OP_WRITE), &
&             op_arg_dat(d_ndux,-1,OP_ID,1,'real(8)',OP_READ))
    call op_par_loop_2(a_eq_b, s_nodes, &
&             op_arg_dat(d_ndy,-1,OP_ID,1,'real(8)',OP_WRITE), &
&             op_arg_dat(d_ndvy,-1,OP_ID,1,'real(8)',OP_READ))
!     DO iNd=1,nNod
!       ndx(iNd)=ndux(iNd)
!       ndy(iNd)=ndvy(iNd)
!     ENDDO

    ! Timing data
    t1=get_time()
    t1=t1-t0
    bookleaf_times%time_in_alegetfvol=bookleaf_times%time_in_alegetfvol+&
&                                     t1

  END SUBROUTINE alegetfvol

  SUBROUTINE fvol(nshape,nnod,nel,nel2,cut,ielnd,ndx0,ndy0,ndx1,ndy1,   &
&                 rdelv, &
&                d_ndx0,d_ndy0,d_ndx1,d_ndy1,d_rdelv)

    USE kinds_mod,ONLY: ink,rlk
    USE op2_bookleaf
    USE common_kernels
    USE ale_getfvol_kernels

    ! Argument list
    INTEGER(KIND=ink),                       INTENT(IN) :: nshape,nnod, &
&                                                          nel,nel2
    REAL(KIND=rlk),                          INTENT(IN) :: cut
    INTEGER(KIND=ink),DIMENSION(nshape,nel2),INTENT(IN) :: ielnd
    REAL(KIND=rlk),   DIMENSION(nnod),       INTENT(IN) :: ndx0,ndy0,   &
&                                                          ndx1,ndy1
    REAL(KIND=rlk),   DIMENSION(nshape,nel2),INTENT(OUT):: rdelv
    type(op_dat) :: d_ndx0,d_ndy0,d_ndx1,d_ndy1,d_rdelv

    ! Local
    INTEGER(KIND=ink) :: iel,ii,jj,jp,n1,n2
    REAL(KIND=rlk)    :: x1,x2,x3,x4,y1,y2,y3,y4,a1,a3,b1,b3

    ! initialise
!     rdelv=0.0_rlk
    call op_par_loop_1(set_zero4,s_elements, &
&           op_arg_dat(d_rdelv,-1,OP_ID,4,'real(8)',OP_WRITE))

    ! construct volumes
    call op_par_loop_18(ale_getfvol_vol,s_elements, &
&           op_arg_dat(d_ndx0,1,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx0,2,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx0,3,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx0,4,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy0,1,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy0,2,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy0,3,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy0,4,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx1,1,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx1,2,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx1,3,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndx1,4,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy1,1,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy1,2,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy1,3,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_ndy1,4,m_el2node,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rdelv,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_gbl(cut,1,'real(8)',OP_READ))

!     DO iel=1,nel
!       DO jj=1,nshape
!         jp=MOD(jj,nshape)+1_ink
!         n1=ielnd(jj,iel)
!         n2=ielnd(jp,iel)
!         x1=ndx0(n1)
!         x2=ndx0(n2)
!         y1=ndy0(n1)
!         y2=ndy0(n2)
!         x3=ndx1(n2)
!         x4=ndx1(n1)
!         y3=ndy1(n2)
!         y4=ndy1(n1)
!         a1=0.25_rlk*(-x1+x2+x3-x4)
!         a3=0.25_rlk*(-x1-x2+x3+x4)
!         b1=0.25_rlk*(-y1+y2+y3-y4)
!         b3=0.25_rlk*(-y1-y2+y3+y4)
!         rdelv(jj,iel)=4.0_rlk*(a1*b3-a3*b1)
!         IF (rdelv(jj,iel).LT.cut) rdelv(jj,iel)=0.0_rlk
!       ENDDO
!     ENDDO

    END SUBROUTINE fvol

END MODULE ale_getfvol_mod
