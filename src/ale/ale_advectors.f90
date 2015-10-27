
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

MODULE ale_advectors_mod

  USE kinds_mod,ONLY: ink,rlk,lok

  IMPLICIT NONE

  PUBLIC :: flux_c1_VL,flux_n1_VL,sum_flux,update_c1,update_n1

CONTAINS

  SUBROUTINE flux_c1_VL(iD1,iD2,iShape,iLSize,iASize,iElEl,iElSd,       &
&                       rCorner,rDel,rVar,rFlux, &
&                       d_ielel,d_ielsd,d_rCorner,d_rDel,d_rVar,d_rFlux)
      USE op2_bookleaf
      USE ale_advectors_kernels
      USE common_kernels
    ! Argument list
    INTEGER(KIND=ink),                         INTENT(IN)  :: iD1,iD2,  &
&                                                             iShape,   &
&                                                             iLSize,   &
&                                                             iASize
    INTEGER(KIND=ink),DIMENSION(iShape,iASize),INTENT(IN)  :: iElEl,    &
&                                                             iElSd
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(IN)  :: rCorner,  &
&                                                             rDel
    REAL(KIND=rlk),   DIMENSION(iASize),       INTENT(IN)  :: rVar
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(OUT) :: rFlux
    type(op_dat) :: d_ielel,d_ielsd,d_rCorner,d_rDel,d_rVar,d_rFlux

    ! Local
    INTEGER(KIND=ink) :: i1,i2,j1,j2,iEl,iE1,iE2
    REAL(KIND=rlk)    :: r1,r2,r3,r4,w1,w2,w3,w4,w5,w6,w7,w8,rV,rGrad

    ! initialise
    call op_par_loop_1(set_zero4,s_elements, &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_WRITE))
!    rFlux=0.0_rlk

    ! construct flux
    DO i1=iD1,iD2
      call op_par_loop_14(ale_advectors_flux,s_elements, &
&           op_arg_dat(d_rCorner,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 1,m_el2el,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 2,m_el2el,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 3,m_el2el,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 4,m_el2el,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rDel,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_dat(d_ielsd,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_gbl(i1,1,'integer(4)',OP_READ))
    ENDDO

  END SUBROUTINE flux_c1_VL

  SUBROUTINE flux_n1_VL(iShape,iLSize,iASize,iElEl,iElSd,rCorner,rDel,  &
&                       rVar,rFlux, &
&                       d_iElEl,d_iElSd,d_rCorner,d_rDel,  &
&                       d_rVar,d_rFlux)

  USE op2_bookleaf
  USE common_kernels
  USE ale_advectors_kernels

    ! Argument list
    INTEGER(KIND=ink),                         INTENT(IN)  :: iShape,   &
&                                                             iLSize,   &
&                                                             iASize
    INTEGER(KIND=ink),DIMENSION(iShape,iASize),INTENT(IN)  :: iElEl,    &
&                                                             iElSd
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(IN)  :: rCorner,  &
&                                                             rDel,rVar
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(OUT) :: rFlux
    type(op_dat) ::                     d_iElEl,d_iElSd,d_rCorner,d_rDel,  &
&                       d_rVar,d_rFlux
    ! Local
    INTEGER(KIND=ink) :: iFaceL,iFaceR,iCorner,iEl,iElL,iElR,iSdL,iSdR, &
&                        iLNdL,iLNdR,iLNNdL,iLNNdR,ii
    REAL(KIND=rlk)    :: w1,w2,w3,w4,w5,w6,w7,w8,rD,rGrad,rV

    ! initialise
!    rFlux=0.0_rlk
        call op_par_loop_1(set_zero4,s_elements, &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_WRITE))

    ! construct flux
    DO iFaceL=1,2
      DO iCorner=1,2
          call op_par_loop_15(ale_advectors_flux2,s_elements, &
&           op_arg_dat(d_rCorner,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rCorner, 4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar, 4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rDel,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_WRITE), &
&           op_arg_dat(d_ielsd,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_gbl(iFaceL,1,'integer(4)',OP_READ), &
&           op_arg_gbl(iCorner,1,'integer(4)',OP_READ))
      ENDDO
    ENDDO

  END SUBROUTINE flux_n1_VL

  SUBROUTINE update_c1(iD1,iD2,iShape,iLSize,iASize,iElEl,iElSd,rBase0, &
&                      rBase1,rCutOff,rFlux,rTotFlux,rVar, &
&                      d_ielel,d_ielsd,d_rBase0,d_rBase1,d_rCutOff,d_rFlux, &
&                      d_rTotFlux, d_rVar)
      USE op2_bookleaf
      USE ale_advectors_kernels
    ! Argument list
    INTEGER(KIND=ink),                         INTENT(IN)    :: iD1,iD2,&
&                                                               iShape, &
&                                                               iLSize, &
&                                                               iASize
    INTEGER(KIND=ink),DIMENSION(iShape,iASize),INTENT(IN)    :: iElEl,  &
&                                                               iElSd
    REAL(KIND=rlk),   DIMENSION(iASize),       INTENT(IN)    :: rBase0, &
&                                                               rBase1, &
&                                                               rCutOff
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(IN)    :: rFlux
    REAL(KIND=rlk),   DIMENSION(iASize),       INTENT(OUT)   :: rTotFlux
    REAL(KIND=rlk),   DIMENSION(iASize),       INTENT(INOUT) :: rVar
    type(op_dat) ::    d_ielel,d_ielsd,d_rBase0,d_rBase1,d_rCutOff,d_rFlux, &
&                      d_rTotFlux, d_rVar
    ! Local
    INTEGER(KIND=ink) :: iEl

    ! calculate total flux
    CALL sum_flux(iD1,iD2,iShape,iLSize,iASize,iElEl,iElSd,rFlux,       &
&                 rTotFlux, &
&                  d_iElEl,d_iElSd,d_rFlux,       &
&                 d_rTotFlux)

    ! update variable
    !Is this on elements?
    call op_par_loop_5(ale_advectors_update_c1,s_elements, &
&           op_arg_dat(d_rBase0, -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rBase1, -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rCutOff,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rTotFlux,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar,-1,OP_ID,1,'real(8)',OP_WRITE))

  END SUBROUTINE update_c1

  SUBROUTINE update_n1(iShape,iUSize,iESize,iCSize,iNSize,iElNd,iElSrt, &
&                      zparallel,rBase0,rBase1,rCut,zActive,rFlux,      &
&                      rTotFlux,rVar, &
&                     d_rBase0,d_rBase1,d_rCut,d_zActive,d_rFlux,      &
&                      d_rTotFlux_bad,d_rVar)
    USE op2_bookleaf
    USE common_kernels
    USE ale_advectors_kernels
    ! Argument list
    INTEGER(KIND=ink),                         INTENT(IN)    :: iShape, &
&                                                               iUSize, &
&                                                               iESize, &
&                                                               iCSize, &
&                                                               iNSize
    INTEGER(KIND=ink),DIMENSION(iShape,iCSize),INTENT(IN)    :: iElNd
    INTEGER(KIND=ink),DIMENSION(iCSize),       INTENT(IN)    :: iElSrt
    LOGICAL(KIND=lok),                         INTENT(IN)    :: zparallel
    REAL(KIND=rlk),   DIMENSION(iNSize),       INTENT(IN)    :: rBase0, &
&                                                               rBase1, &
&                                                               rCut
    !LOGICAL(KIND=lok),DIMENSION(iNSize),       INTENT(IN)    :: zActive
    INTEGER(KIND=ink),DIMENSION(iNSize),      INTENT(OUT)  :: zActive
    REAL(KIND=rlk),   DIMENSION(iShape,iCSize),INTENT(IN)    :: rFlux
    REAL(KIND=rlk),   DIMENSION(iNSize),       INTENT(OUT)   :: rTotFlux
    REAL(KIND=rlk),   DIMENSION(iNSize),       INTENT(INOUT) :: rVar
    type(op_dat) ::    d_rBase0,d_rBase1,d_rCut,d_zActive,d_rFlux,      &
&                      d_rTotFlux_bad,d_rVar
    ! Local
    INTEGER(KIND=ink) :: iEl,iNd,ii,jj

    ! construct total flux
    call op_par_loop_1(set_zero1,s_nodes, &
&           op_arg_dat(d_rscratch18,-1,OP_ID,1,'real(8)',OP_WRITE))
    call op_par_loop_5(ale_advectors_totflux,s_elements, &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rscratch18,1,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_rscratch18,2,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_rscratch18,3,m_el2node,1,'real(8)',OP_INC), &
&           op_arg_dat(d_rscratch18,4,m_el2node,1,'real(8)',OP_INC))

    ! update variable
    call op_par_loop_6(ale_advectors_update_n1,s_nodes, &
&           op_arg_dat(d_rBase0, -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rBase1, -1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rCut,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rscratch18,-1,OP_ID,1,'real(8)',OP_READ), &
&           op_arg_dat(d_rVar,-1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_zActive,-1,OP_ID,1,'integer(4)',OP_READ))

  END SUBROUTINE update_n1

  SUBROUTINE sum_flux(iD1,iD2,iShape,iLSize,iASize,iElEl,iElSd,  &
&                     rFlux,rTotFlux, &
&                     d_ielel,d_ielsd,d_rFlux,d_rTotFlux)
    USE op2_bookleaf
    USE common_kernels
    USE ale_advectors_kernels

    ! Argument list
    INTEGER(KIND=ink),                         INTENT(IN)    :: iD1,iD2,&
&                                                               iShape, &
&                                                               iLSize, &
&                                                               iASize
    INTEGER(KIND=ink),DIMENSION(iShape,iASize),INTENT(IN)    :: iElEl,  &
&                                                               iElSd
    REAL(KIND=rlk),   DIMENSION(iShape,iASize),INTENT(IN)    :: rFlux
    REAL(KIND=rlk),   DIMENSION(iASize),       INTENT(OUT)   :: rTotFlux
    type(op_dat) :: d_ielel,d_ielsd,d_rFlux,d_rTotFlux
    ! Local
    INTEGER(KIND=ink) :: i1,i2,j1,j2,iEl,iE1,iE2,kk
    REAL(KIND=rlk)    :: w1,w2

    call op_par_loop_1(set_zero1,s_elements, &
&           op_arg_dat(d_rTotFlux,-1,OP_ID,1,'real(8)',OP_WRITE))

    DO i1=iD1,iD2
    call op_par_loop_10(ale_advectors_sumflux,s_elements, &
&           op_arg_dat(d_rFlux,-1,OP_ID,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,1,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,2,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,3,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rFlux,4,m_el2el,4,'real(8)',OP_READ), &
&           op_arg_dat(d_rTotFlux,-1,OP_ID,1,'real(8)',OP_WRITE), &
&           op_arg_dat(d_ielsd,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_dat(d_ielel,-1,OP_ID,4,'integer(4)',OP_READ), &
&           op_arg_dat(d_elidx,-1,OP_ID,1,'integer(4)',OP_READ), &
&           op_arg_gbl(i1,1,'integer(4)',OP_READ))
    ENDDO

  END SUBROUTINE sum_flux

END MODULE ale_advectors_mod
