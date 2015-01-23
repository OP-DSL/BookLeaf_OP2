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


MODULE getforce_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getforce_pres(pre,elfx,elfy,a1,a3,b1,b3)

    USE kinds_mod,ONLY: rlk

    implicit none

    REAL(KIND=rlk), DIMENSION(4), INTENT(OUT) :: elfx,elfy
    REAL(KIND=rlk), INTENT(IN) :: pre,a1,a3,b1,b3

    REAL(KIND=rlk) :: w1

    w1=pre
    elfx(1)=w1*(-b3+b1)
    elfx(2)=w1*( b3+b1)
    elfx(3)=w1*( b3-b1)
    elfx(4)=w1*(-b3-b1)
    elfy(1)=w1*( a3-a1)
    elfy(2)=w1*(-a3-a1)
    elfy(3)=w1*(-a3+a1)
    elfy(4)=w1*( a3+a1)

  END SUBROUTINE getforce_pres

  SUBROUTINE getforce_visc(elfx,elfy,qx,qy)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(4), INTENT(INOUT) :: elfx,elfy
    REAL(KIND=rlk), DIMENSION(4), INTENT(IN) :: qx,qy

    INTEGER(KIND=ink) :: jj,jp

    DO jj=1,N_SHAPE
      jp=jj+1_ink
      IF (jp.GT.4_ink) jp=1_ink
      elfx(jj)=elfx(jj)+qx(jj)
      elfx(jp)=elfx(jp)-qx(jj)
      elfy(jj)=elfy(jj)+qy(jj)
      elfy(jp)=elfy(jp)-qy(jj)
    ENDDO

  END SUBROUTINE getforce_visc

END MODULE getforce_kernels