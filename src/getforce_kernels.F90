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

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: elfx,elfy
    REAL(KIND=8), INTENT(IN) :: pre,a1,a3,b1,b3

    REAL(KIND=8) :: w1

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

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(INOUT) :: elfx,elfy
    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: qx,qy

    INTEGER(KIND=4) :: jj,jp

    DO jj=1,4
      jp=jj+1
      IF (jp.GT.4) jp=1
      elfx(jj)=elfx(jj)+qx(jj)
      elfx(jp)=elfx(jp)-qx(jj)
      elfy(jj)=elfy(jj)+qy(jj)
      elfy(jp)=elfy(jp)-qy(jj)
    ENDDO

  END SUBROUTINE getforce_visc

END MODULE getforce_kernels