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


MODULE geometry_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE geometry_calc(a1,a2,a3,b1,b2,b3,elx,ely,cnwt,elvol)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: ONEBYNINE

    implicit none

    REAL(KIND=rlk), DIMENSION(4), INTENT(IN) :: elx,ely
    REAL(KIND=rlk), DIMENSION(4), INTENT(OUT) :: cnwt
    REAL(KIND=rlk), INTENT(OUT) :: a1,a2,a3,b1,b2,b3,elvol

    a1=0.25_rlk*(-elx(1)+elx(2)+elx(3)-elx(4))
    a2=0.25_rlk*( elx(1)-elx(2)+elx(3)-elx(4))
    a3=0.25_rlk*(-elx(1)-elx(2)+elx(3)+elx(4))
    b1=0.25_rlk*(-ely(1)+ely(2)+ely(3)-ely(4))
    b2=0.25_rlk*( ely(1)-ely(2)+ely(3)-ely(4))
    b3=0.25_rlk*(-ely(1)-ely(2)+ely(3)+ely(4))
    cnwt(1)=ONEBYNINE*                            &
&               ((3.0_rlk*b3-b2)*(3.0_rlk*a1-a2)  &
&               -(3.0_rlk*a3-a2)*(3.0_rlk*b1-b2))
    cnwt(2)=ONEBYNINE*                            &
&               ((3.0_rlk*b3+b2)*(3.0_rlk*a1-a2)  &
&               -(3.0_rlk*a3+a2)*(3.0_rlk*b1-b2))
    cnwt(3)=ONEBYNINE*                            &
&               ((3.0_rlk*b3+b2)*(3.0_rlk*a1+a2)  &
                -(3.0_rlk*a3+a2)*(3.0_rlk*b1+b2))
    cnwt(4)=ONEBYNINE*                            &
&               ((3.0_rlk*b3-b2)*(3.0_rlk*a1+a2)  &
                -(3.0_rlk*a3-a2)*(3.0_rlk*b1+b2))
    elvol=4.0_rlk*(a1*b3-a3*b1)

  END SUBROUTINE geometry_calc

  SUBROUTINE geometry_min(elvol,minval)

    USE kinds_mod,ONLY: rlk

    implicit none

    INTEGER(KIND=ink), INTENT(INOUT) :: minval
    REAL(KIND=rlk), INTENT(IN) :: elvol

    IF (elvol.LT.0.0_rlk) THEN
      minval = MIN(minval,0_ink)
    ELSE
      minval = MIN(minval,1_ink)
    ENDIF

  END SUBROUTINE geometry_min

END MODULE geometry_kernels