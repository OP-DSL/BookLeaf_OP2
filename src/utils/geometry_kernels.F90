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

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: elx,ely
    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: cnwt
    REAL(KIND=8), INTENT(OUT) :: a1,a2,a3,b1,b2,b3,elvol

    a1=0.25_8*(-elx(1)+elx(2)+elx(3)-elx(4))
    a2=0.25_8*( elx(1)-elx(2)+elx(3)-elx(4))
    a3=0.25_8*(-elx(1)-elx(2)+elx(3)+elx(4))
    b1=0.25_8*(-ely(1)+ely(2)+ely(3)-ely(4))
    b2=0.25_8*( ely(1)-ely(2)+ely(3)-ely(4))
    b3=0.25_8*(-ely(1)-ely(2)+ely(3)+ely(4))
    cnwt(1)=1.0_8/9.0_8 *                            &
&               ((3.0_8*b3-b2)*(3.0_8*a1-a2)  &
&               -(3.0_8*a3-a2)*(3.0_8*b1-b2))
    cnwt(2)=1.0_8/9.0_8 *                            &
&               ((3.0_8*b3+b2)*(3.0_8*a1-a2)  &
&               -(3.0_8*a3+a2)*(3.0_8*b1-b2))
    cnwt(3)=1.0_8/9.0_8 *                            &
&               ((3.0_8*b3+b2)*(3.0_8*a1+a2)  &
                -(3.0_8*a3+a2)*(3.0_8*b1+b2))
    cnwt(4)=1.0_8/9.0_8 *                            &
&               ((3.0_8*b3-b2)*(3.0_8*a1+a2)  &
                -(3.0_8*a3-a2)*(3.0_8*b1+b2))
    elvol=4.0_8*(a1*b3-a3*b1)

  END SUBROUTINE geometry_calc

  SUBROUTINE geometry_min(elvol,minval)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    INTEGER(KIND=4), INTENT(INOUT) :: minval
    REAL(KIND=8), INTENT(IN) :: elvol

    IF (elvol.LT.0.0_8) THEN
      minval = MIN(minval,0)
    ELSE
      minval = MIN(minval,1)
    ENDIF

  END SUBROUTINE geometry_min

END MODULE geometry_kernels