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


MODULE common_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE gather_fun(a1,a2,a3,a4,b)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(IN) :: a1,a2,a3,a4
    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: b

    b(1) = a1
    b(2) = a2
    b(3) = a3
    b(4) = a4

  END SUBROUTINE gather_fun


  SUBROUTINE a_eq_b_over_c(a,b,c)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(OUT) :: a
    REAL(KIND=8), INTENT(IN) :: b,c

    a = b / c

  END SUBROUTINE a_eq_b_over_c

    SUBROUTINE a_eq_b(a,b)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(OUT) :: a
    REAL(KIND=8), INTENT(IN) :: b

    a = b

  END SUBROUTINE a_eq_b

      SUBROUTINE a_eq_b_4(a,b)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: a
    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: b

    a = b

  END SUBROUTINE a_eq_b_4

  SUBROUTINE a_eq_b_times_c(a,b,c)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(OUT) :: a
    REAL(KIND=8), INTENT(IN) :: b,c

    a = b*c

  END SUBROUTINE a_eq_b_times_c

  SUBROUTINE set_zero1(a)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(OUT) :: a

    a = 0.0_8

  END SUBROUTINE set_zero1

  SUBROUTINE set_zero4(a)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: a

!asdf

    a(1) = 0.0_8
    a(2) = 0.0_8
    a(3) = 0.0_8
    a(4) = 0.0_8

  END SUBROUTINE set_zero4

END MODULE common_kernels
