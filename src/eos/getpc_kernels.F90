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


MODULE getpc_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getpc_update(im,eos_type,eos_param,rho,ein,pre,csqrd)

    USE kinds_mod,ONLY: ink,rlk
    USE bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(OUT) :: pre,csqrd
    REAL(KIND=8), INTENT(IN) :: rho,ein
    INTEGER(KIND=4), INTENT(IN) :: im
    REAL(KIND=8),DIMENSION(6,100),INTENT(IN) :: eos_param
    INTEGER(KIND=4), DIMENSION(100), INTENT(IN) :: eos_type

    REAL(KIND=8) :: t1,t2,t3,t4,t5

    ! getpre
    IF (eos_type(im).EQ.0) THEN ! VOID
      pre=eos_param(1,im)
    ELSEIF (eos_type(im).EQ.1) THEN ! IDEAL GAS
      pre=ein*rho*(eos_param(1,im)-1.0_8)
    ELSEIF (eos_type(im).EQ.2) THEN ! TAIT
      t1=rho/eos_param(3,im)
      pre=eos_param(1,im)*(t1**eos_param(2,im)-1.0_8)
      pre=MAX(pre,eos_param(4,im))
    ELSEIF (eos_type(im).EQ.3) THEN ! JWL
      t1=eos_param(4,im)*eos_param(6,im)/rho
      t2=eos_param(5,im)*eos_param(6,im)/rho
      t3=eos_param(1,im)*rho*ein
      t4=(1.0_8-eos_param(1,im)/t1)*eos_param(2,im)*EXP(-t1)
      t5=(1.0_8-eos_param(1,im)/t2)*eos_param(3,im)*EXP(-t2)
      pre=t3+t4+t5
    ELSE
      pre=-1.0_8
    ENDIF
    IF (ABS(pre).LT.pcut) pre=0.0_8

    ! getcc
    IF (eos_type(im).EQ.0) THEN ! VOID
      csqrd=1.0D-6
    ELSEIF (eos_type(im).EQ.1) THEN ! IDEAL GAS
      csqrd=eos_param(1,im)*(eos_param(1,im)-1.0_8)*ein
    ELSEIF (eos_type(im).EQ.2) THEN ! TAIT
      t1=rho/eos_param(3,im)
      t2=eos_param(2,im)-1.0_8
      csqrd=(eos_param(1,im)*eos_param(2,im))/eos_param(3,im)
      csqrd=csqrd*t1**t2
    ELSEIF (eos_type(im).EQ.3) THEN ! JWL
      t1=eos_param(6,im)/rho
      t2=pre !getpre(im,rho,ein)
      t3=eos_param(4,im)*t1
      t4=eos_param(1,im)/eos_param(4,im)+eos_param(1,im)*t1-t3*t1
      t4=t4*eos_param(2,im)*EXP(-t3)
      t3=eos_param(5,im)*t1
      t5=eos_param(1,im)/eos_param(5,im)+eos_param(1,im)*t1-t3*t1
      t5=t5*eos_param(3,im)*EXP(-t3)
      csqrd=eos_param(1,im)*t2/rho+eos_param(1,im)*ein-t4-t5
    ELSE
      csqrd=-1.0_8
    ENDIF

  END SUBROUTINE getpc_update

  SUBROUTINE getpc_merge(csqrd)

    USE kinds_mod,ONLY: ink,rlk
    USE integers_mod,ONLY: nshape

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: csqrd

    IF (csqrd.LT.0.0_8) THEN
      csqrd = 0.0_8
    ENDIF

  END SUBROUTINE getpc_merge

END MODULE getpc_kernels
