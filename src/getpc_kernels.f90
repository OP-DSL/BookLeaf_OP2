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
    USE reals_mod,ONLY: pcut
    use parameters_mod, ONLY: LI

    implicit none

    REAL(KIND=rlk), INTENT(OUT) :: pre,csqrd
    REAL(KIND=rlk), INTENT(IN) :: rho,ein
    INTEGER(kind=ink), INTENT(IN) :: im
    REAL(KIND=rlk),DIMENSION(6,LI),INTENT(IN) :: eos_param
    INTEGER(KIND=ink), DIMENSION(LI), INTENT(IN) :: eos_type

    REAL(KIND=rlk) :: t1,t2,t3,t4,t5

    !getpre
    SELECT CASE(eos_type(im))
      CASE(0_ink) ! VOID
        pre=eos_param(1,im)
      CASE(1_ink) ! IDEAL GAS
        pre=ein*rho*(eos_param(1,im)-1.0_rlk)
      CASE(2_ink) ! TAIT
        t1=rho/eos_param(3,im)
        pre=eos_param(1,im)*(t1**eos_param(2,im)-1.0_rlk)
        pre=MAX(pre,eos_param(4,im))
      CASE(3_ink) ! JWL
        t1=eos_param(4,im)*eos_param(6,im)/rho
        t2=eos_param(5,im)*eos_param(6,im)/rho
        t3=eos_param(1,im)*rho*ein
        t4=(1.0_rlk-eos_param(1,im)/t1)*eos_param(2,im)*EXP(-t1)
        t5=(1.0_rlk-eos_param(1,im)/t2)*eos_param(3,im)*EXP(-t2)
        pre=t3+t4+t5
      CASE DEFAULT
        pre=-1_rlk
    END SELECT
    IF (ABS(pre).LT.pcut) pre=0.0_rlk

    !getcc
    SELECT CASE(eos_type(im))
      CASE(0_ink) ! VOID
        csqrd=1.0e-6_rlk
      CASE(1_ink) ! IDEAL GAS
        csqrd=eos_param(1,im)*(eos_param(1,im)-1.0_rlk)*ein
      CASE(2_ink) ! TAIT
        t1=rho/eos_param(3,im)
        t2=eos_param(2,im)-1.0_rlk
        csqrd=(eos_param(1,im)*eos_param(2,im))/eos_param(3,im)
        csqrd=csqrd*t1**t2
      CASE(3_ink) ! JWL
        t1=eos_param(6,im)/rho
        t2=pre !getpre(im,rho,ein)
        t3=eos_param(4,im)*t1
        t4=eos_param(1,im)/eos_param(4,im)+eos_param(1,im)*t1-t3*t1
        t4=t4*eos_param(2,im)*EXP(-t3)
        t3=eos_param(5,im)*t1
        t5=eos_param(1,im)/eos_param(5,im)+eos_param(1,im)*t1-t3*t1
        t5=t5*eos_param(3,im)*EXP(-t3)
        csqrd=eos_param(1,im)*t2/rho+eos_param(1,im)*ein-t4-t5
      CASE DEFAULT
        csqrd=-1_rlk
    END SELECT

  END SUBROUTINE getpc_update

  SUBROUTINE getpc_merge(csqrd)

    USE kinds_mod,ONLY: rlk
    USE integers_mod,ONLY: nshape

    implicit none

    REAL(KIND=rlk), INTENT(INOUT) :: csqrd

    IF (csqrd.LT.0.0_rlk) THEN
      csqrd = 0.0_rlk
    ENDIF

  END SUBROUTINE getpc_merge

END MODULE getpc_kernels