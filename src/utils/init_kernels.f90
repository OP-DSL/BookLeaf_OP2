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

MODULE init_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE init_dem(im,mat_rho,mat_ein,rho,ein,elmass,elvol,cnmass,cnwt)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: LI, N_SHAPE

    implicit none

    REAL(KIND=rlk), INTENT(OUT) :: rho,ein,elmass
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: cnmass
    REAL(KIND=rlk), DIMENSION(LI), INTENT(IN) :: mat_rho,mat_ein
    REAL(KIND=rlk), INTENT(IN) :: elvol
    INTEGER(KIND=ink), INTENT(IN) :: im
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: cnwt

    rho=mat_rho(im)
    ein=mat_ein(im)
    elmass=rho*elvol
    cnmass(1) = rho*cnwt(1)
    cnmass(2) = rho*cnwt(2)
    cnmass(3) = rho*cnwt(3)
    cnmass(4) = rho*cnwt(4)


  END SUBROUTINE init_dem


  SUBROUTINE init_subz_pm(elx,ely,rho,spmass)

    USE kinds_mod,ONLY: rlk,ink
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: spmass
    REAL(KIND=rlk), INTENT(IN) :: rho
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: elx,ely

    REAL(KIND=rlk) :: x1,x2,x3,x4,y1,y2,y3,y4,w1,w2,w3,w4
    INTEGER(KIND=ink) :: j1,j2

    x3=0.25_rlk*(elx(1)+elx(2)+elx(3)+elx(4))
    y3=0.25_rlk*(ely(1)+ely(2)+ely(3)+ely(4))
    DO j1=1,N_SHAPE
      x1=elx(j1)
      y1=ely(j1)
      j2=MOD(j1,N_SHAPE)+1_ink
      x2=0.5_rlk*(x1+elx(j2))
      y2=0.5_rlk*(y1+ely(j2))
      j2=MOD(j1+2,N_SHAPE)+1_ink
      x4=0.5_rlk*(x1+elx(j2))
      y4=0.5_rlk*(y1+ely(j2))
      !# Axi-symmetric alternative
      w1=0.25_rlk*(-x1+x2+x3-x4)
      w2=0.25_rlk*(-x1-x2+x3+x4)
      w3=0.25_rlk*(-y1+y2+y3-y4)
      w4=0.25_rlk*(-y1-y2+y3+y4)
      spmass(j1)=4.0_rlk*rho*(w1*w4-w2*w3)
    ENDDO


  END SUBROUTINE init_subz_pm

END MODULE init_kernels
