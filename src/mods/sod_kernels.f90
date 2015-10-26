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


MODULE sod_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE sod_midpoint(ndx,x1,x2)

    USE kinds_mod,ONLY: rlk,ink

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: ndx
    REAL(KIND=rlk), INTENT(INOUT) :: x1,x2

    IF (ndx.LT.x1) x1=ndx
    IF (ndx.GT.x2) x2=ndx

  END SUBROUTINE sod_midpoint

  SUBROUTINE sod_reset(ndx1,ndx2,ndx3,ndx4,xmid,ielmat,rho, &
&                      pre,ein,elmass,cnmass,eos_param,elvol,cnwt)

    USE kinds_mod,ONLY: rlk,ink
    USE parameters_mod,ONLY: LI, N_SHAPE

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: ndx1,ndx2,ndx3,ndx4,elvol
    REAL(KIND=rlk), INTENT(INOUT) :: xmid
    REAL(KIND=rlk), INTENT(OUT) :: rho,pre,ein,elmass
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: cnmass
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: cnwt
    REAL(KIND=rlk), DIMENSION(6,LI), INTENT(IN) :: eos_param
    INTEGER(KIND=ink), INTENT(OUT) :: ielmat


    IF ((0.25_rlk*(ndx1+ndx2+ndx3+ndx4)).LT.xmid) THEN
      ielmat=1_ink
      rho=1.0_rlk
      pre=1.0_rlk
    ELSE
      ielmat=2_ink
      rho=0.125_rlk
      pre=0.1_rlk
    ENDIF
    ein=pre/(rho*(eos_param(1,ielmat)-1.0_rlk))
    elmass=rho*elvol
    cnmass(1)=rho*cnwt(1)
    cnmass(2)=rho*cnwt(2)
    cnmass(3)=rho*cnwt(3)
    cnmass(4)=rho*cnwt(4)

  END SUBROUTINE sod_reset

  SUBROUTINE sod_subz(ndx1,ndx2,ndx3,ndx4,ndy1,ndy2,ndy3,ndy4,rho,spmass)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: ndx1,ndx2,ndx3,ndx4,ndy1,ndy2,ndy3,ndy4,rho
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: spmass

    REAL(KIND=rlk), DIMENSION(N_SHAPE) :: ndx,ndy
    REAL(KIND=rlk) :: x1,x2,x3,x4,y1,y2,y3,y4,w1,w2,w3,w4
    INTEGER(KIND=ink) :: ii,inod

    ndx(1) = ndx1
    ndx(2) = ndx2
    ndx(3) = ndx3
    ndx(4) = ndx4
    ndy(1) = ndy1
    ndy(2) = ndy2
    ndy(3) = ndy3
    ndy(4) = ndy4

    x3=0.25_rlk*(ndx(1)+ndx(2)+ndx(3)+ndx(4))
    y3=0.25_rlk*(ndy(1)+ndy(2)+ndy(3)+ndy(4))
    DO inod=1,4
      x1=ndx(inod)
      y1=ndy(inod)
      ii=MOD(inod,4)+1_ink
      x2=0.5_rlk*(x1+ndx(ii))
      y2=0.5_rlk*(y1+ndy(ii))
      ii=MOD(inod+2,4)+1_ink
      x4=0.5_rlk*(x1+ndx(ii))
      y4=0.5_rlk*(y1+ndy(ii))
      w1=-x1+x2+x3-x4
      w2=-x1-x2+x3+x4
      w3=-y1+y2+y3-y4
      w4=-y1-y2+y3+y4
      spmass(inod)=rho*(w1*w4-w2*w3)
    ENDDO

  END SUBROUTINE sod_subz

END MODULE sod_kernels
