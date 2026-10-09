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

MODULE ale_getfvol_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE ale_getfvol_newpos(ndux,ndvy,ndx,ndy,dt)
    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(IN) :: ndx,ndy
    REAL(KIND=8), INTENT(INOUT) :: ndux,ndvy
    REAL(KIND=8) :: dt

    ndux=ndx+dt*ndux
    ndvy=ndy+dt*ndvy


  END SUBROUTINE ale_getfvol_newpos

  SUBROUTINE ale_getfvol_min1(ndux,ndvy)
    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: ndux,ndvy

    ndux=-1.0_8*ndux
    ndvy=-1.0_8*ndvy


  END SUBROUTINE ale_getfvol_min1

  SUBROUTINE ale_getfvol_vol(ndx01,ndx02,ndx03,ndx04, &
    & ndy01,ndy02,ndy03,ndy04,ndx11,ndx12,ndx13,ndx14, &
    & ndy11,ndy12,ndy13,ndy14,rdelv,cut)
    USE kinds_mod,ONLY: ink,rlk,ink

    implicit none

    REAL(KIND=8), INTENT(IN) :: ndx01,ndx02,ndx03,ndx04, &
    & ndy01,ndy02,ndy03,ndy04,ndx11,ndx12,ndx13,ndx14, &
    & ndy11,ndy12,ndy13,ndy14
    REAL(KIND=8) :: cut
    REAL(KIND=8), DIMENSION(4), INTENT(OUT) :: rDelv
    INTEGER(KIND=4) jj
    REAL(KIND=8) :: x1,x2,x3,x4,y1,y2,y3,y4,a1,a3,b1,b3

    !1.
    jj=1
    x1=ndx01
    x2=ndx02
    y1=ndy01
    y2=ndy02
    x3=ndx12
    x4=ndx11
    y3=ndy12
    y4=ndy11
    a1=0.25_8*(-x1+x2+x3-x4)
    a3=0.25_8*(-x1-x2+x3+x4)
    b1=0.25_8*(-y1+y2+y3-y4)
    b3=0.25_8*(-y1-y2+y3+y4)
    rdelv(jj)=4.0_8*(a1*b3-a3*b1)
    IF (rdelv(jj).LT.cut) rdelv(jj)=0.0_8
    !2.
    jj=2
    x1=ndx02
    x2=ndx03
    y1=ndy02
    y2=ndy03
    x3=ndx13
    x4=ndx12
    y3=ndy13
    y4=ndy12
    a1=0.25_8*(-x1+x2+x3-x4)
    a3=0.25_8*(-x1-x2+x3+x4)
    b1=0.25_8*(-y1+y2+y3-y4)
    b3=0.25_8*(-y1-y2+y3+y4)
    rdelv(jj)=4.0_8*(a1*b3-a3*b1)
    IF (rdelv(jj).LT.cut) rdelv(jj)=0.0_8
    !3.
    jj=3
    x1=ndx03
    x2=ndx04
    y1=ndy03
    y2=ndy04
    x3=ndx14
    x4=ndx13
    y3=ndy14
    y4=ndy13
    a1=0.25_8*(-x1+x2+x3-x4)
    a3=0.25_8*(-x1-x2+x3+x4)
    b1=0.25_8*(-y1+y2+y3-y4)
    b3=0.25_8*(-y1-y2+y3+y4)
    rdelv(jj)=4.0_8*(a1*b3-a3*b1)
    IF (rdelv(jj).LT.cut) rdelv(jj)=0.0_8
    !4.
    jj=4
    x1=ndx04
    x2=ndx01
    y1=ndy04
    y2=ndy01
    x3=ndx11
    x4=ndx14
    y3=ndy11
    y4=ndy14
    a1=0.25_8*(-x1+x2+x3-x4)
    a3=0.25_8*(-x1-x2+x3+x4)
    b1=0.25_8*(-y1+y2+y3-y4)
    b3=0.25_8*(-y1-y2+y3+y4)
    rdelv(jj)=4.0_8*(a1*b3-a3*b1)
    IF (rdelv(jj).LT.cut) rdelv(jj)=0.0_8

  END SUBROUTINE ale_getfvol_vol


END MODULE ale_getfvol_kernels
