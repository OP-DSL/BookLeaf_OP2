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


MODULE getein_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getein_update(elfx,elfy,elu,elv,elmass,ein,ein_out,dt)

    USE kinds_mod,ONLY: rlk
    USE reals_mod,    ONLY: zerocut

    implicit none

    REAL(KIND=rlk), DIMENSION(4), INTENT(IN) :: elfx,elfy,elu,elv
    REAL(KIND=rlk), INTENT(IN) :: elmass,ein
    REAL(KIND=rlk), INTENT(INOUT) :: dt
    REAL(KIND=rlk), INTENT(OUT) :: ein_out

    REAL(KIND=rlk) :: w1

    w1=elfx(1)*elu(1)+elfy(1)*elv(1)+                 &
&   elfx(2)*elu(2)+elfy(2)*elv(2)+                 &
&   elfx(3)*elu(3)+elfy(3)*elv(3)+                 &
&   elfx(4)*elu(4)+elfy(4)*elv(4)
      !# Missing code here that can't be merged
    w1=-w1/MAX(elmass,zerocut)
    ein_out=ein+w1*dt

  END SUBROUTINE getein_update


END MODULE getein_kernels
