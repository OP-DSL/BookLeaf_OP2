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

MODULE gethg_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

    SUBROUTINE gethg_force(elu,elv,kappareg,rho,area,elfx,elfy,dt)
      USE kinds_mod,ONLY: rlk
      USE parameters_mod,ONLY: N_SHAPE

      implicit none
      REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: elu,elv
      REAL(KIND=rlk), INTENT(IN) :: kappareg,rho,area,dt
      REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(INOUT) :: elfx,elfy
      REAL(KIND=rlk) :: w1,w2,w3

      w2=elu(1)-elu(2)+elu(3)-elu(4)
      w3=elv(1)-elv(2)+elv(3)-elv(4)
      w1=-kappareg*rho*area/dt
      w2=w1*w2
      w3=w1*w3
      elfx(1)=elfx(1)+w2
      elfx(2)=elfx(2)-w2
      elfx(3)=elfx(3)+w2
      elfx(4)=elfx(4)-w2
      elfy(1)=elfy(1)+w3
      elfy(2)=elfy(2)-w3
      elfy(3)=elfy(3)+w3
      elfy(4)=elfy(4)-w3
    END SUBROUTINE gethg_force

  END MODULE gethg_kernels