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


MODULE write_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

!   SUBROUTINE write_init(reg_vol,reg_ie,reg_pre, &
! &              reg_pmx,reg_pmn,reg_mass,reg_ke,reg_dmx,reg_dmn)

!     USE kinds_mod,ONLY: ink,rlk

!     implicit none

!     REAL(KIND=8), INTENT(OUT) :: reg_vol,reg_ie,reg_pre, &
! &               reg_pmx,reg_pmn,reg_mass,reg_ke,reg_dmx,reg_dmn

!     reg_vol=0.0_8
!     reg_ie=0.0_8
!     reg_pre=0.0_8
!     reg_pmx=-HUGE(1.0_8)
!     reg_pmn=HUGE(1.0_8)
!     reg_mass=0.0_8
!     reg_ke=0.0_8
!     reg_dmx=-HUGE(1.0_8)
!     reg_dmn=HUGE(1.0_8)

!   END SUBROUTINE write_init

    SUBROUTINE write_regvalues(elvol,elmass,ein,pre, &    !4
&               rho,cnmass,cnwt,ndu1,ndu2,ndu3,ndu4,ndv1, & !8
&               ndv2,ndv3,ndv4,reg_vol,reg_ie, reg_pre, & !6
&       reg_pmx,reg_pmn,reg_mass,reg_ke, reg_dmx,reg_dmn, & !6
&       tot_mom_u,tot_mom_v,ireg,myreg) !4:tab

    USE kinds_mod,ONLY: ink,rlk,ink
    USE op2_bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(IN) :: elvol,elmass,ein,pre, &
&               rho,ndu1,ndu2,ndu3,ndu4,ndv1,ndv2,ndv3,ndv4
    REAL(KIND=8), INTENT(INOUT) :: reg_vol,reg_ie, reg_pre, &
&       reg_pmx,reg_pmn,reg_mass,reg_ke, reg_dmx,reg_dmn, &
&       tot_mom_u, tot_mom_v
    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: cnwt,cnmass
    INTEGER(KIND=4), INTENT(IN) :: ireg,myreg
    REAL(KIND=8) :: c1,w2,w3,w4
    REAL(KIND=8), DIMENSION(4) :: ndua, ndva
    INTEGER(KIND=4) :: ii

    IF (ireg.EQ.(myreg+1)) THEN
      ndua(1) = ndu1
      ndua(2) = ndu2
      ndua(3) = ndu3
      ndua(4) = ndu4
      ndva(1) = ndv1
      ndva(2) = ndv2
      ndva(3) = ndv3
      ndva(4) = ndv4


      ! Condition
      c1=dencut*elvol
      ! Scatter element contributions to region
      reg_vol=reg_vol+elvol
      IF (elmass.GT.c1) THEN
        w2=elmass
        reg_mass=reg_mass+w2
        w3=ein
        w3=w3*w2
        reg_ie=reg_ie+w3
        w4=pre
        w3=w2*w4
        reg_pre=reg_pre+w3
        IF (w4.GT.reg_pmx) reg_pmx=w4
        IF (w4.LT.reg_pmn) reg_pmn=w4
        w4=rho
        IF (w4.GT.reg_dmx) reg_dmx=w4
        IF (w4.LT.reg_dmn) reg_dmn=w4
      ENDIF
      DO ii=1,4
        w2=ndua(ii)
        w3=ndva(ii)
        IF (elmass.GT.c1) THEN
          reg_ke=0.5_8*cnmass(ii)*(w2*w2+w3*w3)+reg_ke
        ENDIF
        w4=rho*cnwt(ii)
        tot_mom_u=tot_mom_u+w2*w4
        tot_mom_v=tot_mom_v+w3*w4
      ENDDO
    ENDIF

  END SUBROUTINE write_regvalues

END MODULE write_kernels
