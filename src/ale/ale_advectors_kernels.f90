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

MODULE ale_advectors_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE ale_advectors_flux(rCorner,rCorner1,rCorner2,rCorner3,rCorner4, &
&   rVar,rVar1,rVar2,rVar3,rVar4,rDel,rFlux,ielsd,i1)
    USE kinds_mod,ONLY: rlk, ink
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: rCorner,rCorner1,rCorner2,rCorner3,rCorner4,rDel
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: rFlux
    REAL(KIND=rlk), INTENT(IN) :: rVar,rVar1,rVar2,rVar3,rVar4
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: ielsd
    INTEGER(KIND=ink), INTENT(IN) :: i1

    INTEGER(KIND=ink) :: i2,j1,j2
    REAL(KIND=rlk) :: r1,r2,r3,r4,w1,w2,w3,w4,w5,w6,w7,w8,rV,rGrad

    i2 = i1+2_ink


!    iE2=iElEl(i2)
    j2=iElSd(i2)
    j1=MOD(i2,N_SHAPE)+1_ink
    r3=rCorner(i2)+rCorner(j1)
    j1=MOD(j2,N_SHAPE)+1_ink
    IF (i2.EQ.3) THEN
      w5=r3+rCorner3(j2)+rCorner3(j1)
    ELSEIF (i2.EQ.4) THEN
      w5=r3+rCorner4(j2)+rCorner4(j1)
    ENDIF
!    iE1=iElEl(i1)
    j2=iElSd(i1)
    j1=i1+1_ink
    r4=rCorner(i1)+rCorner(j1)
    j1=MOD(j2,N_SHAPE)+1_ink
    IF (i1.EQ.1) THEN
      w6=r4+rCorner1(j2)+rCorner1(j1)
    ELSEIF (i1.EQ.2) THEN
      w6=r4+rCorner2(j2)+rCorner2(j1)
    ENDIF
    rV=rVar
    r1=rDel(i1)
    r2=rDel(i2)
    IF (i2.EQ.3) THEN
      w1=rV-rVar3
    ELSEIF (i2.EQ.4) THEN
      w1=rV-rVar4
    ENDIF
    IF (i1.EQ.1) THEN
      w2=rVar1-rV
    ELSEIF (i1.EQ.2) THEN
      w2=rVar2-rV
    ENDIF
    w3=ABS(w1)
    w4=ABS(w2)
    w7=SIGN(1.0_rlk,w2)
    w8=(w4*w6*w6+w3*w5*w5)/(w5*w6*(w5+w6))
    rGrad=w7*MIN(ABS(w8),w3/w5,w4/w6)
    IF (w1*w2.LE.0.0_rlk) rGrad=0.0_rlk
    r1=r1*(rV+rGrad*(r3-0.5_rlk*r1))
    r2=r2*(rV-rGrad*(r4-0.5_rlk*r2))
    rFlux(i1)=r1
    rFlux(i2)=r2

  END SUBROUTINE ale_advectors_flux

  SUBROUTINE ale_advectors_flux2(rCorner,rCorner1,rCorner2,rCorner3,rCorner4, &
&   rVar,rVar1,rVar2,rVar3,rVar4,rDel,rFlux,ielsd,iFaceL,iCorner)
    USE kinds_mod,ONLY: rlk, ink
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: rCorner,rCorner1,rCorner2,rCorner3,rCorner4,rDel
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: rFlux
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: rVar,rVar1,rVar2,rVar3,rVar4
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: ielsd
    INTEGER(KIND=ink), INTENT(IN) :: iFaceL,iCorner

    INTEGER(KIND=ink) :: iFaceR,iSdL,iSdR, &
&                        iLNdL,iLNdR,iLNNdL,iLNNdR,ii
    REAL(KIND=rlk)    :: w1,w2,w3,w4,w5,w6,w7,w8,rD,rGrad,rV

    iFaceR = iFaceL+2_ink
    iLNdL=iFaceL+iCorner-1_ink
    iLNdR=MOD(iFaceR-iCorner+1,N_SHAPE)+1_ink

    rD=0.0_rlk
    iSdL=ielsd(iFaceL)
    iSdR=ielsd(iFaceR)
    ii=iFaceL+2_ink*(iCorner-1_ink)
    IF (rDel(ii).GT.0.0_rlk) THEN
      iLNNdL=MOD(iSdL+iCorner,N_SHAPE)+1_ink
      iLNNdR=MOD(iSdL-iCorner+1,N_SHAPE)+1_ink
      rV=rVar(iLndL)
      rD=rCorner(iLNdL)-0.5_rlk*rDel(ii)
      w5=rCorner(iLNdL)+rCorner(iLNdR)
      IF (iFaceL.EQ.1_ink) THEN
        w6=rCorner1(iLNNdL)+rCorner1(iLNNdR)
        w1=rV-rVar1(iLNNdL)
      ELSEIF (iFaceL.EQ.1_ink) THEN
        w6=rCorner2(iLNNdL)+rCorner2(iLNNdR)
        w1=rV-rVar2(iLNNdL)
      ENDIF

      w2=rVar(iLNdR)-rV
      w3=ABS(w1)
      w4=ABS(w2)
      w7=SIGN(1.0_rlk,w2)
      w8=(w4*w6*w6+w3*w5*w5)/(w5*w6*(w5+w6))
      rGrad=w7*MIN(ABS(w8),w3/w5,w4/w6)
      IF (w1*w2.LE.0.0_rlk) rGrad=0.0_rlk
      rD=rDel(ii)*(rV+rGrad*rD)
    ENDIF
    IF (rDel(ii).LT.0.0_rlk) THEN
      iLNNdL=MOD(iSdR+iCorner-2,N_SHAPE)+1_ink
      iLNNdR=MODULO(iSdR-iCorner-1,N_SHAPE)+1_ink
      rV=rVar(iLNdR)
      rD=rCorner(iLNdR)+0.5_rlk*rDel(ii)
      w5=rCorner(iLNdL)+rCorner(iLNdR)
      IF (iFaceR.EQ.3_ink) THEN
        w6=rCorner3(iLNNdL)+rCorner3(iLNNdR)  ! nel2
        w1=rV-rVar(iLNdL)
        w2=rVar3(iLNNdR)-rV   ! nel2
      ELSEIF (iFaceR.EQ.4_ink) THEN
        w6=rCorner4(iLNNdL)+rCorner4(iLNNdR)  ! nel2
        w1=rV-rVar(iLNdL)
        w2=rVar4(iLNNdR)-rV   ! nel2
      ENDIF
      w3=ABS(w1)
      w4=ABS(w2)
      w7=SIGN(1.0_rlk,w2)
      w8=(w4*w6*w6+w3*w5*w5)/(w5*w6*(w5+w6))
      rGrad=-w7*MIN(ABS(w8),w3/w5,w4/w6)
      IF (w1*w2.LE.0.0_rlk) rGrad=0.0_rlk
      rD=rDel(ii)*(rV+rGrad*rD)
    ENDIF
    rFlux(iLNdL)=rFlux(iLNdL)-rD ! nel1
    rFlux(iLNdR)=rFlux(iLNdR)+rD


  END SUBROUTINE ale_advectors_flux2

  SUBROUTINE ale_advectors_update_c1(rBase0,rBase1,rCutOff,rTotFlux,rVar)
    USE kinds_mod,ONLY: rlk

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: rBase0,rBase1,rCutOff,rTotFlux
    REAL(KIND=rlk), INTENT(OUT) :: rVar

    IF (rBase1.GT.rCutOff) THEN
        rVar=(rVar*rBase0+rTotFlux)/rBase1
      ENDIF
  END SUBROUTINE ale_advectors_update_c1

  SUBROUTINE ale_advectors_totflux(rFlux,rTotFlux1,rTotFlux2,rTotFlux3,rTotFlux4)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: rFlux
    REAL(KIND=rlk), INTENT(INOUT) :: rTotFlux1,rTotFlux2,rTotFlux3,rTotFlux4

    rTotFlux1 = rTotFlux1 + rFlux(1)
    rTotFlux2 = rTotFlux2 + rFlux(2)
    rTotFlux3 = rTotFlux3 + rFlux(3)
    rTotFlux4 = rTotFlux4 + rFlux(4)

  END SUBROUTINE ale_advectors_totflux

  SUBROUTINE ale_advectors_update_n1(rBase0,rBase1,rCut,rTotFlux,rVar,zActive)
    USE kinds_mod,ONLY: rlk

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: rBase0,rBase1,rCut,rTotFlux
    REAL(KIND=rlk), INTENT(OUT) :: rVar
    INTEGER(KIND=ink), INTENT(IN) :: zActive

    IF ((zActive.EQ.1_ink).AND.(rBase1.GT.rCut)) THEN
        rVar=(rVar*rBase0+rTotFlux)/rBase1
    ENDIF

  END SUBROUTINE ale_advectors_update_n1

  SUBROUTINE ale_advectors_sumflux(rFlux,rFlux1,rFlux2,rFlux3,rFlux4,rTotFlux,iElSd,iElEl,iel,i1)
    USE kinds_mod,ONLY: rlk, ink
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: rFlux,rFlux1,rFlux2,rFlux3,rFlux4
    REAL(KIND=rlk), INTENT(OUT) :: rTotFlux
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: iElSd,iElEl
    INTEGER(KIND=ink), INTENT(IN) :: i1,iel

    INTEGER(KIND=ink) :: i2,j1,j2,ie1,ie2
    REAL(KIND=rlk) :: w1,w2

    i2 = i1+2_ink
    iE1=iElEl(i1)
    iE2=iElEl(i2)
    j1=iElSd(i1)
    j2=iElSd(i2)
    IF (i1.EQ.1_ink) THEN
      w1=rFlux1(j1)
      w2=rFlux3(j2)
    ELSEIF (i1.EQ.2_ink) THEN
      w1=rFlux2(j1)
      w2=rFlux4(j2)
    ENDIF

    IF (iE1.EQ.iEl) w1=0.0_rlk
    IF (iE2.EQ.iEl) w2=0.0_rlk
    rTotFlux=rTotFlux-rFlux(i1)-rFlux(i2)+w1+w2
  END SUBROUTINE ale_advectors_sumflux

END MODULE ale_advectors_kernels