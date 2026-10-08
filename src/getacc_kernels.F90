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


MODULE getacc_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getacc_scatter(cnmass,rho,cnwt,cnfx,cnfy, &
&                           ndmass1,ndmass2,ndmass3,ndmass4, &
&                           ndarea1,ndarea2,ndarea3,ndarea4, &
&                           ndub1,ndub2,ndub3,ndub4, &
&                           ndvb1,ndvb2,ndvb3,ndvb4)

    USE kinds_mod,ONLY: ink,rlk
    USE op2_bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(IN) :: rho
    REAL(KIND=8), INTENT(INOUT) :: ndarea1,ndarea2,ndarea3,ndarea4, &
&                           ndub1,ndub2,ndub3,ndub4, &
&                           ndvb1,ndvb2,ndvb3,ndvb4, &
&                           ndmass1,ndmass2,ndmass3,ndmass4
    REAL(KIND=8),DIMENSION(4),INTENT(IN) :: cnmass,cnwt,cnfx,cnfy
    INTEGER(KIND=4) :: jj,ii

    jj=1
    IF (cnmass(jj).GT.zerocut) THEN
      ndmass1=ndmass1+cnmass(jj)
    ELSE
      ii=jj-1
      IF (ii.EQ.0) ii=4
      IF (cnmass(ii).GT.zerocut) THEN
        ndmass1=ndmass1+cnmass(ii)
      ELSE
        ndmass1=ndmass1+rho*cnwt(jj)
      ENDIF
    ENDIF
    ndarea1=ndarea1+cnwt(jj)
    ndub1=ndub1+cnfx(jj)
    ndvb1=ndvb1+cnfy(jj)

    jj=2
    IF (cnmass(jj).GT.zerocut) THEN
      ndmass2=ndmass2+cnmass(jj)
    ELSE
      ii=jj-1
      IF (ii.EQ.0) ii=4
      IF (cnmass(ii).GT.zerocut) THEN
        ndmass2=ndmass2+cnmass(ii)
      ELSE
        ndmass2=ndmass2+rho*cnwt(jj)
      ENDIF
    ENDIF
    ndarea2=ndarea2+cnwt(jj)
    ndub2=ndub2+cnfx(jj)
    ndvb2=ndvb2+cnfy(jj)

    jj=3
    IF (cnmass(jj).GT.zerocut) THEN
      ndmass3=ndmass3+cnmass(jj)
    ELSE
      ii=jj-1
      IF (ii.EQ.0) ii=4
      IF (cnmass(ii).GT.zerocut) THEN
        ndmass3=ndmass3+cnmass(ii)
      ELSE
        ndmass3=ndmass3+rho*cnwt(jj)
      ENDIF
    ENDIF
    ndarea3=ndarea3+cnwt(jj)
    ndub3=ndub3+cnfx(jj)
    ndvb3=ndvb3+cnfy(jj)

    jj=4
    IF (cnmass(jj).GT.zerocut) THEN
      ndmass4=ndmass4+cnmass(jj)
    ELSE
      ii=jj-1
      IF (ii.EQ.0) ii=4
      IF (cnmass(ii).GT.zerocut) THEN
        ndmass4=ndmass4+cnmass(ii)
      ELSE
        ndmass4=ndmass4+rho*cnwt(jj)
      ENDIF
    ENDIF
    ndarea4=ndarea4+cnwt(jj)
    ndub4=ndub4+cnfx(jj)
    ndvb4=ndvb4+cnfy(jj)

  END SUBROUTINE getacc_scatter

  SUBROUTINE getacc_scatter2(cnmass,rho,cnwt,cnfx,cnfy, &
&                           ndmass1, &
&                           ndarea1, &
&                           ndub1, &
&                           ndvb1,jj)

    USE kinds_mod,ONLY: ink,rlk
    USE op2_bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(IN) :: rho
    REAL(KIND=8), INTENT(INOUT) :: ndarea1, &
&                           ndub1, &
&                           ndvb1, &
&                           ndmass1
    REAL(KIND=8),DIMENSION(4),INTENT(IN) :: cnmass,cnwt,cnfx,cnfy
    INTEGER(KIND=4), INTENT(IN) :: jj
    INTEGER(KIND=4) :: ii

    IF (cnmass(jj).GT.zerocut) THEN
      ndmass1=ndmass1+cnmass(jj)
    ELSE
      ii=jj-1
      IF (ii.EQ.0) ii=4
      IF (cnmass(ii).GT.zerocut) THEN
        ndmass1=ndmass1+cnmass(ii)
      ELSE
        ndmass1=ndmass1+rho*cnwt(jj)
      ENDIF
    ENDIF
    ndarea1=ndarea1+cnwt(jj)
    ndub1=ndub1+cnfx(jj)
    ndvb1=ndvb1+cnfy(jj)

  END SUBROUTINE getacc_scatter2

  SUBROUTINE getacc_accel(ndarea,ndmass,ndub,ndvb)

    USE kinds_mod,ONLY: ink,rlk
    USE op2_bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: ndub,ndvb,ndmass
    REAL(KIND=8), INTENT(IN) :: ndarea

    REAL(KIND=8) :: w1

    w1=dencut*ndarea
    IF (ndmass.GT.w1) THEN
      ndub=ndub/ndmass
      ndvb=ndvb/ndmass
    ELSE
      ndub=0.0_8
      ndvb=0.0_8
      ndmass=MAX(zerocut,w1)
    ENDIF

  END SUBROUTINE getacc_accel

  SUBROUTINE getacc_bc(ndub,ndvb,indtype)

    USE kinds_mod,ONLY: ink,rlk,ink
    USE op2_bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: ndub,ndvb
    INTEGER(KIND=4), INTENT(IN) :: indtype

    REAL(KIND=8) :: w2

    IF (indtype.EQ.-1) THEN
      ndub=0.0_8
    ELSEIF (indtype.EQ.-2) THEN
      ndvb=0.0_8
    ELSEIF (indtype.EQ.-3) THEN
      ndub=0.0_8
      ndvb=0.0_8
    ENDIF
    w2=ndub*ndub+ndvb*ndvb
    IF (w2.LT.(accut*accut)) THEN
      ndub=0.0_8
      ndvb=0.0_8
    ENDIF

  END SUBROUTINE getacc_bc

  SUBROUTINE getacc_avgvel(ndub,ndvb,ndu,ndv,dt,dt05)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: ndub,ndvb,ndu,ndv
    REAL(KIND=8), INTENT(INOUT) :: dt,dt05

    REAL(KIND=8) :: w1,w2

    w1=ndu
    w2=ndv
    ndu=w1+dt*ndub
    ndv=w2+dt*ndvb
    ndub=w1+dt05*ndub
    ndvb=w2+dt05*ndvb
  END SUBROUTINE getacc_avgvel

  SUBROUTINE getacc_updpos(ndub,ndvb,ndx,ndy,dt)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(INOUT) :: ndy,ndx,dt
    REAL(KIND=8), INTENT(IN) :: ndub,ndvb

    ndx=ndx+dt*ndub
    ndy=ndy+dt*ndvb

  END SUBROUTINE getacc_updpos


END MODULE getacc_kernels
