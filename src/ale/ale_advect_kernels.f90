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

MODULE ale_advect_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE ale_advect_update(elv,elm,elr,elvpr,elmpr,elrpr, &
&                              cutv,cutm,totv,totm)

    USE kinds_mod,ONLY: rlk
    USE reals_mod,    ONLY: zerocut, dencut

    implicit none

    REAL(KIND=rlk), INTENT(INOUT) :: elvpr,elmpr,elrpr,cutv,cutm,elv,elm,elr
    REAL(KIND=rlk), INTENT(IN) :: totv,totm



    ! store basis variables
      elvpr=elv
      elmpr=elm
      elrpr=elr
      ! construct cut-off's
      cutv=zerocut
      cutm=elvpr*dencut
      ! volume
      elv=elv+totv
      ! mass
      elm=elm+totm
      ! density
      elr=elm/elv

  END SUBROUTINE ale_advect_update

  SUBROUTINE ale_advect_prevolmass(elv0ndm1,elv1,cnm0,cnm1, &
    &      ndv00,ndv10,ndm00,ndv01,ndv11,ndm01,ndv02,ndv12,ndm02,ndv03,ndv13,ndm03)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: elv0ndm1,elv1
    REAL(KIND=rlk), INTENT(INOUT) :: ndv00,ndv10,ndm00,ndv01,ndv11,ndm01,ndv02,ndv12,ndm02,ndv03,ndv13,ndm03
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: cnm0
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: cnm1
    REAL(KIND=rlk) :: w1,w2,w3

    w1=0.25_rlk*elv0ndm1
    w2=0.25_rlk*elv1
    !1st
    w3=cnm1(1)
    cnm0(1)=w3
    ndv00=ndv00+w1
    ndv10=ndv10+w2
    ndm00=ndm00+w3
    !2nd
    w3=cnm1(2)
    cnm0(2)=w3
    ndv01=ndv01+w1
    ndv11=ndv11+w2
    ndm01=ndm01+w3
    !3rd
    w3=cnm1(3)
    cnm0(3)=w3
    ndv02=ndv02+w1
    ndv12=ndv12+w2
    ndm02=ndm02+w3
    !4th
    w3=cnm1(4)
    cnm0(4)=w3
    ndv03=ndv03+w1
    ndv13=ndv13+w2
    ndm03=ndm03+w3

  END SUBROUTINE ale_advect_prevolmass

  SUBROUTINE ale_advect_volmass(delv,delv1,delv2,delv3, delv4, &
&       delm,delm1,delm2,delm3,delm4,dndv,dndm,flux,ielsd,ielel,iel,i1)
    USE kinds_mod,ONLY: rlk, ink
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: delv,delv1,delv2,delv3, delv4, &
&       delm,delm1,delm2,delm3,delm4
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(INOUT) :: flux
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: dndv,dndm
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: ielel,ielsd
    INTEGER(KIND=ink), INTENT(IN) :: iel,i1

    INTEGER(KIND=ink) :: i2,ie1,ie2,is1,is2
    REAL(KIND=rlk)  :: w1,w2,w3,w4

    i2 = i1+2_ink

    ie1=ielel(i1)
    ie2=ielel(i2)
    is1=ielsd(i1)
    is2=ielsd(i2)
    IF (i1.EQ.1) THEN !i1 may only be 1 or 2, so that i2<=4
      w1=delv1(is1)
      w2=delv3(is2)
      w3=delm1(is1)
      w4=delm3(is2)
    ELSEIF (i1.EQ.2) THEN
      w1=delv2(is1)
      w2=delv4(is2)
      w3=delm2(is1)
      w4=delm4(is2)
    ENDIF

    IF (ie1.EQ.iel) THEN !i.e. boundary
      w1=0.0_rlk
      w3=0.0_rlk
    ENDIF
    IF (ie2.EQ.iel) THEN !i.e. boundary
      w2=0.0_rlk
      w4=0.0_rlk
    ENDIF
    w1=w1-delv(i1)
    w2=w2-delv(i2)
    w1=0.25_rlk*(w1-w2)
    dndv(i1)=w1
    dndv(i2)=w1
    w1=w3-delm(i1)
    w2=w4-delm(i2)
    w3=0.25_rlk*(w1-w2)
    dndm(i1)=w3
    dndm(i2)=w3
    w3=0.25_rlk*(w1+w2)
    flux(1)=flux(1)+w3
    flux(2)=flux(2)+w3
    flux(3)=flux(3)+w3
    flux(4)=flux(4)+w3
  END SUBROUTINE ale_advect_volmass

  SUBROUTINE ale_advect_postmass(cnm1,flux, &
&             elv0ndm11,elv0ndm12,elv0ndm13,elv0ndm14)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: flux
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(INOUT) :: cnm1
    REAL(KIND=rlk), INTENT(INOUT) :: elv0ndm11,elv0ndm12,elv0ndm13,elv0ndm14

    cnm1(1)=cnm1(1)+flux(1)
    cnm1(2)=cnm1(2)+flux(2)
    cnm1(3)=cnm1(3)+flux(3)
    cnm1(4)=cnm1(4)+flux(4)
    elv0ndm11=elv0ndm11+flux(1)
    elv0ndm12=elv0ndm12+flux(2)
    elv0ndm13=elv0ndm13+flux(3)
    elv0ndm14=elv0ndm14+flux(4)

  END SUBROUTINE ale_advect_postmass

  SUBROUTINE ale_advect_cutoff(cutv,cutm,ndv0,cut)
    USE kinds_mod,ONLY: rlk
    USE reals_mod,    ONLY: dencut

    implicit none

    REAL(KIND=rlk), INTENT(IN) :: ndv0,cut
    REAL(KIND=rlk), INTENT(OUT) :: cutv,cutm

    cutv=cut
    cutm=dencut*ndv0
  END SUBROUTINE ale_advect_cutoff

  SUBROUTINE ale_advect_markactive(indstatus,indtype,zactive)
    USE kinds_mod,ONLY: ink

    implicit none

    INTEGER(KIND=ink), INTENT(IN) :: indstatus,indtype
    INTEGER(KIND=ink), INTENT(OUT) :: zactive

    IF ((indstatus.GT.0_ink).AND.(indtype.NE.-1_ink).AND.   &
&       (indtype.NE.-3_ink)) THEN
      zactive=1_ink !.TRUE._lok
    ELSE
      zactive=0_ink !.FALSE._lok
    ENDIF
  END SUBROUTINE ale_advect_markactive

  SUBROUTINE ale_advect_markactive2(indstatus,indtype,zactive)
    USE kinds_mod,ONLY: ink

    implicit none

    INTEGER(KIND=ink), INTENT(IN) :: indstatus,indtype
    INTEGER(KIND=ink), INTENT(OUT) :: zactive

    IF ((indstatus.GT.0_ink).AND.(indtype.NE.-2_ink).AND.   &
&         (indtype.NE.-3_ink)) THEN
      zactive=1_ink !.TRUE._lok
    ELSE
      zactive=0_ink !.FALSE._lok
    ENDIF
  END SUBROUTINE ale_advect_markactive2


  SUBROUTINE getq_christiensen_q(elx,ely,elu,elv,qx,qy,qq)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE
    USE reals_mod,    ONLY: zerocut

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: elx,ely,elu,elv
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: qx,qy
    REAL(KIND=rlk), INTENT(OUT) :: qq
    INTEGER(KIND=ink) :: iside, ins
    REAL(KIND=rlk) :: w1,w2,w3,w4,w5,w6,w7,w8,den,uhat,vhat,xhat,yhat

    DO iside=1,N_SHAPE
      ins=MOD(iside,N_SHAPE)+1_ink

      w1=elx(iside)
      w2=elx(ins)
      w3=0.5_rlk*(w1+w2)
      w1=w2-w1
      w2=0.25_rlk*(elx(1)+elx(2)+elx(3)+elx(4))
      w4=ely(iside)
      w5=ely(ins)
      w6=0.5_rlk*(w4+w5)
      w4=w5-w4
      w5=0.25_rlk*(ely(1)+ely(2)+ely(3)+ely(4))
      w7=SQRT((w2-w3)*(w2-w3)+(w5-w6)*(w5-w6))
      w8=SQRT(w1*w1+w4*w4)
      den=1.0_rlk/w7
      xhat=(w5-w6)*den
      yhat=(w3-w2)*den
      den=1.0_rlk/w8
      w1=w1*den
      w2=w4*den
      w3=xhat*w1+yhat*w2
      den=-SIGN(1.0_rlk,w3)*w7
      xhat=xhat*den
      yhat=yhat*den
      uhat=elu(ins)-elu(iside)
      vhat=elv(ins)-elv(iside)
      w5=SQRT((uhat*uhat)+(vhat*vhat))
      w6=uhat*xhat+vhat*yhat
      den=w6/MAX(w5,zerocut)
      qx(iside)=qx(iside)*uhat*den
      qy(iside)=qy(iside)*vhat*den
      IF ((w5.LE.zerocut).OR.(w6.LE.zerocut).OR.(w7.LE.zerocut).OR.   &
&           (w8.LE.zerocut)) THEN
        qx(iside)=0.0_rlk
        qy(iside)=0.0_rlk
      ENDIF
      qq=qq+0.25_rlk*SQRT(qx(iside)*qx(iside)+      &
&               qy(iside)*qy(iside))

    ENDDO

  END SUBROUTINE getq_christiensen_q

END MODULE ale_advect_kernels
