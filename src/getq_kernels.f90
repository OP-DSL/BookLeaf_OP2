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

MODULE getq_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getq_gradcon(du,dv,dx,dy,elu,elv,elx,ely)

    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE
    USE reals_mod,    ONLY: zerocut

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: du,dv,dx,dy
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: elu,elv,elx,ely


    du(1)=elu(2)-elu(1)
    du(2)=elu(3)-elu(2)
    du(3)=elu(4)-elu(3)
    du(4)=elu(1)-elu(4)
    dv(1)=elv(2)-elv(1)
    dv(2)=elv(3)-elv(2)
    dv(3)=elv(4)-elv(3)
    dv(4)=elv(1)-elv(4)
    dx(1)=elx(2)-elx(1)
    dx(2)=elx(3)-elx(2)
    dx(3)=elx(4)-elx(3)
    dx(4)=elx(1)-elx(4)
    dy(1)=ely(2)-ely(1)
    dy(2)=ely(3)-ely(2)
    dy(3)=ely(4)-ely(3)
    dy(4)=ely(1)-ely(4)

  END SUBROUTINE getq_gradcon

  SUBROUTINE getq_christiensen1(du,dv,dx,dy, &
&             du1,du2,du3,du4,dv1,dv2,dv3,dv4,&
&             dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4,&
&            ielsd,scratch,iside)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE
    USE reals_mod,    ONLY: zerocut

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: du,dv,dx,dy, &
&             du1,du2,du3,du4,dv1,dv2,dv3,dv4,&
&             dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: scratch
    INTEGER(KIND=ink), INTENT(INOUT) :: iside
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: ielsd
    INTEGER(KIND=ink) :: is1,is2,ins
    REAL(KIND=rlk) :: w1,w2,w3,w4,den,uhat,vhat,xhat,yhat

    is1=MOD(iside+2_ink,N_SHAPE)+1_ink
    is2=iside+1_ink

    ! edge 1
    w1=du(is1)
    w2=dv(is1)
    w3=dx(is1)
    w4=dy(is1)
    den=SQRT(w1*w1+w2*w2)
    den=1.0_rlk/MAX(den,zerocut)
    uhat=w1*den
    vhat=w2*den
    den=SQRT(w3*w3+w4*w4)
    den=1.0_rlk/MAX(den,zerocut)
    xhat=w3*den
    yhat=w4*den
    den=w3*xhat+w4*yhat
    w1=(w1*uhat+w2*vhat)/SIGN(MAX(ABS(den),zerocut),den)
    w1=1.0_rlk/SIGN(MAX(ABS(w1),zerocut),w1)
    ins=ielsd(iside)
    ins=MOD(ins,N_SHAPE)+1_ink
    IF (iside.eq.1_ink) THEN
      den=dx1(ins)*xhat+dy1(ins)*yhat
      w2=(du1(ins)*uhat+dv1(ins)*vhat)/                        &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ELSE
      den=dx2(ins)*xhat+dy2(ins)*yhat
      w2=(du2(ins)*uhat+dv2(ins)*vhat)/                        &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ENDIF
    scratch(1)=w2*w1
    ins=ielsd(iside+2_ink)
    ins=MOD(ins+2_ink,N_SHAPE)+1_ink
    IF (iside.eq.1_ink) THEN
      den=dx3(ins)*xhat+dy3(ins)*yhat
      w3=(du3(ins)*uhat+dv3(ins)*vhat)/                        &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ELSE
      den=dx4(ins)*xhat+dy4(ins)*yhat
      w3=(du4(ins)*uhat+dv4(ins)*vhat)/                        &
&          SIGN(MAX(ABS(den),zerocut),den)
    ENDIF
    scratch(2)=w3*w1
    ! edge 2
    w1=du(is2)
    w2=dv(is2)
    w3=dx(is2)
    w4=dy(is2)
    den=SQRT(w1*w1+w2*w2)
    den=1.0_rlk/MAX(den,zerocut)
    uhat=w1*den
    vhat=w2*den
    den=SQRT(w3*w3+w4*w4)
    den=1.0_rlk/MAX(den,zerocut)
    xhat=w3*den
    yhat=w4*den
    den=w3*xhat+w4*yhat
    w1=(w1*uhat+w2*vhat)/SIGN(MAX(ABS(den),zerocut),den)
    w1=1.0_rlk/SIGN(MAX(ABS(w1),zerocut),w1)
    ins=ielsd(iside)
    ins=MOD(ins+2_ink,N_SHAPE)+1_ink
    IF (iside.eq.1_ink) THEN
      den=dx1(ins)*xhat+dy1(ins)*yhat
      w2=(du1(ins)*uhat+dv1(ins)*vhat)/                         &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ELSE
      den=dx2(ins)*xhat+dy2(ins)*yhat
      w2=(du2(ins)*uhat+dv2(ins)*vhat)/                         &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ENDIF
    scratch(3)=w2*w1
    ins=ielsd(iside+2_ink)
    ins=MOD(ins,N_SHAPE)+1_ink
    IF (iside.eq.1_ink) THEN
      den=dx3(ins)*xhat+dy3(ins)*yhat
      w3=(du3(ins)*uhat+dv3(ins)*vhat)/                         &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ELSE
      den=dx4(ins)*xhat+dy4(ins)*yhat
      w3=(du4(ins)*uhat+dv4(ins)*vhat)/                         &
  &        SIGN(MAX(ABS(den),zerocut),den)
    ENDIF
    scratch(4)=w3*w1

  END SUBROUTINE getq_christiensen1

  SUBROUTINE getq_christiensen_bc(ielel,scratch,indtype1,indtype2, &
&                               indtype3,indtype4,iside)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: scratch
    INTEGER(KIND=ink), INTENT(INOUT) :: iside
    INTEGER(KIND=ink), INTENT(IN) :: indtype1, indtype2, &
&                                    indtype3, indtype4
    INTEGER(KIND=ink), DIMENSION(N_SHAPE), INTENT(IN) :: ielel
    INTEGER(KIND=ink) :: in1,in2,ins,ic1,ic2
    INTEGER(kind=ink), DIMENSION(4) :: indtype


    indtype(1) = indtype1
    indtype(2) = indtype2
    indtype(3) = indtype3
    indtype(4) = indtype4

    ins=iside+2_ink
    in1=ielel(iside)
    in2=ielel(ins)
    IF (in1.EQ.0_ink) THEN
      ic1=iside
      ic2=MOD(iside,N_SHAPE)+1_ink
      IF (((indtype(ic1).LT.0_ink).AND.(indtype(ic2).LT.0_ink)).AND.&
&          (in2.NE.0_ink)) THEN
        scratch(1)=1.0_rlk
        scratch(3)=1.0_rlk
      ELSE
        scratch(1)=0.0_rlk
        scratch(3)=0.0_rlk
      ENDIF
    ENDIF
    IF (in2.EQ.0_ink) THEN
      ic1=ins
      ic2=MOD(ins,N_SHAPE)+1_ink
      IF (((indtype(ic1).LT.0_ink).AND.(indtype(ic2).LT.0_ink)).AND.&
&          (in1.NE.0_ink)) THEN
        scratch(2)=1.0_rlk
        scratch(4)=1.0_rlk
      ELSE
        scratch(2)=0.0_rlk
        scratch(4)=0.0_rlk
      ENDIF
    ENDIF

  END SUBROUTINE getq_christiensen_bc

  SUBROUTINE getq_christiensen_limiter(du,dv,qx,qy,csqrd,rho,scratch,cq1,cq2,iside)
    USE kinds_mod,ONLY: rlk
    USE parameters_mod,ONLY: N_SHAPE

    implicit none

    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(IN) :: du,dv,scratch
    REAL(KIND=rlk), INTENT(IN) :: csqrd,rho
    REAL(KIND=rlk), INTENT(INOUT) :: cq1,cq2
    REAL(KIND=rlk), DIMENSION(N_SHAPE), INTENT(OUT) :: qx,qy
    INTEGER(KIND=ink), INTENT(INOUT) :: iside
    INTEGER(KIND=ink) :: is1,is2
    REAL(KIND=rlk) :: w1,w2,w3,w4

    is1=MOD(iside+2_ink,N_SHAPE)+1_ink
    is2=iside+1_ink

    w1=cq1*SQRT(csqrd)
    w2=scratch(1)
    w3=scratch(2)
    w2=MIN(0.5_rlk*(w2+w3),2.0_rlk*w2,2.0_rlk*w3,1.0_rlk)
    w2=MAX(0.0_rlk,w2)
    w3=du(is1)
    w4=dv(is1)
    w3=SQRT(w3*w3+w4*w4)
    w3=(1.0_rlk-w2)*rho*(w1+cq2*w3)
    qx(is1)=w3
    qy(is1)=w3
    w2=scratch(3)
    w3=scratch(4)
    w2=MIN(0.5_rlk*(w2+w3),2.0_rlk*w2,2.0_rlk*w3,1.0_rlk)
    w2=MAX(0.0_rlk,w2)
    w3=du(is2)
    w4=dv(is2)
    w3=SQRT(w3*w3+w4*w4)
    w3=(1.0_rlk-w2)*rho*(w1+cq2*w3)
    qx(is2)=w3
    qy(is2)=w3

  END SUBROUTINE getq_christiensen_limiter

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

END MODULE getq_kernels
