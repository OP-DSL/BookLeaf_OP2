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


MODULE getdt_kernels

  USE kinds_mod,ONLY: ink,rlk

  CONTAINS

  SUBROUTINE getdt_cfl(rscratch11,rscratch12,rho,csqrd,qq,elx,ely,ielreg,zdtnotreg,zmidlength)

    USE kinds_mod,ONLY: ink,rlk,ink
!    USE geometry_mod,    ONLY: dlm,dln
    USE bookleaf_consts

    implicit none

    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: elx,ely
    REAL(KIND=8), INTENT(IN) :: rho,csqrd,qq
    REAL(KIND=8), INTENT(OUT) :: rscratch11,rscratch12
    INTEGER(KIND=4), DIMENSION(100) :: zdtnotreg,zmidlength !need to inline this
    INTEGER(KIND=4) :: ielreg

    REAL(KIND=8) :: w1,w2,w3
    !For dlm and dln
    REAL(KIND=8)                              :: x1,x2,y1,y2
    REAL(KIND=8),DIMENSION(4)            :: res

    IF (zdtnotreg(ielreg+1).EQ.1) THEN
      rscratch11=dt_max
      rscratch12=2.2250738585072014D-308
    ELSE
      w1=MAX(rho,zcut)
      w2=MAX(ccut,csqrd)+2.0_8*qq/w1
      IF (zmidlength(ielreg+1).EQ.1) THEN
        x1=elx(1)+elx(2)
        x2=elx(3)+elx(4)
        y1=ely(1)+ely(2)
        y2=ely(3)+ely(4)
        x1=0.5_8*(x1-x2)
        y1=0.5_8*(y1-y2)
        res(1)=x1*x1+y1*y1
        x1=elx(3)+elx(2)
        x2=elx(1)+elx(4)
        y1=ely(3)+ely(2)
        y2=ely(1)+ely(4)
        x1=0.5_8*(x1-x2)
        y1=0.5_8*(y1-y2)
        res(2)=x1*x1+y1*y1
        res(3)=res(1)
        res(4)=res(2)
        w1=MIN(res(1),res(2),res(3),res(4))
      ELSE
        w3=(ely(3)-ely(4))* &
    &      (ely(3)-ely(4))+ &
    &      (elx(3)-elx(4))* &
    &      (elx(3)-elx(4))
        !denom(elx(3),ely(3),elx(4),ely(4))
        IF (w3.LT.zcut) THEN
          res(1)=(0.5_8*(elx(1)+elx(2)) &
    &           -elx(3))*(0.5_8*(elx(1)+elx(2))-elx(3))+ &
    &         (0.5_8*(ely(1)+ely(2))-ely(3))* &
    &         (0.5_8*(ely(1)+ely(2))-ely(3))
          !distpp(elx(1),ely(1),elx(2),ely(2),elx(3),ely(3))
        ELSE
          res(1)= 0.5_8*(ely(3)-ely(4))* &
    &             (elx(1)+elx(2))+0.5_8* &
    &             (ely(1)+ely(2))*(elx(4)-elx(3))+ &
    &              ely(4)*elx(3)-ely(3)*elx(4)
          res(1) = res(1)*res(1)/w3
    !      distpl(elx(1),ely(1),elx(2),ely(2),elx(3),ely(3),elx(4),   &
    !&                   ely(4))/w3
        ENDIF
        w3=(ely(4)-ely(1))* &
    &      (ely(4)-ely(1))+ &
    &      (elx(4)-elx(1))* &
    &      (elx(4)-elx(1))
        !denom(elx(4),ely(4),elx(1),ely(1))
        IF (w3.LT.zcut) THEN
          res(2)=(0.5_8*(elx(2)+elx(3)) &
    &           -elx(4))*(0.5_8*(elx(2)+elx(3))-elx(4))+ &
    &         (0.5_8*(ely(2)+ely(3))-ely(4))* &
    &         (0.5_8*(ely(2)+ely(3))-ely(4))
          !distpp(elx(2),ely(2),elx(3),ely(3),elx(4),ely(4))
        ELSE
          res(2)= 0.5_8*(ely(4)-ely(1))* &
    &             (elx(2)+elx(3))+0.5_8* &
    &             (ely(2)+ely(3))*(elx(1)-elx(4))+ &
    &              ely(1)*elx(4)-ely(4)*elx(1)
          res(2) = res(2)*res(2)/w3
    !      distpl(elx(2),ely(2),elx(3),ely(3),elx(4),ely(4),elx(1),   &
    !&                   ely(1))/w3
        ENDIF
        w3=(ely(1)-ely(2))* &
    &      (ely(1)-ely(2))+ &
    &      (elx(1)-elx(2))* &
    &      (elx(1)-elx(2))
        !denom(elx(1),ely(1),elx(2),ely(2))
        IF (w3.LT.zcut) THEN
          res(3)=(0.5_8*(elx(3)+elx(4)) &
    &           -elx(1))*(0.5_8*(elx(3)+elx(4))-elx(1))+ &
    &         (0.5_8*(ely(3)+ely(4))-ely(1))* &
    &         (0.5_8*(ely(3)+ely(4))-ely(1))
          !distpp(elx(3),ely(3),elx(4),ely(4),elx(1),ely(1))
        ELSE
          res(3)= 0.5_8*(ely(1)-ely(2))* &
    &             (elx(3)+elx(4))+0.5_8* &
    &             (ely(3)+ely(4))*(elx(2)-elx(1))+ &
    &              ely(2)*elx(1)-ely(1)*elx(2)
          res(3) = res(3)*res(3)/w3
    !      distpl(elx(3),ely(3),elx(4),ely(4),elx(1),ely(1),elx(2),   &
    !&                   ely(2))/w3
        ENDIF
        w3=(ely(2)-ely(3))* &
    &      (ely(2)-ely(3))+ &
    &      (elx(2)-elx(3))* &
    &      (elx(2)-elx(3))
        !denom(elx(2),ely(2),elx(3),ely(3))
        IF (w3.LT.zcut) THEN
          res(4)=(0.5_8*(elx(4)+elx(1)) &
    &           -elx(2))*(0.5_8*(elx(4)+elx(1))-elx(2))+ &
    &         (0.5_8*(ely(4)+ely(1))-ely(2))* &
    &         (0.5_8*(ely(4)+ely(1))-ely(2))
          !distpp(elx(4),ely(4),elx(1),ely(1),elx(2),ely(2))
        ELSE
          res(4)=0.5_8*(ely(2)-ely(3))*( &
    &             elx(4)+elx(1))+0.5_8*( &
    &             ely(4)+ely(1))*(elx(3)-elx(2))+ &
    &             ely(3)*elx(2)-ely(2)*elx(3)
          res(4) = res(4)*res(4)/w3
    !      distpl(elx(4),ely(4),elx(1),ely(1),elx(2),ely(2),elx(3),   &
    !&                   ely(3))/w3
        ENDIF
        w1=MIN(res(1),res(2),res(3),res(4))
      ENDIF
      rscratch11=w1/w2
      rscratch12=w1
    ENDIF

  END SUBROUTINE getdt_cfl

  SUBROUTINE getdt_minval(rscratch11,w1)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(IN) :: rscratch11
    REAL(KIND=8), INTENT(INOUT) :: w1

    w1 = MIN(w1,rscratch11)

  END SUBROUTINE getdt_minval

  SUBROUTINE getdt_minloc(rscratch11,elidx,w1,ii)

    USE kinds_mod,ONLY: ink,rlk

    implicit none

    REAL(KIND=8), INTENT(IN) :: rscratch11
    INTEGER(KIND=4), INTENT(INOUT) :: ii
    INTEGER(KIND=4), INTENT(IN) :: elidx
    REAL(KIND=8), INTENT(INOUT) :: w1

    IF ((w1.EQ.rscratch11).AND.(ii.GT.elidx)) THEN
      ii = elidx
    ELSE
      ii = MIN(ii,2147483647)
    ENDIF

  END SUBROUTINE getdt_minloc

  SUBROUTINE getdt_div(elu,elv,a1,a3,b1,b3,elvol,w2_s,w2)

    USE kinds_mod,ONLY: ink,rlk,ink

    implicit none

    REAL(KIND=8), INTENT(IN) :: a1,a3,b1,b3,elvol
    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: elu,elv
    REAL(KIND=8), INTENT(INOUT) :: w2,w2_s
!    INTEGER(KIND=4), INTENT(IN) :: iel
!    INTEGER(KIND=4), INTENT(INOUT) :: ii

    REAL(KIND=8) :: w1

    w1=elu(1)*(-b3+b1)+elv(1)*( a3-a1)+   &
&        elu(2)*( b3+b1)+elv(2)*(-a3-a1)+   &
&        elu(3)*( b3-b1)+elv(3)*(-a3+a1)+   &
&        elu(4)*(-b3-b1)+elv(4)*( a3+a1)
    w1=ABS(w1)/elvol
    w2_s=w1
    IF (w1.GT.w2) THEN
      w2=w1
    ENDIF

  END SUBROUTINE getdt_div

  SUBROUTINE getdt_div_maxloc(iel,w2_s,w2,ii)
    USE kinds_mod,ONLY: ink,rlk,ink
    REAL(KIND=8), INTENT(IN) :: w2_s
    REAL(KIND=8), INTENT(INOUT) :: w2
    INTEGER(KIND=4), INTENT(IN) :: iel
    INTEGER(KIND=4), INTENT(INOUT) :: ii

    IF (w2_s.EQ.w2) THEN
      ii=iel
    ENDIF
  END SUBROUTINE getdt_div_maxloc

  SUBROUTINE getdt_ale_zeul(elu,elv,rscratch12,w2_s,w2)

    USE kinds_mod,ONLY: ink,rlk,ink
    USE bookleaf_consts

    implicit none

    REAL(KIND=8), INTENT(IN) :: rscratch12
    REAL(KIND=8), DIMENSION(4), INTENT(IN) :: elu,elv
    REAL(KIND=8), INTENT(INOUT) :: w2,w2_s
!    INTEGER(KIND=4), INTENT(IN) :: iel
!    INTEGER(KIND=4), INTENT(INOUT) :: ii

    REAL(KIND=8) :: w1

          w1=MAX(elu(1)*elu(1)+elv(1)*elv(1),           &
&                elu(2)*elu(2)+elv(2)*elv(2),           &
&                elu(3)*elu(3)+elv(3)*elv(3),           &
&                elu(4)*elu(4)+elv(4)*elv(4))
    w1=rscratch12/MAX(w1,zerocut)
    w2_s=w1
    IF (w1.LT.w2) THEN
      w2=w1
    ENDIF

  END SUBROUTINE getdt_ale_zeul

  SUBROUTINE getdt_ale_zeul_minloc(iel,w2_s,w2,ii)
    USE kinds_mod,ONLY: ink,rlk,ink
    implicit none
    REAL(KIND=8), INTENT(IN) :: w2_s
    REAL(KIND=8), INTENT(INOUT) :: w2
    INTEGER(KIND=4), INTENT(IN) :: iel
    INTEGER(KIND=4), INTENT(INOUT) :: ii

    IF (w2_s.EQ.w2) THEN
      ii=iel
    ENDIF
  END SUBROUTINE getdt_ale_zeul_minloc

  SUBROUTINE getdt_mindt_reg(elidx,ielreg,idx,reg)
    USE kinds_mod,ONLY: ink,rlk,ink
    implicit none
    INTEGER(KIND=4), INTENT(IN) :: elidx,ielreg
    INTEGER(KIND=4) :: idx
    INTEGER(KIND=4), INTENT(INOUT) :: reg
    IF (idx.EQ.elidx) then
      reg = ielreg+1
    ENDIF
  END SUBROUTINE getdt_mindt_reg
END MODULE getdt_kernels
