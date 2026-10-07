! Constants visible to BookLeaf kernels.
!
! Normal BookLeaf builds obtain these configuration values from reals_mod.
! During OP2 translation they are renamed to op2_const_* and must instead be
! supplied by the backend-generated op2_consts module.  Keeping that choice
! here lets kernel sources remain backend-neutral.
MODULE bookleaf_consts

#ifdef OP2_TRANSLATOR
  USE op2_consts, ONLY: op2_const_dt_min, op2_const_dt_initial,              &
&                      op2_const_dt_max, op2_const_cfl_sf,                   &
&                      op2_const_div_sf, op2_const_dt_g,                     &
&                      op2_const_ccut, op2_const_zcut,                       &
&                      op2_const_zerocut, op2_const_pcut,                    &
&                      op2_const_dencut, op2_const_accut
  ! Retained original kernel modules still use the application names, while
  ! extracted kernels use their renamed op2_const_* counterparts.  Export
  ! both views of each generated value.
  USE op2_consts, ONLY: dt_min => op2_const_dt_min,                          &
&                      dt_initial => op2_const_dt_initial,                   &
&                      dt_max => op2_const_dt_max, cfl_sf => op2_const_cfl_sf,&
&                      div_sf => op2_const_div_sf, dt_g => op2_const_dt_g,   &
&                      ccut => op2_const_ccut, zcut => op2_const_zcut,       &
&                      zerocut => op2_const_zerocut, pcut => op2_const_pcut, &
&                      dencut => op2_const_dencut, accut => op2_const_accut
#else
  USE reals_mod, ONLY: dt_min, dt_initial, dt_max, cfl_sf, div_sf, dt_g,    &
&                      ccut, zcut, zerocut, pcut, dencut, accut
#endif

  IMPLICIT NONE

END MODULE bookleaf_consts
