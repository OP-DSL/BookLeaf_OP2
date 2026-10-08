
MODULE op2_bookleaf_consts

  USE op2_bookleaf_api

#ifdef OP2_TRANSLATOR
  USE op2_consts
#endif

  USE reals_mod, ONLY: dt_min, dt_initial, dt_max, cfl_sf, div_sf, dt_g, &
&                      ccut, zcut, zerocut, pcut, dencut, accut

  CONTAINS

  SUBROUTINE bookleaf_op2_init_const

    USE op2_bookleaf

    implicit none

    call op_decl_const(dt_min, 1, 'real(8)')
    call op_decl_const(dt_initial, 1, 'real(8)')
    call op_decl_const(dt_max, 1, 'real(8)')
    call op_decl_const(cfl_sf, 1, 'real(8)')
    call op_decl_const(div_sf, 1, 'real(8)')
    call op_decl_const(dt_g, 1, 'real(8)')
    call op_decl_const(ccut, 1, 'real(8)')
    call op_decl_const(zcut, 1, 'real(8)')
    call op_decl_const(zerocut, 1, 'real(8)')
    call op_decl_const(pcut, 1, 'real(8)')
    call op_decl_const(dencut, 1, 'real(8)')
    call op_decl_const(accut, 1, 'real(8)')

  END SUBROUTINE bookleaf_op2_init_const

END MODULE op2_bookleaf_consts
