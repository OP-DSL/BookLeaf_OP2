! Central Fortran OP2 API surface for BookLeaf.
!
! Application procedures import this module when they need OP2 types,
! argument constructors, loop interfaces, runtime routines, or HDF5 I/O.
! Keeping the external OP2 module imports here avoids repeating backend
! implementation modules throughout BookLeaf sources.
MODULE op2_bookleaf_api

  USE OP2_Fortran_Declarations
  USE OP2_Fortran_Reference
  USE OP2_Fortran_RT_Support
  USE OP2_Fortran_hdf5_Declarations

  IMPLICIT NONE

END MODULE op2_bookleaf_api
