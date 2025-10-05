MODULE precision_kind_mod
  IMPLICIT NONE
  PRIVATE

  ! Define a kind for double precision floating-point numbers for portability.
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(P=15, R=300)

  PUBLIC :: dp ! Expose the kind parameter
END MODULE precision_kind_mod
