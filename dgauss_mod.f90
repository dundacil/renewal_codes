MODULE dgauss_mod
  USE precision_kind_mod, ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: DGAUSS

CONTAINS

  RECURSIVE PURE FUNCTION DGAUSS(F, A, B, EPS) RESULT(integral_val)
    !! Adaptive recursive 8-point Gaussian quadrature.
    !! Thread- and MPI-safe (no shared state).
    USE precision_kind_mod, ONLY: dp
    IMPLICIT NONE

    INTERFACE
      PURE FUNCTION F(X) RESULT(res)
        USE precision_kind_mod, ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: X
        REAL(dp) :: res
      END FUNCTION F
    END INTERFACE

    REAL(dp), INTENT(IN) :: A, B, EPS
    REAL(dp) :: integral_val
    REAL(dp) :: S8, S16, H, C

    ! Compute central point and half-width
    H = (B - A) / 2.0_dp
    C = (B + A) / 2.0_dp

    ! Coarse 8-point integration
    S8 = DGAUSS_8(F, A, B)

    ! Finer estimate by splitting in half
    S16 = DGAUSS_8(F, A, C) + DGAUSS_8(F, C, B)

    ! Adaptive recursion
    IF (ABS(S16 - S8) < EPS) THEN
       integral_val = S16
    ELSE
       integral_val = DGAUSS(F, A, C, EPS / 2.0_dp) + DGAUSS(F, C, B, EPS / 2.0_dp)
    END IF

  END FUNCTION DGAUSS


  PURE FUNCTION DGAUSS_8(F, A, B) RESULT(integral_val)
    !! Single 8-point Gauss-Legendre quadrature.
    USE precision_kind_mod, ONLY: dp
    IMPLICIT NONE

    INTERFACE
      PURE FUNCTION F(X) RESULT(res)
        USE precision_kind_mod, ONLY: dp
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: X
        REAL(dp) :: res
      END FUNCTION F
    END INTERFACE

    REAL(dp), INTENT(IN) :: A, B
    REAL(dp) :: integral_val
    REAL(dp) :: H, C
    INTEGER :: N

    ! abscissas and weights (8-point Gauss–Legendre)
    REAL(dp), PARAMETER :: XG(8) = [ &
         0.0950125098376374_dp, 0.2816035507792589_dp, &
         0.4580167776572274_dp, 0.6178762444026438_dp, &
         0.7554044083550030_dp, 0.8656312023878318_dp, &
         0.9445750230732326_dp, 0.9894009349916499_dp ]

    REAL(dp), PARAMETER :: WG(8) = [ &
         0.1894506104550685_dp, 0.1826034150449236_dp, &
         0.1691565193950025_dp, 0.1495959888165767_dp, &
         0.1246289712555339_dp, 0.0951585116824928_dp, &
         0.0622535239386479_dp, 0.0271524594117541_dp ]

    H = (B - A) / 2.0_dp
    C = (B + A) / 2.0_dp

    integral_val = 0.0_dp
    DO N = 1, 8
       integral_val = integral_val + WG(N) * (F(C + H * XG(N)) + F(C - H * XG(N)))
    END DO
    integral_val = integral_val * H

  END FUNCTION DGAUSS_8

END MODULE dgauss_mod

