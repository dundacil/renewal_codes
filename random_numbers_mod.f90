! =============================================================================
! MODULE: random_numbers_mod
!
! DESCRIPTION:
!   This module provides functions for generating random numbers from various
!   distributions. The code is a modernization of traditional FORTRAN 77
!   implementations.
!
! PUBLIC FUNCTIONS:
!   - ran2:   Generates a uniform deviate in (0, 1).
!   - rnor:   Generates a standard normal (Gaussian) deviate.
!   - power:  Generates a deviate from a power-law distribution.
!   - levy:   Generates a deviate from a stable Lévy distribution.
! =============================================================================
MODULE random_numbers_mod
  use precision_kind_mod, ONLY: dp  ! (Line 1)
  IMPLICIT NONE                    ! (Line 2)
  
  ! Expose the function for external use
  PUBLIC :: ran2, rnor, power, levy

CONTAINS

  ! ===========================================================================
  FUNCTION ran2(idum) RESULT(random_val)
    ! Modernized Fortran 90 version of the ran2 generator.

    IMPLICIT NONE

    ! --- Arguments ---
    INTEGER, INTENT(INOUT) :: idum

    ! --- Function Result ---
    REAL(dp) :: random_val

    ! --- Parameters from the original source ---
    INTEGER, PARAMETER :: IM1 = 2147483563
    INTEGER, PARAMETER :: IM2 = 2147483399
    INTEGER, PARAMETER :: IMM1 = IM1 - 1
    INTEGER, PARAMETER :: IA1 = 40014
    INTEGER, PARAMETER :: IA2 = 40692
    INTEGER, PARAMETER :: IQ1 = 53668
    INTEGER, PARAMETER :: IQ2 = 52774
    INTEGER, PARAMETER :: IR1 = 12211
    INTEGER, PARAMETER :: IR2 = 3791
    INTEGER, PARAMETER :: NTAB = 32
    INTEGER, PARAMETER :: NDIV = 1 + IMM1 / NTAB
    REAL(dp), PARAMETER :: AM = 1.0_dp / IM1
    REAL(dp), PARAMETER :: EPS = 1.2e-14_dp
    REAL(dp), PARAMETER :: RNMX = 1.0_dp - EPS

    ! --- Local Variables ---
    ! Using the SAVE attribute with initialization replaces the old DATA statement.
    INTEGER, SAVE :: idum2 = 123456789
    INTEGER, DIMENSION(NTAB), SAVE :: iv = 0
    INTEGER, SAVE :: iy = 0
    INTEGER :: j, k

    ! > Long period (> 2e18) random number generator of L’Ecuyer with Bays-Durham shuffle
    ! > and added safeguards. [cite: 1]
    ! > Returns a uniform random deviate between 0.0 and 1.0 (exclusive
    ! > of the endpoint values). [cite: 2]
    ! > RNMX should approximate the largest floating value that is less than 1. [cite: 4]
    
    ! Initialization block.
    ! > Call with idum a negative integer to initialize; thereafter, do not
    ! > alter idum between successive deviates in a sequence. [cite: 3]
    IF (idum <= 0) THEN
      ! > Be sure to prevent idum = 0. [cite: 5]
      idum = MAX(-idum, 1)
      idum2 = idum

      ! > Load the shuffle table (after 8 warm-ups). [cite: 6]
      DO j = NTAB + 8, 1, -1
        k = idum / IQ1
        idum = IA1 * (idum - k * IQ1) - k * IR1
        IF (idum < 0) idum = idum + IM1
        IF (j <= NTAB) iv(j) = idum
      END DO
      iy = iv(1)
    END IF

    ! > Start here when not initializing. [cite: 7]
    ! > Compute idum=mod(IA1*idum,IM1) without overflows by Schrage’s method. [cite: 7]
    k = idum / IQ1
    idum = IA1 * (idum - k * IQ1) - k * IR1
    IF (idum < 0) idum = idum + IM1
    
    ! > Compute idum2=mod(IA2*idum2,IM2) likewise. [cite: 8]
    k = idum2 / IQ2
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2
    IF (idum2 < 0) idum2 = idum2 + IM2

    ! > Will be in the range 1:NTAB. [cite: 9]
    j = 1 + iy / NDIV
    
    ! > Here idum is shuffled, idum and idum2 are combined to generate output. [cite: 10]
    iy = iv(j) - idum2
    iv(j) = idum
    IF (iy < 1) iy = iy + IMM1

    ! > Because users don’t expect endpoint values. [cite: 11]
    random_val = MIN(AM * iy, RNMX)

  END FUNCTION ran2
  ! ===========================================================================


  ! ===========================================================================
  FUNCTION rnor(idum) RESULT(normal_deviate)
    ! Generates a normally distributed deviate with zero mean and unit variance
    ! using the Box-Muller method. It calls ran2 for uniform deviates.
    
    IMPLICIT NONE

    ! --- Arguments ---
    INTEGER, INTENT(INOUT) :: idum

    ! --- Function Result ---
    REAL(dp) :: normal_deviate

    ! --- Parameters ---
    REAL(dp), PARAMETER :: TWopi = 8.0_dp * ATAN(1.0_dp)

    ! --- Local Variables ---
    ! A saved value is used to return the second of two generated deviates
    ! on the next call, improving efficiency.
    INTEGER, SAVE :: iv_saved = 0
    REAL(dp), SAVE :: r_saved, t_saved
    REAL(dp) :: u1, u2

    ! Generate two new deviates if no value is saved or if re-initialization is requested.
    IF (iv_saved == 0 .OR. idum <= 0) THEN
      ! Get two uniform random numbers from ran2.
      u1 = ran2(idum)
      u2 = ran2(idum)
      
      r_saved = SQRT(-2.0_dp * LOG(u1))
      ! > = twopi*u2 [cite: 13]
      t_saved = TWopi * u2
      
      ! Return the first deviate.
      normal_deviate = r_saved * COS(t_saved)
      iv_saved = 1 ! Mark that a value is now saved.
    ELSE
      ! Return the saved deviate.
      normal_deviate = r_saved * SIN(t_saved)
      iv_saved = 0 ! Mark that the saved value has been used.
    END IF

  END FUNCTION rnor
  ! ===========================================================================

  ! ===========================================================================
  FUNCTION power(idum, alpha, t0) RESULT(power_val)
    ! Generates a random number according to a power-law distribution.
    ! The distribution is P(x) ~ (a-1)*T^(a-1)/(T+x)^a, where the input
    ! argument is alpha = -1/(a-1).
    ! Example: to get a power law a=2, alpha = -1/(2-1) = -1

    IMPLICIT NONE

    ! --- Arguments ---
    INTEGER, INTENT(INOUT) :: idum
    REAL(dp), INTENT(IN)   :: alpha, t0

    ! --- Function Result ---
    REAL(dp) :: power_val

    ! --- Local Variables ---
    REAL(dp) :: y

    ! The ran2 function is available directly as it's in the same module.
    y = ran2(idum)

    power_val = t0 * (y**alpha - 1.0_dp)

  END FUNCTION power
  ! ===========================================================================

  ! ===========================================================================
  FUNCTION levy(alpha, sigma, beta, mu, idum) RESULT(levy_val)
    ! Generates a random number from a stable Lévy distribution S_alpha(sigma, beta, mu).
    ! Based on the algorithm from R. Weron, Stat. & Prob. Letters 28 (1996) 165-171.

    IMPLICIT NONE

    ! --- Arguments ---
    REAL(dp), INTENT(IN)   :: alpha, sigma, beta, mu
    INTEGER, INTENT(INOUT) :: idum

    ! --- Function Result ---
    REAL(dp) :: levy_val

    ! --- Parameters ---
    REAL(dp), PARAMETER :: pihalf = 2.0_dp * ATAN(1.0_dp)

    ! --- Saved variables for one-time initialization ---
    LOGICAL, SAVE   :: is_initialized = .FALSE.
    REAL(dp), SAVE  :: balbe, salbe, scale, alphainv

    ! --- Local Variables ---
    REAL(dp) :: v, w

    ! One-time initialization block.
    IF (.NOT. is_initialized) THEN
      IF (alpha /= 1.0_dp) THEN
        balbe = ATAN(beta * TAN(pihalf * alpha)) / alpha
        salbe = (1.0_dp + (beta * TAN(pihalf * alpha))**2)**(0.5_dp / alpha)
      END IF
      ! Note: The 'scale' calculation from the source is specific to the alpha=1 case.
      scale = sigma * beta * LOG(sigma) / pihalf + mu
      alphainv = 1.0_dp / alpha
      is_initialized = .TRUE.
      ! The original code had a "print *" here, which is removed as it's
      ! not standard practice for a library function.
    END IF

    ! Main algorithm.
    v = pihalf - 2.0_dp * pihalf * ran2(idum)
    w = -LOG(ran2(idum))

    IF (alpha == 1.0_dp) THEN
      levy_val = ((pihalf + beta * v) * TAN(v) - &
                 beta * LOG(w * COS(v) / (pihalf + beta * v))) / pihalf
      levy_val = sigma * levy_val + scale
    ELSE
      levy_val = salbe * SIN(alpha * (v + balbe)) / (COS(v))**alphainv * &
                (COS(v - alpha * (v + balbe)) / w)**((1.0_dp - alpha) * alphainv)
      levy_val = sigma * levy_val + mu
    END IF

  END FUNCTION levy
  
  
END MODULE random_numbers_mod



