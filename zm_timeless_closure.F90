module zm_timeless_closure

  ! An alternate closure condition for the Zhang-McFarlane deep convection.
  !
  ! Original author: Sean Patrick Santos <SeanPatrickSantos@gmail.com>

  use shr_kind_mod

  implicit none
#ifdef TEST_MODE
  public
#else
  private
#endif

  interface weight
     module procedure weight
     module procedure weight_1d
  end interface weight

contains

  ! Weighting function for CAPE relaxation (versus large-scale CAPE tendency).
  ! Precondition: x > 0
  pure function weight(x) result(w)
    ! x = deltat / tau
    real(shr_kind_r8), intent(in) :: x
    real(shr_kind_r8) :: w
    w = (1._shr_kind_r8 - exp(-x)) / x
  end function weight

  ! Vector version of above function.
  pure function weight_1d(x) result(w)
    ! x = deltat / tau
    real(shr_kind_r8), intent(in) :: x(:)
    real(shr_kind_r8) :: w(size(x))
    w = (1._shr_kind_r8 - exp(-x)) / x
  end function weight_1d

  ! Rate of CAPE consumption when convection is ongoing.
  ! Specifically, this is used when a_crit < a_prev < a_star.
  pure function cape_consumption_ongoing(tau, a_crit, deltat, weight, a_prev, &
       a_star) result(dadtc)
    ! Convective relaxation timescale (s)
    real(shr_kind_r8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! Weighting of relaxation term vs. large-scale tendency terms
    real(shr_kind_r8), intent(in) :: weight
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc(size(a_prev))

    real(shr_kind_r8) :: cape_relaxation(size(a_prev))
    real(shr_kind_r8) :: dadtls(size(a_prev))

    cape_relaxation = (a_prev - a_crit) / tau
    dadtls = (a_star - a_prev) / deltat
    dadtc = cape_relaxation*weight + dadtls*(1._shr_kind_r8 - weight)
  end function cape_consumption_ongoing

end module zm_timeless_closure
