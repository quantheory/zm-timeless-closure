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

  public cape_consumption_rate

  interface weight
     module procedure weight
     module procedure weight_1d
  end interface weight

contains

  ! Weighting function for CAPE relaxation (versus large-scale CAPE tendency).
  !
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
  !
  ! This is used when a_crit < a_prev < a_star, or when a_crit < a_star < a_prev
  ! and convection is confirmed to last the entire time step.
  !
  ! Precondition: tau > 0 and deltat > 0
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

  ! Rate of CAPE consumption when convection starts mid-timestep.
  !
  ! This is used when a_prev < a_crit < a_star, guaranteeing that convection
  ! will begin during this time step.
  !
  ! Precondition: tau > 0, deltat > 0, a_star > a_crit, and a_star > a_prev
  pure function cape_consumption_starting(tau, a_crit, deltat, a_prev, a_star) &
       result(dadtc)
    ! Convective relaxation timescale (s)
    real(shr_kind_r8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc(size(a_prev))

    ! "1 - f" is the fraction of the time step when convection occurs.
    real(shr_kind_r8) :: omf(size(a_prev))
    real(shr_kind_r8) :: weights(size(a_prev))

    omf = (a_star - a_crit) / (a_star - a_prev)
    weights = weight(omf * deltat / tau)
    dadtc = (a_star - a_crit) * (1._shr_kind_r8 - weights) / deltat
  end function cape_consumption_starting

  ! Fraction of the time step when CAPE will be depleted.
  !
  ! This is used when a_crit < a_prev and a_star < a_prev, to determine the
  ! point at which CAPE will be completely depleted, as a fraction of the host
  ! model coupling time step. If f >= 1, convection will last the entire time
  ! step.
  !
  ! Preconditions: tau > 0, deltat > 0, a_crit < a_prev, and a_star < a_prev
  pure function end_time_frac(tau, a_crit, deltat, a_prev, a_star) result(f)
    ! Convective relaxation timescale (s)
    real(shr_kind_r8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Fraction of time step when convection occurs
    real(shr_kind_r8) :: f(size(a_prev))

    f = log(1 + (a_prev - a_crit) * deltat / ((a_prev - a_star) * tau)) * &
         tau / deltat
  end function end_time_frac

  ! Rate of CAPE consumption when convection ends mid-timestep.
  !
  ! This is used when a_crit < a_prev, a_star < a_prev, and convection has been
  ! confirmed not to last the entire time step.
  !
  ! Precondition: deltat > 0, a_crit < a_prev, and a_star < a_prev
  pure function cape_consumption_ending(a_crit, deltat, f, a_prev, a_star) &
       result(dadtc)
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! Fraction of time step when convection occurs
    real(shr_kind_r8), intent(in) :: f(:)
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc(size(a_prev))

    dadtc = (f * a_star + (1._shr_kind_r8 - f) * a_prev - a_crit) / deltat
  end function cape_consumption_ending

  ! Combined formula for rate of CAPE consumption.
  !
  ! This is the function that should actually be called by the ZM deep
  ! convection. It detects whether consumption is starting, ending, or persists
  ! throughout an entire time step, and calls the corresponding function.
  !
  ! Precondition: tau > 0 and deltat > 0
  pure function cape_consumption_rate(tau, a_crit, deltat, a_prev, a_star) &
       result(dadtc)
    ! Convective relaxation timescale (s)
    real(shr_kind_r8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc(size(a_prev))

    ! Local copy of averaging weight.
    real(shr_kind_r8) :: weight_loc
    ! Fraction of time step with convection.
    real(shr_kind_r8) :: f(size(a_prev))

    weight_loc = weight(deltat/tau)

    where (a_prev <= a_crit)
       where (a_star <= a_crit)
          dadtc = 0._shr_kind_r8
       elsewhere
          dadtc = cape_consumption_starting(tau, a_crit, deltat, a_prev, a_star)
       end where
    elsewhere
       where (a_prev <= a_star)
          f = 1.
       elsewhere
          f = end_time_frac(tau, a_crit, deltat, a_prev, a_star)
       end where
       where (f >= 1.)
          dadtc = cape_consumption_ongoing(tau, a_crit, deltat, weight_loc, &
               a_prev, a_star)
       elsewhere
          dadtc = cape_consumption_ending(a_crit, deltat, f, a_prev, a_star)
       end where
    end where
  end function cape_consumption_rate

end module zm_timeless_closure
