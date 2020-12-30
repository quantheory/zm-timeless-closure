module zm_timeless_closure_ffi

  ! Interface to call zm_timeless_closure from other languages such as C or
  ! Python.
  !
  ! Original author: Sean Patrick Santos <SeanPatrickSantos@gmail.com>

  use, intrinsic :: iso_c_binding

  use zm_timeless_closure

  implicit none

contains

#ifdef TEST_MODE

  ! Weighting function for CAPE relaxation (versus large-scale CAPE tendency).
  !
  ! Precondition: x > 0
  pure subroutine weight_ffi(x, w) bind(c, name="zmtc_weight")
    ! x = deltat / tau
    real(8), intent(in) :: x
    real(8), intent(out) :: w
    w = weight(x)
  end subroutine weight_ffi

  ! Vector version of above function.
  pure subroutine weight_1d_ffi(x, w) bind(c, name="zmtc_weight_1d")
    ! x = deltat / tau
    real(8), intent(in) :: x(:)
    real(8), intent(out) :: w(size(x))
    w = weight(x)
  end subroutine weight_1d_ffi

  ! Rate of CAPE consumption when convection is ongoing.
  !
  ! This is used when a_crit < a_prev < a_star, or when a_crit < a_star < a_prev
  ! and convection is confirmed to last the entire time step.
  !
  ! Precondition: tau > 0 and deltat > 0
  pure subroutine cape_consumption_ongoing_ffi(tau, a_crit, deltat, weight, &
       a_prev, a_star, dadtc) bind(c, name="zmtc_cape_consumption_ongoing")
    ! Convective relaxation timescale (s)
    real(8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(8), intent(in) :: deltat
    ! Weighting of relaxation term vs. large-scale tendency terms
    real(8), intent(in) :: weight
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(8), intent(out) :: dadtc(size(a_prev))

    dadtc = cape_consumption_ongoing(tau, a_crit, deltat, weight, a_prev, &
         a_star)
  end subroutine cape_consumption_ongoing_ffi

  ! Rate of CAPE consumption when convection starts mid-timestep.
  !
  ! This is used when a_prev < a_crit < a_star, guaranteeing that convection
  ! will begin during this time step.
  !
  ! Precondition: tau > 0, deltat > 0, a_star > a_crit, and a_star > a_prev
  pure subroutine cape_consumption_starting_ffi(tau, a_crit, deltat, a_prev, &
       a_star, dadtc) bind(c, name="zmtc_cape_consumption_starting")
    ! Convective relaxation timescale (s)
    real(8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(8), intent(out) :: dadtc(size(a_prev))

    dadtc = cape_consumption_starting(tau, a_crit, deltat, a_prev, a_star)
  end subroutine cape_consumption_starting_ffi

  ! Fraction of the time step when CAPE will be depleted.
  !
  ! This is used when a_crit < a_prev and a_star < a_prev, to determine the
  ! point at which CAPE will be completely depleted, as a fraction of the host
  ! model coupling time step. If f >= 1, convection will last the entire time
  ! step.
  !
  ! Preconditions: tau > 0, deltat > 0, a_crit < a_prev, and a_star < a_prev
  pure subroutine end_time_frac_ffi(tau, a_crit, deltat, a_prev, a_star, f) &
       bind(c, name="zmtc_end_time_frac")
    ! Convective relaxation timescale (s)
    real(8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(8), intent(in) :: a_star(:)
    ! Fraction of time step when convection occurs
    real(8), intent(out) :: f(size(a_prev))

    f = end_time_frac(tau, a_crit, deltat, a_prev, a_star)
  end subroutine end_time_frac_ffi

  ! Rate of CAPE consumption when convection ends mid-timestep.
  !
  ! This is used when a_crit < a_prev, a_star < a_prev, and convection has been
  ! confirmed not to last the entire time step.
  !
  ! Precondition: deltat > 0, a_crit < a_prev, and a_star < a_prev
  pure subroutine cape_consumption_ending_ffi(a_crit, deltat, f, a_prev, &
       a_star, dadtc) bind(c, name="zmtc_cape_consumption_ending")
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(8), intent(in) :: deltat
    ! Fraction of time step when convection occurs
    real(8), intent(in) :: f(:)
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(8), intent(out) :: dadtc(size(a_prev))

    dadtc = cape_consumption_ending(a_crit, deltat, f, a_prev, a_star)
  end subroutine cape_consumption_ending_ffi

#endif
  !endif TEST_MODE

  ! Combined formula for rate of CAPE consumption.
  !
  ! This is the function that should actually be called by the ZM deep
  ! convection. It detects whether consumption is starting, ending, or persists
  ! throughout an entire time step, and calls the corresponding function.
  !
  ! Precondition: tau > 0 and deltat > 0
  pure subroutine cape_consumption_rate_ffi(tau, a_crit, deltat, a_prev, &
       a_star, dadtc) bind(c, name="zmtc_cape_consumption_rate")
    ! Convective relaxation timescale (s)
    real(8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(8), intent(in) :: deltat
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(8), intent(out) :: dadtc(size(a_prev))

    dadtc = cape_consumption_rate(tau, a_crit, deltat, a_prev, a_star)
  end subroutine cape_consumption_rate_ffi

end module zm_timeless_closure_ffi
