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

  public new_zmtc
  public zmtc_t
  public cape_consumption_rate

  ! Structure containing closure-related settings.
  type zmtc_t
    ! Convective relaxation timescale (s)
    real(shr_kind_r8) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8) :: deltat
    ! Weighting of relaxation term vs. large-scale tendency terms
    real(shr_kind_r8) :: weight
  end type zmtc_t

contains

  ! Constructor for zmtc setting struct.
  !
  ! Precondition for all functions: tau > 0 and deltat > 0
  function new_zmtc(tau, a_crit, deltat) result(zmtc)
    ! Convective relaxation timescale (s)
    real(shr_kind_r8), intent(in) :: tau
    ! Threshold CAPE needed for convection to trigger (J/kg)
    real(shr_kind_r8), intent(in) :: a_crit
    ! Coupling time step to host model (s)
    real(shr_kind_r8), intent(in) :: deltat
    ! New ZM timeless closure option struct
    type(zmtc_t) :: zmtc

    zmtc%tau = tau
    zmtc%a_crit = a_crit
    zmtc%deltat = deltat
    zmtc%weight = weight(deltat/tau)
  end function new_zmtc

  ! Weighting function for CAPE relaxation (versus large-scale CAPE tendency).
  !
  ! Precondition: dt_ratio > 0
  elemental function weight(dt_ratio) result(w)
    ! Ratio of host model time step to relaxation timescale.
    real(shr_kind_r8), intent(in) :: dt_ratio
    real(shr_kind_r8) :: w
    w = (1._shr_kind_r8 - exp(-dt_ratio)) / dt_ratio
  end function weight

  ! Rate of CAPE consumption when convection is ongoing.
  !
  ! This is used when a_crit < a_prev < a_star, or when a_crit < a_star < a_prev
  ! and convection is confirmed to last the entire time step.
  elemental function cape_consumption_ongoing(zmtc, a_prev, a_star) result(dadtc)
    ! ZM timeless closure settings struct.
    type(zmtc_t), intent(in) :: zmtc
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc

    real(shr_kind_r8) :: cape_relaxation
    real(shr_kind_r8) :: dadtls

    cape_relaxation = (a_prev - zmtc%a_crit) / zmtc%tau
    dadtls = (a_star - a_prev) / zmtc%deltat
    dadtc = cape_relaxation * zmtc%weight + &
         dadtls * (1._shr_kind_r8 - zmtc%weight)
  end function cape_consumption_ongoing

  ! Rate of CAPE consumption when convection starts mid-timestep.
  !
  ! This is used when a_prev < a_crit < a_star, guaranteeing that convection
  ! will begin during this time step.
  !
  ! Precondition: a_star > a_crit and a_star > a_prev
  elemental function cape_consumption_starting(zmtc, a_prev, a_star) &
       result(dadtc)
    ! ZM timeless closure settings struct.
    type(zmtc_t), intent(in) :: zmtc
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc

    ! "1 - f" is the fraction of the time step when convection occurs.
    real(shr_kind_r8) :: omf
    real(shr_kind_r8) :: weight_omf

    omf = (a_star - zmtc%a_crit) / (a_star - a_prev)
    weight_omf = weight(omf * zmtc%deltat / zmtc%tau)
    dadtc = (a_star - zmtc%a_crit) * (1._shr_kind_r8 - weight_omf) / zmtc%deltat
  end function cape_consumption_starting

  ! Fraction of the time step when CAPE will be depleted.
  !
  ! This is used when a_crit < a_prev and a_star < a_prev, to determine the
  ! point at which CAPE will be completely depleted, as a fraction of the host
  ! model coupling time step. If f >= 1, convection will last the entire time
  ! step.
  !
  ! Preconditions: a_crit < a_prev and a_star < a_prev
  elemental function end_time_frac(zmtc, a_prev, a_star) result(f)
    ! ZM timeless closure settings struct.
    type(zmtc_t), intent(in) :: zmtc
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star
    ! Fraction of time step when convection occurs
    real(shr_kind_r8) :: f

    real(shr_kind_r8) :: dt_ratio

    dt_ratio = zmtc%deltat / zmtc%tau
    ! First order approximation to f.
    f = (a_prev - zmtc%a_crit) / (a_prev - a_star)
    ! Correct to get exact fraction.
    f = log(1 + f * dt_ratio) / dt_ratio
  end function end_time_frac

  ! Rate of CAPE consumption when convection ends mid-timestep.
  !
  ! This is used when a_crit < a_prev, a_star < a_prev, and convection has been
  ! confirmed not to last the entire time step.
  !
  ! Precondition: a_crit < a_prev and a_star < a_prev
  elemental function cape_consumption_ending(zmtc, f, a_prev, a_star) &
       result(dadtc)
    ! ZM timeless closure settings struct.
    type(zmtc_t), intent(in) :: zmtc
    ! Fraction of time step when convection occurs
    real(shr_kind_r8), intent(in) :: f
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc

    dadtc = (f * a_star + (1._shr_kind_r8 - f) * a_prev - zmtc%a_crit) / &
         zmtc%deltat
  end function cape_consumption_ending

  ! Combined formula for rate of CAPE consumption.
  !
  ! This is the function that should actually be called by the ZM deep
  ! convection. It detects whether consumption is starting, ending, or persists
  ! throughout an entire time step, and calls the corresponding function.
  pure function cape_consumption_rate(zmtc, a_prev, a_star) result(dadtc)
    ! ZM timeless closure settings struct.
    type(zmtc_t), intent(in) :: zmtc
    ! Final CAPE after convection occurred in the previous time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_prev(:)
    ! CAPE after large-scale processes act in this time step (J/kg)
    real(shr_kind_r8), intent(in) :: a_star(:)
    ! Magnitude of CAPE consumption by convection (result is positive) (J/kg/s)
    real(shr_kind_r8) :: dadtc(size(a_prev))

    ! Fraction of time step with convection.
    real(shr_kind_r8) :: f(size(a_prev))

    where (a_prev <= zmtc%a_crit)
       where (a_star <= zmtc%a_crit)
          dadtc = 0._shr_kind_r8
       elsewhere
          dadtc = cape_consumption_starting(zmtc, a_prev, a_star)
       end where
    elsewhere
       where (a_prev <= a_star)
          f = 1.
       elsewhere
          f = end_time_frac(zmtc, a_prev, a_star)
       end where
       where (f >= 1.)
          dadtc = cape_consumption_ongoing(zmtc, a_prev, a_star)
       elsewhere
          dadtc = cape_consumption_ending(zmtc, f, a_prev, a_star)
       end where
    end where
  end function cape_consumption_rate

end module zm_timeless_closure
