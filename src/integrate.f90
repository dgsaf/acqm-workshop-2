!>
module integrate

  ! debug compilation
  ! - <STDERR>: file unit for stderr output;
  ! - <DEBUG_INTEGRATE>: verbosity of debug statements;
  !   - 0: none;
  !   - 1: control flow (entering and exiting subroutines/functions);
  !   - 2: allocation, scalar assignment, short array assignment, etc;
  !   - 3: longer array assignment;
  !   - 4: internal variable assignment, allocation, inspecting loops, etc;
  ! - <PREFIX>: prefix every debug statement with this string;
  ! - <ERR>: prefix every error debug statement with this string.
#define STDERR 0
#define DEBUG_INTEGRATE 0
#define PREFIX "[debug] "
#define ERR "[error] "
#define TOL 1.0D-10

  implicit none

contains

  function integrate_trapezoid (n, x, f) result (integral)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n), f(n)
    double precision :: w(n)
    double precision :: integral
    integer :: ii
    integer :: i_err

#if (DEBUG_INTEGRATE >= 1)
    write (STDERR, *) PREFIX, "function integrate_trapezoid()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (n < 2) then
      i_err = 1

#if (DEBUG_INTEGRATE >= 2)
      write (STDERR, *) PREFIX, ERR, "<n> < 2"
#endif

    end if

    ! handle invalid arguments
    if (i_err /= 0) then

#if (DEBUG_INTEGRATE >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

      integral = 0.0d0

#if (DEBUG_INTEGRATE >= 1)
      write (STDERR, *) PREFIX, ERR, "arguments are invalid"
      write (STDERR, *) PREFIX, ERR, "returning <integral> = 0"
      write (STDERR, *) PREFIX, ERR, "exiting function integrate_trapezoid()"
#endif

      return
    end if

#if (DEBUG_INTEGRATE >= 2)
    write (STDERR, *) PREFIX, "arguments are valid"
#endif

#if (DEBUG_INTEGRATE >= 4)
    write (STDERR, *) PREFIX, "<index>, <w>"
#endif

    ! calculate weights
    w(1) = (x(2) - x(1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 4)
    write (STDERR, *) PREFIX, 1, w(1)
#endif

    do ii = 2, n-1
      w(ii) = (x(ii+1) - x(ii-1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 4)
      write (STDERR, *) PREFIX, ii, w(ii)
#endif

    end do

    w(n) = (x(n) - x(n-1)) / 2.0d0

#if (DEBUG_INTEGRATE >= 4)
    write (STDERR, *) PREFIX, n, w(n)
#endif

    ! calculate integral
    integral = 0.0d0

    do ii = 1, n
      integral = integral + (w(ii) * f(ii))
    end do

#if (DEBUG_INTEGRATE >= 2)
    write (STDERR, *) PREFIX, "<integral> = ", integral
#endif

#if (DEBUG_INTEGRATE >= 1)
    write (STDERR, *) PREFIX, "end function integrate_trapezoid()"
#endif

    return

  end function integrate_trapezoid

end module integrate
