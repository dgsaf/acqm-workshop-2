!>
module laguerre

  ! debug compilation
  ! - <STDERR>: file unit for stderr output;
  ! - <DEBUG>: verbosity of debug statements;
  !   - 0: none;
  !   - 1: control flow (entering and exiting subroutines/functions);
  !   - 2: allocation, scalar assignment, short array assignment, etc;
  !   - 3: longer array assignment;
  !   - 4: internal variable assignment, allocation, inspecting loops, etc;
  ! - <PREFIX>: prefix every debug statement with this string.
  ! - <ERR>: prefix every error debug statement with this string.
#define STDERR 0
#define DEBUG 4
#define PREFIX "[debug] "
#define ERR "[error] "

  implicit none

  ! t_basis
  !
  ! Laguerre basis, consisting of basis states
  ! > {|phi_{i}>} for i = 1, ..., <n_basis>
  ! with coordinate-space representation
  ! > phi_{i}(r, theta, phi) = (varphi_{k_{i}, l_{i}}(r) / r)
  ! >                          * Y_{l_{i}, m_{i}}(theta, phi)
  ! where
  ! > varphi_{k, l}(r) = sqrt(alpha * (k - 1)! / (k + l) * (k + 2*l)!)
  ! >                    * (2*alpha*r)^{l+1}
  ! >                    * exp(-alpha*r)
  ! >                    * L_{k - 1}^{2*l + 1}(2*alpha*r)
  ! where L_{i}^{j} are the generalised Laguerre polynomials.
  !
  ! We construct the basis for given:
  ! - <m>: the magnetic quantum number, a conserved symmetry of the one-electron
  !   homonuclear-diatomic-molecule system;
  ! - <parity>: the parity quantum number, <parity> = (-1)^<l>, a conserved
  !   symmetry of the one-electron homonuclear-diatomic-molecule system;
  ! - <l_max>: maximum angular quantum number, <l>, considered in the basis;
  ! - <n_basis_l>: number of basis functions per <l> (from 0 to <l_max>);
  ! - <alpha_l>: value of <alpha> per <l> (from 0 to <l_max>).
  !
  ! As a result of this construction, we also have:
  ! - <n_basis>: the total number of basis states;
  ! - <k_list>: value of <k> for each basis state, k_{i} = k_list(i);
  ! - <l_list>: value of <l> for each basis state, l_{i} = l_list(i).
  !
  ! Furthermore, given a radial grid, <r_grid>, of length <n_r>, we also store:
  ! - <r_grid>: the radial grid;
  ! - <n_r>: the number of points in the radial grid;
  ! - <radial>: the radial basis functions, varphi_{k_{i}, l_{i}}(r), calculated
  !   on the radial grid points.
  !
  ! For the basis to be valid, we must have that:
  ! - <parity> = -1, or +1;
  ! - <l_max> > |<m>|, as we must have <m> in {-<l>, ..., <l>} for all basis
  !   states;
  ! - <n_basis_l(:)> >= 0, as we cannot have a negative amount of basis states
  !   for any value of <l>;
  ! - <alpha_l(:)> > 0.0, since alpha relates to radial variable, and must be
  !   strictly positive;
  ! - <n_basis> >= 1, since we require at least one basis state in a basis.
  !
  ! We note that the basis is ordered in the following way
  ! > {|phi_{1,m,m}>,     |phi_{2,m,m}>,   ...   |phi_{n_basis_l(m),m,m}>,
  ! >  |phi_{1,m+1,m}>,   |phi_{2,m+1,m}>, ...   |phi_{n_basis_l(m+1), m+1, m}>,
  ! >  ...,
  ! >  |phi_{1,l_max,m}>, |phi_{2,l_max,m}>, ... |phi_{n_basis_l(l_max),m+1,m}>}
  ! where
  ! > |phi_{i}> = |phi_{k_{i}, l_{i}, m_{i}}>.
  type t_basis
    integer :: m
    integer :: parity

    integer :: l_max
    integer , allocatable :: n_basis_l(:)
    double precision , allocatable :: alpha_l(:)

    integer :: n_basis
    integer , allocatable :: k_list(:)
    integer , allocatable :: l_list(:)

    integer :: n_r
    double precision , allocatable :: r_grid(:)
    double precision , allocatable :: radial(:, :)
  end type t_basis

contains

  ! setup_basis
  !
  ! For given <basis>, <m>, <parity>, <l_max>, <n_basis_l>, <alpha_l>,
  ! calculates the following:
  ! - <n_basis>;
  ! - <k_list>;
  ! - <l_list>.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0: indicates successful execution;
  ! - 1: indicates invalid arguments.
  subroutine setup_basis (basis, m, parity, l_max, n_basis_l, alpha_l, i_err)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: l_max
    integer , intent(in) :: n_basis_l(0:l_max)
    double precision , intent(in) :: alpha_l(0:l_max)
    integer , intent(in) :: m
    integer , intent(in) :: parity
    integer , intent(out) :: i_err
    integer :: ii, kk, ll

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "subroutine setup_basis()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (abs(parity) /= 1) then
      i_err = 1

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "abs(<parity>) /= 1"
#endif

    end if

    if (l_max < abs(m)) then
      i_err = 1

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<l_max> < abs(<m>)"
#endif

    end if

    if (l_max < 0) then
      i_err = 1

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<l_max> < 0"
#endif

    else
      if (any(n_basis_l(:) < 0)) then
        i_err = 1

#if (DEBUG >= 2)
        write (STDERR, *) PREFIX, ERR, "any(<n_basis_l(:)> < 0)"
#endif

      end if

      if (any(alpha_l(:) < 1.0D-8)) then
        i_err = 1

#if (DEBUG >= 2)
        write (STDERR, *) PREFIX, ERR, "any(<alpha_l(:)> < 1.0D-8)"
#endif

      end if
    end if

    ! handle invalid arguments
    if (i_err /= 0) then

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "arguments are invalid"
#endif

      ! set scalar variables to erroneous values
      basis%m = m
      basis%parity = parity
      basis%l_max = -1
      basis%n_basis = 0
      basis%n_r = 0

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, &
          "<basis> array variables will be left un-allocated"
      write (STDERR, *) PREFIX, ERR, &
          "<basis> scalar variables set to erroneous values"
#endif

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<basis%m> = ", basis%m
      write (STDERR, *) PREFIX, ERR, "<basis%parity> = ", basis%parity
      write (STDERR, *) PREFIX, ERR, "<basis%l_max> = ", basis%l_max
      write (STDERR, *) PREFIX, ERR, "<basis%n_basis> = ", basis%n_basis
      write (STDERR, *) PREFIX, ERR, "<basis%n_r> = ", basis%n_r
#endif

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "exiting subroutine setup_basis()"
#endif

      return
    end if

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "arguments are valid"
#endif

    ! set <m>, <parity>, <l_max>
    basis%m = m
    basis%parity = parity
    basis%l_max = l_max

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "<basis%m> = ", basis%m
    write (STDERR, *) PREFIX, "<basis%parity> = ", basis%parity
    write (STDERR, *) PREFIX, "<basis%l_max> = ", basis%l_max
#endif

    ! allocate <n_basis_l>, <alpha_l>
    allocate(basis%n_basis_l(0:l_max))
    allocate(basis%alpha_l(0:l_max))

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "allocated <basis%n_basis_l>"
    write (STDERR, *) PREFIX, "allocated <basis%alpha_l>"
#endif

    ! calculate <n_basis>, ignoring <n_basis_l> for <l> < |<m>|, and/or for
    ! (-1)**<l> /= <parity>
    basis%n_basis = 0

    do ll = 0, l_max

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<l> = ", ll
#endif

      if ((ll < abs(m)) .or. ((-1)**ll /= parity)) then
        basis%n_basis_l(ll) = 0

#if (DEBUG >= 4)
        write (STDERR, *) PREFIX, "ignoring due to basis symmetry"
#endif

      else
        basis%n_basis_l(ll) = n_basis_l(ll)
        basis%n_basis = basis%n_basis + n_basis_l(ll)
      end if
    end do

    ! set <alpha_l>
    basis%alpha_l(:) = alpha_l(:)

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "<l>, <basis%n_basis_l>, <basis%alpha_l>"
    do ll = 0, basis%l_max
      write (STDERR, *) PREFIX, ll, basis%n_basis_l(ll), basis%alpha_l(ll)
    end do
    write (STDERR, *) PREFIX, "<basis%n_basis> = ", basis%n_basis
#endif

    ! if basis is empty, deallocate memory and return an error code
    if (basis%n_basis == 0) then
      deallocate(basis%n_basis_l)
      deallocate(basis%alpha_l)

      i_err = 1

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "basis is empty"
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
      write (STDERR, *) PREFIX, ERR, "exiting subroutine setup_basis()"
#endif

      return
    end if

    ! allocate <k_list>, <l_list>
    allocate(basis%k_list(basis%n_basis))
    allocate(basis%l_list(basis%n_basis))

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "allocated <basis%k_list>"
    write (STDERR, *) PREFIX, "allocated <basis%l_list>"
#endif

    ! set <k_list>, <l_list>
    ii = 0

    do ll = 0, l_max
      if (basis%n_basis_l(ll) >= 1) then
        do kk = 1, basis%n_basis_l(ll)
          basis%k_list(ii+kk) = kk
          basis%l_list(ii+kk) = ll
        end do
      end if
      ii = ii + basis%n_basis_l(ll)
    end do

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "<i>, <basis%k_list>, <basis%l_list>"
    do ii = 1, basis%n_basis
      write (STDERR, *) PREFIX, ii, basis%k_list(ii), basis%l_list(ii)
    end do
#endif

    ! set <n_r> to zero, indicating radial basis functions have not yet been
    ! plotted on a grid
    basis%n_r = 0

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "<basis%n_r> = ", basis%n_r
#endif

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "end subroutine setup_basis()"
#endif

  end subroutine setup_basis

  ! is_valid
  !
  ! For given <basis>, determines if the basis variables are valid. A <basis> is
  ! invalid if:
  ! - abs(<basis%parity>) /= 1;
  ! - <basis%l_max> < abs(<basis%m>);
  ! - <basis%l_max> < 0;
  ! - <basis%n_basis> < 1;
  ! - <basis%l_max> < 0;
  ! - any(<basis%n_basis_l(:)> < 0);
  ! - any(<basis%alpha_l(:)> < 1.0D-8).
  logical function is_valid (basis)
    type(t_basis) , intent(inout) :: basis

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "function is_valid()"
#endif

    ! check if <basis> is valid
    is_valid = .true.

    if (basis%n_basis < 1) then
      is_valid = .false.

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<basis%n_basis> < 1"
#endif

    end if

    if (abs(basis%parity) /= 1) then
      is_valid = .false.

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "abs(<basis%parity>) /= 1"
#endif

    end if

    if (basis%l_max < abs(basis%m)) then
      is_valid = .false.

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<basis%l_max> < abs(<basis%m>)"
#endif

    end if

    if (basis%l_max < 0) then
      is_valid = .false.

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<basis%l_max> < 0"
#endif

    else
      if (any(basis%n_basis_l(:) < 0)) then
        is_valid = .false.

#if (DEBUG >= 2)
        write (STDERR, *) PREFIX, ERR, "any(<basis%n_basis_l(:)> < 0)"
#endif

      end if

      if (any(basis%alpha_l(:) < 1.0D-8)) then
        is_valid = .false.

#if (DEBUG >= 2)
        write (STDERR, *) PREFIX, ERR, "any(<basis%alpha_l(:)> < 1.0D-8)"
#endif

      end if
    end if

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "end function is_valid()"
#endif

    return
  end function is_valid

  ! setup_radial
  !
  ! For given <basis>, <n_r>, <r_grid>, calculates the radial functions
  ! > varphi_{k_{i}, l_{i}}(r) for i = 1, ..., <n_basis>
  ! on the radial values specified in the grid.
  !
  ! Requires the following variables in <basis> to have already been setup:
  ! - <l_max>;
  ! - <n_basis_l>;
  ! - <alpha_l>;
  ! - <n_basis>.
  ! That is, a call to setup_basis() should have already been made.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine setup_radial (basis, n_r, r_grid, i_err)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    integer , intent(out) :: i_err
    double precision , allocatable :: norm(:)
    double precision :: alpha_grid(n_r)
    double precision :: alpha
    integer :: n_b_l, offset
    integer :: jj, kk, ll

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "subroutine setup_radial()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. is_valid(basis)) then
      i_err = 1
    end if

    if (n_r < 1) then
      i_err = 1

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<n_r> < 1"
#endif

    end if

    ! handle invalid arguments
    if (i_err /= 0) then

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "arguments are invalid"
      write (STDERR, *) PREFIX, ERR, &
          "<basis> array variables will be left un-allocated"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine setup_radial()"
#endif

      return
    end if

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "arguments are valid"
#endif

    ! set <n_r>
    basis%n_r = n_r

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "<basis%n_r> = ", basis%n_r
#endif

    ! allocate <r_grid>, <radial>
    allocate(basis%r_grid(n_r))
    allocate(basis%radial(n_r, basis%n_basis))

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "allocated <basis%r_grid>"
    write (STDERR, *) PREFIX, "allocated <basis%radial>"
#endif

    ! basis: set <r_grid>
    basis%r_grid(:) = r_grid(:)

#if (DEBUG >= 3)
    write (STDERR, *) PREFIX, "<radial index>, <basis%r_grid>"
    do jj = 1, basis%n_r
      write (STDERR, *) PREFIX, jj, basis%r_grid(jj)
    end do
#endif

    ! basis function index offset, incremented by n_basis_l(ll) after each loop
    offset = 0

    ! loop over <l>, basis: set <radial>
    do ll = 0, basis%l_max

#if (DEBUG >= 3)
      write (STDERR, *) PREFIX, "<l> = ", ll
      write (STDERR, *) PREFIX, "<offset> = ", offset
#endif

      ! in-line <n_basis_l>, <alpha> for current <l>
      n_b_l = basis%n_basis_l(ll)
      alpha = basis%alpha_l(ll)

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<n_b_l> = ", n_b_l
      write (STDERR, *) PREFIX, "<alpha> = ", alpha
#endif

      ! if there are no basis states for this value of <l>, then cycle
      if (n_b_l == 0) then

#if (DEBUG >= 4)
        write (STDERR, *) PREFIX, "no basis states for this value of <l>"
        write (STDERR, *) PREFIX, "cycling"
#endif

        cycle
      end if

      ! allocate normalisation constant array
      allocate(norm(n_b_l))

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "allocated <norm>"
#endif

      ! recurrence relation for basis normalisation constants
#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<k>, <norm>"
#endif

      if (n_b_l >= 1) then
        norm(1) = sqrt(alpha / dble((ll+1) * gamma(dble((2*ll)+2))))

#if (DEBUG >= 4)
        write (STDERR, *) PREFIX, 1, norm(1)
#endif

      end if

      if (n_b_l >= 2) then
        do kk = 2, n_b_l
          norm(kk) = norm(kk-1) * sqrt(dble((kk-1) * (kk-1+ll)) &
              / dble((kk+ll) * (kk+(2*ll))))

#if (DEBUG >= 4)
          write (STDERR, *) PREFIX, kk, norm(kk)
#endif

        end do
      end if

      ! in-line <alpha_grid> for current <l>
      alpha_grid(:) = alpha * r_grid(:)

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<radial index>, <alpha_grid>"
      do jj = 1, basis%n_r
        write (STDERR, *) PREFIX, jj, alpha_grid(jj)
      end do
#endif

      ! recurrence relation for basis functions
      if (n_b_l >= 1) then
        basis%radial(:, offset+1) = ((2.0d0 * alpha_grid(:)) ** (ll+1)) &
            * exp(-alpha_grid(:))
      end if

      if (n_b_l >= 2) then
        basis%radial(:, offset+2) = 2.0d0 * (dble(ll+1) - alpha_grid(:)) &
            * basis%radial(:, offset+1)
      end if

      if (n_b_l >= 3) then
        do kk = 3, n_b_l
          basis%radial(:, offset+kk) = &
              ((2.0d0 * (dble(kk-1+ll) - alpha_grid(:)) &
              * basis%radial(:, offset+kk-1)) &
              - dble(kk+(2*ll)-1) * basis%radial(:, offset+kk-2)) &
              / dble(kk-1)
        end do
      end if

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<radial index>, <basis%radial> (un-normalised)"
      do jj = 1, basis%n_r
        write (STDERR, *) PREFIX, jj, basis%radial(jj, offset+1:offset+n_b_l)
      end do
#endif

      ! scale basis radial functions by normalisation constants
      if (n_b_l >= 1) then
        do kk = 1, n_b_l
          basis%radial(:, offset+kk) = basis%radial(:, offset+kk) * norm(kk)
        end do
      end if

#if (DEBUG >= 3)
      write (STDERR, *) PREFIX, "<radial index>, <basis%radial>"
      do jj = 1, basis%n_r
        write (STDERR, *) PREFIX, jj, basis%radial(jj, offset+1:offset+n_b_l)
      end do
#endif

      ! increment basis index offset
      offset = offset + n_b_l

      ! deallocate norm
      deallocate(norm)

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "deallocated <norm>"
#endif

    end do

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "end subroutine setup_basis()"
#endif

  end subroutine setup_radial

  ! overlap
  !
  ! For given <basis>, calculate the overlap matrix elements
  ! > B_{i, j} = < phi_{i} | phi_{j} > for i, j = 1, ..., <n_basis>
  ! using analytic properties of the Laguerre basis.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine overlap (basis, B, i_err)
    type(t_basis) , intent(inout) :: basis
    double precision , intent(out) :: B(basis%n_basis, basis%n_basis)
    integer , intent(out) :: i_err
    integer :: n_b_l, offset
    integer :: kk, ll

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "subroutine overlap()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. is_valid(basis)) then
      i_err = 1
    end if

    ! handle invalid arguments
    if (i_err /= 0) then

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "arguments are invalid"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine overlap()"
#endif

      return
    end if

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "arguments are valid"
#endif

    ! initialise <B>
    B(:, :) = 0.0d0

    ! calculate tri-diagonal overlap matrix elements
    offset = 0

    do ll = 0, basis%l_max

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<l> = ", ll
      write (STDERR, *) PREFIX, "<offset> = ", offset
#endif

      n_b_l = basis%n_basis_l(ll)

      ! if there are no basis states for this value of <l>, then cycle
      if (n_b_l == 0) then

#if (DEBUG >= 4)
        write (STDERR, *) PREFIX, "no basis states for this value of <l>"
        write (STDERR, *) PREFIX, "cycling"
#endif

        cycle
      end if

      ! calculate tri-diagonal overlap matrix elements for current <l>
      if (n_b_l >= 1) then

#if (DEBUG >= 3)
        write (STDERR, *) PREFIX, "<i>, B(i, i), B(i+1, i) = B(i, i+1)"
#endif

        if (n_b_l >= 2) then
          do kk = 1, n_b_l-1
            B(offset+kk, offset+kk) = 1.0d0

            B(offset+kk, offset+kk+1) = - 0.5d0 * sqrt(1 - &
                (dble(ll * (ll + 1)) / dble((kk + ll) * (kk + ll + 1))))

            B(offset+kk+1, offset+kk) = B(offset+kk, offset+kk+1)

#if (DEBUG >= 3)
            write (STDERR, *) PREFIX, &
                offset+kk, B(offset+kk, offset+kk), B(offset+kk+1, offset+kk)
#endif

          end do
        end if

        ! last term (not covered by loop)
        B(offset+n_b_l, offset+n_b_l) = 1.0d0

#if (DEBUG >= 3)
        write (STDERR, *) PREFIX, &
            offset+n_b_l, B(offset+n_b_l, offset+n_b_l)
#endif

      end if

      ! increment basis index offset
      offset = offset + n_b_l
    end do

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "end subroutine overlap()"
#endif

  end subroutine overlap

  ! kinetic
  !
  ! For given <basis>, calculate the overlap matrix elements
  ! > K_{i, j} = < phi_{i} | K | phi_{j} > for i, j = 1, ..., <n_basis>
  ! using analytic properties of the Laguerre basis.
  !
  ! Returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine kinetic (basis, K, i_err)
    type(t_basis) , intent(inout) :: basis
    double precision , intent(out) :: K(basis%n_basis, basis%n_basis)
    integer , intent(out) :: i_err
    double precision :: alpha
    integer :: n_b_l, offset
    integer :: kk, ll

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "subroutine kinetic()"
#endif

    ! check if arguments are valid
    i_err = 0

    if (.not. is_valid(basis)) then
      i_err = 1
    end if

    ! handle invalid arguments
    if (i_err /= 0) then

#if (DEBUG >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG >= 1)
      write (STDERR, *) PREFIX, ERR, "arguments are invalid"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine kinetic()"
#endif

      return
    end if

#if (DEBUG >= 2)
    write (STDERR, *) PREFIX, "arguments are valid"
#endif

    ! initialise <K>
    K(:, :) = 0.0d0

    ! calculate tri-diagonal kinetic matrix elements
    offset = 0

    do ll = 0, basis%l_max

#if (DEBUG >= 4)
      write (STDERR, *) PREFIX, "<l> = ", ll
      write (STDERR, *) PREFIX, "<offset> = ", offset
#endif

      n_b_l = basis%n_basis_l(ll)
      alpha = basis%alpha_l(ll)

      ! if there are no basis states for this value of <l>, then cycle
      if (n_b_l == 0) then

#if (DEBUG >= 4)
        write (STDERR, *) PREFIX, "no basis states for this value of <l>"
        write (STDERR, *) PREFIX, "cycling"
#endif

        cycle
      end if

      ! calculate tri-diagonal kinetic matrix elements for current <l>
      if (n_b_l >= 1) then

#if (DEBUG >= 3)
        write (STDERR, *) PREFIX, "<i>, K(i, i), K(i+1, i) = K(i, i+1)"
#endif

        if (n_b_l >= 2) then
          do kk = 1, n_b_l-1
            K(offset+kk, offset+kk) = 0.5d0 * (alpha ** 2)

            K(offset+kk, offset+kk+1) = (alpha ** 2) * 0.25d0 * sqrt(1 - &
                (dble(ll * (ll + 1)) / dble((kk + ll) * (kk + ll + 1))))

            K(offset+kk+1, offset+kk) = K(offset+kk, offset+kk+1)

#if (DEBUG >= 3)
            write (STDERR, *) PREFIX, &
                offset+kk, K(offset+kk, offset+kk), K(offset+kk+1, offset+kk)
#endif

          end do
        end if

        ! last term (not covered by loop)
        K(offset+n_b_l, offset+n_b_l) = 0.5d0 * (alpha ** 2)

#if (DEBUG >= 3)
        write (STDERR, *) PREFIX, &
            offset+n_b_l, K(offset+n_b_l, offset+n_b_l)
#endif

      end if

      ! increment basis index offset
      offset = offset + n_b_l
    end do

#if (DEBUG >= 1)
    write (STDERR, *) PREFIX, "end subroutine kinetic()"
#endif

  end subroutine kinetic

end module laguerre
