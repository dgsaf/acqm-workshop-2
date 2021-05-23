!>
module laguerre

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
  ! - <n_basis_l>: number of basis functions per <l>;
  ! - <alpha_l>: value of <alpha> per <l>.
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
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine setup_basis (basis, m, parity, l_max, n_basis_l, alpha_l, i_err)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: l_max
    integer , intent(in) :: n_basis_l(l_max)
    integer , intent(in) :: alpha_l(l_max)
    integer , intent(in) :: m
    integer , intent(in) :: parity
    integer , intent(out) :: i_err
    integer :: ii, kk, ll

    ! check if arguments are valid
    i_err = 0

    if (.not. (abs(parity) == 1) &
        .or. (l_max < abs(m)) &
        .or. (any(n_basis_l(:) < 0)) &
        .or. (any(basis%alpha_l(:) < 1.0D-8))) then
      i_err = 1
      return
    end if

    ! set <m>, <parity>, <l_max>
    basis%m = m
    basis%parity = parity
    basis%l_max = l_max

    ! allocate <n_basis_l>, <alpha_l>
    allocate(basis%n_basis_l(l_max))
    allocate(basis%alpha_l(l_max))

    ! set <alpha_l>
    basis%alpha_l(:) = alpha_l(:)

    ! calculate <n_basis>, and ignore <n_basis_l> for <ll> < |<m>|
    basis%n_basis = 0

    do ll = 0, l_max
      if (ll < abs(m)) then
        basis%n_basis_l(ll) = 0
      else
        basis%n_basis_l(ll) = n_basis_l(ll)
        basis%n_basis = basis%n_basis + n_basis_l(ll)
      end if
    end do

    ! if basis is empty, deallocate memory and return an error code
    if (basis%n_basis == 0) then
      deallocate(basis%n_basis_l)
      deallocate(basis%alpha_l)

      i_err = 1
      return
    end if

    ! allocate <k_list>, <l_list>
    allocate(basis%k_list(basis%n_basis))
    allocate(basis%l_list(basis%n_basis))

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

  end subroutine setup_basis

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
    integer , intent(in) :: r_grid(n_r)
    integer , intent(out) :: i_err
    double precision , allocatable :: norm(:)
    double precision :: alpha_grid(n_r)
    double precision :: alpha
    integer :: n_b_l
    integer :: ii, kk, ll

    ! check if arguments are valid
    i_err = 0

    if ((basis%n_basis < 1) &
        .or. (any(basis%n_basis_l(:) < 0)) &
        .or. (any(basis%alpha_l(:) < 0.0d0)) &
        .or. (n_r < 1)) then
      i_err = 1
      return
    end if

    ! basis: set <n_r>, allocate <r_grid>, <radial>
    basis%n_r = n_r
    allocate(basis%r_grid(n_r))
    allocate(basis%radial(n_r, basis%n_basis))

    ! basis: set <r_grid>
    basis%r_grid(:) = r_grid(:)

    ! basis function index offset, incremented by n_basis_l(ll) after each loop
    ii = 0

    ! loop over <l>, basis: set <radial>
    do ll = 0, basis%l_max
      ! in-line <n_basis_l>, <alpha> for current <l>
      n_b_l = basis%n_basis_l(ll)
      alpha = basis%alpha_l(ll)

      ! if there are no basis states for this value of <l>, then cycle
      if (n_b_l == 0) then
        cycle
      end if

      ! allocate normalisation constant array
      allocate(norm(n_b_l))

      ! recurrence relation for basis normalisation constants
      if (n_b_l >= 1) then
        norm(1) = sqrt(alpha / dble((ll+1) * gamma(dble((2*ll)+2))))
      end if

      if (n_b_l >= 2) then
        do kk = 2, n_b_l
          norm(kk) = norm(kk-1) * sqrt(dble((kk-1) * (kk-1+ll)) &
              / dble((kk+ll) * (kk+(2*ll))))
        end do
      end if

      ! in-line <alpha_grid> for current <l>
      alpha_grid(:) = alpha * r_grid(:)

      ! recurrence relation for basis functions
      if (n_b_l >= 1) then
        basis%radial(:, ii+1) = ((2.0d0 * alpha_grid(:)) ** (ll+1)) &
            * exp(-alpha_grid(:))
      end if

      if (n_b_l >= 2) then
        basis%radial(:, ii+2) = 2.0d0 * (dble(ll+1) - alpha_grid(:)) &
            * basis%radial(:, ii+1)
      end if

      if (n_b_l >= 3) then
        do kk = 3, n_b_l
          basis%radial(:, ii+kk) = &
              ((2.0d0 * (dble(kk-1+ll) - alpha_grid(:)) &
              * basis%radial(:, ii+kk-1)) &
              - dble(kk+(2*ll)-1) * basis%radial(:, ii+kk-2)) &
              / dble(kk-1)
        end do
      end if

      ! scale basis radial functions by normalisation constants
      if (n_b_l >= 1) then
        do kk = 1, n_b_l
          basis%radial(:, ii+kk) = basis%radial(:, ii+kk) * norm(kk)
        end do
      end if

      ! increment basis index offset
      ii = ii + n_b_l

      ! deallocate norm
      deallocate(norm)
    end do

  end subroutine setup_radial

end module laguerre
