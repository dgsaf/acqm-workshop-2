!>
module laguerre

  implicit none

  ! t_basis
  !
  ! Laguerre basis, consisting of basis states
  ! > {| phi_{i} >} for i = 1, ..., <n_basis>
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
  ! - <l_max>: maximum angular quantum number, <l>, considered in the basis;
  ! - <n_basis_l>: number of basis functions per <l>;
  ! - <alpha_l>: value of <alpha> per <l>;
  ! - <m>: the magnetic quantum number, a conserved symmetry of the one-electron
  !   homonuclear-diatomic-molecule system;
  ! - <parity>: the parity quantum number, <parity> = (-1)^<l>, a conserved
  !   symmetry of the one-electron homonuclear-diatomic-molecule system.
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
  type t_basis
    integer :: l_max
    integer , allocatable :: n_basis_l(:)
    double precision , allocatable :: alpha_l(:)
    integer :: m
    integer :: parity

    integer :: n_basis
    integer , allocatable :: k_list(:)
    integer , allocatable :: l_list(:)

    integer :: n_r
    double precision , allocatable :: r_grid(:, :)
    double precision , allocatable :: radial(:, :)
  end type t_basis

contains

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
  !
  ! Also returns an error code, <i_err>, where:
  ! - 0 indicates successful execution;
  ! - 1 indicates invalid arguments.
  subroutine setup_radial (basis, n_r, r_grid, i_err)
    type(t_basis) , intent(inout) :: basis
    integer , intent(in) :: n_r
    integer , intent(in) :: r_grid(n_r)
    integer , intent(out) :: i_err
    double precision , allocatable :: norm(:)
    double precision :: alpha_grid(n_r)
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
      ! in-line <n_basis_l> for current <l>
      n_b_l = basis%n_basis_l(ll)

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
      alpha_grid(:) = basis%alpha_l(ll) * r_grid(:)

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
        do kk = 1, n_basis
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
