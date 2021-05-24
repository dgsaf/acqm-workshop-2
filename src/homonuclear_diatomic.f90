!>
program homonuclear_diatomic

  !
  use io
  use laguerre

  implicit none

  ! global flags

  ! one-electron homonuclear-diatomtic-molecule parameters
  integer :: m
  integer :: parity
  integer :: l_max
  integer , allocatable :: n_basis_l(:)
  double precision , allocatable :: alpha_l(:)

  ! radial grid variables
  integer :: n_r
  double precision :: d_r, r_max
  double precision , allocatable :: r_grid(:)

  ! basis variables
  type(t_basis) :: basis

  ! matrix variables

  ! local variables
  integer :: i_err
  integer :: ii

!> program execution
  ! molecule parameters
  m = 2
  parity = 1
  l_max = 5

  allocate(n_basis_l(0:l_max))
  allocate(alpha_l(0:l_max))

  n_basis_l(:) = 5
  alpha_l(:) = 1.0d0

  call setup_basis(basis, m, parity, l_max, n_basis_l, alpha_l, i_err)

  ! radial grid parameters
  d_r = 1
  r_max = 10.0
  n_r = ceiling(r_max / d_r) + 1
  allocate(r_grid(n_r))

  do ii = 1, n_r
    r_grid(ii) = d_r * dble(ii - 1)
  end do

  call setup_radial(basis, n_r, r_grid, i_err)

contains

end program homonuclear_diatomic
