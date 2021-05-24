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
  double precision , allocatable :: B(:, :)
  double precision , allocatable :: K(:, :)
  double precision , allocatable :: V(:, :)

  ! local variables
  integer :: i_err
  integer :: ii

!> program execution
  ! molecule parameters
  m = 0
  parity = 1
  l_max = 3

  allocate(n_basis_l(0:l_max))
  allocate(alpha_l(0:l_max))

  n_basis_l(:) = 5
  alpha_l(:) = 1.0d0

  call setup_basis(basis, m, parity, l_max, n_basis_l, alpha_l, i_err)

  ! radial grid parameters
  d_r = 1.0d-1
  r_max = 1.0d2
  n_r = ceiling(r_max / d_r) + 1
  allocate(r_grid(n_r))

  do ii = 1, n_r
    r_grid(ii) = d_r * dble(ii - 1)
  end do

  call setup_radial(basis, n_r, r_grid, i_err)

  ! matrices
  allocate(B(basis%n_basis, basis%n_basis))
  allocate(K(basis%n_basis, basis%n_basis))
  allocate(V(basis%n_basis, basis%n_basis))

  call overlap(basis, B, i_err)
  call display_matrix(basis%n_basis, basis%n_basis, B)
  call overlap_numeric(basis, B, i_err)
  call display_matrix(basis%n_basis, basis%n_basis, B)

  call kinetic(basis, K, i_err)
  call display_matrix(basis%n_basis, basis%n_basis, K)

  call potential_e_n(basis, 1, 1.0d-5, 10, V, i_err)
  call display_matrix(basis%n_basis, basis%n_basis, V)

contains

end program homonuclear_diatomic
