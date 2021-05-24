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

  ! basis variables
  type(t_basis) :: basis

  ! matrix variables

  ! local variables
  integer :: i_err

  ! program execution
  m = 0
  parity = 1
  l_max = 5

  allocate(n_basis_l(0:l_max))
  allocate(alpha_l(0:l_max))

  n_basis_l(:) = 5
  alpha_l(:) = 1.0d0

  call setup_basis(basis, m, parity, l_max, n_basis_l, alpha_l, i_err)

contains

end program homonuclear_diatomic
