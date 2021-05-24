!>
program potential_curves

  ! debug compilation
  ! - <STDERR>: file unit for stderr output;
  ! - <DEBUG_POTENTIAL_CURVES>: verbosity of debug statements;
  !   - 0: none;
  !   - 1: control flow (entering and exiting subroutines/functions);
  !   - 2: allocation, scalar assignment, short array assignment, etc;
  !   - 3: longer array assignment;
  !   - 4: internal variable assignment, allocation, inspecting loops, etc;
  ! - <PREFIX>: prefix every debug statement with this string.
  ! - <ERR>: prefix every error debug statement with this string.
#define STDERR 0
#define DEBUG_POTENTIAL_CURVES 4
#define PREFIX "[debug] "
#define ERR "[error] "
#define DISPLAY_BASIS 0
#define DISPLAY_VECTOR 0
#define DISPLAY_MATRIX 0
#define TOL 1.0D-10

  use io
  use laguerre

  implicit none

  ! global flags

  ! radial grid variables
  integer :: n_r
  double precision :: d_r, r_max
  double precision , allocatable :: r_grid(:)

  ! basis parameters
  integer :: m
  integer :: parity
  integer :: l_max
  integer , allocatable :: n_basis_l(:)
  double precision , allocatable :: alpha_l(:)

  ! one-electron homonuclear-diatomtic-molecule parameters
  integer :: nuclei_charge
  integer :: lambda_max

  ! axial distance grid parameters
  integer :: n_rz
  double precision :: d_rz, rz_max
  double precision , allocatable :: rz_grid(:)

  ! basis variables
  type(t_basis) :: basis

  ! matrix variables
  double precision , allocatable :: B(:, :), K(:, :), V(:, :), H(:, :)
  double precision , allocatable :: eigen_values(:), eigen_vectors(:, :)

  ! local variables
  integer :: i_err
  integer :: ii

!> program execution

#if (DEBUG_POTENTIAL_CURVES >= 1)
  write (STDERR, *) PREFIX, "program potential_curves"
#endif

  ! read parameters from command line arguments
  call read_input(m, parity, l_max, n_basis_l, alpha_l, nuclei_charge, &
      lambda_max, d_r, r_max, d_rz, rz_max)

  ! check if <r_grid>, <rz_grid> parameters are valid
  if (d_r < TOL) then
    i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<d_r> < TOL"
#endif

  end if

  if (r_max < TOL) then
    i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<r_max> < TOL"
#endif

  end if

  if (d_rz < TOL) then
    i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<d_rz> < TOL"
#endif

  end if

  if (rz_max < TOL) then
    i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<rz_max> < TOL"
#endif

  end if

  ! handle invalid parameters
  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "grid parameters are invalid"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

#if (DEBUG_LAGUERRE >= 2)
    write (STDERR, *) PREFIX, "radial parameters are valid"
#endif

  ! set <n_r>, <n_rz>
  n_r = ceiling(r_max / d_r) + 1
  n_rz = ceiling(rz_max / d_rz) + 1

  ! allocate and set <r_grid>, <rz_grid>
  allocate(r_grid(n_r))
  allocate(rz_grid(n_rz))

  do ii = 1, n_r
    r_grid(ii) = d_r * (ii - 1)
  end do

  do ii = 1, n_rz
    rz_grid(ii) = d_rz * (ii - 1)
  end do

  ! setup <basis>
  call setup_basis(basis, m, parity, l_max, n_basis_l, alpha_l, i_err)

  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "setup_basis() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

  ! setup <basis> radial variables
  call setup_radial(basis, n_r, r_grid, i_err)

  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "setup_radial() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

  ! display basis
#if (DISPLAY_BASIS)
  write (STDERR, *) "<basis%radial>"
  call display_basis(basis%n_r, basis%r_grid, basis%n_basis, basis%radial)
#endif

  ! allocate <B>, <K>, <V>, <H>
  allocate(B(basis%n_basis, basis%n_basis))
  allocate(K(basis%n_basis, basis%n_basis))
  allocate(V(basis%n_basis, basis%n_basis))
  allocate(H(basis%n_basis, basis%n_basis))

  ! calculate overlap matrix
  call overlap(basis, B, i_err)

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<B>"
  call display_matrix(basis%n_basis, basis%n_basis, B)
#endif

  ! calculate kinetic matrix
  call kinetic(basis, K, i_err)

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<K>"
  call display_matrix(basis%n_basis, basis%n_basis, K)
#endif

  ! loop over <rz_grid>
  do ii = 1, n_rz
    ! calculate potential matrix
    call potential_e_n(basis, nuclei_charge, lambda_max, rz_grid(ii), V, i_err)

#if (DISPLAY_MATRIX)
    write (STDERR, *) "<V>"
    call display_matrix(basis%n_basis, basis%n_basis, V)
#endif

    ! calculate electronic hamiltonian
    H(:, :) = K(:, :) + V(:, :) + (1.0d0 / rz_grid(ii))

#if (DISPLAY_MATRIX)
    write (STDERR, *) "<H>"
    call display_matrix(basis%n_basis, basis%n_basis, H)
#endif

    ! solve electronic eigenvalue equation
    call diagonalise(basis, B, H, eigen_values, eigen_vectors, i_err)

#if (DISPLAY_VECTOR)
    write (STDERR, *) "<eigen_values>"
    call display_vector(basis%n_basis, eigen_values)
#endif

#if (DISPLAY_MATRIX)
    write (STDERR, *) "<eigen_vectors>"
    call display_matrix(basis%n_basis, basis%n_basis, eigen_vectors)
#endif

  end do

contains

  ! diagonalise
  !
  ! Note that since the call to rsg modifies the matrices it is given, we send
  ! it copies of B, H.
  subroutine diagonalise (basis, B, H, eigen_values, eigen_vectors, i_err)
    type(t_basis) , intent(in) :: basis
    double precision , intent(in) :: B(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: H(basis%n_basis, basis%n_basis)
    double precision , intent(out) :: eigen_values(basis%n_basis)
    double precision , intent(out) :: &
        eigen_vectors(basis%n_basis, basis%n_basis)
    integer , intent(out) :: i_err
    double precision :: B_copy(basis%n_basis, basis%n_basis)
    double precision :: H_copy(basis%n_basis, basis%n_basis)

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "subroutine diagonalise()"
#endif

    ! create copies of B, H matrices to use in call to rsg subroutine
    B_copy(:, :) = B(:, :)
    H_copy(:, :) = H(:, :)

    ! solve eigenvalue matrix equation
    eigen_values(:) = 0.0d0
    eigen_vectors(:, :) = 0.0d0

    call rsg(basis%n_basis, basis%n_basis, H_copy, B_copy, eigen_values, 1, &
        eigen_vectors, i_err)

    if (i_err /= 0) then

#if (DEBUG_LAGUERRE >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

      write (*, *) "rsg() failed"
      write (*, *) "exiting subroutine diagonalise()"

      return
    end if

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end subroutine diagonalise()"
#endif

  end subroutine diagonalise

end program potential_curves
