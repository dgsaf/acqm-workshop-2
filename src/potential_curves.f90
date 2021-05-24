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
#define DISPLAY_VECTOR 1
#define DISPLAY_MATRIX 1
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

  integer :: n_basis_l_const
  double precision :: alpha_l_const

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
  call read_input(m, parity, l_max, n_basis_l_const, alpha_l_const, &
      nuclei_charge, lambda_max, d_r, r_max, d_rz, rz_max, i_err)

  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "read_input() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

  ! set <n_basis_l>, <alpha_l>
  if ((i_err == 0) .and. (l_max >= 0)) then
    allocate(n_basis_l(0:l_max))
    allocate(alpha_l(0:l_max))

    n_basis_l(:) = n_basis_l_const
    alpha_l(:) = alpha_l_const
  end if

  ! set <n_r>, <n_rz>
  n_r = ceiling(r_max / d_r) + 1
  n_rz = ceiling(rz_max / d_rz) + 1

  ! allocate and set <r_grid>, <rz_grid>
  allocate(r_grid(n_r))
  allocate(rz_grid(n_rz))

#if (DEBUG_POTENTIAL_CURVES >= 4)
  write (STDERR, *) PREFIX, "<radial index>, <r_grid>"
#endif

  do ii = 1, n_r
    r_grid(ii) = d_r * (ii - 1)

#if (DEBUG_POTENTIAL_CURVES >= 4)
    write (STDERR, *) PREFIX, ii, r_grid(ii)
#endif

  end do

#if (DEBUG_POTENTIAL_CURVES >= 4)
  write (STDERR, *) PREFIX, "<radial index>, <r_grid>"
#endif

  do ii = 1, n_rz
    rz_grid(ii) = d_rz * (ii - 1)

#if (DEBUG_POTENTIAL_CURVES >= 4)
    write (STDERR, *) PREFIX, ii, rz_grid(ii)
#endif

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

  ! allocate <eigen_values>, <eigen_vectors>
  allocate(eigen_values(basis%n_basis))
  allocate(eigen_vectors(basis%n_basis, basis%n_basis))

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

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "<rz> = ", rz_grid(ii)
#endif

    ! calculate potential matrix
    call potential_e_n(basis, nuclei_charge, lambda_max, rz_grid(ii), V, i_err)

#if (DISPLAY_MATRIX)
    write (STDERR, *) "<V>"
    call display_matrix(basis%n_basis, basis%n_basis, V)
#endif

    ! calculate electronic hamiltonian, not including 1/R term; that is,
    ! > <H_elec> = <H> + (1.0d0 / rz_grid(ii))
    H(:, :) = K(:, :) + V(:, :)

#if (DISPLAY_MATRIX)
    write (STDERR, *) "<H>"
    call display_matrix(basis%n_basis, basis%n_basis, H)
#endif

    ! solve electronic eigenvalue equation
    call diagonalise(basis, B, H, eigen_values, eigen_vectors, i_err)

    ! shift energies by the 1/R term
    eigen_values(:) = eigen_values(:) + (dble(nuclei_charge ** 2) / rz_grid(ii))

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

#if (DEBUG_POTENTIAL_CURVES >= 2)
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

  ! read_input
  subroutine read_input (m, parity, l_max, n_basis_l_const, alpha_l_const, &
      nuclei_charge, lambda_max, d_r, r_max, d_rz, rz_max, i_err)
    integer , intent(out) :: m, parity, l_max, n_basis_l_const, &
        nuclei_charge, lambda_max
    double precision , intent(out) :: alpha_l_const
    double precision , intent(out) :: d_r, r_max, d_rz, rz_max
    integer , intent(out) :: i_err
    integer :: num_args
    character(len=50) :: arg

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "subroutine read_input()"
#endif

    ! check if command-line arguments are valid
    i_err = 0

    num_args = command_argument_count()

    if (num_args < 11) then
      i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, ERR, "<num_args> < 11"

      write (STDERR, *) PREFIX, ERR, "required arguments are: ",  &
          "<m> <parity> <l_max> <n_basis_l_const> <alpha_l_const> ", &
          "<nuclei_charge> <lambda_max> ", &
          "<d_r> <r_max> <d_rz> <rz_max>"
#endif

    else

      ! read <m>
      call get_command_argument(1, arg)
      read (arg, *) m

      ! read <parity>
      call get_command_argument(2, arg)
      read (arg, *) parity

      if (abs(parity) /= 1) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "abs(<parity>) /= 1"
#endif

      end if

      ! read <l_max>
      call get_command_argument(3, arg)
      read (arg, *) l_max

      if (l_max < abs(m)) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<l_max> < abs(<m>)"
#endif

      end if

      if (l_max < 0) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<l_max> < 0"
#endif

      end if

      ! read <n_basis_l_const>
      call get_command_argument(4, arg)
      read (arg, *) n_basis_l_const

      if (n_basis_l_const < 0) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<n_basis_l_const> < 0)"
#endif

      end if

      ! read <alpha_l_const>
      call get_command_argument(5, arg)
      read (arg, *) alpha_l_const

      if (alpha_l_const < TOL) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<alpha_l_const> < TOL)"
#endif

      end if

      ! read <nuclei_charge>
      call get_command_argument(6, arg)
      read (arg, *) nuclei_charge

      if (nuclei_charge < 1) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<nuclei_charge> < 1)"
#endif

      end if

      ! read <lambda_max>
      call get_command_argument(7, arg)
      read (arg, *) lambda_max

      if (lambda_max < 0) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<nuclei_charge> < 0)"
#endif

      end if

      ! read <d_r>
      call get_command_argument(8, arg)
      read (arg, *) d_r

      if (d_r < TOL) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<d_r> < TOL"
#endif

      end if

      ! read <r_max>
      call get_command_argument(9, arg)
      read (arg, *) r_max

      if (r_max < TOL) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<r_max> < TOL"
#endif

      end if

      ! read <d_rz>
      call get_command_argument(10, arg)
      read (arg, *) d_rz

      if (d_rz < TOL) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<d_rz> < TOL"
#endif

      end if

      ! read <rz_max>
      call get_command_argument(11, arg)
      read (arg, *) rz_max

      if (rz_max < TOL) then
        i_err = 1

#if (DEBUG_POTENTIAL_CURVES >= 2)
        write (STDERR, *) PREFIX, ERR, "<rz_max> < TOL"
#endif

      end if

    end if

    ! handle invalid command-line arguments
    if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
      write (STDERR, *) PREFIX, ERR, "command-line arguments are invalid"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine read_input()"
#endif

      return
    end if

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "command-line arguments are valid"
#endif

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "<m> = ", m
    write (STDERR, *) PREFIX, "<parity> = ", parity
    write (STDERR, *) PREFIX, "<l_max> = ", l_max
    write (STDERR, *) PREFIX, "<n_basis_l_const> = ", n_basis_l_const
    write (STDERR, *) PREFIX, "<alpha_l_const> = ", alpha_l_const
    write (STDERR, *) PREFIX, "<nuclei_charge> = ", nuclei_charge
    write (STDERR, *) PREFIX, "<lambda_max> = ", lambda_max
    write (STDERR, *) PREFIX, "<d_r> = ", d_r
    write (STDERR, *) PREFIX, "<r_max> = ", r_max
    write (STDERR, *) PREFIX, "<d_rz> = ", d_rz
    write (STDERR, *) PREFIX, "<rz_max> = ", rz_max
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end subroutine read_input()"
#endif

  end subroutine read_input

end program potential_curves
