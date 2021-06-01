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
  ! - <PREFIX>: prefix every debug statement with this string;
  ! - <ERR>: prefix every error debug statement with this string;
  ! - <TOL>: double precision tolerance value;
  ! - <DISPLAY_BASIS>: flag if radial basis functions should be displayed;
  ! - <DISPLAY_VECTOR>: flag if vectors should be displayed;
  ! - <DISPLAY_MATRIX>: flag if matrices should be displayed.
#define STDERR 0
#define DEBUG_POTENTIAL_CURVES 2
#define PREFIX "[debug] "
#define ERR "[error] "
#define TOL 1.0D-10
#define DISPLAY_BASIS 0
#define DISPLAY_VECTOR 1
#define DISPLAY_MATRIX 1
#define DP_FORMAT "(f10.4)"

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
  character(len=2000) :: parameter_dir, axial_dir
  logical :: output_exists
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
  n_rz = ceiling(rz_max / d_rz)

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
    rz_grid(ii) = d_rz * ii

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

  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "overlap() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<B>"
  call display_matrix(basis%n_basis, basis%n_basis, B)
#endif

  ! calculate kinetic matrix
  call kinetic(basis, K, i_err)

  if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, ERR, "kinetic() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program potential_curves"
#endif

    call exit(i_err)
  end if

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<K>"
  call display_matrix(basis%n_basis, basis%n_basis, K)
#endif

  ! construct path of parameter directory
  parameter_dir = parameter_directory(basis, n_basis_l_const, alpha_l_const, &
      nuclei_charge, lambda_max)

#if (DEBUG_POTENTIAL_CURVES >= 2)
  write (STDERR, *) PREFIX, "<parameter_dir> = ", trim(adjustl(parameter_dir))
#endif

  ! loop over <rz_grid>
  do ii = 1, n_rz

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "<rz> = ", rz_grid(ii)
#endif

    ! check if data already exists for this calculation
    axial_dir = axial_directory(parameter_dir, rz_grid(ii))

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "<axial_dir> = ", trim(adjustl(axial_dir))
#endif

    inquire(FILE=trim(adjustl(axial_dir))//"complete.txt", EXIST=output_exists)

    if (output_exists) then
#if (DEBUG_POTENTIAL_CURVES >= 1)
      write (STDERR, *) PREFIX, &
          "<axial_dir> already exists for these parameters at this <rz>"
      write (STDERR, *) PREFIX, "cycling"
#endif

      cycle
    end if

    ! calculate potential matrix
    call potential_e_n(basis, nuclei_charge, lambda_max, rz_grid(ii), V, i_err)

    if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
      write (STDERR, *) PREFIX, ERR, "potential_e_n() failed"
      write (STDERR, *) PREFIX, ERR, "cycling"
#endif

      cycle
    end if

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

    if (i_err /= 0) then

#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
      write (STDERR, *) PREFIX, ERR, "diagonalise() failed"
      write (STDERR, *) PREFIX, ERR, "cycling"
#endif

      cycle
    end if

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

    call write_output(basis%n_basis, B, K, V, H, eigen_values, eigen_vectors, &
        parameter_dir, axial_dir)

  end do

#if (DEBUG_POTENTIAL_CURVES >= 1)
  write (STDERR, *) PREFIX, "end program potential_curves"
#endif

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

  ! write_output
  !
  ! Note: we assume by construction that
  ! > <axial_dir> = <parameter_dir>"rz-<rz>/"
  subroutine write_output (n_basis, B, K, V, H, eigen_values, eigen_vectors, &
      parameter_dir, axial_dir)
    integer , intent(in) :: n_basis
    double precision , intent(in) :: B(n_basis, n_basis)
    double precision , intent(in) :: K(n_basis, n_basis)
    double precision , intent(in) :: V(n_basis, n_basis)
    double precision , intent(in) :: H(n_basis, n_basis)
    double precision , intent(in) :: eigen_values(n_basis)
    double precision , intent(in) :: eigen_vectors(n_basis, n_basis)
    character(len=*) , intent(in) :: parameter_dir, axial_dir
    logical :: matrix_exists

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "subroutine write_output()"
#endif

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "<parameter_dir> = ", trim(adjustl(parameter_dir))
    write (STDERR, *) PREFIX, "<axial_dir> = ", trim(adjustl(axial_dir))
#endif

    call execute_command_line("mkdir -p "//trim(adjustl(parameter_dir)))
    call execute_command_line("mkdir -p "//trim(adjustl(axial_dir)))

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "made <parameter_dir>"
    write (STDERR, *) PREFIX, "made <axial_dir>"
#endif

    ! write <B>, <K> to file
    inquire(FILE=trim(adjustl(parameter_dir))//"B.txt", EXIST=matrix_exists)
    if (matrix_exists) then
#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, "<B> already written in <parameter_dir>"
#endif
    else
      call write_matrix(n_basis, n_basis, B, &
          trim(adjustl(parameter_dir))//"B.txt")
#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, "written <B> in <parameter_dir>"
#endif
    end if

    inquire(FILE=trim(adjustl(parameter_dir))//"K.txt", EXIST=matrix_exists)
    if (matrix_exists) then
#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, "<K> already written in <parameter_dir>"
#endif
    else
      call write_matrix(n_basis, n_basis, K, &
          trim(adjustl(parameter_dir))//"K.txt")
#if (DEBUG_POTENTIAL_CURVES >= 2)
      write (STDERR, *) PREFIX, "written <K> in <parameter_dir>"
#endif
    end if

    ! write <V>, <H> to file
    call write_matrix(n_basis, n_basis, V, trim(adjustl(axial_dir))//"V.txt")

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "written <V> in <axial_dir>"
#endif

    call write_matrix(n_basis, n_basis, H, trim(adjustl(axial_dir))//"H.txt")

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "written <H> in <axial_dir>"
#endif

    ! write <eigen_values>, <eigen_vectors> to file
    call write_vector(n_basis, eigen_values, &
        trim(adjustl(axial_dir))//"eigen_values.txt")

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "written <eigen_values> in <axial_dir>"
#endif

    call write_matrix(n_basis, n_basis, eigen_vectors, &
        trim(adjustl(axial_dir))//"eigen_vectors.txt")

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "written <eigen_vectors> in <axial_dir>"
#endif

    ! touch complete file to register that the calculations for this set of
    ! parameters has been completely written to file
    call execute_command_line( &
        "touch "//trim(adjustl(axial_dir))//"complete.txt")

#if (DEBUG_POTENTIAL_CURVES >= 2)
    write (STDERR, *) PREFIX, "touched <complete.txt> in <axial_dir>"
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end subroutine write_output()"
#endif

  end subroutine write_output

  ! parameter_directory
  !
  ! For given <basis>, <n_basis_l_const>, <alpha_l_const>, <nuclei_charge>,
  ! <lambda_max>, construct the name of the directory that output for this
  ! calculation will be written to.
  function parameter_directory (basis, n_basis_l_const, alpha_l_const, &
      nuclei_charge, lambda_max) result (dir)
    type(t_basis) , intent(in) :: basis
    integer , intent(in) :: n_basis_l_const, nuclei_charge, lambda_max
    double precision , intent(in) :: alpha_l_const
    character(len=2000) :: dir
    character(len=100) :: str_m, str_parity, str_l_max, str_n_basis_l_const, &
        str_alpha_l_const, str_nuclei_charge, str_lambda_max, &
        str_d_r, str_r_max

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "function parameter_directory()"
#endif

    ! construct output directory for given parameters
    write (str_m, *) basis%m
    write (str_parity, *) basis%parity
    write (str_l_max, *) basis%l_max
    write (str_n_basis_l_const, *) n_basis_l_const
    write (str_alpha_l_const, DP_FORMAT) alpha_l_const
    write (str_nuclei_charge, *) nuclei_charge
    write (str_lambda_max, *) lambda_max
    write (str_d_r, DP_FORMAT) basis%r_grid(2) - basis%r_grid(1)
    write (str_r_max, DP_FORMAT) basis%r_grid(basis%n_r)

    write (dir, *) &
        "output/", &
        "m-", trim(adjustl(str_m)), ".", &
        "parity-", trim(adjustl(str_parity)), ".", &
        "l_max-", trim(adjustl(str_l_max)), ".", &
        "n_basis_l_const-", trim(adjustl(str_n_basis_l_const)), ".", &
        "alpha_l_const-", trim(adjustl(str_alpha_l_const)), ".", &
        "nuclei_charge-", trim(adjustl(str_nuclei_charge)), ".", &
        "lambda_max-", trim(adjustl(str_lambda_max)), ".", &
        "d_r-", trim(adjustl(str_d_r)), ".", &
        "r_max-", trim(adjustl(str_r_max)), "/"

#if (DEBUG_POTENTIAL_CURVES >= 3)
    write (STDERR, *) PREFIX, "<dir> = ", trim(adjustl(dir))
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end function parameter_directory()"
#endif

  end function parameter_directory

  ! axial_directory
  !
  ! For given <dir>, <rz>, construct the path of the sub-directory
  ! that output for the calculation at this axial distance, <rz>, will be
  ! written to.
  !
  ! Also used to check if the calculation has already been performed before.
  function axial_directory (dir, rz) result (sub_dir)
    character(len=*) , intent(in) :: dir
    double precision , intent(in) :: rz
    character(len=2000) :: sub_dir
    character(len=100) :: str_rz

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "function axial_directory()"
#endif

    ! construct output directory for given parameters
    write (str_rz, DP_FORMAT) rz

    write (sub_dir, *) &
        trim(adjustl(dir)), "rz-", trim(adjustl(str_rz)), "/"

#if (DEBUG_POTENTIAL_CURVES >= 3)
    write (STDERR, *) PREFIX, "<dir> = ", trim(adjustl(sub_dir))
#endif

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end function axial_directory()"
#endif

  end function axial_directory

end program potential_curves
