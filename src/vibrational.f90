!>
program vibrational

  ! debug compilation
  ! - <STDERR>: file unit for stderr output;
  ! - <DEBUG_VIBRATIONAL>: verbosity of debug statements;
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
#define DEBUG_VIBRATIONAL 0
#define PREFIX "[debug] "
#define ERR "[error] "
#define TOL 1.0D-10
#define DISPLAY_BASIS 0
#define DISPLAY_VECTOR 0
#define DISPLAY_MATRIX 0
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
  integer , parameter :: m = 0
  integer , parameter :: parity = 1
  integer , parameter :: l_max = 0
  integer , allocatable :: n_basis_l(:)
  double precision , allocatable :: alpha_l(:)

  integer :: n_basis_l_const
  double precision :: alpha_l_const

  ! one-electron homonuclear-diatomtic-molecule parameters
  double precision :: reduced_mass
  double precision , allocatable :: v_grid(:)

  ! basis variables
  type(t_basis) :: basis

  ! matrix variables
  double precision , allocatable :: B(:, :), K(:, :), V(:, :), H(:, :)
  double precision , allocatable :: eigen_values(:), eigen_vectors(:, :)

  ! eigen_radial radial functions
  double precision , allocatable :: eigen_radial(:, :)

  ! local variables
  character(len=2000) :: parameter_dir
  integer :: i_err
  integer :: ii

!> program execution

#if (DEBUG_VIBRATIONAL >= 1)
  write (STDERR, *) PREFIX, "program vibrational"
#endif

  ! read parameters from command line arguments
  call read_input(n_basis_l_const, alpha_l_const, reduced_mass, &
      d_r, r_max, i_err)

  if (i_err /= 0) then

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "read_input() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
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

  ! set <n_r>
  n_r = ceiling(r_max / d_r) + 1

  ! allocate and set <r_grid>
  allocate(r_grid(n_r))

#if (DEBUG_VIBRATIONAL >= 4)
  write (STDERR, *) PREFIX, "<radial index>, <r_grid>"
#endif

  do ii = 1, n_r
    r_grid(ii) = d_r * (ii - 1)

#if (DEBUG_VIBRATIONAL >= 4)
    write (STDERR, *) PREFIX, ii, r_grid(ii)
#endif

  end do

  ! setup <basis>
  call setup_basis(basis, m, parity, l_max, n_basis_l, alpha_l, i_err)

  if (i_err /= 0) then

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "setup_basis() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
#endif

    call exit(i_err)
  end if

  ! setup <basis> radial variables
  call setup_radial(basis, n_r, r_grid, i_err)

  if (i_err /= 0) then

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "setup_radial() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
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

  ! allocate <eigen_values>, <eigen_vectors>, <eigen_radial>
  allocate(eigen_values(basis%n_basis))
  allocate(eigen_vectors(basis%n_basis, basis%n_basis))
  allocate(eigen_radial(basis%n_r, basis%n_basis))

  ! calculate overlap matrix
  call overlap(basis, B, i_err)

  if (i_err /= 0) then
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "overlap() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
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
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "kinetic() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
#endif
    call exit(i_err)
  end if

  ! scale kinetic matrix for reduced mass
  K(:, :) = K(:, :) / reduced_mass

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<K>"
  call display_matrix(basis%n_basis, basis%n_basis, K)
#endif

  ! construct path of output directory
  parameter_dir = parameter_directory(basis, n_basis_l_const, alpha_l_const, &
      reduced_mass)

  ! allocate <v_grid>, and interpolated from potential energy curvy onto
  ! <r_grid>
  allocate(v_grid(n_r))

  call interpolate_potential(n_r, r_grid, "analytic_data/PEC.1ssg", v_grid, &
      i_err)

  if (i_err /= 0) then
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "interpolate_potential() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
#endif
    call exit(i_err)
  end if

  ! calculate potential matrix
  call potential_spherical(basis, v_grid, V, i_err)

  if (i_err /= 0) then
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "potential_spherical() failed"
    write (STDERR, *) PREFIX, ERR, "exiting program vibrational"
#endif
    call exit(i_err)
  end if

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<V>"
  call display_matrix(basis%n_basis, basis%n_basis, V)
#endif

  ! calculate vibrational hamiltonian
  H(:, :) = K(:, :) + V(:, :)

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<H>"
  call display_matrix(basis%n_basis, basis%n_basis, H)
#endif

  ! solve electronic eigenvalue equation
  call diagonalise(basis, B, H, eigen_values, eigen_vectors, eigen_radial, i_err)

  if (i_err /= 0) then
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, ERR, "diagonalise() failed"
#endif
    call exit(i_err)
  end if

#if (DISPLAY_VECTOR)
  write (STDERR, *) "<eigen_values>"
  call display_vector(basis%n_basis, eigen_values)
#endif

#if (DISPLAY_MATRIX)
  write (STDERR, *) "<eigen_vectors>"
  call display_matrix(basis%n_basis, basis%n_basis, eigen_vectors)
#endif

#if (DISPLAY_BASIS)
  write (STDERR, *) "<eigen_radial>"
  call display_matrix(basis%n_r, basis%r_grid, basis%n_basis, eigen_radial)
#endif

  ! write output
  call write_output(basis, B, K, V, H, eigen_values, eigen_vectors, &
      eigen_radial, parameter_dir)

#if (DEBUG_VIBRATIONAL >= 1)
  write (STDERR, *) PREFIX, "end program vibrational"
#endif

contains

  ! diagonalise
  !
  ! Note that since the call to rsg modifies the matrices it is given, we send
  ! it copies of B, H.
  subroutine diagonalise (basis, B, H, eigen_values, eigen_vectors, &
      eigen_radial, i_err)
    type(t_basis) , intent(in) :: basis
    double precision , intent(in) :: B(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: H(basis%n_basis, basis%n_basis)
    double precision , intent(out) :: eigen_values(basis%n_basis)
    double precision , intent(out) :: &
        eigen_vectors(basis%n_basis, basis%n_basis)
    double precision , intent(out) :: eigen_radial(basis%n_r, basis%n_basis)
    integer , intent(out) :: i_err
    double precision :: B_copy(basis%n_basis, basis%n_basis)
    double precision :: H_copy(basis%n_basis, basis%n_basis)
    double precision :: temp_sum
    integer :: ii, jj, kk

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

    ! calculate eigen-basis
    eigen_radial(:, :) = 0.0d0

    do jj = 1, basis%n_basis
      do ii = 1, basis%n_r
        temp_sum = 0.0d0

        do kk = 1, basis%n_basis
          temp_sum = temp_sum + (eigen_vectors(kk, jj) * basis%radial(ii, kk))
        end do

        eigen_radial(ii, jj) = temp_sum
      end do
    end do

#if (DEBUG_POTENTIAL_CURVES >= 1)
    write (STDERR, *) PREFIX, "end subroutine diagonalise()"
#endif

  end subroutine diagonalise

  ! read_input
  subroutine read_input (n_basis_l_const, alpha_l_const, reduced_mass, &
      d_r, r_max, i_err)
    integer , intent(out) :: n_basis_l_const
    double precision , intent(out) :: alpha_l_const, reduced_mass
    double precision , intent(out) :: d_r, r_max
    integer , intent(out) :: i_err
    integer :: num_args
    character(len=50) :: arg

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "subroutine read_input()"
#endif

    ! check if command-line arguments are valid
    i_err = 0

    num_args = command_argument_count()

    if (num_args < 5) then
      i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "<num_args> < 5"
      write (STDERR, *) PREFIX, ERR, "required arguments are: ",  &
          "<n_basis_l_const> <alpha_l_const> <reduced_mass> ", &
          "<d_r> <r_max>"
#endif
    else

      ! read <n_basis_l_const>
      call get_command_argument(1, arg)
      read (arg, *) n_basis_l_const

      if (n_basis_l_const < 0) then
        i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
        write (STDERR, *) PREFIX, ERR, "<n_basis_l_const> < 0)"
#endif
      end if

      ! read <alpha_l_const>
      call get_command_argument(2, arg)
      read (arg, *) alpha_l_const

      if (alpha_l_const < TOL) then
        i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
        write (STDERR, *) PREFIX, ERR, "<alpha_l_const> < TOL)"
#endif
      end if

      ! read <reduced_mass>
      call get_command_argument(3, arg)
      read (arg, *) reduced_mass

      if (reduced_mass < TOL) then
        i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
        write (STDERR, *) PREFIX, ERR, "<reduced_mass> < TOL)"
#endif
      end if

      ! read <d_r>
      call get_command_argument(4, arg)
      read (arg, *) d_r

      if (d_r < TOL) then
        i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
        write (STDERR, *) PREFIX, ERR, "<d_r> < TOL"
#endif
      end if

      ! read <r_max>
      call get_command_argument(5, arg)
      read (arg, *) r_max

      if (r_max < TOL) then
        i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
        write (STDERR, *) PREFIX, ERR, "<r_max> < TOL"
#endif
      end if
    end if

    ! handle invalid command-line arguments
    if (i_err /= 0) then
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
      write (STDERR, *) PREFIX, ERR, "command-line arguments are invalid"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine read_input()"
#endif
      return
    end if

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "command-line arguments are valid"
#endif

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "<n_basis_l_const> = ", n_basis_l_const
    write (STDERR, *) PREFIX, "<alpha_l_const> = ", alpha_l_const
    write (STDERR, *) PREFIX, "<reduced_mass> = ", reduced_mass
    write (STDERR, *) PREFIX, "<d_r> = ", d_r
    write (STDERR, *) PREFIX, "<r_max> = ", r_max
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "end subroutine read_input()"
#endif

  end subroutine read_input

  ! read_grid
  !
  ! For given <n_r>, <r_grid>, <filename>, read the potential from
  ! "<filename>" and interpolate the potential onto <r_grid>.
  subroutine interpolate_potential (n_r, r_grid, filename, v_grid, i_err)
    integer , intent(in) :: n_r
    double precision , intent(in) :: r_grid(n_r)
    character(len=*) , intent(in) :: filename
    double precision , intent(out) :: v_grid(n_r)
    integer , intent(out) :: i_err
    double precision , allocatable :: r_grid_raw(:), v_grid_raw(:)
    double precision :: r, v
    integer :: n_r_raw
    integer :: fileunit
    integer :: ii
    integer :: io_status

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "subroutine interpolate_potential()"
#endif

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="read")

    n_r_raw = 0
    io_status = 0
    do while (io_status == 0)
      read (fileunit, *, iostat=io_status) r, v
      if (io_status == 0) then
        n_r_raw = n_r_raw + 1
      end if
    end do

#if (DEBUG_VIBRATIONAL >= 3)
    write (STDERR, *) PREFIX, "<n_r_raw> = ", n_r_raw
#endif

    ! handle invalid file read
    if (n_r_raw == 0) then
      i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "<n_r_raw> == 0"
#endif
    end if

    ! handle invalid file read
    if (i_err /= 0) then
      i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
      write (STDERR, *) PREFIX, ERR, "file could not be processed"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine interpolate_potential()"
#endif
      return
    end if

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "file has valid lines"
#endif

    ! re-open file at beginning
    close (fileunit)
    open (unit=fileunit, file=trim(adjustl(filename)), action="read")

    ! allocate <r_grid_raw>, <v_grid_raw>
    allocate(r_grid_raw(n_r_raw))
    allocate(v_grid_raw(n_r_raw))

    ! read potential file
#if (DEBUG_VIBRATIONAL >= 4)
    write (STDERR, *) PREFIX, "<ii>, <r_grid_raw(ii)>, <v_grid_raw(ii)>"
#endif
    ii = 1
    io_status = 0
    do while ((ii <= n_r_raw) .and. (io_status == 0))
      read (fileunit, *, iostat=io_status) r_grid_raw(ii), v_grid_raw(ii)
      if (io_status == 0) then
#if (DEBUG_VIBRATIONAL >= 4)
        write (STDERR, *) PREFIX, ii, r_grid_raw(ii), v_grid_raw(ii)
#endif
        ii = ii + 1
      end if
    end do

    ! handle invalid file read
    if (ii <= n_r_raw) then
      i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "less than <n_r_raw> lines read successfully"
#endif
    end if

    ! handle invalid file read
    if (i_err /= 0) then
      i_err = 1
#if (DEBUG_VIBRATIONAL >= 2)
      write (STDERR, *) PREFIX, ERR, "<i_err> = ", i_err
#endif
#if (DEBUG_VIBRATIONAL >= 1)
      write (STDERR, *) PREFIX, ERR, "file could not be processed"
      write (STDERR, *) PREFIX, ERR, "exiting subroutine interpolate_potential()"
#endif
      return
    end if

    ! close file
    close (fileunit)

    ! interpolate <v_grid> on <r_grid> from <v_grid_raw> on <r_grid_raw>
    call INTRPL(n_r_raw, r_grid_raw, v_grid_raw, n_r, r_grid, v_grid)

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "end subroutine interpolate_potential()"
#endif

  end subroutine interpolate_potential

  ! parameter_directory
  !
  ! For given <basis>, <n_basis_l_const>, <alpha_l_const>, <reduced_mass>,
  ! construct the name of the directory that output for this calculation will be
  ! written to.
  function parameter_directory (basis, n_basis_l_const, alpha_l_const, &
      reduced_mass) result (dir)
    type(t_basis) , intent(in) :: basis
    integer , intent(in) :: n_basis_l_const
    double precision , intent(in) :: alpha_l_const, reduced_mass
    character(len=2000) :: dir
    character(len=100) :: str_n_basis_l_const,  str_alpha_l_const, &
        str_reduced_mass, str_d_r, str_r_max

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "function parameter_directory()"
#endif

    ! construct output directory for given parameters
    write (str_n_basis_l_const, *) n_basis_l_const
    write (str_alpha_l_const, DP_FORMAT) alpha_l_const
    write (str_reduced_mass, DP_FORMAT) reduced_mass
    write (str_d_r, DP_FORMAT) basis%r_grid(2) - basis%r_grid(1)
    write (str_r_max, DP_FORMAT) basis%r_grid(basis%n_r)

    write (dir, *) &
        "output/vib/", &
        "n_basis_l_const-", trim(adjustl(str_n_basis_l_const)), ".", &
        "alpha_l_const-", trim(adjustl(str_alpha_l_const)), ".", &
        "reduced_mass-", trim(adjustl(str_reduced_mass)), ".", &
        "d_r-", trim(adjustl(str_d_r)), ".", &
        "r_max-", trim(adjustl(str_r_max)), "/"

#if (DEBUG_VIBRATIONAL >= 3)
    write (STDERR, *) PREFIX, "<dir> = ", trim(adjustl(dir))
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "end function parameter_directory()"
#endif

  end function parameter_directory

  ! write_output
  subroutine write_output (basis, B, K, V, H, eigen_values, eigen_vectors, &
      eigen_radial, dir)
    type(t_basis) , intent(in) :: basis
    double precision , intent(in) :: B(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: K(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: V(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: H(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: eigen_values(basis%n_basis)
    double precision , intent(in) :: eigen_vectors(basis%n_basis, basis%n_basis)
    double precision , intent(in) :: eigen_radial(basis%n_r, basis%n_basis)
    character(len=*) , intent(in) :: dir

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "subroutine write_output()"
#endif

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "<dir> = ", trim(adjustl(dir))
#endif

    call execute_command_line("mkdir -p "//trim(adjustl(dir)))

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "made <dir>"
#endif

    ! write <B>, <K> to file
    call write_matrix(basis%n_basis, basis%n_basis, B, trim(adjustl(dir))//"B.txt")
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <B> in <dir>"
#endif

    call write_matrix(basis%n_basis, basis%n_basis, K, trim(adjustl(dir))//"K.txt")
#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <K> in <dir>"
#endif

    ! write <V>, <H> to file
    call write_matrix(basis%n_basis, basis%n_basis, V, trim(adjustl(dir))//"V.txt")

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <V> in <dir>"
#endif

    call write_matrix(basis%n_basis, basis%n_basis, H, trim(adjustl(dir))//"H.txt")

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <H> in <dir>"
#endif

    ! write <eigen_values>, <eigen_vectors> to file
    call write_vector(basis%n_basis, eigen_values, &
        trim(adjustl(dir))//"eigen_values.txt")

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <eigen_values> in <dir>"
#endif

    call write_matrix(basis%n_basis, basis%n_basis, eigen_vectors, &
        trim(adjustl(dir))//"eigen_vectors.txt")

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <eigen_vectors> in <dir>"
#endif

    call write_basis(basis%n_r, basis%r_grid, basis%n_basis, eigen_radial, &
        trim(adjustl(dir))//"eigen_radial.txt")

#if (DEBUG_VIBRATIONAL >= 2)
    write (STDERR, *) PREFIX, "written <eigen_radial> in <dir>"
#endif

#if (DEBUG_VIBRATIONAL >= 1)
    write (STDERR, *) PREFIX, "end subroutine write_output()"
#endif

  end subroutine write_output

end program vibrational
