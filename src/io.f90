!>
module io

  implicit none

  ! global flags
  integer , parameter :: dp_default = 4

contains

  ! write_vector
  subroutine write_vector (n, x, filename)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write vector to file
    do ii = 1, n
      write (fileunit, *) x(ii)
    end do

    ! close file
    close (fileunit)

  end subroutine write_vector

  ! write_matrix
  subroutine write_matrix (n_rows, n_cols, A, filename)
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: A(n_rows, n_cols)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write matrix to file
    do ii = 1, n_rows
      write (fileunit, *) A(ii, :)
    end do

    ! close file
    close (fileunit)

  end subroutine write_matrix

  ! write_basis
  subroutine write_basis (n_r, r_grid, n_basis, basis, filename)
    integer , intent(in) :: n_r, n_basis
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: basis(n_r, n_basis)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write r_grid, and basis functions to file
    do ii = 1, n_r
      write (fileunit, *) r_grid(ii), " ", basis(ii, :)
    end do

    ! close file
    close (fileunit)

  end subroutine write_basis

  ! display_vector
  subroutine display_vector (n, x, dp)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    integer , optional , intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = dp_default
    end if

    w = max(ceiling(log10(maxval(abs(x(:))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out vector elements
    do ii = 1, n
      ! if x(ii) will be written as "0.00..0", replace with " .     "
      if (abs(x(ii)) > (10.0**(-d))) then
        write (*, fmt) x(ii)
      else
        write (*, "(a, a, a)") str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
      end if
    end do

  end subroutine display_vector

  ! display_matrix
  subroutine display_matrix (n_rows, n_cols, A, dp)
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: A(n_rows, n_cols)
    integer , optional, intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii, jj

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = dp_default
    end if

    w = max(ceiling(log10(maxval(abs(A(:, :))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out matrix elements
    do ii = 1, n_rows
      do jj = 1, n_cols
        ! if A(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(A(ii, jj)) > (10.0**(-d))) then
          write (*, fmt, advance="no") A(ii, jj)
        else
          write (*, "(a, a, a)", advance="no") &
              str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
        end if
      end do
      write (*, *)
    end do

  end subroutine display_matrix

  ! display_basis
  subroutine display_basis (n_r, r_grid, n_basis, basis, dp)
    integer , intent(in) :: n_r, n_basis
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: basis(n_r, n_basis)
    integer , optional , intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii, jj

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = dp_default
    end if

    w = max(ceiling(log10(maxval(abs(basis(:, :))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out radial grid and radial basis values
    do ii = 1, n_r
      write (*, fmt, advance="no") r_grid(ii)
      do jj = 1, n_basis
        ! if basis(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(basis(ii, jj)) > (10.0**(-d))) then
          write (*, fmt, advance="no") basis(ii, jj)
        else
          write (*, "(a, a, a)", advance="no") &
              str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
        end if
      end do
      write (*, *)
    end do

  end subroutine display_basis

end module io
