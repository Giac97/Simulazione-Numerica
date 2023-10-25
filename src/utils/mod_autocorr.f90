module mod_autocorr
    use iso_fortran_env, only: int32, int64, real64
    implicit none

contains


subroutine calculate_autocorrelation(U, tau_max, chi)
    implicit none
    real(real64), intent(in)    ::  U(:)
    integer(int64), intent(in) ::  tau_max
    real(real64), intent(out)   ::  chi(:)

    integer(int64)     ::  tau, tp 

    real(real64)       ::  sum_a, sum_b, sum_c, sum_d, sum_e
    real(real64)       ::  frac_a, frac_b, num, den
    integer(int64)     ::  t_max

    t_max = size(U)

    do tau = 1, tau_max
        sum_a = 0.
        sum_b = 0.
        sum_c = 0.
        sum_d = 0.
        sum_e = 0.

        do tp = 1, t_max - tau
            sum_a = sum_a + U(tp) * U(tp + tau)
            sum_b = sum_b + U(tp)
            sum_c = sum_c + U(tp + tau)
        enddo

        do tp = 1, t_max
            sum_d = sum_d + U(tp) * U(tp)
            sum_e = sum_e + U(tp)
        enddo

        frac_a = 1. / real(t_max - tau, real64)
        frac_b = 1. / real(t_max, real64)
        
        num = frac_a * sum_a - frac_a * sum_b * frac_a * sum_c
        den = frac_b * sum_d - frac_b * sum_e * frac_b * sum_e

        chi(tau) = num / den 

    enddo 

    
end subroutine calculate_autocorrelation

subroutine CountLines(filename, num_lines, status)
    use iso_fortran_env, only: real64, int64
  implicit none
    character(len=*), intent(in) :: filename
  integer(int64), intent(out) :: num_lines
  integer(int64), intent(out) :: status
  integer(int64) :: unit, i
  character(1000) :: line

  ! Initialize variables
  num_lines = 0
  status = 0

  ! Attempt to open the file
  open(unit=unit, file=filename, status="OLD", action="READ", iostat=status)

  if (status /= 0) then
    return
  end if

  ! Count the number of lines in the file
  do
    read(unit, '(A)', iostat=status) line
    if (status /= 0) exit
    num_lines = num_lines + 1
  end do

  ! Close the file
  close(unit, status='keep')

end subroutine CountLines

subroutine readFile(filename, num_lines, t, U)
    use iso_fortran_env, only: int64, real64
    implicit none
    character(len=*), intent(in)    ::  filename
    integer(int64), intent(in)      ::  num_lines
    integer(int64), intent(out)     ::  t(num_lines)
    real(real64), intent(out)       ::  U(num_lines)
    
    real(real64)       ::  U_p, err_i
    integer(int64)     ::  i
    
    open(20, file=filename, status="OLD", action="READ")
    do i = 1, num_lines
        read(20, *) t(i), U(i), U_p, err_i
    enddo
    close(20)
end subroutine readFile

end module mod_autocorr
