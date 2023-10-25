program autocorrelation
    use iso_fortran_env, only: int64, real64, int32
    use mod_autocorr, only: calculate_autocorrelation, readFile, CountLines
    implicit none
    integer(int32)      ::  i
    integer (int64)     ::   j
    integer(int64)      ::  n_vals, status, tau_max, limit


    real(real64), allocatable   ::   U(:)
    integer(int64), allocatable ::  t(:)
    real(real64), allocatable   ::  t_r(:)
    real(real64), allocatable   ::  chi(:)

    character(256)      ::  arg, fname
    fname = "output_epot.dat"
    tau_max = 200
    limit = 0
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        
        if (trim(arg) == "-help") then
            print*, ""
            print*, "use -taumax or --tm= and an integer value to se the value of tau_max"
            print*, "-fname to enter the name of the file"
        goto 200
        
        else if (trim(arg) == "-taumax" .or. trim(arg) == "--tm=") then
            if (i < command_argument_count()) then
                call get_command_argument(i + 1,arg)
                
                read(arg, *) tau_max
            
            else
                print *, "tau_max requires a value"
                stop
            endif
        else if (trim(arg) == "-fname") then
            if (i < command_argument_count()) then
                call get_command_argument(i + 1, arg)
                read(arg, *) fname
            else
                print *, "File name requires a value"
                stop
            endif
        else if (trim(arg) == "-limit") then
            if (i < command_argument_count()) then
                call get_command_argument(i + 1, arg)
                read(arg, *) limit
            else
                print *, "Limit requires a a value"
                stop
            endif
        endif


    enddo
    call CountLines(fname, n_vals, status) 
    if (limit .ne. 0) then
        n_vals = limit
    endif
    allocate(t(n_vals))
    allocate(U(n_vals))
    allocate(t_r(n_vals))
!    tau_max = 200
    allocate(chi(tau_max))
    call readFile(fname, n_vals, t, U)

    t_r = dble(t)
    
    call calculate_autocorrelation(U, tau_max, chi)

    do i = 1, tau_max
        print *,chi(i) / chi(1)
    enddo
    deallocate(t)
    deallocate(t_r)
    deallocate(U)
    deallocate(chi)
200 continue




end program autocorrelation



