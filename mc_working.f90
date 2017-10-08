program monte_carlo
    use types
    use mcrt_functions
    
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), length_travelled(:, :)
    ! real(dp) :: max_tau, min_tau
    real(dp) :: pos(3), nhat(3)
    real(dp) :: leng, r_max = 200.0*AU, r, tau
    integer :: num_packets, len_tau, j, i, u

    write(*,*) 'Enter number of num_packets'
    read(*,*) num_packets

    ! write(*,*) 'Enter min_tau, max_tau, length of tau'
    ! read(*,*) min_tau, max_tau, len_tau

    len_tau = 7
    
    allocate(tau_list(7))
    allocate(total_scatter(len_tau, num_packets))
    allocate(mean_arr(len_tau))
    allocate(length_travelled(len_tau, num_packets))

    ! write(*,*) 'You have selected: ', num_packets, ' packets, with ', min_tau, max_tau, len_tau, ' as the parameters for tau'
    write(*,*) 'Program running!'

    ! call linspace(min_tau, max_tau, tau_list)
    tau_list = [0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]


    do j = 1, size(tau_list)
        tau = tau_list(j)

        do i = 1, num_packets
            r = 0.0_dp
            call zeros(pos)
            call zeros(nhat)
            
            length_travelled(j, i) = 0

            do while ((r .le. r_max) .and. (r .ge. 0.0_dp))
                call mc_emit(nhat)
                
                leng = mc_gen_L(tau,r_max)

                call mc_update(pos, nhat, leng)
                total_scatter(j,i) = total_scatter(j,i) + 1
                length_travelled(j, i) = length_travelled(j, i) + leng

                r = norm(pos)
            end do

        end do
        write(*,*) j
    end do

    do i = 1, size(tau_list)
        mean_arr(i) = mean(total_scatter(i,:))
    end do
    
    open(newunit = u, file = "Data/total_scatters.txt", status = "replace")
    do i = 1, size(tau_list)
        write(u,*) tau_list(i),',' , mean_arr(i)
    end do

    open(newunit = u, file = "Data/length_travelled.txt", status = "replace")
    do j = 1, size(tau_list)
        do i = 1, num_packets
            write(u,*) tau_list(j) ,',' , length_travelled(j,i)/c
        end do
    end do    

end program
    