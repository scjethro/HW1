program monte_carlo
    use types
    use mcrt_functions
    
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:)
    real(dp) :: pos(3), nhat(3)
    real(dp) :: leng, r_max = 1.0, r, tau, max_tau, min_tau
    integer :: num_packets, len_tau, j, i, u

    write(*,*) 'Enter number of num_packets'
    read(*,*) num_packets

    write(*,*) 'Enter min_tau, max_tau, length of tau'
    read(*,*) min_tau, max_tau, len_tau

    allocate(tau_list(len_tau))
    allocate(total_scatter(len_tau, num_packets))
    allocate(mean_arr(len_tau))

    write(*,*) 'You have selected: ', num_packets, ' packets, with ', min_tau, max_tau, len_tau, ' as the parameters for tau'
    write(*,*) 'Program running!'

    call linspace(min_tau, max_tau, tau_list)
    
    do j = 1, len_tau
        tau = tau_list(j)

        do i = 1, num_packets
            r = 0.0_dp
            call zeros(pos)
            call zeros(nhat)
            
            do while ((r .le. 1.0_dp) .and. (r .ge. 0.0_dp))
                call mc_emit(nhat)
                
                leng = mc_gen_L(tau,r_max)

                call mc_update(pos, nhat, leng)
                total_scatter(j,i) = total_scatter(j,i) + 1
                r = norm(pos)
            end do
        end do
        write(*,*) j
    end do

    do i = 1, len_tau
        mean_arr(i) = mean(total_scatter(i,:))
    end do
    
    open(newunit = u, file = "output_mc.txt", status = "replace")
    do i = 1, len_tau
        write(u,*) tau_list(i),',' , mean_arr(i)
    end do    
    
    deallocate(tau_list)
    deallocate(total_scatter)
    deallocate(mean_arr)

end program
    