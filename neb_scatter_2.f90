program monte_carlo
    use types
    use mcrt_functions
    
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), length_travelled(:, :)
    ! real(dp) :: max_tau, min_tau
    real(dp) :: pos(3), nhat(3), old_pos(3)
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

        do i = 1, num_packets
            r = 0.0_dp
            call zeros(pos)
            call zeros(nhat)
            call zeros(old_pos)
            
            length_travelled(j, i) = 0.0

            do while ((r < r_max) .and. (r .ge. 0.0_dp))
                
                tau = tau_list(j)
                
                old_pos = pos
                
                if (total_scatter(j,i) .eq. 0) then
                    leng = mc_gen_first_L(tau,r_max)
                else
                    leng = mc_gen_L(tau,r_max)
                end if

                call mc_emit(nhat)

                call mc_update(pos, nhat, leng)
                
                r = norm(pos)
                
                if (r < r_max) then
                    total_scatter(j,i) = total_scatter(j,i) + 1
                    length_travelled(j, i) = length_travelled(j, i) + leng
                
                else if (r .ge. r_max) then
                    leng = edge_length(old_pos, nhat, r_max)
                    length_travelled(j, i) = length_travelled(j, i) + leng

                end if
                
            end do

        end do
        write(*,*) j
    end do

    do i = 1, size(tau_list)
        mean_arr(i) = mean(total_scatter(i,:))
    end do
    
    open(newunit = u, file = "Data_nb2/total_scatters.txt", status = "replace")
    do i = 1, size(tau_list)
        write(u,*) tau_list(i),',' , mean_arr(i)
    end do

    open(newunit = u, file = "Data_nb2/length_travelled.txt", status = "replace")
    do j = 1, size(tau_list)
        do i = 1, num_packets
            write(u,*) tau_list(j) ,',' , length_travelled(j,i)/c
        end do
    end do    

end program
    