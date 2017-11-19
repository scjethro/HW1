program monte_carlo_Q1
    use types
    use mcrt_functions

    implicit none
    
    ! delcaration statements
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), length_travelled(:, :)
    real(dp) :: pos(3), nhat(3), old_pos(3)
    real(dp) :: leng, r_max = 200.0*AU, r, tau
    integer :: num_packets, len_tau, j, i, u

    ! ask user for input of the number of packets they would like to run
    ! write(*,*) 'Enter number of num_packets'
    ! read(*,*) num_packets
    num_packets = int(1e5)

    ! memory allocation of different variables now that number of packets is known
    len_tau = 20
    allocate(tau_list(len_tau))
    allocate(total_scatter(len_tau, num_packets))
    allocate(mean_arr(len_tau))
    allocate(length_travelled(len_tau, num_packets))

    write(*,*) 'Program running!'

    ! giving a fixed list of tau values to conform to the question
    ! tau_list = [0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    call linspace(0.1_dp, 50.0_dp, tau_list)

    ! begin to loop over the tau values
    do j = 1, size(tau_list)
        ! for ease of use, define a value of tau as just being the current tau that is running
        tau = tau_list(j)

        ! loop over every packet in the list
        do i = 1, num_packets
            ! set the packets location to the origin
            r = 0.0_dp
            ! generate zero arrays that can be populated, as Kenny does in the Grid code
            call zeros(pos)
            call zeros(nhat)
            call zeros(old_pos)
            
            ! set the distance travelled to 0
            length_travelled(j, i) = 0

            ! while the particle remains in the sphere
            do while ((r .le. r_max) .and. (r .ge. 0.0_dp))
                ! store the old position
                old_pos = pos
                ! emit a particle and create a directino
                call mc_emit(nhat)
                ! generate an appropriate length for the particle to travel
                leng = mc_gen_L(tau,r_max)
                ! update the particle location
                call mc_update(pos, nhat, leng)
                ! calculate the length of the r vector
                r = norm(pos)

                if (r < r_max) then
                    ! if the particle remains in the sphere then add a scatter and add the length
                    total_scatter(j,i) = total_scatter(j,i) + 1
                    length_travelled(j, i) = length_travelled(j, i) + leng
                else if (r >= r_max) then
                    ! if the particle has left the sphere then don't add a scatter as it has shot straight out
                    leng = edge_length(old_pos, nhat, r_max)
                    ! increase the length that the particle has travelled based on the direction it was travelling and the location of the spherical surface.
                    length_travelled(j, i) = length_travelled(j, i) + leng
                end if 

            end do

        end do
        ! use this to indicate which tau value we are evaluating, helps to see how far it is along
        write(*,*) j
    end do

    ! loop over the list of tau values and calculate the mean number of scatters
    do i = 1, size(tau_list)
        mean_arr(i) = mean(total_scatter(i,:))
    end do
    
    ! open a file ready for writing and then loop over the data and write it out
    open(newunit = u, file = "Q1/Data_nb1/total_scatters.txt", status = "replace")
    do i = 1, size(tau_list)
        write(u,*) tau_list(i),',' , mean_arr(i)
    end do

    ! open a file ready for writing and then loop over the data and write it out
    open(newunit = u, file = "Q1/Data_nb1/length_travelled.txt", status = "replace")
    do j = 1, size(tau_list)
        do i = 1, num_packets
            write(u,*) tau_list(j) ,',' , length_travelled(j,i)/c
        end do
    end do    

end program