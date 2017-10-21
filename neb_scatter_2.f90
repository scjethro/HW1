program monte_carlo
    use types
    use mcrt_functions
    
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), length_travelled(:, :), images(:,:,:)
    real(dp) :: pos(3), nhat(3), old_pos(3)
    real(dp) :: r_max = 200.0*AU, min_wt = 1e-3, albedo = 1.0_dp, luminosity = 1.0, p = 0.1, etheta = 0, ephi = 0 ! CORRECT THESE
    real(dp) :: xim, yim, pw_wt, leng, r, tau, ran,  wt
    integer :: num_packets, len_tau, j, i, u, k, xbins = 500, ybins = 500

    ! Parameters for the MC run

    ! write(*,*) 'Enter number of num_packets'
    ! read(*,*) num_packets
    num_packets = 100000

    ! leftover from configurable tau values

    ! write(*,*) 'Enter min_tau, max_tau, length of tau'
    ! read(*,*) min_tau, max_tau, len_tau
    
    ! configure the tau that we want to see

    len_tau = 4
    ! call linspace(min_tau, max_tau, tau_list)
    ! tau_list = [0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    tau_list = [0.1, 1.0, 5.0, 10.0]
    
    ! allocating the memory for things
    ! allocate( tau_list( len_tau ) )
    allocate( total_scatter( len_tau, num_packets ) )
    allocate( mean_arr( len_tau ) )
    allocate( length_travelled( len_tau, num_packets ) )
    allocate( images( len_tau, xbins, ybins ) )

    ! begin program
    write(*,*) 'Program running!'

    do j = 1, size(tau_list)

        do i = 1, num_packets

            wt = luminosity/(num_packets)
            r = 0.0_dp
            tau = tau_list(j)
            call zeros(pos)
            call zeros(old_pos)
            call mc_emit(nhat)

            length_travelled(j, i) = 0.0

            do while ((r < r_max) .and. (r .ge. 0.0_dp))
                ! provided we are within the sphere
                old_pos = pos 
                ! store the old position

                ! if it is the first scatter than you force it to scatter
                if (total_scatter(j,i) .eq. 0) then
                    leng = mc_gen_first_L(tau,r_max)
                    wt = wt*(1-exp(-tau))
                else
                    leng = mc_gen_L(tau,r_max)
                end if

                ! NEE Image Binning MAKE INTO SUBROUNTINE

                pw_wt = gen_pw_wt(wt, pos, nhat, r_max, tau)

                call image_calculate(pos, etheta, ephi, xim, yim)

                call image_bin(xim, yim, xbins, ybins, images, pw_wt, r_max, j)

                ! normal MC stuff

                ! create an nhat vector
                call mc_emit(nhat)
                ! make it walk along the nhat vector a certain length leng
                call mc_update(pos, nhat, leng)
                
                ! check if it will remain in the sphere
                r = norm(pos)

                if (r < r_max) then
                    ! update the length travelled if we remain in the sphere
                    length_travelled(j, i) = length_travelled(j, i) + leng
                    
                    ! calculate if we want to scatter further based on the albedo
                    if (ran .lt. albedo) then
                        ! update the weight
                        wt = albedo * wt
                        ! let it scatter
                        total_scatter(j,i) = total_scatter(j,i) + 1

                        ! then do the killing of packets using the formula in the notes
                        if (wt <= min_wt) then
                            ! generate a random number
                            call RANDOM_NUMBER(ran)
                            
                            ! either kill it or update the weight
                            if (ran .lt. p) then
                                wt = wt/p
                            else
                                continue
                            end if
                        end if  
                    ! if we don't want to scatter then just let it continue without doing anything
                    else
                        continue
                    end if
                else if (r .ge. r_max) then
                    ! if it's outside the sphere then use the geometry to allow us to calculate the "true" length that it's travelled
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

    open(newunit = u, file = "Q2/Data_nb2/image_data.txt", status = "replace")
    do k = 1, len_tau
        do j = 1, xbins
            do i = 1, ybins
                write(u,*) tau_list(k), ',', j, ',' ,i, ',', images(k,j,i)
            end do
        end do    
    end do
    
    ! open(newunit = u, file = "Data_nb2/total_scatters.txt", status = "replace")
    ! do i = 1, size(tau_list)
    !     write(u,*) tau_list(i),',' , mean_arr(i)
    ! end do

    ! open(newunit = u, file = "Data_nb2/length_travelled.txt", status = "replace")
    ! do j = 1, size(tau_list)
    !     do i = 1, num_packets
    !         write(u,*) tau_list(j) ,',' , length_travelled(j,i)/c
    !     end do
    ! end do    

end program
    