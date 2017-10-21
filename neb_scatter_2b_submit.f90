program monte_carlo
    use types
    use mcrt_functions

    implicit none
    ! allocate variables and declare things
    real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), length_travelled(:, :), image_grid(:,:,:), flux_grid(:,:,:)
    real(dp) :: pos(3), nhat(3), old_pos(3), nhat_view(3)
    real(dp) :: r_max = 200.0*AU, min_wt = 1e-3, albedo = 1.0_dp, luminosity = 0.1*lsolar, p = 0.1, etheta = 0.0, ephi = 0.0, distance = 140.0*pc ! CORRECT THESE
    real(dp) :: xim, yim, pw_wt, leng, r, tau, ran,  wt, flux_wt, flux_multiplier, intensity_multiplier, area_pp
    integer :: num_packets, len_tau, j, i, u, k, xbins = 500, ybins = 500

    ! Parameters for the MC run
    num_packets = 10000000

    ! define the flux multiplier
    flux_multiplier = luminosity/( num_packets * distance**2 )

    ! define the intensity multiplier
    area_pp = (r_max/xbins)**2
    intensity_multiplier = luminosity/( distance**2 * num_packets * ( area_pp/(distance**2) ) )

    ! configure the tau that we want to see

    tau_list = [0.1, 5.0, 10.0, 50.0]
    len_tau = size(tau_list)
    
    ! allocating the memory for things
    allocate( total_scatter( len_tau, num_packets ) )
    allocate( mean_arr( len_tau ) )
    allocate( length_travelled( len_tau, num_packets ) )
    allocate( image_grid( len_tau, xbins, ybins ) )
    allocate( flux_grid( len_tau, xbins, ybins ) )

    ! creating an image plane and an image nhat vector that can be used later
    call initialise_image_plane(nhat_view, etheta, ephi)

    write(*,*) 'Program running!'

    do j = 1, size(tau_list)
        ! loop over tau

        do i = 1, num_packets
            ! loop over packets

            ! set location to zero again
            r = 0.0_dp
            ! define tau to make it easier when referencing
            tau = tau_list(j)
            ! make sure everything is zero
            call zeros(pos)
            call zeros(old_pos)
            call mc_emit(nhat)
            ! give the packet an initial weight
            wt = (1-exp(-tau))
            ! set the length travelled to zero
            length_travelled(j, i) = 0.0

            do while ((r < r_max) .and. (r .ge. 0.0_dp))
                ! provided we are within the sphere
                old_pos = pos 
                ! store the old position

                ! if it is the first scatter than you force it to scatter
                if (total_scatter(j,i) .eq. 0) then
                    leng = mc_gen_first_L(tau,r_max)
                else
                    leng = mc_gen_L(tau,r_max)
                end if

                ! NEE image binning
                ! generate the peeling off weight and the flux weight as described in the lectures
                pw_wt = gen_pw_wt(wt, pos, nhat_view, r_max, tau)
                flux_wt = gen_flux_wt(wt, pos, nhat_view, r_max, tau)

                ! calculate the image locations for the packets
                call image_calculate(pos, etheta, ephi, xim, yim)
                ! bin them in their respective arrays according to their weight
                call image_bin(xim, yim, xbins, ybins, image_grid, pw_wt, r_max, j)
                call image_bin(xim, yim, xbins, ybins, flux_grid, flux_wt, r_max, j)

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

                    ! could include the russian roulette here if Kenny wanted
                else if (r .ge. r_max) then
                    ! if it's outside the sphere then use the geometry to allow us to calculate the "true" length that it's travelled
                    leng = edge_length(old_pos, nhat, r_max)
                    length_travelled(j, i) = length_travelled(j, i) + leng
                end if
            end do
        end do
        write(*,*) j
    end do

    open(newunit = u, file = "Q2/Data_nb2_b/image_data.txt", status = "replace")
    do k = 1, len_tau
        do j = 1, xbins
            do i = 1, ybins
                write(u,*) tau_list(k), ',', j, ',' ,i, ',', image_grid(k,j,i)
            end do
        end do    
    end do

    open(newunit = u, file = "Q2/Data_nb2_b/flux_data.txt", status = "replace")
    do k = 1, len_tau
        do j = 1, xbins
            do i = 1, ybins
                write(u,*) tau_list(k), ',', j, ',' ,i, ',', flux_grid(k, j, i)*flux_multiplier
            end do
        end do    
    end do

    open(newunit = u, file = "Q2/Data_nb2_b/intensity_data.txt", status = "replace")
    do k = 1, len_tau
        do j = 1, xbins
            do i = 1, ybins
                write(u,*) tau_list(k), ',', j, ',' ,i, ',', flux_grid(k, j, i)*intensity_multiplier
            end do
        end do    
    end do
    
end program
    