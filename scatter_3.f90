program monte_carlo
      use types
      use mcrt_functions

      implicit none
      ! allocate variables and declare things
      real(dp), allocatable :: tau_list(:), total_scatter(:, :), mean_arr(:), angle_bins(:,:,:), position(:,:), ini_pos(:,:)
      real(dp) :: pos(3), nhat(3), old_pos(3)
      real(dp) :: r_max = 1.0, min_wt = 1e-3, albedo = 1.0_dp, p = 0.1
      real(dp) :: leng, r, tau, ran,  wt, exit_theta, exit_phi
      integer :: num_packets, len_tau, j, i, nbins = 200, u, phi_loc, theta_loc, k

      ! Parameters for the MC run
      num_packets = 1000000

      ! configure the tau that we want to see

      tau_list = [0.1, 1.0, 5.0, 10.0, 20.0, 50.0]
      len_tau = size(tau_list)

      ! allocating the memory for things
      allocate( total_scatter( len_tau, num_packets ) )
      allocate( mean_arr( len_tau ) )
      allocate( angle_bins( len_tau, nbins, nbins ) )
      allocate( position( num_packets, 3 ) )
      allocate( ini_pos( num_packets, 3 ) )

      ! begin program
      write(*,*) 'Program running!'

      do j = 1, size(tau_list)
      
            do i = 1, num_packets
                  wt = 1.0
                  ! r = 0.0_dp
                  tau = tau_list(j)

                  call blob_emit_bias(nhat)
                  call blob_initial_pos(pos, r_max) ! try to fix the separatrix issue and revert to old method?

                  leng = blob_initial_length(pos, nhat, r_max)
                  call mc_update(pos, nhat, leng)

                  ! force the particle to appear inside the sphere to trigger the MC simulation
                  ! this was done because of rounding errors 
                  r = 0.0_dp
                  do while ((r < r_max) .and. (r >= 0.0_dp))
                        ! provided we are within the sphere
                        old_pos = pos 
                        ! store the old position

                        ! generate the length
                        leng = mc_gen_L(tau, r_max)

                        ! create an nhat vector
                        if (r .ne. 0.0_dp) then ! convert into bool
                              call mc_emit(nhat)
                        end if

                        ! make it walk along the nhat vector a certain length leng
                        call mc_update(pos, nhat, leng)

                        ! if (i == 1) then
                        !       write(*,*) pos
                        ! end if
                      
                        ! check if it will remain in the sphere
                        r = norm(pos)

                        if (r < r_max) then
                              ! calculate if we want to scatter further based on the albedo
                              ! implemented the russian roulette style giving packets an extra life before speaking to Kenny
                              ! he says it's not needed in this HW but I've included it because my simulations were run with it
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
                                    ! we might want to consider the fact that its then scattered straight through
                                    continue
                              end if

                        else if (r >= r_max) then
                              ! if it's outside the sphere then use the geometry to allow us to calculate the "true" length that it's travelled

                              leng = edge_length(old_pos, nhat, r_max)

                              ! now walk the particle to that location on the sphere, not to its final location
                              call mc_update(old_pos, nhat, leng)

                              ! calculate the exit angles in the far field based on the direction that the particle is going along
                              call exit_angles(nhat, exit_theta,  exit_phi)

                              ! generate a location in phi and theta as described in report
                              phi_loc = int( exit_phi/(pi)*nbins/2 + nbins/2 )
                              theta_loc = int( cos(exit_theta) * nbins/2 + nbins/2 +1)
                              ! bin the in the array 
                              angle_bins(j, phi_loc, theta_loc) = angle_bins(j, phi_loc, theta_loc) + 1
                        end if
                  end do
            end do
            write(*,*) tau_list(j)
      end do

      open(newunit = u, file = "Q3/Data/angle_binning.txt", status = "replace")
      do k = 1, len_tau
              do j = 1, nbins
                  do i = 1, nbins
                        write(u,*) tau_list(k), ',', j,',',i , ',',angle_bins(k, j, i)
                  end do
            end do    
      end do
end program