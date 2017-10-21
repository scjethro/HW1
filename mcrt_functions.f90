module mcrt_functions
    use types
    implicit none


    contains
      
        real(dp) function sample_exp()
        ! Sample from the exponential distribution using a random number in the range 0 -> 1
            real(dp) :: rand
 
            call RANDOM_NUMBER(rand)
            sample_exp = -log(rand)
        end function
        
        real(dp) function mc_gen_L(tau_max, r_max) result(tau_out)
            ! generate an optical depth and then scale it according to the maximum tau and the max_radius
            real(dp) :: tau, tau_max, r_max
            
            tau = sample_exp()
            tau_out = tau/tau_max*r_max
        end function

        real(dp) function mc_gen_first_L(tau_e, r_max) result(tau_out)
            ! generate the first optical depth in the case that you want forced first scattering
            real(dp) :: tau, tau_e, r_max
            
            ! we use a separate function to calculate the depth to ensure that the value is bounded in the range 0 -> tau_edge
            tau = mc_first_depth(tau_e)
            tau_out = tau/tau_e*r_max
        end function
        
        real(dp) function norm(vector)
            ! calculate the magnitude/norm of a given vector
            real(dp), intent(in) :: vector(:)
            real(dp) :: total
            integer :: i
            
            total = 0.0_dp
            ! note this is not the most elegant solution but it works
            do i = 1,size(vector)
                total = total + vector(i)**2
            end do
            
            norm = sqrt(total)
        end function
        
        real(dp) function mean(array)
            ! calculate the mean of a 1D array that is passed
            real(dp), intent(in) :: array(:)
            real(dp) :: total
            
            total = sum(array)
            
            mean = total/size(array)
        end function

        real(dp) function variance(array)
            ! calculate the variance of a 1D array that is passed
            real(dp) :: array(:), mean_val, total
            integer :: i

            mean_val = mean(array)
            total = 0.0_dp
            do i = 1, size(array)
                total = (array(i) - mean_val)**2
            end do

            variance = total/(size(array)-1)
        end function

        real(dp) function mc_first_depth(tau_edge) result(tau)
            ! generate the correct optical depth for forced first scatter where we want to ensure the value
            ! is bounded on 0 -> tau_edge
            real(dp) , intent(in) :: tau_edge
            real(dp) :: rand
            call RANDOM_NUMBER(rand)

            tau = -log(1-rand*(1-exp(-tau_edge)))
        end function

        real(dp) function edge_length(pos, nhat, r_max)
            ! calculate the distance from a point within a sphere along some direction nhat to the surface
            real(dp), intent(in) :: pos(3), nhat(3), r_max
            real(dp) :: dot, r
            dot = dot_product(pos, nhat)
            r = norm(pos)

            edge_length = -dot + sqrt(dot**2 - (r**2 - r_max**2))
        end function

        real(dp) function blob_initial_length(pos, nhat, r_max)
            ! for Q3 - generate the length a packet should advance to hit the surface of the blob
            real(dp), intent(in) :: pos(3), nhat(3), r_max
            real(dp) :: dot, r

            dot = dot_product(pos, nhat)
            r = norm(pos)

            blob_initial_length = -dot - sqrt(dot**2 - (r**2 - r_max**2))
        end function

        real(dp) function gen_pw_wt(wt, pos, nhat_view, r_max, tau)
            ! generating a peel weight for the next event estimator that can be used when binning to ensure the packets are given the correct weighting
            real(dp), intent(in) :: pos(:), nhat_view(:), r_max, tau, wt
            real(dp) :: tau_peel

            ! using the equation that we discussed during lectures
            tau_peel = tau*edge_length(pos, nhat_view, r_max)/r_max
            gen_pw_wt = 1.0/(4*pi)*exp(-tau_peel)*wt
        end function

        real(dp) function gen_flux_wt(wt, pos, nhat, r_max, tau)
            ! performing the same operation as above but in this case we don't want the 1/4pi factor so we correct for it
            real(dp), intent(in) :: pos(:), nhat(:), r_max, tau, wt
            real(dp) :: pw_wt

            pw_wt = gen_pw_wt(wt, pos, nhat, r_max, tau)
            gen_flux_wt = 4*pi*pw_wt
        end function

        real(dp) function trapz_int(y_array, x_array) result(integral)
            ! a trapezoidal integration routine that takes an x and y array and returns the trapezoidal numerical integral
            ! unused at the moment
            real(dp), intent(in) :: y_array(:), x_array(:)
            integer :: N, i
            
            if ( size(y_array) == size(x_array) ) then
                N = size(y_array)
            else 
                stop 'You have passed 2 arrays of different size, this is not allowed'
            end if

            integral = 0.0
            
            do i = 1, N-1
                integral = integral + 0.5*(y_array(i+1) + y_array(i))*(x_array(i+1) - x_array(i))
            end do
        end function

        subroutine image_calculate(position, theta, phi, x, y)
            ! calculate the image positions using the expressions discussed in class to ensure that we can bin the packets correctly
            real(dp), intent(in) :: position(3), theta, phi
            real(dp), intent(out) :: x, y

            x = position(2)*cos(phi) - position(1)*sin(phi)
            y = position(3)*sin(theta) - position(2)*cos(theta)*sin(phi) - position(1)*cos(theta)*cos(phi)
        end subroutine

        subroutine exit_angles(nhat, theta, phi)
            ! calculate the azimuth and altitude angles of a specific vector nhat
            real(dp), intent(in) :: nhat(3)
            real(dp), intent(out) :: theta, phi

            theta = acos(nhat(3)/norm(nhat))
            phi = atan2(nhat(2), nhat(1))
        end subroutine

        subroutine mc_emit(n)
            ! the emission subroutine for MC packets, creates two random numbers and uses them to generate two values
            ! theta and phi which can be used to create an nhat vector
            real(dp), intent(inout) :: n(:)
            real(dp) :: random_num1, random_num2, theta, phi
            
            if (size(n) /= 3) stop 'Incorrect size of array'
        
            call RANDOM_NUMBER(random_num1)
            call RANDOM_NUMBER(random_num2)
            
            theta = acos(2*random_num1-1)
            phi = 2*pi*random_num2
            
            n(1) = sin(theta)*cos(phi)
            n(2) = sin(theta)*sin(phi)
            n(3) = cos(theta)
        end subroutine
        
        elemental subroutine mc_update(old_pos, n, length)
            ! update the location of the packets by walking them along the distance that has been calculated for them
            real(dp), intent(out) :: old_pos
            real(dp), intent(in) :: n, length
            
            old_pos = old_pos + n*length
        end subroutine

        subroutine blob_initial_pos(pos, r_max)
            ! generating a plane of uniform illumination for the blob, setting it at an arbitrary illumination location (along the x axis)
            real(dp) :: pos(3), rand1, radius, phi, r_max

            call RANDOM_NUMBER(rand1)
            call RANDOM_NUMBER(radius)

            phi = 2*pi*rand1
            pos(1) = -r_max
            pos(2) = sqrt(radius)*sin(phi)
            pos(3) = sqrt(radius)*cos(phi)
        end subroutine

        subroutine blob_emit_bias(n)
            ! ensuring that the particles are emitted parallel towards the target, this subroutine only generates xhat unit vectors
            real(dp), intent(inout) :: n(:)
            n(1) = 1.0
            n(2) = 0.0
            n(3) = 0.0
        end subroutine

        subroutine linspace(start, finish, r)
            ! a subroutine that I wrote to replicate the function np.linspace() in NumPy, useful I've needed a linear space of numbers
            ! simply calculates the difference in start and finish then creates an appropriate number of elements to fill that space
            real(dp), intent(in) :: start, finish
            real(dp), intent(out) :: r(:)
            real(dp) :: step
            integer :: num, i
            
            num = size(r)
            step = (finish - start)/num
            
            do i = 1, num
                r(i) = start + (i-1)*step
            end do
        end subroutine linspace

        elemental subroutine zeros(array)
            ! a subrountine to generate an array of zeros so that you don't worry about the garbage they're filled with.
            ! pretty much useless since you can do array(:,:) = 0 and set it all at once
            real(dp), intent(out) :: array
            array = 0
        end subroutine

        subroutine image_bin(xim, yim, xbins, ybins, image_grid, pw_wt, r_max, tau)
            ! bin the data depending on locations within the x, y space so that an image can be created
            real(dp), intent(in) :: xim, yim, pw_wt, r_max
            integer :: xloc, yloc, tau, xbins, ybins
            real(dp) :: image_grid(:,:,:)

            xloc = int((xim/r_max)*xbins/2 + xbins/2 + 1)
            yloc = int((yim/r_max)*xbins/2 + ybins/2 + 1)

            image_grid(tau, xloc, yloc) = image_grid(tau, xloc, yloc) + pw_wt
        end subroutine

        subroutine initialise_image_plane(nhat_view, etheta, ephi)
            ! create an image plane with a specific vector nhat so that images can be taken onto that plane.
            ! uses minus signs to ensure that the nhat vector points as it should towards the blob in Q3
            real(dp), intent(in) :: etheta, ephi
            real(dp), intent(inout) :: nhat_view(3)
            
            nhat_view(1) = -sin(etheta)*cos(ephi)
            nhat_view(2) = -sin(etheta)*sin(ephi)
            nhat_view(3) = -cos(etheta)
        end subroutine
        
end module mcrt_functions