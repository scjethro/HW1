module mcrt_functions
    use types
    implicit none


    contains
    
        real(dp) function sample_exp()
            real(dp) :: rand
 
            call RANDOM_NUMBER(rand)
            sample_exp = -log(rand)
        end function
        
        real(dp) function mc_gen_L(tau_max, r_max) result(tau_out)
            real(dp) :: tau, tau_max, r_max
            
            tau = sample_exp()
            tau_out = tau/tau_max*r_max
        end function

        real(dp) function mc_gen_first_L(tau_e, r_max) result(tau_out)
            real(dp) :: tau, tau_e, r_max
            
            tau = mc_first_depth(tau_e)
            tau_out = tau/tau_e*r_max
        end function
        
        real(dp) function norm(vector)
            real(dp), intent(in) :: vector(:)
            real(dp) :: total
            integer :: i
            
            total = 0.0_dp
  
            do i = 1,size(vector)
                total = total + vector(i)**2
            end do
            
            norm = sqrt(total)
        end function
        
        real(dp) function mean(array)
            real(dp), intent(in) :: array(:)
            real(dp) :: total
            
            total = sum(array)
            
            mean = total/size(array)
        end function

        real(dp) function variance(array)
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
            real(dp) , intent(in) :: tau_edge
            real(dp) :: rand
            
            call RANDOM_NUMBER(rand)
            tau = -log(1-rand*(1-exp(-tau_edge)))
        end function

        real(dp) function edge_length(pos, nhat, r_max)
            real(dp), intent(in) :: pos(3), nhat(3), r_max
            real(dp) :: dot, r

            dot = dot_product(pos, nhat)
            r = norm(pos)
            edge_length = -dot + sqrt(dot**2 - (r**2 - r_max**2))
        end function

        subroutine image_calculate(position, theta, phi, x, y)
            real(dp), intent(in) :: position(3), theta, phi
            real(dp), intent(out) :: x, y

            x = position(2)*cos(phi) - position(1)*sin(phi)
            y = position(3)*sin(theta) - position(2)*cos(theta)*sin(phi) - position(1)*cos(theta)*cos(phi)

        end subroutine

        subroutine exit_angles(nhat, theta, phi)
            real(dp), intent(in) :: nhat(3)
            real(dp), intent(out) :: theta, phi

            theta = atan(nhat(2)/nhat(1))
            phi = atan(sqrt(nhat(1)**2 + nhat(2)**2)/nhat(1))
        end subroutine

        subroutine mc_emit(n)
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
            real(dp), intent(out) :: old_pos
            real(dp), intent(in) :: n, length
            
            old_pos = old_pos + n*length
        end subroutine

        subroutine linspace(start, finish, r)
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
            real(dp), intent(out) :: array
            array = 0
        end subroutine

        subroutine image_bin(xim, yim, xbins, ybins, image_grid, pw_wt, r_max, tau)
            real(dp), intent(in) :: xim, yim, pw_wt, r_max
            integer :: xloc, yloc, tau, xbins, ybins
            real(dp) :: image_grid(:,:,:)

            xloc = int((xim/r_max)*xbins/2 + xbins/2 + 1)
            yloc = int((yim/r_max)*xbins/2 + ybins/2 + 1)

            image_grid(tau, xloc, yloc) = image_grid(tau, xloc, yloc) + pw_wt

        end subroutine
        
end module mcrt_functions