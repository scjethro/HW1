module types
    implicit none
    ! declaration of various parameters I want to use throughout my code to make things simpler

    ! define the parameter dp to represent real*8/double precision to make the code look cleaner
    integer, parameter :: dp = kind(0.d0)
    
    ! then define another set of paramters that will be useful throughout
    real(dp), parameter :: pi = 4*atan(1.0)
    real(dp), parameter :: c = 299792458.0
    real(dp), parameter :: AU = 1.496e11
    real(dp), parameter :: lsolar = 3.846e26
    real(dp), parameter :: pc = 3.086e16

end module