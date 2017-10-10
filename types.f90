module types
    implicit none
    
    integer, parameter :: dp = kind(0.d0)
    
    real(dp), parameter :: pi = 4*atan(1.0)
    real(dp), parameter :: c = 299792458.0
    real(dp), parameter :: AU = 1.496e11
    real(dp), parameter :: lsolar = 3.846e26
    real(dp), parameter :: pc = 3.086e16

end module
    