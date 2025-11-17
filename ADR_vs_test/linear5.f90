subroutine Linear_5thi( fh , fp )
    implicit none
    real :: gamma0 , gamma1 , gamma2
    real :: beta(3)
    real :: tau5
    real :: epsilon
    real :: wt0 , wt1 , wt2
    real :: wsum
    real :: f30 , f31 , f32
    real :: fh
    real :: fp(5)
    real :: CT
    integer ::flag
       

    fh = 3.0/128.0*fp(1) - 5.0/32.0*fp(2) + 45.0/64.0*fp(3) + 15.0/32.0*fp(4) - 5.0/128.0*fp(5)  



    endsubroutine