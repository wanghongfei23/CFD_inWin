!!===============================================!!
!!                 差分：显式                     !!
!!===============================================!!
Subroutine Diff_E4 ( Hx, Stencil, du_SP )
implicit none
integer,parameter:: ODiff = 4
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
real(kind=8):: CoeDif_IN(ODiff)
    CoeDif_IN(1:4) = (/ 1./24,  -9./8,  9./8,  -1./24  /)
    du_SP =  sum( CoeDif_IN(:) * Stencil(:) )/Hx
End Subroutine


Subroutine Diff_E6 ( Hx, Stencil, du_SP )
implicit none
integer,parameter:: ODiff = 6
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
real(kind=8):: CoeDif_IN(ODiff)
    CoeDif_IN(1:6) = (/ -3./640,  25./384,  -75./64,  75./64,  -25./384,  3./640  /)
    du_SP =  sum( CoeDif_IN(:) * Stencil(:) )/Hx
End Subroutine


Subroutine Diff_E8 ( Hx, Stencil, du_SP )
implicit none
integer,parameter:: ODiff = 8
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
real(kind=8):: CoeDif_IN(ODiff)
    CoeDif_IN(1:8) = (/  5./7168,  -49./5120,  245./3072,  -1225./1024,  1225./1024,  -245./3072,  49./5120,  -5./7168  /)
    du_SP =  sum( CoeDif_IN(:) * Stencil(:) )/Hx
End Subroutine


Subroutine Diff_E10 ( Hx, Stencil, du_SP )
implicit none
integer,parameter:: ODiff = 10
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
real(kind=8):: CoeDif_IN(ODiff)
    CoeDif_IN(1:10) = (/  -35./294912,  405./229376,  -567./40960,  735./8192,  -19845./16384,  19845./16384,  -735./8192,  567./40960,  -405./229376,  35./294912  /)
    du_SP =  sum( CoeDif_IN(:) * Stencil(:) )/Hx
End Subroutine


Subroutine Diff_E12 ( Hx, Stencil, du_SP )
implicit none
integer,parameter:: ODiff = 12
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
real(kind=8):: CoeDif_IN(ODiff)
    CoeDif_IN(1:12) = (/  63./2883584,  -847./2359296,  5445./1835008,  -22869./1310720,  12705./131072,  -160083./131072,  160083./131072,  -12705./131072,  22869./1310720,  -5445./1835008,  847./2359296,  -63./2883584  /)
    du_SP =  sum( CoeDif_IN(:) * Stencil(:) )/Hx
End Subroutine