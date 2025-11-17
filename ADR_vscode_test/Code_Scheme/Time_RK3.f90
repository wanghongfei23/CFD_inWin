Subroutine TVD_RK3 ( T0, Tt )
use CaseSetup
use Variables_Zone1
implicit none
real,Intent(IN):: T0, Tt
real:: Ht
real,dimension(0:Nx):: uk
real,dimension(0:Nx):: RHS_1, RHS_2, RHS_3

    Ht = Tt - T0
    !====================================> TVD-RK3
    !---> Step 1
        uk = u0
        call RHS_ghostNode ( uk, RHS_1 )

    !---> Step 2
        uk = u0 + RHS_1 *Ht
        call RHS_ghostNode ( uk, RHS_2 )

    !---> Step 3
        uk = u0 + (RHS_1 + RHS_2) *Ht/4
        call RHS_ghostNode ( uk, RHS_3 )
    !====================================>

    ut = u0 + ( RHS_1 + RHS_2 + 4*RHS_3 ) *Ht/6
End Subroutine

