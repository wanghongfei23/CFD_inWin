
Subroutine PP_limiter_Flux_1D_Euler ( Nx, lambda, qPrim, barF_SP, hatF_PP, hatF_HO )
use CaseSetup, only: nEq
implicit none
!integer,parameter:: Limiter_PP_Type = 1     !> Form_1 (traditional,   CFL<=0.5)
integer,parameter:: Limiter_PP_Type = 2     !> Form_2 (LF type only,  CFL<=1.0)
!integer,parameter:: Limiter_PP_Type = 3     !> Form_3 (HLLC included, CFL<=0.5)
integer,INTENT(IN ):: Nx
real,INTENT(IN   ):: lambda
real,INTENT(IN   )::   qPrim(nEq,Nx)
real,INTENT(IN   ):: barF_SP(nEq,Nx)
real,INTENT(IN   ):: hatF_PP(nEq,Nx+1)
real,INTENT(INOUT):: hatF_HO(nEq,Nx+1)
!> Var
real,dimension(nEq,Nx):: qCons
real,dimension(nEq,Nx):: WL_PP, WR_PP
real,dimension(nEq,Nx):: WL_HO, WR_HO
!> Tool 
integer:: iFP, iSP
real:: theta, theta_R, theta_L
real:: rho_PP, pp_PP
real:: rho_HO, pp_HO
real,parameter:: lim_eps = 1E-13
    
    !> 1st-order PP Splitting
    Do iSP = 1,Nx
        !call PrimToCons ( qPrim(:,iSP), qCons(:,iSP) )

        IF( Limiter_PP_Type == 1 )Then
            WL_PP(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_PP(:,iSP  )  !> Form_1 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2}
            WR_PP(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_PP(:,iSP+1)  !> Form_1 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2}
        ElseIF( Limiter_PP_Type == 2 )Then
            WL_PP(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_PP(:,iSP  ) - lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2} - lam*F_{i}
            WR_PP(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_PP(:,iSP+1) + lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2} + lam*F_{i}
        ElseIF( Limiter_PP_Type == 3 )Then
            WL_PP(:,iSP) = qCons(:,iSP) + 2*lambda*( hatF_PP(:,iSP  ) - barF_SP(:,iSP) )  !> Form_3 ==>  W^L_{i} = U_{i} + 2*lam*( F_{i-1/2} - F_{i} )
            WR_PP(:,iSP) = qCons(:,iSP) - 2*lambda*( hatF_PP(:,iSP+1) - barF_SP(:,iSP) )  !> Form_3 ==>  W^R_{i} = U_{i} - 2*lam*( F_{i+1/2} - F_{i} )
        Else
            write(*,*) 'Split Form NOT Selected!';  stop
        EndIF
    End do
    

    !>---------------------------------------->
    Do iSP = 1,Nx
        IF( Limiter_PP_Type == 1 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_HO(:,iSP  )  !> Form_1 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2}
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_HO(:,iSP+1)  !> Form_1 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2}
        ElseIF( Limiter_PP_Type == 2 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_HO(:,iSP  ) - lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2} - lam*F_{i}
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_HO(:,iSP+1) + lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2} + lam*F_{i}
        ElseIF( Limiter_PP_Type == 3 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*( hatF_HO(:,iSP  ) - barF_SP(:,iSP) )  !> Form_3 ==>  W^L_{i} = U_{i} + 2*lam*( F_{i-1/2} - F_{i} )
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*( hatF_HO(:,iSP+1) - barF_SP(:,iSP) )  !> Form_3 ==>  W^R_{i} = U_{i} - 2*lam*( F_{i+1/2} - F_{i} )
        Else
            write(*,*) 'Split Form NOT Selected!';  stop
        EndIF
    End do
    
    Do iFP = 1,Nx+1
        theta_R = 1.0
        theta_L = 1.0
        
        !> Check WR_{iFP}
        IF( iFP /= 1 )Then  !> Note: qConsR_{0} does NOT exist
            rho_PP = WR_PP(1,iFP-1)   !> F_{i-1/2} --> W^R_{i-1}
            rho_HO = WR_HO(1,iFP-1)   !> F_{i-1/2} --> W^R_{i-1}
        
            IF( rho_PP < 0.0 )Then
                theta_R = 0.0;      write(*,*) 'Warning: Negative PP_density in Flux!'; write(*,*) rho_PP; write(*,*)''
            ElseIF( rho_HO < lim_eps )Then
                theta_R = ( lim_eps - rho_PP  )/( rho_HO - rho_PP )
            EndIF
        EndIF
        
        !> Check WL_{iFP}
        IF( iFP /= Nx+1 )Then   !> Note: qConsL_{Nx+1} does NOT exist
            rho_PP = WL_PP(1,iFP)   !> F_{i-1/2} --> W^L_{i}
            rho_HO = WL_HO(1,iFP)   !> F_{i-1/2} --> W^L_{i}
        
            IF( rho_PP < 0.0 )Then
                theta_L = 0.0;      write(*,*) 'Warning: Negative PP_density in Flux!'; write(*,*) rho_PP; write(*,*)''
            ElseIF( rho_HO < lim_eps )Then
                theta_L = ( lim_eps - rho_PP  )/( rho_HO - rho_PP )
            EndIF
        EndIF
        
        !> Fix density
        theta = min( theta_R, theta_L )
        hatF_HO(:,iFP) =  (1.0-theta)*hatF_PP(:,iFP)  +  theta*hatF_HO(:,iFP)
    End Do
    
    
    !>---------------------------------------->
    Do iSP = 1,Nx
        IF( Limiter_PP_Type == 1 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_HO(:,iSP  )  !> Form_1 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2}
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_HO(:,iSP+1)  !> Form_1 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2}
        ElseIF( Limiter_PP_Type == 2 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*hatF_HO(:,iSP  ) - lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^L_{i} = U_{i} + 2*lam*F_{i-1/2} - lam*F_{i}
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*hatF_HO(:,iSP+1) + lambda*barF_SP(:,iSP)  !> Form_2 ==>  W^R_{i} = U_{i} - 2*lam*F_{i+1/2} + lam*F_{i}
        ElseIF( Limiter_PP_Type == 3 )Then
            WL_HO(:,iSP) = qCons(:,iSP) + 2*lambda*( hatF_HO(:,iSP  ) - barF_SP(:,iSP) )  !> Form_3 ==>  W^L_{i} = U_{i} + 2*lam*( F_{i-1/2} - F_{i} )
            WR_HO(:,iSP) = qCons(:,iSP) - 2*lambda*( hatF_HO(:,iSP+1) - barF_SP(:,iSP) )  !> Form_3 ==>  W^R_{i} = U_{i} - 2*lam*( F_{i+1/2} - F_{i} )
        Else
            write(*,*) 'Split Form NOT Selected!';  stop
        EndIF
    End do
    
    
    Do iFP = 1,Nx+1
        theta_R = 1.0
        theta_L = 1.0
        
        !> Check WR_{iFP}
        IF( iFP /= 1 )Then  !> Note: qConsR_{0} does NOT exist
            call Cal_pressure_1D ( WR_PP(:,iFP-1), pp_PP )   !> F_{i-1/2} --> W^R_{i-1}
            call Cal_pressure_1D ( WR_HO(:,iFP-1), pp_HO )   !> F_{i-1/2} --> W^R_{i-1}
        
            IF( pp_PP < 0.0 )Then
                theta_R = 0.0;      !write(*,*) 'Warning: Negative PP_pressure_R in Flux!  index = ',iFP-1; write(*,*) WR_PP(:,iFP-1); write(*,*) pp_PP; write(*,*)''
            ElseIF( pp_HO < lim_eps )Then
                theta_R = ( lim_eps - pp_PP  )/( pp_HO - pp_PP )
            EndIF
        EndIF
        
        !> Check WL_{iFP}
        IF( iFP /= Nx+1 )Then   !> Note: qConsL_{Nx+1} does NOT exist
            call Cal_pressure_1D ( WL_PP(:,iFP), pp_PP )   !> F_{i-1/2} --> W^L_{i}
            call Cal_pressure_1D ( WL_HO(:,iFP), pp_HO )   !> F_{i-1/2} --> W^L_{i}
        
            IF( pp_PP < 0.0 )Then
                theta_L = 0.0;      !write(*,*) 'Warning: Negative PP_pressure_L in Flux!  index = ',iFP; write(*,*) WL_PP(:,iFP); write(*,*) pp_PP; write(*,*)''
            ElseIF( pp_HO < lim_eps )Then
                theta_L = ( lim_eps - pp_PP  )/( pp_HO - pp_PP )
            EndIF
        EndIF
        
        !> Fix density
        theta = min( theta_R, theta_L )
        hatF_HO(:,iFP) =  (1.0-theta)*hatF_PP(:,iFP)  +  theta*hatF_HO(:,iFP)
    End Do
End Subroutine
!------------------------------------------------------------------------------------------<








!--------------------------------------------------------------
Subroutine Cal_pressure_1D ( cons, pp )
use Parameters, only: vga   !< (GAMMA - 1.0)
implicit none
real,INTENT(IN ):: cons(3)
real,INTENT(OUT):: pp
!real:: E,rho,rhoU
!    rho  = cons(1)
!    rhoU = cons(2)
!    E    = cons(3)
    pp = vga*( cons(3) - 0.5*cons(2)**2/cons(1) )
End Subroutine 
