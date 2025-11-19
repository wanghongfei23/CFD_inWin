
Subroutine PP_limiter_ConsVar_1D_Euler( Nx, qL_PP, qL_HO )
use CaseSetup, only: nEq
implicit none
integer,INTENT(IN ):: Nx
real(kind=8),INTENT(IN   ):: qL_PP(nEq,Nx+1)    !> 1st-order interpolation
real(kind=8),INTENT(INOUT):: qL_HO(nEq,Nx+1)    !> High-order interpolation
!> Tool 
integer:: iFP
real(kind=8):: theta
real(kind=8):: rho_PP, pp_PP
real(kind=8):: rho_HO, pp_HO
real(kind=8),parameter:: lim_eps = 1E-13


    !> �����ܶ�
    !>---------------------------------------->
    Do iFP = 1,Nx+1
        rho_PP = qL_PP(1,iFP)
        rho_HO = qL_HO(1,iFP)
        
        IF( rho_PP < 0.0 )Then
            theta = 0.0; write(*,*) 'Warning: Negative PP_density in Cons!'; write(*,*) rho_PP; write(*,*)''
        ElseIF( rho_HO < lim_eps )Then
            theta = ( lim_eps - rho_PP  )/( rho_HO - rho_PP )
        Else
            theta = 1.0
        EndIF
        
        qL_HO(:,iFP) =  (1.0-theta)*qL_PP(:,iFP)  +  theta*qL_HO(:,iFP)
    End Do
    !<----------------------------------------<
    
    
    !> ����ѹ��
    !>---------------------------------------->
    Do iFP = 1,Nx+1
        pp_PP = qL_PP(nEq,iFP)
        pp_HO = qL_HO(nEq,iFP)
        
        IF( pp_PP < 0.0 )Then
            theta = 0.0; write(*,*) 'Warning: Negative PP_pressure in Cons!'; write(*,*) pp_PP; write(*,*)''
        ElseIF( pp_HO < lim_eps )Then
            theta = ( lim_eps - pp_PP  )/( pp_HO - pp_PP )
        Else
            theta = 1.0
        EndIF
        
        qL_HO(:,iFP) =  (1.0-theta)*qL_PP(:,iFP)  +  theta*qL_HO(:,iFP)
    End Do
    !<----------------------------------------<
        
End Subroutine
!!======================================================================================