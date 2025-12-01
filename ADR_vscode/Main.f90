!! ----------------------------------------------!!
Module CaseSetup
implicit none
    ! ====
    ! 王鸿飞开关
    logical:: useWCNS_MR        = .false.    !> 是否使用WCNS-MR
    logical:: useWCNS_JS        = .false.    !> 是否使用WCNS-JS
    logical:: useWCNS_Z         = .false.    !> 是否使用WCNS-Z 
    logical:: useTCNS           = .false.    !> 是否使用TCNS 
    logical:: useTCNS_A         = .false.    !> 是否使用TCNS_A 
    logical:: useTCNS_S         = .false.    !> 是否使用TCNS_S 
    logical:: useTCNS_ASF102    = .true.    !> 是否使用TCNS_ASF102 
    ! ====
    
    integer,parameter:: nEq = 1
    integer:: Order != 11
    integer:: OInt  != Order
    integer:: ODiff != Order+1
End Module
!! ----------------------------------------------!!
!! ----------------------------------------------!!
Module Variables_Zone1
implicit none
    integer:: Nx
    real(kind=8):: Hx
    real(kind=8),allocatable,dimension( : ) :: SP_X, FP_X
    real(kind=8),allocatable,dimension( : ) :: u0, ut
    real(kind=8),allocatable,dimension(:,:) :: ADR
End Module
!! ----------------------------------------------!!

program ADR_Analysis
use CaseSetup
use Variables_Zone1
use Parameters, only: PI
implicit none
integer:: iSP, nn, rr
real(kind=8),parameter:: T0 = 0.0
real(kind=8),parameter:: Tt = 1.E-7 !> tau
real(kind=8),parameter:: Length = 1.0
real(kind=8):: lambda_n   !> 波长
real(kind=8)::    phi_n   !> 波数
real(kind=8):: Re0, Im0
real(kind=8):: Ret, Imt
real(kind=8):: dx, dt
complex*16 :: cmpLn, cmpS0, cmpSt
!----------------------------------------
    Nx = 499 ! 原有的
    ! Nx = 199
    ! Nx = 99
    Hx = Length/(Nx)
    
    Allocate(  SP_X( 0:Nx )   )
    Allocate(    u0( 0:Nx )   )
    Allocate(    ut( 0:Nx )   )
    Allocate( ADR( 3,0:Nx/2 ) )
    
    !> 网格
    Do iSP = 0,Nx
        SP_X(iSP) = iSP*Hx
    End Do
    
    Do rr = 2,6
        Order = 2*rr - 1
        OInt  = Order
        ODiff = Order+1
        !================================> ADR
        Do nn = 1,Nx/2
            print *, "nn = ",nn
            print *, "Nx = ",Nx
        
            !> 波长
            lambda_n = 1.d0/nn
               phi_n = 2.d0*pi*nn/nx
    
            !> 初始值
            Do iSP = 0,Nx
                u0(iSP) = cos( 2*PI/lambda_n *SP_X(iSP) )
            End Do
        
            !> 傅里叶变换 (t=0)
            Re0 = 0.d0
            Im0 = 0.d0
            do iSP = 0,Nx
                Re0 = Re0 + u0(iSP)*cos(dble(iSP)*phi_n)
                Im0 = Im0 - u0(iSP)*sin(dble(iSP)*phi_n)
            enddo
            Re0 = Re0/dble(Nx)
            Im0 = Im0/dble(Nx)
        
        
            !========================> t=0 -> t=tau
            call TVD_RK3 ( T0, Tt )
            !========================>
        
        
            !> 傅里叶变换 (t=tau)
            Ret = 0.d0
            Imt = 0.d0
            do iSP = 0,Nx
                Ret = Ret + ut(iSP)*cos(dble(iSP)*phi_n)
                Imt = Imt - ut(iSP)*sin(dble(iSP)*phi_n)
            enddo
            Ret = Ret/dble(Nx)
            Imt = Imt/dble(Nx)

            !> 计算ADR
            cmpS0 = cmplx(Re0,Im0,kind=8)
            cmpSt = cmplx(Ret,Imt,kind=8)
            cmpLn = cmpSt/cmpS0
            cmpLn = zlog(cmpLn)  !<自然对数
            dx = Hx
            dt = Tt
            ADR(1,nn) =   phi_n
            ADR(2,nn) = - AIMAG(cmpLn) * dx/dt
            ADR(3,nn) =    REAL(cmpLn) * dx/dt
        End Do
        !================================>
        call Output
    End Do
stop
End


! 王鸿飞
subroutine Output
use CaseSetup
use Variables_Zone1
implicit none
character(len=50 ):: C1,C2,C3
character(len=400):: file_Name
integer:: nn

    IF( useWCNS_MR )Then
        C1 = 'WCNS-MR-O'
    ElseIF( useWCNS_JS )Then
        C1 = 'WCNS-JS-O'
    ElseIF( useWCNS_Z )Then
        C1 = 'WCNS-Z-O'
    ElseIF( useTCNS )Then
        C1 = 'TCNS-Z-O'
    ElseIF( useTCNS_A )Then
        C1 = 'TCNS-A-O'
    ElseIF( useTCNS_S )Then
        C1 = 'TCNS-S-O'
    ElseIF( useTCNS_ASF102 )Then
        C1 = 'TCNS-ASF102-O'
    EndIF


    write(C2, '(I8)') Order
    
    write(C3, '(I8)') Nx
    
    
    file_Name = './ADR_N'//trim(adjustl(C3))//'_'//trim(adjustl(C1))//trim(adjustl(C2))//'_lam=1.plt'
    
    open( 99, file=trim(file_Name) ) 
        write(99,*) 'VARIABLES="omega" "Re" "Im" '
        Do nn = 1, nx/2
            write(99,*) ADR(:,nn)
        End Do
    close(99)
End subroutine
