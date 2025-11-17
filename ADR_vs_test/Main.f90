Module CaseSetup
implicit none
    !=======================================
    logical:: useWCNS_MR        = .false.    !> �Ƿ�ʹ��WCNS-MR
    logical:: useWCNS_JS        = .false.    !> �Ƿ�ʹ��WCNS-JS
    logical:: useWCNS_Z         = .false.    !> �Ƿ�ʹ��WCNS-Z 
    logical:: useTCNS_Cong      = .false.    !> �Ƿ�ʹ��TCNS-Cong only 5th order
    logical:: useTCNS           = .false.    !> �Ƿ�ʹ��TCNS only 5th order
    !=======================================
    
    integer,parameter:: nEq = 1
    integer:: Order != 11
    integer:: OInt  != Order
    integer:: ODiff != Order+1
End Module


!!------------------------!!
Module Variables_Zone1
implicit none
    integer:: Nx
    real:: Hx
    real,allocatable,dimension( : ) :: SP_X, FP_X
    real,allocatable,dimension( : ) :: u0, ut
    real,allocatable,dimension(:,:) :: ADR
End Module
!!--------------------------------------------------------------------------------------




!!======================================================================================
program ADR_Analysis
use CaseSetup
use Variables_Zone1
use Parameters, only: PI
implicit none
integer:: iSP, nn, rr
real,parameter:: T0 = 0.0
real,parameter:: Tt = 1.E-5 !> tau
real,parameter:: Length = 1.0
real:: lambda_n   !> Wave length
real::    phi_n   !> Wave number
real:: Re0, Im0
real:: Ret, Imt
real:: dx, dt
complex*16 :: cmpLn, cmpS0, cmpSt
!----------------------------------------
    Nx = 499
    Hx = Length/(Nx)
    
    Allocate(  SP_X( 0:Nx )   )
    Allocate(    u0( 0:Nx )   )
    Allocate(    ut( 0:Nx )   )
    Allocate( ADR( 3,0:Nx/2 ) )
    
    !> Mesh
    Do iSP = 0,Nx
        SP_X(iSP) = iSP*Hx
    End Do
    
    
    !Do rr = 3,3
        Order = 5
        OInt  = Order
        ODiff = Order+1
        !================================> ADR
        Do nn = 1,Nx/2
            print *, "nn = ",nn
            print *, "Nx = ",Nx
        
            !> Wavelength
            lambda_n = 1.d0/nn
               phi_n = 2.d0*pi*nn/nx
    
            !> Initial values
            Do iSP = 0,Nx
                u0(iSP) = cos( 2*PI/lambda_n *SP_X(iSP) )
            End Do
        
            !> Fourier Transformation (t=0)
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
        
        
            !> Fourier Transformation (t=tau)
            Ret = 0.d0
            Imt = 0.d0
            do iSP = 0,Nx
	            Ret = Ret + ut(iSP)*cos(dble(iSP)*phi_n)
                Imt = Imt - ut(iSP)*sin(dble(iSP)*phi_n)
	        enddo
            Ret = Ret/dble(Nx)
            Imt = Imt/dble(Nx)

            !> Compute ADR
            cmpS0 = cmplx(Re0,Im0)
            cmpSt = cmplx(Ret,Imt)
            cmpLn = cmpSt/cmpS0
            cmpLn = zlog(cmpLn)  !<��Ȼ����
            dx = Hx
            dt = Tt
            ADR(1,nn) =   phi_n
            ADR(2,nn) = - AIMAG(cmpLn) * dx/dt
            ADR(3,nn) =    REAL(cmpLn) * dx/dt
        End Do
        !================================>
        call Output
    !End Do
    
stop
End
    
    
    

!&=========================================================================|
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
    EndIF

    write(C2, '(I8)') Order
    
    write(C3, '(I8)') Nx
    
    
    file_Name = './ADR_N_lam=1.plt'
    
    open( 99, file=trim(file_Name) ) 
        write(99,*) 'VARIABLES="omega" "Re" "Im" '
        Do nn = 1, nx/2
            write(99,*) ADR(:,nn)
        End Do
    close(99)
End subroutine