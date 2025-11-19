Subroutine RHS_ghostNode ( u_In, RHS )
use CaseSetup, only: Order
use Variables_Zone1, only: Nx, Hx
implicit none
real(kind=8),INTENT(IN ),dimension(0:Nx):: u_In
real(kind=8),INTENT(OUT),dimension(0:Nx)::  RHS
!integer,parameter:: span = (Order+1)/2
!integer,parameter:: numGNSP = Order
!integer,parameter:: numGNFP = (Order-1)/2
!real(kind=8),dimension( -numGNSP+0: Nx+numGNSP  )::  u_GN  !> SP
!real(kind=8),dimension( -numGNFP+0: Nx+numGNFP+1):: uL_GN  !> FP
integer:: span, numGNSP, numGNFP
real(kind=8),allocatable,dimension(:)::  u_GN, uL_GN
!> Tool 
integer:: iSP, iSPs,iSPe
integer:: iFP, iFPs,iFPe
real(kind=8):: primL, primR, Stencil(Order+1)
real(kind=8):: du_SP
!!--------------------------------------------------------
span = (Order+1)/2
numGNSP = Order
numGNFP = (Order-1)/2
Allocate(  u_GN(-numGNSP+0: Nx+numGNSP  ) )
Allocate( uL_GN(-numGNFP+0: Nx+numGNFP+1) )

    Do iSP = 1,numGNSP
        u_GN(0-iSP) = u_In(Nx-iSP+1)
    End do

    Do iSP = 0,Nx
        u_GN(iSP) = u_In(iSP)
    End do

    Do iSP = 1,numGNSP
        u_GN(Nx+iSP) = u_In(Nx+iSP-1)
    End do

    !===> 插值
    Do iFP = -numGNFP+0, Nx+numGNFP+1
        iSPs = iFP - span
        iSPe = iFP + span-1
        Stencil(:) = u_GN( iSPs:iSPe )
        
        call NLI_ghostNode_my_AS ( Stencil, primL, primR )
        ! call NLI_ghostNode ( Stencil, primL, primR )
        
        uL_GN(iFP) = primL
    End do
    
    !===> 差分
    Do iSP = 0,Nx
        iFPs = iSP - span+1
        iFPe = iSP + span
        Stencil(:) = uL_GN( iFPs:iFPe )
        
        call Diff_ghostNode ( Hx, Stencil, du_SP )
        
        RHS(iSP) = -du_SP
    End do
End Subroutine


subroutine NLI_ghostNode ( Stencil, primL, primR )
use CaseSetup, only: OInt
implicit none
integer,parameter:: nVar = 1
real(kind=8),INTENT(IN ):: Stencil(OInt+1)
real(kind=8),INTENT(OUT):: primL, primR
real(kind=8):: StencilL(nVar,OInt)
real(kind=8):: StencilR(nVar,OInt)
    !----------------------------------------
    StencilL(1,:) = Stencil(1:OInt)
    StencilR(1,:) = 0.0
    primR = 0.0

    IF(     OInt == 3 )Then
        call WCNS_MR_O3 ( nVar, StencilL, primL )
    ElseIF( OInt == 5 )Then
        call WCNS_MR_O5 ( nVar, StencilL, primL )
    ElseIF( OInt == 7 )Then
        call WCNS_MR_O7 ( nVar, StencilL, primL )
    ElseIF( OInt == 9 )Then
        call WCNS_MR_O9 ( nVar, StencilL, primL )
    ElseIF( OInt == 11 )Then
        call WCNS_MR_O11( nVar, StencilL, primL )
    EndIF
End subroutine



subroutine NLI_ghostNode_my_AS ( Stencil, primL, primR )
use CaseSetup, only: OInt
implicit none
integer,parameter:: nVar = 1
real(kind=8),INTENT(IN ):: Stencil(OInt+1)
real(kind=8),INTENT(OUT):: primL, primR
real(kind=8):: StencilL(nVar,OInt)
real(kind=8):: StencilR(nVar,OInt)
    !----------------------------------------
    StencilL(1,:) = Stencil(1:OInt)
    StencilR(1,:) = 0.0
    primR = 0.0

    IF(     OInt == 3 )Then
        call TCNS_ASF102_O3 ( nVar, StencilL, primL )
    ElseIF( OInt == 5 )Then
        call TCNS_ASF102_O5 ( nVar, StencilL, primL )
    ElseIF( OInt == 7 )Then
        call TCNS_ASF102_O7 ( nVar, StencilL, primL )
    ElseIF( OInt == 9 )Then
        call TCNS_ASF102_O9 ( nVar, StencilL, primL )
    ElseIF( OInt == 11 )Then
        call TCNS_ASF102_O11( nVar, StencilL, primL )
    EndIF
End subroutine





subroutine Diff_ghostNode ( Hx, Stencil, du_SP )
use CaseSetup, only: ODiff
implicit none
real(kind=8),INTENT(IN ):: Hx
real(kind=8),INTENT(IN ):: Stencil(ODiff)
real(kind=8),INTENT(OUT):: du_SP
    IF(     ODiff == 4 )Then
        call Diff_E4 ( Hx, Stencil, du_SP )
    ElseIF( ODiff == 6 )Then
        call Diff_E6 ( Hx, Stencil, du_SP )
    ElseIF( ODiff == 8 )Then
        call Diff_E8 ( Hx, Stencil, du_SP )
    ElseIF( ODiff == 10 )Then
        call Diff_E10( Hx, Stencil, du_SP )
    ElseIF( ODiff == 12 )Then
        call Diff_E12( Hx, Stencil, du_SP )
    EndIF
End subroutine