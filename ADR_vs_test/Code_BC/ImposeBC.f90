
!&===================================================&!
!&         Traditional SD Boundary Method            &!
!&===================================================&!
Subroutine SetBoundaryValueByType (BC_Side, BC_Type, lineID, qBnd)
use Variables
use Coordinates
use Coefficients,only: CoeInt_BNR
implicit none
integer,INTENT(IN):: BC_Side, BC_Type, lineID
real,INTENT(OUT):: qBnd(4)
integer:: iVar,oppSide
real:: stencil(nEq,Order)
real:: q_IN(nEq), DM(2)
    !======> Far-Field
    IF( BC_Type == FarField .or. BC_Type == Inflow .or. BC_Type == Outflow )Then
        call Far_Field( qBnd(1),qBnd(2),qBnd(3),qBnd(4) )
      
    !======> Slip Wall
    ElseIF( BC_Type == SlipWall )Then
        call BoundaryFluxPointDM (BC_Side, lineID, DM)        !> Get Boundary FP Normal
        call BoundaryStencilPrim (BC_Side, lineID, Stencil)   !> Get Downwind SP Interpolation Stencil
        Do iVar=1,nEq
            q_IN(iVar) = sum( CoeInt_BNR(Order-1,:) * Stencil(iVar,:) )   !> Interpolate for 'q_IN'
        End Do
        call BC_SlipWall_vecN(DM,q_IN, qBnd)    !> Compute 'qBnd'    
      
    !======> Periodic
    ElseIF( BC_Type == Periodic )Then
        call OppositeSide (BC_Side, oppSide)
        call BoundaryStencilPrim (oppSide, lineID, Stencil)   !> Get Downwind SP Interpolation Stencil
        Do iVar=1,nEq
            qBnd(iVar) = sum( CoeInt_BNR(Order-1,:) * Stencil(iVar,:) )   !> Interpolate for 'qBnd'
        End Do
      
    !======> NOT Included
    Else
        write (*,*) 'BC_Type NOT Included in Traditional SD Boundary Method Yet!!!  BC_Type = ', BC_Type
        stop
    EndIF
End Subroutine
    

!> 取边界逆风插值模板：距边界最近的Order个SP值
!&===================================================&!
Subroutine BoundaryStencilPrim (BC_Side, lineID, Stencil)
use Variables,only: nEq,Nx,Ny, Order, qPrim, WestSide,EastSide,SouthSide,NorthSide
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: Stencil(nEq,Order)
    !> West: stencil SP={ 1, 2, ..., P }
    IF( BC_Side == WestSide )Then
        stencil(:, 1:Order)  =  qPrim(:, 1:Order,       lineID)
        call FlipArray(nEq,Order, stencil)
    !> East: stencil SP={ Nx-P+1, Nx-P+2, ..., Nx }
    ElseIF( BC_Side == EastSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, Nx-Order+1:Nx, lineID)
        
    !> South: stencil SP={ 1, 2, ..., P }
    ElseIF( BC_Side == SouthSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, lineID, 1:Order)
        call FlipArray(nEq,Order, stencil)
    !> North: stencil SP={ Ny-P+1, Ny-P+2, ..., Ny }
    ElseIF( BC_Side == NorthSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, lineID, Ny-Order+1:Ny)
    EndIF
End Subroutine
    
    
!> 取边界FPx/FPy处法向
!&===================================================&!
Subroutine BoundaryFluxPointDM (BC_Side, lineID, DM)
use Variables,only: nEq,Nx,Ny, WestSide,EastSide,SouthSide,NorthSide
use Coordinates,only: FPx_JDMx,FPx_JDMy, FPy_JDMx,FPy_JDMy
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: DM(2)
    !> West: Left Boundary FPx
    IF( BC_Side == WestSide )Then
        DM(1)  =  FPx_JDMx( 1, lineID )
        DM(2)  =  FPx_JDMy( 1, lineID )
    !> East: Right Boundary FPx
    ElseIF( BC_Side == EastSide )Then   
        DM(1)  =  FPx_JDMx( Nx+1, lineID )
        DM(2)  =  FPx_JDMy( Nx+1, lineID )
        
    !> South: Lower Boundary FPy
    ElseIF( BC_Side == SouthSide )Then   
        DM(1)  =  FPy_JDMx( lineID, 1 )
        DM(2)  =  FPy_JDMy( lineID, 1 )
    !> North: Upper Boundary FPy
    ElseIF( BC_Side == NorthSide )Then   
        DM(1)  =  FPy_JDMx( lineID, Ny+1 )
        DM(2)  =  FPy_JDMy( lineID, Ny+1 )
    EndIF
End Subroutine
    
    
!> 取对侧边界（用于周期性边界）
!&===================================================&!
Subroutine OppositeSide (SideIn, SideOut)
use Variables, only: WestSide,EastSide,SouthSide,NorthSide
implicit none
integer,INTENT(IN ):: SideIn
integer,INTENT(OUT):: SideOut
    IF( SideIN == WestSide )Then
        SideOut = EastSide
    ElseIF( SideIN == EastSide )Then   
        SideOut = WestSide
    ElseIF( SideIN == SouthSide )Then   
        SideOut = NorthSide
    ElseIF( SideIN == NorthSide )Then   
        SideOut = SouthSide
    EndIF
End Subroutine