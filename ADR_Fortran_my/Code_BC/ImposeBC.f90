!&===================================================&!
!&         传统SD边界方法            &!
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
    !======> 远场
    IF( BC_Type == FarField .or. BC_Type == Inflow .or. BC_Type == Outflow )Then
        call Far_Field( qBnd(1),qBnd(2),qBnd(3),qBnd(4) )
      
    !======> 滑移壁面
    ElseIF( BC_Type == SlipWall )Then
        call BoundaryFluxPointDM (BC_Side, lineID, DM)        !> 获取边界FP法向量
        call BoundaryStencilPrim (BC_Side, lineID, Stencil)   !> 获取下风向SP插值模板
        Do iVar=1,nEq
            q_IN(iVar) = sum( CoeInt_BNR(Order-1,:) * Stencil(iVar,:) )   !> 插值得到'q_IN'
        End Do
        call BC_SlipWall_vecN(DM,q_IN, qBnd)    !> 计算'qBnd'    
      
    !======> 周期性
    ElseIF( BC_Type == Periodic )Then
        call OppositeSide (BC_Side, oppSide)
        call BoundaryStencilPrim (oppSide, lineID, Stencil)   !> 获取下风向SP插值模板
        Do iVar=1,nEq
            qBnd(iVar) = sum( CoeInt_BNR(Order-1,:) * Stencil(iVar,:) )   !> 插值得到'qBnd'
        End Do
      
    !======> 未包含
    Else
        write (*,*) '传统SD边界方法中未包含的BC_Type!!!  BC_Type = ', BC_Type
        stop
    EndIF
End Subroutine
    

!> 取边界插值模板值：取边界附近的Order个SP值
!&===================================================&!
Subroutine BoundaryStencilPrim (BC_Side, lineID, Stencil)
use Variables,only: nEq,Nx,Ny, Order, qPrim, WestSide,EastSide,SouthSide,NorthSide
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: Stencil(nEq,Order)
    !> 西边: 模板SP={ 1, 2, ..., P }
    IF( BC_Side == WestSide )Then
        stencil(:, 1:Order)  =  qPrim(:, 1:Order,       lineID)
        call FlipArray(nEq,Order, stencil)
    !> 东边: 模板SP={ Nx-P+1, Nx-P+2, ..., Nx }
    ElseIF( BC_Side == EastSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, Nx-Order+1:Nx, lineID)
        
    !> 南边: 模板SP={ 1, 2, ..., P }
    ElseIF( BC_Side == SouthSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, lineID, 1:Order)
        call FlipArray(nEq,Order, stencil)
    !> 北边: 模板SP={ Ny-P+1, Ny-P+2, ..., Ny }
    ElseIF( BC_Side == NorthSide )Then   
        stencil(:, 1:Order)  =  qPrim(:, lineID, Ny-Order+1:Ny)
    EndIF
End Subroutine
    
    
!> 取边界FPx/FPy的法向量
!&===================================================&!
Subroutine BoundaryFluxPointDM (BC_Side, lineID, DM)
use Variables,only: nEq,Nx,Ny, WestSide,EastSide,SouthSide,NorthSide
use Coordinates,only: FPx_JDMx,FPx_JDMy, FPy_JDMx,FPy_JDMy
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: DM(2)
    !> 西边: 左边界FPx
    IF( BC_Side == WestSide )Then
        DM(1)  =  FPx_JDMx( 1, lineID )
        DM(2)  =  FPx_JDMy( 1, lineID )
    !> 东边: 右边界FPx
    ElseIF( BC_Side == EastSide )Then   
        DM(1)  =  FPx_JDMx( Nx+1, lineID )
        DM(2)  =  FPx_JDMy( Nx+1, lineID )
        
    !> 南边: 下边界FPy
    ElseIF( BC_Side == SouthSide )Then   
        DM(1)  =  FPy_JDMx( lineID, 1 )
        DM(2)  =  FPy_JDMy( lineID, 1 )
    !> 北边: 上边界FPy
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