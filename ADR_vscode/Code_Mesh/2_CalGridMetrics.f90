!subroutine CalGridMetrics
!use Variables, only: Nx,Ny, useSCMM
!use Coordinates
!use Coefficients, only: QWeights_X,QWeights_Y
!implicit none
!integer i,j
!real(kind=8):: Row_SP_X(Nx),    Row_SP_Y(Nx),     Xxi_j(Nx),  Yxi_j(Nx)
!real(kind=8):: Row_FP_X(Nx+1),  Row_FP_Y(Nx+1)
!real(kind=8):: Col_SP_X(Ny),    Col_SP_Y(Ny),     Xeta_i(Ny), Yeta_i(Ny)
!real(kind=8):: Col_FP_X(Ny+1),  Col_FP_Y(Ny+1)
!real(kind=8):: XYeta_xi(Nx,Ny),  XetaY_xi(Nx,Ny),  XYxi_eta(Nx,Ny),  XxiY_eta(Nx,Ny) !> 用于SCMM
!real(kind=8):: XYeta_xi_j(Nx),   XetaY_xi_j(Nx),   XYxi_eta_i(Ny),   XxiY_eta_i(Ny)  !> 用于SCMM
!real(kind=8):: Area
!
!    !==================================!
!    !           SP度量                !
!    !==================================!
!    Allocate( SP_Xxi(Nx,Ny),SP_Yxi(Nx,Ny), SP_Xeta(Nx,Ny),SP_Yeta(Nx,Ny) ,SP_Jij(Nx,Ny) )  
!
!    Do j=1,Ny
!        Row_FP_X(1:Nx+1) = FPx_X(1:Nx+1,j)
!        Row_FP_Y(1:Nx+1) = FPx_Y(1:Nx+1,j)
!        call Difference(1,Nx,hx, Row_FP_X, Xxi_j)
!        call Difference(1,Nx,hx, Row_FP_Y, Yxi_j)
!        SP_Xxi(1:Nx,j) = Xxi_j(1:Nx)
!        SP_Yxi(1:Nx,j) = Yxi_j(1:Nx)
!    End Do
!
!    Do i=1,Nx
!        Col_FP_X(1:Ny+1) = FPy_X(i,1:Ny+1)
!        Col_FP_Y(1:Ny+1) = FPy_Y(i,1:Ny+1)
!        call Difference(1,Ny,hy, Col_FP_X, Xeta_i)
!        call Difference(1,Ny,hy, Col_FP_Y, Yeta_i)
!        SP_Xeta(i,1:Ny) = Xeta_i(1:Ny)
!        SP_Yeta(i,1:Ny) = Yeta_i(1:Ny)
!    End Do
!    
!
!    !==================================!
!    !          FPx度量                !
!    !==================================!
!    Allocate( FPx_Xeta(Nx+1,Ny), FPx_Yeta(Nx+1,Ny) )
!    Allocate( FPx_JDMx(Nx+1,Ny), FPx_JDMy(Nx+1,Ny) )
!
!    Do i=1,Nx+1
!        Col_FP_X(1:Ny+1) = FP0_X(i,1:Ny+1)
!        Col_FP_Y(1:Ny+1) = FP0_Y(i,1:Ny+1)
!        call Difference(1,Ny,hy, Col_FP_X, Xeta_i)
!        call Difference(1,Ny,hy, Col_FP_Y, Yeta_i)
!        FPx_Xeta(i,1:Ny) = Xeta_i(1:Ny)
!        FPx_Yeta(i,1:Ny) = Yeta_i(1:Ny)
!        
!        FPx_JDMx(i,1:Ny) =  Yeta_i(1:Ny)
!        FPx_JDMy(i,1:Ny) = -Xeta_i(1:Ny)
!    End Do
!    
!
!    !==================================!
!    !          FPy度量                !
!    !==================================!
!    allocate(  FPy_Xxi(Nx,Ny+1),  FPy_Yxi(Nx,Ny+1) )
!    Allocate( FPy_JDMx(Nx,Ny+1), FPy_JDMy(Nx,Ny+1) )
!
!    Do j=1,Ny+1
!        Row_FP_X(1:Nx+1) = FP0_X(1:Nx+1,j)
!        Row_FP_Y(1:Nx+1) = FP0_Y(1:Nx+1,j)
!        call Difference(1,Nx,hx, Row_FP_X, Xxi_j)
!        call Difference(1,Nx,hx, Row_FP_Y, Yxi_j)
!        FPy_Xxi(1:Nx,j) = Xxi_j(1:Nx)
!        FPy_Yxi(1:Nx,j) = Yxi_j(1:Nx)
!        
!        FPy_JDMx(1:Nx,j) = -Yxi_j(1:Nx)
!        FPy_JDMy(1:Nx,j) =  Xxi_j(1:Nx)
!    End Do
!
!
!
!    !==================================!
!    !             雅可比矩阵            !
!    !==================================!
!    IF( useSCMM )Then
!        write(*,*) "SCMM开启."
!        Do j=1,Ny
!        Do i=1,Nx+1
!            Row_FP_X(i) = FPx_X(i,j)*FPx_Yeta(i,j)
!            Row_FP_Y(i) = FPx_Y(i,j)*FPx_Xeta(i,j)
!        End Do
!        call Difference(1,Nx,hx, Row_FP_X, XYeta_xi_j)
!        call Difference(1,Nx,hx, Row_FP_Y, XetaY_xi_j)
!        Do i=1,Nx
!            XYeta_xi(i,j) = XYeta_xi_j(i)
!            XetaY_xi(i,j) = XetaY_xi_j(i)
!        End Do
!        End Do
!    
!        Do i=1,Nx
!        Do j=1,Ny+1
!            Col_FP_X(j) = FPy_X(i,j)*FPy_Yxi(i,j)
!            Col_FP_Y(j) = FPy_Y(i,j)*FPy_Xxi(i,j)
!        End Do
!        call Difference(1,Ny,hy, Col_FP_X, XYxi_eta_i)
!        call Difference(1,Ny,hy, Col_FP_Y, XxiY_eta_i)
!        Do j=1,Ny
!            XYxi_eta(i,j) = XYxi_eta_i(j)
!            XxiY_eta(i,j) = XxiY_eta_i(j)
!        End Do
!        End Do
!
!        Do j=1,Ny
!        Do i=1,Nx
!            SP_Jij(i,j) = 0.5* ( XYeta_xi(i,j)  -  XYxi_eta(i,j)  +  XxiY_eta(i,j)  -  XetaY_xi(i,j) )
!        End Do
!        End Do
!    Else    
!        write(*,*) "SCMM关闭."
!        Do j=1,Ny
!        Do i=1,Nx
!            SP_Jij(i,j) = SP_Xxi(i,j)*SP_Yeta(i,j) - SP_Yxi(i,j)*SP_Xeta(i,j)
!        End Do
!        End Do 
!    End IF
!    
!    
!    !------------------------------------
!    Area = 0.0
!    Do j=1,Ny
!    Do i=1,Nx
!        Area = Area + abs(SP_Jij(i,j)) *hx*QWeights_X(i) *hy*QWeights_Y(j)
!    End Do 
!    End Do
!    write(*,*) "物理面积 = ", Area
!End Subroutine