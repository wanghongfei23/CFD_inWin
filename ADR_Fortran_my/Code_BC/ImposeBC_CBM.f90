!&===================================================&!
!&         特征边界方法            &!
!&===================================================&!  
subroutine SetBoundaryValueByType_CBM (BC_Side, BC_Type, lineID, qBnd_CBM)
use Mat_Manipulations
use Variables
use Coordinates
use Coefficients,only: CoeInt_BNR
implicit none
integer,INTENT(IN):: BC_Side, BC_Type, lineID
real,INTENT(OUT):: qBnd_CBM(nEq)
real:: qPrim_Inf(nEq) !> 远场qPrim在X-Y坐标系中
real:: qPrim_Stencil(nEq,Order) !> qPrim在X-Y坐标系中
real:: hatW_SP(nEq,Order), hatW_Bnd(nEq), hatW_Inf(nEq)   !> hatW: qPrim在N-T坐标系中
real:: hatQ_SP(nEq,Order), hatQ_Bnd(nEq), hatQ_Inf(nEq)   !> hatQ: qCons在N-T坐标系中
real:: hatV_SP(nEq,Order), hatV_Int(nEq), hatV_Inf(nEq)   !> hatV: qChar在N-T坐标系中
integer:: eIDX,iSP,iVar
real:: DM(2), mod, VNx,VNy, VTx,VTy
real:: LEV(4,4)
real:: matRot(2,2), matIRot(2,2)
!---------------------------------------------------!
        
    
    !======> 情况：远场
    IF( BC_Type == FarField )Then
        call Far_Field (qBnd_CBM(1), qBnd_CBM(2), qBnd_CBM(3), qBnd_CBM(4))
        return
    !======> 情况：周期性
    ElseIF( BC_Type == Periodic )Then
        call CBM_Periodic (BC_Side, lineID, qBnd_CBM)
        return
    EndIF
     

    
    !==============================================================================================================!
    !======> 获取插值模板
    call BoundaryFluxPointDM (BC_Side, lineID, DM)              !> 获取边界FP法向量
    call BoundaryStencilPrim (BC_Side, lineID, qPrim_Stencil)   !> 获取下风向SP插值模板
        
    !======> 获取边界单位内法向量和切向量 --> 旋转矩阵
    IF( BC_Side == SouthSide .or. BC_Side == WestSide )Then
        eIDX = 1
    ElseIF( BC_Side == NorthSide .or. BC_Side == EastSide )Then
        eIDX = -1
    End IF
    mod = sqrt(DM(1)**2 + DM(2)**2)
    VNx = eIDX * DM(1)/mod   !> 旋转变换矩阵的单位向量
    VNy = eIDX * DM(2)/mod
    VTx = -VNy
    VTy =  VNx
    matRot(1,1:2) = (/ VNx, VNy /)  !>  Rot = [ cos sin]
    matRot(2,1:2) = (/ VTx, VTy /)  !         [-sin cos]
    matIRot(1,1:2) = (/  VNx, -VNy /)  !> IRot = [cos -sin]
    matIRot(2,1:2) = (/ -VTx,  VTy /)  !         [sin  cos]
    
    !======> 旋转SP原始变量到N-T坐标系，然后转换到特征变量系
    Do iSP=1,Order
        hatW_SP(:,iSP) = CoorSysRotation( qPrim_Stencil(:,iSP), matRot )    !> qPrim_{iSP} --> hatW_{iSP}
        call PrimToCons( hatW_SP(:,iSP), hatQ_SP(:,iSP) )                   !>  hatW_{iSP} --> hatQ_{iSP}
    End Do
    
    !======> 计算N-T坐标系下Nearest_SP处的左特征向量矩阵
    call GetLeftEigenVectors_Cons_ortho( hatW_SP(:,Order), LEV )    !> hatW_{Nearest_SP} --> LEV
            
    !======> 计算N-T坐标系下的所有SP特征变量值
    Do iSP=1,Order
        hatV_SP(:,iSP) =  MatMulVec (LEV, hatQ_SP(:,iSP), nEq, nEq)    !> hatV_{iSP} = LEV * hatQ_{iSP}
    End Do
        
    !======> 从各SP特征变量值插值得到边界(Interpolated)特征变量值
    Do iVar=1,4
        hatV_Int(iVar) = sum( CoeInt_BNR(Order-1,:) * hatV_SP(iVar,:) )    !> hatV_{iSP} --> hatV_{Bnd}
    End Do
    
    !======> 取X-Y坐标系下FarField原始变量值，然后转换到N-T坐标系下的FarField特征变量值
    call Far_Field (qPrim_Inf(1), qPrim_Inf(2), qPrim_Inf(3), qPrim_Inf(4))   !> qPrim_Infty
    hatW_Inf =  CoorSysRotation( qPrim_Inf, matRot )                          !> qPrim_Infty --> hatW_Infty
    call PrimToCons (hatW_Inf, hatQ_Inf)                                      !>  hatW_Infty --> hatQ_Infty
    hatV_Inf =  MatMulVec (LEV, hatQ_Inf, nEq, nEq)                           !>  hatQ_Infty --> hatV_Infty
    !==============================================================================================================!
    
    
    !======> 情况：无滑移壁面
    IF( BC_Type == SlipWall )Then
        call CBM_SlipWall (LEV, hatV_Int, hatW_Bnd)
     
    !======> 情况：流入边界
    ElseIF( BC_Type == Inflow )Then
        IF( BC_Side /= WestSide )Then
            write(*,*) '流出边界不在西边! BC_Side =',BC_Side;    stop
        End IF
        call CBM_Inflow (LEV, hatV_Int, hatV_Inf, hatW_Bnd)
    !======> 情况：流出边界
    ElseIF( BC_Type == Outflow )Then
        IF( BC_Side /= EastSide )Then
            write(*,*) '流出边界不在东边! BC_Side =',BC_Side;    stop
        End IF
        call CBM_Outflow (LEV, hatV_Int, hatV_Inf, hatW_Bnd)
     
    !======> 情况：未包含
    Else
        write (*,*) 'CBM中未包含的BC_Type!!!  BC_Type = ', BC_Type
        stop
    EndIF
        
    
    !======> 反旋转到原始变量的X-Y坐标系
    qBnd_CBM = CoorSysRotation( hatW_Bnd, matIRot )   !> hatW_{Bnd} --> qPrim_{Bnd}
    
End subroutine
    
    
    
!&===================================================&!
!&           计算边界处的CBM通量            &!
!&===================================================&! 
Subroutine CBM_BoundaryFlux (BC_Side, lineID, Flux)
use Variables,only: nEq, WBC,EBC,SBC,NBC, WestSide,EastSide,SouthSide,NorthSide, FarField,SlipWall,Inflow,Outflow
use Coordinates,only: WB_Type,EB_Type,SB_Type,NB_Type
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: Flux(nEq)
real:: DM(2), qPrim_Bnd(nEq)
integer:: BC_Type

    IF( BC_Type == FarField )Then
        return
    EndIF
    
    !> 西边
    IF( BC_Side == WestSide )Then
        call BoundaryFluxPointDM (BC_Side, lineID, DM)
        qPrim_Bnd = WBC(:,lineID)
        BC_Type = WB_Type(lineID)
    !> 东边
    ElseIF( BC_Side == EastSide )Then   
        call BoundaryFluxPointDM (BC_Side, lineID, DM)
        qPrim_Bnd = EBC(:,lineID)
        BC_Type = EB_Type(lineID)
    !> 南边
    ElseIF( BC_Side == SouthSide )Then   
        call BoundaryFluxPointDM (BC_Side, lineID, DM)
        qPrim_Bnd = SBC(:,lineID)
        BC_Type = SB_Type(lineID)
    !> 北边
    ElseIF( BC_Side == NorthSide )Then   
        call BoundaryFluxPointDM (BC_Side, lineID, DM)
        qPrim_Bnd = NBC(:,lineID)
        BC_Type = NB_Type(lineID)
    EndIF
    
    !===> 转换为通量
    IF( BC_Type == SlipWall .or. BC_Type == Inflow .or. BC_Type == Outflow )Then
        call PrimToFlux(qPrim_Bnd, DM, Flux)
    EndIF
End Subroutine
    


!&===================================================&!!&===================================================&!
Subroutine CBM_Inflow (LEV, hatV_Int, hatV_Inf, hatW_Bnd)
use Mat_Manipulations
use Variables,only: nEq
implicit none
real,INTENT(IN ):: LEV(nEq,nEq)
real,INTENT(IN ):: hatV_Int(nEq)    !> hatV: qChar在N-T坐标系中
real,INTENT(IN ):: hatV_Inf(nEq)    !> hatV: qChar在N-T坐标系中
real,INTENT(OUT):: hatW_Bnd(nEq)    !> hatW: qPrim在N-T坐标系中
real:: hatQ_Bnd(nEq)                !> hatQ: qCons在N-T坐标系中
real:: matL44(nEq,nEq), vecR4(nEq), matIL44(nEq,nEq)
    !======> 组装LHS系数矩阵 & RHS向量
    matL44 = LEV
    vecR4(1) = hatV_Int(1)  !> 亚音速:      qn-c<0 --> q_Int
    vecR4(2) = hatV_Inf(2)  !> 亚音速流入: qn>0 --> q_Infty
    vecR4(3) = hatV_Inf(3)  !> 亚音速流入: qn>0 --> q_Infty
    vecR4(4) = hatV_Inf(4)  !> 亚音速:      qn+c>0 --> q_Infty
    !======> 求解线性方程组得到特征变量
    call Mat_Inversion (matL44, 4, matIL44)
    hatQ_Bnd = MatMulVec (matIL44, vecR4, 4,4)
    !======> 转换得到原始变量
    call ConsToPrim( hatQ_Bnd, hatW_Bnd )   !> hatQ_{Bnd} --> hatW_{Bnd}
End Subroutine
   
!&===================================================&!!&===================================================&!  
Subroutine CBM_Outflow (LEV, hatV_Int, hatV_Inf, hatW_Bnd)
use Mat_Manipulations
use Variables,only: nEq
implicit none
real,INTENT(IN ):: LEV(nEq,nEq)
real,INTENT(IN ):: hatV_Int(nEq)    !> hatV: qChar在N-T坐标系中
real,INTENT(IN ):: hatV_Inf(nEq)    !> hatV: qChar在N-T坐标系中
real,INTENT(OUT):: hatW_Bnd(nEq)    !> hatW: qPrim在N-T坐标系中
real:: hatQ_Bnd(nEq)                !> hatQ: qCons在N-T坐标系中
real:: matL44(nEq,nEq), vecR4(nEq), matIL44(nEq,nEq)
    !======> 组装LHS系数矩阵 & RHS向量
    matL44 = LEV
    vecR4(1) = hatV_Int(1)  !> 亚音速:       qn-c<0 --> q_Int
    vecR4(2) = hatV_Int(2)  !> 亚音速流出: qn<0 --> q_Int
    vecR4(3) = hatV_Int(3)  !> 亚音速流出: qn<0 --> q_Int
    vecR4(4) = hatV_Inf(4)  !> 亚音速:       qn+c>0 --> q_Infty
    !======> 求解线性方程组得到特征变量
    call Mat_Inversion (matL44, 4, matIL44)
    hatQ_Bnd = MatMulVec (matIL44, vecR4, 4,4)
    !======> 转换得到原始变量
    call ConsToPrim( hatQ_Bnd, hatW_Bnd )   !> hatQ_{Bnd} --> hatW_{Bnd}
End Subroutine
    
    
    
    
!&===================================================&!!&===================================================&!
Subroutine CBM_SlipWall (LEV, hatV_Int, hatW_Bnd)
use Mat_Manipulations
use Variables,only: nEq
implicit none
real,INTENT(IN ):: LEV(nEq,nEq)
real,INTENT(IN ):: hatV_Int(nEq)    !> hatV: qChar在N-T坐标系中
real,INTENT(OUT):: hatW_Bnd(nEq)    !> hatW: qPrim在N-T坐标系中
real:: hatQ_Bnd(nEq)                !> hatQ: qCons在N-T坐标系中
real:: matL33(3,3), vecR3(3), vecL3(3), matIL33(3,3)
integer:: iRow
    !======> 组装LHS系数矩阵 & RHS向量
    Do iRow=1,3
        matL33(iRow,:) = (/ LEV(iRow,1), LEV(iRow,3), LEV(iRow,4) /)
        vecR3 (iRow)   = hatV_Int(iRow)
    End Do
    !======> 求解线性方程组得到特征变量
    matIL33 = GetMatInversion_N33( matL33 )    !> 直接按系数3*3矩阵求逆
    vecL3 =  MatMulVec (matIL33, vecR3, 3,3)
    !======> 转换得到原始变量
    hatQ_Bnd(1) = vecL3(1)
    hatQ_Bnd(2) =   0.0   !> 无穿透条件
    hatQ_Bnd(3) = vecL3(2)
    hatQ_Bnd(4) = vecL3(3)
    call ConsToPrim( hatQ_Bnd, hatW_Bnd )   !> hatQ_{Bnd} --> hatW_{Bnd}
    hatW_Bnd(2) =   0.0   !> 无穿透条件
End Subroutine
    
    
    
    
!&===================================================&!!&===================================================&!
Subroutine CBM_Periodic (BC_Side, lineID, qBnd)
use Variables,only: nEq,Order
use Coefficients,only: CoeInt_BNR
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(OUT):: qBnd(nEq)
integer:: oppSide, iVar
real:: Stencil(nEq,Order) !> qPrim在X-Y坐标系中

    call OppositeSide (BC_Side, oppSide)
    call BoundaryStencilPrim (oppSide, lineID, Stencil)   !> 获取下风向SP插值模板
    
    Do iVar=1,nEq
        qBnd(iVar) = sum( CoeInt_BNR(Order-1,:) * Stencil(iVar,:) )   !> 插值得到'qBnd'
    End Do
End Subroutine










!&===================================================&!!&===================================================&!!&===================================================&!
subroutine Check_CBM(BC_Side, lineID, qBnd)
USE lapack95
use Mat_Manipulations
use Variables
use Coordinates
use Coefficients,only: CoeInt_BNR
implicit none
integer,INTENT(IN):: BC_Side, lineID
real,INTENT(IN):: qBnd(4)

integer:: i,j, eIDX,iSP,iVar,iRow, flag
integer:: info, ipiv2(2), ipiv3(3), ipiv4(40)    !> For LAPACK Functions
real:: DM(2)
real:: mod, VNx,VNy, VTx,VTy, matRot(2,2), matIRot(2,2)
real:: LEV(4,4)
real:: REV(4,4)
real:: qPrim_Stencil(nEq,Order) !> qPrim在X-Y坐标系中
real:: qBnd_CBM(nEq)
real:: hatW_SP(nEq,Order), hatW_Bnd(nEq)   !> hat_W: qPrim在N-T坐标系中
real:: hatQ_SP(nEq,Order), hatQ_Bnd(nEq)   !> hat_Q: qCons在N-T坐标系中
real:: hatV_SP(nEq,Order), hatV_Int(nEq)   !> hat_V: qChar在N-T坐标系中

real:: wr3(3), wi3(3), vr3(3,3), vl3(3,3)
real:: wr4(4), wi4(4), vr4(4,4), vl4(4,4)
real:: mat33L(3,3), mat33L_I(3,3), vec3R(3), vec3L(3)
real:: mat44L(4,4), mat44L_I(4,4), vec4R(4)
real:: mat4E(4,4)


    call BoundaryFluxPointDM (BC_Side, lineID, DM)              !> 获取边界FP法向量
    call BoundaryStencilPrim (BC_Side, lineID, qPrim_Stencil)   !> 获取下风向SP插值模板
        
        
    !=================================================> 计算单位法向量和切向量
    IF( BC_Side == SouthSide .or. BC_Side == WestSide )Then
        eIDX = 1
    ElseIF( BC_Side == NorthSide .or. BC_Side == EastSide )Then
        eIDX = -1
    End IF
    mod = sqrt(DM(1)**2 + DM(2)**2)
    VNx = eIDX * DM(1)/mod
    VNy = eIDX * DM(2)/mod
    VTx = -VNy
    VTy =  VNx
    matRot(1,1:2) = (/ VNx, VNy /)  !>  Rot = [ cos sin]  旋转矩阵
    matRot(2,1:2) = (/ VTx, VTy /)  !         [-sin cos]
    matIRot(1,1:2) = (/  VNx, -VNy /)  !> IRot = [cos -sin]  反旋转矩阵
    matIRot(2,1:2) = (/ -VTx,  VTy /)  !         [sin  cos]
    
                
    !=================================================> 旋转原始变量到N-T坐标系，然后转换到特征变量系
    Do iSP=1,Order
        hatW_SP(:,iSP) = CoorSysRotation( qPrim_Stencil(:,iSP), matRot )    !> qPrim_{iSP} --> hatW_{iSP}
        call PrimToCons( hatW_SP(:,iSP), hatQ_SP(:,iSP) )    !> hatW_{iSP} --> hatQ_{iSP}
    End Do
    
    !=================================================> 计算N-T坐标系下特征结构 hatW_{Nearest_SP} --> LEV
    DM = (/ 1.0, 0.0 /)
    call GetLeftEigenVectors_Cons( hatW_SP(:,Order), DM, LEV )
    !call GetLeftEigenVectors_Cons_ortho( hatW_SP(:,Order), LEV )
    
            mat44L = 0.0;
            mat44L(2,2) = qPrim_Stencil(2,Order)*VNx + qPrim_Stencil(3,Order)*VNy;
            mat44L(3,3) = mat44L(2,2);
            mat44L(1,1) = mat44L(2,2) - sqrt(1.4*qPrim_Stencil(4,Order)/qPrim_Stencil(1,Order));
            mat44L(4,4) = mat44L(2,2) + sqrt(1.4*qPrim_Stencil(4,Order)/qPrim_Stencil(1,Order));
            write(*,*) "---> 特征结构_An = ";    write(*,*) (mat44L(i,:),i=1,4);     write(*,*)''
            
            !call GetRightEigenVectors_Cons( hatW_SP(:,Order), DM, REV )
            call GetRightEigenVectors_Cons_ortho( hatW_SP(:,Order), REV )
            write(*,*) "---> LEV = ";    write(*,*) (LEV(i,:),i=1,4);     write(*,*)''
            write(*,*) "---> REV = ";    write(*,*) (REV(i,:),i=1,4);     write(*,*)''
            flag = Mat_CheckInverse( LEV, REV, 4 )
            
            mat44L_I = MATMUL(REV,MATMUL(mat44L,LEV))
            wr4 = GetMatEigenvalues_MKL(mat44L_I, 4)
            write(*,*) "特征值 LEV *Lambda_n *REV = ";    write(*,*) wr4;     write(*,*)''
        
    !=================================================> 计算N-T坐标系下的特征变量 hatV_{iSP} = LEV * hatQ_{iSP}
    Do iSP=1,Order
        !hatV_SP(:,iSP) =  MatMulVec( LEV, hatQ_SP(:,iSP), 4,4 )
        Do iVar=1,3
            hatV_SP(iVar,iSP) =  sum( LEV(iVar,:) * hatQ_SP(:,iSP) )
        End Do
    End Do
            
         
    !=================================================> 插值得到边界特征变量 hatV_{Bnd} 前3个分量 v1_{Bnd},v2_{Bnd},v3_{Bnd}
    Do iVar=1,3
        hatV_Int(iVar) = sum( CoeInt_BNR(Order-1,:) * hatV_SP(iVar,:) )
    End Do
    
    !=================================================> 组装LHS系数矩阵 & RHS向量，然后求解特征变量 hatQ_{Bnd}
    Do iRow=1,3
        mat33L(iRow,:) = (/ LEV(iRow,1), LEV(iRow,3), LEV(iRow,4) /)
        vec3R (iRow)   = hatV_Int(iRow)
    End Do
    
    mat33L_I = GetMatInversion_N33( mat33L )    !> 3*3矩阵求逆
    vec3L =  MatMulVec( mat33L_I, vec3R, 3,3 )
    !vec3L = SolveLinearSystem_LS_MKL( mat33L, vec3R, 3 )   !> Solve Ax=B --> x
    
    hatQ_Bnd(1) = vec3L(1)
    hatQ_Bnd(2) = 0.0   !> 强制无穿透条件
    hatQ_Bnd(3) = vec3L(2)
    hatQ_Bnd(4) = vec3L(3)
    call ConsToPrim( hatQ_Bnd, hatW_Bnd )
    
    
    !=================================================> 反旋转到X-Y坐标系得到原始变量 hatW_{Bnd} --> qPrim_{Bnd}
    hatW_Bnd(2) = 0.0   !> 强制无穿透条件
    qBnd_CBM = CoorSysRotation( hatW_Bnd, matIRot )
        
        
        write(*,*) "CBM qBnd_1 = ",qBnd_CBM(1), "   ref-CBM = ",abs( qBnd(1) - qBnd_CBM(1) )
        write(*,*) "CBM qBnd_2 = ",qBnd_CBM(2), "   ref-CBM = ",abs( qBnd(2) - qBnd_CBM(2) )
        write(*,*) "CBM qBnd_3 = ",qBnd_CBM(3), "   ref-CBM = ",abs( qBnd(3) - qBnd_CBM(3) )
        write(*,*) "CBM qBnd_4 = ",qBnd_CBM(4), "   ref-CBM = ",abs( qBnd(4) - qBnd_CBM(4) )
        write(*,*) qBnd_CBM(4)/qBnd_CBM(1)**(1.4) - 1.0
        write(*,*) "========================================================================"
    
    !stop

End subroutine