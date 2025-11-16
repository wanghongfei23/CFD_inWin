
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
!Subroutine Cal_beta_S1 ( xiL, xiR, beta1 )
!use CaseSetup, only: Order
!implicit none
!real,INTENT(IN ):: xiL, xiR
!real,INTENT(OUT):: beta1
!real,parameter:: epsilon = 1E-40
!integer,parameter:: coeP = (Order-1)/2
!real:: r01, r11
!real:: sig0, sig1, sigSum
!real:: xi0, xi1
!
!    !beta1 = min( xiL**2 , xiR**2 )    !> New beta_1
!    !return
!    
!    xi0 = xiL**2
!    xi1 = xiR**2
!    
!    IF( xi0 >= xi1 )Then
!        r01 = 1./11
!    Else
!        r01 = 10./11
!    EndIF
!    r11 = 1 - r01
!    
!    sig0 = r01 * ( 1 + abs(xi0 - xi1)**coeP/(xi0 + epsilon) )    !> Tau --> Method 1
!    sig1 = r11 * ( 1 + abs(xi0 - xi1)**coeP/(xi1 + epsilon) )
!    sigSum = sig0 + sig1
!    beta1 = ( (sig0*xiL + sig1*xiR)/sigSum )**2
!End Subroutine
    
    
Subroutine Cal_beta_S1 ( xiL, xiR, beta1 )
implicit none
real,INTENT(IN ):: xiL, xiR
real,INTENT(OUT):: beta1
    beta1 = min( xiL**2 , xiR**2 )    !> New beta_1
    return
End Subroutine
    
    
    
    
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
Subroutine WCNS_MR_O3 ( nVar, var_In, var_Out )
use Parameters
implicit none
integer,parameter:: Order = 3            !> Accuracy
integer,parameter:: nSub = (Order+1)/2   !> Number of sub-stencils
integer,parameter:: nDrv =  Order-1      !> Highest derivative order for SI
integer,INTENT(IN):: nVar
real,INTENT(IN ):: var_In(nVar,Order)    !> Stencil
real,INTENT(OUT):: var_Out(nVar)         !> qL
real,parameter:: epsilon = 1E-10    !> WENO-MR
integer:: iVar, iDrv, iSub, coeP
real:: varS1(nVar)  , qL_S1, pL_S1, beta_S1
real:: varS2(nVar,3), qL_S2, pL_S2, beta_S2
real:: Drv_qS(nDrv,nSub)
real:: Drv_pS(nDrv,nSub)
real:: Tau, alpha(nSub), alphaSum, w_NL(nSub)
real:: xiL(nVar), xiR(nVar)

real:: MR_IntS2(3), MR_SI_S2(3,2)
!-------------------------------------------------
    MR_IntS2(1:3    ) = (/  -1./8 ,   3./4,  3./8  /)   !> S2 Interpolation
    MR_SI_S2(1:3, 1 ) = (/  -1./2 ,    0.0,  1./2  /)   !> S2 SI_u'
    MR_SI_S2(1:3, 2 ) = (/    1.0 ,   -2.0,   1.0  /)   !> S2 SI_u''
!-------------------------------------------------
    
    
    !===========================> Interior FP
    varS1(:)      = var_In(:,  2 )    !> S1 = {    n2 |    }
    varS2(:, 1:3) = var_In(:, 1:3)    !> S2 = { n1 n2 | n3 }
    xiL(:) = varS2(:,2) - varS2(:,1)
    xiR(:) = varS2(:,3) - varS2(:,2)
    
    
    !======================================================>!======================================================>
    Do iVar=1,nVar
        !====================================> poly_Q
        qL_S1 = varS1(iVar)                                !> Sub-stencil Interpolation: S1
        qL_S2 = sum( MR_IntS2(:) * varS2(iVar,:) )    !> Sub-stencil Interpolation: S2
        
        Drv_qS = 0.0
        Do iDrv = 1,2
            Drv_qS(iDrv,S2) = sum( MR_SI_S2(:,iDrv) * varS2(iVar,:) )  !> poly_Q Drv: S2
        End Do
        
        !====================================> poly_P
        pL_S1 = qL_S1                                     
        pL_S2 = 1./r22 * ( qL_S2 - r12*pL_S1 )            
        
        !> poly_P derivatives
        Drv_pS = 0.0
        Drv_pS(:,S2) = 1./r22 * ( Drv_qS(:,S2) )
        
        !> Lambda: scale derivatives
        Do iSub = 1,nSub
        Do iDrv = 1,2*(iSub-1)
            coeP = iDrv-1
            Drv_pS(iDrv,iSub) = lambda**coeP *Drv_pS(iDrv,iSub)**2
        End Do
        End Do
        
        !====================================> Beta & Tau
        !> Smoothness Indicators
        call Cal_beta_S1 ( xiL(iVar), xiR(iVar), beta_S1 )
        beta_S2 = sum( Drv_pS(1:2 ,S2) )
        
        !> Tau
        coeP = 2    !> coeP=1(Zhu) --> 2nd-order
        tau = ( abs(beta_S2-beta_S1) )**coeP
        
        !====================================> Weights & Combination
        alpha(S1) = r12 * ( 1 + ( tau/(epsilon + beta_S1) ) )
        alpha(S2) = r22 * ( 1 + ( tau/(epsilon + beta_S2) ) )
        
        alphaSum = sum(alpha)
        w_NL = alpha/alphaSum
        var_Out(iVar) = w_NL(1)*pL_S1  +  w_NL(2)*pL_S2
    End Do
End Subroutine
        
        
        
        
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
Subroutine WCNS_MR_O5 ( nVar, var_In, var_Out )
use Parameters
implicit none
integer,parameter:: Order = 5            !> Accuracy
integer,parameter:: nSub = (Order+1)/2   !> Number of sub-stencils
integer,parameter:: nDrv =  Order-1      !> Highest derivative order for SI
integer,INTENT(IN ):: nVar
real,INTENT(IN ):: var_In(nVar,Order) !> optimal模板
real,INTENT(OUT):: var_Out(nVar)         !> qL
real,parameter:: epsilon = 1E-10    !> WENO-MR
integer:: iVar, iDrv, iSub, coeP
real:: varS1(nVar)  , qL_S1, pL_S1, beta_S1
real:: varS2(nVar,3), qL_S2, pL_S2, beta_S2
real:: varS3(nVar,5), qL_S3, pL_S3, beta_S3
real:: Drv_qS(nDrv,nSub)
real:: Drv_pS(nDrv,nSub)
real:: Tau, alpha(nSub), alphaSum, w_NL(nSub)
real:: xiL(nVar), xiR(nVar)

real:: MR_IntS2(3), MR_SI_S2(3,2)
real:: MR_IntS3(5), MR_SI_S3(5,4)
!-------------------------------------------------
    MR_IntS2(1:3   ) = (/  -1./8,  3./4,  3./8  /)   !> S2 Interpolation
    MR_SI_S2(1:3, 1) = (/  -1./2,   0.0,  1./2  /)  !> S2 SI_u'
    MR_SI_S2(1:3, 2) = (/    1.0,  -2.0,   1.0  /)  !> S2 SI_u''
                  
    MR_IntS3(1:5   ) = (/  3./128,  -5./32,  45./64,  15./32,  -5./128  /)   !> S3 Interpolation
    MR_SI_S3(1:5, 1) = (/   1./12,  -2./3 ,    0.0 ,   2./3 ,  -1./12   /)   !> S3 SI_u'
    MR_SI_S3(1:5, 2) = (/  -1./12,   4./3 ,  -5./2 ,   4./3 ,  -1./12   /)   !> S3 SI_u''
    MR_SI_S3(1:5, 3) = (/   -1./2,    1.0 ,    0.0 ,   -1.0 ,   1./2    /)   !> S3 SI_u'''
    MR_SI_S3(1:5, 4) = (/     1.0,   -4.0 ,    6.0 ,   -4.0 ,    1.0    /)   !> S3 SI_u''''
!-------------------------------------------------


    !===========================> Interior FP
    varS1(:)      = var_In(:,  3 )    !> S1 = {      n3 |      }
    varS2(:, 1:3) = var_In(:, 2:4)    !> S2 = {   n2 n3 | n4   }   
    varS3(:, 1:5) = var_In(:, 1:5)    !> S3 = {n1 n2 n3 | n4,n5}
    xiL(:) = varS2(:,2) - varS2(:,1)
    xiR(:) = varS2(:,3) - varS2(:,2)
    
    
    !======================================================>!======================================================>
    Do iVar=1,nVar
        !====================================> poly_Q
        qL_S1 = varS1(iVar)                                !> Sub-stencil Interpolation: S1
        qL_S2 = sum( MR_IntS2(:) * varS2(iVar,:) )    !> Sub-stencil Interpolation: S2
        qL_S3 = sum( MR_IntS3(:) * varS3(iVar,:) )    !> Sub-stencil Interpolation: S3
        
        Drv_qS = 0.0
        Do iDrv = 1,2
            Drv_qS(iDrv,S2) = sum( MR_SI_S2(:,iDrv) * varS2(iVar,:) )  !> poly_Q Drv: S2
        End Do
        Do iDrv = 1,4
            Drv_qS(iDrv,S3) = sum( MR_SI_S3(:,iDrv) * varS3(iVar,:) )  !> poly_Q Drv: S3
        End Do
        
        !====================================> poly_P
        pL_S1 = qL_S1                                     
        pL_S2 = 1./r22 * ( qL_S2 - r12*pL_S1 )            
        pL_S3 = 1./r33 * ( qL_S3 - r13*pL_S1 - r23*pL_S2 )
        
        !> poly_P derivatives
        Drv_pS = 0.0
        Drv_pS(:,S2) = 1./r22 * ( Drv_qS(:,S2) )
        Drv_pS(:,S3) = 1./r33 * ( Drv_qS(:,S3) - r23*Drv_pS(:,S2) )
        
        !> Lambda: scale derivatives
        Do iSub = 1,nSub
        Do iDrv = 1,2*(iSub-1)
            coeP = iDrv-1
            Drv_pS(iDrv,iSub) = lambda**coeP *Drv_pS(iDrv,iSub)**2
        End Do
        End Do
        
        !====================================> Beta & Tau
        !> Smoothness Indicators
        call Cal_beta_S1 ( xiL(iVar), xiR(iVar), beta_S1 )
        beta_S2 = sum( Drv_pS(1:2 ,S2) )
        beta_S3 = sum( Drv_pS(1:4 ,S3) )
        
        !> Tau
        coeP = nSub-1
        tau = ( ( abs(beta_S3-beta_S1) + abs(beta_S3-beta_S2) )/(nSub-1) )**coeP
        
        !====================================> Weights & Combination
        alpha(S1) = r13 * ( 1 + tau/(epsilon + beta_S1) )
        alpha(S2) = r23 * ( 1 + tau/(epsilon + beta_S2) )
        alpha(S3) = r33 * ( 1 + tau/(epsilon + beta_S3) )
           
        alphaSum = sum(alpha)
        w_NL = alpha/alphaSum
        var_Out(iVar) = w_NL(1)*pL_S1  +  w_NL(2)*pL_S2  +  w_NL(3)*pL_S3
    End Do
End Subroutine
        
        
        
        
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
Subroutine WCNS_MR_O7 ( nVar, var_In, var_Out )
use Parameters
implicit none
integer,parameter:: Order = 7            !> Accuracy
integer,parameter:: nSub = (Order+1)/2   !> Number of sub-stencils
integer,parameter:: nDrv =  Order-1      !> Highest derivative order for SI
integer,INTENT(IN ):: nVar
real,INTENT(IN ):: var_In(nVar,Order) !> optimal模板
real,INTENT(OUT):: var_Out(nVar)         !> qL
real,parameter:: epsilon = 1E-10    !> WENO-MR
integer:: iVar, iDrv, iSub, coeP
real:: varS1(nVar)  , qL_S1, pL_S1, beta_S1
real:: varS2(nVar,3), qL_S2, pL_S2, beta_S2
real:: varS3(nVar,5), qL_S3, pL_S3, beta_S3
real:: varS4(nVar,7), qL_S4, pL_S4, beta_S4
real:: Drv_qS(nDrv,nSub)
real:: Drv_pS(nDrv,nSub)
real:: Tau, alpha(nSub), alphaSum, w_NL(nSub)
real:: xiL(nVar), xiR(nVar)

real:: MR_IntS2(3), MR_SI_S2(3,2)
real:: MR_IntS3(5), MR_SI_S3(5,4)
real:: MR_IntS4(7), MR_SI_S4(7,6)
!-------------------------------------------------
    MR_IntS2(1:3   ) = (/  -1./8,  3./4,  3./8  /)   !> S2 Interpolation
    MR_SI_S2(1:3, 1) = (/  -1./2,   0.0,  1./2  /)  !> S2 SI_u'
    MR_SI_S2(1:3, 2) = (/    1.0,  -2.0,   1.0  /)  !> S2 SI_u''
                  
    MR_IntS3(1:5   ) = (/  3./128,  -5./32,  45./64,  15./32,  -5./128  /)   !> S3 Interpolation
    MR_SI_S3(1:5, 1) = (/   1./12,  -2./3 ,    0.0 ,   2./3 ,  -1./12   /)   !> S3 SI_u'
    MR_SI_S3(1:5, 2) = (/  -1./12,   4./3 ,  -5./2 ,   4./3 ,  -1./12   /)   !> S3 SI_u''
    MR_SI_S3(1:5, 3) = (/   -1./2,    1.0 ,    0.0 ,   -1.0 ,   1./2    /)   !> S3 SI_u'''
    MR_SI_S3(1:5, 4) = (/     1.0,   -4.0 ,    6.0 ,   -4.0 ,    1.0    /)   !> S3 SI_u''''
    
    MR_IntS4(1:7   ) = (/  -5./1024,  21./512,  -175./1024,  175./256,  525./1024,  -35./512,   7./1024  /)    !> S4 Interpolation
    MR_SI_S4(1:7, 1) = (/  -1./60  ,   3./20 ,    -3./4   ,      0.0 ,    3./4   ,   -3./20 ,   1./60    /)    !> S4 SI_u'
    MR_SI_S4(1:7, 2) = (/   1./90  ,  -3./20 ,     3./2   ,  -49./18 ,    3./2   ,   -3./20 ,   1./90    /)    !> S4 SI_u''
    MR_SI_S4(1:7, 3) = (/   1./8   ,   -1.0  ,    13./8   ,      0.0 ,  -13./8   ,      1.0 ,  -1./8     /)    !> S4 SI_u'''
    MR_SI_S4(1:7, 4) = (/  -1./6   ,    2.0  ,   -13./2   ,   28./3  ,  -13./2   ,      2.0 ,  -1./6     /)    !> S4 SI_u''''
    MR_SI_S4(1:7, 5) = (/  -1./2   ,    2.0  ,    -5./2   ,      0.0 ,    5./2   ,     -2.0 ,   1./2     /)    !> S4 SI_u'''''
    MR_SI_S4(1:7, 6) = (/    1.0   ,   -6.0  ,     15.0   ,    -20.0 ,    15.0   ,     -6.0 ,    1.0     /)    !> S4 SI_u''''''
!-------------------------------------------------


    !===========================> Interior FP
        varS1(:)      = var_In(:,  4 )    !> S1 = {          n4 |          }
        varS2(:, 1:3) = var_In(:, 3:5)    !> S2 = {       n3 n4 | n5       }   
        varS3(:, 1:5) = var_In(:, 2:6)    !> S3 = {    n2 n3 n4 | n5 n6    }
        varS4(:, 1:7) = var_In(:,  : )    !> S4 = { n1 n2 n3 n4 | n5 n6 n7 }
        
        xiL(:) = varS2(:,2) - varS2(:,1)
        xiR(:) = varS2(:,3) - varS2(:,2)
    
    
    !======================================================>!======================================================>
    Do iVar=1,nVar
        !====================================> poly_Q
        qL_S1 = varS1(iVar)                                !> Sub-stencil Interpolation: S1
        qL_S2 = sum( MR_IntS2(:) * varS2(iVar,:) )    !> Sub-stencil Interpolation: S2
        qL_S3 = sum( MR_IntS3(:) * varS3(iVar,:) )    !> Sub-stencil Interpolation: S3
        qL_S4 = sum( MR_IntS4(:) * varS4(iVar,:) )    !> Sub-stencil Interpolation: S4
        
        Drv_qS = 0.0
        Do iDrv = 1,2
            Drv_qS(iDrv,S2) = sum( MR_SI_S2(:,iDrv) * varS2(iVar,:) )  !> poly_Q Drv: S2
        End Do
        Do iDrv = 1,4
            Drv_qS(iDrv,S3) = sum( MR_SI_S3(:,iDrv) * varS3(iVar,:) )  !> poly_Q Drv: S3
        End Do
        Do iDrv = 1,6
            Drv_qS(iDrv,S4) = sum( MR_SI_S4(:,iDrv) * varS4(iVar,:) )  !> poly_Q Drv: S4
        End Do
        
        !====================================> poly_P
        pL_S1 = qL_S1                                     
        pL_S2 = 1./r22 * ( qL_S2 - r12*pL_S1 )            
        pL_S3 = 1./r33 * ( qL_S3 - r13*pL_S1 - r23*pL_S2 )
        pL_S4 = 1./r44 * ( qL_S4 - r14*pL_S1 - r24*pL_S2 - r34*pL_S3 )
        
        !> poly_P derivatives
        Drv_pS = 0.0
        Drv_pS(:,S2) = 1./r22 * ( Drv_qS(:,S2) )
        Drv_pS(:,S3) = 1./r33 * ( Drv_qS(:,S3) - r23*Drv_pS(:,S2) )
        Drv_pS(:,S4) = 1./r44 * ( Drv_qS(:,S4) - r24*Drv_pS(:,S2) - r34*Drv_pS(:,S3) )
        
        !> Lambda: scale derivatives
        Do iSub = 1,nSub
        Do iDrv = 1,2*(iSub-1)
            coeP = iDrv-1
            Drv_pS(iDrv,iSub) = lambda**coeP *Drv_pS(iDrv,iSub)**2
        End Do
        End Do
        
        !====================================> Beta & Tau
        !> Smoothness Indicators
        call Cal_beta_S1 ( xiL(iVar), xiR(iVar), beta_S1 )
        beta_S2 = sum( Drv_pS(1:2 ,S2) )
        beta_S3 = sum( Drv_pS(1:4 ,S3) )
        beta_S4 = sum( Drv_pS(1:6 ,S4) )
        
        !> Tau
        coeP = nSub-1
        tau = ( ( abs(beta_S4-beta_S1) + abs(beta_S4-beta_S2) + abs(beta_S4-beta_S3) )/(nSub-1) )**coeP
        
        !====================================> Weights & Combination
        alpha(S1) = r14 * ( 1 + tau/(epsilon + beta_S1) )
        alpha(S2) = r24 * ( 1 + tau/(epsilon + beta_S2) )
        alpha(S3) = r34 * ( 1 + tau/(epsilon + beta_S3) )
        alpha(S4) = r44 * ( 1 + tau/(epsilon + beta_S4) )
        
        alphaSum = sum(alpha)
        w_NL = alpha/alphaSum
        var_Out(iVar) = w_NL(1)*pL_S1  +  w_NL(2)*pL_S2  +  w_NL(3)*pL_S3  +  w_NL(4)*pL_S4
    End Do
End Subroutine
        
        
        
        
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
Subroutine WCNS_MR_O9 ( nVar, var_In, var_Out )
use Parameters
implicit none
integer,parameter:: Order = 9            !> Accuracy
integer,parameter:: nSub = (Order+1)/2   !> Number of sub-stencils
integer,parameter:: nDrv =  Order-1      !> Highest derivative order for SI
integer,INTENT(IN ):: nVar
real,INTENT(IN ):: var_In(nVar,Order) !> optimal模板
real,INTENT(OUT):: var_Out(nVar)         !> qL
real,parameter:: epsilon = 1E-10    !> WENO-MR
integer:: iVar, iDrv, iSub, coeP
real:: varS1(nVar)  , qL_S1, pL_S1, beta_S1
real:: varS2(nVar,3), qL_S2, pL_S2, beta_S2
real:: varS3(nVar,5), qL_S3, pL_S3, beta_S3
real:: varS4(nVar,7), qL_S4, pL_S4, beta_S4
real:: varS5(nVar,9), qL_S5, pL_S5, beta_S5
real:: Drv_qS(nDrv,nSub)
real:: Drv_pS(nDrv,nSub)
real:: Tau, alpha(nSub), alphaSum, w_NL(nSub)
real:: xiL(nVar), xiR(nVar)

real:: MR_IntS2(3), MR_SI_S2(3,2)
real:: MR_IntS3(5), MR_SI_S3(5,4)
real:: MR_IntS4(7), MR_SI_S4(7,6)
real:: MR_IntS5(9), MR_SI_S5(9,8)
!-------------------------------------------------
    MR_IntS2(1:3   ) = (/  -1./8,  3./4,  3./8  /)   !> S2 Interpolation
    MR_SI_S2(1:3, 1) = (/  -1./2,   0.0,  1./2  /)  !> S2 SI_u'
    MR_SI_S2(1:3, 2) = (/    1.0,  -2.0,   1.0  /)  !> S2 SI_u''
                  
    MR_IntS3(1:5   ) = (/  3./128,  -5./32,  45./64,  15./32,  -5./128  /)   !> S3 Interpolation
    MR_SI_S3(1:5, 1) = (/   1./12,  -2./3 ,    0.0 ,   2./3 ,  -1./12   /)   !> S3 SI_u'
    MR_SI_S3(1:5, 2) = (/  -1./12,   4./3 ,  -5./2 ,   4./3 ,  -1./12   /)   !> S3 SI_u''
    MR_SI_S3(1:5, 3) = (/   -1./2,    1.0 ,    0.0 ,   -1.0 ,   1./2    /)   !> S3 SI_u'''
    MR_SI_S3(1:5, 4) = (/     1.0,   -4.0 ,    6.0 ,   -4.0 ,    1.0    /)   !> S3 SI_u''''
    
    MR_IntS4(1:7   ) = (/  -5./1024,  21./512,  -175./1024,  175./256,  525./1024,  -35./512,   7./1024  /)    !> S4 Interpolation
    MR_SI_S4(1:7, 1) = (/  -1./60  ,   3./20 ,    -3./4   ,      0.0 ,    3./4   ,   -3./20 ,   1./60    /)    !> S4 SI_u'
    MR_SI_S4(1:7, 2) = (/   1./90  ,  -3./20 ,     3./2   ,  -49./18 ,    3./2   ,   -3./20 ,   1./90    /)    !> S4 SI_u''
    MR_SI_S4(1:7, 3) = (/   1./8   ,   -1.0  ,    13./8   ,      0.0 ,  -13./8   ,      1.0 ,  -1./8     /)    !> S4 SI_u'''
    MR_SI_S4(1:7, 4) = (/  -1./6   ,    2.0  ,   -13./2   ,   28./3  ,  -13./2   ,      2.0 ,  -1./6     /)    !> S4 SI_u''''
    MR_SI_S4(1:7, 5) = (/  -1./2   ,    2.0  ,    -5./2   ,      0.0 ,    5./2   ,     -2.0 ,   1./2     /)    !> S4 SI_u'''''
    MR_SI_S4(1:7, 6) = (/    1.0   ,   -6.0  ,     15.0   ,    -20.0 ,    15.0   ,     -6.0 ,    1.0     /)    !> S4 SI_u''''''
    
    MR_IntS5(1:9   ) = (/  35./32768,  -45./4096,  441./8192,  -735./4096,  11025./16384,  2205./4096,  -735./8192,  63./4096,  -45./32768  /)    !> S5 Interpolation
    MR_SI_S5(1:9, 1) = (/   1./280,  -4./105,     1./5  ,    -4./5 ,       0.0,     4./5 ,   -1./5  ,   4./105,  -1./280                    /)    !> S5 SI_u'
    MR_SI_S5(1:9, 2) = (/  -1./560,   8./315,    -1./5  ,     8./5 ,  -205./72,     8./5 ,   -1./5  ,   8./315,  -1./560                    /)    !> S5 SI_u''
    MR_SI_S5(1:9, 3) = (/  -7./240,   3./10 ,  -169./120,    61./30,       0.0,   -61./30,  169./120,  -3./10 ,   7./240                    /)    !> S5 SI_u'''
    MR_SI_S5(1:9, 4) = (/   7./240,  -2./5  ,   169./60 ,  -122./15,    91./8 ,  -122./15,  169./60 ,  -2./5  ,   7./240                    /)    !> S5 SI_u''''
    MR_SI_S5(1:9, 5) = (/   1./6  ,  -3./2  ,    13./3  ,   -29./6 ,       0.0,    29./6 ,  -13./3  ,   3./2  ,  -1./6                      /)    !> S5 SI_u'''''
    MR_SI_S5(1:9, 6) = (/  -1./4  ,      3.0,      -13.0,      29.0,   -75./2 ,      29.0,     -13.0,      3.0,  -1./4                      /)    !> S5 SI_u''''''
    MR_SI_S5(1:9, 7) = (/  -1./2  ,      3.0,       -7.0,       7.0,       0.0,      -7.0,       7.0,     -3.0,   1./2                      /)    !> S5 SI_u'''''''
    MR_SI_S5(1:9, 8) = (/      1.0,     -8.0,       28.0,     -56.0,      70.0,     -56.0,      28.0,     -8.0,    1.0                      /)    !> S5 SI_u''''''''
!-------------------------------------------------


    !===========================> Interior FP
        varS1(:)      = var_In(:,  5 )    !> S1 = {             n5 |             }
        varS2(:, 1:3) = var_In(:, 4:6)    !> S2 = {          n4 n5 | n6          }    --> xiL,xiR
        varS3(:, 1:5) = var_In(:, 3:7)    !> S3 = {       n3 n4 n5 | n6 n7       }
        varS4(:, 1:7) = var_In(:, 2:8)    !> S4 = {    n2 n3 n4 n5 | n6 n7 n8    }
        varS5(:, 1:9) = var_In(:,  : )    !> S5 = { n1 n2 n3 n4 n5 | n6 n7 n8 n9 }
        
        xiL(:) = varS2(:,2) - varS2(:,1)
        xiR(:) = varS2(:,3) - varS2(:,2)
    
    
    !======================================================>!======================================================>
    Do iVar=1,nVar
        !====================================> poly_Q
        qL_S1 = varS1(iVar)                                !> Sub-stencil Interpolation: S1
        qL_S2 = sum( MR_IntS2(:) * varS2(iVar,:) )    !> Sub-stencil Interpolation: S2
        qL_S3 = sum( MR_IntS3(:) * varS3(iVar,:) )    !> Sub-stencil Interpolation: S3
        qL_S4 = sum( MR_IntS4(:) * varS4(iVar,:) )    !> Sub-stencil Interpolation: S4
        qL_S5 = sum( MR_IntS5(:) * varS5(iVar,:) )    !> Sub-stencil Interpolation: S5
        
        Drv_qS = 0.0
        Do iDrv = 1,2
            Drv_qS(iDrv,S2) = sum( MR_SI_S2(:,iDrv) * varS2(iVar,:) )  !> poly_Q Drv: S2
        End Do
        Do iDrv = 1,4
            Drv_qS(iDrv,S3) = sum( MR_SI_S3(:,iDrv) * varS3(iVar,:) )  !> poly_Q Drv: S3
        End Do
        Do iDrv = 1,6
            Drv_qS(iDrv,S4) = sum( MR_SI_S4(:,iDrv) * varS4(iVar,:) )  !> poly_Q Drv: S4
        End Do
        Do iDrv = 1,8
            Drv_qS(iDrv,S5) = sum( MR_SI_S5(:,iDrv) * varS5(iVar,:) )  !> poly_Q Drv: S5
        End Do
        
        !====================================> poly_P
        pL_S1 = qL_S1                                     
        pL_S2 = 1./r22 * ( qL_S2 - r12*pL_S1 )            
        pL_S3 = 1./r33 * ( qL_S3 - r13*pL_S1 - r23*pL_S2 )
        pL_S4 = 1./r44 * ( qL_S4 - r14*pL_S1 - r24*pL_S2 - r34*pL_S3 )
        pL_S5 = 1./r55 * ( qL_S5 - r15*pL_S1 - r25*pL_S2 - r35*pL_S3 - r45*pL_S4 )
        
        !> poly_P derivatives
        Drv_pS = 0.0
        Drv_pS(:,S2) = 1./r22 * ( Drv_qS(:,S2) )
        Drv_pS(:,S3) = 1./r33 * ( Drv_qS(:,S3) - r23*Drv_pS(:,S2) )
        Drv_pS(:,S4) = 1./r44 * ( Drv_qS(:,S4) - r24*Drv_pS(:,S2) - r34*Drv_pS(:,S3) )
        Drv_pS(:,S5) = 1./r55 * ( Drv_qS(:,S5) - r25*Drv_pS(:,S2) - r35*Drv_pS(:,S3) - r45*Drv_pS(:,S4) )
        
        !> Lambda: scale derivatives
        Do iSub = 1,nSub
        Do iDrv = 1,2*(iSub-1)
            coeP = iDrv-1
            Drv_pS(iDrv,iSub) = lambda**coeP *Drv_pS(iDrv,iSub)**2
        End Do
        End Do
        
        !====================================> Beta & Tau
        !> Smoothness Indicators
        call Cal_beta_S1 ( xiL(iVar), xiR(iVar), beta_S1 )
        beta_S2 = sum( Drv_pS(1:2 ,S2) )
        beta_S3 = sum( Drv_pS(1:4 ,S3) )
        beta_S4 = sum( Drv_pS(1:6 ,S4) )
        beta_S5 = sum( Drv_pS(1:8 ,S5) )
        
        !> Tau
        coeP = nSub-1
        tau = ( ( abs(beta_S5-beta_S1) + abs(beta_S5-beta_S2) + abs(beta_S5-beta_S3) + abs(beta_S5-beta_S4) )/(nSub-1) )**coeP
        
        !====================================> Weights & Combination
        alpha(S1) = r15 * ( 1 + tau/(epsilon + beta_S1) )
        alpha(S2) = r25 * ( 1 + tau/(epsilon + beta_S2) )
        alpha(S3) = r35 * ( 1 + tau/(epsilon + beta_S3) )
        alpha(S4) = r45 * ( 1 + tau/(epsilon + beta_S4) )
        alpha(S5) = r55 * ( 1 + tau/(epsilon + beta_S5) )
        
        alphaSum = sum(alpha)
        w_NL = alpha/alphaSum
        var_Out(iVar) = w_NL(1)*pL_S1  +  w_NL(2)*pL_S2  +  w_NL(3)*pL_S3  +  w_NL(4)*pL_S4  +  w_NL(5)*pL_S5
    End Do
End Subroutine
        
        
        
        
!!==============================================================!!
!!                           WCNS-MR                            !!
!!==============================================================!!
Subroutine WCNS_MR_O11 ( nVar, var_In, var_Out )
use Parameters
implicit none
integer,parameter:: Order = 11           !> Accuracy
integer,parameter:: nSub = (Order+1)/2   !> Number of sub-stencils
integer,parameter:: nDrv =  Order-1      !> Highest derivative order for SI
integer,INTENT(IN ):: nVar
real,INTENT(IN ):: var_In(nVar,Order) !> optimal模板
real,INTENT(OUT):: var_Out(nVar)         !> qL
real,parameter:: epsilon = 1E-10    !> WENO-MR
integer:: iVar, iDrv, iSub, coeP
real:: varS1(nVar)   , qL_S1, pL_S1, beta_S1
real:: varS2(nVar, 3), qL_S2, pL_S2, beta_S2
real:: varS3(nVar, 5), qL_S3, pL_S3, beta_S3
real:: varS4(nVar, 7), qL_S4, pL_S4, beta_S4
real:: varS5(nVar, 9), qL_S5, pL_S5, beta_S5
real:: varS6(nVar,11), qL_S6, pL_S6, beta_S6
real:: Drv_qS(nDrv,nSub)
real:: Drv_pS(nDrv,nSub)
real:: Tau, alpha(nSub), alphaSum, w_NL(nSub)
real:: xiL(nVar), xiR(nVar)

real:: MR_IntS2( 3), MR_SI_S2( 3, 2)
real:: MR_IntS3( 5), MR_SI_S3( 5, 4)
real:: MR_IntS4( 7), MR_SI_S4( 7, 6)
real:: MR_IntS5( 9), MR_SI_S5( 9, 8)
real:: MR_IntS6(11), MR_SI_S6(11,10)
!-------------------------------------------------
    MR_IntS2(1:3   ) = (/  -1./8,  3./4,  3./8  /)   !> S2 Interpolation
    MR_SI_S2(1:3, 1) = (/  -1./2,   0.0,  1./2  /)  !> S2 SI_u'
    MR_SI_S2(1:3, 2) = (/    1.0,  -2.0,   1.0  /)  !> S2 SI_u''
                  
    MR_IntS3(1:5   ) = (/  3./128,  -5./32,  45./64,  15./32,  -5./128  /)   !> S3 Interpolation
    MR_SI_S3(1:5, 1) = (/   1./12,  -2./3 ,    0.0 ,   2./3 ,  -1./12   /)   !> S3 SI_u'
    MR_SI_S3(1:5, 2) = (/  -1./12,   4./3 ,  -5./2 ,   4./3 ,  -1./12   /)   !> S3 SI_u''
    MR_SI_S3(1:5, 3) = (/   -1./2,    1.0 ,    0.0 ,   -1.0 ,   1./2    /)   !> S3 SI_u'''
    MR_SI_S3(1:5, 4) = (/     1.0,   -4.0 ,    6.0 ,   -4.0 ,    1.0    /)   !> S3 SI_u''''
    
    MR_IntS4(1:7   ) = (/  -5./1024,  21./512,  -175./1024,  175./256,  525./1024,  -35./512,   7./1024  /)    !> S4 Interpolation
    MR_SI_S4(1:7, 1) = (/  -1./60  ,   3./20 ,    -3./4   ,      0.0 ,    3./4   ,   -3./20 ,   1./60    /)    !> S4 SI_u'
    MR_SI_S4(1:7, 2) = (/   1./90  ,  -3./20 ,     3./2   ,  -49./18 ,    3./2   ,   -3./20 ,   1./90    /)    !> S4 SI_u''
    MR_SI_S4(1:7, 3) = (/   1./8   ,   -1.0  ,    13./8   ,      0.0 ,  -13./8   ,      1.0 ,  -1./8     /)    !> S4 SI_u'''
    MR_SI_S4(1:7, 4) = (/  -1./6   ,    2.0  ,   -13./2   ,   28./3  ,  -13./2   ,      2.0 ,  -1./6     /)    !> S4 SI_u''''
    MR_SI_S4(1:7, 5) = (/  -1./2   ,    2.0  ,    -5./2   ,      0.0 ,    5./2   ,     -2.0 ,   1./2     /)    !> S4 SI_u'''''
    MR_SI_S4(1:7, 6) = (/    1.0   ,   -6.0  ,     15.0   ,    -20.0 ,    15.0   ,     -6.0 ,    1.0     /)    !> S4 SI_u''''''
    
    MR_IntS5(1:9   ) = (/  35./32768,  -45./4096,  441./8192,  -735./4096,  11025./16384,  2205./4096,  -735./8192,  63./4096,  -45./32768  /)    !> S5 Interpolation
    MR_SI_S5(1:9, 1) = (/   1./280,  -4./105,     1./5  ,    -4./5 ,       0.0,     4./5 ,   -1./5  ,   4./105,  -1./280                    /)    !> S5 SI_u'
    MR_SI_S5(1:9, 2) = (/  -1./560,   8./315,    -1./5  ,     8./5 ,  -205./72,     8./5 ,   -1./5  ,   8./315,  -1./560                    /)    !> S5 SI_u''
    MR_SI_S5(1:9, 3) = (/  -7./240,   3./10 ,  -169./120,    61./30,       0.0,   -61./30,  169./120,  -3./10 ,   7./240                    /)    !> S5 SI_u'''
    MR_SI_S5(1:9, 4) = (/   7./240,  -2./5  ,   169./60 ,  -122./15,    91./8 ,  -122./15,  169./60 ,  -2./5  ,   7./240                    /)    !> S5 SI_u''''
    MR_SI_S5(1:9, 5) = (/   1./6  ,  -3./2  ,    13./3  ,   -29./6 ,       0.0,    29./6 ,  -13./3  ,   3./2  ,  -1./6                      /)    !> S5 SI_u'''''
    MR_SI_S5(1:9, 6) = (/  -1./4  ,      3.0,      -13.0,      29.0,   -75./2 ,      29.0,     -13.0,      3.0,  -1./4                      /)    !> S5 SI_u''''''
    MR_SI_S5(1:9, 7) = (/  -1./2  ,      3.0,       -7.0,       7.0,       0.0,      -7.0,       7.0,     -3.0,   1./2                      /)    !> S5 SI_u'''''''
    MR_SI_S5(1:9, 8) = (/      1.0,     -8.0,       28.0,     -56.0,      70.0,     -56.0,      28.0,     -8.0,    1.0                      /)    !> S5 SI_u''''''''
    
    MR_IntS6(1:11    ) = (/  -63./262144,  385./131072,  -4455./262144,  2079./32768,  -24255./131072,  43659./65536,  72765./131072,  -3465./32768,  6237./262144,  -495./131072,  77./262144  /)    !> S6 Interpolation
    MR_SI_S6(1:11,  1) = (/   -1./1260,      5./504  ,    -5./84  ,      5./21  ,     -5./6  ,      0.0    ,      5./6  ,    -5./21  ,     5./84  ,    -5./504  ,    1./1260  /)    !> S6 SI_u'
    MR_SI_S6(1:11,  2) = (/    1./3150,     -5./1008 ,     5./126 ,     -5./21  ,      5./3  ,  -5269./1800,      5./3  ,    -5./21  ,     5./126 ,    -5./1008 ,    1./3150  /)    !> S6 SI_u''
    MR_SI_S6(1:11,  3) = (/   41./6048,  -1261./15120,   541./1120,  -4369./2520,   1669./720,      0.0    ,  -1669./720,  4369./2520,  -541./1120,  1261./15120,  -41./6048  /)    !> S6 SI_u'''
    MR_SI_S6(1:11,  4) = (/  -41./7560,   1261./15120,  -541./840 ,   4369./1260,  -1669./180,   1529./120 ,  -1669./180,  4369./1260,  -541./840 ,  1261./15120,  -41./7560  /)    !> S6 SI_u''''
    MR_SI_S6(1:11,  5) = (/  -13./288 ,     19./36   ,   -87./32  ,     13./2   ,   -323./48 ,      0.0    ,    323./48 ,   -13./2   ,    87./32  ,   -19./36   ,   13./288   /)    !> S6 SI_u'''''
    MR_SI_S6(1:11,  6) = (/   13./240 ,    -19./24   ,    87./16  ,    -39./2   ,    323./8  ,  -1023./20  ,    323./8  ,   -39./2   ,    87./16  ,   -19./24   ,   13./240   /)    !> S6 SI_u''''''
    MR_SI_S6(1:11,  7) = (/    5./24  ,    -13./6    ,    69./8   ,     -17.0   ,     63./4  ,      0.0    ,    -63./4  ,     17.0   ,   -69./8   ,    13./6    ,   -5./24    /)    !> S6 SI_u'''''''
    MR_SI_S6(1:11,  8) = (/   -1./3   ,     13./3    ,    -23.0   ,      68.0   ,    -126.0  ,    154.0    ,    -126.0  ,     68.0   ,    -23.0   ,    13./3    ,   -1./3     /)    !> S6 SI_u''''''''
    MR_SI_S6(1:11,  9) = (/   -1./2   ,       4.0    ,   -27./2   ,      24.0   ,     -21.0  ,      0.0    ,      21.0  ,    -24.0   ,    27./2   ,     -4.0    ,    1./2     /)    !> S6 SI_u'''''''''
    MR_SI_S6(1:11, 10) = (/     1.0   ,     -10.0    ,     45.0   ,    -120.0   ,     210.0  ,   -252.0    ,     210.0  ,   -120.0   ,     45.0   ,    -10.0    ,     1.0     /)    !> S6 SI_u''''''''''
!-------------------------------------------------


    !===========================> Interior FP
        varS1(:)       = var_In(:,  6  )    !> S1 = {                n6 |                  }
        varS2(:, 1:3 ) = var_In(:, 5:7 )    !> S2 = {             n5 n6 | n7               }    --> xiL,xiR
        varS3(:, 1:5 ) = var_In(:, 4:8 )    !> S3 = {          n4 n5 n6 | n7 n8            }
        varS4(:, 1:7 ) = var_In(:, 3:9 )    !> S4 = {       n3 n4 n5 n6 | n7 n8 n9         }
        varS5(:, 1:9 ) = var_In(:, 2:10)    !> S5 = {    n2 n3 n4 n5 n6 | n7 n8 n9 n10     }
        varS6(:, 1:11) = var_In(:,  :  )    !> S6 = { n1 n2 n3 n4 n5 n6 | n7 n8 n9 n10 n11 }
        
        xiL(:) = varS2(:,2) - varS2(:,1)
        xiR(:) = varS2(:,3) - varS2(:,2)
    
    
    !======================================================>!======================================================>
    Do iVar=1,nVar
        !====================================> poly_Q
        qL_S1 = varS1(iVar)                                !> Sub-stencil Interpolation: S1
        qL_S2 = sum( MR_IntS2(:) * varS2(iVar,:) )    !> Sub-stencil Interpolation: S2
        qL_S3 = sum( MR_IntS3(:) * varS3(iVar,:) )    !> Sub-stencil Interpolation: S3
        qL_S4 = sum( MR_IntS4(:) * varS4(iVar,:) )    !> Sub-stencil Interpolation: S4
        qL_S5 = sum( MR_IntS5(:) * varS5(iVar,:) )    !> Sub-stencil Interpolation: S5
        qL_S6 = sum( MR_IntS6(:) * varS6(iVar,:) )    !> Sub-stencil Interpolation: S6
        
        Drv_qS = 0.0
        Do iDrv = 1,2
            Drv_qS(iDrv,S2) = sum( MR_SI_S2(:,iDrv) * varS2(iVar,:) )  !> poly_Q Drv: S2
        End Do
        Do iDrv = 1,4
            Drv_qS(iDrv,S3) = sum( MR_SI_S3(:,iDrv) * varS3(iVar,:) )  !> poly_Q Drv: S3
        End Do
        Do iDrv = 1,6
            Drv_qS(iDrv,S4) = sum( MR_SI_S4(:,iDrv) * varS4(iVar,:) )  !> poly_Q Drv: S4
        End Do
        Do iDrv = 1,8
            Drv_qS(iDrv,S5) = sum( MR_SI_S5(:,iDrv) * varS5(iVar,:) )  !> poly_Q Drv: S5
        End Do
        Do iDrv = 1,10
            Drv_qS(iDrv,S6) = sum( MR_SI_S6(:,iDrv) * varS6(iVar,:) )  !> poly_Q Drv: S6
        End Do
        
        !====================================> poly_P
        pL_S1 = qL_S1                                     
        pL_S2 = 1./r22 * ( qL_S2 - r12*pL_S1 )            
        pL_S3 = 1./r33 * ( qL_S3 - r13*pL_S1 - r23*pL_S2 )
        pL_S4 = 1./r44 * ( qL_S4 - r14*pL_S1 - r24*pL_S2 - r34*pL_S3 )
        pL_S5 = 1./r55 * ( qL_S5 - r15*pL_S1 - r25*pL_S2 - r35*pL_S3 - r45*pL_S4 )
        pL_S6 = 1./r66 * ( qL_S6 - r16*pL_S1 - r26*pL_S2 - r36*pL_S3 - r46*pL_S4 - r56*pL_S5 )
        
        !> poly_P derivatives
        Drv_pS = 0.0
        Drv_pS(:,S2) = 1./r22 * ( Drv_qS(:,S2) )
        Drv_pS(:,S3) = 1./r33 * ( Drv_qS(:,S3) - r23*Drv_pS(:,S2) )
        Drv_pS(:,S4) = 1./r44 * ( Drv_qS(:,S4) - r24*Drv_pS(:,S2) - r34*Drv_pS(:,S3) )
        Drv_pS(:,S5) = 1./r55 * ( Drv_qS(:,S5) - r25*Drv_pS(:,S2) - r35*Drv_pS(:,S3) - r45*Drv_pS(:,S4) )
        Drv_pS(:,S6) = 1./r66 * ( Drv_qS(:,S6) - r26*Drv_pS(:,S2) - r36*Drv_pS(:,S3) - r46*Drv_pS(:,S4) - r56*Drv_pS(:,S5) )
        
        !> Lambda: scale derivatives
        Do iSub = 1,nSub
        Do iDrv = 1,2*(iSub-1)
            coeP = iDrv-1
            Drv_pS(iDrv,iSub) = lambda**coeP *Drv_pS(iDrv,iSub)**2
        End Do
        End Do
        
        !====================================> Beta & Tau
        !> Smoothness Indicators
        call Cal_beta_S1 ( xiL(iVar), xiR(iVar), beta_S1 )
        beta_S2 = sum( Drv_pS(1:2 ,S2) )
        beta_S3 = sum( Drv_pS(1:4 ,S3) )
        beta_S4 = sum( Drv_pS(1:6 ,S4) )
        beta_S5 = sum( Drv_pS(1:8 ,S5) )
        beta_S6 = sum( Drv_pS(1:10,S6) )
        
        !> Tau
        coeP = nSub-1
        tau = ( ( abs(beta_S6-beta_S1) + abs(beta_S6-beta_S2) + abs(beta_S6-beta_S3) + abs(beta_S6-beta_S4) + abs(beta_S6-beta_S5) )/(nSub-1) )**coeP
        
        !====================================> Weights & Combination
        alpha(S1) = r16 * ( 1 + tau/(epsilon + beta_S1) )
        alpha(S2) = r26 * ( 1 + tau/(epsilon + beta_S2) )
        alpha(S3) = r36 * ( 1 + tau/(epsilon + beta_S3) )
        alpha(S4) = r46 * ( 1 + tau/(epsilon + beta_S4) )
        alpha(S5) = r56 * ( 1 + tau/(epsilon + beta_S5) )
        alpha(S6) = r66 * ( 1 + tau/(epsilon + beta_S6) )
        
        
        alphaSum = sum(alpha)
        w_NL = alpha/alphaSum
        var_Out(iVar) = w_NL(1)*pL_S1  +  w_NL(2)*pL_S2  +  w_NL(3)*pL_S3  +  w_NL(4)*pL_S4  +  w_NL(5)*pL_S5  +  w_NL(6)*pL_S6
    End Do
End Subroutine