        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct  2 21:12:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RHS_GHOSTNODE__genmod
          INTERFACE 
            SUBROUTINE RHS_GHOSTNODE(U_IN,RHS)
              USE VARIABLES_ZONE1, ONLY :                               &
     &          NX,                                                     &
     &          HX
              USE CASESETUP, ONLY :                                     &
     &          ORDER
              REAL(KIND=4), INTENT(IN) :: U_IN(0:NX)
              REAL(KIND=4), INTENT(OUT) :: RHS(0:NX)
            END SUBROUTINE RHS_GHOSTNODE
          END INTERFACE 
        END MODULE RHS_GHOSTNODE__genmod
