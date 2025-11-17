        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct  2 21:12:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIFF_GHOSTNODE__genmod
          INTERFACE 
            SUBROUTINE DIFF_GHOSTNODE(HX,STENCIL,DU_SP)
              USE CASESETUP, ONLY :                                     &
     &          ODIFF
              REAL(KIND=4), INTENT(IN) :: HX
              REAL(KIND=4), INTENT(IN) :: STENCIL(ODIFF)
              REAL(KIND=4), INTENT(OUT) :: DU_SP
            END SUBROUTINE DIFF_GHOSTNODE
          END INTERFACE 
        END MODULE DIFF_GHOSTNODE__genmod
