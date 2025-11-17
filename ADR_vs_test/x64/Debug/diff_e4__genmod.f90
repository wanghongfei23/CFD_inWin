        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct  2 21:12:22 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIFF_E4__genmod
          INTERFACE 
            SUBROUTINE DIFF_E4(HX,STENCIL,DU_SP)
              REAL(KIND=4), INTENT(IN) :: HX
              REAL(KIND=4), INTENT(IN) :: STENCIL(4)
              REAL(KIND=4), INTENT(OUT) :: DU_SP
            END SUBROUTINE DIFF_E4
          END INTERFACE 
        END MODULE DIFF_E4__genmod
