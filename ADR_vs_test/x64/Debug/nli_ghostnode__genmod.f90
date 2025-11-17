        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct  2 21:12:23 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NLI_GHOSTNODE__genmod
          INTERFACE 
            SUBROUTINE NLI_GHOSTNODE(STENCIL,PRIML,PRIMR)
              USE CASESETUP, ONLY :                                     &
     &          OINT
              REAL(KIND=4), INTENT(IN) :: STENCIL(OINT+1)
              REAL(KIND=4), INTENT(OUT) :: PRIML
              REAL(KIND=4), INTENT(OUT) :: PRIMR
            END SUBROUTINE NLI_GHOSTNODE
          END INTERFACE 
        END MODULE NLI_GHOSTNODE__genmod
