
!&=============================================================================================&!
!&>                                        Parameters                                          &!
!&=============================================================================================&! 
Module Parameters
implicit none
    real,parameter:: machine_zero = 1E-12
    real,parameter:: PI = 3.141592653589793238462643
    real,parameter:: gamma = 1.4
    real,parameter:: oga   = 1.0/gamma   !< 1.0/GAMMA
    real,parameter:: vga   = 0.4         !< (GAMMA - 1.0)
    real,parameter:: ovga  = 2.5         !< 1.0/(GAMMA - 1.0)
    integer,parameter:: ZERO = 0
    integer,parameter::  One = 1
    
    integer,parameter:: Keep =  1
    integer,parameter:: Flip = -1

    integer,parameter:: interiorFP = 1
    integer:: L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,  R10,R9,R8,R7,R6,R5,R4,R3,R2,R1
    integer,parameter:: S1 = 1
    integer,parameter:: S2 = 2
    integer,parameter:: S3 = 3
    integer,parameter:: S4 = 4
    integer,parameter:: S5 = 5
    integer,parameter:: S6 = 6
    integer,parameter,dimension(10) :: Drv = (/ 1,2,3,4,5,6,7,8,9,10 /)
    
    real,parameter:: r12 =  1./11
    real,parameter:: r22 = 10./11
    
    real,parameter:: r13 =   1./111
    real,parameter:: r23 =  10./111
    real,parameter:: r33 = 100./111
    
    real,parameter:: r14 =    1./1111
    real,parameter:: r24 =   10./1111
    real,parameter:: r34 =  100./1111
    real,parameter:: r44 = 1000./1111
    
    real,parameter:: r15 =     1./11111
    real,parameter:: r25 =    10./11111
    real,parameter:: r35 =   100./11111
    real,parameter:: r45 =  1000./11111
    real,parameter:: r55 = 10000./11111
    
    real,parameter:: r16 =      1./111111
    real,parameter:: r26 =     10./111111
    real,parameter:: r36 =    100./111111
    real,parameter:: r46 =   1000./111111
    real,parameter:: r56 =  10000./111111
    real,parameter:: r66 = 100000./111111
    
    real,parameter:: lambda = 1.0
    
End module