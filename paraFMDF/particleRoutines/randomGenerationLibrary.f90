!--------------------------------------------------------------------------------------------------
!> \file
!! Random number generation for use with Monte-Carlo Methods
!
! Computational Multiphase Transport Laboratory
! Department of Mechanical Engineering
! University of Illinois at Chicago
! 842 W. Taylor Street
! Lab ERF 2022
! Chicago, Il
!--------------------------------------------------------------------------------------------------
!> Library used for generating random numbers for Monte Carlo simulation
    module randomGenerationLibrary
        implicit none
    
        contains
        
        !--------------------------------------------------------------------------------------
        !> Generates a uniformly random number over the interval [0,1]
        function uniformRandom() result(randomNumber)
        !This function is intended to generate a uniform set of random numbers between 0 and 1
    
            implicit none
        
            integer                            :: seedSize
            integer                            :: clock,i
            integer, dimension(:), allocatable :: seed
            double precision                   :: randomNumber
            
            !Initialize the seed
            call RANDOM_SEED(size=seedSize)
            allocate(seed(seedSize))
            call SYSTEM_CLOCK(COUNT=clock)
            seed = clock + 37 * (/ (i-1,i=1,seedSize) /)
            call RANDOM_SEED(GET=seed)
            deallocate(seed)
        
            !Generate the uniform random number
            call RANDOM_NUMBER(randomNumber)
        end function uniformRandom
        !--------------------------------------------------------------------------------------
        !> Generates a uniformly random number over the interval [-1,1]
        function uniformZeroedRandom() result(randomNumber)
        !This function is intended to generate a uniform set of random numbers between -1 and 1
    
            implicit none
        
            integer                            :: seedSize
            integer                            :: clock,i
            integer, dimension(:), allocatable :: seed
            double precision                   :: randomNumber,randy
            
            !Initialize the seed
            call RANDOM_SEED(size=seedSize)
            allocate(seed(seedSize))
            call SYSTEM_CLOCK(COUNT=clock)
            seed = clock + 37 * (/ (i-1,i=1,seedSize) /)
            call RANDOM_SEED(GET=seed)
            deallocate(seed)
        
            !Generate the uniform random number between 0 and 1
            call RANDOM_NUMBER(randy)
            
            !Scale it between -1 and 1
            randomNumber = (randy*2.0) - 1.0
            
        end function uniformZeroedRandom
        !--------------------------------------------------------------------------------------
        !> Generates a gaussian random number between [0,1]
        function gaussianRandom() result(gaussianNumber)
        !This function uses uniformRandom to generate a gaussian random number between 0 and 1
        
            implicit none
        
            double precision :: gaussianNumber, ran1, ran2
            double precision :: pi
            
            parameter (pi = 4.0*atan(1.0))
            
            !Scale the uniform number to be gaussian
            ran1 = sqrt(2.0*(-alog(1 - real(uniformRandom()))))
            ran2 = 2.0*pi*uniformRandom()
            gaussianNumber = ran1*cos(ran2)
        end function gaussianRandom
    end module randomGenerationLibrary
!--------------------------------------------------------------------------------------------------