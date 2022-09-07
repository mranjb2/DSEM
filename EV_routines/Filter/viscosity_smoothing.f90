!> @brief Localize and smooth the entropy viscosity by averaging the entropy viscosity at neighbouring grid points
!!
!! Average muhg on Gauss-Gauss points.
      SUBROUTINE entropy_viscosity_smooting(d,ngrid)
!     ---------------------------------------------------------------------
!     Smoothens muhg with FD type filter on spectral grid
!     ---------------------------------------------------------------------
      USE domain_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
      USE mpi
!
          IMPLICIT NONE
!
          TYPE (domain) :: d(ngp)
!
          INTEGER          :: ngrid,n,m,l,id
          DOUBLE PRECISION :: TWELVETH, SIXTH, TENTH
!     ---------------------------------------------------------------------
!     Constants in fd type filters
!     ---------------------------------------------------------------------
          SIXTH     = 1.0/6.0d0
          TENTH     = 0.10d0
          TWELVETH  = 1.0/12.0d0
!
          DO id = 1,ngrid
!     ---------------------------------------------------------------------
!     Corners points
!     ---------------------------------------------------------------------
             n = 1
             m = 1
             l = 1
             d(id)%muhg(1,1,1) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n+1,m,l)   +  &
                                              d(id)%muhg(n,m+1,l)   +  &
                                              d(id)%muhg(n,m,l+1))
             n = d(id)%ncg(1)
             m = 1
             l = 1
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n-1,m,l)   +  &
                                              d(id)%muhg(n,m+1,l)   +  &
                                              d(id)%muhg(n,m,l+1)) 

             n = 1
             m = d(id)%ncg(2)
             l = 1 
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n+1,m,l)   +  &
                                              d(id)%muhg(n,m-1,l)   +  &
                                              d(id)%muhg(n,m,l+1))
                                              
             n = 1 
             m = 1
             l = d(id)%ncg(3)
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n+1,m,l)   +  &
                                              d(id)%muhg(n,m+1,l)   +  &
                                              d(id)%muhg(n,m,l-1)) 

             n = 1
             m = d(id)%ncg(2)
             l = d(id)%ncg(3) 
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n+1,m,l)   +  &
                                              d(id)%muhg(n,m-1,l)   +  &
                                              d(id)%muhg(n,m,l-1))                                              
                                              
             n = d(id)%ncg(1) 
             m = 1
             l = d(id)%ncg(3)
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n-1,m,l)   +  &
                                              d(id)%muhg(n,m+1,l)   +  &
                                              d(id)%muhg(n,m,l-1)) 

             n = d(id)%ncg(1)
             m = d(id)%ncg(2)
             l = 1 
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n-1,m,l)   +  &
                                              d(id)%muhg(n,m-1,l)   +  &
                                              d(id)%muhg(n,m,l+1))                                                 
                                              
                                              
             n = d(id)%ncg(1)
             m = d(id)%ncg(2)
             l = d(id)%ncg(3) 
             d(id)%muhg(n,m,l) = SIXTH*(3.0d0*d(id)%muhg(n  ,m,l)   +  &
                                              d(id)%muhg(n-1,m,l)   +  &
                                              d(id)%muhg(n,m-1,l)   +  &
                                              d(id)%muhg(n,m,l-1))  
!     ---------------------------------------------------------------------
!     Side Plates
!     ---------------------------------------------------------------------
             n = 1
             DO m=2,d(id)%ncg(2)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0* d(id)%muhg(n  ,m,l)   +  &
                                                  d(id)%muhg(n+1,m,l)   +  &
                                                  d(id)%muhg(n,m-1,l)   +  &
                                                  d(id)%muhg(n,m+1,l)   +  &
                                                  d(id)%muhg(n,m,l+1)   +  &
                                                  d(id)%muhg(n,m,l-1))
                END DO
             ENDDO

             m = 1
             DO n=2,d(id)%ncg(1)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0* d(id)%muhg(n  ,m,l)   +  &
                                                  d(id)%muhg(n,m+1,l)   +  &
                                                  d(id)%muhg(n-1,m,l)   +  &
                                                  d(id)%muhg(n+1,m,l)   +  &
                                                  d(id)%muhg(n,m,l+1)   +  &
                                                  d(id)%muhg(n,m,l-1))
                END DO
             ENDDO

             l = 1
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0* d(id)%muhg(n  ,m,l)   +  &
                                                  d(id)%muhg(n,m,l+1)   +  &
                                                  d(id)%muhg(n,m-1,l)   +  &
                                                  d(id)%muhg(n,m+1,l)   +  &
                                                  d(id)%muhg(n-1,m,l)   +  &
                                                  d(id)%muhg(n+1,m,l))
                END DO
             ENDDO

             n = d(id)%ncg(1)
             DO m=2,d(id)%ncg(2)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0*d(id)%muhg(n  ,m,l)    +  &
                                                  d(id)%muhg(n-1,m,l)   +  &
                                                  d(id)%muhg(n,m-1,l)   +  &
                                                  d(id)%muhg(n,m+1,l)   +  &
                                                  d(id)%muhg(n,m,l+1)   +  &
                                                  d(id)%muhg(n,m,l-1))
                END DO
             ENDDO

             m = d(id)%ncg(2)
             DO n=2,d(id)%ncg(1)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0*d(id)%muhg(n  ,m,l)    +  &
                                                  d(id)%muhg(n,m-1,l)   +  &
                                                  d(id)%muhg(n-1,m,l)   +  &
                                                  d(id)%muhg(n+1,m,l)   +  &
                                                  d(id)%muhg(n,m,l+1)   +  &
                                                  d(id)%muhg(n,m,l-1))
                END DO
             ENDDO

             l = d(id)%ncg(3)
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                d(id)%muhg(n,m,l) = TENTH*(5.0d0*d(id)%muhg(n  ,m,l)    +  &
                                                  d(id)%muhg(n,m,l-1)   +  &
                                                  d(id)%muhg(n,m-1,l)   +  &
                                                  d(id)%muhg(n,m+1,l)   +  &
                                                  d(id)%muhg(n-1,m,l)   +  &
                                                  d(id)%muhg(n+1,m,l))
                END DO
             ENDDO




!     ---------------------------------------------------------------------
!     Interior
!     ---------------------------------------------------------------------
               
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                   DO l=2,d(id)%ncg(3)-1
                   d(id)%muhg(n,m,l) = TWELVETH*(6.0d0*d(id)%muhg(n  ,m,l) +  &
                                                       d(id)%muhg(n+1,m,l) +  &
                                                       d(id)%muhg(n-1,m,l) +  &
                                                       d(id)%muhg(n,m+1,l) +  &
                                                       d(id)%muhg(n,m-1,l) +  &
                                                       d(id)%muhg(n,m,l+1) +  &
                                                       d(id)%muhg(n,m,l-1) )
                   ENDDO
                ENDDO
             ENDDO

          END DO
!
          RETURN                                   

         END SUBROUTINE entropy_viscosity_smooting
!
!///////////////////////////////////////////////////////////////////////////////////////////////