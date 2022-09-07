!>\file EViscosity_GL.f90
!! This file contains subroutines to calculate the entropy viscosity.

!> @brief Determine entropy viscosity based on the entropy redidual on the Gauss-Gauss points.
!!
!! \f$ Dh_1 = \frac{\partial S}{\partial t} + \frac{\partial (Fs_x)}{\partial x} + \frac{\partial (Fs_y)}{\partial y} \f$:
!!
!! \f$ Dh_2 = \frac{S}{\rho} (\frac{\partial \rho}{\partial t} + \nabla . (\rho u)) \f$
!!
!! Determine the time derivative with a backward differencing, \f$ \frac{\partial S}{\partial t} \f$  and \f$ \frac{\partial \rho}{\partial t} \f$.
!!
!! Determine the derivative of the entropy flux in x and y,
!! \f$ \frac{\partial Fs_x}{\partial x} \f$ and \f$ \frac{\partial Fs_y}{\partial y} \f$
!!
!! \f$ \nabla . (\rho u) \f$
!!
!! Determine entropy viscosity based on \f$ Dh_1 \f$ and \f$ Dh_2 \f$ and
!! an average grid spacing that pre-processed in Dist.
      SUBROUTINE Entropy_Viscosity(d,ngrid) 
!                    
      USE domain_definition
      USE physics
      USE mpi
      USE input_data

!
      IMPLICIT NONE 
!
      TYPE (domain)                                :: d(ngp)
!
!
      INTEGER                                      :: ngrid,nmort,id,n,m,l
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)        :: rho,umag,T					!< density, velocity magnitude, and temperature
      DOUBLE PRECISION                             :: delta_g, h					!< grid spacing
      DOUBLE PRECISION                             :: biS						    	!< EV residual
      DOUBLE PRECISION                             :: mur       					!< The calculated viscosity
      DOUBLE PRECISION                             :: u, v, w, p					!< velocity components and pressure
      DOUBLE PRECISION                             :: cbeta, ducros				!< \f$ C_{\beta} \f$ and Ducros sensor
      DOUBLE PRECISION                             :: muE,rhomax,UTmax

      DO id = 1,ngrid
!     ------------------------------------------------------
!     convert conservative to primitive variables
!     ------------------------------------------------------
      rhomax = 0.0d0
      UTmax  = 0.0d0

         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               DO n = 1,d(id)%ncg(1)
!               
                  rho(n,m,l)   =  d(id)%Q(n,m,l,1)*d(id)%jacob(n,m,l)
                  u            =  d(id)%Q(n,m,l,2)/d(id)%Q(n,m,l,1)
                  v            =  d(id)%Q(n,m,l,3)/d(id)%Q(n,m,l,1)
                  w            =  d(id)%Q(n,m,l,4)/d(id)%Q(n,m,l,1)
                  umag(n,m,l)  =  sqrt(u**2 + v**2 + w**2)
                  p            = (d(id)%Q(n,m,l,5)*d(id)%jacob(n,m,l) -    &
                                  0.5d0*rho(n,m,l)*umag(n,m,l))*(gamma-1.0d0) 
                  T(n,m,l)     =  p*gamma/rho(n,m,l)
!
!        ---------------------------------------------------
!        calculate maximum density and wave velocity in element
!        ---------------------------------------------------
                  rhomax       = max (rhomax,rho(n,m,l))
                  UTmax        = max (UTmax, umag(n,m,l)+ sqrt(gamma*T(n,m,l))) 
!        ---------------------------------------------------     
               ENDDO
            ENDDO
         ENDDO
    
!     LAD viscosity
!

         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               DO n = 1,d(id)%ncg(1)
                  rho(n,m,l)         = d(id)%Q     (n,m,l,1)*d(id)%jacob(n,m,l)
                  h                  = delta_g     (d(id),n,m,l)
                  biS                = d(id)%RHS_ent(n,m,l)
                  ducros             = d(id)%ducros(n,m,l)
                  mur                = Cmu*re*rho(n,m,l)*h**2.d0*abs(biS)
!                  mumax              = Cmax*rhomax * h * UTmax

                  if (ducrosSensor) mur = mur * ducros
                  
                  d(id)%muhg(n,m,l)  = min(max(mur,0.d0),mumax)
                  
               ENDDO
            ENDDO
         ENDDO
         
         IF ( d(id)%bcond(4) /= 'interface') THEN
            d(id)%muhg (d(id)%ncg(1),:,:) = 0.0d0
         END IF 
  
      END DO    ! END OF ngrid DO LOOP.   
!

      RETURN 
!
      END SUBROUTINE Entropy_Viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
!> @brief This routine calculates the local distance between Gauss points in the physical space.
      SUBROUTINE Dist(d,ngrid)
!
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
          DOUBLE PRECISION :: SIXTH, THIRD, FIFTH
!     ---------------------------------------------------------------------
!     Constants 
!     ---------------------------------------------------------------------
          THIRD  = 1.0/3.0d0
          FIFTH  = 0.20d0
          SIXTH  = 1.0/6.0d0
!
          DO id = 1,ngrid
!     ---------------------------------------------------------------------
!     Corners points
!     ---------------------------------------------------------------------
             n = 1
             m = 1
             l = 1
             d(id)%hx(n,m,l) = THIRD*(      - d(id)%xg(1,n  ,m,l) &
                                            - d(id)%xg(2,n  ,m,l) &
                                            - d(id)%xg(3,n  ,m,l) &                                              
                                            + d(id)%xg(1,n+1,m,l) &
                                            + d(id)%xg(2,n,m+1,l) &
                                            + d(id)%xg(3,n,m,l+1))
             n = d(id)%ncg(1)
             m = 1
             l = 1
             d(id)%hx(n,m,l) = THIRD*(      + d(id)%xg(1,n  ,m,l) &
                                            - d(id)%xg(2,n  ,m,l) &
                                            - d(id)%xg(3,n  ,m,l) &                                              
                                            - d(id)%xg(1,n-1,m,l) &
                                            + d(id)%xg(2,n,m+1,l) &
                                            + d(id)%xg(3,n,m,l+1))

             n = 1
             m = d(id)%ncg(2)
             l = 1
             d(id)%hx(n,m,l) = THIRD*(      - d(id)%xg(1,n  ,m,l) &
                                            + d(id)%xg(2,n  ,m,l) &
                                            - d(id)%xg(3,n  ,m,l) &                                              
                                            + d(id)%xg(1,n+1,m,l) &
                                            - d(id)%xg(2,n,m-1,l) &
                                            + d(id)%xg(3,n,m,l+1))

             n = 1 
             m = 1
             l = d(id)%ncg(3)
             d(id)%hx(n,m,l) = THIRD*(      - d(id)%xg(1,n  ,m,l) &
                                            - d(id)%xg(2,n  ,m,l) &
                                            + d(id)%xg(3,n  ,m,l) &                                              
                                            + d(id)%xg(1,n+1,m,l) &
                                            + d(id)%xg(2,n,m+1,l) &
                                            - d(id)%xg(3,n,m,l-1))

             n = 1
             m = d(id)%ncg(2)
             l = d(id)%ncg(3) 
             d(id)%hx(n,m,l) = THIRD*(      - d(id)%xg(1,n  ,m,l) &
                                            + d(id)%xg(2,n  ,m,l) &
                                            + d(id)%xg(3,n  ,m,l) &                                              
                                            + d(id)%xg(1,n+1,m,l) &
                                            - d(id)%xg(2,n,m-1,l) &
                                            - d(id)%xg(3,n,m,l-1))                                            
                                              
             n = d(id)%ncg(1) 
             m = 1
             l = d(id)%ncg(3)
             d(id)%hx(n,m,l) = THIRD*(      + d(id)%xg(1,n  ,m,l) &
                                            - d(id)%xg(2,n  ,m,l) &
                                            + d(id)%xg(3,n  ,m,l) &                                              
                                            - d(id)%xg(1,n-1,m,l) &
                                            + d(id)%xg(2,n,m+1,l) &
                                            - d(id)%xg(3,n,m,l-1)) 

             n = d(id)%ncg(1)
             m = d(id)%ncg(2)
             l = 1 
             d(id)%hx(n,m,l) = THIRD*(      + d(id)%xg(1,n  ,m,l) &
                                            + d(id)%xg(2,n  ,m,l) &
                                            - d(id)%xg(3,n  ,m,l) &                                              
                                            - d(id)%xg(1,n-1,m,l) &
                                            - d(id)%xg(2,n,m-1,l) &
                                            + d(id)%xg(3,n,m,l+1))                                               
                                              
                                              
             n = d(id)%ncg(1)
             m = d(id)%ncg(2)
             l = d(id)%ncg(3) 
             d(id)%hx(n,m,l) = THIRD*(      + d(id)%xg(1,n  ,m,l) &
                                            + d(id)%xg(2,n  ,m,l) &
                                            + d(id)%xg(3,n  ,m,l) &                                              
                                            - d(id)%xg(1,n-1,m,l) &
                                            - d(id)%xg(2,n,m-1,l) &
                                            - d(id)%xg(3,n,m,l-1))  
!     ---------------------------------------------------------------------
!     Side Plates
!     ---------------------------------------------------------------------
             n = 1
             DO m=2,d(id)%ncg(2)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%hx(n,m,l) = FIFTH*(        - d(id)%xg(1,n  ,m,l) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(2,n,m-1,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 + d(id)%xg(3,n,m,l+1) &
                                                 - d(id)%xg(3,n,m,l-1))
                END DO
             ENDDO

             m = 1
             DO n=2,d(id)%ncg(1)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%hx(n,m,l) = FIFTH*(        - d(id)%xg(2,n  ,m,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 + d(id)%xg(3,n,m,l+1) &
                                                 - d(id)%xg(3,n,m,l-1))
                END DO
             ENDDO

             l = 1
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                d(id)%hx(n,m,l) = FIFTH*(        + d(id)%xg(3,n,m,l+1) &
                                                 - d(id)%xg(3,n,  m,l) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 - d(id)%xg(2,n,m-1,l))
                END DO
             ENDDO

             n = d(id)%ncg(1)
             DO m=2,d(id)%ncg(2)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%hx(n,m,l) = FIFTH*(        + d(id)%xg(1,n  ,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 - d(id)%xg(2,n,m-1,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 + d(id)%xg(3,n,m,l+1) &
                                                 - d(id)%xg(3,n,m,l-1))
                END DO
             ENDDO

             m = d(id)%ncg(2)
             DO n=2,d(id)%ncg(1)-1
                DO l=2,d(id)%ncg(3)-1
                d(id)%hx(n,m,l) = FIFTH*(        + d(id)%xg(2,n  ,m,l) &
                                                 - d(id)%xg(2,n,m-1,l) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 + d(id)%xg(3,n,m,l+1) &
                                                 - d(id)%xg(3,n,m,l-1))
                END DO
             ENDDO

             l = d(id)%ncg(3)
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                d(id)%hx(n,m,l) = FIFTH*(        - d(id)%xg(3,n,m,l-1) &
                                                 + d(id)%xg(3,n,  m,l) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 - d(id)%xg(2,n,m-1,l))
                END DO
             ENDDO




!     ---------------------------------------------------------------------
!     Interior
!     ---------------------------------------------------------------------
               
             DO n=2,d(id)%ncg(1)-1
                DO m=2,d(id)%ncg(2)-1
                   DO l=2,d(id)%ncg(3)-1
                   d(id)%hx(n,m,l) = SIXTH*(     - d(id)%xg(3,n,m,l-1) &
                                                 + d(id)%xg(3,n,m,l+1) &
                                                 + d(id)%xg(1,n+1,m,l) &
                                                 - d(id)%xg(1,n-1,m,l) &
                                                 + d(id)%xg(2,n,m+1,l) &
                                                 - d(id)%xg(2,n,m-1,l))
                   ENDDO
                ENDDO
             ENDDO

          END DO
!
          RETURN  
      END SUBROUTINE Dist    
!
!////////////////////////////////////////////////////////////////////////////////
!
!> @brief This routine calculates the filtered rate of strain.
subroutine Strain(d,ngridl)
    use domain_definition
    use physics
    implicit none
    
    type(domain),dimension(ngridl)              :: d            !< Domain
    integer,intent(in)                          :: ngridl       !< Number of grid points
    integer                                     :: id, eqs, x, y, z, nv
    double precision, dimension(nx,ny,nz,neq)   :: Qx, Qy, Qz
    double precision                            :: filteredStrain
    
    !We loop over the elements of the grid
    do id=1,ngridl-1
        
        ! Now calculate the solution derivatives for each of the equations (to be used)
        do nv=1,neq
            call Gauss_Deriv(d(id)%gmetg,d(id)%Qlgg(:,:,:,nv),d(id)%Qglg(:,:,:,nv),     &
                         d(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv),   &
                         d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz                        )

        end do
        
        ! finally calculate the sgs viscosity
        do z=1,d(id)%ncg(3)
            do y=1,d(id)%ncg(2)
                do x=1,d(id)%ncg(1)
                    d(id)%Strain(x,y,z) = filteredStrain(d(id),Qx,Qy,Qz,x,y,z)


                end do
            end do
        end do
    end do
    
end subroutine Strain