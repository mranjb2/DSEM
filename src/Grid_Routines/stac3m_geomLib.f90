!> \file stac3m_geomLib.f90
!! Contains routines pertaining to setting up Gauss and Lobatto points.

!> @brief Set the computational points within each element
!!
!! Set up the spaces of gauss-lobatto points and gauss points in computational domain    
      SUBROUTINE set_points(d)
!
!.........................................................
!     set the computational space points on this subdomain
!.........................................................
!
      USE domain_definition
      USE constants
!
      implicit double precision (a-h,o-z) 
!
      TYPE(domain)          :: d
      INTEGER, DIMENSION(3) :: na
!
!     --------------------------
!     set up gauss-lobatto grids                                            
!     --------------------------
!
      d%cx = 0.0d0
      na = d%nc - 1
      DO k = 1,3
         DO n = 1,d%nc(k) 
            xx = cos((n - 1)*pi/na(k)) 
            d%cx(n,k) = 0.5d0*(1.d0 - xx)
         END DO
      END DO
!
!     ------------------
!     set up gauss grids                                                    
!     ------------------
!
      d%cxg = 0.0d0
      DO k = 1,3
         DO n = 1,d%ncg(k)
            xx = cos((2.d0*n-1.d0)*pi/(2.d0*(d%ncg(k) - 1.d0) + 2.d0)) 
            d%cxg(n,k) = 0.5d0*(1.d0 - xx)
         END DO
      END DO
!
      RETURN
      END SUBROUTINE set_points
!
!> @brief Set the computational points within each mortar
!!
!! Set up the computational spaces of gauss points in every mortar                            
      SUBROUTINE mortar_gauss_points(ncg,cxg,weight)
!
!.........................................................
!     compute the computational space gauss points
!.........................................................
!
      USE size
      USE constants
!
      implicit double precision (a-h,o-z) 
!
      INTEGER, DIMENSION(2)               :: ncg
      DOUBLE PRECISION, DIMENSION(nmax,2) :: cxg,weight
!
      weight = 0.0d0
      DO k = 1,2
         DO n = 1,ncg(k)
            xx = cos((2.d0*n-1.d0)*pi/(2.d0*(ncg(k) - 1.d0) + 2.d0)) 
            cxg(n,k) = 0.5d0*(1.d0 - xx)
         END DO
      END DO
!
      RETURN
      END SUBROUTINE mortar_gauss_points
