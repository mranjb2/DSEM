!///////////////////////////////////////////////////////////////////////////////////////////////////////
!////////											////////
!////////	SolnDerivs.f90									////////
!////////											////////
!////////	contains:									////////
!////////											////////
!////////	      SUBROUTINE Gauss_Deriv(gmet,wlgg,wglg,wggl,wx,wy,wz,ncg,dmx,dmy,dmz  	////////
!////////											////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the gauss point derivatives
!!
!
      SUBROUTINE Gauss_Deriv(gmet,wlgg,wglg,wggl,wx,wy,wz,ncg,  &
     &                       dmx,dmy,dmz)                                   
!
!     compute the gauss point derivatives                               
!
!     DATE: 5/25/95                                                     
!     ROUTINES CALLED:dx2dg                                             
!                     dy2dg                                             
!     APPLICABILITY: parabolic problems
!
!
      USE size
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      DOUBLE PRECISION wlgg(nx,ny,nz),wglg(nx,ny,nz),wggl(nx,ny,nz)
      DOUBLE PRECISION wx(nx,ny,nz),wy(nx,ny,nz), wz(nx,ny,nz) 
      DOUBLE PRECISION gmet(3,3,nx,ny,nz)
      DOUBLE PRECISION dmx(nx,nx),dmy(ny,ny),dmz(nz,nz) 
      DOUBLE PRECISION w1x(nx,ny,nz),w2y(nx,ny,nz),w3z(nx,ny,nz) 
      INTEGER          ncg(3)
!
! Initializations
!
      wx = 0.0d0 
      wy = 0.0d0
      wz = 0.0d0
!
! Compute derivatives
!

      CALL dx3dga(wlgg,ncg(:)+1,w1x,ncg(:),-0.5d0,dmx)
      CALL dy3dga(wglg,ncg(:)+1,w2y,ncg(:),-0.5d0,dmy)
      CALL dz3dga(wggl,ncg(:)+1,w3z,ncg(:),-0.5d0,dmz)

!
      DO m = 1,ncg(1) 
         DO n = 1,ncg(2) 
           DO k = 1,ncg(3) 
              wx(n,m,k) = (gmet(1,1,n,m,k)*w1x(n,m,k) + &
                          gmet(2,1,n,m,k)*w2y(n,m,k) + &
                          gmet(3,1,n,m,k)*w3z(n,m,k))
           END DO 
         END DO 
      END DO 

!
      DO m = 1,ncg(1) 
         DO n = 1,ncg(2) 
           DO k = 1,ncg(3) 
              wy(n,m,k) = (gmet(1,2,n,m,k)*w1x(n,m,k) + &
                          gmet(2,2,n,m,k)*w2y(n,m,k) + &
                          gmet(3,2,n,m,k)*w3z(n,m,k))
           END DO 
         END DO 
      END DO 

!
      DO m = 1,ncg(1) 
         DO n = 1,ncg(2) 
           DO k = 1,ncg(3) 
              wz(n,m,k) = (gmet(1,3,n,m,k)*w1x(n,m,k) + &
                          gmet(2,3,n,m,k)*w2y(n,m,k) + &
                          gmet(3,3,n,m,k)*w3z(n,m,k))
           END DO 
         END DO 
      END DO 

!
      RETURN 
      END                                           
