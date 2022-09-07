!> \file matrices.f90
!! Spectral matrix routines
!///////////////////////////////////////////////////////////////////////
!////////							////////
!////////	matrices.f90					////////
!////////							////////
!////////	contains:					////////
!////////							////////
!////////	   SUBROUTINE set_matrices(d)			////////
!////////	   SUBROUTINE dermat(ncol,no,xo,nnew,xnew,b)	////////
!////////	   SUBROUTINE intrpmat(ncol,no,xo,nnew,xnew,b)	////////
!////////							////////
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets up arrays for chebyshev derivatives.
!!
!!     check to see if these matrices have been set up already
!!     for x, y, and z direction. If it is, point to that array. Otherwise,
!!     add to the list and point to the new one
      SUBROUTINE set_matrices(d) 
!
!.......................................................................
!     date: 2/13/95                                                     
!     routines called: dermat                                           
!                      intrpmat                                         
!     applicability: all stac3m codes
!
!     sets up arrays for chebyshev derivatives.                            
!.......................................................................
!
      USE domain_definition
      USE Order_Matrices
      USE mpi_par

      IMPLICIT none

      TYPE (domain) :: d
      INTEGER       :: k, listloc
      LOGICAL       :: use_copy
!
!     ----------------------------------------------------------
!     check to see if these matrices have been set up already
!     for x direction. If it is, point to that array. Otherwise,
!     add to the list and point to the new one
!     ----------------------------------------------------------
!
      use_copy = .false.
      DO k = 1,num_ordersx
         IF ( order_mapx(k) == d%nc(1) )     THEN
            use_copy = .true.
            listloc = k
            EXIT
         END IF
      END DO

      IF ( use_copy )     THEN
         d%dmx => dmx(:,:,listloc)
         d%bx  => bx(:,:,listloc)
      ELSE
         num_ordersx = num_ordersx + 1
         IF( num_ordersx > max_orders )     THEN
            WRITE(6,*) 'Too many different element orders. Increase max_orders in size'
            STOP 'Too many element orders'
         END IF
         order_mapx(num_ordersx) = d%nc(1)
         CALL dermat  (nx,d%nc(1),d%cx(:,1),d%ncg(1),d%cxg(:,1),dmx(:,:,num_ordersx)) 
         CALL intrpmat(nx,d%ncg(1),d%cxg(:,1),d%nc(1),d%cx(:,1),bx(:,:,num_ordersx)) 

         d%dmx => dmx(:,:,num_ordersx)
         d%bx  => bx (:,:,num_ordersx)
      END IF
!
!     ----------------------------------------------------------
!     check to see if these matrices have been set up already
!     for y direction. If it is, point to that array. Otherwise,
!     add to the list and point to the new one
!     ----------------------------------------------------------
!
      use_copy = .false.
      DO k = 1,num_ordersy
         IF ( order_mapy(k) == d%nc(2) )     THEN
            use_copy = .true.
            listloc = k
            EXIT
         END IF
      END DO

      IF ( use_copy )     THEN
         d%dmy => dmy(:,:,listloc)
         d%by  => by(:,:,listloc)
      ELSE
         num_ordersy = num_ordersy + 1
         IF( num_ordersy > max_orders )     THEN
            WRITE(6,*) 'Too many different element orders. Increase max_orders in size'
            STOP 'Too many element orders'
         END IF
         order_mapy(num_ordersy) = d%nc(2)
         CALL dermat  (ny,d%nc(2),d%cx(:,2),d%ncg(2),d%cxg(:,2),dmy(:,:,num_ordersy)) 
         CALL intrpmat(ny,d%ncg(2),d%cxg(:,2),d%nc(2),d%cx(:,2),by(:,:,num_ordersy)) 
         d%dmy => dmy(:,:,num_ordersy)
         d%by  => by (:,:,num_ordersy)
      END IF
!
!     ----------------------------------------------------------
!     check to see if these matrices have been set up already
!     for z direction. If it is, point to that array. Otherwise,
!     add to the list and point to the new one
!     ----------------------------------------------------------
!
      use_copy = .false.
      DO k = 1,num_ordersz
         IF ( order_mapz(k) == d%nc(3) )     THEN
            use_copy = .true.
            listloc = k
            EXIT
         END IF
      END DO

      IF ( use_copy )     THEN
         d%dmz => dmz(:,:,listloc)
         d%bz  => bz(:,:,listloc)
      ELSE
         num_ordersz = num_ordersz + 1
         IF( num_ordersz > max_orders )     THEN
            WRITE(6,*) 'Too many different element orders. Increase max_orders in size'
            STOP 'Too many element orders'
         END IF
         order_mapz(num_ordersz) = d%nc(3)
         CALL dermat  (nz,d%nc(3),d%cx(:,3),d%ncg(3),d%cxg(:,3),dmz(:,:,num_ordersz)) 
         CALL intrpmat(nz,d%ncg(3),d%cxg(:,3),d%nc(3),d%cx(:,3),bz(:,:,num_ordersz)) 
         d%dmz => dmz(:,:,num_ordersz)
         d%bz  => bz (:,:,num_ordersz)
      END IF
      
!
      RETURN 
      END SUBROUTINE set_matrices
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the differentiation matrix.
!!
!
      SUBROUTINE dermat(ncol,no,xo,nnew,xnew,b) 
!
!.......................................................................
!     compute the differentiation matrix.                                 
!
!     date: 2/13/95                                                     
!     routines called: polynder.f                                       
!     dependencies: constants                                                
!
!     applicability:all                                                 
!.......................................................................
!
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DIMENSION xo(*),xnew(*) 
      DOUBLE PRECISION b(ncol,*),polynder 
!
!     ------------------------------------------------------
!     computes derivative matrix with fix for rounding error            
!     ------------------------------------------------------
!
      DO k = 1,nnew 
         sum = 0.d0 
         DO j = 1,no 
            IF ( j .EQ. k ) CYCLE 
            b(k,j) = polynder(j,xnew(k),no,xo) 
            IF ( dabs(b(k,j)) .LT. eps )      b(k,j) = 0.d0 
            sum = sum - b(k,j) 
            b(k,k) = sum 
         end do 
      END DO 
!
      RETURN 
      END SUBROUTINE dermat
!                                                                       
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the interpolation matrix.   
!!
!                                                                       
      SUBROUTINE intrpmat(ncol,no,xo,nnew,xnew,b) 
!                                                                       
!.......................................................................
!     date: 2/13/95                                                     
!     routines called: none                                             
!     includes: none                                                    
!     applicability:all                                                 
!                                                                       
!     compute the interpolation matrix.                                 
!.......................................................................
!                                                                       
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DIMENSION xo(*),xnew(*) 
      DOUBLE PRECISION b(ncol,*),polyn 
!                                                                       
      DO k = 1,nnew 
         DO j = 1,no 
            b(k,j) = polyn(j,xnew(k),no,xo)
         END DO
      END DO
!                                                                       
      RETURN 
      END SUBROUTINE intrpmat
