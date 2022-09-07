!>\file Prolong_gtol.f90
!! This file contains the subroutine to send EV from Gauss grid to Lobatto grid.

!> @brief This routines sends the filtered artificial viscosity calculated
!! on the Gauss-Gauss points to Lobatto points. The viscous flux
!! then will be calculated on Lobatto points.
      SUBROUTINE EV_Prolong(d,ngrid)      
!
      USE size
      USE domain_definition
      USE mpi
!
      IMPLICIT NONE
!
      TYPE (domain)   :: d(ngp)
!
      INTEGER         :: ngrid,id,m,n,l
! 
      
      DO id= 1,ngrid
!
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               CALL interp(nx,d(id)%bx,                            &
                              d(id)%muhg  (:,m,l)   ,d(id)%ncg(1), &
                              d(id)%muhlgg(:,m,l)   ,d(id)%nc(1)) 
            END DO
         END DO
!
         DO l = 1,d(id)%ncg(3)
            DO n = 1,d(id)%ncg(1)
               CALL interp(ny,d(id)%by,                            &
                              d(id)%muhg  (n,:,l)   ,d(id)%ncg(2), &
                              d(id)%muhglg(n,:,l)   ,d(id)%nc(2))

            END DO
         END DO
!
         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,                            &
                              d(id)%muhg  (n,m,:)   ,d(id)%ncg(3), &
                              d(id)%muhggl(n,m,:)   ,d(id)%nc(3))
            END DO
         END DO         
!
      END DO 
!
      RETURN
!
      END SUBROUTINE EV_Prolong
!
!//////////////////////////////////////////////////////////////////////////////