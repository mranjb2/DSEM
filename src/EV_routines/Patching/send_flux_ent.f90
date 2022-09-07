!> \file
!! EV flux routines
!///////////////////////////////////////////////////////////////////////////////////////////////
!                                                                                       ////////
!                                                                                       ////////
!        SUBROUTINE send_flux_to_faces_ent(mrtr,d,imort)                                ////////
!                                                                                       ////////
!             Send calculated entropy flux on mortars to the side Gauss-Lobatto points  ////////
!                                                                                       ////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!> @brief Send calculated entropy flux on mortars to the side Gauss-Lobatto points.
      SUBROUTINE send_flux_to_faces_ent(mrtr,d,mortar_flux,imort)
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
!
      IMPLICIT NONE
!      
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
!
      INTEGER                                    :: idms,ids,nprocm,idm,idml,nprocs,idsl  
      INTEGER                                    :: mortar_side,iface,i,j,ifacem,imort    
      DOUBLE PRECISION, DIMENSION(nmax,nmax,2)   :: mortar_flux,flux
      INTEGER, DIMENSION(0:3,2)                  :: orient_map = &
                                                 RESHAPE((/1,2,1,2,2,1,2,1/),(/4,2/))
!
!     --------------------------------
!     send the flux to the master face
!     --------------------------------
!
      idm    = mrtr%id(1)
      ids    = mrtr%id(2)
      nprocm = ngridmpi(idm,2)
      idml   = ngridmpi(idm,1)
      IF (ids /= 0) THEN
        nprocs = ngridmpi(ids,2)
        idsl   = ngridmpi(ids,1)
      ENDIF
    

      mortar_side = 1
      iface = mrtr%iface(mortar_side)
           
         DO j = 1,mrtr%len(2,mortar_side)
            DO i = 1,mrtr%len(1,mortar_side)
               flux(i,j,:) = mortar_flux(i,j,:)
            END DO
         END DO       
! 

      SELECT CASE (iface) 
      
                                                                   
         CASE (1)
           
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%ge(i,1,j) = flux(i,j,1)
                     d(idml)%se(i,1,j) = flux(i,j,2)
                  END DO
               END DO           

         CASE (2)
           
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%ge(i,d(idml)%nc(2),j) = flux(i,j,1)
                     d(idml)%se(i,d(idml)%nc(2),j) = flux(i,j,2)
                  END DO
               END DO            

         CASE (3)
            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%he(i,j,1) = flux(i,j,1)
                     d(idml)%te(i,j,1) = flux(i,j,2)
                  END DO
               END DO           

         CASE (4)
            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%fe(d(idml)%nc(1),i,j) = flux(i,j,1)
                     d(idml)%re(d(idml)%nc(1),i,j) = flux(i,j,2)
                  END DO
               END DO           

         CASE (5)
          
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%he(i,j,d(idml)%nc(3)) = flux(i,j,1)
                     d(idml)%te(i,j,d(idml)%nc(3)) = flux(i,j,2)
                  END DO
               END DO           

         CASE (6)

               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%fe(1,i,j) = flux(i,j,1)
                     d(idml)%re(1,i,j) = flux(i,j,2)
                  END DO
               END DO

      END SELECT
!
!     ----------------------------------------------------------
!     Map and rotate the mortar fluxes onto the slave spaces
!     ----------------------------------------------------------
!
      mortar_side = 2

      ifacem= mrtr%iface(1)
      IF ( ids == 0 .OR. &
           d(idml)%bcond(ifacem) == 'periodm' .OR. &
           d(idml)%bcond(ifacem) == 'periods' )     RETURN  
      IF (nmortedge(imort,1) == 1) RETURN
      iface = mrtr%iface(mortar_side)

           
      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
               flux(i,j,:) = mortar_flux(i,j,:)
         END DO
      END DO
          

      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

      SELECT CASE (iface)

         CASE (1)

            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%ge(i,1,j) = mrtr%sign*flux(i,j,1)
                     d(idsl)%se(i,1,j) = mrtr%sign*flux(i,j,2)
                     END DO
                  END DO
               

         CASE (2)

            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%ge(i,d(idsl)%nc(2),j) = mrtr%sign*flux(i,j,1)
                     d(idsl)%se(i,d(idsl)%nc(2),j) = mrtr%sign*flux(i,j,2)
                  END DO
               END DO
            

         CASE (3)

           
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%he(i,j,1) = mrtr%sign*flux(i,j,1)
                     d(idsl)%te(i,j,1) = mrtr%sign*flux(i,j,2)
                  END DO
               END DO
            

         CASE (4)

            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%fe(d(idsl)%nc(1),i,j) = mrtr%sign*flux(i,j,1)
                     d(idsl)%re(d(idsl)%nc(1),i,j) = mrtr%sign*flux(i,j,2)
                  END DO
               END DO
            

         CASE (5)

            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%he(i,j,d(idsl)%nc(3)) = mrtr%sign*flux(i,j,1)
                     d(idsl)%te(i,j,d(idsl)%nc(3)) = mrtr%sign*flux(i,j,2)
                  END DO
               END DO
            

         CASE (6)

            
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%fe(1,i,j) = mrtr%sign*flux(i,j,1)
                     d(idsl)%re(1,i,j) = mrtr%sign*flux(i,j,2)
                  END DO
               END DO
           

         END SELECT
!
      RETURN
      END SUBROUTINE send_flux_to_faces_ent
!
