!> @file
!> This file cantains subroutines that send the fluxes to faces and mortars
!///////////////////////////////////////////////////////////////////////////////
!////////								////////
!////////	send_flux.f90						////////
!////////								////////
!////////	contains:						////////
!////////								////////
!////////	      SUBROUTINE send_flux_to_faces(mrtr,d,mortar_flux)	
!////////             SUBROUTINE send_flux_to_mortars(mrtr,r,s,t,imort,d)
!////////             SUBROUTINE send_flux_to_mortars_mpi(mrtr,r,s,t,imort,d)
!////////	      SUBROUTINE send_flux_to_faces_and_add(mrtr,d,mortar_flux)	
!////////	      SUBROUTINE send_flux_to_faces_and_add_mpi(mrtr,d,mortar_flux,imort)	
!////////								////////
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends fluxes to the face.
!!
!!     Once the mortar fluxes are computed, send them back to the 
!!     subdomains
!
      SUBROUTINE send_flux_to_faces(mrtr,d,mortar_flux,imort)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
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

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    flux(i,j,nv) = mortar_flux(i,j,nv)
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,1,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,d(idml)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(d(idml)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,d(idml)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(1,i,j,nv) = flux(i,j,nv)
                  END DO
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
           d(idml)%bcond(ifacem) == 'periods' )     RETURN  ! this is a boundary mortar
      IF (nmortedge(imort,1) == 1) RETURN
      iface = mrtr%iface(mortar_side)

      IF( mrtr%conforming(1,2) .AND. &
          mrtr%conforming(2,2) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%lenmortar(2)
                 DO i = 1,mrtr%lenmortar(1)
                    flux(i,j,nv) = mortar_flux(i,j,nv)
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF

      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%g(i,1,j,nv) = mrtr%sign*flux(i,j,nv)
                     END DO
                  END DO
               END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%g(i,d(idsl)%nc(2),j,nv) = mrtr%sign*flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%h(i,j,1,nv) = mrtr%sign*flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%f(d(idsl)%nc(1),i,j,nv) = mrtr%sign*flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%h(i,j,d(idsl)%nc(3),nv) = mrtr%sign*flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%f(1,i,j,nv) = mrtr%sign*flux(i,j,nv)
                  END DO
               END DO
            END DO

         END SELECT
!
      RETURN
      END SUBROUTINE send_flux_to_faces
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends fluxes to the masters.
!!
!
      SUBROUTINE send_flux_to_mortars(mrtr,imort,d)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
      !DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t
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

!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idml,i,1,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idml,i,d(idml)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idml,i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idml,d(idml)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idml,i,j,d(idml)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idml,1,i,j,nv)  
                  END DO
               END DO
            END DO

      END SELECT

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    mrtr%fv(i,j,nv,1) = flux(i,j,nv)/mrtr%norm(i,j) 
                   
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF
!
!     ----------------------------------------------------------
!     Map and rotate the mortar fluxes onto the slave spaces
!     ----------------------------------------------------------
!
      mortar_side = 2

      iface = mrtr%iface(mortar_side)
      IF ( ids == 0 ) RETURN  ! this is a boundary mortar


      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idsl,i,1,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idsl,i,d(idsl)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idsl,i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idsl,d(idsl)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idsl,i,j,d(idsl)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idsl,1,i,j,nv)  
                  END DO
               END DO
            END DO

      END SELECT

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    mrtr%fv(i,j,nv,2) = flux(i,j,nv)/mrtr%norm(i,j) 
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF

!
      RETURN
      END SUBROUTINE send_flux_to_mortars
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends fluxes to the mortars in different cores.
!!
!
      SUBROUTINE send_flux_to_mortars_mpi(mrtr,imort,d)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
      USE User_Data
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
      !DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t
      INTEGER, DIMENSION(0:3,2)                  :: orient_map = &
                                                 RESHAPE((/1,2,1,2,2,1,2,1/),(/4,2/))
!
!     --------------------------------
!     send the flux to the master face
!     --------------------------------
!
    idm    = nmortmpi(imort,2)
    ids    = nmortmpi(imort,3)
    nprocm = ngridmpi(idm,2)
    nprocs = ngridmpi(ids,2)
    idml   = ngridmpi(idm,1)
    idsl   = ngridmpi(ids,1)
!
!     ----------------------------------------------------------
!     Map and rotate the mortar fluxes onto the slave spaces
!     ----------------------------------------------------------
!
     IF ( myid == nprocm) THEN
      mortar_side = 1
      iface = mrtr%iface(mortar_side)

!
      SELECT CASE (iface)
         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idml,i,1,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = s(idml,i,d(idml)%nc(2),j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idml,i,j,1,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idml,d(idml)%nc(1),i,j,nv)
                  END DO
               END DO
            END DO
         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = t(idml,i,j,d(idml)%nc(3),nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = r(idml,1,i,j,nv)
                  END DO
               END DO
            END DO

      END SELECT

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    mrtr%fv(i,j,nv,1) = flux(i,j,nv)/mrtr%norm(i,j)
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF
   END IF

   IF (myid == nprocs) THEN
      mortar_side = 2

      iface = nslavempi(imort,1)
      nlen1 = nslavempi(imort,2)
      nlen2 = nslavempi(imort,3)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = s(idsl,i,1,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = s(idsl,i,d(idsl)%nc(2),j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = t(idsl,i,j,1,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = r(idsl,d(idsl)%nc(1),i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = t(idsl,i,j,d(idsl)%nc(3),nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%fv(i,j,nv,1) = r(idsl,1,i,j,nv)
                  END DO
               END DO
            END DO

      END SELECT

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    mrtr%fv(i,j,nv,1) = mrtr%fv(i,j,nv,1)/mrtr%norm(i,j)   
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF

       CALL MPI_ISEND (                                          &
            mrtr%fv(1,1,1,1),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocm,imort,comm1d,req(nreq),ierr)

     ! CALL send_q_mpi(nprocm,imort,mrtr%fv(:,:,:,1))
     ! CALL send_q_mpi(nprocm,imort,flux(:,:,:))

    ELSEIF (myid == nprocm) THEN

       CALL MPI_IRECV (                                          &
            mrtr%fv(1,1,1,2),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocs,imort,comm1d,req(nreq),ierr)
       !CALL recv_q_mpi(nprocs,imort,mrtr%fv(:,:,:,2))

    END IF

!
      RETURN
      END SUBROUTINE send_flux_to_mortars_mpi
!
!//////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine adds fluxes in the face.
!!
!
      SUBROUTINE send_flux_to_faces_add(mrtr,d,mortar_flux,imort)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
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

      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    flux(i,j,nv) = mortar_flux(i,j,nv)
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF
    IF (myid ==1 .and. ids ==0) THEN
     ! write(*,*) 'chick idml,iface',idm,idml,iface
     ! write(*,*) 'chick flux',flux(2,2,:)
    ENDIF
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,1,j,nv) = d(idml)%g(i,1,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,d(idml)%nc(2),j,nv) = d(idml)%g(i,d(idml)%nc(2),j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,1,nv) = d(idml)%h(i,j,1,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(d(idml)%nc(1),i,j,nv) = d(idml)%f(d(idml)%nc(1),i,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,d(idml)%nc(3),nv) = d(idml)%h(i,j,d(idml)%nc(3),nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(1,i,j,nv) = d(idml)%f(1,i,j,nv) - flux(i,j,nv)/re
                  END DO
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
           d(idml)%bcond(ifacem) == 'periods' )     RETURN  ! this is a boundary mortar


      iface = mrtr%iface(mortar_side)
      IF( mrtr%conforming(1,2) .AND. &
          mrtr%conforming(2,2) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%lenmortar(2)
                 DO i = 1,mrtr%lenmortar(1)
                   IF ( mrtr%sign == -1) THEN
             	     flux(i,j,nv) = -mortar_flux(i,j,nv)
            	   ELSE
                     flux(i,j,nv) = mortar_flux(i,j,nv)
                   ENDIF
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF

      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%g(i,1,j,nv) = d(idsl)%g(i,1,j,nv) - flux(i,j,nv)/re
                     END DO
                  END DO
               END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%g(i,d(idsl)%nc(2),j,nv) = d(idsl)%g(i,d(idsl)%nc(2),j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%h(i,j,1,nv) = d(idsl)%h(i,j,1,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%f(d(idsl)%nc(1),i,j,nv) = d(idsl)%f(d(idsl)%nc(1),i,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%h(i,j,d(idsl)%nc(3),nv) = d(idsl)%h(i,j,d(idsl)%nc(3),nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%f(1,i,j,nv) = d(idsl)%f(1,i,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         END SELECT
!
      RETURN
      END SUBROUTINE send_flux_to_faces_add
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine adds fluxes in the face in different cores.
!!
!
!
      SUBROUTINE send_flux_to_faces_add_mpi(mrtr,d,mortar_flux,imort)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
      INTEGER, DIMENSION(0:3,2)                  :: orient_map = &
                                                 RESHAPE((/1,2,1,2,2,1,2,1/),(/4,2/))
!
!     --------------------------------
!     send the flux to the master face
!     --------------------------------
!
    mortar_side = 1
    idm    = nmortmpi(imort,2)
    ids    = nmortmpi(imort,3)
    nprocm = ngridmpi(idm,2)
    nprocs = ngridmpi(ids,2)
    idml   = ngridmpi(idm,1)
    idsl   = ngridmpi(ids,1)

    IF (myid == nprocm) THEN
      iface = mrtr%iface(mortar_side)


      IF( mrtr%conforming(1,1) .AND.           &
          mrtr%conforming(2,1) )          THEN

           DO nv = 1,neq
              DO j = 1,mrtr%len(2,mortar_side)
                 DO i = 1,mrtr%len(1,mortar_side)
                    flux(i,j,nv) = mortar_flux(i,j,nv)
                 END DO
              END DO
           END DO

      ELSE
         CALL map_to_face_space(mortar_side,mrtr,mortar_flux,flux)
      END IF
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,1,j,nv) = d(idml)%g(i,1,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%g(i,d(idml)%nc(2),j,nv) = d(idml)%g(i,d(idml)%nc(2),j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,1,nv) = d(idml)%h(i,j,1,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(d(idml)%nc(1),i,j,nv) = d(idml)%f(d(idml)%nc(1),i,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%h(i,j,d(idml)%nc(3),nv) = d(idml)%h(i,j,d(idml)%nc(3),nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%f(1,i,j,nv) = d(idml)%f(1,i,j,nv) - flux(i,j,nv)/re
                  END DO
               END DO
            END DO

      END SELECT
    ENDIF
 

!
!     ----------------------------------------------------------
!     Map and rotate the mortar fluxes onto the slave spaces
!     ----------------------------------------------------------
!
     
    IF (myid == nprocm) THEN
      mortar_side = 2
      ifacem= mrtr%iface(1)
      IF ( ids == 0 .OR. &
           d(idml)%bcond(ifacem) == 'periodm' .OR. &
           d(idml)%bcond(ifacem) == 'periods' )     RETURN  ! this is a boundary mortar

      DO nv = 1,neq
        DO j = 1,mrtr%lenmortar(2)
          DO i = 1,mrtr%lenmortar(1)
            IF ( mrtr%sign == -1) THEN
             mrtr%fv(i,j,nv,2) = -mortar_flux(i,j,nv)
            ELSE
             mrtr%fv(i,j,nv,2) = mortar_flux(i,j,nv)
            ENDIF
          END DO
        END DO
      END DO

      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

       CALL MPI_ISEND (                                          &
            mrtr%fv(1,1,1,2),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocs,imort,comm1d,req(nreq),ierr)

      ! CALL send_q_mpi(nprocs,imort,mrtr%fv(:,:,:,2))

    ELSEIF (myid == nprocs) THEN

      mortar_side = 2
      iface = nslavempi(imort,1)

      IF ( ids == 0 .OR. &
           d(idsl)%bcond(iface) == 'periodm' .OR. &
           d(idsl)%bcond(iface) == 'periods' )     RETURN  ! this is a boundary mortar

       CALL MPI_IRECV (                                          &
            mrtr%fv(1,1,1,1),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocm,imort,comm1d,req(nreq),ierr)
      ! CALL recv_q_mpi(nprocm,imort,mrtr%fv(:,:,:,1))
    ENDIF
!
      RETURN
      END SUBROUTINE send_flux_to_faces_add_mpi
!
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine adds fluxes over different cores.
!!
!
!
      SUBROUTINE add_flux_mpi(mrtr,d,imort)
!
!......................................................................
!     date: 10/24/02
!
!     stac3m version
!
!     add the fluxes that are sent with non-blocking 
!     message passing
!
!......................................................................

      USE domain_definition
      USE mortar_definition
      USE physics
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
      INTEGER, DIMENSION(0:3,2)                  :: orient_map = &
                                                 RESHAPE((/1,2,1,2,2,1,2,1/),(/4,2/))
      idm    = nmortmpi(imort,2)
      ids    = nmortmpi(imort,3)
      nprocm = ngridmpi(idm,2)
      nprocs = ngridmpi(ids,2)
      idml   = ngridmpi(idm,1)
      idsl   = ngridmpi(ids,1)
      iface = nslavempi(imort,1)
      nlen1 = nslavempi(imort,2)
      nlen2 = nslavempi(imort,3)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%g(i,1,j,nv) = d(idsl)%g(i,1,j,nv) - mrtr%fv(i,j,nv,1)/re
                     END DO
                  END DO
               END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%g(i,d(idsl)%nc(2),j,nv) = d(idsl)%g(i,d(idsl)%nc(2),j,nv) - mrtr%fv(i,j,nv,1)/re
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%h(i,j,1,nv) = d(idsl)%h(i,j,1,nv) - mrtr%fv(i,j,nv,1)/re
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%f(d(idsl)%nc(1),i,j,nv) = d(idsl)%f(d(idsl)%nc(1),i,j,nv) - mrtr%fv(i,j,nv,1)/re
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%h(i,j,d(idsl)%nc(3),nv) = d(idsl)%h(i,j,d(idsl)%nc(3),nv) - mrtr%fv(i,j,nv,1)/re
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%f(1,i,j,nv) = d(idsl)%f(1,i,j,nv) - mrtr%fv(i,j,nv,1)/re
                  END DO
               END DO
            END DO

         END SELECT
      RETURN
      END SUBROUTINE add_flux_mpi
