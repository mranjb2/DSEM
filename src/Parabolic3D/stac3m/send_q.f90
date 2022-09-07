!> @file
!> This file cantains subroutines that send the q to faces and mortars
!///////////////////////////////////////////////////////////////////////////////
!////////								////////
!////////	send_q.f90						////////
!////////								////////
!////////	contains:						////////
!////////								////////
!////////	      SUBROUTINE send_q_to_faces(mrtr,d,mortar_q)	////////
!////////	      SUBROUTINE send_q_to_faces_mpi(mrtr,d,mortar_q,imort)	////////
!////////	      SUBROUTINE send_q_to_mortars(mrtr,d)      	////////
!////////	      SUBROUTINE send_q_to_mortars_mpi(mrtr,d,imort)	////////
!////////								////////
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends q to the boundary face.
!!
!     Once the mortar fluxes are computed, send them back to the 
!     subdomains
!
      SUBROUTINE send_q_to_faces(mrtr,d,flux,imort)
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
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_q,flux

      idm    = mrtr%id(1)
      ids    = mrtr%id(2)
      nprocm = ngridmpi(idm,2)
      idml   = ngridmpi(idm,1)
      IF (ids /= 0) THEN
        nprocs = ngridmpi(ids,2)
        idsl   = ngridmpi(ids,1)
      ENDIF
                                                                              
      IF (mrtr%id(2) ==0) THEN  ! boundary face
!
!     --------------------------------
!     send the flux to the boundary face
!     --------------------------------
!
      mortar_side = 1
      iface = mrtr%iface(mortar_side)
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,1,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)
            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,d(idml)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(d(idml)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,d(idml)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(1,i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

      END SELECT
      
      ELSE
!
!     --------------------------------
!     send the flux to the master face
!     --------------------------------
!
      mortar_side = 1
      iface = mrtr%iface(mortar_side)

!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,1,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,d(idml)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(d(idml)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,d(idml)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(1,i,j,nv) = flux(i,j,nv)
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

      IF (nmortedge(imort,1) ==1) RETURN

      iface = mrtr%iface(mortar_side)


      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qglg(i,1,j,nv) = flux(i,j,nv)
                  END DO
              END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qglg(i,d(idsl)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qggl(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qlgg(d(idsl)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qggl(i,j,d(idsl)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idsl)%Qlgg(1,i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         END SELECT
!
       END IF
      RETURN
      END SUBROUTINE send_q_to_faces
!
!//////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends q to the boundary face in different cores.
!!
!
      SUBROUTINE send_q_to_faces_mpi(mrtr,d,flux,imort)
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
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_q,flux

!
!     --------------------------------
!     send the flux to the boundary face
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
!
      iface = mrtr%iface(mortar_side)
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,1,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)
            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qglg(i,d(idml)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(d(idml)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qggl(i,j,d(idml)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     d(idml)%Qlgg(1,i,j,nv) = flux(i,j,nv)
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
!
      iface = mrtr%iface(mortar_side)
!
      IF ( mrtr%orient /= 0 )     CALL rotate_flux_to_face(mrtr,flux)

      CALL send_q_mpi(nprocs,imort,flux(:,:,:))                                

    ELSEIF (myid == nprocs) THEN


      mortar_side = 2
      iface = nslavempi(imort,1)
      nlen1 = nslavempi(imort,2)
      nlen2 = nslavempi(imort,3)
      IF ( ids == 0 .OR. &
           d(idsl)%bcond(iface) == 'periodm' .OR. &
           d(idsl)%bcond(iface) == 'periods' )     RETURN  ! this is a boundary mortar

      CALL recv_q_mpi(nprocm,imort,flux(:,:,:))                           

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qglg(i,1,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qglg(i,d(idsl)%nc(2),j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qggl(i,j,1,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qlgg(d(idsl)%nc(1),i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qggl(i,j,d(idsl)%nc(3),nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     d(idsl)%Qlgg(1,i,j,nv) = flux(i,j,nv)
                  END DO
               END DO
            END DO

         END SELECT
!
       END IF
      RETURN
      END SUBROUTINE send_q_to_faces_mpi
!//////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends the q to the mortars.
!!
!  Once the fluxes are interpolated, send them to the mortars
! 
      SUBROUTINE send_q_to_mortars(mrtr,d,imort)
!
!......................................................................
!     date: 12/11/01
!
!     stac3m version
!
!     Once the fluxes are interpolated, send them to the mortars
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
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux

      idm    = nmortmpi(imort,2)
      ids    = nmortmpi(imort,3) 
      nprocm = ngridmpi(idm,2)
      idml   = ngridmpi(idm,1)
      IF (ids /= 0 ) THEN
        nprocs = ngridmpi(ids,2)
        idsl   = ngridmpi(ids,1)
      ENDIF

      IF (mrtr%id(2) ==0) THEN  ! boundary face
!
!     --------------------------------
!     send the flux to the boundary face
!     --------------------------------
!
      mortar_side = 1
      iface = mrtr%iface(mortar_side)
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qglg(i,1,j,nv)  
                  END DO
               END DO
            END DO

         CASE (2)
            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qglg(i,d(idml)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qggl(i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qlgg(d(idml)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qggl(i,j,d(idml)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qlgg(1,i,j,nv)  
                  END DO
               END DO
            END DO
      END SELECT

      DO nv = 1,neq
         DO j = 1,mrtr%len(2,mortar_side)
            DO i = 1,mrtr%len(1,mortar_side)
                 mrtr%Q(i,j,nv,1) = flux(i,j,nv) 
            END DO
         END DO
      END DO

      
      ELSE
!
!     --------------------------------
!     send the flux to the master face
!     --------------------------------
!
      mortar_side = 1
      iface = mrtr%iface(mortar_side)
!
      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qglg(i,1,j,nv)  
                  END DO
               END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qglg(i,d(idml)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qggl(i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qlgg(d(idml)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qggl(i,j,d(idml)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idml)%Qlgg(1,i,j,nv)  
                  END DO
               END DO
            END DO
      END SELECT

       DO nv = 1,neq
          DO j = 1,mrtr%len(2,mortar_side)
             DO i = 1,mrtr%len(1,mortar_side)
                 mrtr%Q(i,j,nv,1) = flux(i,j,nv) 
             END DO
          END DO
       END DO

!
!     ----------------------------------------------------------
!     Map and rotate the mortar fluxes onto the slave spaces
!     ----------------------------------------------------------
!
      mortar_side = 2
     
      IF ( ids == 0 )     RETURN  ! this is a boundary mortar
      iface = mrtr%iface(mortar_side)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qglg(i,1,j,nv)  
                  END DO
              END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qglg(i,d(idsl)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qggl(i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qlgg(d(idsl)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qggl(i,j,d(idsl)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,mortar_side)
                  DO i = 1,mrtr%len(1,mortar_side)
                     flux(i,j,nv) = d(idsl)%Qlgg(1,i,j,nv)  
                  END DO
               END DO
            END DO
         END SELECT

      DO nv = 1,neq
         DO j = 1,mrtr%len(2,mortar_side)
            DO i = 1,mrtr%len(1,mortar_side)
                 mrtr%Q(i,j,nv,2) = flux(i,j,nv)
            END DO
         END DO
      END DO

!
       END IF
      RETURN
      END SUBROUTINE send_q_to_mortars
!
!//////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sends the q to the mortars in different cores.
!!
!
      SUBROUTINE send_q_to_mortars_mpi(mrtr,d,imort)
!
!......................................................................
!     date: 12/11/01
!
!     stac3m version
!
!     Once the fluxes are interpolated, send them to the mortars
!
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux


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

!
    IF (myid == nprocs) THEN
!
! determine  Q on slave mortar
!
      mortar_side = 2

      iface = nslavempi(imort,1)
      nlen1 = nslavempi(imort,2)
      nlen2 = nslavempi(imort,3)

      SELECT CASE (iface)

         CASE (1)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qglg(i,1,j,nv)  
                  END DO
              END DO
            END DO

         CASE (2)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qglg(i,d(idsl)%nc(2),j,nv)  
                  END DO
               END DO
            END DO

         CASE (3)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qggl(i,j,1,nv)  
                  END DO
               END DO
            END DO

         CASE (4)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qlgg(d(idsl)%nc(1),i,j,nv)  
                  END DO
               END DO
            END DO

         CASE (5)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qggl(i,j,d(idsl)%nc(3),nv)  
                  END DO
               END DO
            END DO

         CASE (6)

            DO nv = 1,neq
               DO j = 1,nlen2
                  DO i = 1,nlen1
                     mrtr%Q(i,j,nv,1) = d(idsl)%Qlgg(1,i,j,nv)  
                  END DO
               END DO
            END DO
         END SELECT

       CALL MPI_ISEND (                                          &
            mrtr%Q(1,1,1,1),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocm,imort,comm1d,req(nreq),ierr)

    ELSEIF (myid == nprocm) THEN

       CALL MPI_IRECV (                                          &
            mrtr%Q(1,1,1,2),nmax*nmax*neq, MPI_DOUBLE_PRECISION,           &
            nprocs,imort,comm1d,req(nreq),ierr)

    ENDIF

!
      RETURN
      END SUBROUTINE send_q_to_mortars_mpi
