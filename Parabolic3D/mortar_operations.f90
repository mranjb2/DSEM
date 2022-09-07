!> @file
!> This file cantains subroutines that compute either interface conditions or boundary conditions
!> along the mortars, and send the mortar fluxes to the domain sides.
!///////////////////////////////////////////////////////////////////////////////////////////////
!////////										////////
!////////	mortar_operations.f90							////////
!////////										////////
!////////	contains:								////////
!////////										////////
!////////	      SUBROUTINE mortar_fluxes(mrtr,d,time)				////////
!////////	      SUBROUTINE mortar_fluxes_mpi(mrtr,d,time)				////////
!////////	      SUBROUTINE interface_flux(mrtr,flux,property1,property2,num_props)////////
!////////	      SUBROUTINE map_to_mortar_space(iside,mrtr)			////////
!////////	      SUBROUTINE map_to_face_space(iside,mrtr,mortar_flux,flux)		////////
!////////	      SUBROUTINE rotate_to_mortar(mrtr)					////////
!////////	      SUBROUTINE rotate_flux_to_face(mrtr,flux)				////////
!////////	      SUBROUTINE write_matrix(nc,u)					////////
!////////	      SUBROUTINE V_Avg_Soln(mrtr,Qbound)                                //////// 
!////////										////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the mortar fluxes and sends the mortar fluxes to the domain faces.
!!
      SUBROUTINE mortar_fluxes_mpi(mrtr,d,time,imort)
!
!......................................................................
!     Date: 11/24/98
!
!     Compute either interface conditions or boundary conditions
!     along the mortars, and send the mortar fluxes to the domain
!     sides. The interface flux routine is found here. The boundary
!     flux routine is found in the problem file.
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE Material_Properties
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
 
      idm    = nmortmpi(imort,2)
      nprocm = ngridmpi(idm,2)

      IF (myid == nprocm) THEN
         id    = mrtr%id(1)
         idml  = ngridmpi(id,1)
         iface = mrtr%iface(1)
         IF (d(idml)%bcond(iface) == 'periods' .or. &
             d(idml)%bcond(iface) == 'periodm')     THEN     ! this is a boundary mortar
    
            CALL boundary_flux(time,mrtr,d(idml)%bcond(iface),&
                            0.0d0,0,flux)
         ELSE
         ! write(*,*) 'chick mortar_Qm',mrtr%Q(2,2,:,1)
         ! write(*,*) 'chick mortar_Qs',mrtr%Q(2,2,:,2)
         ! write(*,*) 'chick id', mrtr%id(1),mrtr%id(2),imort,idm
            CALL interface_flux(mrtr,flux, &
                 0.0d0,&
                 0.0d0,0)
         ENDIF
      ENDIF
!
!     ------------------------------------------
!     send the mortar fluxes to the domain faces
!     ------------------------------------------
!
      CALL send_flux_to_faces(mrtr,d,flux,imort)
!
      RETURN
      END SUBROUTINE mortar_fluxes_mpi
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the mortar fluxes and sends the mortar fluxes to the domain faces.
!!
!
!
      SUBROUTINE mortar_fluxes(mrtr,d,time,imort)
!
!......................................................................
!     Date: 11/24/98
!
!     Compute either interface conditions or boundary conditions
!     along the mortars, and send the mortar fluxes to the domain
!     sides. The interface flux routine is found here. The boundary
!     flux routine is found in the problem file.
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE Material_Properties
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      id    = mrtr%id(1)
      idml  = ngridmpi(id,1)
      iface = mrtr%iface(1)
      IF ( mrtr%id(2) == 0 .or. d(idml)%bcond(iface) == 'periods' .or. &
           d(idml)%bcond(iface) == 'periodm')     THEN     ! this is a boundary mortar
    
         CALL boundary_flux(time,mrtr,d(idml)%bcond(iface),&
                            0.0d0,0,flux)

      ELSE                               ! this is an interface mortar

         CALL interface_flux(mrtr,flux, &
              0.0d0,&
              0.0d0,0)

      END IF
!
!     ------------------------------------------
!     send the mortar fluxes to the domain faces
!     ------------------------------------------
!

      CALL send_flux_to_faces(mrtr,d,flux,imort)
!
      RETURN
      END SUBROUTINE mortar_fluxes
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the flux along the face using the riemann solver.
!!
!!     Both conforming and nonconforming interfaces are handled here.
!!     The procedure here is first to rotate the face data to orient 
!!     itself with the mortar. Next the data is projected if 
!!     necessary to match the mortar space dimensions. The fluxes
!!     are computed with the riemann solver after that. The fluxes are then
!!     projected back to the face space dimensions. Finally, the fluxes
!!     are rotated back to match the slave face.
!
      SUBROUTINE interface_flux(mrtr,flux,property1,property2,num_props)
!
!......................................................................
!     Date: 11/24/98
!     
!     Compute the flux along the face using the riemann solver
!
!     Both conforming and nonconforming interfaces are handled here.
!     The procedure here is first to rotate the face data to orient 
!     itself with the mortar. Next the data is projected if 
!     necessary to match the mortar space dimensions. The fluxes
!     are computed with the riemann solver after that. The fluxes are then
!     projected back to the face space dimensions. Finally, the fluxes
!     are rotated back to match the slave face.
!......................................................................
!
      USE mortar_definition
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(neq+mxprops)   :: Ql,Qr
      DOUBLE PRECISION, DIMENSION(neq)           :: f
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1,property2
!
!     ------------------------------------------------------
!     Rotate and project the solutions onto the mortar space
!     ------------------------------------------------------
!
      IF ( mrtr%orient /= 0 )     CALL rotate_to_mortar(mrtr)
      DO iside = 1,2
         IF ( .NOT.mrtr%conforming(1,iside)   .OR. &
              .NOT.mrtr%conforming(2,iside) )      &
            CALL map_to_mortar_space(iside,mrtr)
      END DO
!
!     ---------------------------------------
!     Solve the Riemann problem on the mortar
!     ---------------------------------------
!
      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO np = 1,num_props  ! add material properties to end of vector
               Ql(neq+np) = property2(np)
               Qr(neq+np) = property1(np)
            END DO

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                 DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,2)
                     Qr(nv) = mrtr%Q(i,j,nv,1)
                  END DO

                  CALL riemann(mrtr%n_hat(:,i,j),mrtr%norm(i,j),Ql,Qr,f)

                  DO nv = 1,neq
                     flux(i,j,nv) = f(nv)
                  END DO
	       END DO
            END DO 

         CASE (2,4,5)

            DO np = 1,num_props  ! add material properties to end of vector
               Ql(neq+np) = property1(np)
               Qr(neq+np) = property2(np)
            END DO

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                 DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO

                  CALL riemann(mrtr%n_hat(:,i,j),mrtr%norm(i,j),Ql,Qr,f)

                  DO nv = 1,neq
                     flux(i,j,nv) = f(nv)
                  END DO
	       END DO
            END DO 

      END SELECT
!
      RETURN
      END SUBROUTINE interface_flux
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine interpolates the face solutions onto the mortar space.
!!
!!     There are three cases to consider: conforming in one of the two directions
!!     or conforming in no direction at all. Iside defines either the
!!     master (1) or slave (2) side of the mortar that is to be
!!     interpolated.
!
      SUBROUTINE map_to_mortar_space(iside,mrtr)
!
!......................................................................
!     Interpolate the face solutions onto the mortar space. There are 
!     three cases to consider: conforming in one of the two directions
!     or conforming in no direction at all. Iside defines either the
!     master (1) or slave (2) side of the mortar that is to be
!     interpolated.
!......................................................................
!
      USE mortar_definition
      USE Edge_Mappings
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: temp
      DOUBLE PRECISION, DIMENSION(nmax)          :: Qtemp,vtemp
      DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: prolong1,prolong2
!
      IF ( iside == 1 )     THEN
         prolong1 => mrtr%prolong1_dir1
         prolong2 => mrtr%prolong1_dir2
         idir = 1
         jdir = 2
      ELSE
         prolong1 => mrtr%prolong2_dir1
         prolong2 => mrtr%prolong2_dir2
         SELECT CASE (mrtr%orient)
            CASE (DFLT,B1F2,B1B2,F1B2)
               idir = 1
               jdir = 2
            CASE (F2F1,B2F1,B2B1,F2B1)
               idir = 2
               jdir = 1
         END SELECT
      END IF
!
!     -------------------------------------------------------------
!     At least one direction is non-conforming. Check to see if one
!     turns out to be conforming.
!     -------------------------------------------------------------
!
      IF      ( mrtr%conforming(1,iside) ) THEN ! direction 1 is conforming

         DO nv = 1,neq
            DO i = 1,mrtr%lenmortar(1)
               DO j = 1,mrtr%len(jdir,iside)
                  Qtemp(j) = mrtr%Q(i,j,nv,iside)
               END DO
               CALL interp(nmax,prolong2,Qtemp,mrtr%len(jdir,iside),&
                           vtemp,mrtr%lenmortar(2))
               DO j = 1,mrtr%lenmortar(2)
                  mrtr%Q(i,j,nv,iside) = vtemp(j)
               END DO
            END DO
         END DO
         
      ELSE IF ( mrtr%conforming(2,iside) ) THEN ! direction 2 is conforming

         DO nv = 1,neq
            DO j = 1,mrtr%lenmortar(2)
               CALL interp(nmax,prolong1,mrtr%Q(:,j,nv,iside),mrtr%len(idir,iside),&
                           vtemp,mrtr%lenmortar(1))
               DO i = 1,mrtr%lenmortar(1)
                  mrtr%Q(i,j,nv,iside) = vtemp(i)
               END DO
            END DO
         END DO
         
      ELSE ! fully two direction interpolation needed

         DO nv = 1,neq
            DO j = 1,mrtr%len(jdir,iside)
               CALL interp(nmax,prolong1,mrtr%Q(:,j,nv,iside),mrtr%len(idir,iside),&
                           vtemp,mrtr%lenmortar(1))
               DO i = 1,mrtr%lenmortar(1)
                  temp(i,j) = vtemp(i)
               END DO
            END DO
            DO i = 1,mrtr%lenmortar(1)
               DO j = 1,mrtr%len(jdir,iside)
                  Qtemp(j) = temp(i,j)
               END DO
               CALL interp(nmax,prolong2,Qtemp,mrtr%len(jdir,iside),&
                           vtemp,mrtr%lenmortar(2))
               DO j = 1,mrtr%lenmortar(2)
                  mrtr%Q(i,j,nv,iside) = vtemp(j)
               END DO
            END DO
         END DO

      END IF
!
      RETURN
!
      END SUBROUTINE map_to_mortar_space
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine maps the mortar flux onto the face space.
!!
!!  There are three cases to consider: conforming in one of the two directions
!!  or conforming in no direction at all.
!
      SUBROUTINE map_to_face_space(iside,mrtr,mortar_flux,flux)
!
!......................................................................
!     Map the mortar flux onto the face space. There are three
!     cases to consider: conforming in one of the two directions
!     or conforming in no direction at all.
!......................................................................
!
      USE mortar_definition
      USE Edge_Mappings
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: temp
      DOUBLE PRECISION, DIMENSION(nmax)          :: Ftemp,vtemp
      DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: restrict1,restrict2
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: mortar_flux,flux
!
      IF ( iside == 1 )     THEN
         restrict1 => mrtr%restrict1_dir1
         restrict2 => mrtr%restrict1_dir2
         idir = 1
         jdir = 2
      ELSE
         restrict1 => mrtr%restrict2_dir1
         restrict2 => mrtr%restrict2_dir2
         SELECT CASE (mrtr%orient)
            CASE (DFLT,B1F2,B1B2,F1B2)
               idir = 1
               jdir = 2
            CASE (F2F1,B2F1,B2B1,F2B1)
               idir = 2
               jdir = 1
         END SELECT
      END IF
!
!     -------------------------------------------------------------
!     At least one direction is non-conforming. Check to see if one
!     turns out to be conforming.
!     -------------------------------------------------------------
!
      IF      ( mrtr%conforming(1,iside) ) THEN ! direction 1 is conforming

         DO nv = 1,neq
            DO i = 1,mrtr%lenmortar(1)
               DO j = 1,mrtr%lenmortar(2)
                  Ftemp(j) = mortar_flux(i,j,nv)
               END DO
               CALL interp(nmax,restrict2,Ftemp,mrtr%lenmortar(2),&
                           vtemp,mrtr%len(jdir,iside))
               DO j = 1,mrtr%len(jdir,iside)
                  flux(i,j,nv) = vtemp(j)
               END DO
            END DO
         END DO
         
      ELSE IF ( mrtr%conforming(2,iside) ) THEN ! direction 2 is conforming

         DO nv = 1,neq
            DO j = 1,mrtr%lenmortar(2)
               CALL interp(nmax,restrict1,mortar_flux(:,j,nv),mrtr%lenmortar(1),&
                           vtemp,mrtr%len(idir,iside))
               DO i = 1,mrtr%len(idir,iside)
                  flux(i,j,nv) = vtemp(i)
               END DO
            END DO
         END DO
         
      ELSE ! fully two direction interpolation needed

         DO nv = 1,neq
            DO j = 1,mrtr%lenmortar(2)
               CALL interp(nmax,restrict1,mortar_flux(:,j,nv),mrtr%lenmortar(1),&
                           vtemp,mrtr%len(idir,iside))
               DO i = 1,mrtr%len(idir,iside)
                  temp(i,j) = vtemp(i)
               END DO
            END DO
            DO i = 1,mrtr%len(idir,iside)
               DO j = 1,mrtr%lenmortar(2)
                  Ftemp(j) = temp(i,j)
               END DO
               CALL interp(nmax,restrict2,Ftemp,mrtr%lenmortar(2),&
                           vtemp,mrtr%len(jdir,iside))
               DO j = 1,mrtr%len(jdir,iside)
                  flux(i,j,nv) = vtemp(j)
               END DO
            END DO
         END DO

      END IF
!
      RETURN
!
      END SUBROUTINE map_to_face_space
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine rotates the data to align it along the mortar.
!!
!!     If a slave face is rotated with respect to a mortar, then this
!!     routine is called to rotate the data to align it along the mortar.
!!     It is assumed that copying the data and writing it back onto the
!!     mortar is more efficient than doing loop index modifications on the
!!     fly when setting up the left/right sides before solving the riemann
!!     problem. If this is not true, then this routine can be eliminated, 
!!     and the mapping can be done point by point.
!
      SUBROUTINE rotate_to_mortar(mrtr)
!
!......................................................................
!     If a slave face is rotated with respect to a mortar, then this
!     routine is called to rotate the data to align it along the mortar.
!     It is assumed that copying the data and writing it back onto the
!     mortar is more efficient than doing loop index modifications on the
!     fly when setting up the left/right sides before solving the riemann
!     problem. If this is not true, then this routine can be eliminated, 
!     and the mapping can be done point by point.
!......................................................................
!
      USE mortar_definition
      USE Edge_Mappings
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: temp
!
      SELECT CASE ( mrtr%orient )

         CASE ( F1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(i,l) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(1,2)
                  DO i = 1,mrtr%len(2,2)
                     mrtr%Q(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( B1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(k,l) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     mrtr%Q(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE( B1F2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(k,j) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(1,2)
                  DO i = 1,mrtr%len(2,2)
                     mrtr%Q(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( F2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     temp(j,i) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%Q(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( B2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(k,i) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%Q(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( B2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(k,l) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%Q(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( F2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(j,l) = mrtr%Q(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%Q(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

      END SELECT
!
      RETURN
      END SUBROUTINE rotate_to_mortar
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine rotates the flux to align it along the face.
!!
!!     If an element face is rotated with respect to a mortar, then this
!!     routine is called to rotate the flux to align it along the face.
!!     It is assumed that copying the data and writing it back onto the
!!     mortar is efficient.
!
      SUBROUTINE rotate_flux_to_face(mrtr,flux)
!
!......................................................................
!     If an element face is rotated with respect to a mortar, then this
!     routine is called to rotate the flux to align it along the face.
!     It is assumed that copying the data and writing it back onto the
!     mortar is efficient.
!......................................................................
!
      USE mortar_definition
      USE Edge_Mappings
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: temp
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
!
      SELECT CASE ( mrtr%orient )

         CASE ( F1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(i,j) = flux(i,l,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( B1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(i,j) = flux(k,l,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE( B1F2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(i,j) = flux(k,j,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( F2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     temp(i,j) = flux(j,i,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( B2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(i,j) = flux(k,i,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( B2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(i,j) = flux(k,l,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( F2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(i,j) = flux(j,l,nv)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     flux(i,j,nv) = temp(i,j)
                  END DO
               END DO
            END DO

      END SELECT
!
      RETURN
      END SUBROUTINE rotate_flux_to_face
!
      SUBROUTINE write_matrix(nc,u)
      use size
      integer,dimension(2) :: nc
      double precision, dimension(nmax,nmax) :: u

      write(6,*) '**'
      do j = 1,nc(2)
         write(6,*) (u(i,j),i=1,nc(1))
      end do
      write(6,*) '**'

      return
      end subroutine write_matrix

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine averages the master and slave solution and set 
!> the master and slave Q to this average.
!!
!
      SUBROUTINE V_Avg_Soln(mrtr,Qbound)
!
!......................................................................
!     Date: 03/14/01
!
!     Average the master and slave solution and set the master and
!     slave Q to this average.
!  
!     applicability: Navier-Stokes
!......................................................................

      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Q
!-----
!
!
!     ------------------------------------------------------
!     Rotate and project the solutions onto the mortar space
!     ------------------------------------------------------
!

!      IF ( mrtr%orient /= 0 )     CALL rotate_to_mortar(mrtr)
!     Rotation is not needed, it is already performed in the Euler code

      DO iside = 1,2
         IF ( .NOT.mrtr%conforming(1,iside)   .OR. &
              .NOT.mrtr%conforming(2,iside) )      &
            CALL map_to_mortar_space(iside,mrtr)
      END DO
!
!
!     ---------------------------------------
!     Average the solution on the mortar
!     ---------------------------------------
!

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
                     Qm     = mrtr%Q(i,j,nv,1)
                     Qs     = mrtr%Q(i,j,nv,2)
                     Q(nv)  = (Qm + Qs)*0.5d0
                  END DO

                  mrtr%Q(i,j,1:neq,2)  = Q(1:neq)
                  mrtr%Q(i,j,1:neq,1)  = Q(1:neq)
                  Qbound(i,j,1:neq)    = Q(1:neq)

               END DO
            END DO

      RETURN
      END SUBROUTINE V_Avg_Soln
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine rotates the viscid fluxes data to align it along the mortar.
!!
!!     If a slave face is rotated with respect to a mortar, then this
!!     routine is called to rotate the data to align it along the mortar.
!!     It is assumed that copying the data and writing it back onto the
!!     mortar is more efficient than doing loop index modifications on the
!!     fly when setting up the left/right sides before solving the riemann
!!     problem. If this is not true, then this routine can be eliminated, 
!!     and the mapping can be done point by point.
!
      SUBROUTINE rotate_fv_to_mortar(mrtr)
!
!......................................................................
!     If a slave face is rotated with respect to a mortar, then this
!     routine is called to rotate the data to align it along the mortar.
!     It is assumed that copying the data and writing it back onto the
!     mortar is more efficient than doing loop index modifications on the
!     fly when setting up the left/right sides before solving the riemann
!     problem. If this is not true, then this routine can be eliminated, 
!     and the mapping can be done point by point.
!......................................................................
!
      USE mortar_definition
      USE Edge_Mappings
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: temp
!
      SELECT CASE ( mrtr%orient )

         CASE ( F1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(i,l) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(1,2)
                  DO i = 1,mrtr%len(2,2)
                     mrtr%fv(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( B1B2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  l = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(k,l) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     mrtr%fv(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE( B1F2 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     k = mrtr%len(1,2) - i + 1
                     temp(k,j) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO j = 1,mrtr%len(1,2)
                  DO i = 1,mrtr%len(2,2)
                     mrtr%fv(i,j,nv,2) = temp(i,j)
                  END DO
               END DO
            END DO

         CASE ( F2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     temp(j,i) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%fv(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( B2F1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     temp(k,i) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%fv(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( B2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  k = mrtr%len(2,2) - j + 1
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(k,l) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%fv(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

         CASE ( F2B1 )

            DO nv = 1,neq
               DO j = 1,mrtr%len(2,2)
                  DO i = 1,mrtr%len(1,2)
                     l = mrtr%len(1,2) - i + 1
                     temp(j,l) = mrtr%fv(i,j,nv,2)
                  END DO
               END DO

               DO i = 1,mrtr%len(1,2)
                  DO j = 1,mrtr%len(2,2)
                     mrtr%fv(j,i,nv,2) = temp(j,i)
                  END DO
               END DO
            END DO

      END SELECT
!
      RETURN
      END SUBROUTINE rotate_fv_to_mortar
!
