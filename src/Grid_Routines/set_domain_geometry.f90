!
!> @file
!! This file contains the subroutines that set up the geometry for the domain
!!
!///////////////////////////////////////////////////////////////////////
!> @brief
!! This subroutine sets up the grid and metric terms for each element/subdomain
      SUBROUTINE set_geometry(nodes,d,mrtr,num_surf,surface_array)
!
!.......................................................................
!     Set up the grid and metric terms for each element/subdomain
!     This version is for the STAC^3M codes, and compute the metric terms 
!     at the Lobatto/Gauss points.
!.......................................................................
!
      USE Domain_Definition
      USE Mortar_Definition
      USE FE_Data_Types

      IMPLICIT  none

      DOUBLE PRECISION, DIMENSION(nmax,2) :: cxf,weight
      DOUBLE PRECISION , DIMENSION(3)     :: cx
      DOUBLE PRECISION , DIMENSION(3,3)   :: grad_x,grad_cx
      DOUBLE PRECISION                    :: J,sqr
      DOUBLE PRECISION, DIMENSION(6)      :: end_vals = (/0.0d0,1.0d0,0.0d0,1.0d0,1.0d0,0.0d0/)
      INTEGER                             :: isurf,iface,num_surf,ic,nxy,nYX
      INTEGER                             :: l,m,n,mrtr_no,mrtr_side
      logical                             :: all_straight
!
      TYPE (node_vector)                  :: nodes
      TYPE(face_interp), DIMENSION(6)     :: face_data
      TYPE (surface), DIMENSION(max_surf) :: surface_array
      TYPE (domain)                       :: d
      TYPE (mortar),DIMENSION(nmp)        :: mrtr
!
!     ----------------------
!     set corner information
!     ----------------------
!
      DO ic = 1,8
         d%corner(:,ic) = nodes%node(d%node(ic))%x
      END DO
!
!     -----------------------
!     set boundary conditions
!     -----------------------
!
      DO iface = 1,6
         isurf = d%which_surface(iface)
         IF ( isurf /= 0 )     d%bcond(iface) = surface_array(isurf)%bcond
      END DO
!
!     ----------------------------------------------------------------------------
!>     Set the interpolation data for the faces. If the element is a hex8 and there
!>     is no modification to any side, then we can skip setting the faces and use
!>     the plain hex8 routines. If any side is modified, then the general routine
!>     must be used.
!     ----------------------------------------------------------------------------
!
      all_straight = .FALSE.
      IF( d%type == 8 )     THEN
         all_straight = .TRUE.
         DO iface = 1,6
            isurf = d%which_surface(iface)
            IF (isurf /= 0 )     THEN
               IF( surface_array(isurf)%surface_type == 'equation' ) then
                  all_straight = .FALSE.
                  EXIT
               END IF
            END IF
         END DO
      END IF
!
!     --------------------------------------------------------
!     compute the interior collocation points and metric terms
!     --------------------------------------------------------
!
!                                   ------------------
      IF( all_straight )     THEN ! all straight sides
!                                   ------------------
!
!        -----------------------------------------------------
!        jacobian and grid points are computed at cell centers
!        -----------------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL hex8_map(d%corner,cx,d%xg(:,n,m,l))
                  CALL grad_hex8_map(d%corner,cx,grad_x)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%jacob(n,m,l) = 1.d0/J
                  DO nYX = 1,3
                   DO nxy = 1,3
                     d%gmetg(nYX,nxy,n,m,l) = grad_cx(nYX,nxy) 
                   ENDDO
                  ENDDO
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!        lobatto-gauss-gauss face quantities (faces 4&6)
!        -----------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%nc(1)
                  cx(1) = d%cx(n,1)
                  CALL grad_hex8_map(d%corner,cx,grad_x)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(1,:,n,m,l) = J*grad_cx(1,:)
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!        gauss-lobatto-gauss face quantities (faces 1&2)
!        -----------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%nc(2)
               cx(2) = d%cx(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL grad_hex8_map(d%corner,cx,grad_x)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(2,:,n,m,l) = J*grad_cx(2,:)
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!       gauss-gauss-lobatto face quantities (faces 3&5)
!        -----------------------------------------------
!
         DO l = 1,d%nc(3)
            cx(3) = d%cx(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL grad_hex8_map(d%corner,cx,grad_x)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(3,:,n,m,l) = J*grad_cx(3,:)
               END DO
            END DO
         END DO
!
!            ------------------------
      ELSE ! at least one curved side
!            ------------------------

         DO iface = 1,6
            IF ( d%which_surface(iface) == 0 )     THEN
               CALL set_face_from_element(iface,d,nodes,face_data(iface))
            ELSE
               isurf = d%which_surface(iface)
               IF (surface_array(isurf)%surface_type == 'equation')     THEN
                  CALL set_face_from_surface(iface,d,surface_array(isurf),face_data(iface))
               ELSE
                  CALL set_face_from_element(iface,d,nodes,face_data(iface))
               END IF
            END IF
         END DO
!
!        -----------------------------------------------------
!        jacobian and grid points are computed at cell centers
!        -----------------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL general_hex_map(cx,d%xg(:,n,m,l),d%corner,face_data)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%jacob(n,m,l) = 1.0d0/J
                  DO nYX = 1,3
                   DO nxy = 1,3
                     d%gmetg(nYX,nxy,n,m,l) = grad_cx(nYX,nxy) 
                   ENDDO
                  ENDDO
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!        lobatto-gauss-gauss face quantities (faces 4&6)
!        -----------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%nc(1)
                  cx(1) = d%cx(n,1)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(1,:,n,m,l) = J*grad_cx(1,:)
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!        gauss-lobatto-gauss face quantities (faces 1&2)
!        -----------------------------------------------
!
         DO l = 1,d%ncg(3)
            cx(3) = d%cxg(l,3)
            DO m = 1,d%nc(2)
               cx(2) = d%cx(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(2,:,n,m,l) = J*grad_cx(2,:)
               END DO
            END DO
         END DO
!
!        -----------------------------------------------
!       gauss-gauss-lobatto face quantities (faces 3&5)
!        -----------------------------------------------
!
         DO l = 1,d%nc(3)
            cx(3) = d%cx(l,3)
            DO m = 1,d%ncg(2)
               cx(2) = d%cxg(m,2)
               DO n = 1,d%ncg(1)
                  cx(1) = d%cxg(n,1)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
                  CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
                  d%gmet(3,:,n,m,l) = J*grad_cx(3,:)
               END DO
            END DO
         END DO

      END IF
!
!     -----------------------
!     compute mortar geometry
!     -----------------------
!
!     ----------
!     face 1 & 2
!     ----------
!
      DO iface = 1,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(2) = end_vals(iface)
         DO l = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(3) = cxf(l,2)
            DO n = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(1) = cxf(n,1)
               IF ( all_straight )     THEN
                  CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,n,l))
                  CALL grad_hex8_map(d%corner,cx,grad_x)
               ELSE
                  CALL general_hex_map(cx,mrtr(mrtr_no)%xg(:,n,l),d%corner,face_data)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
               END IF
               CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
               sqr = SQRT(grad_cx(2,1)**2 + grad_cx(2,2)**2 + grad_cx(2,3)**2)
               mrtr(mrtr_no)%n_hat(1,n,l)  = grad_cx(2,1)/sqr
               mrtr(mrtr_no)%n_hat(2,n,l)  = grad_cx(2,2)/sqr
               mrtr(mrtr_no)%n_hat(3,n,l)  = grad_cx(2,3)/sqr
               mrtr(mrtr_no)%norm(n,l) = J*sqr
            END DO
         END DO
      END DO
!
!     ----------
!     face 3 & 5
!     ----------
!
      DO iface = 3,5,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(3) = end_vals(iface)
         DO m = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(2) = cxf(m,2)
            DO n = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(1) = cxf(n,1)
               IF ( all_straight )     THEN
                  CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,n,m))
                  CALL grad_hex8_map(d%corner,cx,grad_x)
               ELSE
                  CALL general_hex_map(cx,mrtr(mrtr_no)%xg(:,n,m),d%corner,face_data)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
               END IF
               CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
               sqr = SQRT(grad_cx(3,1)**2 + grad_cx(3,2)**2 + grad_cx(3,3)**2)
               mrtr(mrtr_no)%n_hat(1,n,m)  = grad_cx(3,1)/sqr
               mrtr(mrtr_no)%n_hat(2,n,m)  = grad_cx(3,2)/sqr
               mrtr(mrtr_no)%n_hat(3,n,m)  = grad_cx(3,3)/sqr
               mrtr(mrtr_no)%norm(n,m) = J*sqr
            END DO
         END DO
      END DO
!
!     ----------
!     face 4 & 6
!     ----------
!
      DO iface = 4,6,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(1) = end_vals(iface)
         DO l = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(3) = cxf(l,2)
            DO m = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(2) = cxf(m,1)
               IF ( all_straight )     THEN
                  CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,m,l))
                  CALL grad_hex8_map(d%corner,cx,grad_x)
               ELSE
                  CALL general_hex_map(cx,mrtr(mrtr_no)%xg(:,m,l),d%corner,face_data)
                  CALL general_hex_grad(cx,grad_x,d%corner,face_data)
               END IF
               CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
               sqr = SQRT(grad_cx(1,1)**2 + grad_cx(1,2)**2 + grad_cx(1,3)**2)
               mrtr(mrtr_no)%n_hat(1,m,l)  = grad_cx(1,1)/sqr
               mrtr(mrtr_no)%n_hat(2,m,l)  = grad_cx(1,2)/sqr
               mrtr(mrtr_no)%n_hat(3,m,l)  = grad_cx(1,3)/sqr
               mrtr(mrtr_no)%norm(m,l) = J*sqr
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE set_geometry
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!! This subroutine sets up the interpolation data that defines a face. In this case,
!! we take the information from the element point definition.
      SUBROUTINE set_face_from_element(iface,d,nodes,face_data)
!
!.......................................................................
!     Set up the interpolation data that defines a face. In this case,
!     we take the information from the element point definition.
!.......................................................................
!
      USE FE_Data_Types
      USE Domain_Definition
      USE Element_Library

      TYPE (domain)      :: d !< element
      TYPE (face_interp) :: face_data  !< face data
      TYPE (node_vector) :: nodes !< node vector type
!
      SELECT CASE (d%type)

         CASE(8)
            face_data%num_knots = (/2,2/)
            face_data%knots(1:2,1) = cx_hex8
            face_data%knots(1:2,2) = cx_hex8
            DO j = 1,2
               DO i = 1,2
                  inode = d%node(Intrp_Map_Hex8(i,j,iface))
                  face_data%points(:,i,j) = nodes%node(inode)%x
               END DO
            END DO

         CASE(26)
            face_data%num_knots = (/3,3/)
            face_data%knots(1:3,1) = cx_hex26
            face_data%knots(1:3,2) = cx_hex26
            DO j = 1,3
               DO i = 1,3
                  inode = d%node(Intrp_Map_Hex26(i,j,iface))
                  face_data%points(:,i,j) = nodes%node(inode)%x
               END DO
            END DO

         CASE(64)
            !face_data%num_knots = (/4,4/)
            !face_data%knots(1:4,1) = cx_hex64
            !face_data%knots(1:4,2) = cx_hex64
            !DO j = 1,4
            !   DO i = 1,4
            !      inode = d%node(Intrp_Map_Hex64(i,j,iface))
            !      face_data%points(:,i,j) = nodes%node(inode)%x
            !   END DO
            !END DO

         CASE DEFAULT
            WRITE(6,*) 'Unsupported element type :',d%type
      END SELECT
!
      RETURN
      END SUBROUTINE set_face_from_element
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!! This subroutine sets up the faces on the surface.
!!
!!     if the face is attached to a curved surface, evaluate the interpolation
!!    points from that surface. This requires the finding of where the face
!!     attaches to the surface. It is assumed here that the surface on which the
!!     face lies is known, and it is just a matter of finding the parametric
!!     coordinates of the face on the surface. Note, the first four nodes of a
!!     face are the same for all the supported element types.
!
      SUBROUTINE set_face_from_surface(face_id,d,the_surface,face_data)
!
!.......................................................................
!     if the face is attached to a curved surface, evaluate the interpolation
!     points from that surface. This requires the finding of where the face
!     attaches to the surface. It is assumed here that the surface on which the
!     face lies is known, and it is just a matter of finding the parametric
!     coordinates of the face on the surface. Note, the first four nodes of a
!     face are the same for all the supported element types.

!     @@@@@ not implemented @@@@@@@@@
!.......................................................................
!
      USE FE_Data_Types
      USE Domain_Definition
      USE Element_Library

      TYPE (domain)     :: d
      TYPE(face_interp) :: face_data
      TYPE (surface)    :: the_surface

      DOUBLE PRECISION, DIMENSION(3,4) :: face_point
      INTEGER                          :: the_node,face_id
!
!     ------------------------------------------------------------------------
!     get the four nodes for this face (hex8 is sufficient to get the corners)
!     ------------------------------------------------------------------------
!
      DO k = 1,4
         the_node = d%node(Face_Map_hex8(k,face_id))
         face_point(:,k) = d%corner(:,the_node)
      END DO
!
!     ------------------------------------------------------------------------
!     find the face-parameter location of the fout face nodes
!     ------------------------------------------------------------------------
!
      
      
!
      RETURN
      END SUBROUTINE set_face_from_surface
