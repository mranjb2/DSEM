!> \file refine_geometry.f90
!! Routines for refining geometry in cases where we wish to switch polynomial order or Lobatto subdomains are small
!
!///////////////////////////////////////////////////////////////////////////////
!////////                                                               ////////
!////////       refine_geometry.f90                                     ////////
!////////                                                               ////////
!////////       contains:                                               ////////
!////////             SUBROUTINE set_no_points(d,mrtr) 			////////
!////////             SUBROUTINE set_geometry_refine(d,mrtr)            ////////
!////////             SUBROUTINE reset_pointers(d,mrtr)		        ////////
!////////                                                               ////////
!///////////////////////////////////////////////////////////////////////////////
!> @brief Adjust the number of Gauss points per domain.
!!
!! Define the length of each mortar
!
      SUBROUTINE set_no_points(d,mrtr,nmort)
!
!.......................................................................
!     adjust the # of Gauss/Gauss points per domain
!
!     DATE: 01/18/02

      USE domain_definition 
      USE mortar_definition
      USE mpi_par
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)


      ncgnew = resolution

      DO id = 1,ngridl 
        d(id)%ncg(:) = ncgnew
        d(id)%nc(:)  = ncgnew + 1 
      ENDDO

      DO j = 1,nmortl
        mrtr(j)%lenmortar(1) = ncgnew 
        mrtr(j)%lenmortar(2) = ncgnew 
        mrtr(j)%len(1,1)     = ncgnew 
        mrtr(j)%len(1,2)     = ncgnew 
        mrtr(j)%len(2,1)     = ncgnew 
        mrtr(j)%len(2,2)     = ncgnew 
      ENDDO

      DO j = 1,nmort
        nslavempi(j,2)       = ncgnew
        nslavempi(j,3)       = ncgnew
      ENDDO

      RETURN
     END SUBROUTINE set_no_points
      
!
!> @brief Refine the geometry
!!
!! This operation is only valid for straight sides. First jacobian and grid 
!! points are computed at cell centers. Then the quantities of lobatto-gauss-gauss 
!! face, gauss-lobatto-gauss face, and gauss-gauss-lobatto face are computed. Finally 
!! the mortar geometries of face 1, 2, 3, 4, 5, and 6 are computed.
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE set_geometry_refine(d,mrtr)
!
!.......................................................................
!     set metric terms (valid only for straight sides)
!
!     DATE: 01/18/02

      USE domain_definition 
      USE mortar_definition
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d
      TYPE (mortar) :: mrtr(nmp)

      DOUBLE PRECISION, DIMENSION(nmax,2) :: cxf,weight
      DOUBLE PRECISION , DIMENSION(3)     :: cx
      DOUBLE PRECISION , DIMENSION(3,3)   :: grad_x,grad_cx
      DOUBLE PRECISION                    :: J,sqr
      DOUBLE PRECISION, DIMENSION(6)      :: end_vals = (/0.0d0,1.0d0,0.0d0,1.0d0,1.0d0,0.0d0/)
      INTEGER                             :: isurf,iface,num_surf,ic,nxy,nYX
      INTEGER                             :: l,m,n,mrtr_no,mrtr_side
!

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
         mrtr_no   = nmortmpi(mrtr_no,1)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(2) = end_vals(iface)
         DO l = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(3) = cxf(l,2)
            DO n = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(1) = cxf(n,1)
               CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,n,l))
               CALL grad_hex8_map(d%corner,cx,grad_x)
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
         mrtr_no   = nmortmpi(mrtr_no,1)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(3) = end_vals(iface)
         DO m = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(2) = cxf(m,2)
            DO n = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(1) = cxf(n,1)
               CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,n,m))
               CALL grad_hex8_map(d%corner,cx,grad_x)
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
         mrtr_no   = nmortmpi(mrtr_no,1)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         CALL mortar_gauss_points(mrtr(mrtr_no)%lenmortar,cxf,weight)
         cx(1) = end_vals(iface)
         DO l = 1,mrtr(mrtr_no)%lenmortar(2)
            cx(3) = cxf(l,2)
            DO m = 1,mrtr(mrtr_no)%lenmortar(1)
               cx(2) = cxf(m,1)
               CALL hex8_map(d%corner,cx,mrtr(mrtr_no)%xg(:,m,l))
               CALL grad_hex8_map(d%corner,cx,grad_x)
               CALL Invert_Metric_Tensor(grad_x,grad_cx,J)
               sqr = SQRT(grad_cx(1,1)**2 + grad_cx(1,2)**2 + grad_cx(1,3)**2)
               mrtr(mrtr_no)%n_hat(1,m,l)  = grad_cx(1,1)/sqr
               mrtr(mrtr_no)%n_hat(2,m,l)  = grad_cx(1,2)/sqr
               mrtr(mrtr_no)%n_hat(3,m,l)  = grad_cx(1,3)/sqr
               mrtr(mrtr_no)%norm(m,l) = J*sqr
            END DO
         END DO
      END DO

      RETURN
     END SUBROUTINE set_geometry_refine

!
!> @brief Reset the pointers
!!
!! Reset pointers for each domain
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE reset_pointers(d,mrtr)
!
!.......................................................................
!     reset pointers
!
!     DATE: 01/18/02

      USE domain_definition 
      USE mortar_definition
      USE Order_Matrices
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)
!
!
!     --------------
!     reset pointers
!     --------------
!
      DO id = 1,ngridl
         DO k = 1,num_ordersx
            IF ( order_mapx(k) == d(id)%nc(1) )     THEN
               d(id)%dmx => dmx(:,:,k)
               d(id)%bx  => bx (:,:,k)
               EXIT
            END IF
         END DO
      END DO
      DO id = 1,ngridl
         DO k = 1,num_ordersy
            IF ( order_mapy(k) == d(id)%nc(2) )     THEN
               d(id)%dmy => dmy(:,:,k)
               d(id)%by  => by (:,:,k)
               EXIT
            END IF
         END DO
      END DO
      DO id = 1,ngridl
         DO k = 1,num_ordersz
            IF ( order_mapz(k) == d(id)%nc(3) )     THEN
               d(id)%dmz => dmz(:,:,k)
               d(id)%bz  => bz (:,:,k)
               EXIT
            END IF
         END DO
      END DO

      DO j = 1,nmortl
         CALL restart_mortar_matrices(mrtr(j))
      END DO


      RETURN
     END SUBROUTINE reset_pointers
