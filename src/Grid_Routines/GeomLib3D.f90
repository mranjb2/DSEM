!> \file GeomLib3D.f90 
!! Mapping library for finite element geometries
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
!////////												////////
!////////	geomlib3d.f90										////////
!////////												////////
!////////	contains:										////////
!////////												////////
!////////	      SUBROUTINE quad4_map(corner,cx,x)							////////
!////////	      SUBROUTINE hex8_map(corner,cx,x)							////////
!////////	      SUBROUTINE grad_hex8_map(corner,cx,grad_x)					////////
!////////	      SUBROUTINE general_hex_map(cx,x,corner_point,face_data)				////////
!////////	      SUBROUTINE compute_hex_map(cx,x,face,edge,corner)					////////
!////////	      SUBROUTINE compute_face_point(face_data,cx,p)					////////
!////////	      SUBROUTINE compute_Poly_2D(face_data,cx,p)					////////
!////////	      SUBROUTINE general_hex_grad(cx,grad_x,corner_point,face_data)			////////
!////////	      SUBROUTINE compute_grad_hex_map(cx,grad_x,face,face_der,edge,edge_der,corner)	////////
!////////	      SUBROUTINE compute_face_derivative(face_data,cx,grad)				////////
!////////	      SUBROUTINE compute_Poly_2D_der(face_data,cx,grad)					////////
!////////	      SUBROUTINE Invert_Metric_Tensor(grad_x,grad_cx,J)					////////
!////////												////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
!> @brief Compute the transfinite mapping for a quad4 plate element
!!
!!     corner = the 4 corners in standard fe format
!!
!!     cx     = computational space variable (X,Y)
!!
!!     x      = physical space variable (x,y)
      SUBROUTINE quad4_map(corner,cx,x)
!
!
      IMPLICIT none
!
      DOUBLE PRECISION, DIMENSION(3,4),INTENT(IN) :: corner !< four corners of element in standard FE format
      DOUBLE PRECISION, DIMENSION(2),INTENT(IN)   :: cx !< computational space variable (X,Y)
      DOUBLE PRECISION, DIMENSION(3),INTENT(OUT)  :: x !< physical space variable (x,y)
!
      INTEGER                                     :: j
!
      DO j = 1,3
         x(j)  =   corner(j,1)*(1.d0 - cx(1))*(1.d0 - cx(2)) &
                 + corner(j,2)* cx(1)        *(1.d0 - cx(2)) &
                 + corner(j,3)* cx(1)        * cx(2)         &
                 + corner(j,4)*(1.d0 - cx(1))* cx(2)
      END DO
!
      RETURN

      END SUBROUTINE quad4_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the transfinite mapping for a hex8 element
!!
!!     corner = the 8 corners in standard fe format
!!
!!     cx     = computational space variable (X,Y,Z)
!!
!!     x      = physical space variable (x,y,z)
      SUBROUTINE hex8_map(corner,cx,x)
!
!
      IMPLICIT none
!
      DOUBLE PRECISION, DIMENSION(3,8),INTENT(IN) :: corner !< 8 corners in standard FE format
      DOUBLE PRECISION, DIMENSION(3),INTENT(IN)   :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION, DIMENSION(3),INTENT(OUT)  :: x !< physical space variable (x,y,z)
!
      INTEGER                                     :: j
!
      DO j = 1,3
         x(j)  =   corner(j,1)*(1.d0 - cx(1))*(1.d0 - cx(2))*(1.d0 - cx(3)) &
                 + corner(j,2)* cx(1)        *(1.d0 - cx(2))*(1.d0 - cx(3)) &
                 + corner(j,3)* cx(1)        * cx(2)        *(1.d0 - cx(3)) &
                 + corner(j,4)*(1.d0 - cx(1))* cx(2)        *(1.d0 - cx(3)) &
                 + corner(j,5)*(1.d0 - cx(1))*(1.d0 - cx(2))* cx(3)         &
                 + corner(j,6)* cx(1)        *(1.d0 - cx(2))* cx(3)         &
                 + corner(j,7)* cx(1)        * cx(2)        * cx(3)         &
                 + corner(j,8)*(1.d0 - cx(1))* cx(2)        * cx(3)
      END DO
!
      RETURN

      END SUBROUTINE hex8_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the gradient of the transfinite mapping for a hex8 element.
!> The gradient is stored as d x_i/d X_j
!!
!!     corner = the 8 corners in standard fe format
!!
!!     cx     = computational space variable (X,Y,Z)
!!
!!     grad_x = physical space gradient
!!
!!              grad_x(:,1) = (x_X,y_X,z_X), etc.
!
!
      SUBROUTINE grad_hex8_map(corner,cx,grad_x)
!
!
      IMPLICIT none
!
      DOUBLE PRECISION, DIMENSION(3,8),INTENT(IN)  :: corner !< 8 corners in standard fe format
      DOUBLE PRECISION, DIMENSION(3),INTENT(IN)    :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION, DIMENSION(3,3),INTENT(OUT) :: grad_x !< physical space gradient
      INTEGER                                      :: i
!
      DO i = 1,3
         grad_x(i,1) =  -corner(i,1)*(1.d0 - cx(2))*(1.d0 - cx(3)) &
                       + corner(i,2)*(1.d0 - cx(2))*(1.d0 - cx(3)) &
                       + corner(i,3)* cx(2)        *(1.d0 - cx(3)) &
                       - corner(i,4)* cx(2)        *(1.d0 - cx(3)) &
                       - corner(i,5)*(1.d0 - cx(2))* cx(3)         &
                       + corner(i,6)*(1.d0 - cx(2))* cx(3)         &
                       + corner(i,7)* cx(2)        * cx(3)         &
                       - corner(i,8)* cx(2)        * cx(3)
!
         grad_x(i,2) =  -corner(i,1)*(1.d0 - cx(1))*(1.d0 - cx(3)) &
                       - corner(i,2)* cx(1)        *(1.d0 - cx(3)) &
                       + corner(i,3)* cx(1)        *(1.d0 - cx(3)) &
                       + corner(i,4)*(1.d0 - cx(1))*(1.d0 - cx(3)) &
                       - corner(i,5)*(1.d0 - cx(1))* cx(3)         &
                       - corner(i,6)* cx(1)        * cx(3)         &
                       + corner(i,7)* cx(1)        * cx(3)         &
                       + corner(i,8)*(1.d0 - cx(1))* cx(3)
!
         grad_x(i,3) =  -corner(i,1)*(1.d0 - cx(1))*(1.d0 - cx(2)) &
                       - corner(i,2)* cx(1)        *(1.d0 - cx(2)) &
                       - corner(i,3)* cx(1)        * cx(2)         &
                       - corner(i,4)*(1.d0 - cx(1))* cx(2)         &
                       + corner(i,5)*(1.d0 - cx(1))*(1.d0 - cx(2)) &
                       + corner(i,6)* cx(1)        *(1.d0 - cx(2)) &
                       + corner(i,7)* cx(1)        * cx(2)         &
                       + corner(i,8)*(1.d0 - cx(1))* cx(2)
      END DO
!
      RETURN

      END SUBROUTINE grad_hex8_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Given the six faces and four corners of a hex element, and a
!>   computational point cx, return the physical space location x
!!
!
      SUBROUTINE general_hex_map(cx,x,corner_point,face_data)
!
!
      USE FE_Data_Types

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION , DIMENSION(3)   :: cx !< computational space variable (X,Y,Z)
	  DOUBLE PRECISION , DIMENSION(3)   :: x !< physical space variable (x,y,z)
      DOUBLE PRECISION , DIMENSION(3,6) :: face_point !< face points
      DOUBLE PRECISION , DIMENSION(3,12):: edge_point !< edge coordinates
      DOUBLE PRECISION , DIMENSION(3,8) :: corner_point !< corner coordinates
      TYPE(face_interp), DIMENSION(6)   :: face_data   !< face strings
!
!     compute edges
!
      CALL compute_face_point(face_data(1),(/cx(1),0.d0/) ,edge_point(:,1))
      CALL compute_face_point(face_data(1),(/1.0d0,cx(3)/),edge_point(:,2))
      CALL compute_face_point(face_data(1),(/cx(1),1.0d0/),edge_point(:,3))
      CALL compute_face_point(face_data(1),(/0.0d0,cx(3)/),edge_point(:,4))
      CALL compute_face_point(face_data(2),(/cx(1),0.d0/) ,edge_point(:,5))
      CALL compute_face_point(face_data(2),(/1.0d0,cx(3)/),edge_point(:,6))
      CALL compute_face_point(face_data(2),(/cx(1),1.0d0/),edge_point(:,7))
      CALL compute_face_point(face_data(2),(/0.0d0,cx(3)/),edge_point(:,8))
      CALL compute_face_point(face_data(4),(/cx(2),0.d0/) ,edge_point(:,10))
      CALL compute_face_point(face_data(6),(/cx(2),0.d0/) ,edge_point(:,9))
      CALL compute_face_point(face_data(4),(/cx(2),1.0d0/),edge_point(:,11))
      CALL compute_face_point(face_data(6),(/cx(2),1.0d0/),edge_point(:,12))
!
!     compute faces
!
      CALL compute_face_point(face_data(1),(/cx(1),cx(3)/),face_point(:,1))
      CALL compute_face_point(face_data(2),(/cx(1),cx(3)/),face_point(:,2))
      CALL compute_face_point(face_data(3),(/cx(1),cx(2)/),face_point(:,3))
      CALL compute_face_point(face_data(4),(/cx(2),cx(3)/),face_point(:,4))
      CALL compute_face_point(face_data(5),(/cx(1),cx(2)/),face_point(:,5))
      CALL compute_face_point(face_data(6),(/cx(2),cx(3)/),face_point(:,6))
!
!     compute the mapping
!
      CALL compute_hex_map(cx,x,face_point,edge_point,corner_point)
!
      RETURN
      END SUBROUTINE general_hex_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the transfinite mapping for a general hex element
!! 
!! cx          = computational space variable (X,Y,Z)
!!
!! x           = resultant physical space variable (x,y,z)
!!
!! face        = face interpolant location for appropriate local face coordinate
!!
!! edge        = edge interpolant location for appropriate edge coordinate
!!
!! corner      = The 8 corners in standard fe format
!
      SUBROUTINE compute_hex_map(cx,x,face,edge,corner)
!
!
      IMPLICIT none
!
!     input...
!
      DOUBLE PRECISION, DIMENSION(3)   ,INTENT(IN) :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION, DIMENSION(3,6) ,INTENT(IN) :: face !< face interpolant location for appropriate local face coordinate
      DOUBLE PRECISION, DIMENSION(3,12),INTENT(IN) :: edge !< edge interpolant location for appropriate edge coordinate
      DOUBLE PRECISION, DIMENSION(3,8) ,INTENT(IN) :: corner !< 8 corners in standard fe format
!
!     output...
!
      DOUBLE PRECISION, DIMENSION(3),INTENT(OUT)  :: x !< resultant physical space variable (x,y,z)
!
!     local...
!
      INTEGER                                     :: j,iface
!-----
      DO j = 1,3
!
!        face contributions
!

         x(j) =  face(j,6)*(1.d0 - cx(1)) + face(j,4)*cx(1) &
               + face(j,1)*(1.d0 - cx(2)) + face(j,2)*cx(2) &
               + face(j,3)*(1.d0 - cx(3)) + face(j,5)*cx(3)
!
!        edge contributions
!
         x(j) = x(j) - edge(j,1) *(1.d0 - cx(2))*(1.d0 - cx(3)) &
                     - edge(j,3) *(1.d0 - cx(2))*        cx(3)  &
                     - edge(j,5) *        cx(2) *(1.d0 - cx(3)) &
                     - edge(j,7) *        cx(2) *        cx(3)  &
                     - edge(j,9) *(1.d0 - cx(1))*(1.d0 - cx(3)) &
                     - edge(j,12)*(1.d0 - cx(1))*        cx(3)  &
                     - edge(j,10)*(1.d0 - cx(3))*        cx(1)  &
                     - edge(j,11)*        cx(1) *        cx(3)  &
                     - edge(j,4) *(1.d0 - cx(1))*(1.d0 - cx(2)) &
                     - edge(j,8) *(1.d0 - cx(1))*        cx(2)  &
                     - edge(j,2) *        cx(1) *(1.d0 - cx(2)) &
                     - edge(j,6) *        cx(1) *        cx(2)
!
!        corner contributions
!
         x(j) = x(j) + corner(j,1)*(1.d0 -  cx(1))*(1.d0 -  cx(2))*(1.d0 -  cx(3)) &
                     + corner(j,5)*(1.d0 -  cx(1))*(1.d0 -  cx(2))*         cx(3)  &
                     + corner(j,4)*(1.d0 -  cx(1))*         cx(2) *(1.d0 -  cx(3)) &
                     + corner(j,8)*(1.d0 -  cx(1))*         cx(2) *         cx(3)  &
                     + corner(j,2)*         cx(1) *(1.d0 -  cx(2))*(1.d0 -  cx(3)) &
                     + corner(j,6)*         cx(1) *(1.d0 -  cx(2))*         cx(3)  &
                     + corner(j,3)*         cx(1) *         cx(2) *(1.d0 -  cx(3)) &
                     + corner(j,7)*         cx(1) *         cx(2) *         cx(3)
      END DO
!
      RETURN
      END SUBROUTINE compute_hex_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the nodes on a subdomain face according by 
!> interpolating the face_data
!!
!! Use bi-linear mapping if it is a flat side (for speed), 
!! otherwise use general lagrange interpolation
      SUBROUTINE compute_face_point(face_data,cx,p)
!
!
      USE FE_Data_Types

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE(face_interp)                :: face_data
      DOUBLE PRECISION , DIMENSION(2)  :: cx !< computational space (X,Y)
      DOUBLE PRECISION , DIMENSION(3)  :: p !< polynomial
!
!     use bi-linear mapping of this is a flat side (for speed) otherwise
!     do general lagrange interpolation
!
      IF( face_data%num_knots(1) == 2 .AND. face_data%num_knots(2) == 2)     THEN
         DO j = 1,3
            p(j)  =   face_data%points(j,1,1)*(1.d0 - cx(1))*(1.d0 - cx(2)) &
                    + face_data%points(j,2,1)* cx(1)        *(1.d0 - cx(2)) &
                    + face_data%points(j,2,2)* cx(1)        * cx(2)         &
                    + face_data%points(j,1,2)*(1.d0 - cx(1))* cx(2)
         END DO
      ELSE
         CALL compute_Poly_2D(face_data,cx,p)
      END IF
      
      RETURN
      END SUBROUTINE compute_face_point
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the lagrange interpolant through the points p(i,j)
!> at the location (x,y)
!!
!
      SUBROUTINE compute_Poly_2D(face_data,cx,p)
!
!
      USE FE_Data_Types

      IMPLICIT none
!
      TYPE(face_interp)                 :: face_data !< face data
      DOUBLE PRECISION, DIMENSION(2)    :: cx !< computational space variable
      DOUBLE PRECISION, DIMENSION(3)    :: p !< polynomial
      DOUBLE PRECISION, DIMENSION(nmax) :: l_i
      DOUBLE PRECISION                  :: l_j
      DOUBLE PRECISION                  :: polyn
      INTEGER                           :: i,j,k
      EXTERNAL polyn
!
      DO i = 1,face_data%num_knots(1)
         l_i(i) = polyn(i,cx(1),face_data%num_knots(1),face_data%knots(:,1))
      END DO
!
      p = 0.0d0
      DO j = 1,face_data%num_knots(2)
         l_j = polyn(j,cx(2),face_data%num_knots(2),face_data%knots(:,2))
         DO i = 1,face_data%num_knots(1)
            DO k = 1,3
               p(k) = p(k) + face_data%points(k,i,j)*l_i(i)*l_j
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE compute_Poly_2D
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Given the six faces and four corners of a hex element, and a
!> computational point cx, return the jacobian matrix grad_x = d x_i/d X_j
!!
!!     cx           = computational space variable (X,Y,Z)
!!
!!     grad_x       = physical space gradient
!!
!!                    grad_x(:,1) = (x_X,y_X,z_X), etc.
!!
!!     corner_point = the 8 corners in standard fe format
!!
!!     face_data    = interpolant data for the 6 faces
!
      SUBROUTINE general_hex_grad(cx,grad_x,corner_point,face_data)
!
!
      USE FE_Data_Types

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION , DIMENSION(3)      :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION , DIMENSION(3,3)    :: grad_x !< gradient of space variable (x_X,y_X,z_X) in i,j,k
      DOUBLE PRECISION , DIMENSION(3,2)    :: grad_2d !< two dimenional gradient
      DOUBLE PRECISION , DIMENSION(3,6)    :: face_point
      DOUBLE PRECISION , DIMENSION(3,2,6)  :: face_der
      DOUBLE PRECISION , DIMENSION(3,12)   :: edge_point
      DOUBLE PRECISION , DIMENSION(3,12)   :: edge_der
      DOUBLE PRECISION , DIMENSION(3,8)    :: corner_point
      TYPE(face_interp), DIMENSION(6)      :: face_data !< interpolant data for the 6 faces
      DOUBLE PRECISION                     :: eps
!
!     compute edges
!
      CALL compute_face_point(face_data(1),(/cx(1),0.d0/) ,edge_point(:,1))
      CALL compute_face_point(face_data(1),(/1.0d0,cx(3)/),edge_point(:,2))
      CALL compute_face_point(face_data(1),(/cx(1),1.0d0/),edge_point(:,3))
      CALL compute_face_point(face_data(1),(/0.0d0,cx(3)/),edge_point(:,4))
      CALL compute_face_point(face_data(2),(/cx(1),0.d0/) ,edge_point(:,5))
      CALL compute_face_point(face_data(2),(/1.0d0,cx(3)/),edge_point(:,6))
      CALL compute_face_point(face_data(2),(/cx(1),1.0d0/),edge_point(:,7))
      CALL compute_face_point(face_data(2),(/0.0d0,cx(3)/),edge_point(:,8))
      CALL compute_face_point(face_data(4),(/cx(2),0.d0/) ,edge_point(:,10))
      CALL compute_face_point(face_data(6),(/cx(2),0.d0/) ,edge_point(:,9))
      CALL compute_face_point(face_data(4),(/cx(2),1.0d0/),edge_point(:,11))
      CALL compute_face_point(face_data(6),(/cx(2),1.0d0/),edge_point(:,12))
!
!     compute faces
!
      CALL compute_face_point(face_data(1),(/cx(1),cx(3)/),face_point(:,1))
      CALL compute_face_point(face_data(2),(/cx(1),cx(3)/),face_point(:,2))
      CALL compute_face_point(face_data(3),(/cx(1),cx(2)/),face_point(:,3))
      CALL compute_face_point(face_data(4),(/cx(2),cx(3)/),face_point(:,4))
      CALL compute_face_point(face_data(5),(/cx(1),cx(2)/),face_point(:,5))
      CALL compute_face_point(face_data(6),(/cx(2),cx(3)/),face_point(:,6))
!
!     compute edge derivatives
!
      CALL compute_face_derivative(face_data(1),(/cx(1),0.d0/) ,grad_2d)
      edge_der(:,1) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(1),(/1.0d0,cx(3)/),grad_2d)
      edge_der(:,2) = grad_2d(:,2)
      CALL compute_face_derivative(face_data(1),(/cx(1),1.0d0/),grad_2d)
      edge_der(:,3) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(1),(/0.0d0,cx(3)/),grad_2d)
      edge_der(:,4) = grad_2d(:,2)
      CALL compute_face_derivative(face_data(2),(/cx(1),0.d0/) ,grad_2d)
      edge_der(:,5) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(2),(/1.0d0,cx(3)/),grad_2d)
      edge_der(:,6) = grad_2d(:,2)
      CALL compute_face_derivative(face_data(2),(/cx(1),1.0d0/),grad_2d)
      edge_der(:,7) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(2),(/0.0d0,cx(3)/),grad_2d)
      edge_der(:,8) = grad_2d(:,2)
      CALL compute_face_derivative(face_data(6),(/cx(2),0.d0/) ,grad_2d)
      edge_der(:,9) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(4),(/cx(2),0.d0/) ,grad_2d)
      edge_der(:,10) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(4),(/cx(2),1.0d0/),grad_2d)
      edge_der(:,11) = grad_2d(:,1)
      CALL compute_face_derivative(face_data(6),(/cx(2),1.0d0/),grad_2d)
      edge_der(:,12) = grad_2d(:,1)
!
!     compute face derivatives
!
      CALL compute_face_derivative(face_data(1),(/cx(1),cx(3)/),face_der(:,:,1))
      CALL compute_face_derivative(face_data(2),(/cx(1),cx(3)/),face_der(:,:,2))
      CALL compute_face_derivative(face_data(3),(/cx(1),cx(2)/),face_der(:,:,3))
      CALL compute_face_derivative(face_data(4),(/cx(2),cx(3)/),face_der(:,:,4))
      CALL compute_face_derivative(face_data(5),(/cx(1),cx(2)/),face_der(:,:,5))
      CALL compute_face_derivative(face_data(6),(/cx(2),cx(3)/),face_der(:,:,6))
!
!     compute the mapping
!
      CALL compute_grad_hex_map(cx,grad_x,face_point,face_der, &
                                edge_point,edge_der,corner_point)
!
!     zero out rounded quantities
!
      eps = 100.d0*EPSILON(eps)
      DO i = 1,3
         DO j = 1,3
            IF( ABS(grad_x(i,j)) <= eps )     grad_x(i,j) = 0.0D0
         END DO
      END DO
!
      RETURN
      END SUBROUTINE general_hex_grad
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the gradient for a general hex element
!!
!!     cx          = computational space variable (X,Y,Z)
!!
!!     grad_x      = resultant gradient (x_X,x_Y,x_Z)
!!
!!
!!     face        = face interpolant location for appropriate local face coordinate
!!
!!     face_der    = face interpolant derivatives for appropriate local face coordinate
!!
!!     edge        = edge interpolant location for appropriate edge coordinate
!!
!!     edge_der    = edge interpolant derivative for appropriate edge coordinate
!!
!!     corner      = The 8 corners in standard fe format
!
      SUBROUTINE compute_grad_hex_map(cx,grad_x,face,face_der,edge,edge_der,corner)
!
!
      IMPLICIT none
!
!     input...
!
      DOUBLE PRECISION, DIMENSION(3)    ,INTENT(IN) :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION, DIMENSION(3,6)  ,INTENT(IN) :: face !< face interpolant location for appropriate local face coordinate
      DOUBLE PRECISION, DIMENSION(3,2,6),INTENT(IN) :: face_der !< face interpolant derivatives for appropriate local face coordinate
      DOUBLE PRECISION, DIMENSION(3,12) ,INTENT(IN) :: edge !< edge interpolant location for appropriate edge coordinate
      DOUBLE PRECISION, DIMENSION(3,12) ,INTENT(IN) :: edge_der !< edge interpolant derivative for appropriate edge coordinate
      DOUBLE PRECISION, DIMENSION(3,8)  ,INTENT(IN) :: corner !< 8 corners in standard fe format
!
!     output...
!
      DOUBLE PRECISION, DIMENSION(3,3),INTENT(OUT)  :: grad_x !<  resultant gradient (x_X,x_Y,x_Z)
!
!     local...
!
      INTEGER                                     :: j,iface
!-----
      DO j = 1,3
!
!        ------------
!        X derivative
!        ------------
!
         grad_x(j,1) =  -face(j,6) + face(j,4)                                  &
                       + face_der(j,1,1)*(1.d0 - cx(2)) + face_der(j,1,2)*cx(2) &
                       + face_der(j,1,3)*(1.d0 - cx(3)) + face_der(j,1,5)*cx(3)
!
         grad_x(j,1) = grad_x(j,1) - edge_der(j,1)*(1.d0 - cx(2))*(1.d0 - cx(3)) &
                                   - edge_der(j,3)*(1.d0 - cx(2))*        cx(3)  &
                                   - edge_der(j,5)*        cx(2) *(1.d0 - cx(3)) &
                                   - edge_der(j,7)*        cx(2) *        cx(3)  &
                                   + edge(j,9) *(1.d0 - cx(3)) &
                                   + edge(j,12)*        cx(3)  &
                                   - edge(j,10)*(1.d0 - cx(3)) &
                                   - edge(j,11)*        cx(3)  &
                                   + edge(j,4)*(1.d0 - cx(2))  &
                                   + edge(j,8)*        cx(2)   &
                                   - edge(j,2)*(1.d0 - cx(2))  &
                                   - edge(j,6)*        cx(2)
!
         grad_x(j,1) = grad_x(j,1) - corner(j,1)*(1.d0 -  cx(2))*(1.d0 -  cx(3)) &
                                   - corner(j,5)*(1.d0 -  cx(2))*         cx(3)  &
                                   - corner(j,4)*         cx(2) *(1.d0 -  cx(3)) &
                                   - corner(j,8)*         cx(2) *         cx(3)  &
                                   + corner(j,2)*(1.d0 -  cx(2))*(1.d0 -  cx(3)) &
                                   + corner(j,6)*(1.d0 -  cx(2))*         cx(3)  &
                                   + corner(j,3)*         cx(2) *(1.d0 -  cx(3)) &
                                   + corner(j,7)*         cx(2) *         cx(3)
!        ------------
!        Y derivative
!        ------------
!
         grad_x(j,2) =  face_der(j,1,6)*(1.d0 - cx(1)) + face_der(j,1,4)*cx(1) &
                      - face(j,1) + face(j,2) &
                      + face_der(j,2,3)*(1.d0 - cx(3)) + face_der(j,2,5)*cx(3)
!
         grad_x(j,2) = grad_x(j,2) + edge(j,1)*(1.d0 - cx(3)) &
                                   + edge(j,3)*        cx(3)  &
                                   - edge(j,5)*(1.d0 - cx(3)) &
                                   - edge(j,7)*        cx(3)  &
                                   - edge_der(j,9) *(1.d0 - cx(1))*(1.d0 - cx(3)) &
                                   - edge_der(j,12)*(1.d0 - cx(1))*        cx(3)  &
                                   - edge_der(j,10)*(1.d0 - cx(3))*        cx(1)  &
                                   - edge_der(j,11)*        cx(1) *        cx(3)  &
                                   + edge(j,4) *(1.d0 - cx(1)) &
                                   - edge(j,8) *(1.d0 - cx(1)) &
                                   + edge(j,2) *        cx(1)  &
                                   - edge(j,6) *        cx(1)
!
         grad_x(j,2) = grad_x(j,2) - corner(j,1)*(1.d0 -  cx(1))*(1.d0 -  cx(3)) &
                                   - corner(j,5)*(1.d0 -  cx(1))*         cx(3)  &
                                   + corner(j,4)*(1.d0 -  cx(1))*(1.d0 -  cx(3)) &
                                   + corner(j,8)*(1.d0 -  cx(1))*         cx(3)  &
                                   - corner(j,2)*         cx(1) *(1.d0 -  cx(3)) &
                                   - corner(j,6)*         cx(1) *         cx(3)  &
                                   + corner(j,3)*         cx(1) *(1.d0 -  cx(3)) &
                                   + corner(j,7)*         cx(1) *         cx(3)
!        ------------
!        Z derivative
!        ------------
!
         grad_x(j,3) =  face_der(j,2,6)*(1.d0 - cx(1)) + face_der(j,2,4)*cx(1) &
                      + face_der(j,2,1)*(1.d0 - cx(2)) + face_der(j,2,2)*cx(2) &
                      - face(j,3) + face(j,5)
!
         grad_x(j,3) = grad_x(j,3) + edge(j,1) *(1.d0 - cx(2))  &
                                   - edge(j,3) *(1.d0 - cx(2))  &
                                   + edge(j,5) *        cx(2)   &
                                   - edge(j,7) *        cx(2)   &
                                   + edge(j,9) *(1.d0 - cx(1))  &
                                   - edge(j,12)*(1.d0 - cx(1))  &
                                   + edge(j,10)*        cx(1)   &
                                   - edge(j,11)*        cx(1)   &
                                   - edge_der(j,4) *(1.d0 - cx(1))*(1.d0 - cx(2)) &
                                   - edge_der(j,8) *(1.d0 - cx(1))*        cx(2)  &
                                   - edge_der(j,2) *        cx(1) *(1.d0 - cx(2)) &
                                   - edge_der(j,6) *        cx(1) *        cx(2)
!
         grad_x(j,3) = grad_x(j,3) - corner(j,1)*(1.d0 -  cx(1))*(1.d0 -  cx(2))  &
                                   + corner(j,5)*(1.d0 -  cx(1))*(1.d0 -  cx(2))  &
                                   - corner(j,4)*(1.d0 -  cx(1))*         cx(2)   &
                                   + corner(j,8)*(1.d0 -  cx(1))*         cx(2)   &
                                   - corner(j,2)*         cx(1) *(1.d0 -  cx(2))  &
                                   + corner(j,6)*         cx(1) *(1.d0 -  cx(2))  &
                                   - corner(j,3)*         cx(1) *         cx(2)   &
                                   + corner(j,7)*         cx(1) *         cx(2)
      END DO
!
      RETURN
      END SUBROUTINE compute_grad_hex_map
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the derivative in the two local coordinate directions
!> on a subdomain face according by interpolating the face_data
!!
!
      SUBROUTINE compute_face_derivative(face_data,cx,grad)
!
!
      USE FE_Data_Types

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      TYPE(face_interp)                  :: face_data !< interpolant data for the 6 faces
      DOUBLE PRECISION , DIMENSION(2)    :: cx !< computational space variable (X,Y,Z)
      DOUBLE PRECISION , DIMENSION(3,2)  :: grad !< gradient 
!
!!     use bi-linear mapping of this is a flat side (for speed) otherwise
!!     do general lagrange interpolation
!
      IF( face_data%num_knots(1) == 2 .AND. face_data%num_knots(2) == 2)     THEN
         DO j = 1,3
            grad(j,1)  =  -face_data%points(j,1,1)*(1.d0 - cx(2)) &
                         + face_data%points(j,2,1)*(1.d0 - cx(2)) &
                         + face_data%points(j,2,2)       * cx(2)  &
                         - face_data%points(j,1,2)       * cx(2)

            grad(j,2)  =  -face_data%points(j,1,1)*(1.d0 - cx(1)) &
                         - face_data%points(j,2,1)*        cx(1)  &
                         + face_data%points(j,2,2)*        cx(1)  &
                         + face_data%points(j,1,2)*(1.d0 - cx(1))
         END DO
      ELSE
         CALL compute_Poly_2D_der(face_data,cx,grad)
      END IF

      RETURN
      END SUBROUTINE compute_face_derivative
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the lagrange interpolant through the points p(i,j)
!>  at the location (x,y)
!!
!
      SUBROUTINE compute_Poly_2D_der(face_data,cx,grad)
!
!
      USE FE_Data_Types

      IMPLICIT none
!
      TYPE(face_interp)                 :: face_data !< interpolant data 
      DOUBLE PRECISION, DIMENSION(2)    :: cx !< computational space variable
      DOUBLE PRECISION, DIMENSION(3,2)  :: grad
      DOUBLE PRECISION, DIMENSION(nmax) :: l_i
      DOUBLE PRECISION                  :: l_j
      DOUBLE PRECISION                  :: polynder,polyn
      INTEGER                           :: i,j,k
      EXTERNAL polynder,polyn
!
      grad = 0.0d0
!
!     first direction derivative
!
      DO i = 1,face_data%num_knots(1)
         l_i(i) = polynder(i,cx(1),face_data%num_knots(1),face_data%knots(:,1))
      END DO
!
      DO j = 1,face_data%num_knots(2)
         l_j = polyn(j,cx(2),face_data%num_knots(2),face_data%knots(:,2))
         DO i = 1,face_data%num_knots(1)
            DO k = 1,3
               grad(k,1) = grad(k,1) + face_data%points(k,i,j)*l_i(i)*l_j
            END DO
         END DO
      END DO
!
!     second direction derivative
!
      DO i = 1,face_data%num_knots(1)
         l_i(i) = polyn(i,cx(1),face_data%num_knots(1),face_data%knots(:,1))
      END DO
!
      DO j = 1,face_data%num_knots(2)
         l_j = polynder(j,cx(2),face_data%num_knots(2),face_data%knots(:,2))
         DO i = 1,face_data%num_knots(1)
            DO k = 1,3
               grad(k,2) = grad(k,2) + face_data%points(k,i,j)*l_i(i)*l_j
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE compute_Poly_2D_der
!
!     ///////////////////////////////////////////////////////////////////////
!> @brief Compute the inverse of the metric tensor and the jacobian.
!>  The gradient is stored as d x_i/d X_j
!!
!!     grad_x  = physical space gradient
!!
!!               grad_x(:,1) = (x_X,y_X,z_X), etc.
!!
!!     grad_cx = computational space gradient
!!
!!               grad_cx(:,1) = (X_x,Y_x,Z_x), etc.
!!
!!     J       = transformation jacobian
!!
      SUBROUTINE Invert_Metric_Tensor(grad_x,grad_cx,J)
!
!     ......................................................................
!
      IMPLICIT none

      DOUBLE PRECISION, DIMENSION(3,3),INTENT(IN)  :: grad_x !< physical space gradient
      DOUBLE PRECISION, DIMENSION(3,3),INTENT(OUT) :: grad_cx !< computational space gradient
      DOUBLE PRECISION                             :: J !< jacobian
	  DOUBLE PRECISION                             :: small_value
      INTEGER                                      :: n,m
!
      small_value = 100*EPSILON(small_value)
!
      J =    grad_x(1,1)*grad_x(2,2)*grad_x(3,3) &
           - grad_x(1,2)*grad_x(2,1)*grad_x(3,3) &
           - grad_x(1,1)*grad_x(2,3)*grad_x(3,2) &
           + grad_x(1,3)*grad_x(2,1)*grad_x(3,2) &
           + grad_x(1,2)*grad_x(2,3)*grad_x(3,1) &
           - grad_x(1,3)*grad_x(2,2)*grad_x(3,1)

      grad_cx(1,1) =  ( grad_x(2,2)*grad_x(3,3)    & 
                     -  grad_x(2,3)*grad_x(3,2))/J
      grad_cx(2,1) =  (-grad_x(2,1)*grad_x(3,3)    & 
                     +  grad_x(2,3)*grad_x(3,1))/J
      grad_cx(3,1) =  ( grad_x(2,1)*grad_x(3,2)    & 
                     -  grad_x(3,1)*grad_x(2,2))/J

      grad_cx(1,2) =  (-grad_x(1,2)*grad_x(3,3)    & 
                     +  grad_x(1,3)*grad_x(3,2))/J
      grad_cx(2,2) =  ( grad_x(1,1)*grad_x(3,3)    & 
                     -  grad_x(1,3)*grad_x(3,1))/J
      grad_cx(3,2) =  (-grad_x(1,1)*grad_x(3,2)    & 
                     +  grad_x(1,2)*grad_x(3,1))/J

      grad_cx(1,3) =  ( grad_x(1,2)*grad_x(2,3)    & 
                     -  grad_x(2,2)*grad_x(1,3))/J
      grad_cx(2,3) =  (-grad_x(1,1)*grad_x(2,3)    & 
                     +  grad_x(1,3)*grad_x(2,1))/J
      grad_cx(3,3) =  ( grad_x(1,1)*grad_x(2,2)    & 
                     -  grad_x(1,2)*grad_x(2,1))/J
!
      DO m = 1,3
         DO n = 1,3
            IF ( ABS(grad_cx(n,m)) <= small_value )     grad_cx(n,m) = 0.0d0
         END DO
      END DO
!
      RETURN
      END SUBROUTINE Invert_Metric_Tensor
