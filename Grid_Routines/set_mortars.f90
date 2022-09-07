!> @file
!! This file contains the subroutines that generate and set up the mortars
!!
!///////////////////////////////////////////////////////////////////////////////
!////////								 ////////
!////////       set_mortars.f90						 ////////
!////////								 ////////
!////////       contains:						 ////////
!////////								 ////////
!////////             SUBROUTINE generate_mortars(d,ngrid,mrtr,nmort)    ////////
!////////             SUBROUTINE set_mortar(mrtr,mrtr_no,d,the_face)     ////////
!////////								 ////////
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine generates the mortars and element connections.
!!
!! First, set up the element connections through the mortars. Second, 
!! generate faces for the mortars, count total number. Finally, 
!! generate the mortars and element connections.
      SUBROUTINE generate_mortars(d,ngrid,mrtr,nmort)
!
!.......................................................................
!     set up the element connections through the mortars
!.......................................................................
!
      USE face_operations
      USE Mortar_Definition
      USE File_Units
!
      TYPE (domain),DIMENSION(ngp)        :: d !< element list
      TYPE (mortar),DIMENSION(nmp)        :: mrtr !< mortar list
      TYPE (face_list), POINTER           :: face_list_head,face_list_tail,the_face
!
!     --------------------------------------------------
!     generate faces for the mortars, count total number
!     --------------------------------------------------
!
      CALL make_list_of_faces(d,ngrid,face_list_head,face_list_tail)
      CALL count_faces(face_list_head,nmort)
!
!     ------------------------------------------------
!     now generate the mortars and element connections
!     ------------------------------------------------
!
      IF(.not.ASSOCIATED(face_list_head))     THEN
         WRITE(ioutunit,*) 'the face list is empty'
         STOP 'ERROR: no faces in the face list'
      ELSE
         the_face => face_list_head
         mrtr_number = 1
         DO
            CALL set_mortar(mrtr(mrtr_number),mrtr_number,d,the_face)
            the_face => the_face%next
            if(.not.ASSOCIATED(the_face))     EXIT
            mrtr_number = mrtr_number + 1
         END DO
      END IF
!
!     -----------------------------------
!     the face lists are no longer needed
!     -----------------------------------
!
      CALL deallocate_face_list(face_list_head)
!
      END SUBROUTINE generate_mortars
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets up the mortar information, such as mapping to the right
!> interface of the element, setting up data for mortar.
!!
!!     Mapping arrays for transformations follow. Given the face number,
!!     edgemap returns the coordinate directions of the face. For example
!!     face 1 defines a square in the (X,Z) plane, so
!!     edgemap(1,1) = 1 and edgemap(2,1) = 3.
!!
!!     Next, sign_map returns the
!!     sign of the dot product of the normals of two faces that form
!!     a mortar. For example, if two faces of the same face number form
!!     an interface, then their normals must point in opposite directions.
!!     In the case that faces 1 of two elements form an interface, then
!!     sign_map(1,1) = -1.
!
       SUBROUTINE set_mortar(mrtr,mrtr_no,d,the_face)
!.......................................................................
!     set the mortar information
!.......................................................................
!
      USE FE_Data_Types
      USE Mortar_Definition
      USE domain_definition
      USE Edge_Mappings
!
      TYPE (domain),DIMENSION(ngp)  :: d !< element list
      TYPE (mortar)                 :: mrtr !< mortar list
      TYPE (face_list)              :: the_face !< face list
!
!
!--------------------------------------------------------------------------
!     Mapping arrays for transformations follow. Given the face number,
!     edgemap returns the coordinate directions of the face. For example
!     face 1 defines a square in the (X,Z) plane, so
!     edgemap(1,1) = 1 and edgemap(2,1) = 3.
!
!     Next, sign_map returns the
!     sign of the dot product of the normals of two faces that form
!     a mortar. For example, if two faces of the same face number form
!     an interface, then their normals must point in opposite directions.
!     In the case that faces 1 of two elements form an interface, then
!     sign_map(1,1) = -1.
!
!--------------------------------------------------------------------------
!

      INTEGER, DIMENSION(2,6) :: edge_map = RESHAPE((/1,3,1,3,1,2,2,3,1,2,2,3/),(/2,6/))
      INTEGER, DIMENSION(6,6) :: sign_map = RESHAPE((/-1, 1,-1, 1, 1,-1,  &
                                                       1,-1, 1,-1,-1, 1, &
                                                      -1, 1,-1, 1, 1,-1, &
                                                        1,-1, 1,-1,-1, 1, &
                                                       1,-1, 1,-1,-1, 1, &
                                                       -1, 1,-1, 1, 1,-1/),(/6,6/))
!
!    ---------------
!     set mortar data
!     ---------------
!
       idl = the_face%face_data%from_element(1)
       idr = the_face%face_data%from_element(2)
       ifl = the_face%face_data%from_side(1)
       ifr = the_face%face_data%from_side(2)
!
       mrtr%id    = the_face%face_data%from_element
       mrtr%iface = the_face%face_data%from_side
!
       mrtr%len(1,1) = d(idl)%ncg(edge_map(1,ifl))
       mrtr%len(2,1) = d(idl)%ncg(edge_map(2,ifl))
!
!     -------------------------------------------------
!     determine whether the mortar is conforming or not
!     -------------------------------------------------
!
      IF (idr == 0 )     THEN ! this is a boundary mortar (always conforming)

         mrtr%lenmortar(1) = mrtr%len(1,1)
         mrtr%lenmortar(2) = mrtr%len(2,1)
         mrtr%len(1,2)     = 0
         mrtr%len(2,2)     = 0
         mrtr%conforming   = .true.
         mrtr%orient       = 0
         mrtr%sign         = 1

      ELSE                   ! this is an interface mortar

         mrtr%sign          = sign_map(ifl,ifr)
         mrtr%orient        = d(idr)%orientation(ifr)
         mrtr%len(1,2)      = d(idr)%ncg(edge_map(1,ifr))
         mrtr%len(2,2)      = d(idr)%ncg(edge_map(2,ifr))
         SELECT CASE (mrtr%orient)
            CASE (DFLT,B1F2,B1B2,F1B2)
               j1 = 1
               j2 = 2
            CASE (F2F1,B2F1,B2B1,F2B1)
               j1 = 2
               j2 = 1
         END SELECT
!
!        --------------------------------------------------------------
!        choose mortar dimension so that outflow condition is satisfied
!        --------------------------------------------------------------
!
         mrtr%lenmortar(1) = max(mrtr%len(1,1),mrtr%len(j1,2))
         mrtr%lenmortar(2) = max(mrtr%len(2,1),mrtr%len(j2,2))
!
!        ------
!        master
!        ------
!
         IF ( mrtr%lenmortar(1) == mrtr%len(1,1) ) THEN
            mrtr%conforming(1,1) = .true.
         ELSE
            mrtr%conforming(1,1) = .false.
         END IF
         IF ( mrtr%lenmortar(2) == mrtr%len(2,1) ) THEN
            mrtr%conforming(2,1) = .true.
         ELSE
            mrtr%conforming(2,1) = .false.
         END IF
!
!        ------
!        slave
!        ------
!
         IF ( mrtr%lenmortar(1) == mrtr%len(j1,2) ) THEN
              mrtr%conforming(1,2) = .true.
         ELSE
              mrtr%conforming(1,2) = .false.
         END IF
         IF ( mrtr%lenmortar(2) == mrtr%len(j2,2) ) THEN
              mrtr%conforming(2,2) = .true.
         ELSE
              mrtr%conforming(2,2) = .false.
         END IF

      END IF
!
!     ---------------
!     set element data
!     ---------------
!
      d(idl)%mortar(1,ifl) = mrtr_no
      d(idl)%mortar(2,ifl) = 1
      IF( idr == 0 )     THEN
         d(idl)%ibtype(ifl) = 0
      ELSE
         d(idl)%ibtype(ifl)   = 1
         d(idr)%ibtype(ifr)   = 1
         d(idl)%bcond(ifl)    = 'interface'
         d(idr)%mortar(1,ifr) = mrtr_no
         d(idr)%mortar(2,ifr) = 2
         d(idr)%bcond(ifr)    = 'interface'
      END IF
!
      RETURN
      END SUBROUTINE set_mortar


!
!///////////////////////////////////////////////////////////////////////
!> @brief
!! This subroutine computes whether the normal of the master and slave side
!! have the same sign.
!!     sj = dot(n_hat_master,n_hat_slave)
!!     if sj positive then they have the same sign, if negative
!!     then they have opposite sign.
! 
!     date: 11/21/01
!     programmer: Guus Jacobs
!///////////////////////////////////////////////////////////////////////

      SUBROUTINE set_mortar_nsign(mrtr,d)

      USE mortar_definition
      USE domain_definition
      USE Edge_Mappings
 
      IMPLICIT DOUBLE PRECISION(a-h,o-z)

      TYPE(domain) :: d(ngp)
      TYPE(mortar) :: mrtr
  
      DOUBLE PRECISION :: n_hat(3,2)

      idm = mrtr%id(1)
      ids = mrtr%id(2)
      ifm = mrtr%iface(1)
      ifs = mrtr%iface(2)


    IF (d(idm)%ibtype(ifm) == 1) THEN
!    Compute the normal along the master and slave domain faces
 
      
     SELECT CASE(ifm) 
       CASE(1)
          sqr         = sqrt(d(idm)%gmet(2,1,1,1,1)**2+ d(idm)%gmet(2,2,1,1,1)**2 +   &
                        d(idm)%gmet(2,3,1,1,1)**2)

          n_hat(1,1)  = d(idm)%gmet(2,1,1,1,1)/sqr 
          n_hat(2,1)  = d(idm)%gmet(2,2,1,1,1)/sqr 
          n_hat(3,1)  = d(idm)%gmet(2,3,1,1,1)/sqr 
       CASE(2)
          nc          = d(idm)%nc(2)
          sqr         = sqrt(d(idm)%gmet(2,1,1,nc,1)**2+ d(idm)%gmet(2,2,1,nc,1)**2 +   &
                        d(idm)%gmet(2,3,1,nc,1)**2)
          n_hat(1,1)  = d(idm)%gmet(2,1,1,nc,1)/sqr 
          n_hat(2,1)  = d(idm)%gmet(2,2,1,nc,1)/sqr 
          n_hat(3,1)  = d(idm)%gmet(2,3,1,nc,1)/sqr 
       CASE(3)
          sqr         = sqrt(d(idm)%gmet(3,1,1,1,1)**2+ d(idm)%gmet(3,2,1,1,1)**2 +   &
                        d(idm)%gmet(3,3,1,1,1)**2)
          n_hat(1,1)  = d(idm)%gmet(3,1,1,1,1)/sqr 
          n_hat(2,1)  = d(idm)%gmet(3,2,1,1,1)/sqr 
          n_hat(3,1)  = d(idm)%gmet(3,3,1,1,1)/sqr 
       CASE(4)
          nc          = d(idm)%nc(1)
          sqr         = sqrt(d(idm)%gmet(1,1,nc,1,1)**2+ d(idm)%gmet(1,2,nc,1,1)**2 +   &
                        d(idm)%gmet(1,3,nc,1,1)**2)
          n_hat(1,1)  = d(idm)%gmet(1,1,nc,1,1)/sqr 
          n_hat(2,1)  = d(idm)%gmet(1,2,nc,1,1)/sqr 
          n_hat(3,1)  = d(idm)%gmet(1,3,nc,1,1)/sqr 
       CASE(5)
          nc          = d(idm)%nc(3)
          sqr         = sqrt(d(idm)%gmet(3,1,1,1,nc)**2+ d(idm)%gmet(3,2,1,1,nc)**2 +   &
                        d(idm)%gmet(3,3,1,1,nc)**2)
          n_hat(1,1)  = d(idm)%gmet(3,1,1,1,nc)/sqr 
          n_hat(2,1)  = d(idm)%gmet(3,2,1,1,nc)/sqr 
          n_hat(3,1)  = d(idm)%gmet(3,3,1,nc,1)/sqr 
       CASE(6)
          sqr         = sqrt(d(idm)%gmet(1,1,1,1,1)**2+ d(idm)%gmet(1,2,1,1,1)**2 +   &
                        d(idm)%gmet(1,3,1,1,1)**2)
          n_hat(1,1)  = d(idm)%gmet(1,1,1,1,1)/sqr 
          n_hat(2,1)  = d(idm)%gmet(1,2,1,1,1)/sqr 
          n_hat(3,1)  = d(idm)%gmet(1,3,1,1,1)/sqr 
     END SELECT


     SELECT CASE(ifs) 
       CASE(1)
          sqr         = sqrt(d(ids)%gmet(2,1,1,1,1)**2+ d(ids)%gmet(2,2,1,1,1)**2 +   &
                        d(ids)%gmet(2,3,1,1,1)**2)
          n_hat(1,2)  = d(ids)%gmet(2,1,1,1,1)/sqr 
          n_hat(2,2)  = d(ids)%gmet(2,2,1,1,1)/sqr 
          n_hat(3,2)  = d(ids)%gmet(2,3,1,1,1)/sqr 
       CASE(2)
          nc          = d(ids)%nc(2)
          sqr         = sqrt(d(ids)%gmet(2,1,1,nc,1)**2+ d(ids)%gmet(2,2,1,nc,1)**2 +   &
                        d(ids)%gmet(2,3,1,nc,1)**2)
          n_hat(1,2)  = d(ids)%gmet(2,1,1,nc,1)/sqr 
          n_hat(2,2)  = d(ids)%gmet(2,2,1,nc,1)/sqr 
          n_hat(3,2)  = d(ids)%gmet(2,3,1,nc,1)/sqr 
       CASE(3)
          sqr         = sqrt(d(ids)%gmet(3,1,1,1,1)**2+ d(ids)%gmet(3,2,1,1,1)**2 +   &
                        d(ids)%gmet(3,3,1,1,1)**2)
          n_hat(1,2)  = d(ids)%gmet(3,1,1,1,1)/sqr 
          n_hat(2,2)  = d(ids)%gmet(3,2,1,1,1)/sqr 
          n_hat(3,2)  = d(ids)%gmet(3,3,1,1,1)/sqr 
       CASE(4)
          nc          = d(ids)%nc(1)
          sqr         = sqrt(d(ids)%gmet(1,1,nc,1,1)**2+ d(ids)%gmet(1,2,nc,1,1)**2 +   &
                        d(ids)%gmet(1,3,nc,1,1)**2)
          n_hat(1,2)  = d(ids)%gmet(1,1,nc,1,1)/sqr 
          n_hat(2,2)  = d(ids)%gmet(1,2,nc,1,1)/sqr 
          n_hat(3,2)  = d(ids)%gmet(1,3,nc,1,1)/sqr 
       CASE(5)
          nc          = d(ids)%nc(3)
          sqr         = sqrt(d(ids)%gmet(3,1,1,1,nc)**2+ d(ids)%gmet(3,2,1,1,nc)**2 +   &
                        d(ids)%gmet(3,3,1,1,nc)**2)
          n_hat(1,2)  = d(ids)%gmet(3,1,1,1,nc)/sqr 
          n_hat(2,2)  = d(ids)%gmet(3,2,1,1,nc)/sqr 
          n_hat(3,2)  = d(ids)%gmet(3,3,1,nc,1)/sqr 
       CASE(6)
          sqr         = sqrt(d(ids)%gmet(1,1,1,1,1)**2+ d(ids)%gmet(1,2,1,1,1)**2 +   &
                        d(ids)%gmet(1,3,1,1,1)**2)
          n_hat(1,2)  = d(ids)%gmet(1,1,1,1,1)/sqr 
          n_hat(2,2)  = d(ids)%gmet(1,2,1,1,1)/sqr 
          n_hat(3,2)  = d(ids)%gmet(1,3,1,1,1)/sqr 
     END SELECT

!
!    Compute sj and nsign
!

     sj = 0
     DO i=1,3
       sj = sj + n_hat(i,1)*n_hat(i,2) 
     ENDDO
!
     sgn = SIGN(1.0d0,sj)
     mrtr%nsign = NINT(sgn)
   ELSE 
     mrtr%nsign = 1
   ENDIF

!  SELECT CASE(mrtr%orient)
!    CASE(B1F2,F1B2,B2F1,F2B1)
!      mrtr%nsign = -mrtr%nsign
!  END SELECT
!     write(*,*) 'mrtr%nsign',mrtr%nsign,mrtr%sign

 
      RETURN
     END SUBROUTINE set_mortar_nsign













