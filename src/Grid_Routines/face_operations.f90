!> \file face_operations.f90
!! Contains all finite element operations related to the element faces
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!////////											////////
!////////	face_operations.f90								////////
!////////											////////
!////////	contains:									////////
!////////											////////
!////////	      MODULE face_operations							////////
!////////	      SUBROUTINE make_list_of_faces(d,num_elements,face_list_head,face_list_tail)////////
!////////	      SUBROUTINE initialize_face_list(head,tail)				////////
!////////	      SUBROUTINE initialize_first_face(new_face,head,tail)			////////
!////////	      SUBROUTINE make_face(id,face_num,new_face)				////////
!////////	      SUBROUTINE set_face_data(id,face_num,new_face,the_element)		////////
!////////	      SUBROUTINE Insert_to_face_List(element_id,element_face,new_face,head,tail)////////
!////////	      SUBROUTINE Delete_face(face_to_delete)					////////
!////////	      SUBROUTINE Insert_face(new_face,previous,next)				////////
!////////	      SUBROUTINE append_face(new_face,tail)					////////
!////////	      SUBROUTINE prepend_face(new_face,head)					////////
!////////	      LOGICAL FUNCTION faces_match(target_nodes,test_nodes,orientation)		////////
!////////	      SUBROUTINE count_faces(head,number_of_faces)				////////
!////////	      SUBROUTINE gather_boundary_faces(boundary_list_head,boundary_list_tail,&	////////
!////////	      SUBROUTINE print_face_list(ioutunit,head)					////////
!////////	      SUBROUTINE print_face_list_backwards(ioutunit,tail)			////////
!////////	      SUBROUTINE deallocate_face_list(head)					////////
!////////											////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!> This module contains the routines necessary to create and traverse
!! the doubly linked list for the faces linked list
      MODULE face_operations
!
      USE Domain_Definition
      USE FE_Data_Types
      USE Element_Library
      USE Edge_Mappings
!
      PUBLIC  :: make_list_of_faces,print_face_list,Delete_face,count_faces
!
      PRIVATE :: make_face,initialize_face_list,initialize_first_face,Insert_to_face_List
      PRIVATE :: Insert_face,append_face,prepend_face
!
!     contains the routines necessary to create and traverse the doubly linked list
!     for the faces linked list
!
!
!     --------
      CONTAINS
!     --------
!
!> @brief This routine takes the element list and makes a list of the faces.
!!
!! Setting links from the element to the face and back from the face to the element. These are done by
!! initializing the list, and then setting up the first six faces from the first element. 
!! After that, check the face list to see if the face to be added is already in the list.
      SUBROUTINE make_list_of_faces(d,num_elements,face_list_head,face_list_tail)
!
!     this routine takes the element list and makes a list of the faces.
!     Links to the face and back from the face to the element are set.
!
      TYPE (domain),DIMENSION(ngp) :: d
      TYPE (face_list), POINTER    :: face_list_head,face_list_tail,new_face
      INTEGER                      :: type, orientation
!
!-----
!
!     initialize the list, and then set up the first six faces from the
!     first element. After that, check the face list to see if the face
!     to be added is already in the list.
!
      CALL initialize_face_list(face_list_head,face_list_tail)
!
!     first element
!
      id    = 1
      iface = 1
      CALL make_face(id,iface,new_face)
      CALL initialize_first_face(new_face,face_list_head,face_list_tail)
      CALL set_face_data(id,iface,new_face,d(id))

      DO iface = 2,6
         CALL make_face(id,iface,new_face)
         CALL set_face_data(id,iface,new_face,d(id))
         CALL Insert_to_face_List(id,iface,new_face,face_list_head, &
                                  face_list_tail,orientation)
         d(id)%orientation(iface) = orientation
      END DO
!
!     subsequent elements
!
      DO id = 2, num_elements   ! for each element
         DO iface = 1,6                  ! for each side of a hex element
            CALL make_face(id,iface,new_face)
            CALL set_face_data(id,iface,new_face,d(id))
            CALL Insert_to_face_List(id,iface,new_face,face_list_head, &
                                     face_list_tail,orientation)
            d(id)%orientation(iface) = orientation
          END DO
      END DO

      RETURN   
!
      END SUBROUTINE make_list_of_faces
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Initialize the face list
!!
!! Nullify the pointers named head and tail of the face lists
!
      SUBROUTINE initialize_face_list(head,tail)
!
!     nullify head and tail
!
      TYPE (face_list), POINTER :: head,tail
!
!-----
!
      NULLIFY (head,tail)

      RETURN
      END SUBROUTINE initialize_face_list
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Add the first element to the list
      SUBROUTINE initialize_first_face(new_face,head,tail)
!
!     add the first element to the list
!
      TYPE (face_list), POINTER :: head,tail
      TYPE (face_list), POINTER :: new_face
!
!---- 
!
      NULLIFY(new_face%previous,new_face%next)
      head => new_face
      tail => new_face

      RETURN
      END SUBROUTINE initialize_first_face
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Make element face
!!
!! Allocate space for face list and zero out the list
      SUBROUTINE make_face(id,face_num,new_face)
!
!     make a face record and zero out fields
!
      TYPE (face_list), POINTER :: new_face
      INTEGER                   :: face_num,id
!
!-----
!
      ALLOCATE(new_face, STAT = ierr)
      IF( ierr /= 0)     THEN
         WRITE(6,*) 'cannot allocate new face for element ',id,' face ',face_num
         STOP 'face allocation error'
      END IF

      ALLOCATE(new_face%face_data, STAT = ier)
      IF( ierr /= 0)     THEN
         WRITE(6,*) 'cannot allocate new face data for element ', &
                     id,' face ',face_num
         STOP 'face data allocation error'
      END IF

      new_face%face_data%key            = 0
      new_face%face_data%node           = 0
      new_face%face_data%from_element   = 0
      new_face%face_data%from_side      = 0

      NULLIFY (new_face%previous,new_face%next)

      RETURN
      END SUBROUTINE make_face
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine sets the face nodes used in the face and the connectivity between face and element
      SUBROUTINE set_face_data(id,face_num,new_face,the_element)
!
!     set the face nodes used in this face and the face/element
!     connectivity
!
      TYPE (face_list), POINTER :: new_face
      TYPE (domain)             :: the_element
      INTEGER                   :: face_num,id, min_node
!
!     face to element connectivity
!
      new_face%face_data%from_element(1) = id
      new_face%face_data%from_side(1)    = face_num
!
!     find the nodes associated with this face in standard order
!
      SELECT CASE(the_element%type)

         CASE(8) ! hex8 element - 1st order
            DO k = 1,4
               new_face%face_data%node(k) = the_element%node(Face_Map_hex8(k,face_num))
            END DO

         CASE(26) ! hex26 - second order
            DO k = 1,9
               new_face%face_data%node(k) = the_element%node(Face_Map_hex26(k,face_num))
            END DO

         CASE(64) ! hex64 - third order
            DO k = 1,16
               new_face%face_data%node(k) = the_element%node(Face_Map_hex64(k,face_num))
            END DO

      END SELECT

      min_node = MINVAL(new_face%face_data%node(1:4))
      new_face%face_data%key = min_node
!
      RETURN
      END SUBROUTINE set_face_data
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine works as inserting a face into the list
!!
!! Insert a face into the list. If the face already exists in the
!! list, delete this one and point to the one that's already in there.
      SUBROUTINE Insert_to_face_List(element_id,element_face,new_face, &
                                     head,tail,orientation)
!
!     insert a face into the list. If the face already exists in the
!     list, delete this one and point to the one that's already in there
!
      TYPE (face_list), POINTER :: head,tail,current,next,temp
      TYPE (face_list), POINTER :: new_face
      INTEGER                   :: element_id,element_face, current_key, orientation
!-----
!
      key = new_face%face_data%key
      orientation = 0

      IF(key < head%face_data%key )          THEN ! insert face before head of list
         CALL prepend_face(new_face,head)
         new_face%face_data%from_element(1) = element_id
         new_face%face_data%from_side(1)    = element_face
         RETURN

      ELSE IF (key > tail%face_data%key)     THEN ! insert face at end of list
         CALL append_face(new_face,tail)
         new_face%face_data%from_element(1) = element_id
         new_face%face_data%from_side(1)    = element_face
         RETURN

      ELSE                                        ! find location in list and insert
!
!        search through list from beginning until
!        we go past all key values less than the
!        target key
!
         current => head
         DO
            IF ( key <= current%face_data%key ) EXIT
            current => current%next
            IF( .NOT.ASSOCIATED(current) )     THEN
               WRITE(6,*) 'end of list hit ', key,element_face
               STOP 'end of face list hit'
            END IF
         END DO
!
!        if this is a new key, insert the face and we are done
!     
         IF ( key < current%face_data%key )     THEN ! new unique key
            current => current%previous
            next => current%next
            CALL Insert_face(new_face,current,next)          ! this face with this key 
            new_face%face_data%from_element(1) = element_id  ! isn't in list so add it
            new_face%face_data%from_side(1)    = element_face
            RETURN
         END IF
!
!        this key is already in the list. search until the last
!        face with this key. if a face matches along the way, point to that face.
!        otherwise add to the list
!
         DO
!
            IF ( faces_match(current%face_data%node(1:4), &
                             new_face%face_data%node(1:4),orientation) )     THEN
               temp => new_face
               new_face  => current
               DEALLOCATE(temp)
               new_face%face_data%from_element(2) = element_id
               new_face%face_data%from_side(2)    = element_face
               RETURN
            END IF

            IF ( .NOT.ASSOCIATED(current%next) )     THEN ! we're at end of list
               CALL append_face(new_face,tail)
               new_face%face_data%from_element(1) = element_id
               new_face%face_data%from_side(1)    = element_face
               RETURN
            END IF

            IF ( current%next%face_data%key > key )   THEN ! last one with this key
               next => current%next
               CALL Insert_face(new_face,current,next)          ! this face with this key 
               new_face%face_data%from_element(1) = element_id  ! isn't in list so add it
               new_face%face_data%from_side(1)    = element_face
               RETURN
            END IF

            current => current%next

         END DO

      END IF

      STOP 'funny termination of Insert_to_face_List'
      END SUBROUTINE Insert_to_face_List
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine functions as deleting a face in the list
      SUBROUTINE Delete_face(face_to_delete)
!
!     delete a face in the list.
!
      TYPE (face_list), POINTER :: previous,next
      TYPE (face_list), POINTER :: face_to_delete
!-----
!
      previous => face_to_delete%previous
      next     => face_to_delete%next
      previous%next => next
      next%previous => previous
      DEALLOCATE(face_to_delete%face_data)
      DEALLOCATE(face_to_delete)

      RETURN
      END SUBROUTINE Delete_face
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine functions as inserting a face between previous and next one
      SUBROUTINE Insert_face(new_face,previous,next)
!
!     insert a face between previous and next
!
      TYPE (face_list), POINTER :: previous,next
      TYPE (face_list), POINTER :: new_face
      INTEGER              :: key
!-----
!
      previous%next => new_face
      new_face%next => next

      next%previous     => new_face
      new_face%previous => previous

      RETURN
      END SUBROUTINE Insert_face
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine functions as adding a new element to the end of list
      SUBROUTINE append_face(new_face,tail)
!
!     add a new element to the end of list
!
      TYPE (face_list), POINTER :: tail
      TYPE (face_list), POINTER :: new_face
!-----
!
      tail%next         => new_face
      new_face%previous => tail
      NULLIFY(new_face%next)
      tail => new_face

      RETURN
   END SUBROUTINE append_face
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine works as adding a new element to the head of the list
      SUBROUTINE prepend_face(new_face,head)
!
!     add a new element to the head of the list
!
      TYPE (face_list), POINTER :: head
      TYPE (face_list), POINTER :: new_face
!-----
!
      new_face%next => head
      head%previous => new_face
      NULLIFY(new_face%previous)
      head => new_face

      RETURN
      END SUBROUTINE prepend_face
!
!///////////////////////////////////////////////////////////////////////
!> @brief Test to see if two faces match
!!
!! Two faces match if the four nodes are the same. Note that if one of 
!! them fails to match, then the face fails to match.Orientation gives 
!! the relative orientation of the test face to the target face. In this 
!! routine, orientation is measured in 90 degree increments: angle = orientation*pi/2
      LOGICAL FUNCTION faces_match(target_nodes,test_nodes,orientation)
!
!......................................................................
!     Test to see if two faces match. Two faces match iff the four 
!     nodes are the same. Note that if one of them fails to match,
!     then the face fails to match.Orientation gives the relative  
!     orientation of the test face to the target face. In this routine,
!     orientation is measured in 90 degree increments: angle = orientation*pi/2
!......................................................................
!
      INTEGER, DIMENSION(4), INTENT(IN) :: target_nodes,test_nodes
      INTEGER, INTENT(OUT)              :: orientation
      LOGICAL                           :: hit
      INTEGER, DIMENSION(4)             :: c
!
      faces_match       = .false.

      n = 0
      DO i = 1,4
         hit = .false.
         DO k = 1,4
            IF( target_nodes(i) == test_nodes(k) )     then
               hit = .true.
               c(i) = k
               n = n + 1
               EXIT
            END IF
         END DO
         IF (.NOT.hit)     RETURN
      END DO
!
      faces_match = .true.
!
!     -------------------------
!     determine the orientation
!     -------------------------
!
      is1 = side_map(c(1),c(2))
      is2 = side_map(c(1),c(4))
      IF (is1 == ILLEG .OR. is2 == ILLEG )     THEN
         write(6,*) 'Illegal node mapping'
         write(6,*)  target_nodes
         write(6,*)  test_nodes
         STOP 'Illegal orientation of element faces'
      END IF
      orientation = orient_map(is1,is2)
      IF(orientation == ILLEG )     THEN
         write(6,*) 'Illegal orientation of element faces'
         write(6,*)  target_nodes
         write(6,*)  test_nodes
         STOP 'Illegal orientation of element faces'
      END IF

      RETURN
      END FUNCTION faces_match
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine is to count the number of faces in the list
      SUBROUTINE count_faces(head,number_of_faces)
!
!     count the number of faces in the list
!
      TYPE (face_list), POINTER :: the_face,head
      INTEGER                   :: number_of_faces
!-----
!
      number_of_faces = 0
      IF(.not.ASSOCIATED(head))     THEN
         number_of_faces = 0
      ELSE
         the_face => head
         DO
            number_of_faces = number_of_faces + 1
            the_face => the_face%next
            if(.not.ASSOCIATED(the_face))     EXIT
         END DO
      END IF

      RETURN
      END SUBROUTINE count_faces
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine is to gather boundary faces in a list
!!
!! this routine takes an face list and creates from it a list of
!! faces that occur along external boundaries. These faces may
!! occur on different boundaries and will not necessarily be
!! in any kind of order.
      SUBROUTINE gather_boundary_faces(boundary_list_head,boundary_list_tail,&
                                       face_list_head)
!
!     this routine takes an face list and creates from it a list of
!     faces that occur along external boundaries. These faces may
!     occur on different boundaries and will not necessarily be
!     in any kind of order
!
!     faces are determined to be on a boundary if the from_domain(2)
!     field = 0.
!
      TYPE (face_list), POINTER :: boundary_list_head,boundary_list_tail,&
                                   face_list_head,new_face,face
!------
!
!     initialize the new list
!
      CALL initialize_face_list(boundary_list_head,boundary_list_tail)
!
!     add faces to the new list
!
      face => face_list_head
      DO
         IF (face%face_data%from_element(2) == 0)      THEN

            ALLOCATE(new_face)

            IF(.NOT.ASSOCIATED(boundary_list_head))    THEN
               new_face%face_data => face%face_data
               CALL initialize_first_face(new_face,boundary_list_head,boundary_list_tail)
            ELSE
               CALL append_face(new_face,boundary_list_tail)
               new_face%face_data => face%face_data
            END IF

         END IF
         face => face%next
         IF(.NOT.ASSOCIATED(face))     EXIT
      END DO
!
      RETURN   
      END SUBROUTINE gather_boundary_faces
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine is used to print a face list from first to last
      SUBROUTINE print_face_list(ioutunit,head)
!
!     print a list
!
      TYPE (face_list), POINTER :: the_face,head
!-----
!
      IF(.not.ASSOCIATED(head))     THEN
         WRITE(ioutunit,*) 'the list is empty'
      ELSE
         the_face => head
         DO
            WRITE(ioutunit,*) the_face%face_data%key,the_face%face_data%node,&
                              the_face%face_data%from_element,&
                              the_face%face_data%from_side
            the_face => the_face%next
            if(.not.ASSOCIATED(the_face))     EXIT
         END DO
      END IF

      RETURN
      END SUBROUTINE print_face_list
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine is used to print a face list from last to first
      SUBROUTINE print_face_list_backwards(ioutunit,tail)
!
!     print a list from last to first
!
      TYPE (face_list), POINTER :: the_face,tail
!-----
!
      IF(.not.ASSOCIATED(tail))     THEN
         WRITE(ioutunit,*) 'the list is empty'
      ELSE
         the_face => tail
         DO
            WRITE(ioutunit,*) the_face%face_data%key,the_face%face_data%node,&
                              the_face%face_data%from_element,&
                              the_face%face_data%from_side
            the_face => the_face%previous
            if(.not.ASSOCIATED(the_face))     EXIT
         END DO
      END IF

      RETURN
      END SUBROUTINE print_face_list_backwards
!
!///////////////////////////////////////////////////////////////////////
!> @brief This subroutine is used to free up the memory taken by a list
      SUBROUTINE deallocate_face_list(head)
!
!     free up the memory taken by a list
!
      TYPE (face_list), POINTER :: the_face,next_face,head
!-----
!
      IF(.not.ASSOCIATED(head))     THEN
         WRITE(ioutunit,*) 'the list is empty'
      ELSE
         the_face => head
         DO
            next_face => the_face%next
            DEALLOCATE(the_face)
            the_face => next_face
            if(.not.ASSOCIATED(the_face))     EXIT
         END DO
      END IF

      RETURN
      END SUBROUTINE deallocate_face_list
!
!///////////////////////////////////////////////////////////////////////
!
      END MODULE face_operations
