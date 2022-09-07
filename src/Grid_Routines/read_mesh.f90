!> @file 
!! This file contains the routines necessary to read the mesh file
!! in order to get the information such as nodes, elements, and boundary
!! conditions.
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!////////											////////
!////////	read_mesh.f90									////////
!////////											////////
!////////	contains:									////////
!////////											////////
!////////	      SUBROUTINE read_mesh_file(iunit,d,num_elements,nodes)			////////
!////////	      SUBROUTINE read_specification_file(iunit,num_surf,surfaces)		////////
!////////	      SUBROUTINE read_neutral_file(iunit,d,num_elements,nodes)			////////
!////////	      SUBROUTINE read_p_header_card(iunit,it,id,iv,kc,n1,n2,n3,n4,n5)		////////
!////////	      SUBROUTINE read_p_summary(iunit,file_date,file_time,version)		////////
!////////	      SUBROUTINE read_p_node_data_c1(iunit,x,y,z)				////////
!////////	      SUBROUTINE read_p_node_data_c2(iunit,icf,gtype,ndf,config,cid,pspc)	////////
!////////	      SUBROUTINE read_p_element_data_c1(iunit,nodes,config,pid,ceid,th1,th2,th3)////////
!////////	      SUBROUTINE read_p_element_data_c2(iunit,lnodes)				////////
!////////	      SUBROUTINE read_p_element_data_c3(iunit,adata)				////////
!////////	      SUBROUTINE read_nodes(iunit,nodes)					////////
!////////	      SUBROUTINE read_elements(iunit,d,num_elements)				////////
!////////	      SUBROUTINE read_which_boundary_ids(iunit,d)				////////
!////////	      SUBROUTINE write_err_msg(msg_n,input,maximum,where)			////////
!////////											////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief
!> Store keywords needed to navigate input file
!!
!
      MODULE global_infile_keywords ! store keywords needed to navigate input file
         SAVE
         INTEGER                           :: num_keywords = 7
         CHARACTER (LEN=10), DIMENSION(7)  :: keyword_list = &
                                            (/ 'surfaces  ', &
                                               'attachment', &
                                               'nodes     ', &
                                               'patran    ', &
                                               'elements  ', &
                                               'simple    ', &
                                               'boundary  '  &
                                            /)
      END MODULE global_infile_keywords
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read simple mesh 1.1 file for getting the information
!> of nodes, elements, and boundaries.
!!
!
     SUBROUTINE read_mesh_file(iunit,d,num_elements,nodes)

!     .......................................................................
!
!     read the mesh file. For the format, see the individual read routines
!
!     .......................................................................

      USE FE_Data_Types
      USE Domain_Definition
      USE global_infile_keywords
      USE File_Units
!
      TYPE (domain),DIMENSION(ngp)  :: d
      TYPE (node_vector)            :: nodes
!
      CHARACTER (LEN = 132) :: input_line
      CHARACTER (LEN = 132) :: get_string_value
      CHARACTER (LEN = 10)  :: keyword_in_, the_word
      LOGICAL               :: eof
      CHARACTER (LEN = 3)   :: version_no
!
!
      version_no = 'nil'
      DO
         CALL skip_comments(input_line,iunit,eof)
            IF(eof)     EXIT
            the_word = keyword_in_(input_line,keyword_list,num_keywords,kword)

            SELECT CASE(the_word)

               CASE('patran')
                  CALL read_neutral_file(iunit,d,num_elements,nodes)

               CASE('nodes')
                  CALL read_nodes(iunit,nodes)

               CASE('elements')
                  CALL read_elements(iunit,d,num_elements)

               CASE('boundary')
                  CALL read_which_boundary_ids(iunit,d)

               CASE('simple')
                  version_no = get_string_value(input_line)

               CASE DEFAULT
                  WRITE(err_unit,*) 'unrecognized keyword in input file: ',the_word
                  WRITE(err_unit,*) input_line
         END SELECT
      END DO
!
      RETURN
!
      END SUBROUTINE read_mesh_file
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read specific mesh file
!!
!
      SUBROUTINE read_specification_file(iunit,num_surf,surfaces)

!     .......................................................................
!     read the surfaces from the input file
!     .......................................................................

      USE FE_Data_Types
      USE Material_Properties
      USE File_Units

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (surface), DIMENSION(max_surf)       :: surfaces
!
      CHARACTER (LEN = 132) :: input_line
      CHARACTER (LEN = 132) :: get_string_value
      INTEGER               :: get_int_value
      LOGICAL               :: eof
!
!     ------------------------------------
!     read in surface/boundary information
!     ------------------------------------
!
      n = 0
      DO     ! until an "end all" line is read

         CALL skip_comments(input_line,iunit,eof)
         IF(eof)     EXIT
         IF(input_line(1:7) == 'end all')   EXIT

         n = n + 1
         if( n > max_surf ) call write_err_msg(8,n,max_surf,0)

         surfaces(n)%surface_id   = get_int_value(input_line)
         READ(iunit,'(a132)')       input_line
         surfaces(n)%bcond        = get_string_value(input_line)
         READ(iunit,'(a132)')       input_line
         surfaces(n)%surface_type = get_string_value(input_line)

         SELECT CASE(surfaces(n)%surface_type) ! see what type of surface

            CASE ('equation') ! define surface by a parametric equation

               DO k = 1,3     ! read in equations for (x,y,z)
                  READ(iunit,'(a132)') input_line
                  surfaces(n)%equation(k) = TRIM(input_line)
               END DO
               READ(iunit,'(a132)') input_line   ! skip over the end surface ... line

            CASE ('element  ') ! define surface from the element data

               READ(iunit,'(a132)') input_line   ! skip over the end surface ... line
!
            CASE DEFAULT

               WRITE(err_unit,*) 'surface type ',surfaces(n)%surface_type,' not supported'
               STOP 'unsupported surface type found in geometry file'

         END SELECT
      END DO
      num_surf = n
!
!     ------------------------------------
!     read in material information
!     ------------------------------------
!
      num_prop = 0
      CALL skip_comments(input_line,iunit,eof)
      IF(eof)     RETURN
      num_prop = get_int_value(input_line)
      DO
         CALL skip_comments(input_line,iunit,eof)
         IF(eof)     EXIT
         IF(input_line(1:7) == 'end all')   EXIT

         mid = get_int_value(input_line)
         IF (mid > max_materials )   call write_err_msg(11,mid,max_materials,0)  
         READ(iunit,*) (material_property(j,mid), j = 1,num_prop)
      END DO

      RETURN
      END SUBROUTINE read_specification_file
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to  read a patran neutral file and store 
!> information needed by the code currently supports packets 1,2,26, and 99.
!!
!
      SUBROUTINE read_neutral_file(iunit,d,num_elements,nodes)
!
      USE Method_Data
      USE Domain_Definition
      USE FE_Data_Types
      IMPLICIT none
!
      TYPE (domain),DIMENSION(ngp) :: d !< list of elements
      TYPE (node_vector)           :: nodes    !< list of nodes
!
      INTEGER               :: iunit,it,id,iv,kc,n1,n2,n3,n4,n5
      INTEGER               :: icf,ndf,config,cid,pspc(6)
      INTEGER               :: nnodes,pid,ceid,lnodes(10)
      INTEGER               :: offset_node,offset_element,k,j,num_elements
      CHARACTER(LEN = 4)    :: file_date(3), file_time(2), version(3)
      CHARACTER (LEN = 132) :: input_line
      CHARACTER(LEN=1)      :: gtype
      DOUBLE PRECISION      :: x,y,z
      DOUBLE PRECISION      :: th1,th2,th3
      LOGICAL               :: first_node, first_element

!-----
      first_node     = .true.
      first_element  = .true.
      offset_node    = 0
      offset_element = 0
!
      DO
         CALL read_p_header_card(iunit,it,id,iv,kc,n1,n2,n3,n4,n5)

         SELECT CASE(it)

            CASE(1)                                  ! node cards

               CALL read_p_node_data_c1(iunit,x,y,z)
               IF(first_node)     THEN
                  offset_node = id - 1
                  first_node  = .false.
               END IF
               nodes%node(id - offset_node)%x(1) = x
               nodes%node(id - offset_node)%x(2) = y
               nodes%node(id - offset_node)%x(3) = z
               CALL read_p_node_data_c2(iunit,icf,gtype,ndf,config,cid,pspc)

            CASE(2)                                  ! element cards

               IF(first_element)     THEN
                  offset_element = id - 1
                  first_element = .false.
               END IF
               CALL read_p_element_data_c1(iunit,nnodes,config,pid,ceid,th1,th2,th3)
               CALL read_p_element_data_c2(iunit,lnodes)
               DO k = 1,nnodes
                  d(id - offset_element)%node(k) = lnodes(k) - offset_node
               END DO
               d(id - offset_element)%type  = iv

               IF(method == 'dg    ')     THEN
                  d(id - offset_element)%nc(1) = n3 + 1
                  d(id - offset_element)%nc(2) = n4 + 1
                  d(id - offset_element)%nc(3) = n5 + 1
                  d(id - offset_element)%ncg(1) = n3 + 1
                  d(id - offset_element)%ncg(2) = n4 + 1
                  d(id - offset_element)%ncg(3) = n5 + 1
               ELSE
                  d(id - offset_element)%nc(1) = n3 + 2
                  d(id - offset_element)%nc(2) = n4 + 2
                  d(id - offset_element)%nc(3) = n5 + 2
                  d(id - offset_element)%ncg(1) = n3 + 1
                  d(id - offset_element)%ncg(2) = n4 + 1
                  d(id - offset_element)%ncg(3) = n5 + 1
               END IF

               DO j = 1,n1                           ! ignore adata associated values
                  READ(iunit,*) input_line
               END DO

            CASE(26)                                 ! summary cards

               nodes%num_nodes       = n1
               num_elements = n2
               IF(n1 > max_nodes)     CALL write_err_msg(1,n1,max_nodes,0)
               IF(n2 > ngp)  CALL write_err_msg(2,n2,ngp,0)
               CALL read_p_summary(iunit,file_date,file_time,version)

            CASE(99)  ! end card

               RETURN

            CASE DEFAULT                             ! ignore all other information

               DO j = 1,kc
                  READ(iunit,*) input_line
               END DO

         END SELECT
      END DO

      RETURN

      END SUBROUTINE read_neutral_file
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to  read a patran packet header card
!!
!
      SUBROUTINE read_p_header_card(iunit,it,id,iv,kc,n1,n2,n3,n4,n5)

!     .........................................................................
!     read a patran packet header card
!     .........................................................................

      INTEGER, INTENT(OUT) :: iunit,it,id,iv,kc,n1,n2,n3,n4,n5
!-----
      READ(iunit,fmt='(i2,8i8)') it,id,iv,kc,n1,n2,n3,n4,n5
!
      RETURN
      END SUBROUTINE read_p_header_card
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to  read a patran packet summary card (packet 26)
!!
!
      SUBROUTINE read_p_summary(iunit,file_date,file_time,version)


      CHARACTER(LEN = 4) :: file_date(3), file_time(2), version(3)
!----
      READ(iunit,fmt='(3a4,2a4,3a4)') file_date,file_time,version
!
      RETURN
      END SUBROUTINE read_p_summary
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read a patran packet node card (packet 1, card 1)
!!
!
      SUBROUTINE read_p_node_data_c1(iunit,x,y,z)


      DOUBLE PRECISION :: x,y,z
!----
      READ(iunit,fmt = '(3e16.9)') x,y,z
!
      RETURN
      END SUBROUTINE read_p_node_data_c1
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read a patran packet node card (packet 1, card 2)
!!
!
      SUBROUTINE read_p_node_data_c2(iunit,icf,gtype,ndf,config,cid,pspc)



      INTEGER          :: icf,ndf,config,cid,pspc(6)
      CHARACTER(LEN=1) :: gtype
!----
      READ(iunit,fmt = '(i1,1a1,i8,i8,i8,2x,6i1)') icf,gype,ndf,config,cid,&
                                                (pspc(j),j=1,6)
!
      RETURN
      END SUBROUTINE read_p_node_data_c2
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read a patran packet element card (packet 2, card 1)
!!
!
      SUBROUTINE read_p_element_data_c1(iunit,nodes,config,pid,ceid,th1,th2,th3)


      INTEGER          :: iunit,nodes,config,pid,ceid
      DOUBLE PRECISION :: th1,th2,th3
!----
      READ(iunit,fmt='(i8,i8,i8,i8,3e16.9)') nodes,config,pid,ceid,th1,th2,th3
!
      RETURN
      END SUBROUTINE read_p_element_data_c1
!
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read a patran packet element card (packet 2, card 2)
!!
!
      SUBROUTINE read_p_element_data_c2(iunit,lnodes)

!     .........................................................................
!     read a patran packet element card (packet 2, card 2)
!     .........................................................................

      INTEGER          :: iunit,lnodes(10)
!-----
      READ(iunit,fmt='(10i8)') (lnodes(j),j=1,4)
!
      RETURN
      END SUBROUTINE read_p_element_data_c2
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read a patran packet element card (packet 2, card 3)
!!
!
      SUBROUTINE read_p_element_data_c3(iunit,adata)

!     .........................................................................
!     read a patran packet element card (packet 2, card 3)
!     .........................................................................

      INTEGER                        :: iunit
      DOUBLE PRECISION, DIMENSION(5) :: adata
!-----
      READ(iunit,fmt='(5e16.9)') (adata(j),j=1,5)
!
      RETURN
      END SUBROUTINE read_p_element_data_c3
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read the nodes from the input file in simple node list format
!!
!
      SUBROUTINE read_nodes(iunit,nodes)

!
      USE FE_Data_Types
      TYPE (node_vector) :: nodes
!
      CHARACTER (LEN = 132) :: input_line
      DOUBLE PRECISION      :: z
!-----
      j = 0
      DO     ! until an "end nodes" line is read

         READ(iunit,'(a132)') input_line
         IF(input_line(1:3) == 'end')   EXIT
   
         j = j + 1
         if(j > max_nodes) call write_err_msg(1,j,max_nodes,0)
         read(input_line,*) node_n,(nodes%node(node_n)%x(k),k=1,3)
		 !write(*,*) node_n,nodes%node(node_n)%x(1)

      END DO
      nodes%num_nodes = j

      RETURN
      END SUBROUTINE read_nodes
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read the elements from the input file in the simple element format
!!
!
      SUBROUTINE read_elements(iunit,d,num_elements)


!
      USE Method_Data
      USE Domain_Definition
      USE FE_Data_Types
      USE Domain_Definition

      TYPE (domain),DIMENSION(ngp) :: d
      CHARACTER (LEN = 256)        :: input_line
      INTEGER                      :: element_id,etype
      CHARACTER (LEN = 4)          :: topology
!-----

      j = 0
      DO     ! until an "end elements" line is read

         READ(iunit,'(a256)') input_line
         IF(input_line(1:3) == 'end')   EXIT
   
         j = j + 1
         if(j > ngp) call write_err_msg(2,j,max_nodes,0)

         read(input_line,*) element_id,topology,etype,&
                            d(element_id)%material_id,&
                           (d(element_id)%nc(k),k=1,3)

         read(iunit,*) (d(element_id)%node(k),k=1,etype)
         d(element_id)%type = etype

         IF( method == 'dg    ' )     THEN
            DO k = 1,3
               d(element_id)%nc(k)   = d(element_id)%nc(k) + 1
               d(element_id)%ncg(k)  = d(element_id)%nc(k)
            END DO
         ELSE
            DO k = 1,3
               d(element_id)%ncg(k)  = d(element_id)%nc(k) + 1
               d(element_id)%nc(k)   = d(element_id)%ncg(k) + 1
            END DO
         END IF

      END DO
      num_elements = j

      RETURN
      END SUBROUTINE read_elements
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to read the surface from which the element 
!> face boundary condition comes from when the input file is in the 
!> simple element format. 
!!
!
      SUBROUTINE read_which_boundary_ids(iunit,d)

!
      USE FE_Data_Types
      USE Domain_Definition

      TYPE (domain),DIMENSION(ngp) :: d
      INTEGER                      :: element_id,face_id
!
      CHARACTER (LEN = 132) :: input_line
!-----
      DO id = 1,ngp
         DO face_id = 1,6
            d(id)%which_surface(face_id) = 0
         END DO
      END DO

      DO     ! until an "end boundary" line is read

         READ(iunit,'(a132)') input_line
         IF(input_line(1:3) == 'end')   EXIT
   
         READ(input_line,*) element_id,face_id, &
         d(element_id)%which_surface(face_id)

      END DO

      RETURN
      END SUBROUTINE read_which_boundary_ids
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine is used to write out an error message and stop the 
!> computation if the input is bigger than the maximum.
!!
!
      SUBROUTINE write_err_msg(msg_n,input,maximum,where)
!
!
      USE File_Units
      INTEGER,INTENT(in) :: msg_n,input,maximum,where
!
!---
      SELECT CASE(msg_n)
        CASE(1)
           WRITE(err_unit,*) 'specification file has too many nodes!'
        CASE(2)
           WRITE(err_unit,*) 'specification file has too many elements/subdomains!' 
        CASE(3)
           WRITE(err_unit,*) 'specification ',where,' has too many points!' 
        CASE(4)
           WRITE(err_unit,*) 'specification file has too many surfaces!'
        CASE(5)
           WRITE(err_unit,*) 'specification file has too many nodes along surface ',where
        case(6)
           WRITE(err_unit,*) 'specification file has too many nodes with boundary condition ', where
        CASE(7)
           WRITE(err_unit,*) 'Too many knots for spline on surface ',where
        CASE(8)
           WRITE(err_unit,*) 'specification file has too many faces!'
        CASE(9)
           WRITE(err_unit,*) 'too many points along surface ', where
        CASE(10)
           WRITE(err_unit,*) 'specification file has too many different boundary conditions!'
        CASE(11)
           WRITE(err_unit,*) 'specification file has too many materials!'
      END SELECT
!
      WRITE(err_unit,*) 'code currently dimensioned for ',maximum 
!     WRITE(err_unit,*) 'the data file has ',input
      WRITE(err_unit,*) 'make appropriate changes in size.f and recompile'
      STOP       'An input parameter is larger than dimensioned size'
!
      END SUBROUTINE write_err_msg
