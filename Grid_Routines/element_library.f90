!> \file element_library.f90 Contains routines and data relating to standard finite elements
!
!> @brief Shared data storage for hexahedral elements
!!
!! Contains face and interpolation mappings for HEX8, HEX26, and HEX64 elements
      MODULE Element_Library


      DOUBLE PRECISION, PARAMETER :: oth = 1.d0/3.d0, tth = 2.d0/3.d0
!
!    ------------------------------------------------------------------
!     mapping of element ordered nodes to the 6 faces of a HEX8
!     mapping of the nodes on a face for the face interpolant of a HEX8
!    ------------------------------------------------------------------
!
      INTEGER, DIMENSION(4,6) :: Face_Map_hex8 = &
      RESHAPE((/1,2,6,5, &
                4,3,7,8, &
                1,2,3,4, &
                2,3,7,6, &
                5,6,7,8, &
                1,4,8,5/),(/4,6/)) !< mapping of element ordered nodes to the 6 faces of a HEX8
!
      INTEGER, DIMENSION(2,2,6) :: Intrp_Map_Hex8 = &
      RESHAPE((/1,2,5,6, &
                4,3,8,7, &
                1,2,4,3, &
                2,3,6,7, &
                5,6,8,7, &
                1,4,5,8/),(/2,2,6/)) !< mapping of the nodes on a face for the face interpolant of a HEX8
!
!    ------------------------------------------------------------------
!     mapping of element ordered nodes to the 6 faces of a HEX26
!     mapping of the nodes on a face for the face interpolant of a HEX26
!    ------------------------------------------------------------------
!
      INTEGER, DIMENSION(9,6) :: Face_Map_hex26 = &
      RESHAPE((/1,2,6,5,9,14,17,13,25,            &
                4,3,7,8,11,15,19,16,26,           &
                1,2,3,4,9,10,11,12,21,            &
                2,3,7,6,10,15,18,14,24,           &
                5,6,7,8,17,18,19,20,22,           &
                1,4,8,5,12,16,20,13,23/),(/9,6/)) !< mapping of element ordered nodes to the 6 faces of a HEX26

      INTEGER, DIMENSION(3,3,6) :: Intrp_Map_Hex26 = &
      RESHAPE((/1,9,2,13,25,14,5,17,6,               &
                4,11,3,16,26,15,8,19,7,              &
                1,9,2,12,21,10,4,11,3,               &
                2,10,3,14,24,15,6,18,7,              &
                5,17,6,20,22,18,8,19,7,              &
                1,12,4,13,23,16,5,20,8/),(/3,3,6/)) !< mapping of the nodes on a face for the face interpolant of a HEX26
!
!    ------------------------------------------------------------------
!     mapping of element ordered nodes to the 6 faces of a HEX64
!     mapping of the nodes on a face for the face interpolant of a HEX64
!    ------------------------------------------------------------------
!

      INTEGER, DIMENSION(16,6) :: Face_Map_hex64 = &
      RESHAPE((/1,2,6,5,9,10,18,22,26,25,21,17,0,0,0,0,    &
                4,3,7,8,14,13,19,23,29,30,24,20,0,0,0,0,   &
                1,2,3,4,9,10,11,12,13,14,15,16,0,0,0,0,    &
                2,3,7,6,11,12,19,23,28,27,22,18,0,0,0,0,   &
                5,6,7,8,25,26,27,28,29,30,31,32,0,0,0,0,   &
                1,4,8,5,16,15,20,24,31,32,21,17,0,0,0,0/),(/16,6/)) !< mapping of element ordered nodes to the 6 faces of a HEX64

      INTEGER, DIMENSION(4,4,6) :: Intrp_Map_Hex64 = &
      RESHAPE((/1,2,6,5,9,10,18,22,26,25,21,17,0,0,0,0,    &
                4,3,7,8,14,13,19,23,29,30,24,20,0,0,0,0,   &
                1,2,3,4,9,10,11,12,13,14,15,16,0,0,0,0,    &
                2,3,7,6,11,12,19,23,28,27,22,18,0,0,0,0,   &
                5,6,7,8,25,26,27,28,29,30,31,32,0,0,0,0,   &
                1,4,8,5,16,15,20,24,31,32,21,17,0,0,0,0/),(/4,4,6/)) !< mapping of the nodes on a face for the face interpolant of a HEX64
!
      DOUBLE PRECISION, DIMENSION(2) :: cx_hex8  = (/0.0d0,1.0d0/)
      DOUBLE PRECISION, DIMENSION(3) :: cx_hex26 = (/0.0d0,0.5d0,1.0d0/)
      DOUBLE PRECISION, DIMENSION(4) :: cx_hex64 = (/0.0d0,oth,tth,1.0d0/)
!
      END MODULE element_library
!
!> @brief Edge mapping storage of hexahedral elements
!!
!! Orientation map: This is a 4x4 array that determines how a slave 
!! face is mapped onto a mortar. The value 1 or 2 indicates the index
!! while the sign indicates forward or backward. This array, then assigns
!! the constants defined above to the particular copy operation that must be
!! done. For instance, (index1,index2) of the mortar mates with (index2,-index1)
!! of the slace, then the orientation is assigned the value F2B1 for "forward 2,
!! backward 1".
!!
!! Side Map: Given two corner points, return the index in which an array runs
!! betweeen these two corners. This is shown in the figure below:
      MODULE Edge_Mappings
!
!     -------------------------
!     list of defined constants
!     -------------------------
!
      INTEGER, PARAMETER :: ILLEG = -99, DFLT = 0, B1F2 = 1, B1B2 = 2, &
                            F1B2  = 3  , F2F1 = 4, B2F1 = 5, B2B1 = 6, &
                            F2B1  = 7
!
!     -----------------------------------------------------------------------------
!     Orientation map: This is a 4x4 array that determines how a slave 
!     face is mapped onto a mortar. The value 1 or 2 indicates the index
!     while the sign indicates forward or backward. This array, then assigns
!     the constants defined above to the particular copy operation that must be
!     done. For instance, (index1,index2) of the mortar mates with (index2,-index1)
!     of the slace, then the orientation is assigned the value F2B1 for "forward 2,
!     backward 1".
!     -----------------------------------------------------------------------------
!
      INTEGER, PARAMETER, DIMENSION(-2:2,-2:2) :: orient_map =              &
                          RESHAPE((/ILLEG,B1B2 ,ILLEG,F1B2 ,ILLEG,          &
                                    B2B1 ,ILLEG,ILLEG,ILLEG,F2B1,           &
                                    ILLEG,ILLEG,ILLEG,ILLEG,ILLEG,          &
                                    B2F1 ,ILLEG,ILLEG,ILLEG,F2F1,           &
                                    ILLEG,B1F2 ,ILLEG,DFLT ,ILLEG/),(/5,5/)) !< Determines how a slave face is mapped onto a mortar.
!
!     -----------------------------------------------------------------------------
!     Side Map: Given two corner points, return the index in which an array runs
!     betweeen these two corners. This is shown in the figure below:
!
!               4 -----------3-------------3
!               |                          |
!          /\   |                          |
!          |    |                          |
!               |                          |
!          2    |                          |
!               4                          2
!          x    |                          |
!          e    |                          |
!          d    |                          |
!          n    |                          |
!          i    |                          |
!               |                          |
!               1 -----------1-------------2
!                        index 1 -->
!     -----------------------------------------------------------------------------
!
      INTEGER, PARAMETER, DIMENSION(4,4) :: side_map = &
                                 RESHAPE((/ILLEG, -1   , ILLEG, -2   ,          &
                                           1    , ILLEG, -2   , ILLEG,           &
                                           ILLEG, 2    , ILLEG, 1    ,           &
                                           2    , ILLEG, -1   , ILLEG/),(/4,4/)) 
      END MODULE Edge_Mappings
