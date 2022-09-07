      MODULE size
         SAVE
!
!        --------------------------------------
!        Parameters associated with the physics
!        --------------------------------------
!
         INTEGER, PARAMETER :: neq = 5            ! number of equations in system
         INTEGER, PARAMETER :: max_materials = 1  ! number of materials
         INTEGER, PARAMETER :: mxprops = 1        ! number of material properties
!
!        -----------------------------------------
!        Parameters associated with the grid size
!        ----------------------------------------
!
         INTEGER, PARAMETER          :: ngp  = 10648                  ! number of elements
         INTEGER, PARAMETER          :: nx   = 10, ny = 10, nz = 10      ! dimensions in each direction
         INTEGER, PARAMETER          :: nmax = 10                           ! max of the above
        
         INTEGER, PARAMETER          :: nref = 0                           ! dimension refinement
!
!        ---------------------------------------------------------------
!        The following defines the number of different polynomial orders
!        allowed in the mesh. This is limited to save on matrix storage.
!        ---------------------------------------------------------------
!
         INTEGER, PARAMETER :: max_orders = 3              ! number of different orders
         INTEGER, PARAMETER :: max_mrtr_mat = max_orders*max_orders/2 + 1

         INTEGER, PARAMETER :: nmp        = ngp*6          ! number of mortars
         INTEGER, PARAMETER :: max_nodes  = 20*ngp         ! number of element nodes
         INTEGER, PARAMETER :: max_surf   = 6              ! number of surfaces 
                                                           ! that define physical boundary
!   
!        --------------------------------------
!        Parameters associated with statistics
!        --------------------------------------
!
         INTEGER, PARAMETER :: nav     = 4!4                 ! neq+nav, number of averages
         INTEGER, PARAMETER :: nfluct  =  43               ! number of fluct. quant.
         
!
         INTEGER, PARAMETER :: npart  = 1                ! number of particles
         INTEGER, PARAMETER ::  nfour = 4                !/////////////////

      END MODULE size
