!//////////////////////////////////////////////////////////////////////////////////////
!
         MODULE old_size
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
         INTEGER, PARAMETER :: ngp0  = 10016! number of elements
         INTEGER, PARAMETER :: nx0   = 4  , ny0 = 4 , nz0 = 4  ! dimensions in each direction
         INTEGER, PARAMETER :: nmax0 = 4                   ! max of the above
!        
         INTEGER, PARAMETER :: nref =       0              ! dimension refinement
!
!        ---------------------------------------------------------------
!        The following defines the number of different polynomial orders
!        allowed in the mesh. This is limited to save on matrix storage.
!        ---------------------------------------------------------------
!
         INTEGER, PARAMETER :: max_orders = 3              ! number of different orders
         INTEGER, PARAMETER :: max_mrtr_mat = max_orders*max_orders/2 + 1

         INTEGER, PARAMETER :: nmp        = ngp0*6          ! number of mortars
         INTEGER, PARAMETER :: max_nodes  = 20*ngp0         ! number of element nodes
         INTEGER, PARAMETER :: max_surf   = 6              ! number of surfaces 
                                                           ! that define physical boundary
!   
!        --------------------------------------
!        Parameters associated with statistics
!        --------------------------------------
!
         INTEGER, PARAMETER :: nav     = 4                 ! neq+nav, number of averages
         INTEGER, PARAMETER :: nfluct  =  43               ! number of fluct. quant.
         
!
         INTEGER, PARAMETER :: npart = 1000                ! number of particles

      END MODULE old_size
!
!     //////////////////////////////////////////////////////////////////////////////////////
!
      MODULE old_domain_definition
      USE old_size
      SAVE
!
      TYPE old_domain
!
!       INTEGER :: neq = 5
!       INTEGER :: ngp = 2500
!       INTEGER :: n = 6
!
         INTEGER                                      :: type            ! (8 = hex8,... number of nodes in hex)
         INTEGER, DIMENSION(26)                       :: node
         INTEGER, DIMENSION(6)                        :: orientation 
         INTEGER, DIMENSION(2,6)                      :: mortar          ! mortar # and side for each face
         INTEGER, DIMENSION(6)                        :: which_surface   ! surface a face is attached to

         INTEGER, DIMENSION(3)                        :: nc,ncg          ! grid size

         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: Q               ! dependent variables
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: Qlgg            ! dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: Qglg            ! dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: Qggl            ! dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: f,g,h           ! fluxes
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0,neq)    :: g_Q             ! rk arrays
!
!        material properties are stored bu material id
!
         INTEGER                                 :: material_id
!
!        geometry arrays (gauss points)
!
         DOUBLE PRECISION, DIMENSION(nmax0,3)       :: cxg
         DOUBLE PRECISION, DIMENSION(3,nx0,ny0,nz0)   :: xg
         DOUBLE PRECISION, DIMENSION(nx0,ny0,nz0)     :: jacob ! equals 1/J for speed
         DOUBLE PRECISION, DIMENSION(3,3,nx0,ny0,nz0) :: gmet  ! metric/J
         DOUBLE PRECISION, DIMENSION(3,3,nx0,ny0,nz0) :: gmetg ! metric/J
!
!                        (lobatto points)
!
         DOUBLE PRECISION, DIMENSION(3,8)    :: corner
         DOUBLE PRECISION, DIMENSION(nmax0,3) :: cx
!
!        matrices
!
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmx      ! x derivative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmy      ! y dervative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmz      ! y dervative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: bx       ! x prolongation matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: by       ! y prolongation matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: bz       ! y prolongation matrix
!
!        face arrays
!
         INTEGER, DIMENSION(6)               :: ibtype
         CHARACTER(LEN = 9) , DIMENSION(6)   :: bcond

         CHARACTER(LEN=8)                    :: domain_type
!
      END TYPE
!
      END MODULE old_domain_definition
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE read_soution(ngrid,d)
!
!......................................................................
!     date: 3/28/2014
!
!     read in restart file interpolate the solution Q to the new polynomial order
!......................................................................
!
        USE domain_definition
        USE old_domain_definition
!       USE mortar_definition
!       USE stats_definition
!       USE particle_definition
!       USE Order_Matrices
!       USE Material_Properties
!       USE User_Data
!       USE input_data
!       USE mpi_par
!       USE constants
!       USE physics
!       USE mpi
!
      IMPLICIT none
!      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
       TYPE (domain)        :: d(ngp)
       TYPE (old_domain)    :: da

       INTEGER          :: read_unit,ngrid,nmort,id,i,j,k,l,m,n
       DOUBLE PRECISION :: hx,sum1,sum2,sum3,sum4,sum5,xmapped,ymapped,zmapped
       DOUBLE PRECISION :: polyn,hy(ny),hz(nz)
       INTEGER          :: num_ordersx0,num_ordersy0,num_ordersz0
       INTEGER          :: order_mapx0,order_mapy0,order_mapz0
       DOUBLE PRECISION, DIMENSION(nx0,nx0,3) :: bx0,dmx0
       DOUBLE PRECISION, DIMENSION(ny0,ny0,3) :: by0,dmy0
       DOUBLE PRECISION, DIMENSION(nz0,nz0,3) :: bz0,dmz0
       INTEGER :: ngrid0,nmort0,num_mrtr_matrices0
       INTEGER, DIMENSION(2,3*3/2+1) :: mrtr_mat_map0
       DOUBLE PRECISION, DIMENSION(nmax0,nmax0,3*3/2+1) :: prolong0,restrict0
!
!
	  write(*,*) "Reading from solution.rst ..."

      read_unit = 61
      OPEN (unit=read_unit,file='solution.rst',status='old',form='unformatted')
!      
!     -----------------
!     read in run sizes
!     -----------------
!
      READ(read_unit) ngrid0,nmort0
!
!     -------------------------
!     read in computed matrices
!     -------------------------
!
      READ(read_unit) num_ordersx0,num_ordersy0,num_ordersz0
      READ(read_unit) order_mapx0,order_mapy0,order_mapz0
      READ(read_unit) bx0
      READ(read_unit) by0
      READ(read_unit) bz0
      READ(read_unit) dmx0
      READ(read_unit) dmy0
      READ(read_unit) dmz0

      READ(read_unit) num_mrtr_matrices0
      READ(read_unit) mrtr_mat_map0
      READ(read_unit) prolong0
      READ(read_unit) restrict0
!
!
      DO id = 1,ngrid
        
        !------- Read the domain
          
          CALL zero_domain(da)
          READ(read_unit) da%type            ! (8 = hex8,... number of nodes in hex)
          READ(read_unit) da%node
          READ(read_unit) da%orientation 
          READ(read_unit) da%mortar          ! mortar # and side for each face
          READ(read_unit) da%which_surface   ! surface a face is attached to

          READ(read_unit) da%nc,da%ncg       ! grid size

          READ(read_unit) da%Q               ! dependent variables
          READ(read_unit) da%Qlgg            ! dependent variables
          READ(read_unit) da%Qglg            ! dependent variables
          READ(read_unit) da%Qggl            ! dependent variables
          READ(read_unit) da%f               ! fluxes
          READ(read_unit) da%g               ! fluxes
          READ(read_unit) da%h               ! fluxes
          READ(read_unit) da%g_Q             ! rk arrays
!
          READ(read_unit) da%material_id
!
          READ(read_unit) da%cx,da%cxg
          READ(read_unit) da%xg
          READ(read_unit) da%jacob
          READ(read_unit) da%gmet
          READ(read_unit) da%gmetg
!
          READ(read_unit) da%corner
!
          READ(read_unit) da%ibtype
          READ(read_unit) da%bcond
          READ(read_unit) da%domain_type
          
          !------- Interpolate the domain
          
          DO k = 1,d(id)%ncg(3)
            DO j = 1,d(id)%ncg(2)
              DO i = 1,d(id)%ncg(1)
              
              xmapped = d(id)%cxg(i,1)
              ymapped = d(id)%cxg(j,2)
              zmapped = d(id)%cxg(k,3)
              
              ! Compute Lagrangian interpolating polynomial in 2 directions
              
              hy = 0.0d0
              DO m = 1,da%ncg(2)
                hy(m) = polyn(m,ymapped,da%ncg(2),da%cxg(1:da%ncg(2),1))
              END DO
              
              hz = 0.0d0
              DO l = 1,da%ncg(3)
                hz(l) = polyn(l,zmapped,da%ncg(3),da%cxg(1:da%ncg(3),1))
              END DO
              
              !
              ! Interpolate fluid propterties to the new grid position with spectral interpolation
              !
              
              sum1 = 0.0d0
              sum2 = 0.0d0
              sum3 = 0.0d0
              sum4 = 0.0d0
              sum5 = 0.0d0
              
              do m = 1,da%ncg(1)
                hx = polyn(m,xmapped,da%ncg(1),da%cxg(1:da%ncg(1),1))
                do n = 1,da%ncg(2)
                  do l = 1,da%ncg(3)
                    sum1 = sum1 + da%Q(m,n,l,1)*hx*hy(n)*hz(l)
                    sum2 = sum2 + da%Q(m,n,l,2)*hx*hy(n)*hz(l)
                    sum3 = sum3 + da%Q(m,n,l,3)*hx*hy(n)*hz(l)
                    sum4 = sum4 + da%Q(m,n,l,4)*hx*hy(n)*hz(l)
                    sum5 = sum5 + da%Q(m,n,l,5)*hx*hy(n)*hz(l)
                  end do
                end do
              end do
              
              d(id)%Q(i,j,k,1) = sum1
              d(id)%Q(i,j,k,2) = sum2
              d(id)%Q(i,j,k,3) = sum3
              d(id)%Q(i,j,k,4) = sum4
              d(id)%Q(i,j,k,5) = sum5
              
              END DO
            END DO
          END DO
            
        END DO
!
        CLOSE(read_unit)
!
      RETURN
      END SUBROUTINE
!

