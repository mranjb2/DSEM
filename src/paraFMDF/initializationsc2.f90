!> \file initializationsc2.f90
!! Variables initializations, main modules contained here
!///////////////////////////////////////////////////////////////////////
!///////                            ////////
!///////                            ////////
!///////   initializations :NS STAC3M CODE VERSION      ////////
!///////                            ////////
!///////      contains: Method_Data             ////////
!///////                constants               ////////
!///////                rk_coefs                ////////
!///////                domain_definition           ////////
!///////                boundary_definition         ////////
!///////                mortar_definition           ////////
!///////                            ////////
!///////////////////////////////////////////////////////////////////////
!> Defines numerical method
      MODULE Method_Data
         CHARACTER (LEN = 6) :: method = 'stac3m'  ! for STAC3M
      END MODULE Method_Data
!
!     //////////////////////////////////////////////////////////////////////////////////////
!> Common numerical constants
      MODULE constants
!
         DOUBLE PRECISION, PARAMETER :: pi    = 3.141592653589793238462643d0 !< \f$ \pi \f$
         DOUBLE PRECISION, PARAMETER :: twopi = 2.d0*pi !< \f$ 2 \pi \f$
         DOUBLE PRECISION            :: eps !< Convergence criteria
         DOUBLE PRECISION            :: big !< A very large machine number
!
      END MODULE constants
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Runge Kutta coefficient storage
      MODULE rk_coefs
!
         INTEGER                         :: kord !< RK order
         DOUBLE PRECISION, DIMENSION (5) :: ark,brk,crk
         DOUBLE PRECISION                :: dtconst !< timestepping constant
!
      END MODULE rk_coefs
!
!     //////////////////////////////////////////////////////////////////////////////////////
!> Main definition of grid/elements/domain
      MODULE domain_definition
         USE size
         SAVE
!
!     date 11/9/1998
!> Element datatype
      TYPE domain
!
         INTEGER                                      :: type            !< (8 = hex8,... number of nodes in hex)
         INTEGER, DIMENSION(26)                       :: node			!< nodes
         INTEGER, DIMENSION(6)                        :: orientation 	!< element orientation
         INTEGER, DIMENSION(2,6)                      :: mortar          !< mortar # and side for each face
         INTEGER, DIMENSION(6)                        :: which_surface   !< surface a face is attached to

         INTEGER, DIMENSION(3)                        :: nc,ncg          !< grid size

         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: Q,Q_old,Q_old2  !< dependent variables
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: Qlgg            !< dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: Qglg            !< dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: Qggl            !< dependent variables on faces
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: f,g,h           !< fluxes
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)    :: g_Q             !< rk arrays
!
!        material properties are stored bu material id
!
         INTEGER                                 :: material_id	!< material if set
!
!        geometry arrays (gauss points)
!
         DOUBLE PRECISION, DIMENSION(nmax,3)       :: cxg !< mappged gauss point location
         DOUBLE PRECISION, DIMENSION(3,nx,ny,nz)   :: xg !< physical gauss point location
         DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: jacob !< equals 1/J for speed
         DOUBLE PRECISION, DIMENSION(3,3,nx,ny,nz) :: gmet  !< metric/J
         DOUBLE PRECISION, DIMENSION(3,3,nx,ny,nz) :: gmetg !< metric/J
!
!                        (lobatto points)
!
         DOUBLE PRECISION, DIMENSION(3,8)    :: corner
         DOUBLE PRECISION, DIMENSION(nmax,3) :: cx
!
!        EV matrices
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: Sh,Sh_old,Sh_old2    !< Entropy and 2-level entropy history
         DOUBLE PRECISION, DIMENSION(nx,ny,nz) :: Sh_lgg,Sh_glg,Sh_ggl !< Entropy on Lobatto grids
         DOUBLE PRECISION, DIMENSION(nx,ny,nz) :: RHS_ent              !< Entropy viscosity method residual
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: fe,re,ge,se,he,te    !< Entropy fluxes
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: rho_old,rho_old2     !< Density history
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: muhg                 !< Entropy viscosity
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: hx                   !< local distance between Gauss points
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: muhlgg,muhglg,muhggl !< Entropy viscosity on Lobatto grids
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: Strain               !< filtered rate of strain
         DOUBLE PRECISION, DIMENSION(nx,nx,nx) :: ducros, rhoSensor    !< Ducros and density-based sensor

!        Smagorinsky
         DOUBLE PRECISION, DIMENSION(nx,ny,nz)        :: nu_t       !< LES viscosity on ggg points
         DOUBLE PRECISION, DIMENSION(nx,ny,nz)        :: nu_t_lgg   !< LES viscosity on lgg points
         DOUBLE PRECISION, DIMENSION(nx,ny,nz)        :: nu_t_glg   !< LES viscosity on glg points
         DOUBLE PRECISION, DIMENSION(nx,ny,nz)        :: nu_t_ggl   !< LES viscosity on ggl points
         double precision, dimension(nx,ny,nz)        :: nu_x       !< LES viscosity derivative on ggg points \f[ \frac{\partial \nu}{\partial x} \f]
         double precision, dimension(nx,ny,nz)        :: nu_y       !< LES viscosity derivative on ggg points \f[ \frac{\partial \nu}{\partial y} \f]
         double precision, dimension(nx,ny,nz)        :: nu_z       !< LES viscosity derivative on ggg points \f[ \frac{\partial \nu}{\partial z} \f]
!        Cs calculation
          ! DOUBLE PRECISION, DIMENSION(9)             :: tt   
          ! DOUBLE PRECISION, DIMENSION(nx,ny,nz,9)    :: Tensor        
!        matrices
!
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmx      !< x derivative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmy      !< y dervative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: dmz      !< y dervative matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: bx       !< x prolongation matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: by       !< y prolongation matrix
         DOUBLE PRECISION, DIMENSION(:,:), POINTER  :: bz       !< y prolongation matrix
!
!        face arrays
!
         INTEGER, DIMENSION(6)                      :: ibtype
         CHARACTER(LEN = 9) , DIMENSION(6)          :: bcond

         CHARACTER(LEN=8)                           :: domain_type
!
!        FMDF Quantities
!
         double precision, dimension(nx,ny,nz,20)        :: s     !< Scalar arrays for FMDF
         double precision, dimension(nx+2,ny+2,nz+2,20)  :: sumT  !< Summation array for numerator of Favre average
         double precision, dimension(nx+2,ny+2,nz+2)     :: sumW  !< Summation array for denominator of Favre average
         double precision, dimension(nx+2,ny+2,nz+2)     :: sumR  !< Summation array for reaction source term
         double precision, dimension(nx,ny,nz)           :: eRho  !< \f$ \langle \rho \rangle_l \f$
         double precision, dimension(nx,ny,nz)           :: react !< Reaction source on particle
         double precision, dimension(nx,ny,nz)           :: DP		!< pressure compressibility term
         double precision, dimension(nx,ny,nz)           :: Po		!< old pressure array
         double precision, dimension(nx,ny,nz)           :: Poo   !< old pressure array
         integer, dimension(nx,ny,nz)                    :: pCount!< Particles attributed to gauss point
         integer                                         :: np    !< number of particles in element

      END TYPE
!
      END MODULE domain_definition
!
!//////////////////////////////////////////////////////////////////////////////////////
!> All numerical matrices
      MODULE Order_Matrices
         USE size
         SAVE
!
!        ---------------------------------------------
!        element derivative and interpolation matrices
!        ---------------------------------------------
!
         INTEGER                                               :: num_ordersx = 0, num_ordersy = 0
         INTEGER                                               :: num_ordersz = 0
         INTEGER, DIMENSION(max_orders)                        :: order_mapx  = 0, order_mapy  = 0
         INTEGER, DIMENSION(max_orders)                        :: order_mapz  = 0
         DOUBLE PRECISION, DIMENSION(nx,nx,max_orders), TARGET :: bx  = 0.0d0
         DOUBLE PRECISION, DIMENSION(ny,ny,max_orders), TARGET :: by  = 0.0d0
         DOUBLE PRECISION, DIMENSION(nz,nz,max_orders), TARGET :: bz  = 0.0d0
         DOUBLE PRECISION, DIMENSION(nx,nx,max_orders), TARGET :: dmx = 0.0d0
         DOUBLE PRECISION, DIMENSION(ny,ny,max_orders), TARGET :: dmy = 0.0d0
         DOUBLE PRECISION, DIMENSION(nz,nz,max_orders), TARGET :: dmz = 0.0d0
!
!        ---------------------------------------------
!        mortar restriction and prolongation matrices
!        ---------------------------------------------
!
         INTEGER                                                     :: num_mrtr_matrices = 0
         INTEGER, DIMENSION(2,max_mrtr_mat)                          :: mrtr_mat_map = 0
         DOUBLE PRECISION, DIMENSION(nmax,nmax,max_mrtr_mat), TARGET :: prolong = 0.0d0,restrict = 0.0d0
!
      END MODULE Order_Matrices
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Material properties if material is set
      MODULE Material_Properties
         USE size
         SAVE
!
!        ---------------------------------------------
!        save material properties here
!        ---------------------------------------------
!
         INTEGER                                               :: num_prop !< number of properties
         DOUBLE PRECISION, DIMENSION(mxprops,max_materials)    :: material_property !< the material properties
!
      END MODULE Material_Properties
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Main mortar definitions stored here
      MODULE mortar_definition
         USE size
         SAVE
!> Interface mortar datatype for FE elements
      TYPE mortar
!
         INTEGER :: lenmortar(2) !< dimension of the mortar
         INTEGER :: len(2,2)     !< dimension of element face (2 directions x 2 sides)
         INTEGER :: id(2)        !< id of the contributing elements
         INTEGER :: iface(2)     !< face # of contributing elements
         INTEGER :: orient       !< orientation of #2 element face (orient*pi/2)
         INTEGER :: sign         !< orientation of normals (+1/-1)
         INTEGER :: nsign        !< orientation of normals of bordering domain faces  (+1/-1)
!
         DOUBLE PRECISION, DIMENSION(3,nmax,nmax)     :: n_hat                         !< face normal
         DOUBLE PRECISION, DIMENSION(nmax,nmax)       :: norm                          !< normalization
         DOUBLE PRECISION, DIMENSION(nmax,nmax,6)     :: S                             !< normalization
         DOUBLE PRECISION, DIMENSION(:,:), POINTER    :: prolong1_dir1,prolong1_dir2   !< prolongation arrays
         DOUBLE PRECISION, DIMENSION(:,:), POINTER    :: prolong2_dir1,prolong2_dir2   !< prolongation arrays
         DOUBLE PRECISION, DIMENSION(:,:), POINTER    :: restrict1_dir1,restrict1_dir2 !< restriction arrays
         DOUBLE PRECISION, DIMENSION(:,:), POINTER    :: restrict2_dir1,restrict2_dir2 !< restriction arrays
         DOUBLE PRECISION, DIMENSION(nmax,nmax,neq,2) :: Q                  !< solutions along mortar
         DOUBLE PRECISION, DIMENSION(nmax,nmax,neq,2) :: fv                 !< viscous flux along mortar
         DOUBLE PRECISION, DIMENSION(3,nmax,nmax)     :: xg                 !< positions along mortar
         LOGICAL         , DIMENSION(2,2)             :: conforming         !< mortar is conforming or not (2 directions x 2 sides)
!
      END TYPE
!
      END MODULE mortar_definition
!
!///////////////////////////////////////////////////////////////////////
!> MPI particle global storage
      MODULE mpi_par
      USE size

      INTEGER                          :: ierr,myid,comm1d,nreq
      INTEGER                          :: numprocs !< number of processors
      INTEGER,ALLOCATABLE              :: stat(:) !< statistics array for particles
      INTEGER,ALLOCATABLE              :: stata(:,:) !< secondary statistics array for particles
      INTEGER,DIMENSION(nmp*8*2)       :: req            
      INTEGER                          :: ngridl !< number of local elements to core
      INTEGER                          :: nmortl,nedge,nmortd
      INTEGER,ALLOCATABLE              :: ngridmpi(:,:)
      INTEGER,ALLOCATABLE              :: ngridedge(:)
      INTEGER,ALLOCATABLE              :: nmortmpi(:,:)
      INTEGER,ALLOCATABLE              :: nslavempi(:,:)
      INTEGER,ALLOCATABLE              :: nmortedge(:,:)
      INTEGER,ALLOCATABLE              :: nmortdouble(:)
      INTEGER,ALLOCATABLE              :: nn_procs(:)   !/////// added by cmparing by Ahmads Code.

      END MODULE mpi_par
!
!//////////////////////////////////////////////////////////////////////////////////////
!> MPI particle global storage
      MODULE mpi_par_part
      USE size
          integer, parameter  :: send_p_length = 33 !< length of particle array for communication
          INTEGER             :: npart_send_tot,npart_send,npart_recv
          INTEGER, PARAMETER  :: npart_mpi = npart/4!npart/10+10
          INTEGER             :: local_part_mpi(npart_mpi*3)
          INTEGER,ALLOCATABLE :: global_part_mpi(:,:)
          DOUBLE PRECISION    :: send_drop(npart_mpi*send_p_length)
          DOUBLE PRECISION    :: recv_drop(npart_mpi*send_p_length)
      END MODULE mpi_par_part
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!> Contains statistics data type
      MODULE stats_definition
         USE size
         SAVE
!
         INTEGER                                        :: nsample      !< number of samples
!> Average statistics data type
      TYPE statfl
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq+nav)     :: Q_av         !< Average properties
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,4)           :: Q_xyz        !< Average properties
!
      END TYPE
!> Fluctating statistics data type
      TYPE statfluct
         DOUBLE PRECISION, DIMENSION(nx,ny,nz,nfluct)      :: Q_fluct      !< Fluctuating properties
!
      END TYPE

!
  END MODULE  stats_definition
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!> Particle definition
      MODULE particle_definition
!
      USE size
      SAVE
!> Particle datatype with all particle quantities
      TYPE particle
!
         INTEGER                           :: onoff  !< take particle into account (onoff=1)
                                                     ! don't take particle into account (onoff=0)
         INTEGER                           :: ngrid  !< subdomain particle is in
!
!     Particle Position, Velocity and Properties
!

         DOUBLE PRECISION, DIMENSION(3)               :: Xp     !< particle coordinate
         DOUBLE PRECISION, DIMENSION(3)               :: Xpmap  !< particle coordinate in mapped space
         DOUBLE PRECISION, DIMENSION(3)               :: Vp     !< particle velocity
         DOUBLE PRECISION                             :: Tp     !< particle temperature
         double precision                             :: T      !< temperature on particle
         DOUBLE PRECISION                             :: Mp     !< particle mass

!
!
!     Particle Position, Velocity and Properties for the Adam-Bashfort time-integration
!

         DOUBLE PRECISION, DIMENSION(3)    :: Xpnm     !< particle coordinate
         DOUBLE PRECISION, DIMENSION(3)    :: Vp1nm    !< particle velocity
         DOUBLE PRECISION                  :: Tpnm     !< particle temperature
         DOUBLE PRECISION                  :: Mpnm     !< particle mass

         DOUBLE PRECISION, DIMENSION(3)    :: Vfp    !< fluid velocity at particle position
         DOUBLE PRECISION                  :: Tfp    !< fluid temperature at particle position
         DOUBLE PRECISION                  :: Rhofp  !< fluid density at particle position
         DOUBLE PRECISION                  :: Yffp   !< fluid vapor fraction at particle position

!       Particle fmdf quantities
        double precision                    :: w         !< particle weight
        double precision                    :: rho       !< particle density
        double precision,dimension(20)      :: scalar    !< scalar particle properties
        double precision                    :: nu_t      !< turbulent viscosity
        double precision,dimension(3)       :: nu_d      !< viscous fluxes
        double precision                    :: omega_m   !< sub-grid mixing frequency
        double precision                    :: P         !< particle pressure
        double precision                    :: rho_a     !< particle analytical density
        double precision                    :: dp        !< material derivative pressure term
        double precision                    :: react     !< particle reaction source term
        double precision                    :: eRho      !< ensemble number density
        double precision,dimension(3,2)     :: b         !< particle basis function
        integer                             :: il,jl,kl  !< particle bounds
        double precision,dimension(3,nx)    :: h         !< spectral basis


!
      END TYPE
!
!       FMDF Pointer particle data type
!> Particle pointer (linked list) datatype
        type ePointers
            type(particle), pointer             :: p !< accessor to particle pointer
        end type
  END MODULE particle_definition
!
!///////////////////////////////////////////////////////////////////////
!> Parallel quantities for particles
      MODULE part_par

      INTEGER            :: IEVAP
      DOUBLE PRECISION   :: SCF,SCO,DA,ZE,CE,RCOEF,YFINIT
      DOUBLE PRECISION   :: PEF,PEO,PD0M,PM0,PSI,CRE,CTAUP,CYFPS,CF2,CF3,CF4,PD0
      INTEGER            :: nprtmax
      INTEGER            :: nprtrepl = 0

      END MODULE part_par
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Finite element data type module
      MODULE FE_Data_Types

      USE size
!> Node data type
      TYPE node
         DOUBLE PRECISION, DIMENSION(3) :: x
         LOGICAL                        :: on_boundary
      END TYPE node
      !> Node vector data type
      TYPE node_vector
         INTEGER                          :: num_nodes
         TYPE(node), DIMENSION(max_nodes) :: node
      END TYPE node_vector
!> Element face data type
      TYPE face
         INTEGER          :: key
         INTEGER          :: node(16)
         INTEGER          :: from_element(2)
         INTEGER          :: from_side(2)
         DOUBLE PRECISION :: start(2),end(2)
      END TYPE face
      !> Face pointer list data type
      TYPE face_list
         TYPE (face)     , POINTER  :: face_data
         TYPE (face_list), POINTER  :: previous,next
      END TYPE
      !> Face pointer
      TYPE face_pointer
         TYPE (face), POINTER :: the_face
      END TYPE
      !> Face interpolation data type
      TYPE face_interp ! data for interpolating a face
         DOUBLE PRECISION, DIMENSION(3,nmax,nmax) :: points
         DOUBLE PRECISION, DIMENSION(nmax,2)      :: knots
         INTEGER         , DIMENSION(2)           :: num_knots
      END TYPE
      !> Surface data type
      TYPE surface
         INTEGER          :: surface_id
         CHARACTER(132)   :: equation(3)   ! equation for each variable
         CHARACTER(LEN=9) :: surface_type  !  = 'equation' or  = 'element'
         CHARACTER(LEN=9) :: bcond         ! one of the approved bc's
       END TYPE surface

       END MODULE FE_Data_Types
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Zeros the computational domain before reading in conditions
      SUBROUTINE zero_domain(d) 
!
!     date: 11/9/98                                                     
!
!     routines called: none
!
      USE domain_definition
!
      TYPE (domain) :: d
!
      d%type = 8            ! (8 = hex8,... number of nodes in hex)
      d%node = 0
      d%orientation = 0
      d%mortar = 0         ! mortar # and side for each face
      d%which_surface = 0  ! surface a face is attached to

      d%nc  = 0
      d%ncg = 0            ! grid size

      d%Q    = 0.0d0       ! dependent variables
      d%Qlgg = 0.0d0       ! dependent variables on faces
      d%Qglg = 0.0d0       ! dependent variables on faces
      d%Qggl = 0.0d0       ! dependent variables on faces
      d%f    = 0.0d0
      d%g    = 0.0d0
      d%h    = 0.0d0          ! fluxes
      d%g_Q  = 0.0d0        ! rk arrays
!
      d%material_id = 0
!
      d%cx = 0.0d0
      d%xg = 0.0d0
      d%jacob = 0.0d0
      d%gmet = 0.0d0
      ! d%Tensor = 0
      ! d%tt = 0
!
      d%corner = 0.0d0
      d%cxg     = 0.0d0
!
!     matrices
!
      NULLIFY(d%dmx)
      NULLIFY(d%dmy)
      NULLIFY(d%dmz)
      NULLIFY(d%bx)
      NULLIFY(d%by)
      NULLIFY(d%bz)
!
      d%ibtype = 1
      d%bcond = 'interface'
      d%domain_type = 'interior'
!
      RETURN
      END SUBROUTINE zero_domain
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE zero_mortar(mrtr) 
!
!     date: 11/9/98                                                     
!
!     routines called: none
!     applicability: all
!
      USE mortar_definition
!
      TYPE (mortar) :: mrtr
!
      mrtr%lenmortar = 0    ! dimension of the mortar
      mrtr%len = 0          ! dimension of element face (2 directions x 2 sides)
      mrtr%id = 0           ! id of the contributing elements
      mrtr%iface  = 0       ! face # of contributing elements
      mrtr%orient = 0       ! orientation of #2 element face (orient*pi/2)
      mrtr%sign = 0         ! orientation of normals (+1/-1)
      mrtr%nsign = 0         ! orientation of normals (+1/-1)
!
      mrtr%n_hat = 0.0
      mrtr%norm  = 0.0       ! geometry arrays along mortar

      NULLIFY(mrtr%prolong1_dir1,mrtr%prolong1_dir2)
      NULLIFY(mrtr%prolong2_dir1,mrtr%prolong2_dir2)
      NULLIFY(mrtr%restrict1_dir1,mrtr%restrict1_dir2)
      NULLIFY(mrtr%restrict2_dir1,mrtr%restrict2_dir2)
      mrtr%Q  = 0.0                      ! solutions along mortar
      mrtr%fv  = 0.0                      ! solutions along mortar
      mrtr%xg  = 0.0                     ! positions along mortar
      mrtr%conforming = .TRUE.           ! mortar is conforming or not (2 directions x 2 sides)
!
      RETURN
!
      END SUBROUTINE zero_mortar
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Zeros all matrices
      SUBROUTINE zero_matrices
!
         USE Order_Matrices

         num_ordersx = 0
     num_ordersy = 0
         num_ordersz = 0
         order_mapx  = 0
     order_mapy  = 0
         order_mapz  = 0
         bx  = 0.0d0
         by  = 0.0d0
         bz  = 0.0d0
         dmx = 0.0d0
         dmy = 0.0d0
         dmz = 0.0d0
         num_mrtr_matrices = 0
         mrtr_mat_map = 0
         prolong = 0.0d0
     restrict = 0.0d0
!
        RETURN
      END SUBROUTINE zero_matrices
!
!///////////////////////////////////////////////////////////////////////
!> Zeros all statistics arrays
      SUBROUTINE zero_stats(stats,statsfluct)
!
!     date: 6/19/01
!
!     routines called: none
!     applicability: all
!
      USE stats_definition
!
      TYPE (statfl)    :: stats
      TYPE (statfluct) :: statsfluct

      stats%Q_av              = 0.0d0
      stats%Q_xyz             = 0.0d0
      statsfluct%Q_fluct      = 0.0d0

      RETURN
!
   END SUBROUTINE zero_stats
!
!///////////////////////////////////////////////////////////////////////
!> Zeros all quantities on particle for entire list
      SUBROUTINE zero_drops(drop)
!
      USE particle_definition
!
      TYPE(particle) :: drop
!
      drop%onoff = 0
      drop%ngrid = 0

!
      drop%Xp    = 0.0d0
      drop%Xpmap = 0.0d0
      drop%Vp    = 0.0d0
      drop%Mp    = 0.0d0
!
      drop%Xpnm  = 0.0d0
      drop%Vp1nm = 0.0d0
      drop%Tpnm  = 0.0d0
      drop%Mpnm  = 0.0d0
!
      drop%Vfp   = 0.0d0
      drop%Tfp   = 0.0d0
      drop%Rhofp = 0.0d0
      drop%Yffp  = 0.0d0
!
      RETURN
!
   END SUBROUTINE zero_drops
!
!///////////////////////////////////////////////////////////////////////
!
