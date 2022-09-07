!> \file input_lib.f90
!! For reading bump file and setting user properties
!///////////////////////////////////////////////////////////////
!////////						////////
!////////	input_lib.f90				////////
!////////						////////
!////////	contains:				////////
!////////						////////
!////////	   MODULE File_Units			////////
!////////	   MODULE keywords			////////
!////////	   MODULE input_data			////////
!////////	   SUBROUTINE get_input_values(iunit)	////////
!////////						////////
!///////////////////////////////////////////////////////////////
!
!
!.......................................................................
!>     Store the file units to be used for i/o
!!
!!     i/o files and units used by main code (units 1-19 are reserved)
!!     These can be changed here according to operating system.
!!
!!     input file           = unit 5
!!     standard output      = unit 6
!!     error output         = unit 6
!!     mesh file            = unit 7
!!     geometry file        = unit 8
!!     restart file         = unit 11
!!     plot file            = unit 12      
!!     fv dat file          = unit 14
!!     fvs dat file         = unit 17
!!     fvp dat file         = unit 18
!!     droplet file         = unit 19
   MODULE File_Units

      INTEGER, PARAMETER  :: mesh_unit = 7,spec_unit = 8
      INTEGER, PARAMETER  :: in_unit   = 5,rst_unit = 11
      INTEGER, PARAMETER  :: iout_unit = 15,plt_unit = 12
      INTEGER, PARAMETER  :: drop_unit = 19
      INTEGER, PARAMETER  :: fvt_unit = 14
      INTEGER, PARAMETER  :: fvs_unit = 17
      INTEGER, PARAMETER  :: stat_unit = 7
      INTEGER, PARAMETER  :: fvp_unit = 18
      INTEGER, PARAMETER  :: err_unit  = 15
      INTEGER, PARAMETER  :: dissipated_unit = 20
      ! INTEGER, PARAMETER  :: average_unit = 21
      INTEGER, PARAMETER  :: rms_unit = 22
      INTEGER, PARAMETER  :: espec_unit = 23
      INTEGER, PARAMETER  :: vcalc_unit = 24
      INTEGER, PARAMETER  :: smagvcalc_unit = 25
      INTEGER, PARAMETER  :: umodal_unit = 30
      INTEGER, PARAMETER  :: vmodal_unit = 31
      INTEGER, PARAMETER  :: wmodal_unit = 32

   END MODULE File_Units
!
!///////////////////////////////////////////////////////////////////////
!> Store the keywords needed to read the input file
   MODULE keywords
!
!.......................................................................
!     store the keywords needed to read the input file
!.......................................................................
!
      SAVE
!
!     ---------------------------
!     keywords for setting values
!     ---------------------------
!
      INTEGER, PARAMETER                :: max_wrd = 18
      INTEGER, PARAMETER                :: max_tgl = 28
      INTEGER                           :: num_keywords = max_wrd
      CHARACTER (LEN=10), DIMENSION(max_wrd) :: keyword_list =  &
                                               (/ 'title     ', &
                                                  'run       ', &
                                                  'mesh      ', &
                                                  'geometry  ', &
                                                  'resolution', &
                                                  'highest   ', &
                                                  'lowest    ', &
                                                  'length    ', &
                                                  'order     ', &
                                                  'source    ', &
                                                  'cfl       ', &
                                                  'final     ', &
                                                  'maximum   ', &
                                                  'frequency ', &
                                                  'interval  ', &
                                                  'droplet   ', &
                                                  'drinter   ', &
                                                  'drfreq    '  &
                                               /)
      LOGICAL, DIMENSION(max_wrd) :: keyword_set = .false.
!
!     --------------------------------------------------------
!     Toggle words take on either one value or the next in the
!     list. for example or logal/global
!     --------------------------------------------------------
!
      INTEGER                           ::num_togglewords = max_tgl
      CHARACTER (len=10), DIMENSION(max_tgl) :: toggle_list =   &
                                               (/               &
                                                  'movie      ', &
                                                  'endplot    ', &
                                                  'local      ', &
                                                  'global     ', &
                                                  'end        ', &
                                                  'periodic   ', &
                                                  'accurate   ', &
                                                  'steady     ', &
                                                  'verbose    ', &
                                                  'compact    ', &
                                                  'statistic  ', &
                                                  'average    ', &
                                                  'rms        ', &
                                                  'drops      ', &
                                                  'nodrops    ', &
                                                  'restfldr   ', &
                                                  'drendpl    ', &
                                                  'drmov      ', &
                                                  'metis      ', &
                                                  'continue   ', &
                                                  'smagorinsk ', &
                                                  'rhoSens    ', &
                                                  'shockSens  ', &
                                                  'shock      ', &
                                                  'ducrosSens ', &
                                                  'smoothEV   ', &
                                                  'fmdf       ', &
                                                  'prntdis    '  &
                                         /)
      LOGICAL, DIMENSION(max_tgl) :: toggle_set = .false.
!
   END MODULE keywords
!
!///////////////////////////////////////////////////////////////////////
!> Storage for data read in from standard input
   MODULE input_data
!
!.......................................................................
!     storage for data read in from standard input
!.......................................................................
!
      SAVE
!
      LOGICAL              :: local = .false., global = .true.       ! local/global time step
      LOGICAL              :: verbose  = .false., compact = .true.   ! do verbose or compact output for element information
      LOGICAL              :: movie    = .false., endplot = .true.   ! write or don't write movie file
      LOGICAL              :: average    = .true., rms = .false.     ! average statistics, rms statistics
      LOGICAL              :: rmscont = .false.                      ! continue averaging
      LOGICAL              :: statistics = .false.                   ! if true, compute statistics

      LOGICAL              :: time_accurate = .false., steady_state = .true.
      LOGICAL              :: defined_resolution = .false.           ! if true, compute the approx order
      LOGICAL              :: save_restarts = .false.                ! if true, save restart files periodically
      LOGICAL              :: drmov    = .false., drendpl = .true.   ! write or don't write movie file
      LOGICAL              :: restfldrop = .false.                   ! if true, restart fluid + insert particles
      LOGICAL              :: drops    = .false., nodrops = .true.   ! solve or don't solve for droplets
      LOGICAL              :: dropper = .false.                      ! if true, restart fluid + insert particles
      LOGICAL              :: usemetis = .false.                     ! if true, create/read metis file
      LOGICAL              :: smagorinsky = .false.                  ! if true, use smagorinsky turbulence model
      LOGICAL              :: rhoSensor = .false.                    ! if true, apply rho sensor to smagorinsky turbulence model
      LOGICAL              :: shockSensor = .false.                  ! if true, apply shock sensor to smagorinsky turbulence model
      LOGICAL              :: shock = .false.                        ! if true, use shock capturing method
      LOGICAL              :: ducrosSensor = .false.                 ! if true, apply Ducros sensor to EV method
      LOGICAL              :: Smoothing = .false.                    ! if true, apply EV smoothing
      logical              :: isfmdf   = .false.                     ! Do we turn on FMDF?
      logical              :: dissipation   = .false.                ! used to compute energy dissipation
!
      INTEGER              :: max_steps = 0                          ! maximum number of steps
      INTEGER              :: nout = 1                               ! interval between calls to do_user_stuff
      INTEGER              :: npout = 1                               ! interval between calls to insert particles
      INTEGER              :: iord = 3                               ! runge-kutta order
      INTEGER              :: min_order = 5                          !
      INTEGER              :: max_order                              !
!
      DOUBLE PRECISION     :: stab = 0.9d0                           ! cfl number
      DOUBLE PRECISION     :: tfinal = 0.0d0, tout = 1.0d0           ! final time, time interval for output
      DOUBLE PRECISION     :: tpout = 1.0d0                          ! time interval for particle insert
      DOUBLE PRECISION     :: char_length = 1.d0                     ! characteristic length scale
      DOUBLE PRECISION     :: resolution                             ! number of points/characteristic length
      CHARACTER (LEN = 80) :: title = 'parabolic computation'
      CHARACTER (LEN = 32) :: filename       = 'none'                ! restart/output file names
      CHARACTER (LEN = 32) :: mesh_filename  = 'none'                ! mesh or restart file name
      CHARACTER (LEN = 32) :: geom_filename  = 'none'                ! geometry file name
      CHARACTER (LEN = 32) :: drop_filename  = 'none'                ! droplet file name

   END MODULE input_data
!
!///////////////////////////////////////////////////////////////////////
!> Read in the input data-file from unit "iunit" and save the values
   SUBROUTINE get_input_values(iunit)
!
!.......................................................................
!     Read in the input data-file from unit "iunit"  and save the values
!     these values are saved in the modules used below
!.......................................................................
!
      USE size
      USE keywords
      USE input_data
      USE physics
      USE user_keywords
      USE physics_keywords
      USE constants
      USE File_Units
!
      IMPLICIT none
!
      DOUBLE PRECISION      :: get_dp_value
      INTEGER               :: get_int_value
      CHARACTER (LEN = 132) :: input_line, get_string_value
      CHARACTER (LEN = 10)  :: keyword_in_, the_word
      CHARACTER (LEN = 7)   :: the_list
      LOGICAL               :: eof
      INTEGER               :: iunit,k,kword,length,n
!
!     -----------------------------------------------------
!     write out non-comment input lines as they are read in
!     -----------------------------------------------------
!
!        WRITE(iout_unit,*) '********************************************'
!        WRITE(iout_unit,*) '*                                          *'
!        WRITE(iout_unit,*) '*             input states                 *'
!        WRITE(iout_unit,*) '*                                          *'
!        WRITE(iout_unit,*) '********************************************'
!        WRITE(iout_unit,*) ' '
!
      DO
         CALL skip_comments(input_line,iunit,eof)
         IF(eof)     EXIT
         the_word = keyword_in_(input_line,keyword_list,num_keywords,kword)
         IF(the_word == 'no keyword')     THEN ! try the toggle list
            the_word = keyword_in_(input_line,toggle_list,num_togglewords,kword)
         END IF
         IF(the_word == 'no keyword')     THEN ! try the physics list
            the_word = keyword_in_(input_line,physics_keyword_list,num_physics_keywords,kword)
            IF(the_word /= 'no keyword') the_list = 'physics'
         END IF
         IF(the_word == 'no keyword')     THEN ! try the user list
            the_word = keyword_in_(input_line,user_keyword_list,num_user_keywords,kword)
            IF(the_word /= 'no keyword') the_list = 'user'
         END IF
!        WRITE(iout_unit,*) trim(input_line)
!
!        ---------------------------------------------
!        first see if it is a general program variable
!        ---------------------------------------------
!
         SELECT CASE(the_word)

            CASE('no keyword')
!              WRITE(iout_unit,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
!              WRITE(iout_unit,*) 'no keyword found in this input line.'
!              WRITE(iout_unit,*) 'check spellings, or put a comment "!" in column one'
!              WRITE(iout_unit,*) 'input line ignored'

            CASE('title')
               title = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('source')
               source_type = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('run')
               filename = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('mesh')
               mesh_filename = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('geometry')
               geom_filename = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('order')
               iord = get_int_value(input_line)
               keyword_set(kword) = .true.

            CASE('maximum')
               max_steps = get_int_value(input_line)
               keyword_set(kword) = .true.

            CASE('highest')
               max_order = get_int_value(input_line)
               keyword_set(kword) = .true.
               defined_resolution = .true.
               max_order          = min(max_order,nx,ny)

            CASE('lowest')
               min_order = get_int_value(input_line)
               keyword_set(kword) = .true.
               defined_resolution = .true.

            CASE('cfl')
               stab = get_dp_value(input_line)
               keyword_set(kword) = .true.

            CASE('length')
               char_length = get_dp_value(input_line)
               keyword_set(kword) = .true.

            CASE('resolution')
               resolution = get_dp_value(input_line)
               keyword_set(kword) = .true.
               defined_resolution = .true.

            CASE('final')
               tfinal = get_dp_value(input_line)
               keyword_set(kword) = .true.

            CASE('interval')
               tout = get_dp_value(input_line)
               keyword_set(kword) = .true.

            CASE('frequency')
               nout = get_int_value(input_line)
               keyword_set(kword) = .true.

            CASE('drinter')
               tpout = get_dp_value(input_line)
               dropper = .true.
               keyword_set(kword) = .true.

            CASE('drfreq')
               npout = get_int_value(input_line)
               dropper = .true.
               keyword_set(kword) = .true.

            CASE('local')
               local  = .true.
               global = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('global')
               global = .true.
               local  = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('end')
               save_restarts = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('periodic')
               save_restarts  = .true.
               toggle_set((kword+1)/2) = .true.

            CASE('accurate')
               time_accurate = .true.
               steady_state  = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('steady')
               steady_state  = .true.
               time_accurate = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('movie')
               movie = .true.
               endplot = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('endplot')
               endplot = .true.
               movie = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('verbose')
               verbose = .true.
               compact = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('statistic')
               statistics  = .true.
               toggle_set((kword+1)/2) = .true.

            CASE('average')
               average  = .true.
               rms      = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('rms')
               average  = .false.
               rms      = .true.
               toggle_set((kword+1)/2) = .true.

            CASE('continue')
               rmscont  = .true.
               toggle_set((kword+1)/2) = .true.

            CASE('droplet')
               drop_filename = get_string_value(input_line)
               keyword_set(kword) = .true.

            CASE('drops')
               drops   = .true.
               nodrops = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('nodrops')
               nodrops = .true.
               drops   = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('drmov')
               drmov = .true.
               drendpl = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('drendpl')
               drendpl = .true.
               drmov = .false.
               toggle_set((kword+1)/2) = .true.

            CASE('restfldr')
               restfldrop  = .true.
               toggle_set((kword+1)/2) = .true.

            CASE('metis')
               usemetis  = .true.
               toggle_set((kword+1)/2) = .true.
               
            CASE('smagorinsk')
               smagorinsky  = .true.
               toggle_set((kword+1)/2) = .true.
               
            CASE('rhoSens')
               rhoSensor  = .true.
               toggle_set((kword+1)/2) = .true.
               
            CASE('shockSens')
               shockSensor  = .true.
               toggle_set((kword+1)/2) = .true.
               
            CASE('shock')
               shock  = .true.
               toggle_set((kword+1)/2) = .true.
             
            CASE('ducrosSens')
               ducrosSensor  = .true.
               toggle_set((kword+1)/2) = .true.
             
            CASE('smoothEV')
               Smoothing  = .true.
               toggle_set((kword+1)/2) = .true.
             
            case('fmdf')
                isfmdf   = .true.
                toggle_set((kword+1)/2) = .true.
            
            case('prntdis')
               dissipation = .true.
               toggle_set((kword+1)/2) =.true.


            CASE default ! it must be one of the user keywords or physics keywords

               IF(the_list == 'physics')     THEN
                   CALL set_physics_values(the_word,input_line,kword)
               ELSE
                   CALL set_user_values(the_word,input_line,kword)
               END IF

         END SELECT
!
      END DO
!
!     -------------------
!     set mode variables:
!     -------------------
!
!     --------------------------------------------------------------------------
!     If the resolution is not defined by any variables, then it is assumed that
!     the number of grid points is determined by the input file
!     --------------------------------------------------------------------------
!
      IF( .NOT.defined_resolution )     THEN
         DO k = 1,num_keywords
            IF (keyword_list(k) == 'lowest    ')       keyword_set(k) = .true.
            IF (keyword_list(k) == 'length    ')       keyword_set(k) = .true.
            IF (keyword_list(k) == 'highest   ')       keyword_set(k) = .true.
            IF (keyword_list(k) == 'resolution')       keyword_set(k) = .true.
         END DO
      END IF
!
!     --------------------------------------------------------------------
!     if computation is time accurate then need  final_time, time interval
!     otherwise, if steady, then need output interval and local/global
!     --------------------------------------------------------------------
!
      IF(time_accurate)     THEN
         DO k = 1,num_keywords
            IF (keyword_list(k) == 'frequency ')       keyword_set(k) = .true.
         END DO
         DO k = 1,num_togglewords
            IF (toggle_list(k)  == 'local     ')       toggle_set((k+1)/2) = .true.
         END DO
      ELSE                          ! steady-state
         DO k = 1,num_keywords
            IF (keyword_list(k) == 'final     ')       keyword_set(k) = .true.
            IF (keyword_list(k) == 'interval  ')       keyword_set(k) = .true.
         END DO
      END IF
!
!     -----------------------------------------------------
!     alert user to which values have only the defaults set
!     -----------------------------------------------------
!
      DO k = 1,num_keywords
         IF(.not.keyword_set(k))     THEN
!          WRITE(iout_unit,*) 'default value for keyword ',TRIM(keyword_list(k)),' used'
         END IF
      END DO
      DO k = 1,num_togglewords/2
         IF(.not.toggle_set(k))     THEN
!          WRITE(iout_unit,*) 'default value for keywords ',TRIM(toggle_list(2*k-1)),&
!                      ' and ',TRIM(toggle_list(2*k)),' used'
         END IF
      END DO
      DO k = 1,num_physics_keywords
         IF(.not.physics_keyword_set(k))     THEN
!          write(iout_unit,*) 'default value for keyword ',TRIM(physics_keyword_list(k)),' used'
         END IF
      END DO
      DO k = 1,num_user_keywords
         IF(.not.user_keyword_set(k))     THEN
!          write(iout_unit,*) 'default value for keyword ',TRIM(user_keyword_list(k)),' used'
         END IF
      END DO
!
      RETURN
   END SUBROUTINE get_input_values
