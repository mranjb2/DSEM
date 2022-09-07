!
       PROGRAM Driver_3D
!
!.......................................................................
!      solve conservative form equations
!
!      creation date: 8/27/98                                            
!      modifications:                                                    
!
!      user manual: see www.scri.fsu.edu/~dak/dak-home.html/research/stac^m
!
!.......................................................................
!
!     -----------------
!     modules in use...
!
!     dimensions:
!     -----------------
!
      USE size
!
!     -----------------
!     Definitions:
!     -----------------
!
      USE domain_definition
      USE mortar_definition
      USE stats_definition
      USE particle_definition
      USE keywords
      USE FE_Data_Types
      USE Method_Data
      USE File_Units
!
!     -----------------
!     Data storage:
!     -----------------
!
      USE constants
      USE physics
      USE rk_coefs
      USE input_data
      USE Order_Matrices
      USE Material_Properties
      USE User_Data
      USE part_par
!
!     -----------------
!     Work modules
!     -----------------
!
      USE mpi_par
      USE mpi_par_part
!
!     -----------------
!     Other definitions
!     -----------------
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (node_vector)                  :: nodes
      TYPE (domain), DIMENSION(ngp)       :: d             ! define array of subdomains 
      TYPE (mortar), DIMENSION(nmp)       :: mrtr          ! define array of mortars
      TYPE (statfl)  , DIMENSION(ngp)     :: stats ! define array of particles
      TYPE (statfluct)  , DIMENSION(ngp)  :: statsfluct ! define array of particles
      TYPE (particle), DIMENSION(npart)   :: drop          ! define array of particles

      TYPE (surface), DIMENSION(max_surf) :: surface_array ! boundary surfaces
!
      DOUBLE PRECISION, DIMENSION(ngp) :: delta_t
      DOUBLE PRECISION                 :: dt
      CHARACTER(LEN=32)                :: rstname,outname,pltname,fvtname,fvsname,fname,fnamea,staname,fvpname
      CHARACTER                        :: date*8,time_str*10,zone*5
      LOGICAL                          :: restart,is_NaN,Is_INF
!
!     ----------------------------------
!     variables for elapsed time logging
!     ----------------------------------
!
      DOUBLE PRECISION                 :: compute_time
      INTEGER,DIMENSION(8)             :: start_time_values,end_time_values,elapsed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MPI INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INCLUDE 'mpif.h'

      CALL MPI_INIT( ierr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      WRITE(*,*) "proc",myid, "of",numprocs,"is alive"

      comm1d   = MPI_COMM_WORLD
      ALLOCATE(stat(MPI_STATUS_SIZE),STAT=ierr)   
      ALLOCATE(stata(MPI_STATUS_SIZE,nmp*8*2),STAT=ierr)             
      ALLOCATE(ngridmpi(ngp*numprocs,2),STAT=ierr)             
      ALLOCATE(ngridedge(ngp*numprocs),STAT=ierr)             
      ALLOCATE(nmortmpi(nmp*numprocs,3),STAT=ierr)             
      ALLOCATE(nslavempi(nmp*numprocs,3),STAT=ierr)             
      ALLOCATE(nmortedge(nmp*numprocs,3),STAT=ierr)             
      ALLOCATE(nmortdouble(nmp*numprocs),STAT=ierr)             
      ALLOCATE(global_part_mpi(npart_mpi*numprocs,3),STAT=ierr)             
      ngridmpi        = 0
      ngridedge       = 0
      nmortmpi        = 0
      nslavempi       = 0
      nmortedge       = 0
      nmortdouble     = 0
      global_part_mpi = 0

!
!     -------------------------
!     compute machine constants
!     -------------------------
!
      eps = EPSILON(1.d0)
      big = HUGE(1.d0)
!
!     -----------------
!     zero out objects                                                             
!     -----------------
!
      k = 0
      ngrid = 0
      nprt  = 0
      time  = 0.0d0
      DO id = 1,ngp
         CALL zero_domain(d(id))
         CALL zero_stats(stats(id),statsfluct(id))
      END DO
      DO j = 1,nmp
         CALL zero_mortar(mrtr(j))
      END DO
      DO ip=1,npart
         CALL zero_drops(drop(ip))
      ENDDO
!
!     -----------------------------------
!     read in data from initial data file
!     -----------------------------------
!
      OPEN (unit=in_unit,file = "bump.in", status='old')  ! replace for unix machines
         CALL get_input_values(in_unit)
         REWIND(in_unit)
         IF (myid == 0 ) THEN
           CALL set_io_files(rstname,outname,pltname,fvtname,fvsname,restart,staname,fvpname)
         ENDIF
         CALL MPI_BCAST(restart,1,MPI_LOGICAL,0,comm1d,ierr)
         CALL MPI_BCAST(outname,32,MPI_CHARACTER,0,comm1d,ierr)
         OPEN (unit=iout_unit,file=outname)
         !IF (myid == 0) CALL write_file(in_unit,iout_unit)
      CLOSE (in_unit)

!
!     ------
!     setup:                             
!     ------
!
      CALL set_up_physics      ! Equation dependent physics parameters
      CALL rkcoefs(iord)       ! Set up runge-kutta coefficients
!
      IF( restart )     THEN   ! read the restart file

        IF (myid==0) THEN
          OPEN (unit=rst_unit,file=filename,status='old',form='unformatted')
        ENDIF

            CALL read_restart(rst_unit,ngrid,d,nmort,mrtr,stats,statsfluct,drop,nprt,time)
       
        IF (myid==0) THEN
          CLOSE(rst_unit)
        ENDIF
        
        IF (defined_resolution) THEN
            CALL set_no_points(d,mrtr,nmort) 
            CALL zero_matrices
            DO id = 1,ngridl
              CALL set_points(d(id))
              CALL set_geometry_refine(d(id),mrtr)
              CALL set_matrices(d(id))
            ENDDO
            
            CALL reset_pointers(d,mrtr)
!
            time = 0.0 
            DO id = 1,ngridl 
               CALL Xinitc(id,d(id),mrtr)
            END DO
!
            DO id = 1,ngridl
               DO nv = 1,neq
                  DO l = 1,d(id)%ncg(3)
                     DO m = 1,d(id)%ncg(2) 
                        DO n = 1,d(id)%ncg(1)
                           d(id)%Q(n,m,l,nv)  = d(id)%Q(n,m,l,nv)/d(id)%jacob(n,m,l)
                        END DO
                     END DO 
                  END DO 
               END DO
            END DO
!
 
        ENDIF
!
        IF (drops) THEN
          IF (restfldrop) THEN
		   write(*,*)'Reading droplet file...'
           OPEN (unit=drop_unit,file=drop_filename,status='old')
           CALL readpart(drop_unit,nprt,drop,ngridl,d)
		   write(*,*)'There are ',nprt,'droplets added to the flow'
		  else
			  write(*,*)'Droplets read from restart file...'
          ENDIF
!
!     ------------------------------------------------------------------
!     initialize particle parameters and find the subdomain the particle
!     is in and map the particle location and velocity to that domain.
!     ------------------------------------------------------------------

           CALL set_part_par()
           DO ip=1,nprt
             IF (drop(ip)%onoff == 1) THEN
               ip_grid = drop(ip)%ngrid
               CALL part_map(d(ip_grid),drop(ip))
               IF (d(ip_grid)%ncg(1) <7) THEN
                 CALL spec_Pinterp(d(ip_grid),drop(ip),myid,ip)
               ELSE
                   CALL Lagr_Pinterp(d(ip_grid),drop(ip),myid,ip)
               ENDIF
               IF (restfldrop) THEN
                 CALL part_initc(drop(ip),nprt)
               ENDIF
             ENDIF
           END DO
           IF (restfldrop) THEN
             CLOSE(drop_unit)
           ENDIF
        END IF
!
        CALL user_setup(d,ngridl,mrtr,nmortl,restart) ! do user defined operations
!
      ELSE                     ! start from scratch
!
!        ------------------------------
!        read in mesh and geometry data
!        ------------------------------
!
        IF (myid==0) THEN
         OPEN (unit=mesh_unit,file=mesh_filename,status='old') 
         OPEN (unit=spec_unit,file=geom_filename,status='old')
            CALL read_mesh_file(mesh_unit,d,ngrid,nodes)
            CALL read_specification_file(spec_unit,num_surf,surface_array)
         CLOSE(mesh_unit)
         CLOSE(spec_unit)

!
!        ----------------
!        generate mortars
!        ----------------
!
         CALL generate_mortars(d,ngrid,mrtr,nmort)
!
!        -----------------------------------------------
!        Set up spectral related arrays, then 
!        compute domain and mortar grid and metric terms
!        -----------------------------------------------
!
         DO id = 1,ngrid            
            CALL set_points(d(id))   ! gauss points/weights, etc.
            CALL set_geometry(nodes,d(id),mrtr,num_surf,surface_array)
         END DO


         CALL user_setup(d,ngrid,mrtr,nmort,restart)

        ENDIF

         CALL grid_partitioning(d,mrtr,ngrid,nmort)

         DO id = 1,ngridl            
            CALL set_matrices(d(id)) ! differentiation and interpolation matrices
         END DO

         DO j = 1,nmortl
            CALL set_mortar_matrices(mrtr(j))
            CALL set_mortar_nsign(mrtr(j),d) 
         END DO
!
!
!        --------------------------------------
!        set up initial conditions on each grid
!        --------------------------------------
!
         time = 0.0 
         DO id = 1,ngridl 
            CALL initc(id,d(id),mrtr)
         END DO
!
!        ----------------------------------------------------------
!        convert initial conditions, Q, to JQ Note that jacob = 1/J
!        ----------------------------------------------------------
!
         DO id = 1,ngridl
            DO nv = 1,neq
               DO l = 1,d(id)%ncg(3)
                  DO m = 1,d(id)%ncg(2) 
                     DO n = 1,d(id)%ncg(1)
                        d(id)%Q(n,m,l,nv)  = d(id)%Q(n,m,l,nv)/d(id)%jacob(n,m,l)
                     END DO
                  END DO 
               END DO 
            END DO
         END DO
!
!        ----------------------------------------------------------
!        Update solutions from restart file (if needed)
!        ----------------------------------------------------------
! 
!          CALL read_soution(ngrid,d)
!
!     --------------------------
!     do preprocessing for MPI
!     --------------------------
!
      CALL preproc_mpi_driver(nmort) 
!
!
!     ------------------------------------
      END IF  ! End restart/startup branch
!     ------------------------------------
!     --------------------------
!     write out grid information                                                         
!     --------------------------
!
IF (myid==0) THEN
      WRITE(iout_unit,9070) 
      WRITE(iout_unit,'(30x,a80)') title 
      WRITE(iout_unit,9070) 
      IF(restart)     THEN 
         WRITE(iout_unit,*) 'This run restarted from file:',filename 
      ENDIF 
      WRITE(iout_unit,*) 'Order of Runge-Kutta used = ',iord 
      WRITE(iout_unit,*)
      IF(verbose)     THEN
         WRITE(iout_unit,fmt='(/,1x,30x,a16,/)')'grid information' 
         DO 1100 id = 1,ngridl 
            WRITE(iout_unit,9000)id
            WRITE(iout_unit,*) '--------------------'
            IF ( method == 'dg    ' )     THEN
               WRITE(iout_unit,*)      'approximation order: ',d(id)%nc-1
            ELSE
               WRITE(iout_unit,*)      'approximation order: ',d(id)%nc-2
            END IF
            WRITE(iout_unit,9020)   (d(id)%corner(:,j),j=1,8) 
            !WRITE(iout_unit,*)      'Boundary conditions: ',((d(id)%bcond(j),'  '),j=1,6)
 1100    CONTINUE 
         IF(nmort.gt.0)     THEN 
            WRITE(iout_unit,fmt='(/,1x,25x,a22,/)') 'connection information' 
            WRITE(iout_unit,9080) 
            WRITE(iout_unit,9090) 
         ENDIF 
         DO m = 1,nmortl 
            WRITE(iout_unit,9050) m,mrtr(m)%id(1),mrtr(m)%iface(1),       &
            mrtr(m)%id(2),mrtr(m)%iface(2),mrtr(m)%orient,                &
            mrtr(m)%conforming, mrtr(m)%sign
         END DO
      END IF
ENDIF
!
!     --------------------------------------------------------------
!     compute time step if time_accurate is requested.
!     Set the number of steps and time step
!     so there is an integer number of steps betweeen movie 
!     output if movies are requested or run information, otherwise.
!     --------------------------------------------------------------
!
      IF ( time_accurate ) THEN
         dt = big 
         DO id = 1,ngridl 
            CALL dtcalc(dt,d(id)) 
         END DO
         IF (numprocs >1) THEN
            CALL gathbcast_dt(dt)
         ENDIF
         dt      = dt*stab
!         nout    = NINT(tout/dt) + 1
         npout   = NINT(tpout/dt) + 1
!         dt      = tout/dble(nout)
         dtmin   = dt
         delta_t = dt
         nsteps  = NINT((tfinal - time)/dt)
         max_steps_time = nsteps
write(*,*) 'nout etc', dt, stab, nout
         IF(max_steps_time > max_steps .and. myid == 0)     THEN
            WRITE(iout_unit,*) 'Final time takes too many steps (',max_steps_time,')'
            WRITE(iout_unit,*) 'Only the requested ',max_steps,' will be taken'
            WRITE(*,*) 'Final time takes too many steps (',max_steps_time,')'
            WRITE(*,*) 'Only the requested ',max_steps,' will be taken'
         END IF
         max_steps = MIN(max_steps_time,max_steps)
      END IF
      nout_userstuff = nout*1000
!
!     --------------
!     time step loop                                                        
!     --------------
!
      IF (myid == 0) THEN
        movie_file_no = 0
        CALL DATE_AND_TIME(date,time_str,zone,start_time_values)
        WRITE(iout_unit,*) 'Computation started: ',date, time_str
        OPEN (unit=stat_unit,file=staname)
      ENDIF

      DO k = 1,max_steps
          kiteration = k

!
	  if(myid==0) write(*,*)'Timestep: ',k, 'Sol:', d(1)%Q(1,1,1,1), "Time=", k*dt
         IF ( steady_state )     THEN ! compute the time step
            dtmin = big
            DO id = 1,ngridl
               dt = big 
               call dtcalc(dt,d(id))
               delta_t(id) = dt*stab
               dtmin = min(dtmin,delta_t(id))
            END DO
            IF ( global )     THEN   ! use the global time step
               dt = big
               DO id = 1,ngridl
                  dt = min(dt,delta_t(id))
               END DO
               IF (numprocs >1) THEN
                  CALL gathbcast_dt(dt)
               ENDIF
               delta_t = dt
               dtmin = dt
            END IF
         END IF
!
         CALL tstep(time,delta_t,d,ngridl,mrtr,nmort,resid) ! take a step

         IF ( Is_NaN(resid) .OR. Is_INF(resid) )     THEN 
              write(*,*) 'NaN appeared'
              OPEN(unit=22,file='NAN.dat')
              write(22,*) 'NaN appeared at node',myid
              WRITE(22,*) 'NAN or INF appeared at step = ',k 
              CLOSE(22)
              CALL MPI_ABORT(comm1d,nnn,ierr)
              WRITE(iout_unit,*) 'NAN or INF appeared at step = ',k 
              CLOSE(iout_unit) 
              STOP
         ENDIF 
!
         IF (drops) THEN
           IF (numprocs > 1) THEN
             CALL gather_part_matrix()
             CALL exchange_particles(drop,nprt)
           ENDIF
           CALL part_integrate(k,delta_t,d,drop,ngridl,nprt,time)
           IF (numprocs > 1) THEN
             CALL MPI_WAITALL(nreq,req(1:nreq),stata(:,1:nreq),ierr)
             CALL part_integrate_mpi(k,delta_t,d,drop,ngridl,nprt,time)
           ENDIF
           CALL WrtPartFile(fvp_unit,k,max_steps,time,drop,movie_file_no)
           IF (dropper .and. MOD(k,npout) ==0) THEN
             OPEN (unit=drop_unit,file=drop_filename,status='old')
             CALL readpart_per(drop_unit,nprt,drop,ngridl,d,delta_t)
             CLOSE(drop_unit)
           ENDIF
         ENDIF
!
         time = time + dt
!        IF( MOD( k,10 ) == 0) CALL do_user_stuff(time,ngridl,d,mrtr,nmort,stat_unit)
         CALL do_user_stuff(time,ngridl,d,mrtr,nmort,stat_unit)
         IF(itimeinflowbc==1) THEN
          DO id=1,ngridl
           CALL inflow_time_bc(id,d(id),mrtr,k,dt)
          END DO
         END IF

         IF( MOD( k,nout ) == 0)   THEN               ! do periodic operations

            IF( MOD( k,nout_userstuff ) == 0)   THEN               ! do periodic operations
              CALL do_user_stuff(time,ngridl,d,mrtr,nmort,stat_unit)
            ENDIF

            IF( save_restarts )     THEN              ! save a restart file for safety
               OPEN(UNIT=rst_unit, FILE = rstname, FORM='unformatted') 
                  CALL write_restart(rst_unit,ngrid,d,nmort,mrtr,stats,statsfluct,drop,time) ! for safety
               CLOSE(rst_unit)
            END IF

            IF( movie )     THEN
              IF (myid==0) THEN
               movie_file_no = movie_file_no + 1
               WRITE(fname, fmt='(a5,i4.4,a4)') 'movie',movie_file_no,'.plt'
               WRITE(fnamea, fmt='(a5,i4.4,a4)') 'movie',movie_file_no,'.fvt'

               OPEN(unit=plt_unit,file=fname) 
               OPEN(unit=fvt_unit,file=fnamea) 
              ENDIF

              CALL wrtdat(plt_unit,fvt_unit,fvs_unit,time,ngrid,d,.true.,stats,statsfluct,.false.)

              IF (myid==0) THEN
               CLOSE(plt_unit)
               CLOSE(fvt_unit)
              ENDIF
            END IF

            IF (statistics .and. MOD(k,nout) == 0) THEN
              CALL Stats_All(d,ngridl,stats,statsfluct,mrtr,nmort)
            END IF
!
         END IF

         IF(time_accurate .AND. time >= tfinal)     go to 1900 
!
      END DO
      k = k - 1 
!
 1900 continue 
!
!     ------------
!     finishing up                                                       
!     ------------
!
      write(*,*) 'chick finish'
!     ----------------------------
!     create metis file if desired
!     ----------------------------
!
      IF ((restart .eqv. .false.) .AND. (usemetis)) THEN
       CALL create_metis_file(d, ngrid, mrtr, nmort) 
      ENDIF
     
!     ------------------
!     get run statistics
!     ------------------
!
      IF (myid == 0) THEN     
        CALL DATE_AND_TIME(date,time_str,zone,end_time_values)
        WRITE(iout_unit,*) 'Computation completed: ',date, time_str
        WRITE(iout_unit,*) ' '
        WRITE(iout_unit,*) 'Run Statistics:'
        compute_time = elapsed_time(start_time_values,end_time_values,elapsed)
        WRITE(iout_unit,*) ' '
        WRITE(iout_unit,*) 'Number of elapsed minutes = ',compute_time
        IF(k > 0)     THEN
          WRITE(iout_unit,*) 'Time steps per hour = ',NINT(k*60/compute_time),' assuming 100% utilization'
        END IF
      ENDIF
!
!     -----------------------------------
!     print results and save restart file                                   
!     -----------------------------------
!
      IF (myid == 0) THEN
        OPEN(unit=rst_unit,file = rstname,form='unformatted') 
      ENDIF
!
      CALL write_restart(rst_unit,ngrid,d,nmort,mrtr,stats,statsfluct,drop,time) ! for safety
!
      IF (myid == 0) THEN
        CLOSE(rst_unit) 
      ENDIF

!
!     --------------------------------
!     convert solutions, JQ, back to Q
!     --------------------------------
!
      DO id = 1,ngridl
         DO nv = 1,neq
            DO l = 1,d(id)%ncg(3)
               DO m = 1,d(id)%ncg(2) 
                  DO n = 1,d(id)%ncg(1)
                     d(id)%Q(n,m,l,nv)  = d(id)%Q(n,m,l,nv)*d(id)%jacob(n,m,l)
                  END DO
               END DO 
            END DO 
         END DO
      END DO
!
!     ------------------------------------
!     do problem specific computations now                              
!     ------------------------------------
!
      CALL cleanup(time,ngridl,d, mrtr,nmort) 
!
      IF( endplot )     THEN
write(*,*) 'endplot is written now'
        IF (myid == 0) THEN
          OPEN(unit=plt_unit,file=pltname) 
          OPEN(unit=fvt_unit,file=fvtname) 
          IF ( statistics ) OPEN(unit=fvs_unit,file=fvsname)
        ENDIF
!
!        uncomment out this subroutine for writing .plt and .fvt files 
!
        CALL wrtdat(plt_unit,fvt_unit,fvs_unit,time,ngrid,d,.false.,stats,statsfluct,statistics)
!
        IF (myid == 0) THEN
          CLOSE(plt_unit)
          CLOSE(fvt_unit)
          IF ( statistics ) CLOSE(fvs_unit)
        ENDIF
      END IF
!
      CLOSE (stat_unit)
      CLOSE (iout_unit)

      CALL MPI_FINALIZE( ierr )

      STOP 'normal termination'
 
 9000 format (/,1x,'grid  ',i4,5x) 
 9020 format (1x,'corners : ',/,3(5x,1pe12.5)) 
 9060 format (12x,i2,9x,i2) 
 9070 format (/,80('*'),/) 
 9080 format (11x,'master',8x,'slave') 
 9090 format (1x,'mortar',2x,'id   side',5x,'id  side',10x,'orientation'&
     &        ,8x,'conforming',10x,'sign')
 9050 format (1x,i5,5x,i5,2x,i5,7x,i5,2x,i2,15x,i2,10x,4(l1,5x),10x,i2)
!
!

      END PROGRAM Driver_3D
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE set_io_files(rstname,outname,pltname,fvtname,fvsname,restart,staname,fvpname)
!
!.......................................................................
!     set input,output,plot and restart files from filename                          
!.......................................................................
!
      USE input_data
      IMPLICIT none
!
      LOGICAL            :: restart
      CHARACTER(LEN=32)  :: rstname,outname,pltname,fvtname,fvsname,staname,fvpname
      CHARACTER(LEN=3)   :: tag
      INTEGER            :: lp
!
      lp = index(filename,'.')

      IF(lp == 0 )     THEN   ! no extension given
         restart = .FALSE.
         outname = TRIM(filename)//'.out'
         pltname = TRIM(filename)//'.plt'
         fvtname = TRIM(filename)//'.fvt'
         fvsname = TRIM(filename)//'.fvs'
         staname = TRIM(filename)//'.sta'
         rstname = TRIM(filename)//'.rst'
         fvpname = TRIM(filename)//'.fvp'

      ELSE ! extension given could be start or restart

         tag = filename(lp+1:lp+3) 
         IF(tag.eq.'rst')     restart = .true. 
         IF(restart)     then 
            rstname = filename(1:lp-1)//'+.rst' 
         ELSE 
            rstname = filename(1:lp)//'rst' 
         ENDIF 
         outname    = filename(1:lp)//'out'
         pltname    = filename(1:lp)//'plt'
         fvtname    = filename(1:lp)//'fvt'
         fvsname    = filename(1:lp)//'fvs'
         staname    = filename(1:lp)//'sta'
         fvpname    = filename(1:lp)//'fvp'

      END IF
!
      RETURN
      END

