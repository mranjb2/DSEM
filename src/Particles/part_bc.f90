!> \file part_bc.f90
!! Particle boundary conditions
!///////////////////////////////////////////////////////////////
!////////						////////
!////////	part_bc.f90		        	////////
!////////                                               ////////
!////////       contains:                               ////////
!////////                                               ////////
!////////	 SUBROUTINE part_bc(d,dr)		////////
!////////	 SUBROUTINE Part_Walladiab_Elast()      ////////
!////////        SUBROUTINE Part_Interf()		////////
!///////////////////////////////////////////////////////////////
!> Applies particle boundary conditions
       SUBROUTINE part_bc(d,drop,ip,ngrid)
!
!      date: 05/03/00
!     
!
!      this subroutine applies interface and boundary condition
!      to the particle in case it moves outside a domain
!
     
      USE size
      USE domain_definition
      USE particle_definition
      USE mpi_par
      USE mpi_par_part

!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)    :: d(ngp)
      TYPE(particle)    :: drop(npart)

      INTEGER           ::ngrid,pgrid
      LOGICAL           :: node_cross 
!
      node_cross = .false.
!
      pgrid = drop(ip)%ngrid
!
!     Find out if particle crosses a face or a corner point
!
      ib1 = 0
      ib2 = 0
      ib3 = 0
      ncorner_cross = 0
!
      IF (drop(ip)%Xpmap(1) < 0.0) ib1 = 6
      IF (drop(ip)%Xpmap(1) > 1.0) ib1 = 4
      IF (drop(ip)%Xpmap(2) < 0.0) ib2 = 1
      IF (drop(ip)%Xpmap(2) > 1.0) ib2 = 2
      IF (drop(ip)%Xpmap(3) < 0.0) ib3 = 3
      IF (drop(ip)%Xpmap(3) > 1.0) ib3 = 5

      IF (ib1/=0 .and. ib2 ==0 .and. ib3==0) THEN
        ib = ib1
        mortar_side = d(pgrid)%mortar(1,ib)
        IF (nmortedge(mortar_side,1) == 1) node_cross = .true.
      ELSEIF (ib2/=0 .and. ib1 ==0 .and. ib3==0) THEN
        ib = ib2
        mortar_side = d(pgrid)%mortar(1,ib)
        IF (nmortedge(mortar_side,1) == 1) node_cross = .true.
      ELSEIF (ib3/=0 .and. ib1 ==0 .and. ib2==0) THEN
        ib = ib3
        mortar_side = d(pgrid)%mortar(1,ib)
        IF (nmortedge(mortar_side,1) == 1) node_cross = .true.
      ELSE
         ncorner_cross = 1
         ms1 = d(pgrid)%mortar(1,ib1)
         ms2 = d(pgrid)%mortar(1,ib2)
         ms3 = d(pgrid)%mortar(1,ib3)
         IF (nmortedge(ms1,1) ==1 .or. nmortedge(ms2,1) ==1 .or. &
             nmortedge(ms3,1) ==1) THEN
            node_cross = .true.
         ENDIF
      ENDIF
!
! Recognize mortar, check if it crosses to another node, if so add entry
! to the local_part_mpi
!
    IF (node_cross) THEN ! this particle crosses nodes
        CALL periodic_part(drop(ip)%Xp(1),drop(ip)%Xp(2),drop(ip)%Xp(3))
        npart_send   = npart_send + 1
        IF (npart_send > npart_mpi) THEN
          write(*,*) 'number of particles send is out of bounds'
        ENDIF
        npos = (npart_send-1)*3+1

        SELECT CASE (ncorner_cross)
          CASE(0)    
            IF (d(pgrid)%bcond(ib) == 'walladiab') drop(ip)%onoff = 0 !@TODO: Check this.
            IF (d(pgrid)%bcond(ib) == 'wallisoth') drop(ip)%onoff = 0 !@TODO: CHECK THIS.
            IF (d(pgrid)%bcond(ib) == 'outflow')   drop(ip)%onoff = 0
            IF (d(pgrid)%bcond(ib) == 'inflow')    drop(ip)%onoff = 0 !@TODO: CHECK THIS
            IF (drop(ip)%onoff == 0) THEN
                npart_send = npart_send -1
                RETURN
            ENDIF
          
            IF (myid == nmortedge(mortar_side,2)) THEN
              local_part_mpi(npos)  = ip          
              local_part_mpi(npos+1)  = myid          
              local_part_mpi(npos+2)  = nmortedge(mortar_side,3)
            ELSE 
              local_part_mpi(npos)  = ip          
              local_part_mpi(npos+1)  = myid          
              local_part_mpi(npos+2)  = nmortedge(mortar_side,2)
            ENDIF
          CASE(1)    
!   this can be improved for general coding
            CALL find_processor(drop(ip)%Xp(1),drop(ip)%Xp(2),drop(ip)%Xp(3),nproc_part)
			write(*,*)'finding processor of drop',ip,'nproc_part',nproc_part,'line 110'
            IF (nproc_part< 0 .or. nproc_part> numprocs) THEN
               drop(ip)%onoff = 0
			   write(*,*)'drop',ip,'turned off line 113'
               RETURN
            ENDIF
            local_part_mpi(npos)  = ip          
            local_part_mpi(npos+1)  = myid          
            local_part_mpi(npos+2)  =  nproc_part
       END SELECT
    ELSE
       SELECT CASE(ncorner_cross) 
      
         CASE(0)    
           IF (d(pgrid)%bcond(ib) == 'interface') then
		   		CALL Part_Interf(drop(ip),ngrid,d)
		   		!write(*,*)'drop',ip,'is at interface line 124'
			end if
           IF (d(pgrid)%bcond(ib) == 'walladiab') drop(ip)%onoff = 0
           IF (d(pgrid)%bcond(ib) == 'wallisoth') drop(ip)%onoff = 0
           IF (d(pgrid)%bcond(ib) == 'outflow')   drop(ip)%onoff = 0
           IF (d(pgrid)%bcond(ib) == 'inflow')    drop(ip)%onoff = 0
          CASE(1)    
		   !write(*,*)'drop',ip,'is at interface line 128'
           CALL Part_Interf(drop(ip),ngrid,d)
       END SELECT

    ENDIF
!        
      RETURN
   END SUBROUTINE part_bc
  
!
!///////////////////////////////////////////////////////////////////////////////
!
!> Determines if particles crossed an interface
      SUBROUTINE Part_Interf(dr,ngrid,d)
      
      USE size       
      USE domain_definition
      USE particle_definition

!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)     :: d(ngp)
      TYPE (particle)   :: dr
      INTEGER           :: ngrid,pgrid

!
!      Find new subdomain
!
       CALL periodic_part(dr%Xp(1),dr%Xp(2),dr%Xp(3))
!
       CALL find_subd(dr,ngrid,d)
       IF (dr%onoff == 1) THEN
        pgrid = dr%ngrid
        CALL part_map(d(pgrid),dr)
       ENDIF

      RETURN
   END SUBROUTINE part_Interf
!
!///////////////////////////////////////////////////////////////////////////////
!> Periodic boundary condition on particles
       SUBROUTINE periodic_part(x1,x2,x3) !added x3 dimension
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
       double precision :: lx,ly,lz
       
       lx = 2.0d0
       ly = 2.0d0
       lz = 2.0d0
!
! set periodic bc
       IF (x1 > lx) THEN
           x1 = x1 - lx
       ENDIF
       IF (x1 < 0.0d0) THEN
           x1 = x1 + lx
       ENDIF
       IF (x2 > ly) THEN
           x2 = x2 - ly
       ENDIF
       IF (x2 < 0.0d0) THEN
           x2 = x2 + ly
       ENDIF
	   !The following code snippet was added to fix periodic particles
	   !write(*,*)'x3',x3
	   if (x3 > lz) then
		   x3 = x3 - lz
	   end if
	   if (x3 < 0.0d0) then
		   x3 = x3 + lz
	   end if
	   !end of the addition
!
      
       RETURN
     END SUBROUTINE periodic_part
