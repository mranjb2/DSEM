!> \file readpart.f90
!! Reading of particle input file
!///////////////////////////////////////////////////////////////////////////////////////
!////////                                                                       ////////
!////////       readpart.f90                                                    ////////
!////////                                                                       ////////
!////////       contains:                                                       ////////
!////////                                                                       ////////
!////////          SUBROUTINE readpart(iunit,nprt,drop)                         ////////
!////////          SUBROUTINE readpart_per(iunit,nprt,nprtdum,drop)                         ////////
!////////                                                                       ////////
!///////////////////////////////////////////////////////////////////////////////////////
!> Reads particles from file
   SUBROUTINE readpart(iunit,nprt,drop,ngrid,d)
!
!     date: 1/12/00
!
      USE size
      USE particle_definition
      USE domain_definition
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (particle), DIMENSION(npart) :: drop 
      TYPE (particle)                   :: dr
      TYPE (domain), DIMENSION(ngp)     :: d
!
!........................................................................................
!  read in the initial coordinates of the particle from the input file format
!
!........................................................................................

      nprt_dum = 0
      READ(iunit,*) nprt
      DO i=1,nprt  
           READ(iunit,*) a,b,c
	   dr%Xp(1)=a
	   dr%Xp(2)=b
	   dr%Xp(3)=c
 
           dr%onoff=1
           CALL find_subd(dr,ngrid,d)
           IF (dr%onoff == 1) THEN

             nprt_dum=nprt_dum+1 
             drop(nprt_dum) = dr
           ENDIF
      END DO
      nprt = nprt_dum
       
       RETURN
     END SUBROUTINE readpart
!
!///////////////////////////////////////////////////////////////////////////////////////
!> Read particle file
   SUBROUTINE readpart_per(iunit,nprt,drop,ngrid,d,dt)
!
!     date: 1/12/00
!
      USE size
      USE particle_definition
      USE domain_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (particle), DIMENSION(npart) :: drop 
      TYPE (particle)                   :: dr
      TYPE (domain), DIMENSION(ngp)     :: d
!
      DOUBLE PRECISION,DIMENSION(ngp)  :: dt
!
!........................................................................................
!  read in the initial coordinates of the particle from the input file format
!
!........................................................................................

      nprt_dum = 1
      READ(iunit,*) nprta
      DO 1000 i=1,nprta  
           READ(iunit,*) a,b,c
	   dr%Xp(1)=a
	   dr%Xp(2)=b
	   dr%Xp(3)=c
 
           dr%onoff=1
           CALL find_subd(dr,ngrid,d)
           IF (dr%onoff == 1) THEN
             DO ip=nprt_dum,npart
               IF (drop(ip)%onoff == 0) THEN
                 drop(ip)       = dr
                 drop(ip)%onoff = 1
                 ip_grid        = drop(ip)%ngrid

                 CALL part_map(d(ip_grid),drop(ip))
                 IF (d(ip_grid)%ncg(1) <7) THEN
                   CALL spec_Pinterp(d(ip_grid),drop(ip))
                 ELSE
                   CALL Lagr_Pinterp(d(ip_grid),drop(ip))
                 ENDIF
                 CALL part_initc(drop(ip),nprt)
!
                 CALL integr_onepart_euler(dt,drop(ip))
	         CALL part_map(d(ip_grid),drop(ip))
                 IF (drop(ip)%Xpmap(1) > 1.0d0 .or. drop(ip)%Xpmap(1)<0.0d0 .or.      & 
                   drop(ip)%Xpmap(2) > 1.0d0 .or. drop(ip)%Xpmap(2)<0.0d0 .or.     &
                   drop(ip)%Xpmap(3) > 1.0d0 .or. drop(ip)%Xpmap(3)<0.0d0 ) THEN 
                   CALL part_bc(d,drop,ip,ngrid)
                 END IF
!   
                 nprt_dum = ip
!
                 GO TO 1000
               ENDIF
               IF (ip==npart) write(*,*) 'out of bounds in readpart.f90'
             ENDDO
           ENDIF
1000      CONTINUE
      IF (nprt_dum > nprt) THEN 
          nprt = nprt_dum
      ENDIF
       
       RETURN
     END SUBROUTINE readpart_per
!

         
