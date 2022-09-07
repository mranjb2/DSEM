!> \file metis_file.f90
!! Creation/reading of METIS files
!///////////////////////////////////////////////////////////////////////
!f///////							////////
!///////          metis_file.f90        			////////
!///////							////////
!///////	  contains: create_metis_file			////////
!///////	            read_metis_file			////////
!///////							////////
!///////////////////////////////////////////////////////////////////////

!......................................................................
!     Creates metis graph file used for partitioning
!......................................................................
!
!> @brief Creates metis graph file used for partitioning.
    SUBROUTINE create_metis_file(d,ngrid,mrtr,nmort)!

      
      USE mpi_par
      USE input_data
      USE domain_definition
      USE mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)
!
      INTEGER               :: id, iface, im, im_side, ndom
      INTEGER               :: lp
      INTEGER, DIMENSION(6) :: ids
      CHARACTER(LEN=32)     :: metisfile
      CHARACTER(LEN=3)      :: tag
	  write(*,*) 'Operation: Creating METIS .graph file'
!
!     open the graph file
!
      lp = index(filename,'.')
      IF(lp == 0 )     THEN   ! no extension given
        metisfile  = TRIM(filename)//'.graph'
      ELSE
        tag        = filename(lp+1:lp+3) 
        metisfile  = filename(1:lp)//'graph'
      ENDIF
 
 
      OPEN(UNIT=21, FILE=metisfile) 
!
!     determine number of vertices and edges (=mortars)
!     and write to file
!
      ndom = 0
      DO im = 1,nmort
        IF (nmortmpi(im,2) /= 0 .and. nmortmpi(im,3) /= 0) ndom = ndom + 1
      ENDDO

      WRITE(21,10) ngrid, ndom
!
!     determine connectivity per domain and write to file
!
      DO id = 1,ngrid
        ndom = 0
        DO iface = 1,6
          im      = d(id)%mortar(1,iface)
          im_side = d(id)%mortar(2,iface) 
          IF (im_side == 1) THEN
            im_side = 2
          ELSE
            im_side = 1
          ENDIF
          IF (nmortmpi(im, im_side+1) /= 0) THEN 
            ndom = ndom + 1
            ids(ndom) = nmortmpi(im, im_side+1)            
          ENDIF
        ENDDO
        WRITE(21, *) (ids(j), j=1,ndom)
      ENDDO

      CLOSE(21)

10    FORMAT(2i10)

      RETURN
     END SUBROUTINE

!......................................................................
!     Reads metis graph file used for partitioning
!......................................................................
!
!> @brief Reads metis graph file used for partitioning.
    SUBROUTINE read_metis_file(id, nno, ngrid)!

      
      USE mpi_par
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      INTEGER               :: lp
      CHARACTER(LEN=32)     :: metisfile
      CHARACTER(LEN=3)      :: tag
!
!     open the graph file
!
	  if(myid==0 .and. id==1) write(*,*) 'Operation: Reading METIS .graph file'
      lp = index(filename,'.')
      IF(lp == 0 )     THEN   ! no extension given
        IF (numprocs < 10) THEN
          write(metisfile,'(A,".graph.part.",I1)')   TRIM(filename), numprocs
        ELSE IF (numprocs < 100) THEN
          write(metisfile,'(A,".graph.part.",I2)')   TRIM(filename), numprocs
        ELSE IF (numprocs < 1000) THEN 
          write(metisfile,'(A,".graph.part.",I3)')   TRIM(filename), numprocs
        ELSE 
          write(metisfile,'(A,".graph.part.",I4)')   TRIM(filename), numprocs
        ENDIF
 
      ELSE
        tag        = filename(lp+1:lp+3) 
        IF (numprocs < 10) THEN
          write(metisfile,'(A,"graph.part.",I1)')   filename(1:lp), numprocs
        ELSE IF (numprocs < 100) THEN
          write(metisfile,'(A,"graph.part.",I2)')   filename(1:lp), numprocs
        ELSE IF (numprocs < 1000) THEN 
          write(metisfile,'(A,"graph.part.",I3)')  filename(1:lp), numprocs
        ELSE 
          write(metisfile,'(A,"graph.part.",I4)')  filename(1:lp), numprocs
        ENDIF
      ENDIF
 
 
      OPEN(UNIT=21, FILE=metisfile) 
!
!  Read the processor number for the domain    
!
      DO idm = 1,ngrid
         READ(21,*) nno
         IF (idm == id) THEN
           EXIT
         ENDIF
      ENDDO

      CLOSE(21)

      RETURN
     END SUBROUTINE
