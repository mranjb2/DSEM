!> \file
!! Writes particle data to VTK
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE WrtPartFile(iunit,k,max_step,time,drop,movie_file_no)
!
      USE particle_definition
      USE input_data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      TYPE( particle ) :: drop(npart)
!
!
      CHARACTER(LEN=32) :: fname
!

      IF (drmov) THEN
         IF( MOD( k,nout ) == 0)   THEN               ! do periodic operations
           IF (myid==0) THEN
             movie_file_no = movie_file_no + 1
             WRITE(fname, fmt='(a8,i5.5,a4)') 'particle',movie_file_no,'.vtk'
             OPEN(unit=iunit,file=fname)
           ENDIF
           CALL WrtPartFilea(iunit,time,drop)
           IF (myid==0) CLOSE(iunit)
         ENDIF
      ELSEIF (drendpl) THEN
         IF (k  > max_step-2) THEN
           IF (myid==0) THEN
             movie_file_no = movie_file_no + 1
             WRITE(fname, fmt='(a8,i5.5,a4)') 'particle',movie_file_no,'.vtk'
             OPEN(unit=iunit,file=fname)
           ENDIF
           CALL WrtPartFilea(iunit,time,drop)
           IF (myid==0) CLOSE(iunit)
         ENDIF  
      ENDIF 

      RETURN
     END SUBROUTINE WrtPartFile
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE WrtPartFilea(iunit,time,drop)
!
      USE particle_definition
      USE mpi_par
      USE mpi
	  use size
	  use part_par
	  use constants
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      TYPE( particle ) :: drop(npart)
      TYPE( particle ) :: dropa(npart)
!
!
!  Determine the total number of particles per processor, reduce later 
!
  
      nprta = 0
      DO i=0,numprocs-1
        IF (myid==i) THEN
          DO ip=1,npart
            IF (drop(ip)%onoff==1) nprta=nprta+1
          ENDDO
        ENDIF
      ENDDO
!
!   Send the droplets from other processors to the root processer and write
!   to file
!
      IF (numprocs>1) THEN

        CALL MPI_REDUCE(nprta,nprtb,1,MPI_INTEGER, MPI_SUM,0, &
             comm1d, ierr)

        DO i=0,numprocs-1
           IF (i==0) THEN
              IF (myid==0) WRITE(iunit,20) nprtb
              IF (myid ==0) CALL wrtpartLoc(iunit,drop,time)
           ELSE
              IF (myid==i) THEN
                CALL send_drop_plot(drop,i+100000,0)
              ELSEIF (myid==0) THEN
               CALL recv_drop_plot(dropa,i+100000,i)
               CALL wrtpartLoc(iunit,dropa,time)
              ENDIF
!                
           ENDIF
        ENDDO
!
       ELSE   ! one processor 
	   !write(*,*)'nprt',nprt
	   !write(*,*)'npart',npart
	   !write(*,*)'nprta',nprta
	   !write(*,*)'npartb',npartb
	   !--- write the header portion of the vtk file
         WRITE(iunit,20) '# vtk DataFile Version 3.0'
         WRITE(iunit,20) 'Three-Dimensional Particle Tracking'
		 write(iunit,20) 'ASCII'
		 write(iunit,20) 'DATASET POLYDATA'
		 write(iunit,20) 'FIELD FieldData 1'
		 write(iunit,20) 'TIME 1 1 double'
		 write(iunit,*) time
		 write(iunit,*)
	   !--- write the point data for the particles	 
		 write(iunit,22) 'POINTS ',nprta,' float'
         CALL wrtpartLoc(iunit,drop,time)
		 write(iunit,*)
	   !--- write some scalar data for the particles
	   	 write(iunit,21) 'POINT_DATA ',nprta
		 write(iunit,20) 'SCALARS temperature float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartTemp(iunit,drop,time)
		 write(iunit,*)
	     write(iunit,20) 'SCALARS nu float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartVis(iunit,drop,time)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS omega float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartOmega(iunit,drop,time)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS scalar1 float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartScalar(iunit,drop,time,1)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS scalar2 float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartScalar(iunit,drop,time,2)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS scalar3 float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartScalar(iunit,drop,time,3)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS scalar4 float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartScalar(iunit,drop,time,4)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS scalar5 float'
		 write(iunit,20) 'LOOKUP_TABLE default'
		 call wrtpartScalar(iunit,drop,time,5)
		 write(iunit,*)
		 write(iunit,20) 'SCALARS weight float'
 		 write(iunit,20) 'LOOKUP_TABLE default'
 		 call wrtpartWeight(iunit,drop,time)
 		 write(iunit,*)
         write(iunit,20) 'SCALARS nu_x float'
         write(iunit,20) 'LOOKUP_TABLE default'
         call wrtpartNu_d(iunit,drop,time,1)
         write(iunit,*)
         write(iunit,20) 'SCALARS nu_y float'
         write(iunit,20) 'LOOKUP_TABLE default'
         call wrtpartNu_d(iunit,drop,time,2)
         write(iunit,*)
         write(iunit,20) 'SCALARS nu_z float'
         write(iunit,20) 'LOOKUP_TABLE default'
         call wrtpartNu_d(iunit,drop,time,3)
         write(iunit,*)
	   !--- write some vector data for particles
	     write(iunit,20) 'VECTORS velocity float'
		 call wrtpartVel(iunit,drop,time)
		 write(iunit,*) 
		
		 
       ENDIF  ! if loop on numprocessors
!
!

20    FORMAT(A)
21    format(A,I9)
22    format(A,I9,A)
      RETURN
      END SUBROUTINE  WrtPartFilea
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE wrtpartLoc(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
      USE particle_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE(particle)                       :: drop(npart) 
      
       
      DO ip=1,npart
          IF (drop(ip)%onoff==1) THEN
          WRITE(iunit,10) drop(ip)%Xp(1),drop(ip)%Xp(2),drop(ip)%Xp(3)
          ENDIF
      END DO

10    FORMAT(6f24.16)
20    FORMAT(i9)

      RETURN
   END SUBROUTINE wrtpartLoc
   
!
!///////////////////////////////////////////////////////////////////////
!
  SUBROUTINE wrtpartTemp(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
     USE particle_definition
!
     IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
     TYPE(particle)                       :: drop(npart) 
      
       
     DO ip=1,npart
         IF (drop(ip)%onoff==1) THEN
         WRITE(iunit,10) drop(ip)%Tfp
         ENDIF
     END DO

10    FORMAT(6f24.16)
20    FORMAT(i9)

     RETURN
  END SUBROUTINE wrtpartTemp
!///////////////////////////////////////////////////////////////////////
!
  SUBROUTINE wrtpartVis(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
     USE particle_definition
!
     IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
     TYPE(particle)                       :: drop(npart) 


     DO ip=1,npart
         IF (drop(ip)%onoff==1) THEN
         WRITE(iunit,10) drop(ip)%nu_t
         ENDIF
     END DO

10    FORMAT(6f24.16)
20    FORMAT(i9)

     RETURN
  END SUBROUTINE wrtpartVis
!///////////////////////////////////////////////////////////////////////
!
SUBROUTINE wrtpartNu_d(iunit,drop,time,n)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
   USE particle_definition
!
   IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
   TYPE(particle)                       :: drop(npart) 


   DO ip=1,npart
       IF (drop(ip)%onoff==1) THEN
       WRITE(iunit,10) drop(ip)%nu_d(n)
       ENDIF
   END DO

10    FORMAT(6f24.16)
20    FORMAT(i9)

   RETURN
END SUBROUTINE wrtpartNu_d
!///////////////////////////////////////////////////////////////////////
!
  SUBROUTINE wrtpartOmega(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
     USE particle_definition
!
     IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
     TYPE(particle)                       :: drop(npart) 


     DO ip=1,npart
         IF (drop(ip)%onoff==1) THEN
         WRITE(iunit,10) drop(ip)%omega_m
         ENDIF
     END DO

10    FORMAT(6f24.16)
20    FORMAT(i9)

     RETURN
  END SUBROUTINE wrtpartOmega
!///////////////////////////////////////////////////////////////////////
!
 SUBROUTINE wrtpartVel(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
    USE particle_definition
!
    IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
    TYPE(particle)                       :: drop(npart) 
      
       
    DO ip=1,npart
        IF (drop(ip)%onoff==1) THEN
        WRITE(iunit,10) drop(ip)%Vp(1),drop(ip)%Vp(2),drop(ip)%Vp(3)
        ENDIF
    END DO

10    FORMAT(6f24.16)
20    FORMAT(i6)

    RETURN
 END SUBROUTINE wrtpartVel
 
  SUBROUTINE wrtpartScalar(iunit,drop,time,which)
 !
 !     date: 01/17/00
 !
 !     Write out the the particle solution for plotting by fieldview 
 !
 !
     USE particle_definition
 !
     IMPLICIT DOUBLE PRECISION (a-h,o-z)
 !
     TYPE(particle)                     :: drop(npart) 
     integer                            :: which


     DO ip=1,npart
         IF (drop(ip)%onoff==1) THEN
         WRITE(iunit,10) drop(ip)%scalar(which)
         ENDIF
     END DO

 10    FORMAT(6f24.16)
 20    FORMAT(i6)

     RETURN
  END SUBROUTINE wrtpartScalar

 SUBROUTINE wrtpartWeight(iunit,drop,time)
!
!     date: 01/17/00
!
!     Write out the the particle solution for plotting by fieldview 
!
!
    USE particle_definition
!
    IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
    TYPE(particle)                     :: drop(npart) 


    DO ip=1,npart
        IF (drop(ip)%onoff==1) THEN
        WRITE(iunit,10) drop(ip)%w
        ENDIF
    END DO

10    FORMAT(6f24.16)
20    FORMAT(i6)

    RETURN
 END SUBROUTINE wrtpartWeight
   
