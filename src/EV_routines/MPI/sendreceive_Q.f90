!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_q_mpi(inode,itag,Qloc)
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
       DOUBLE PRECISION,DIMENSION(nmax,nmax,neq) :: Qloc
!
       CALL MPI_SEND (                                          & 
            Qloc(:,:,:),nmax*nmax*neq, MPI_DOUBLE_PRECISION,   	& 	
            inode,itag,comm1d,ierr)
       RETURN
     END SUBROUTINE send_q_mpi
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_q_mpi(inode,itag,Qloc)
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
       DOUBLE PRECISION,DIMENSION(nmax,nmax,neq) :: Qloc
!
       CALL MPI_RECV (                                          & 
            Qloc(:,:,:),nmax*nmax*neq, MPI_DOUBLE_PRECISION,   	& 	
            inode,itag,comm1d,stat,ierr)
       RETURN
     END SUBROUTINE recv_q_mpi
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_q_mpi_noblock(inode,itag,Qloc)
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
       DOUBLE PRECISION,DIMENSION(nmax,nmax,neq) :: Qloc
!
       CALL MPI_ISEND (                                          & 
            Qloc,nmax*nmax*neq, MPI_DOUBLE_PRECISION,   	& 	
            inode,itag,comm1d,req(nreq),ierr)

      ! write(*,*) 'chick nreq',myid,nreq

       RETURN
     END SUBROUTINE send_q_mpi_noblock
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_q_mpi_noblock(inode,itag,Qloc)
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
       DOUBLE PRECISION,DIMENSION(nmax,nmax,neq) :: Qloc
!
       CALL MPI_IRECV (                                          & 
            Qloc,nmax*nmax*neq, MPI_DOUBLE_PRECISION,   	& 	
            inode,itag,comm1d,req(nreq),ierr)

       RETURN
     END SUBROUTINE recv_q_mpi_noblock
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE wait_send_q_to_mortar()
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)


!       write(*,*) 'chick req',req(1:nreq)
       CALL MPI_WAITALL(nreq,req(1:nreq),stata(:,1:nreq),ierr)

       RETURN
     END SUBROUTINE wait_send_q_to_mortar
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE wait_send_flux_to_face()
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)


       CALL MPI_WAITALL(nreq,req(1:nreq),stata(:,1:nreq),ierr)

       RETURN
     END SUBROUTINE wait_send_flux_to_face
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE wait_send_q_to_face()
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)


       nedge1 = nedge*5
       nedge2 = nedge*6

       CALL MPI_WAITALL(nedge,req(nedge1:nedge2),stata(:,nedge1:nedge2),ierr)

       RETURN
     END SUBROUTINE wait_send_q_to_face
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE wait_send_ent_to_mortar()
!
       USE size
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)


    !   write(*,*) 'chick req',req(1:nreq)
       CALL MPI_WAITALL(nreq,req(1:nreq),stata(:,1:nreq),ierr)

       RETURN
     END SUBROUTINE wait_send_ent_to_mortar
!
!//////////////////////////////////////////////////////////////