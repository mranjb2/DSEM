!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_domain(d,nsend,inode,itag)
      

      USE size
      USE domain_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(domain) :: d
!
!   Send the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 +1
      ncount1 = nx*ny*nz*neq
      ncount2 = nmax*4*neq
      IF (nsend==1) THEN
        CALL MPI_SEND (             		     	            &
             d%type,1, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%node(:),26, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%orientation(:),6, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%mortar(:,:),12, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%which_surface(:),6, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%nc(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%ncg(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%Qlgg,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%Qglg,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
      ELSEIF (nsend==2) THEN
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%Qggl,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%f,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%g,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%h,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%g_Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%material_id,1, MPI_INTEGER,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%cxg,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
      ELSEIF(nsend==3) THEN
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%xg,3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%jacob,nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%gmet,3*3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%gmetg,3*3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%corner,3*8, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%cx,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%ibtype,6, MPI_INTEGER,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%cx,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%bcond,6*9, MPI_CHARACTER,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%domain_type,8, MPI_CHARACTER,inode,itag,     &
             comm1d,ierr)
      ENDIF

      itag = ndum

      RETURN 
      END SUBROUTINE send_domain
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_domain(d,nsend,inode,itag)
      

      USE size
      USE domain_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(domain) :: d
!
!   Receive the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 + 1
      ncount1 = nx*ny*nz*neq
      ncount2 = nmax*4*neq
      IF (nsend==1) THEN
        CALL MPI_RECV (             		     	            &
             d%type,1, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%node(:),26, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%orientation(:),6, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%mortar(:,:),12, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%which_surface(:),6, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%nc(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%ncg(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%Qlgg,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%Qglg,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
      ELSEIF (nsend==2) THEN
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%Qggl,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%f,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%g,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%h,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%g_Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%material_id,1, MPI_INTEGER,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%cxg,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
      ELSEIF(nsend==3) THEN
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%xg,3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%jacob,nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%gmet,3*3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%gmetg,3*3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%corner,3*8, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%cx,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%ibtype,6, MPI_INTEGER,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%cx,nmax*3, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%bcond,6*9, MPI_CHARACTER,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%domain_type,8, MPI_CHARACTER,inode,itag,     &
             comm1d,stat,ierr)
      ENDIF

      itag = ndum

      RETURN
      END SUBROUTINE recv_domain
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_mrtr(mrtr,inode,itag)
      

      USE size
      USE mortar_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(mortar) :: mrtr

!
!   Send the mortar, array by array
!
        ndum = itag
        itag = ngp*50 + (itag-1)*20
        CALL MPI_SEND (             		     	            		&
             mrtr%lenmortar(:),2, MPI_INTEGER,inode,itag, 	                    		& 
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%len(:,:),4, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%id(:),2, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%iface(:),2, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%orient,1, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%sign,1, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%nsign,1, MPI_INTEGER,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%n_hat(:,:,:),3*nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%norm(:,:),nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     		&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%Q(:,:,:,:),nmax*nmax*neq*2, MPI_DOUBLE_PRECISION,inode,itag,     	&
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%fv(:,:,:,:),nmax*nmax*neq*2, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%xg(:,:,:),3*nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            		&
             mrtr%conforming(:,:),4, MPI_LOGICAL,inode,itag,    &
             comm1d,ierr)
        
!
        itag = ndum
      RETURN
      END SUBROUTINE send_mrtr
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_mrtr(mrtr,inode,itag)
      

      USE size
      USE mortar_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(mortar) :: mrtr
!
!   Send the mortar, array by array
!
       ndum = itag
       itag = ngp*50 + (itag-1)*20
        CALL MPI_RECV (             		     	            		&
             mrtr%lenmortar(:),2, MPI_INTEGER,inode,itag, 	                    		& 
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%len(:,:),4, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%id(:),2, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%iface(:),2, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%orient,1, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%sign,1, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%nsign,1, MPI_INTEGER,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%n_hat(:,:,:),3*nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%norm(:,:),nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     		&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%Q(:,:,:,:),nmax*nmax*neq*2, MPI_DOUBLE_PRECISION,inode,itag,     	&
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%fv(:,:,:,:),nmax*nmax*neq*2, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%xg(:,:,:),3*nmax*nmax, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            		&
             mrtr%conforming(:,:),4, MPI_LOGICAL,inode,itag,    &
             comm1d,stat,ierr)
!

         itag = ndum
      RETURN
    END SUBROUTINE recv_mrtr
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_domain_plot(d,nsend,inode,itag,stats,statis)
      

      USE size
      USE stats_definition
      USE domain_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(domain) :: d
      TYPE (statfl)  :: stats 
!
      LOGICAL :: statis
!
!   Send the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 +1+ 30
      ncount1 = nx*ny*nz*neq
      ncount2 = nmax*4*neq
      IF (nsend==1) THEN
        CALL MPI_SEND (             		     	            &
             d%nc(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%ncg(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%xg,3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             d%jacob,nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,ierr)
        IF (statis) THEN
          itag = itag+1
           CALL MPI_SEND (             		     	            &
             stats%Q_av,nx*ny*nz*(neq+nav), MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,ierr)
        ENDIF
      ENDIF

      itag = ndum

      RETURN 
      END SUBROUTINE send_domain_plot
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_domain_plot(d,nsend,inode,itag,stats,statis)
      

      USE size
      USE domain_definition
      USE stats_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE(domain) :: d
      TYPE (statfl)  :: stats 
!
      LOGICAL       :: statis
!
!   Receive the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 + 1+ 30
      ncount1 = nx*ny*nz*neq
      ncount2 = nmax*4*neq
      IF (nsend==1) THEN
        CALL MPI_RECV (             		     	            &
             d%nc(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%ncg(:),3, MPI_INTEGER,inode,itag, 		    &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%Q,ncount1, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%xg,3*nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             d%jacob,nx*ny*nz, MPI_DOUBLE_PRECISION,inode,itag,     &
             comm1d,stat,ierr)
        IF (statis) THEN
          itag = itag+1
          CALL MPI_RECV (             		     	            &
             stats%Q_av,nx*ny*nz*(neq+nav), MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,stat,ierr)
        ENDIF
      ENDIF

      itag = ndum

      RETURN
      END SUBROUTINE recv_domain_plot

!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_stats_plot(stats,nsend,inode,itag)

      USE size
      USE stats_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE (statfl)  :: stats 
!
!   Send stats
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 +1 + 30
      IF (nsend==3) THEN
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             stats%Q_av,nx*ny*nz*(neq+nav), MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,ierr)
      ENDIF

      itag = ndum

      RETURN 
      END SUBROUTINE send_stats_plot
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_stats_plot(stats,nsend,inode,itag)
      

      USE size
      USE stats_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
      TYPE (statfl) :: stats
!
!   Receive the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 + 1+ 30
      IF (nsend==3) THEN
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             stats%Q_av,nx*ny*nz*(neq+nav), MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,stat,ierr)
      ENDIF

      itag = ndum

      RETURN
      END SUBROUTINE recv_stats_plot
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE send_statsfluct_plot(statsfluct,nsend,inode,itag)

      USE size
      USE stats_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)


      TYPE (statfluct)  :: statsfluct 
!
!   Send stats
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 +1 + 30
      IF (nsend==2) THEN
        itag = itag+1
        CALL MPI_SEND (             		     	            &
             statsfluct%Q_fluct,nx*ny*nz*nfluct, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,ierr)
      ENDIF

      itag = ndum

      RETURN 
      END SUBROUTINE send_statsfluct_plot
!
!//////////////////////////////////////////////////////////////
!
      SUBROUTINE recv_statsfluct_plot(statsfluct,nsend,inode,itag)
      

      USE size
      USE stats_definition
      USE mpi_par
      USE mpi
  
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

!
      TYPE (statfluct) :: statsfluct
!
!   Receive the domain, array by array
!
      ndum  = itag
      itag  = (itag -1)*40 + (nsend-1)*10 + 1+ 30
      IF (nsend==2) THEN
        itag = itag+1
        CALL MPI_RECV (             		     	            &
             statsfluct%Q_fluct,nx*ny*nz*nfluct, MPI_DOUBLE_PRECISION,inode,itag,          &
             comm1d,stat,ierr)
      ENDIF

      itag = ndum

      RETURN
      END SUBROUTINE recv_statsfluct_plot
!
!//////////////////////////////////////////////////////////////
!
