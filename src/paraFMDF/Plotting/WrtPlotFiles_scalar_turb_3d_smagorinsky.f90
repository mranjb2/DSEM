!> \file
!! Plotting routines
!///////////////////////////////////////////////////////////////////////
!///////                            ////////
!///////                            ////////
!///////          DISCONTINUOUS GALERKIN VERSION        ////////
!///////                            ////////
!///////      contains: wrtdat              ////////
!///////            fill_plot_arrays            ////////
!///////                            ////////
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE wrtdat(plt_unit,fvt_unit,fvs_unit,time,ngrid,d,scale,stats,statsfluct,statis)
!
!......................................................................
!     date: 11/13/98
!     Write out the solution for plotting by tecplot
!
!     (see physlib for variable mappings)
!......................................................................
!
!
      USE domain_definition
      USE stats_definition
      USE Order_Matrices
      USE physics
      USE constants
      USE mpi_par
      USE input_data
      USE mpi
      USE size
      use turbulence
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
!      INCLUDE "mpif.h"

!
      TYPE(domain)                                          :: d(ngp)
      TYPE(domain)                                          :: da
      TYPE (statfl)  , DIMENSION(ngp)                       :: stats 
      TYPE (statfl)                                         :: statsa, statsb
      TYPE (statfluct)                                      :: statsfluct(ngp)
      TYPE (statfluct)                                      :: statsflucta
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,neq+nav)    :: Q
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,9)          :: Tensor
      double precision,dimension(nx+2,ny+2,nz+2,10)         :: s
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)            :: x,y,z,pCount,nu,r,nux,nuy,nuz
      DOUBLE PRECISION                                      :: p,mach_no,u,v,w,rho,temp
      LOGICAL                                               :: scale,statis
      INTEGER                                               :: plt_unit,fvt_unit,fvs_unit,plot_unit,handshake
      INTEGER,DIMENSION(3)                                  :: ncg
      !CHARACTER, DIMENSION(80) :: title

!-----
    IF (myid ==0) THEN  
      !len = LEN_TRIM(title)
!
!     ----------------------------
!     write out header information
!     ----------------------------
!
      nvar = 6
      WRITE(plt_unit,*) ' TITLE = "test"'
      IF (statistics) THEN
       IF (average) THEN
         WRITE(plt_unit,*) ' VARIABLES = "x","y","z","rho_av","u_fav","v_fav","w_fav", &
                          "p_av","u_rav","v_rav","w_rav","T_av"'
       ElSEIF (rms) THEN
         WRITE(plt_unit,*) ' VARIABLES = "x","y","z","uu","vv","ww","uv","uw","vw","ut","vt","wt"'
       END IF   
      ELSE
        write(plt_unit,*) 'TEXT X=50, Y=50, T="',time,'"'
        IF (isfmdf) THEN
           WRITE(plt_unit,*) ' VARIABLES = "x","y","z","p","u","v","w","rho","T",  &
                           "S1","S2","S3","S4","S5","pCount","react","nu_t","nu_x","nu_y","nu_z","EV","NUT"'
        ELSE IF (nsc > 0) THEN
           WRITE(plt_unit,*) ' VARIABLES = "x","y","z","u","v","w","rho","T",&
                           "react","EV","NUT","T2","F","O","H2O","CO2","N2"'
        ELSE
           WRITE(plt_unit,*) ' VARIABLES = "x","y","z","u","v","w","rho","T",  &
                           "react","EV","NUT","T1","T2","T3","T4","T5","T6","T7","T8","T9"'
        END IF
      END IF

      IF (statis) WRITE(fvs_unit,*) ngrid
    ENDIF
      DO id = 1,ngrid
        idl = ngridmpi(id,1)
        IF (myid == 0 .AND. ngridmpi(id,2) == 0) THEN
          ncg(:) = d(idl)%ncg(:)
        ELSEIF (numprocs >1) THEN
          nsnode = ngridmpi(id,2)                      ! node from which one is sending
          nrnode = 0                                   ! receiving node (master in this case)
          itag   = id*10
          IF (myid /= 0 .and. myid==nsnode) THEN
            CALL MPI_SEND (                                         &
               d(idl)%ncg(:),3,MPI_INTEGER,nrnode,itag,             &
               comm1d,ierr)
          ELSEIF (myid == 0 .and. nsnode /= 0) THEN
            CALL MPI_RECV (                                         &
               ncg(:),3, MPI_INTEGER,nsnode,itag,                   &
               comm1d,stat,ierr)
          ENDIF
        ENDIF
!
        IF (myid == 0) THEN
          !WRITE(plt_unit,*)  ncg(1)+2, ncg(2)+2, ncg(3)+2
          !WRITE(fvt_unit,*)  ncg(1)+2, ncg(2)+2, ncg(3)+2,nvar
          IF (statis) WRITE(fvs_unit,*)  ncg(1)+2, ncg(2)+2, ncg(3)+2,nvar
        ENDIF
      END DO
!
      CALL zero_domain(da) 
      IF (statis) CALL zero_stats(statsa,statsflucta)


      DO id = 1,ngrid
        !write(*,*) 'Write grid with id', id, ' of total grid', ngrid
!
!    Write plt file
!
         idl  = ngridmpi(id,1)
         IF (myid == 0 .and. ngridmpi(id,2) == 0) THEN
            da = d(idl)
            IF (statis) THEN
             statsa = stats(idl)
             IF (rms) THEN 
               statsa%Q_av(:,:,:,1:6) = statsfluct(idl)%Q_fluct(:,:,:,8:13)
               statsa%Q_av(:,:,:,7:8) = statsfluct(idl)%Q_fluct(:,:,:,1:2)
               statsa%Q_av(:,:,:,9) = statsfluct(idl)%Q_fluct(:,:,:,14)
             END IF
            ENDIF
         ELSEIF (numprocs >1) THEN
            nsnode = ngridmpi(id,2)                      ! node from which one is sending
            nrnode = 0                                   ! receiving node (master in this case)
            itaga   = id
            DO ns=1,1
              IF (myid /= 0 .and. myid==nsnode) THEN
                IF (statistics) THEN
                 statsb = stats(idl)
                 IF (rms) THEN 
                  statsb%Q_av(:,:,:,1:6) = statsfluct(idl)%Q_fluct(:,:,:,8:13)
                  statsb%Q_av(:,:,:,7:8) = statsfluct(idl)%Q_fluct(:,:,:,1:2 )
                  statsb%Q_av(:,:,:,9  ) = statsfluct(idl)%Q_fluct(:,:,:,14  )
                 END IF
                END IF
                 CALL send_domain_plot(d(idl),ns,nrnode,itaga,statsb,statis)
              ELSEIF (myid == 0 .and. nsnode /= 0) THEN
                 CALL recv_domain_plot(da,ns,nsnode,itaga,statsa,statis)
              ENDIF
            ENDDO

         ENDIF

 
       IF (myid ==0) THEN

!     --------------
!     reset pointers
!     --------------
!
         DO k = 1,num_ordersx
            IF ( order_mapx(k) == da%nc(1) )     THEN
               da%bx  => bx (:,:,k)
               EXIT
            END IF
         END DO
         DO k = 1,num_ordersy
            IF ( order_mapy(k) == da%nc(2) )     THEN
               da%by  => by (:,:,k)
               EXIT
            END IF
         END DO
         DO k = 1,num_ordersz
            IF ( order_mapz(k) == da%nc(3) )     THEN
               da%bz  => bz (:,:,k)
               EXIT
            END IF
         END DO

         IF (statistics) THEN
          CALL fill_stat_arrays(da,statsa,Q,x,y,z,scale)
         ELSE 
          CALL fill_plot_arrays(da,Q,x,y,z,scale)
          call fill_scal_arrays(da,s,x,y,z,scale)
          call fill_num_arrays(da,pCount,x,y,z,scale)
          call fill_nu_arrays(da,id,nu,nux,nuy,nuz,x,y,z,scale)
          call fill_r_arrays(da,r,x,y,z,scale)
          call fill_tensor_arrays(da,Tensor,x,y,z,scale)
         ENDIF
!
         WRITE (plt_unit,*) ' ZONE I=',da%ncg(1)+2,' J=',da%ncg(2)+2, &
         ' K=',da%ncg(3)+2,' , F = POINT SOLUTIONTIME=',time 

         IF (statistics) THEN
           IF (average) THEN 
             DO l = 1,da%ncg(3) + 2
               DO m = 1,da%ncg(2) + 2
                  DO n = 1,da%ncg(1) + 2
                    rho_av = Q(n,m,l,1)
                    u_fav = Q(n,m,l,2)/Q(n,m,l,1)
                    v_fav = Q(n,m,l,3)/Q(n,m,l,1)  
                    w_fav = Q(n,m,l,4)/Q(n,m,l,1)
                    p_av  = Q(n,m,l,5)
                    u_rav = Q(n,m,l,6)
                    v_rav = Q(n,m,l,7)
                    w_rav = Q(n,m,l,8)
                    t_av  = Q(n,m,l,9)
                    WRITE(plt_unit,50) x(n,m,l), y(n,m,l), z(n,m,l), &
                         rho_av,u_fav,v_fav,w_fav,p_av,u_rav,v_rav,w_rav,t_av
                  END DO
               END DO
             END DO  
           ELSEIF (rms) THEN
             DO l = 1,da%ncg(3) + 2
               DO m = 1,da%ncg(2) + 2
                  DO n = 1,da%ncg(1) + 2
                    uu = Q(n,m,l,1)
                    vv = Q(n,m,l,2)
                    ww = Q(n,m,l,3)
                    uv = Q(n,m,l,4)
                    uw = Q(n,m,l,5)
                    vw = Q(n,m,l,6)
                    ut = Q(n,m,l,7)
                    vt = Q(n,m,l,8)
                    wt = Q(n,m,l,9)
                    WRITE(plt_unit,50) x(n,m,l), y(n,m,l), z(n,m,l), &
                         uu,vv,ww,uv,uw,vw,ut,vt,wt
                  END DO
               END DO
             END DO
           ENDIF 
         ELSE
             DO l = 1,da%ncg(3) + 2
               DO m = 1,da%ncg(2) + 2
                  DO n = 1,da%ncg(1) + 2
                    u    = Q(n,m,l,2)/Q(n,m,l,1)
                    v    = Q(n,m,l,3)/Q(n,m,l,1)
                    w    = Q(n,m,l,4)/Q(n,m,l,1)
                    rho  = Q(n,m,l,1)
                    p    = (gamma-1.d0)*(Q(n,m,l,5) - 0.5d0*rho*(u**2 + v**2 + w**2))
                    temp = p*gamma*mach*mach/rho
                    
                   
                    
                    IF (isfmdf) THEN
                       WRITE(plt_unit,40) x(n,m,l), y(n,m,l), z(n,m,l),        &
                       &                  p, u, v, w, rho, temp, s(n,m,l,1),   &
                       &                  s(n,m,l,2), s(n,m,l,3), s(n,m,l,4),  &
                       &                  s(n,m,l,5), pCount(n,m,l),r(n,m,l),  &
                       &                  nu(n,m,l),nux(n,m,l),nuy(n,m,l),nuz(n,m,l),s(n,m,l,6),s(n,m,l,7)
                    ELSE IF (nsc > 0) THEN
                       WRITE(plt_unit,40) x(n,m,l), y(n,m,l), z(n,m,l),        &
                       &                  u, v, w, rho, temp,r(n,m,l),         &
                       &                  s(n,m,l,6),s(n,m,l,7)
                    ELSE
                       WRITE(plt_unit,40) x(n,m,l), y(n,m,l), z(n,m,l),         &
                       &                  u, v, w, rho, temp,r(n,m,l),          &
                       &                  s(n,m,l,6),s(n,m,l,7),Tensor(l,m,n,1),&
                                          Tensor(l,m,n,2),Tensor(l,m,n,3),Tensor(l,m,n,4),&
                                          Tensor(l,m,n,5),Tensor(l,m,n,6),Tensor(l,m,n,7),&
                                          Tensor(l,m,n,8),Tensor(l,m,n,9)
                    END IF
                  END DO
                END DO
              END DO
         END IF 

!
!     Write fvt file
!

      ENDIF
     
      IF (numprocs >1) THEN
      IF (myid==0) handshake = 1
      CALL MPI_BCAST(handshake,1,MPI_INTEGER,0,comm1d,ierr)
      ENDIF
 
      END DO
30    FORMAT(e20.8)
40    FORMAT(31e20.8)
50    FORMAT(12e20.8)
60    FORMAT(9e20.8)
!
!
      RETURN
   END SUBROUTINE wrtdat
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE fill_plot_arrays(d,Q,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,neq) :: Q
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
      IF ( scale )     THEN ! scale out the jacobians for movies
         DO nv = 1,neq
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Q(n+1,m+1,l+1,nv) = d%Q(n,m,l,nv)*d%jacob(n,m,l)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO nv = 1,neq
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Q(n+1,m+1,l+1,nv) = d%Q(n,m,l,nv)
                  END DO
               END DO
            END DO
         END DO
      END IF

      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      DO nv = 1,neq
         CALL interp_int_to_faces(Q(:,:,:,nv),d%ncg,d%nc,d%bx,d%by,d%bz)
      END DO
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      DO nv = 1,neq
         CALL interp_face_to_edge(Q(:,:,:,nv),d%ncg,d%nc,d%bx,d%bz)
      END DO
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      DO nv = 1,neq
         CALL interp_edge_to_corner(Q(:,:,:,nv),d%ncg,d%nc,d%bx)
      END DO
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END


   SUBROUTINE fill_r_arrays(d,r,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: r
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
         !DO nv = 1,10
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     r(n+1,m+1,l+1) = d%react(n,m,l)
                  END DO
               END DO
            END DO
         !END DO


      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      !DO nv = 1,10
         CALL interp_int_to_faces(r(:,:,:),d%ncg,d%nc,d%bx,d%by,d%bz)
      !END DO
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      !DO nv = 1,10
         CALL interp_face_to_edge(r(:,:,:),d%ncg,d%nc,d%bx,d%bz)
      !END DO
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      !DO nv = 1,10
         CALL interp_edge_to_corner(r(:,:,:),d%ncg,d%nc,d%bx)
      !END DO
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END
                                           
!
!///////////////////////////////////////////////////////////////////////


   SUBROUTINE fill_num_arrays(d,pCount,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: pCount
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
         !DO nv = 1,10
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     pCount(n+1,m+1,l+1) = dble(d%pCount(n,m,l))
                  END DO
               END DO
            END DO
         !END DO


      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      !DO nv = 1,10
         CALL interp_int_to_faces(pCount(:,:,:),d%ncg,d%nc,d%bx,d%by,d%bz)
      !END DO
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      !DO nv = 1,10
         CALL interp_face_to_edge(pCount(:,:,:),d%ncg,d%nc,d%bx,d%bz)
      !END DO
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      !DO nv = 1,10
         CALL interp_edge_to_corner(pCount(:,:,:),d%ncg,d%nc,d%bx)
      !END DO
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END
                                           
!
!///////////////////////////////////////////////////////////////////////

   SUBROUTINE fill_scal_arrays(d,s,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,10)  :: s
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
         DO nv = 1,10
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                  	 IF (nv==6) THEN
                  	    s(n+1,m+1,l+1,nv) = d%muhg(n,m,l)
                  	 ELSE IF (nv==7) THEN
                  	    s(n+1,m+1,l+1,nv) = d%nu_t(n,m,l)
                  	 ELSE
                        s(n+1,m+1,l+1,nv) = d%s(n,m,l,nv)
                     END IF
                  END DO
               END DO
            END DO
         END DO


      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      DO nv = 1,10
         CALL interp_int_to_faces(s(:,:,:,nv),d%ncg,d%nc,d%bx,d%by,d%bz)
      END DO
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      DO nv = 1,10
         CALL interp_face_to_edge(s(:,:,:,nv),d%ncg,d%nc,d%bx,d%bz)
      END DO
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      DO nv = 1,10
         CALL interp_edge_to_corner(s(:,:,:,nv),d%ncg,d%nc,d%bx)
      END DO
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END

   SUBROUTINE fill_nu_arrays(d,id,nu,nux,nuy,nuz,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
      use turbulence
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: nu,nux,nuy,nuz
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
         !DO nv = 1,10
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     nu(n+1,m+1,l+1) = d%nu_t(n,m,l)
                     nux(n+1,m+1,l+1) = d%nu_x(n,m,l)
                     nuy(n+1,m+1,l+1) = d%nu_y(n,m,l)
                     nuz(n+1,m+1,l+1) = d%nu_z(n,m,l)
                  END DO
               END DO
            END DO
         !END DO


      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!

      CALL interp_int_to_faces(nu,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(nux,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(nuy,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(nuz,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!

      CALL interp_face_to_edge(nu,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(nux,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(nuy,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(nuz,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      CALL interp_edge_to_corner(nu,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(nux,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(nuy,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(nuz,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END
                                           
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE fill_stat_arrays(d,stats,Q,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
!
      USE domain_definition
      USE stats_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
      TYPE (statfl)                                    :: stats
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,neq+nav) :: Q
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
      IF ( scale )     THEN ! scale out the jacobians for movies
         DO nv = 1,neq+nav
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Q(n+1,m+1,l+1,nv) = stats%Q_av(n,m,l,nv)/nsample
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO nv = 1,neq+nav
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Q(n+1,m+1,l+1,nv) = stats%Q_av(n,m,l,nv)/nsample
                  END DO
               END DO
            END DO
         END DO
      END IF

      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      DO nv = 1,neq+nav
         CALL interp_int_to_faces(Q(:,:,:,nv),d%ncg,d%nc,d%bx,d%by,d%bz)
      END DO
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      DO nv = 1,neq+nav
         CALL interp_face_to_edge(Q(:,:,:,nv),d%ncg,d%nc,d%bx,d%bz)
      END DO
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      DO nv = 1,neq+nav
         CALL interp_edge_to_corner(Q(:,:,:,nv),d%ncg,d%nc,d%bx)
      END DO
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE fill_tensor_arrays(d,Tensor,x,y,z,scale) 
!
!......................................................................
!     date: 11/13/98                                                      
!
!     Fill the plot arrays that will be written out. Includes gauss point
!     values in the interior and interpolated values at the boundaries.
!......................................................................
    USE domain_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      TYPE (domain)                                  :: d
!
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2,9)   :: Tensor
      DOUBLE PRECISION,DIMENSION(nx+2,ny+2,nz+2)     :: x,y,z
      LOGICAL                                        :: scale
!
!     -------------------------------------------
!     store interior point quantities
!     -------------------------------------------
!
        DO id=1,9
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Tensor(n+1,m+1,l+1,id) = d%Tensor(n,m,l,id)
                  END DO
               END DO
            END DO
        END DO
        

      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 1,d%ncg(1)
               x(n+1,m+1,l+1)   = d%xg(1,n,m,l)
               y(n+1,m+1,l+1)   = d%xg(2,n,m,l)
               z(n+1,m+1,l+1)   = d%xg(3,n,m,l)
            END DO
         END DO 
      END DO
!
!     ----------------------------------------------------------------
!     interpolate values to the sides and put into the plotting arrays
!     ----------------------------------------------------------------
!
      do id=1,9
        CALL interp_int_to_faces(Tensor(:,:,:,id),d%ncg,d%nc,d%bx,d%by,d%bz)
      end do
      CALL interp_int_to_faces(x,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(y,d%ncg,d%nc,d%bx,d%by,d%bz)
      CALL interp_int_to_faces(z,d%ncg,d%nc,d%bx,d%by,d%bz)
!
      do id=1,9
        CALL interp_face_to_edge(Tensor(:,:,:,id),d%ncg,d%nc,d%bx,d%bz)
      end do
      CALL interp_face_to_edge(x,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(y,d%ncg,d%nc,d%bx,d%bz)
      CALL interp_face_to_edge(z,d%ncg,d%nc,d%bx,d%bz)
!
      do id=1,9
        CALL interp_edge_to_corner(Tensor(:,:,:,id),d%ncg,d%nc,d%bx)
      end do
      CALL interp_edge_to_corner(x,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(y,d%ncg,d%nc,d%bx)
      CALL interp_edge_to_corner(z,d%ncg,d%nc,d%bx)

!
      RETURN 
      END
 
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE interp_int_to_faces(u,ncg,nc,bx,by,bz)
!
!......................................................................
!     interpolate from the interior to the faces
!......................................................................
!
      USE size
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      DOUBLE PRECISION, DIMENSION(nx+2,ny+2,nz+2) :: u
      DOUBLE PRECISION, DIMENSION(nmax+2)         :: temp
      INTEGER, DIMENSION(3)                       :: ncg,nc
      DOUBLE PRECISION                            :: bx(nx,nx),by(ny,ny),bz(nz,nz)
!
!     ----------------
!     To faces 6 and 4
!     ----------------
!
      DO l = 2,ncg(3)+1
         DO m = 2,ncg(2)+1
            DO n = 1,ncg(1)
               temp(n) = u(n+1,m,l)
            END DO
            CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
            u(1      ,m,l) = u_left
            u(ncg(1)+2,m,l) = u_right
         END DO
      END DO
!
!     ----------------
!     To faces 1 and 2
!     ----------------
!
      DO l = 2,ncg(3)+1
         DO n = 2,ncg(1)+1
            DO m = 1,ncg(2)
               temp(m) = u(n,m+1,l)
            END DO
            CALL interp_to_ends(ny,by,temp,ncg(2),u_left,u_right,nc(2))
            u(n,1      ,l) = u_left
            u(n,ncg(2)+2,l) = u_right
         END DO
      END DO
!
!     ----------------
!     To faces 3 and 5
!     ----------------
!
      DO m = 2,ncg(2)+1
         DO n = 2,ncg(1)+1
            DO l = 1,ncg(3)
               temp(l) = u(n,m,l+1)
            END DO
            CALL interp_to_ends(nz,bz,temp,ncg(3),u_left,u_right,nc(3))
            u(n,m      ,1) = u_left
            u(n,m,ncg(3)+2) = u_right
         END DO
      END DO
!
      RETURN
      END SUBROUTINE interp_int_to_faces
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE interp_face_to_edge(u,ncg,nc,bx,bz)
!
!......................................................................
!     interpolate from the interior of the faces to the edges
!......................................................................
!
      USE size
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      DOUBLE PRECISION, DIMENSION(nx+2,ny+2,nz+2) :: u
      DOUBLE PRECISION, DIMENSION(nmax+2)         :: temp
      INTEGER, DIMENSION(3)                       :: ncg,nc
      DOUBLE PRECISION                            :: bx(nx,nx),bz(nz,nz)
!
!     -----------------------------
!     left/right edges along face 1
!     -----------------------------
!
      DO l = 2,ncg(3)+1
         DO n = 1,ncg(1)
            temp(n) = u(n+1,1,l)
         END DO
         CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
         u(1      ,1,l)  = u_left
         u(ncg(1)+2,1,l) = u_right
      END DO
!     -----------------------------
!     left/right edges along face 2
!     -----------------------------
!
      DO l = 2,ncg(3)+1
         DO n = 1,ncg(1)
            temp(n) = u(n+1,ncg(2)+2,l)
         END DO
         CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
         u(1      ,ncg(2)+2,l)  = u_left
         u(ncg(1)+2,ncg(2)+2,l) = u_right
      END DO
!
!     -----------------------------
!     top/bottom edges along face 1
!     -----------------------------
!
      DO n = 2,ncg(1)+1
         DO l = 1,ncg(3)
            temp(l) = u(n,1,l+1)
         END DO
         CALL interp_to_ends(nz,bz,temp,ncg(3),u_left,u_right,nc(3))
         u(n,1,1      ) = u_left
         u(n,1,ncg(3)+2) = u_right
      END DO
!
!     -----------------------------
!     top/bottom edges along face 2
!     -----------------------------
!
      DO n = 2,ncg(1)+1
         DO l = 1,ncg(3)
            temp(l) = u(n,ncg(2)+2,l+1)
         END DO
         CALL interp_to_ends(nz,bz,temp,ncg(3),u_left,u_right,nc(3))
         u(n,ncg(2)+2,1      ) = u_left
         u(n,ncg(2)+2,ncg(3)+2) = u_right
      END DO
!
!     -----------------------------
!     top/bottom edges along face 6
!     -----------------------------
!
      DO m = 2,ncg(2)+1
         DO l = 1,ncg(3)
            temp(l) = u(1,m,l+1)
         END DO
         CALL interp_to_ends(nz,bz,temp,ncg(3),u_left,u_right,nc(3))
         u(1,m,1      ) = u_left
         u(1,m,ncg(3)+2) = u_right
      END DO
!
!     -----------------------------
!     top/bottom edges along face 4
!     -----------------------------
!
      DO m = 2,ncg(2)+1
         DO l = 1,ncg(3)
            temp(l) = u(ncg(1)+2,m,l+1)
         END DO
         CALL interp_to_ends(nz,bz,temp,ncg(3),u_left,u_right,nc(3))
         u(ncg(1)+2,m,1      ) = u_left
         u(ncg(1)+2,m,ncg(3)+2) = u_right
      END DO
!
      RETURN
      END SUBROUTINE interp_face_to_edge
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE interp_edge_to_corner(u,ncg,nc,bx)
!
!......................................................................
!     interpolate from the interior of the edges to the corners
!......................................................................
!
      USE size
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      DOUBLE PRECISION, DIMENSION(nx+2,ny+2,nz+2) :: u
      DOUBLE PRECISION, DIMENSION(nmax+2)         :: temp
      INTEGER, DIMENSION(3)                       :: ncg,nc
      DOUBLE PRECISION                            :: bx(nx,nx)
!
!     -----------------------------
!     corners 1 & 2
!     -----------------------------
!
      DO n = 1,ncg(1)
         temp(n) = u(n+1,1,1)
      END DO
      CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
      u(1,1,1      ) = u_left
      u(ncg(1)+2,1,1) = u_right
!
!     -----------------------------
!     corners 3 & 4
!     -----------------------------
!
      DO n = 1,ncg(1)
         temp(n) = u(n+1,ncg(2)+2,1)
      END DO
      CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
      u(1,ncg(2)+2,1      ) = u_left
      u(ncg(1)+2,ncg(2)+2,1) = u_right
!
!     -----------------------------
!     corners 5 & 6
!     -----------------------------
!
      DO n = 1,ncg(1)
         temp(n) = u(n+1,1,ncg(3)+2)
      END DO
      CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
      u(1      ,1,ncg(3)+2) = u_left
      u(ncg(1)+2,1,ncg(3)+2) = u_right
!
!     -----------------------------
!     corners 7 & 8
!     -----------------------------
!
      DO n = 1,ncg(1)
         temp(n) = u(n+1,ncg(2)+2,ncg(3)+2)
      END DO
      CALL interp_to_ends(nx,bx,temp,ncg(1),u_left,u_right,nc(1))
      u(1      ,ncg(2)+2,ncg(3)+2) = u_left
      u(ncg(1)+2,ncg(2)+2,ncg(3)+2) = u_right
!
      RETURN
      END SUBROUTINE interp_edge_to_corner
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE interp_to_ends(ncol,b,u,ncg,uleft,uright,nc) 
!
!......................................................................
!     date: 11/10/98                                                     
!     routines called: none                                             
!     includes: none                                                    
!     applicability:all                                                 
!
!     perform 1D interpolation from gauss points to the ends of the interval
!......................................................................
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DOUBLE PRECISION :: b(ncol,*),u(*),uleft,uright
      INTEGER :: ncg,nc
!
      uleft = 0.0d0 
      DO k = 1,ncg 
         uleft = uleft + b(1,k)*u(k) 
      END DO

      uright = 0.0d0 
      DO k = 1,ncg 
         uright = uright + b(nc,k)*u(k) 
      END DO
!
      RETURN 
      END SUBROUTINE interp_to_ends
