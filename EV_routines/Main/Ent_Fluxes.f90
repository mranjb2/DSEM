!>\file Ent_Fluxes.f90
!! This file contains all subroutines related to the entropy fluxes.
!///////////////////////////////////////////////////////////////////////////////////////////////
!                                                                                       ////////
!                                                                                       ////////
!        SUBROUTINE Patch_Fluxes(d,ngrid)                                               ////////
!                                                                                       ////////
!             Calculate and Patch the Entropy Fluxes on Gauss-Lobatto points:           ////////                       
!                                                                                       ////////
!             Contains:                                                                 ////////
!                 SUBROUTINE Int_fluxes_Ent(d(id))                                      ////////
!                 SUBROUTINE mortar_fluxes_ent(mrtr(jl),d,time,itag,a)                  ////////
!                                                                                       ////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief This routine calculates and patches the Entropy Fluxes on Gauss-Lobatto points.
      SUBROUTINE Patch_Fluxes(d,ngrid,mrtr,nmort) 
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
      USE Order_Matrices
      USE input_data
      USE mpi
!
      IMPLICIT NONE 
      TYPE (domain)                          :: d(ngp)     !< domain
      TYPE (mortar)                          :: mrtr(nmp)  !< mortar
!
      INTEGER                                :: ngrid,nmort,id,j,jl,itag,iface,idm,idml,np       
      DOUBLE PRECISION, DIMENSION(nmax,nmax) :: Sbound
!
!     --------------------------------------------------
!     compute interior point fluxes and save edge values
!     on the mortars on different processors
!     --------------------------------------------------        
      DO id = 1,ngrid 
       IF (ngridedge(id) == 1) THEN
         CALL Int_Fluxes_Ent(d(id))
       ENDIF
      END DO

!
!     --------------------------------------------------
!     compute interior point entropy fluxes and save 
!     edge values on the mortars on equal processors
!     --------------------------------------------------
!
      DO id = 1,ngrid 
         IF (ngridedge(id) == 0) THEN
           CALL Int_Fluxes_Ent(d(id))
         ENDIF
      END DO
!
!
!     --------------------------------------------------
!     compute the mortar fluxes on different processors
!     --------------------------------------------------
!
      DO j = 1,nmort
        IF (myid == nmortedge(j,2) .and. nmortedge(j,1) /= 1) THEN
            jl   = nmortmpi(j,1)
            itag = j
            CALL mortar_fluxes_ent(mrtr(jl),d,itag)  
         ENDIF
      END DO
!
!     --------------------------------------------------
!     compute the mortar fluxes on equal processors
!     --------------------------------------------------
!
      DO j = 1,nmort
         IF (myid == nmortedge(j,2) .and. nmortedge(j,1) == 1) THEN
           jl   = nmortmpi(j,1)
           itag = j
           CALL mortar_fluxes_ent(mrtr(jl),d,itag) 
! 
          ENDIF
      END DO
!
!
      RETURN 
      END SUBROUTINE Patch_Fluxes                                          
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!                                                                                       ////////
!                                                                                       ////////
!        SUBROUTINE Int_Fluxes_Ent(d)                                                   ////////
!                                                                                       ////////
!             Calculate the Entropy Fluxes on Gauss-Lobatto points:                     ////////
!                 Project S onto Gauss-Lobatto points: S_lg                             ////////
!                 Then determine points values:        Fs=u*S_lg                        ////////
!                                                                                       ////////
!             Contains:                                                                 ////////
!                 SUBROUTINE Proj_X(d,ngrid)                                            ////////
!                 SUBROUTINE Proj_Y(d,ngrid)                                            ////////
!                 SUBROUTINE Proj_Z(d,ngrid)                                            ////////
!                 SUBROUTINE feflux(Q,S,f,r)                                            ////////
!                 SUBROUTINE geflux(Q,S,g,s)                                            ////////
!                 SUBROUTINE heflux(Q,S,h,t)                                            ////////
!                                                                                       ////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief This routine calculate the Entropy Fluxes on Gauss-Lobatto points:
!!
!! 1. Project \f$S\f$ onto Gauss-Lobatto points: \f$S_{lg}\f$
!!
!! 2. Determine points values: \f$ F_s = u * S_{lg} \f$
      SUBROUTINE Int_Fluxes_Ent(d) 
!
      USE domain_definition
      USE mortar_definition
      USE Material_Properties
      USE mpi_par
      USE Order_Matrices
      USE input_data
      USE mpi
!
      IMPLICIT NONE 
      TYPE (domain) :: d  !< domain
!
      CALL Proj_X(d) 
!
      CALL Proj_Y(d)
!
      CALL Proj_Z(d)
!
      RETURN
      END SUBROUTINE Int_Fluxes_Ent
!      
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the "X" direction fluxes on the lobatto-gauss-gauss points
!! by first interpolating the solutions to the lobatto-gauss-gauss points
!! and then computing the fluxes. Save the boundary and interface
!! solutions along the faces to be later sent to the mortars. This code
!! block is for between faces 6 and 4.
      SUBROUTINE Proj_X(d)
!                                                        
      USE domain_definition
          
      IMPLICIT NONE 
!      
      TYPE (domain)                             :: d
!
      INTEGER                                   :: ngrid,id,m,n,l,nv  
      DOUBLE PRECISION                          :: Shv,fv,rv,gv,sv,hv,tv         
      DOUBLE PRECISION, DIMENSION(nmax)         :: Shin,Shout
      DOUBLE PRECISION, DIMENSION(ngp)          :: dt    
      DOUBLE PRECISION, DIMENSION(neq)          :: Qv 
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Sh 
! 
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
!          
            DO n = 1,d%nc(1)                         
              Shin(n)                = d%Sh(n,m,l)
            END DO
!            
             CALL interp(nx,d%bx,                         &
                         Shin,d%ncg(1),                   &
                         Shout,d%nc(1))
!             
               DO n = 1,d%nc(1)
                 d%Sh_lgg(n,m,l)     = Shout(n)
              END DO
!              
         END DO
      END DO
!      
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
            DO n = 2,d%nc(1)-1
               DO nv = 1,neq
!               
                  Qv(nv)       = d%Qlgg(n,m,l,nv)
                  Shv          = d%Sh_lgg(n,m,l)
!                   
               END DO
!               
               CALL feflux(Qv,Shv,fv,rv)
               CALL geflux(Qv,Shv,gv,sv)
               CALL heflux(Qv,Shv,gv,sv)
! 
               d%fe(n,m,l)  =  d%gmet(1,1,n,m,l)*fv + &
                               d%gmet(1,2,n,m,l)*gv + &
                               d%gmet(1,3,n,m,l)*hv
                               
               d%re(n,m,l)  =  d%gmet(1,1,n,m,l)*rv + &
                               d%gmet(1,2,n,m,l)*sv + &
                               d%gmet(1,3,n,m,l)*tv 
          END DO
         END DO
       END DO
!      
      RETURN                                                                       
      END SUBROUTINE Proj_X
!      
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the "Y" direction fluxes on the lobatto/gauss points
!! by first interpolating the solutions to the lobatto/gauss points
!! and then computing the fluxes. Save the boundary and interface
!! solutions along the faces to be later sent to the mortars. This block
!! is for between faces 1 and 2.
      SUBROUTINE Proj_Y(d)
!                                                       
      USE domain_definition                                                                                           
    
      IMPLICIT NONE 
      TYPE (domain)                             :: d
!
      INTEGER                                   :: ngrid,id,m,n,l,nv      
      DOUBLE PRECISION, DIMENSION(nmax)         :: Shin,Shout
      DOUBLE PRECISION                          :: time,Shv,fv,rv,gv,sv,hv,tv
      DOUBLE PRECISION, DIMENSION(ngp)          :: dt
      DOUBLE PRECISION, DIMENSION(neq)          :: Qv
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Sh   
! 
      
      DO l = 1,d%ncg(3) 
         DO n = 1,d%ncg(1)
!          
            DO m = 1,d%ncg(2)
                Shin(m) = d%Sh(n,m,l)
            END DO
!            
            CALL interp(ny,d%by,                         &
                        Shin,d%ncg(2),                   &
                        Shout,d%nc(2))               
!               
               DO m = 1,d%nc(2)
                   d%Sh_glg(n,m,l)        = Shout(m)
               END DO
!               
         END DO
      END DO
!
      DO l = 1,d%ncg(3)
         DO m = 2,d%nc(2)-1
            DO n = 1,d%ncg(1)
               DO nv = 1,neq
!               
                  Qv(nv)       = d%Qglg(n,m,l,nv)
                  Shv          = d%Sh_glg(n,m,l) 
!                  
               END DO
!               
               CALL feflux(Qv,Shv,fv,rv)
               CALL geflux(Qv,Shv,gv,sv)
               CALL heflux(Qv,Shv,hv,tv)
! 
               d%ge(n,m,l)     = d%gmet(2,1,n,m,l)*fv + &
                                 d%gmet(2,2,n,m,l)*gv + &
                                 d%gmet(2,3,n,m,l)*hv
                                   
               d%se(n,m,l)     = d%gmet(2,1,n,m,l)*rv + &
                                 d%gmet(2,2,n,m,l)*sv + &
                                 d%gmet(2,3,n,m,l)*tv
!                
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE Proj_Y
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the "Z" direction fluxes on the lobatto/gauss points
!! by first interpolating the solutions to the lobatto/gauss points
!! and then computing the fluxes. Save the boundary and interface
!! solutions along the faces to be later sent to the mortars. This code
!! block is for between faces 3 and 5.
      SUBROUTINE Proj_Z(d)
!                                                       
      USE domain_definition                                                                                          
    
     
      IMPLICIT NONE 
      TYPE (domain)                             :: d
!
      INTEGER                                   :: ngrid,id,m,n,l,nv      
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Q
      DOUBLE PRECISION, DIMENSION(nmax)         :: Qin,Shin,Qout,Shout
      DOUBLE PRECISION                          :: time,Shv,fv,rv,gv,sv,hv,tv
      DOUBLE PRECISION, DIMENSION(ngp)          :: dt
      DOUBLE PRECISION, DIMENSION(neq)          :: Qv
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Sh   
!
      DO m = 1,d%ncg(2) 
         DO n = 1,d%ncg(1)
! 
            DO l = 1,d%ncg(3)
                Shin(l)                = d%Sh(n,m,l)
             END DO
!             
             CALL interp(nz,d%bz,                        &
                         Shin,d%ncg(3),                  &
                         Shout,d%nc(3))
!                         
               DO l = 1,d%nc(3)
                   d%Sh_ggl(n,m,l)     = Shout(l)
               END DO
!               
            END DO
          END DO
!
      DO n = 1,d%ncg(1)
         DO m = 1,d%ncg(2)
            DO l = 2,d%nc(3)-1
!            
               DO nv = 1,neq
                  Qv(nv) = d%Qggl(n,m,l,nv)
                  Shv    = d%Sh_ggl(n,m,l) 
               END DO
!               
               CALL feflux(Qv,Shv,fv,rv)
               CALL geflux(Qv,Shv,gv,sv)
               CALL heflux(Qv,Shv,hv,tv)
               
               d%he(n,m,l) = d%gmet(3,1,n,m,l)*fv + &
                             d%gmet(3,2,n,m,l)*gv + &
                             d%gmet(3,3,n,m,l)*hv  
                             
               d%te(n,m,l) = d%gmet(3,1,n,m,l)*rv + &
                             d%gmet(3,2,n,m,l)*sv + &
                             d%gmet(3,3,n,m,l)*tv
!                 
!               
            END DO
         END DO
      END DO
        
!
      RETURN
      END SUBROUTINE Proj_Z
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the "X" direction entropy fluxes on the lobatto/gauss points.
      SUBROUTINE feflux(Q,S,f,r) 
!
      IMPLICIT NONE
!      
      DOUBLE PRECISION, DIMENSION(5) :: Q
      DOUBLE PRECISION               :: S,f,r,u,rho,rhou
!
      rho  = Q(1)
      rhou = Q(2)
!
      u = rhou/rho 
!
      f = u*S
      r=  rhou      
!
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the "Y" direction entropy fluxes on the lobatto/gauss points.
      SUBROUTINE geflux(Q,S,g,se) 
!
      IMPLICIT NONE
!       
      DOUBLE PRECISION, DIMENSION(5) :: Q
      DOUBLE PRECISION               :: S,g,se,v,rho,rhov
!
      rho  = Q(1)
      rhov = Q(3)
!
      v = rhov/rho 
!
      g = S*v 
      se = rhov
!
      RETURN 
      END
!
!//////////////////////////////////////////////////////////////////////////////       
!
!> @brief Compute the "Z" direction entropy fluxes on the lobatto/gauss points.
      SUBROUTINE heflux(Q,S,h,t) 
!
      IMPLICIT NONE
!      
      DOUBLE PRECISION, DIMENSION(5) :: Q
      DOUBLE PRECISION               :: S,h,t,w,rho,rhow
!
      rho  = Q(1)
      rhow = Q(4)
!
      w = rhow/rho 
!
      h = S*w 
      t = rhow
!
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////