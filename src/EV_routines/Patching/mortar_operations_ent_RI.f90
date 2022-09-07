!>\file mortar_operations_ent_RI.f90
!! This file contains subroutines related to the entropy viscosity method on mortars.

!> @brief Determine entropy interface flux with Lax-Friedrichs.
      SUBROUTINE mortar_fluxes_ent(mrtr,d,imort)
!
      USE domain_definition
      USE mortar_definition
      USE Material_Properties
      USE mpi_par
!
      IMPLICIT NONE
      TYPE (domain)                              :: d(ngp) 
      TYPE (mortar)                              :: mrtr
      INTEGER                                    :: id,idml,iface,imort
      DOUBLE PRECISION, DIMENSION(nmax,nmax,2)   :: fluxB
!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      CALL Interface_RI(mrtr,fluxB)
!
!     ------------------------------------------
!     send the mortar fluxes to the domain faces
!     ------------------------------------------
!
      CALL send_flux_to_faces_ent(mrtr,d,fluxB,imort)
!
      RETURN
      END SUBROUTINE mortar_fluxes_ent
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Solve the Riemann problem on the mortar.
      SUBROUTINE Interface_RI(mrtr,flux)
!
      USE mortar_definition
      
      IMPLICIT NONE
      
      TYPE (mortar)                              :: mrtr
      INTEGER                                    :: i,j,nv
      DOUBLE PRECISION                           :: rhol,rhor,ul,ur,vl,vr,wl,wr
      DOUBLE PRECISION                           :: Sl,Sr,fl,fr,rl,rr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,2)   :: flux
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql,Qr
      DOUBLE PRECISION, DIMENSION(2)             :: f

!     ---------------------------------------
!     Solve the Riemann problem on the mortar
!     ---------------------------------------
!

      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)
!
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
!               
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,2)
                     Qr(nv) = mrtr%Q(i,j,nv,1)
                  END DO                  
!                     
                     CALL Riemann_Ent(Ql,Qr,mrtr%n_hat(:,i,j),mrtr%norm(i,j),f)                 
!                  
                     flux(i,j,1) = f(1)
                     flux(i,j,2) = f(2)
!                     
               END DO
            END DO   
!
           CASE (2,4,5)
!
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
!               
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO
!                   
                     CALL Riemann_Ent(Ql,Qr,mrtr%n_hat(:,i,j),mrtr%norm(i,j),f)          
!                 
                     flux(i,j,1) = f(1)
                     flux(i,j,2) = f(2)
!                     
               END DO
            END DO 
!

      END SELECT
!
      RETURN
      END SUBROUTINE Interface_RI
!
!/////////////////////////////////////////////////////////////////////////////////////////
!                                                                                 ////////
!                                                                                 ////////
!       Riemann_Ent  (Ql,Qr,n_hat,ds,flux)                                        ////////
!                                                                                 ////////
!         Solve the Riemann problem on the mortar to find the Riemann solution    ////////
!         in order to calculate the entropy on the mortar.                        ////////
!         Contains:                                                               ////////
!            feflux(Q,S,fs,rs)                                                    ////////
!            geflux(Q,S,gs,ss)                                                    //////// 
!            heflux(Q,S,hs,ts)                                                    ////////
!                                                                                 ////////
!/////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief Solve the Riemann problem on the mortar to find the Riemann solution
!! in order to calculate the entropy on the mortar.
      SUBROUTINE Riemann_Ent(Ql,Qr,n_hat,ds,flux)                                         
!
      USE physics
!
      IMPLICIT NONE 
      DOUBLE PRECISION, DIMENSION(5) :: Ql,Qr,Q
      DOUBLE PRECISION, DIMENSION(2) :: flux
      DOUBLE PRECISION               :: fs,gs,hs,rs,ss,ts
      DOUBLE PRECISION               :: en1,en2,en3
      DOUBLE PRECISION               :: rho,rhou,rhov,rhow,rhoe
      DOUBLE PRECISION               :: rhon,rhoun,rhovn,rhown,rhoen
      DOUBLE PRECISION               :: qdotn,qdotnn,p,a,pn,an,e0,eplus,eminus
      DOUBLE PRECISION               :: rplus,rminus,u,v,w,S,etos
      DOUBLE PRECISION               :: q_tanx,q_tany,q_tanz,ds
      DOUBLE PRECISION, DIMENSION(3) :: n_hat
!
      en1 = n_hat(1)
      en2 = n_hat(2)
      en3 = n_hat(3)

      rho  = Ql(1)
      rhou = Ql(2)
      rhov = Ql(3)
      rhow = Ql(4)
      rhoe = Ql(5)

      rhon  = Qr(1)
      rhoun = Qr(2)
      rhovn = Qr(3)
      rhown = Qr(4)
      rhoen = Qr(5)

      qdotn  = (en1*rhou + en2*rhov + en3*rhow)/rho 
      qdotnn = (en1*rhoun + en2*rhovn + en3*rhown)/rhon

      p = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2 + rhow**2)) 
      a = sqrt(gamma*p/rho) 
      pn = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2 + rhown**2)) 
      an = sqrt(gamma*pn/rhon)

      e0     = 0.5d0*(qdotn + qdotnn) 
      eplus  = e0 + 0.5d0*(a + an) 
      eminus = e0 - 0.5d0*(a + an) 
!
      IF(e0 .gt. 0.0d0)     THEN 
         etos = p/rho**gamma
         q_tanx = rhou/rho - qdotn*en1
         q_tany = rhov/rho - qdotn*en2
         q_tanz = rhow/rho - qdotn*en3
      else 
         etos = pn/rhon**gamma 
         q_tanx = rhoun/rhon - qdotnn*en1
         q_tany = rhovn/rhon - qdotnn*en2
         q_tanz = rhown/rhon - qdotnn*en3
      endif 
!
      if(eplus .ge. 0.0d0)     then 
         rplus = qdotn + 2.d0/(gamma-1.d0)*a 
      else 
         rplus = qdotnn + 2.d0/(gamma-1.d0)*an 
      endif 
!
      if(eminus .ge. 0.0d0)     then 
         rminus = qdotn - 2.d0/(gamma-1.d0)*a 
      else 
         rminus = qdotnn - 2.d0/(gamma-1.d0)*an 
      endif 
!
      qdotn = 0.5d0*(rplus + rminus) 
      a = 0.25d0*(gamma-1.d0)*(rplus - rminus)

      rho = (a*a/(gamma*etos))**(1./(gamma-1.d0)) 
      p = rho**gamma*etos
      u = q_tanx + en1*qdotn
      v = q_tany + en2*qdotn
      w = q_tanz + en3*qdotn
!
      rhou = rho*u 
      rhov = rho*v
      rhow = rho*w
      rhoe = p/(gamma-1.d0) + 0.5d0*rho*(u**2 + v**2 + w**2)

      Q(1) = rho
      Q(2) = rhou
      Q(3) = rhov
      Q(4) = rhow
      Q(5) = rhoe
!
      S = (rho/(gamma**2-gamma))*log(p/rho**gamma)
!      
      call feflux(Q,S,fs,rs)
      call geflux(Q,S,gs,ss) 
      call heflux(Q,S,hs,ts) 
!
      flux(1) = ds*(en1*fs + en2*gs + en3*hs)
      flux(2) = ds*(en1*rs + en2*ss + en3*ts)
!
      return 
      END
!
!///////////////////////////////////////////////////////////////////////