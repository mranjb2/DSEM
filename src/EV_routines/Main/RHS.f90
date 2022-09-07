!> \file 
!! This file contains the subroutines related to the EV method residual.

!> @brief This routine calculates the residual (D) based on which the entropy viscosity is calculated.
!!
!! The entropy transport equation for Navier-Stokes system of equations reads
!!
!! \f$ \frac{\partial s}{{\partial t}}  + \frac{\partial }{{\partial {x_j}}}\left( { {u_j} \mathop s } \right) - \frac{ \Phi + \Gamma  }{{T}} -{\Lambda} \geq  0 \f$
!!
!! with \f$ s \f$ denoting the fluid's entropy, and
!!
!! \f$ \Phi=\frac{1}{Re_f}\left( { {S_{ij}}-\frac{1}{3}\delta_{ij}{S_{kk}} } \right)\left( { {S_{ij}}-\frac{1}{3}\delta_{ij}{S_{kk}} } \right) \f$
!!
!! \f$ \Gamma=\frac{1}{Re_f Pr_f(\gamma-1) Ma^2_f}\frac{1}{T}\left( { \frac{\partial T}{\partial x_j} \frac{\partial T}{\partial x_j}} \right) \f$
!!
!! \f$ \Lambda=\frac{1}{Re_f Pr_f(\gamma-1) Ma^2_f}\frac{\partial}{\partial x_j}\left(\frac{1}{T} {  \frac{\partial T}{\partial x_j}} \right) \f$
!!
!! Then, the residual is
!! 
!! \f$ D = \frac{\Phi + \Gamma}{T} \f$
      SUBROUTINE RHS_ent_trans(d,ngrid) 
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE mpi_par
      USE User_Data
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d(ngp) 
!
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxlgg, Qylgg,Qzlgg			  	!< Solution derivatives on lgg points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxglg, Qyglg,Qzglg    			!< Solution derivatives on glg points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxggl, Qyggl,Qzggl,Qzzggl  !< Solution derivatives on ggl points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qx, Qy, Qz            			!< Solution derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxx, Qyy, Qzz         			!< Solution second derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: px, py, pz            			!< Pressure derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Tx, Ty, Tz           			!< Temperature derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Txx, Tyy, Tzz   				  	!< Temperature second derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Plgg, Pglg, Pggl				  	!< Pressure on lobato points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Tlgg, Tglg, Tggl 			  	!< Temperature on lobato points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Txlgg, Txglg, Txggl  			!< Temperature x-gradients on lobato points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Tylgg, Tyglg, Tyggl  			!< Temperature y-gradients on lobato points
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: Tzlgg, Tzglg, Tzggl     		!< Temperature z-gradients on lobato points
      DOUBLE PRECISION, DIMENSION(neq)          :: Q,Qx_loc,Qy_loc,Qz_loc,Qzz_loc
      DOUBLE PRECISION, DIMENSION(3)            :: T1_loc,p1_loc,T2_loc
      DOUBLE PRECISION                          :: RHS, ducros, rhoSensor
!
      DO id = 1,ngrid
!        ------------------------------------------------
!        Compute solution gradients at gauss/gauss points
!        ------------------------------------------------
!
         Plgg (:,:,:) = (gamma-1.d0)*(d(id)%Qlgg(:,:,:,5) - 0.5d0*d(id)%Qlgg(:,:,:,1)*     &
                        ((d(id)%Qlgg(:,:,:,2)/d(id)%Qlgg(:,:,:,1))**2 +                    &
                         (d(id)%Qlgg(:,:,:,3)/d(id)%Qlgg(:,:,:,1))**2 +                    &
                         (d(id)%Qlgg(:,:,:,4)/d(id)%Qlgg(:,:,:,1))**2 ))
         
         Pglg (:,:,:) = (gamma-1.d0)*(d(id)%Qglg(:,:,:,5) - 0.5d0*d(id)%Qglg(:,:,:,1)*     &
                        ((d(id)%Qglg(:,:,:,2)/d(id)%Qglg(:,:,:,1))**2 +                    &
                         (d(id)%Qglg(:,:,:,3)/d(id)%Qglg(:,:,:,1))**2 +                    &
                         (d(id)%Qglg(:,:,:,4)/d(id)%Qglg(:,:,:,1))**2 ))
                         
         Pggl (:,:,:) = (gamma-1.d0)*(d(id)%Qggl(:,:,:,5) - 0.5d0*d(id)%Qggl(:,:,:,1)*     &
                        ((d(id)%Qggl(:,:,:,2)/d(id)%Qggl(:,:,:,1))**2 +                    &
                         (d(id)%Qggl(:,:,:,3)/d(id)%Qggl(:,:,:,1))**2 +                    &
                         (d(id)%Qggl(:,:,:,4)/d(id)%Qggl(:,:,:,1))**2 ))
                         
         Tlgg (:,:,:) = Plgg (:,:,:)*gamma/d(id)%Qlgg(:,:,:,1)
         Tglg (:,:,:) = Pglg (:,:,:)*gamma/d(id)%Qglg(:,:,:,1)
         Tggl (:,:,:) = Pggl (:,:,:)*gamma/d(id)%Qggl(:,:,:,1)
         
         DO nv = 1,neq
            CALL Gauss_Deriv(d(id)%gmetg,d(id)%Qlgg(:,:,:,nv),d(id)%Qglg(:,:,:,nv),       &
                             d(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv), &
			     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
         END DO
!

        CALL Gauss_Deriv(d(id)%gmetg,Tlgg(:,:,:),Tglg(:,:,:),Tggl(:,:,:),       &
                         Tx(:,:,:),Ty(:,:,:),Tz(:,:,:),                   &
			             d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
			             
        CALL Gauss_Deriv(d(id)%gmetg,Plgg(:,:,:),Pglg(:,:,:),Pggl(:,:,:),       &
                         px(:,:,:),py(:,:,:),pz(:,:,:),                   &
			             d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)

!        -----------------------------------------
!        Prolong the derivatives to the cell faces
!        -----------------------------------------
!

      DO nv = 1,neq
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
                CALL interp(nx,d(id)%bx,Qx(:,m,l,nv),d(id)%ncg(1),Qxlgg(:,m,l,nv),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Qy(:,m,l,nv),d(id)%ncg(1),Qylgg(:,m,l,nv),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Qz(:,m,l,nv),d(id)%ncg(1),Qzlgg(:,m,l,nv),d(id)%nc(1))
            END DO
         END DO

         DO l = 1,d(id)%ncg(3)
            DO n = 1,d(id)%ncg(1)
                CALL interp(ny,d(id)%by,Qx(n,:,l,nv),d(id)%ncg(2),Qxglg(n,:,l,nv),d(id)%nc(2))
               CALL interp(ny,d(id)%by,Qy(n,:,l,nv),d(id)%ncg(2),Qyglg(n,:,l,nv),d(id)%nc(2))
              CALL interp(ny,d(id)%by,Qz(n,:,l,nv),d(id)%ncg(2),Qzglg(n,:,l,nv),d(id)%nc(2))
            END DO
         END DO

         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,Qx(n,m,:,nv),d(id)%ncg(3),Qxggl(n,m,:,nv),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Qy(n,m,:,nv),d(id)%ncg(3),Qyggl(n,m,:,nv),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Qz(n,m,:,nv),d(id)%ncg(3),Qzggl(n,m,:,nv),d(id)%nc(3))
            END DO
         END DO

      END DO
      
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
                CALL interp(nx,d(id)%bx,Tx(:,m,l),d(id)%ncg(1),Txlgg(:,m,l),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Ty(:,m,l),d(id)%ncg(1),Tylgg(:,m,l),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Tz(:,m,l),d(id)%ncg(1),Tzlgg(:,m,l),d(id)%nc(1))
            END DO
         END DO

         DO l = 1,d(id)%ncg(3)
            DO n = 1,d(id)%ncg(1)
                CALL interp(ny,d(id)%by,Tx(n,:,l),d(id)%ncg(2),Txglg(n,:,l),d(id)%nc(2))
                CALL interp(ny,d(id)%by,Ty(n,:,l),d(id)%ncg(2),Tyglg(n,:,l),d(id)%nc(2))
                CALL interp(ny,d(id)%by,Tz(n,:,l),d(id)%ncg(2),Tzglg(n,:,l),d(id)%nc(2))
            END DO
         END DO

         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,Tx(n,m,:),d(id)%ncg(3),Txggl(n,m,:),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Ty(n,m,:),d(id)%ncg(3),Tyggl(n,m,:),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Tz(n,m,:),d(id)%ncg(3),Tzggl(n,m,:),d(id)%nc(3))
            END DO
         END DO
         
!        ------------------------------------------------
!        Compute solution gradients at gauss/gauss points
!        ------------------------------------------------
!
!!!!! NOTE THAT GMETG IS USED INSTEAD OF GMET
         DO nv = 1,neq
            CALL Gauss_Deriv(d(id)%gmetg,Qxlgg(:,:,:,nv),Qyglg(:,:,:,nv),               &
                             Qzggl(:,:,:,nv),Qxx(:,:,:,nv),Qyy(:,:,:,nv),Qzz(:,:,:,nv), &
			     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
         END DO
            CALL Gauss_Deriv(d(id)%gmetg,Txlgg(:,:,:),Tyglg(:,:,:),            &
                             Tzggl(:,:,:),Txx(:,:,:),Tyy(:,:,:),Tzz(:,:,:), &
			                 d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
      
!
         DO m = 1,d(id)%ncg(2)
           DO n = 1,d(id)%ncg(1)
             DO l = 1,d(id)%ncg(3)
               Q(1:neq)      = d(id)%Q(n,m,l,1:neq)*d(id)%jacob(n,m,l)
               Qx_loc(1:neq) = Qx(n,m,l,1:neq)
               Qy_loc(1:neq) = Qy(n,m,l,1:neq)
               Qz_loc(1:neq) = Qz(n,m,l,1:neq)
               T1_loc(1)     = Tx(n,m,l)
               T1_loc(2)     = Ty(n,m,l)               
               T1_loc(3)     = Tz(n,m,l)               
               p1_loc(1)     = px(n,m,l)
               p1_loc(2)     = py(n,m,l)               
               p1_loc(3)     = pz(n,m,l)
               T2_loc(1)     = Txx(n,m,l)
               T2_loc(2)     = Tyy(n,m,l)               
               T2_loc(3)     = Tzz(n,m,l)  

               CALL calculate_terms(Q,Qx_loc,Qy_loc,Qz_loc,T1_loc,p1_loc,T2_loc,RHS,ducros,rhoSensor)
               
               d(id)%RHS_ent(n,m,l) = RHS
               d(id)%ducros(n,m,l) = ducros
               d(id)%rhoSensor(n,m,l) = rhoSensor
               
              END DO
            END DO
         END DO		                 
			                 
     END DO
     
      RETURN                                                           
      END SUBROUTINE RHS_ent_trans
!////////////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief This routine calculates the residual (D) based on which the entropy viscosity is calculated. It also calculates
!! the Ducros sensor and the turbulence model density-based sensor
!!
!! Ducros sensor:
!!
!! \f$ \Psi = \frac{(\nabla .u)^2}{(\nabla .u)^2 + |\omega|^2 + \epsilon} \f$
!!
!! Density-based sensor:
!!
!! \f$ \theta_2 = \sqrt{
!! \frac{1}{\dfrac{1}{\rho_x^2 + \epsilon}+\dfrac{1}{\rho_y^2 + \epsilon}+1} +
!! \frac{1}{\dfrac{1}{\rho_x^2 + \epsilon}+\dfrac{1}{\rho_z^2 + \epsilon}+1} +
!! \frac{1}{\dfrac{1}{\rho_y^2 + \epsilon}+\dfrac{1}{\rho_z^2 + \epsilon}+1} 
!! }. \f$
      SUBROUTINE calculate_terms(Q,Q_x,Q_y,Q_z,T1,p1,T2,RHS,ducros,rhoSensor)
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Q,Q_x,Q_y,Q_z
      DOUBLE PRECISION, DIMENSION(3) :: T1,p1,T2
      DOUBLE PRECISION :: RHS					        		!< The EV residual
      DOUBLE PRECISION :: ducros, rhoSensor				!< Ducros sensor and the density-based sensor
      DOUBLE PRECISION :: rhoX2, rhoY2, rhoZ2, eps
!
	  eps = 1.d-9
!
      rho  = Q(1)
      rhou = Q(2)
      rhov = Q(3)
      rhow = Q(4)
      rhoe = Q(5)
!
      u = rhou/rho
      v = rhov/rho
      w = rhow/rho
      p = (gamma-1.d0)*(rhoe - 0.5d0*rho*(u**2 + v**2 + w**2))
      T = p*gamma/rho
      
!
      ux = (Q_x(2) - u*Q_x(1))/rho
      uy = (Q_y(2) - u*Q_y(1))/rho
      uz = (Q_z(2) - u*Q_z(1))/rho
      vx = (Q_x(3) - v*Q_x(1))/rho
      vy = (Q_y(3) - v*Q_y(1))/rho
      vz = (Q_z(3) - v*Q_z(1))/rho
      wx = (Q_x(4) - w*Q_x(1))/rho
      wy = (Q_y(4) - w*Q_y(1))/rho
      wz = (Q_z(4) - w*Q_z(1))/rho
!
      px = p1(1)
      py = p1(2)
      pz = p1(3)
!
      tx = T1(1)
      ty = T1(2)
      tz = T1(3)  
!
      txx = T2(1)
      tyy = T2(2)
      tzz = T2(3)  
!
      phi  = (2.0d0*(ux**2+vy**2+wz**2)+((vx+uy)**2+(uz+wx)**2+(vz+wy)**2)-2.0d0*(ux+vy+wz)**2/3.0d0)
      gam  = (tx**2+ty**2+tz**2)/(T*(gamma-1.0d0)*pr) 
      lam  = ( (((Q_x(1)*p-rho*px)/(gamma*p**2))*tx+txx/t) +          &
               (((Q_y(1)*p-rho*py)/(gamma*p**2))*ty+tyy/t) +          &
               (((Q_z(1)*p-rho*pz)/(gamma*p**2))*tz+tzz/t) ) /((gamma-1.0d0)*pr)
!
      div = ux + vy + wz
      cur = sqrt ( (wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2 )
      !dRho = sqrt(Q_x(1)**2.d0+Q_y(1)**2.d0+Q_z(1)**2.d0)
      rhoX2 = Q_x(1)**2.d0
      rhoY2 = Q_y(1)**2.d0
      rhoZ2 = Q_z(1)**2.d0
!
      RHS = (((phi+gam)/T)/re)
      ducros = max(0.d0,min(1.d0,(div**2/(div**2+cur**2+1.0d-9))))
      !rhoSensor = dRho/(dRho+1.d0)
      rhoSensor = sqrt(														&
      						(1.d0/((1.d0/(rhoX2+eps))+(1.d0/(rhoY2+eps))+1.d0))+	&
      						(1.d0/((1.d0/(rhoX2+eps))+(1.d0/(rhoZ2+eps))+1.d0))+	&
      						(1.d0/((1.d0/(rhoY2+eps))+(1.d0/(rhoZ2+eps))+1.d0))		&
      						)

      RETURN
      END
