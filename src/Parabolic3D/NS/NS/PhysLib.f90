!> \file PhysLib.f90
!! Physics routines for the euler-gas dynamics equations
!
!     Physics routines for the euler gas-dynamics equations
!
!     The variable mappings are
!
!              Q(1) = rho
!              Q(2) = rhou
!              Q(3) = rhov
!              Q(4) = rhow
!              Q(5) = rhoe
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the X inviscid flux.
      SUBROUTINE fflux(Q,f) 
!
!......................................................................
!     date: 11/9/98          
!     routines called: none  
!     applicability: Euler equations     
!
!     compute the X inviscid flux
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Q,f
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
!
      f(1) = rhou 
      f(2) = p + rhou*u 
      f(3) = rhou*v 
      f(4) = rhou*w 
      f(5) = u*(rhoe + p) 
!
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the Y inviscid flux.
      SUBROUTINE gflux(Q,g) 
!
!......................................................................
!     date: 11/9/98          
!     routines called: none  
!     applicability: Euler equations     
!
!     compute the Y inviscid flux
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DOUBLE PRECISION, DIMENSION(5) :: Q,g
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
!
      g(1) = rhov 
      g(2) = rhou*v 
      g(3) = p + rhov*v 
      g(4) = rhow*v 
      g(5) = v*(rhoe + p) 
!
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the Z inviscid flux.
      SUBROUTINE hflux(Q,h)
!
!......................................................................
!     date: 11/9/98          
!     routines called: none  
!     applicability: all     
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DOUBLE PRECISION, DIMENSION(5) :: Q,h
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
!
      h(1) = rhow 
      h(2) = rhou*w 
      h(3) = rhov*w
      h(4) = p + rhow*w 
      h(5) = w*(rhoe + p) 
!
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief use normal riemann invariants to solve the riemann problem.
      SUBROUTINE riemann_alt(n_hat,ds,Ql,Qr,flux)                                         
!
!......................................................................
!     date: 11/9/98
!     routines called: fflux, gflux, hflux
!
!     use normal riemann invariants to solve the riemann problem  
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DOUBLE PRECISION, DIMENSION(5) :: Ql,Qr,f,g,Q,flux,h
      DOUBLE PRECISION               :: q_tanx,q_tany,q_tanz
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
      IF(e0 > 0.0d0)     THEN 
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
      if(eplus >= 0.0d0)     then 
         rplus = qdotn + 2.d0/(gamma-1.d0)*a 
      else 
         rplus = qdotnn + 2.d0/(gamma-1.d0)*an 
      endif 
!
      if(eminus >= 0.0d0)     then 
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
      call fflux(Q,f)
      call gflux(Q,g) 
      call hflux(Q,h) 
!
      DO k = 1,5
         flux(k) = ds*(en1*f(k) + en2*g(k) + en3*h(k))
      END DO
!
      return 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the Roe flux with entropy fix. 
      SUBROUTINE riemann(nhat,ds,Qleft,Qright,flux)                                 
!
!......................................................................
!     date: 12/15/99 from the 2d code
!     routines called: none  
!     applicability: all     
!
!     compute the Roe flux with entropy fix.                            
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Qleft,Qright,flux
      DOUBLE PRECISION, DIMENSION(3) :: nhat
!
      rho  = Qleft(1)
      rhou = Qleft(2)
      rhov = Qleft(3)
      rhow = Qleft(4)
      rhoe = Qleft(5)

      rhon  = Qright(1)
      rhoun = Qright(2)
      rhovn = Qright(3)
      rhown = Qright(4)
      rhoen = Qright(5)

      ul = rhou/rho 
      vl = rhov/rho 
      wl = rhow/rho 
      pleft = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2 + rhow**2 )) 
!
      ur = rhoun/rhon 
      vr = rhovn/rhon 
      wr = rhown/rhon 
      pright = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2+ rhown**2)) 
!
      ql = nhat(1)*ul + nhat(2)*vl + nhat(3)*wl
      qr = nhat(1)*ur + nhat(2)*vr + nhat(3)*wr
      hl = 0.5d0*(ul*ul + vl*vl + wl*wl) + gamma/(gamma-1.d0)*pleft/rho 
      hr = 0.5d0*(ur*ur + vr*vr + wr*wr) + gamma/(gamma-1.d0)*pright/rhon 
!
!     ---------------------
!     square root averaging  
!     ---------------------
!
      rtd = sqrt(rho*rhon) 
      betal = rho/(rho + rtd) 
      betar = 1.d0 - betal 
      utd = betal*ul + betar*ur 
      vtd = betal*vl + betar*vr 
      wtd = betal*wl + betar*wr 
      htd = betal*hl + betar*hr 
      atd2 = (gamma-1.d0)*(htd - 0.5d0*(utd*utd + vtd*vtd + wtd*wtd)) 
      atd = sqrt(atd2) 
      qtd = utd*nhat(1) + vtd*nhat(2)  + wtd*nhat(3)
!
      IF(qtd >= 0.0d0)     THEN

         dw1 = 0.5d0*((pright - pleft)/atd2 - (qr - ql)*rtd/atd) 
         sp1 = qtd - atd 
         sp1m = min(sp1,0.0d0) 
         hd1m = ((gamma+1.d0)/4.d0*atd/rtd)*dw1 
         eta1 = max(-abs(sp1) - hd1m,0.0d0) 
         udw1 = dw1*(sp1m - 0.5d0*eta1) 
         rql = rho*ql 
         flux(1) = ds*(rql + udw1) 
         flux(2) = ds*(rql*ul + pleft*nhat(1) + udw1*(utd - atd*nhat(1))) 
         flux(3) = ds*(rql*vl + pleft*nhat(2) + udw1*(vtd - atd*nhat(2))) 
         flux(4) = ds*(rql*wl + pleft*nhat(3) + udw1*(wtd - atd*nhat(3))) 
         flux(5) = ds*(rql*hl + udw1*(htd - qtd*atd)) 

      ELSE 

         dw4 = 0.5d0*((pright - pleft)/atd2 + (qr - ql)*rtd/atd) 
         sp4 = qtd + atd 
         sp4p = max(sp4,0.0d0) 
         hd4 = ((gamma+1.d0)/4.d0*atd/rtd)*dw4 
         eta4 = max(-abs(sp4) + hd4,0.0d0) 
         udw4 = dw4*(sp4p + 0.5d0*eta4) 
         rqr = rhon*qr 
         flux(1) = ds*(rqr - udw4) 
         flux(2) = ds*(rqr*ur + pright*nhat(1) - udw4*(utd + atd*nhat(1))) 
         flux(3) = ds*(rqr*vr + pright*nhat(2) - udw4*(vtd + atd*nhat(2))) 
         flux(4) = ds*(rqr*wr + pright*nhat(3) - udw4*(wtd + atd*nhat(3))) 
         flux(5) = ds*(rqr*hr - udw4*(htd + qtd*atd)) 
      ENDIF
      RETURN 
      END SUBROUTINE riemann

!
!///////////////////////////////////////////////////////////////////////
!
!> @brief use normal riemann invariants to solve the riemann problem.
      SUBROUTINE riemann_solution(n_hat,Ql,Qr,Q)                                         
!
!......................................................................
!     date: 03/12/01
!
!     use normal riemann invariants to solve the riemann problem  
!......................................................................
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DOUBLE PRECISION, DIMENSION(5) :: Ql,Qr,Q
      DOUBLE PRECISION               :: q_tanx,q_tany,q_tanz
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
      IF(e0 > 0.0d0)     THEN 
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
      if(eplus >= 0.0d0)     then 
         rplus = qdotn + 2.d0/(gamma-1.d0)*a 
      else 
         rplus = qdotnn + 2.d0/(gamma-1.d0)*an 
      endif 
!
      if(eminus >= 0.0d0)     then 
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
      return 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the viscous flux in x-direction.
      SUBROUTINE x_vis_flux(Q,Q_x,Q_y,Q_z,f,amu)
!
!     date: 03/14/01
!     routines called: diffmol
!     applicability: Navier-Stokes equations
!
!     compute the viscous flux in x-direction
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Q,Q_x,Q_y,Q_z,f
!     EXTERNAL diffmol
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

      px = (gamma-1.d0)*(Q_x(5) - 0.5d0*Q_x(1)*(u*u + v*v + w*w) - rho*(u*ux + v*vx + w*wx))
      tx = (px - p*Q_x(1)/rho)*mach*mach*gamma/rho
!
      amu    = amu
      akappa = 1.0d0 + (amu-1.0d0)*2.0d-1/(gamma-1.0d0)
      div = ux + vy + wz
      strainxy = uy + vx
      strainxz = wx + uz
!
      f(1) = 0.0d0
      f(2) = 2.d0*amu*(ux - div/3.d0)
      f(3) = amu*strainxy
      f(4) = amu*strainxz
      f(5) = u*f(2) + v*f(3) + w*f(4) + akappa*tx/(pr*mach*mach*(gamma-1.0d0))
!
      RETURN
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the viscous flux in y-direction.
      SUBROUTINE y_vis_flux(Q,Q_x,Q_y,Q_z,g,amu)
!
!     date: 03/14/01
!     routines called: diffmol
!     applicability: Navier-Stokes equations
!
!     compute the viscous flux in y-direction
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Q,Q_x,Q_y,Q_z,g
!     EXTERNAL diffmol
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

      py = (gamma-1.d0)*(Q_y(5) - 0.5d0*Q_y(1)*(u*u + v*v + w*w) - rho*(u*uy + v*vy + w*wy))
      ty = (py - p*Q_y(1)/rho)*mach*mach*gamma/rho
!
      amu      =  amu
      akappa = 1.0d0 + (amu-1.0d0)*2.0d-1/(gamma-1.0d0)
      div      = ux + vy + wz
      strainxy = uy + vx
      strainyz = wy + vz
!
      g(1) = 0.0d0
      g(2) = amu*strainxy
      g(3) = 2.d0*amu*(vy - div/3.d0)
      g(4) = amu*strainyz
      g(5) = v*g(3) + u*g(2) + w*g(4) + akappa*ty/(pr*mach*mach*(gamma-1.0d0))
!
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
!
!> @brief compute the viscous flux in z-direction.
      SUBROUTINE z_vis_flux(Q,Q_x,Q_y,Q_z,h,amu)
!
!     date: 03/14/01
!     routines called: diffmol
!     applicability: Navier-Stokes equations
!
!     compute the viscous flux in z-direction
!
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Q,Q_x,Q_y,Q_z,h
!     EXTERNAL diffmol
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

      pz = (gamma-1.d0)*(Q_z(5) - 0.5d0*Q_z(1)*(u*u + v*v + w*w) - rho*(u*uz + v*vz + w*wz))
      tz = (pz - p*Q_z(1)/rho)*mach*mach*gamma/rho
!
      amu      = amu
      akappa = 1.0d0 + (amu-1.0d0)*2.0d-1/(gamma-1.0d0)
      div      = ux + vy + wz
      strainxz = wx + uz
      strainyz = wy + vz
!
      h(1) = 0.0d0
      h(2) = amu*strainxz
      h(3) = amu*strainyz
      h(4) = 2.d0*amu*(wz - div/3.d0)
      h(5) = v*h(3) + u*h(2) + w*h(4) + akappa*tz/(pr*mach*mach*(gamma-1.0d0))
!
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
!
!> @brief use normal riemann invariants to solve the riemann problem.
      SUBROUTINE riemann_alt_periods(n_hat,ds,Ql,Qr,flux)
!
!......................................................................
!     date: 11/9/98
!     routines called: fflux, gflux, hflux
!
!     use normal riemann invariants to solve the riemann problem
!......................................................................
!
      USE physics
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Ql,Qr,f,g,Q,flux,h
      DOUBLE PRECISION               :: q_tanx,q_tany,q_tanz
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
      IF(e0 > 0.0d0)     THEN
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
      if(eplus >= 0.0d0)     then
         rplus = qdotn + 2.d0/(gamma-1.d0)*a
      else
         rplus = qdotnn + 2.d0/(gamma-1.d0)*an
      endif
!
      if(eminus >= 0.0d0)     then
         rminus = qdotn - 2.d0/(gamma-1.d0)*a
      else
         rminus = qdotnn - 2.d0/(gamma-1.d0)*an
      endif
!
      qdotn = 0.5d0*(rplus + rminus)
      a = 0.25d0*(gamma-1.d0)*(rplus - rminus)
     rho = (a*a/(gamma*etos))**(1./(gamma-1.d0))
      p = rho**gamma*etos
      pold = p
      dpdx = -3.d0/re
      dpdx = dpdxturb
      p = p + dpdx*char_length
      rho = p*rho/pold
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
      call fflux(Q,f)
      call gflux(Q,g)
      call hflux(Q,h)
!
      DO k = 1,5
         flux(k) = ds*(en1*f(k) + en2*g(k) + en3*h(k))
      END DO
!
      return
      END
!
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief use normal riemann invariants to solve the riemann problem.
      SUBROUTINE osher(n_hat,ds,Ql,Qr,flux)
!
!......................................................................
!     date: 11/9/98
!     routines called: fflux, gflux, hflux
!
!     use normal riemann invariants to solve the riemann problem
!......................................................................
!
      USE physics
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(5) :: Ql,Qr,f,g,Q,flux,h
      DOUBLE PRECISION               :: q_tanx,q_tany,q_tanz
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


      p  = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2 + rhow**2))
      a  = sqrt(gamma*p/rho)
      pn = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2 + rhown**2))
      an = sqrt(gamma*pn/rhon)
!
      q_tanx = rhoun/rhon - qdotnn*en1
      q_tany = rhovn/rhon - qdotnn*en2
      q_tanz = rhown/rhon - qdotnn*en3
!
      z    = (gamma-1.0d0)/(2.0d0*gamma)
      Ha   = -qdotnn/an+2.0d0/(gamma-1.0d0)
      Ha   = Ha/(qdotn/a+2.0d0/(gamma-1.0d0))
      pn   = p/Ha**(1/z)
!   compute p*
      pstar = a+an-(qdotnn-qdotn)*(gamma-1.0d0)/2.0d0
      pstar = pstar/(a/p**z+an/pn**z)
      pstar = pstar**(1.0d0/z)
!   compute u*
      ustar = Ha*qdotn/a+qdotnn/an+2.0d0*(Ha-1)/(gamma-1.0d0)
      ustar = ustar/(Ha/a+1.0d0/an)
!  compute rho23
      rho23=pstar*gamma/an**2
!
      rho  = rho23
      p    = pstar
      qdotn= ustar

!
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
      call fflux(Q,f)
      call gflux(Q,g)
      call hflux(Q,h)
!
      DO k = 1,5
         flux(k) = ds*(en1*f(k) + en2*g(k) + en3*h(k))
      END DO
!
      return

      END
!

