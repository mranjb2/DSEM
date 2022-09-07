!> \file part_interp_vis.f90
!! Particle interpolation routines
!////////////////////////////////////////////////////////////////////////////////////////
!/////     SUBROUTINE spec_Pinterp(pgrid,nsurr,msurr,dom,drop,ip)                 ///////
!/////                                                                            ///////
!/////                                                                            ///////
!/////     -----------------------------------------------                        ///////
!/////     interpolate the solution values spectrally to the                      ///////
!/////     particle position (03/15/00)                                           ///////
!/////     -----------------------------------------------                        ///////
!////////////////////////////////////////////////////////////////////////////////////////
!> @brief Spectral particle interpolation
!!
!! Uses spectral basis of p-order equal to element to determine primitive quantites at particle location.
!! This is costly if p-order is high
    SUBROUTINE spec_Pinterp(d,dr,myid,ip)
!
      USE size
      USE physics
      USE domain_definition
      USE particle_definition
        use turbulence

!
      implicit none
      TYPE (domain)     :: d
      TYPE(particle)    :: dr

      DOUBLE PRECISION  :: Qpart(neq)
      double precision  :: nu, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10
      double precision  :: hx, dp
      double precision  :: nu_d(3)
      DOUBLE PRECISION  :: polyn,hy(ny),hz(nz) 
      double precision  :: u,v,w,rho,rhoe,P,T
      LOGICAL           :: Is_NaN
      integer           :: n,m,l, myid, ip

!
! Compute Lagrangian interpolating polynomial in 2 directions
!
      hy=0.0d0
      hz=0.0d0
      DO m=1,d%ncg(2)
         hy(m)=polyn(m,dr%Xpmap(2),d%ncg(2),d%cxg(1:d%ncg(2),2))
      ENDDO
      dr%h(2,:) = hy(:)
      DO l=1,d%ncg(3)
         hz(l)=polyn(l,dr%Xpmap(3),d%ncg(3),d%cxg(1:d%ncg(3),3))
      ENDDO
      dr%h(3,:) = hz(:)
! 
! Interpolate fluid propterties to the particle position with spectral interpolation
!

    sum1 = 0.0d0
    sum2 = 0.0d0
    sum3 = 0.0d0
    sum4 = 0.0d0
    sum5 = 0.0d0
    sum6 = 0.0d0
    sum7 = 0.0d0
    sum8 = 0.0d0
    sum9 = 0.0d0
    sum10 = 0.0d0
    do m=1,d%ncg(1)
        hx = polyn(m,dr%Xpmap(1),d%ncg(1),d%cxg(1:d%ncg(1),1))
        dr%h(1,m) = hx
        do n=1,d%ncg(2)
            do l=1,d%ncg(3)
                sum1 = sum1 + d%Q(m,n,l,1)*hx*hy(n)*hz(l)*d%jacob(m,n,l)
                sum2 = sum2 + d%Q(m,n,l,2)*hx*hy(n)*hz(l)*d%jacob(m,n,l)
                sum3 = sum3 + d%Q(m,n,l,3)*hx*hy(n)*hz(l)*d%jacob(m,n,l)
                sum4 = sum4 + d%Q(m,n,l,4)*hx*hy(n)*hz(l)*d%jacob(m,n,l)
                sum5 = sum5 + d%Q(m,n,l,5)*hx*hy(n)*hz(l)*d%jacob(m,n,l)
                sum6 = sum6 + d%nu_t(m,n,l)*hx*hy(n)*hz(l)
                sum7 = sum7 + d%nu_x(m,n,l)*hx*hy(n)*hz(l)
                sum8 = sum8 + d%nu_y(m,n,l)*hx*hy(n)*hz(l)
                sum9 = sum9 + d%nu_z(m,n,l)*hx*hy(n)*hz(l)
                sum10 = sum10 + d%DP(m,n,l)*hx*hy(n)*hz(l)
            end do
        end do
    end do
    Qpart(1) = sum1
    Qpart(2) = sum2
    Qpart(3) = sum3
    Qpart(4) = sum4
    Qpart(5) = sum5
    nu       = sum6
    nu_d(1)  = sum7
    nu_d(2)  = sum8
    nu_d(3)  = sum9
    dp       = sum10

    !write(*,*)'nu on drop:',nu
    rho = Qpart(1)
    u   = Qpart(2)/Qpart(1)
    v   = Qpart(3)/Qpart(1)
    w   = Qpart(4)/Qpart(1)
    rhoe= Qpart(5)
    P   = (rhoe - rho*0.5d0*(u*u+v*v+w*w))*(gamma-1.0d0)
    T   = P*gamma*mach*mach/rho
    
    dr%Vfp(1)   = u
    dr%Vfp(2)   = v
    dr%Vfp(3)   = w
    dr%P        = P
    dr%Rhofp    = rho
    dr%Tfp      = T
    dr%T        = T
    dr%nu_d     = nu_d
    dr%dp       = dp
    
    


!      dr%Vfp(1)=Qpart(2)/Qpart(1)
!      dr%Vfp(2)=Qpart(3)/Qpart(1)
!      dr%Vfp(3)=Qpart(4)/Qpart(1)
!      dp       =(Qpart(5)-(Qpart(2)**2+Qpart(3)**2+Qpart(4)**2)/(Qpart(1)*2.0d0))*(gamma-1.0d0)
!      dr%P     =dp*d%jacob(1,1,1)!Qpart(1)*dr%Tfp*d%jacob(1,1,1) / 1.40d0
!      dr%Rhofp =Qpart(1)*d%jacob(1,1,1)
!      dr%Tfp   = dr%P * gamma / dr%Rhofp * d%jacob(1,1,1)
!      dr%Tp   = dr%P * gamma * mach * mach / (Qpart(1))/d%jacob(1,1,1)
!      dr%Rhofp =Qpart(1)*d%jacob(1,1,1)
      if (nu.ge.0) then
          dr%nu_t = nu
      else
          dr%nu_t = 0.0d0
      end if
  
       IF ( Is_NaN(dr%Vfp(1)) ) THEN
         write(*,*) dr%Xpmap,myid,ip
         write(*,*) dr%Xp
       ENDIF

      RETURN
    END SUBROUTINE spec_Pinterp

!////////////////////////////////////////////////////////////////////////////////////////      
!/////     SUBROUTINE Lagr_Pinterp(d,dr)                              ///////       
!/////                                                                            ///////       
!/////                                                                            ///////       
!/////     -----------------------------------------------                        ///////       
!/////     interpolate the solution values to the                                 ///////       
!/////     particle position with a sixth order Lagriangian                       ///////       
!/////     polynomial that uses points within a subdomain                         ///////
!/////     -----------------------------------------------                        ///////       
!////////////////////////////////////////////////////////////////////////////////////////
!> @brief Lagrangian interpolant for particles
!!
!! Uses a 6th order Lagrangian polynomial to interpolate quantities to particle location. For high-p cases (over P6)
!! this is significantly faster than the spectral interpolant. However, note that the interpolant is single sided (does not cross element interfaces), and this can cause issues with high gradients near element interfaces.
    SUBROUTINE Lagr_Pinterp(d,dr,myid,ip)         
!
      USE size
      USE constants
      USE domain_definition
      USE particle_definition
      USE physics 
      use turbulence

!
      implicit none
!
      TYPE (domain)    :: d
      TYPE(particle)   :: dr
!
      INTEGER, PARAMETER :: np = 2 ! interpolation order
      DOUBLE PRECISION, DIMENSION(neq) :: Qpart
      DOUBLE PRECISION :: polyn,hy(6),hz(6) 
        double precision :: Xp, Yp, Zp
        double precision :: sum1, sum2, sum3, sum4, sum5, sum6
        double precision :: nu, dp, hx
        double precision :: x,pie, analValS, analValNu, delta_g, percent,y,z
        double precision  :: u,v,w,rho,rhoe,P,T
      INTEGER          :: nsurr, msurr,lsurr
        integer         :: n,m,l, np1, np2, l1, m1, n1
        integer         :: myid, ip
      INTEGER          :: nn(2),mm(2),ll(2)
      LOGICAL          :: Is_NaN
! 
     np1 = np/2-1
     np2 = np/2
!
!    Find the cell the particle is located in
!
      Xp   = ACOS(1.0d0-2.0d0*dr%Xpmap(1))/pi
      Yp   = ACOS(1.0d0-2.0d0*dr%Xpmap(2))/pi
      Zp   = ACOS(1.0d0-2.0d0*dr%Xpmap(3))/pi
!
! Determine first, last and surrounding grid point in x, y and z direction.
!
      nsurr= INT(DBLE(d%ncg(1))*Xp+0.5d0)
!
      IF (nsurr- np1 <= 1)  THEN
         nn(1) = 1 
         nn(2) = np
      ELSEIF (nsurr+np2 >= d%ncg(1)) THEN
         nn(1) = d%ncg(1)-np+1
         nn(2) = d%ncg(1)
      ELSE
         nn(1) = nsurr-np1
         nn(2) = nsurr+np2
      ENDIF
!
      msurr= INT(DBLE(d%ncg(2))*Yp+0.5d0)
!
      IF (msurr-np1 <= 1)  THEN
         mm(1) = 1 
         mm(2) = np
      ELSEIF (msurr+np2 >= d%ncg(2)) THEN
         mm(1) = d%ncg(2)-np+1
         mm(2) = d%ncg(2)
      ELSE
         mm(1) = msurr-np1
         mm(2) = msurr+np2
      ENDIF
!
      lsurr= INT(DBLE(d%ncg(3))*Zp+0.5d0)
!
      IF (lsurr-np1 <= 1)  THEN
         ll(1) = 1 
         ll(2) = np
      ELSEIF (lsurr+np2 >= d%ncg(3)) THEN
         ll(1) = d%ncg(3)-np+1
         ll(2) = d%ncg(3)
      ELSE
         ll(1) = lsurr-np1
         ll(2) = lsurr+np2
      ENDIF
!
! Compute Lagrangian interpolating polynomial in 2 directions
!
!      write(*,*) 'hello4',mm(1),mm(2),msurr,Yp,dr%Xpmap(2),dr%Xp(2)
!      write(*,*) d%corner(2,1),d%corner(2,4)
      hy=0.0d0
      hz=0.0d0
      DO m=1,np
         hy(m)=polyn(m,dr%Xpmap(2),np,d%cxg(mm(1):mm(2),2)) 
      ENDDO
!      write(*,*) 'hello5',ll(1),ll(2),lsurr
      DO l=1,np
         hz(l)=polyn(l,dr%Xpmap(3),np,d%cxg(ll(1):ll(2),3))
      ENDDO
!
! Interpolate to particle position using np^th order 
! Lagrangian interpolation
!
!      write(*,*) 'hello6',nn(1),nn(2),nsurr
!      DO nv=1,neq
!        som  = 0.0d0
!        DO n=1,np
!          hx=polyn(n,dr%Xpmap(1),np,d%cxg(nn(1):nn(2),1))
!          n1 = nn(1)-1+n
!          DO m=1,np
!            m1 = mm(1)-1+m
!            DO l=1,np
!              l1 = ll(1)-1+l
!              som  = som  + d%Q(n1,m1,l1,nv)*hx*hy(m)*hz(l)
!            END DO
!          END DO
!        END DO
!        Qpart(nv)  = som
!      END DO

    sum1 = 0.0d0
    sum2 = 0.0d0
    sum3 = 0.0d0
    sum4 = 0.0d0
    sum5 = 0.0d0
    sum6 = 0.0d0
    do n =1,np
        hx = polyn(n,dr%Xpmap(1),np,d%cxg(nn(1):nn(2),1))
        n1 = nn(1)-1+n
        do m=1,np
            m1=mm(1)-1+m
            do l=1,np
                l1=ll(1)-1+l
                
                sum1 = sum1 + d%Q(n1,m1,l1,1)*hx*hy(m)*hz(l)
                sum2 = sum2 + d%Q(n1,m1,l1,2)*hx*hy(m)*hz(l)
                sum3 = sum3 + d%Q(n1,m1,l1,3)*hx*hy(m)*hz(l)
                sum4 = sum4 + d%Q(n1,m1,l1,4)*hx*hy(m)*hz(l)
                sum5 = sum5 + d%Q(n1,m1,l1,5)*hx*hy(m)*hz(l)
                sum6 = sum6 + nu_t(dr%ngrid,n1,m1,l1)*hx*hy(m)*hz(l)
            end do
        end do
    end do
    Qpart(1) = sum1
    Qpart(2) = sum2
    Qpart(3) = sum3
    Qpart(4) = sum4
    Qpart(5) = sum5
    nu       = sum6
    
!   write(*,*)'nu on drop:',nu
!
!    pie = 3.141592653589790d0
!    x = dr%Xp(1)
!    y = dr%Xp(2)
!    z = dr%Xp(3)
!    analValS = pi * sqrt(0.50d0 * sin(pie*x)*sin(pie*x) + cos(pie*y)*cos(pie*y) + cos(pie*z)*cos(pie*z))
!    
!    write(*,*)'anal val S on drop:',analValS
!    
!    analValNu = (C_s*delta_g(d,1,1,1))*(C_s*delta_g(d,1,1,1))*analValS
!    
!    write(*,*)'anal val nu on drop:',analValNu
!    
!    percent = 100.0d0 * abs(analValNu - nu)/analValNu
!    write(*,*)'Percent Relative Error: ',percent
    



!
! Determine the fluid variables at the droplet position
!
      rho = Qpart(1)
      u   = Qpart(2)/Qpart(1)
      v   = Qpart(3)/Qpart(1)
      w   = Qpart(4)/Qpart(1)
      rhoe= Qpart(5)
      P   = (rhoe - rho*0.5d0*(u*u+v*v+w*w))*(gamma-1.0d0)
      T   = P*gamma*mach*mach/rho

      dr%Vfp(1)   = u
      dr%Vfp(2)   = v
      dr%Vfp(3)   = w
      dr%P        = P
      dr%Rhofp    = rho
      dr%Tfp      = T
      dr%T        = T
    !write(*,*)'dr%Rhofp',dr%Rhofp

       IF ( Is_NaN(dr%Vfp(1)) ) THEN
         write(*,*) dr%Xpmap,myid,ip
         write(*,*) dr%Xp
       ENDIF
 
  
      RETURN
    END SUBROUTINE Lagr_Pinterp

