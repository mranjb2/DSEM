!> \file part_map.f90
!! Particle mapping routines
!///////////////////////////////////////////////////////////////////////////////
!////////                                                               ////////
!////////       part_map.f90                                            ////////
!////////                                                               ////////
!////////       contains:                                               ////////
!////////                                                               ////////
!//////// SUBROUTINE part_map(d,dbnd,dr)                                ////////
!//////// SUBROUTINE part_unmap(d,dr)                                   ////////
!////////                                                               ////////
!///////////////////////////////////////////////////////////////////////////////
!> Calls appropriate mapping routine for element type
      SUBROUTINE part_map(d,dr) 
!                                                                                
!                                                                             
!     Map the particle onto the subdomain it is in                             
!
!
!     date: 05/04/00
!
      
      USE domain_definition
      USE particle_definition

!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)     :: d
      TYPE (particle)   :: dr
!
      !call map_hex(d,dr)
      call map_pos(d,dr)
      
      RETURN
   END SUBROUTINE part_map

!
!///////////////////////////////////////////////////////////////////////////////
!> @brief Maps particle from physical space to mapped space
!!
!! This subroutine assumes the element has straight sides (rectangular with no skew). Very fast, but cannot handle complex geometry.
      SUBROUTINE map_pos(d,dr) 
!                                                                                
!                                                                             
!     Map the particle onto the subdomain it is in                             
!
!
!     date: 05/04/00
!
      
      USE domain_definition
      USE particle_definition

!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)     :: d
      TYPE (particle)   :: dr
!
      xcor1 = d%corner(1,1) 
      xcor2 = d%corner(1,2) 
      ycor1 = d%corner(2,1) 
      ycor4 = d%corner(2,4) 
      zcor1 = d%corner(3,1) 
      zcor5 = d%corner(3,5)
      
      xp    = dr%Xp(1)
      yp    = dr%Xp(2)
      zp    = dr%Xp(3)

!
!     Determine the X, Y and Z coordinate of the particle in mapped space
!
!     For a quadrilateral cell with straight sides
! 
         X     =  (xp-xcor1)/(xcor2-xcor1)
         Y     =  (yp-ycor1)/(ycor4-ycor1)
         Z     =  (zp-zcor1)/(zcor5-zcor1)
!
        dr%Xpmap(1) = X
        dr%Xpmap(2) = Y
        dr%Xpmap(3) = Z

      RETURN
   END SUBROUTINE map_pos

!> @ brief Maps particle from physical space to mapped space in element
!!
!! Note that this subroutine solves the equation for a skewed element that only has straight z-direction sides. It can handle complex geometry but is much slower that the map_pos routine.
   subroutine map_hex(d,p)
       use domain_definition
       use particle_definition
       implicit none
       type(domain)     :: d
       type(particle)   :: p
       double precision :: a1,a2,a3,a4,b1,b2,b3,b4,aa,bb,cc,m,l
       
       a1   =   d%corner(1,1)
       b1   =   d%corner(2,1)
       a2   =   d%corner(1,2) - a1
       b2   =   d%corner(2,2) - b1
       a3   =   d%corner(1,4) - a1
       b3   =   d%corner(2,4) - b1
       a4   =   a1 - d%corner(1,2) - d%corner(1,4) + d%corner(1,3)
       b4   =   b1 - d%corner(2,2) - d%corner(2,4) + d%corner(2,3)
       
       aa = a4*b3-a3*b4
       bb = a4*b1-a1*b4+a2*b3-a3*b2+b4*p%Xp(1)-a4*p%Xp(2)
       cc = a2*b1-a1*b2+b2*p%Xp(1)-a2*p%Xp(2)
       
       if (aa==0) then
           m = -cc/bb
       else
           !write(*,*) 'skewed elem mapping'
           m = (-bb + dsqrt(bb*bb-4.0*aa*cc))/(2.0*aa)
       end if
       
       l = (p%Xp(1)-a1-a3*m)/(a2+a4*m)
       
       p%Xpmap(1) = l
       p%Xpmap(2) = m
       p%Xpmap(3) = (p%Xp(3)-d%corner(3,1))/(d%corner(3,5)-d%corner(3,1))
       
   end subroutine map_hex
