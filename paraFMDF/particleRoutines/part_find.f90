!> \file part_find.f90 
!! Particle searching algorithms
!///////////////////////////////////////////////////////////////////////////////
!////////                                                               ////////
!////////       part_find.f90                                           ////////
!////////                                                               ////////
!////////       contains:                                               ////////
!////////                                                               ////////
!//////// SUBROUTINE find_subd(dombound,dom,drop,ip)                    ////////
!//////// SUBROUTINE find_particle_cell(ncg,Xpmap,nsurr,msurr,lsurr)    ////////
!////////                                                               ////////
!///////////////////////////////////////////////////////////////////////////////
!> Given the particle and an element, this subroutine determines if the particle is bound by a line of the element
   SUBROUTINE find_subd(p,ngrid,d)
!
!     date: 01/14/00

!     applicability: mortar versions, euler/navier-stokes
!
!     
!
!
      USE size
      USE domain_definition
      USE particle_definition
      USE mpi_par

!
      implicit none
      TYPE (domain)    :: d(ngp)
      TYPE(particle)   :: p
      integer, intent(in) :: ngrid
      
      double precision      :: Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Px, Py
      
      integer               :: id
      
      logical               :: xbound,ybound,zbound,dombound

!
!     ---------------------------------------------------------
!     determine if the particle is bound by the domain and in 
!     which subdomain the particle is found
!     ---------------------------------------------------------
!
!   To determine if the particle is bound by a line of the element side we calculate:
!          [x2 - x1    x3 - x1]
!       det[                  ]
!          [y2 - y1    y3 - y1]
!   The sign will tell us if the particle is left or right
    
    Px = p%Xp(1)
    Py = p%Xp(2)
    
    do id=1,ngrid
        Ax = d(id)%corner(1,1)
        Ay = d(id)%corner(2,1)
        Bx = d(id)%corner(1,2)
        By = d(id)%corner(2,2)
        Cx = d(id)%corner(1,3)
        Cy = d(id)%corner(2,3)
        Dx = d(id)%corner(1,4)
        Dy = d(id)%corner(2,4)
        
        if ( (Bx-Ax)*(Py-Ay)-(By-Ay)*(Px-Ax) .ge. 0.0d0 ) then !we're above the bottom
            if ( (Cx-Bx)*(Py-By)-(Cy-By)*(Px-Bx) .ge. 0.0d0 ) then ! we're to the left of the right
                if ( (Dx-Cx)*(Py-Cy)-(Dy-Cy)*(Px-Cx) .ge. 0.0d0 ) then ! we're below the top
                    if ( (Ax-Dx)*(Py-Dy)-(Ay-Dy)*(Px-Dx) .ge. 0.0d0 ) then ! we're to the right of the left
                        p%ngrid = id
                        goto 569
                    end if
                end if
            end if
        end if
    end do
    !if loop is completed without a hit then the particle isn't in the domain anymore
    p%onoff = 0
    !if the particle hit then it jumps here
569 continue    

   END SUBROUTINE find_subd
!
!///////////////////////////////////////////////////////////////////////////////
!> Determines the lobatto cell in which the particle is located.
      SUBROUTINE find_particle_cell(ncg,Xpmap,nsurr)   
!
!     date: 05/16/03
!
!     applicability: mortar versions, euler/navier-stokes, stac3m
!     find the lobatto cell the particle is located in and returns
!     the number of the Guassian grid point in this cell
!
      USE constants
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      DOUBLE PRECISION  :: Xpmap(3)
      INTEGER           :: ncg(3),nsurr(3)
!
      DO i=1,3
         dum      = DBLE(ncg(i))*ACOS(1.0d0-2.0d0*Xpmap(i))/pi
         nsurr(i) = INT(1.0d0+dum)
      ENDDO
!
      RETURN
   END SUBROUTINE find_particle_cell
   
