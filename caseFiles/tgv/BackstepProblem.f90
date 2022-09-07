!
!///////////////////////////////////////////////////////////////////////
!
      MODULE user_keywords
!
!......................................................................
!     store user modifiable keywords needed to read the input file.
!     any number of user keywords can be added to be accessible from
!     the user modifiable routines:
!
!        user_setup
!        do_user_stuff
!        cleanup
!        initc
!        + other user supplied routines called by them
!
!     communication is through the data defined in the module "User_Data"
!     The values in User_Data are set below in the routine "set_user_values"
!......................................................................
!
      SAVE
!
      INTEGER                           :: num_user_keywords = 3
      CHARACTER (LEN=10), DIMENSION(3)  :: user_keyword_list = (/ &
                                                'theta     '    , &
                                                'phi       '    ,  &
                                                'outlet    '      &
                                           /)
      LOGICAL, DIMENSION(3) :: user_keyword_set = .false.
      END MODULE user_keywords
!
!
!///////////////////////////////////////////////////////////////////////
!
      MODULE User_Data
         USE size
         SAVE

         DOUBLE PRECISION          :: pout         ! outlet pressure
         DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t        ! viscous fluxes
         DOUBLE PRECISION :: xxin(38),uuin(38)
         INTEGER          :: itimeinflowbc 

      END MODULE User_Data
!
!///////////////////////////////////////////////////////////////////////
!
      MODULE Skin_Data
         USE size
         SAVE
         DOUBLE PRECISION          :: bmatu1b(nmax,nmax,nmax)      ! integration matrix
         DOUBLE PRECISION          :: bmataa(nmax-1,nmax-1)      ! integration matrix

      END MODULE Skin_Data
!
!///////////////////////////////////////////////////////////////////////
module particleInlets
    use size
    integer                          :: inletCount
    integer, allocatable             :: inletList(:)
    double precision, allocatable    :: cumPart(:)
    
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE set_user_values(the_word,input_line,kword)
!
!......................................................................
!     read in the input data-file and save the values in User_Data
!......................................................................
!
      USE physics
      USE User_Data
      USE user_keywords
      USE constants
!
      DOUBLE PRECISION      :: get_dp_value   ! extracts a double prec number
      INTEGER               :: get_int_value  ! extracts an integer
      CHARACTER (LEN = 132) :: input_line
      CHARACTER (LEN = 10)  :: the_word
!
!     ---------------------------------
!     two user keywords for this problem
!     ---------------------------------
!
      SELECT CASE(the_word)
         CASE('theta')
            theta = get_dp_value(input_line)
            theta = theta*180.d0/pi
            user_keyword_set(kword) = .true.
         CASE('phi')
            phi = get_dp_value(input_line)
            phi = phi*180.d0/pi
            user_keyword_set(kword) = .true.
         CASE('outlet')
            pout = get_dp_value(input_line)
            user_keyword_set(kword) = .true.


      END SELECT
!
      RETURN
      END SUBROUTINE SET_USER_VALUES
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE user_setup(d,ngrid,mrtr,nmort,restart)
!
!......................................................................
!     ALLOWS PROBLEM SPECIFIC COMPUTATIONS TO BE DONE
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE Skin_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)

!
      DOUBLE PRECISION                  :: xwgt(nmax),cxg(nmax)
      INTEGER:: o
      INTEGER,DIMENSION(ngp)            :: npoin
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
      DOUBLE PRECISION                  :: polyn
      INTEGER, PARAMETER                :: legendre = 1
      INTEGER,DIMENSION(3300)           :: mortmat
      LOGICAL                           :: nodouble,restart
      integer                           :: yperiodic, xperiodic, zperiodic


      
IF (.not. restart) THEN
      mortmat = 0
      chick = 1
      IF (chick ==1) THEN
!
!  set periodic boundary conditions
!

!  find the interfaces  at which periodic boundary conditions
!  apply, adjust domain and mortar parameters
      zperiodic = 1
      if(zperiodic == 1) then
    !write(*,*)'hi jon line 142'
            kkk = 0
            do idm = 1,ngrid
                  if(d(idm)%bcond(3) == 'outflow') then
                        d(idm)%bcond(3)  = 'interface'
                        d(idm)%ibtype(3) = 1
                        mrtr_no          = d(idm)%mortar(1,3)
                        d(idm)%which_surface(3) = 0
                        do ids=1,ngrid
                              if (d(ids)%bcond(5)    =='outflow'            .and.  &
                                  d(idm)%corner(1,1) == d(ids)%corner(1,5)  .and.  &
                                  d(idm)%corner(2,1) == d(ids)%corner(2,5))  then
                        !write(*,*)'hello jon line 154'
                        kkk          = kkk + 1
                        mrtr_no_sl   = d(ids)%mortar(1,5)
                        mortmat(kkk) = mrtr_no_sl 
                        
                        mrtr(mrtr_no)%len(:,2)  = mrtr(mrtr_no_sl)%len(:,1)
                        mrtr(mrtr_no)%id(2)     = ids
                        mrtr(mrtr_no)%iface(2)  = 5
                        d(ids)%bcond(5)         = 'interface'
                        d(ids)%ibtype(5)        = 1
                        d(ids)%which_surface(5) = 0
                        d(ids)%mortar(1,5)      = mrtr_no
                        d(ids)%mortar(2,5)      = 2
                    end if
                end do
            end if
        end do
    end if


      yperiodic = 1
      if(yperiodic == 1) then
      !This is periodic in y
      !kkk = 0 !a counter
      DO idm = 1,ngrid
         IF (d(idm)%bcond(1) == 'outflow') THEN  !if inflow on face 1
            d(idm)%bcond(1)  = 'interface'      !Change to interface
            d(idm)%ibtype(1) = 1
            mrtr_no          = d(idm)%mortar(1,1) !Finds the corresponding mortar
            d(idm)%which_surface(1) = 0 
            DO ids = 1,ngrid                      !Over all the elements
               IF (d(ids)%bcond(2) == 'outflow'              .and.  & !is the face on the other side of the element inflow on face 2?
                   d(idm)%corner(1,1) == d(ids)%corner(1,4) .and.  & !do corner points match?
                   d(idm)%corner(3,1) == d(ids)%corner(3,4)) THEN    !matches the x & z -corner 1 is x and corner 3 is z
                                                                     !if 1 and 4 match then same element

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,2)
                   mortmat(kkk) = mrtr_no_sl
                  
                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids 
                   mrtr(mrtr_no)%iface(2) = 2
                   d(ids)%bcond(2)  = 'interface'
                   d(ids)%ibtype(2) = 1
                   d(ids)%which_surface(2) = 0
                   d(ids)%mortar(1,2) = mrtr_no
                   d(ids)%mortar(2,2) = 2
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      end if
!
!  find the interfaces  at which periodic boundary conditions
!  apply, adjust domain and mortar parameters
!
    chicka=1 !periodic in x
    IF (chicka==1) THEN
      DO idm = 1,ngrid
         IF (d(idm)%bcond(6) == 'outflow') THEN
             d(idm)%bcond(6)  = 'interface'
             d(idm)%ibtype(6) = 1
             mrtr_no          = d(idm)%mortar(1,6)
             d(idm)%which_surface(6) = 0

            DO ids = 1,ngrid
               IF (d(ids)%bcond(4) == 'outflow'              .and. &
                   d(idm)%corner(3,1) == d(ids)%corner(3,2) .and.  &
                   d(idm)%corner(2,1) == d(ids)%corner(2,2)) THEN

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,4)
                   mortmat(kkk) = mrtr_no

                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids
                   mrtr(mrtr_no)%iface(2) = 4
                   mrtr(mrtr_no_sl)%len(:,2) = mrtr(mrtr_no)%len(:,1)
                   mrtr(mrtr_no_sl)%id(2)    = idm
                   mrtr(mrtr_no_sl)%iface(2) = 6

                   d(ids)%bcond(4)  = 'interface'
                   d(ids)%ibtype(4) = 1
                   d(ids)%which_surface(4) = 0
                   d(ids)%mortar(1,4) = mrtr_no
                   d(ids)%mortar(2,4) = 1
              ENDIF
            ENDDO
         ENDIF
      ENDDO
    ENDIF


      nmortdouble = kkk
      ENDIF
!
!  Clean up mortars
!
      nmorta = 0
      DO im = 1,nmort
        nodouble = .true.
        DO imd = 1,nmortdouble
           IF (im==mortmat(imd)) THEN
             nodouble = .false.
           ENDIF
        ENDDO
        IF (nodouble) THEN 
          nmorta = nmorta +1
          mrtr(nmorta) = mrtr(im)
          idm          = mrtr(im)%id(1) 
          ids          = mrtr(im)%id(2)
          ifacem       = mrtr(im)%iface(1) 
          ifaces       = mrtr(im)%iface(2) 
          d(idm)%mortar(1,ifacem) = nmorta 
          IF (ids /= 0) THEN
            d(ids)%mortar(1,ifaces) = nmorta 
          ENDIF
        ENDIF
 
      ENDDO
      write(*,*) 'chick nmort old and new due to periodic boundary conditions',nmort,nmorta
      nmort = nmorta
ENDIF
!      ama      = 0.4
!      twall    = 1/ama**2
      itimeinflowbc=0
!
!  Initialize matrices
!
      bmatu1b = 0.0d0
      bmataa  = 0.0d0
      r       = 0.0d0
      s       = 0.0d0
      t       = 0.0d0
!
!  Compute integration matrices
!
      cxg  = 0.0d0
      ncg  = d(1)%ncg(1)
      CALL gaussq(legendre,ncg,0.0d0,0.0d0,0,endpts,work,cxg,xwgt)
!
      DO i = 1,ncg
         DO j = 1,ncg
           DO l = 1,ncg
            sum = 0.0d0
            DO k = 1,ncg
               hl = polyn(l,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
               hj = polyn(j,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
               hi = polyn(i,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
               sum = sum + 0.5d0*xwgt(k)*hi*hj*hl
            END DO
            bmatu1b(i,j,l) = sum
           END DO
         END DO
      END DO
!
      CALL compute_bmata(bmataa,d(1)%cxg(1:ncg,1),ncg,ncg)

!

      RETURN
      END SUBROUTINE user_setup
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE do_user_stuff(time,ngrid,d,mrtr,nmort,stat_unit)
!
!......................................................................
!     allow the user to do periodic computations/io,
!     such as a probe measurement. Frequency of calling is set by the
!     variable "nout" i.e. by the keyword "frequency"
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
      USE size
      USE particle_definition
      USE part_par
      
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)   :: d(ngp)
      TYPE (mortar)   :: mrtr(nmp)
      type (particle) :: drop(1)
!
      INTEGER :: stat_unit
!
!     Write stuff to file with extension .sta
!
!      IF (myid ==0) THEN
!        u = d(1)%Q(3,3,3,2)
!        v = d(1)%Q(3,3,3,3)
!        w = d(1)%Q(3,3,3,4)
!        rho = d(1)%Q(3,3,3,1)
!        p = d(1)%Q(3,3,3,5)*(gamma-1)-0.5*rho*(u*u+v*v+w*w)
!        write(stat_unit,*)  time, p, itimeinflowbc
!      END IF
!      RETURN

!   do i=1,ngrid
!       ncg = d(i)%ncg(1)
!       write(*,*)'n:',i
!       write(*,*)'max:',maxval(d(i)%Q(1:ncg,1:ncg,1:ncg,2)),'min:',minval(d(i)%Q(1:ncg,1:ncg,1:ncg,2))
!   end do
!


    
      END SUBROUTINE do_user_stuff
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE cleanup(time,ngrid,d)
!
!......................................................................
!     perform problem dependent computations
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
!
!     ----------------------------
!     do nothing
!     ----------------------------
!
!
      RETURN
      END SUBROUTINE cleanup
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initc(id,d,mrtr)
!
!......................................................................
!     Set up the initial conditions as a function of x and y
!     these conditions are for a point source solution
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
!   type(bounds)             :: dbnd
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: fv(5)
      INTEGER                      :: neighbor
      double precision              :: u,v,w,gOffset
    
!
!     -------------------------
!     Inputs 
!     -------------------------
!



!
!     ----------------------------
!     compute interior point values
!     ----------------------------
!
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               x = d%xg(1,n,m,l)
               y = d%xg(2,n,m,l)
               z = d%xg(3,n,m,l)
               CALL exactsol(x,y,z,d%Q(n,m,l,1),d%Q(n,m,l,2), &
                             d%Q(n,m,l,3),d%Q(n,m,l,4),d%Q(n,m,l,5),fv(:))
           !d%Q(n,m,l,1) = rho
           !d%Q(n,m,l,2) = u*rho
           !d%Q(n,m,l,3) = v*rho
           !d%Q(n,m,l,4) = w*rho
           !d%Q(n,m,l,5) = press/(gamma-1.0d0)+0.5d0*d%Q(n,m,l,1)*(u**2+v**2+w**2 )
           !fv(:) = 0.0d0
            END DO 
         END DO
      END DO
!
!     --------------------------
!     compute boundary solutions
!     --------------------------
!
!     ----------
!     face 1 & 2
!     ----------
!
      DO iface = 1,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         DO l = 1,d%ncg(3)
            DO n = 1,d%ncg(1)
               x = mrtr(mrtr_no)%xg(1,n,l)
               y = mrtr(mrtr_no)%xg(2,n,l)
               z = mrtr(mrtr_no)%xg(3,n,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(n,l,1,2),mrtr(mrtr_no)%Q(n,l,2,2), &
                                   mrtr(mrtr_no)%Q(n,l,3,2),mrtr(mrtr_no)%Q(n,l,4,2), &
                                   mrtr(mrtr_no)%Q(n,l,5,2),mrtr(mrtr_no)%fv(n,l,:,2))
               !mrtr(mrtr_no)%Q(n,l,1,2) = rho
               !mrtr(mrtr_no)%Q(n,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(n,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(n,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(n,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(n,l,:,2) = 0
            END DO
         END DO
      END DO
!
!     ----------
!     face 3 & 5
!     ----------
!
      DO iface = 3,5,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               x = mrtr(mrtr_no)%xg(1,n,m)
               y = mrtr(mrtr_no)%xg(2,n,m)
               z = mrtr(mrtr_no)%xg(3,n,m)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(n,m,1,2),mrtr(mrtr_no)%Q(n,m,2,2), &
                                   mrtr(mrtr_no)%Q(n,m,3,2),mrtr(mrtr_no)%Q(n,m,4,2), &
                                   mrtr(mrtr_no)%Q(n,m,5,2),mrtr(mrtr_no)%fv(n,m,:,2))

               !mrtr(mrtr_no)%Q(n,m,1,2) = rho
               !mrtr(mrtr_no)%Q(n,m,2,2) = u*rho
               !mrtr(mrtr_no)%Q(n,m,3,2) = v*rho
               !mrtr(mrtr_no)%Q(n,m,4,2) = w*rho
               !mrtr(mrtr_no)%Q(n,m,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(n,m,:,2) = 0.0d0
            END DO
         END DO
      END DO
!
!     ----------
!     face 4 & 6
!     ----------
!
      DO iface = 4,6,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         IF (iface==6) THEN
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               x = mrtr(mrtr_no)%xg(1,m,l)
               y = mrtr(mrtr_no)%xg(2,m,l)
               z = mrtr(mrtr_no)%xg(3,m,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(m,l,1,2),mrtr(mrtr_no)%Q(m,l,2,2), &
                                   mrtr(mrtr_no)%Q(m,l,3,2),mrtr(mrtr_no)%Q(m,l,4,2), &
                                   mrtr(mrtr_no)%Q(m,l,5,2),mrtr(mrtr_no)%fv(m,l,:,2))

               !mrtr(mrtr_no)%Q(m,l,1,2) = rho
               !mrtr(mrtr_no)%Q(m,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(m,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(m,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(m,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(m,l,:,2) = 0.0d0
            END DO
         END DO
         ELSEIF (iface==4) THEN ! outflow boundary
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               x = mrtr(mrtr_no)%xg(1,m,l)
               y = mrtr(mrtr_no)%xg(2,m,l)
               z = mrtr(mrtr_no)%xg(3,m,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(m,l,1,2),mrtr(mrtr_no)%Q(m,l,2,2), &
                                   mrtr(mrtr_no)%Q(m,l,3,2),mrtr(mrtr_no)%Q(m,l,4,2), &
                                   mrtr(mrtr_no)%Q(m,l,5,2),mrtr(mrtr_no)%fv(m,l,:,2))

               !mrtr(mrtr_no)%Q(m,l,1,2) = rho
               !mrtr(mrtr_no)%Q(m,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(m,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(m,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(m,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(m,l,:,2) = 0.0d0
            END DO
         END DO
         ENDIF
      END DO
      
!
      RETURN
      END SUBROUTINE initc
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Xinitc(id,d,mrtr)
!
!......................................................................
!     Interpolate solution values to new grid
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: cxg(nx),Q(nx,ny,nz,neq),polyn
      INTEGER                      :: p
!
      Q = 0.0d0
      cxg = 0.0d0
!
!     Define old Chebyshev grid parameters
!
      ncg1 = 5
      ncg2 = 5
      ncg3 = 5
      cxg = 0.0d0
      DO n = 1,ncg1
         xx = cos((2.d0*DBLE(n)-1.d0)*pi/(2.d0*(DBLE(ncg1) - 1.d0) + 2.d0))
         cxg(n) = 0.5d0*(1.d0 - xx)
      END DO


!
!     Interpolate Q to refined grid within each domain
!
      DO i=1,d%ncg(1)
        DO j=1,d%ncg(2)
          DO k=1,d%ncg(3)
            X =d%cxg(i,1)
            Y =d%cxg(j,2)
            Z =d%cxg(k,3)
!
            DO nv=1,neq
              som  = 0.0d0
              DO m=1,ncg1
                DO n=1,ncg2
                  DO p=1,ncg3
                    hx=polyn(m,X,ncg1,cxg(:))
                    hy=polyn(n,Y,ncg2,cxg(:))
                    hz=polyn(p,Z,ncg3,cxg(:))
                    som  = som  + d%Q(m,n,p,nv)*hx*hy*hz
                  ENDDO
                ENDDO
              ENDDO
!
              Q(i,j,k,nv)=som
!
            ENDDO
!
          ENDDO
        ENDDO
      ENDDO
!
!   Multiply by jacobian: this one is scaled out again in the the main program
!
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               DO nv = 1,neq
                   d%Q(n,m,l,nv) = Q(n,m,l,nv)*d%jacob(n,m,l)
               END DO
            END DO
        END DO
      END DO
!
!
!
      RETURN
      END SUBROUTINE Xinitc
!
!///////////////////////////////////////////////////////////////////////
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE exactsol(x,y,z,rho,rhou,rhov,rhow,rhoe,fv)
!
!......................................................................
!     compute the initial for the channel flow with periodic bc
!......................................................................
!
      USE constants
      USE User_Data
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION :: fv(5)
      DOUBLE PRECISION :: hlow, hupp, uu, vv, ww, p, rho1, som, initvel, Tref
      DOUBLE PRECISION :: xx(16),uuu(16)
      DOUBLE PRECISION :: xx1(38),uuu1(38)
      DOUBLE PRECISION :: polyn, hz
      double precision, parameter :: pie = 3.141592653589793238462643383279502d0
      double precision  :: radius, smooth, Mar
!
   Mar = 0.4d0
   
   Temp = 0.1d0*(sin(pie*x)+sin(pie*y)) + 1.0d0
   uu = Mar*Mar*T
   vv = 0.0d0
   ww = 0.0d0
   P = 1.0d0

    rho     = P*gamma/Temp
    rhou 	= rho*uu
    rhov 	= rho*vv
    rhow 	= rho*ww
    rhoe 	= p/(gamma-1.0d0) + 0.5d0*rho*(uu**2 + vv**2 + ww**2)
    fv(:) = 0.0d0
!
      RETURN
      END SUBROUTINE exactsol
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE compute_source(x,y,z,time,source,delta_t)
!
!......................................................................
!     physical source term
!......................................................................
!
      USE size
      USE physics
      USE User_Data
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(neq) :: source
      double precision                  :: coef, xshift, yshift, scale
!
!         coef    = 10.0d0
!         scale   = 0.020d0
!         xshift  = 1.0d0
!         yshift  = 1.0d0
!
!         source(1) = 0.0d0
!         source(2) = coef * exp( -( (x-xshift)**2 + (y-yshift)**2 ) / scale )
!         source(3) = 0.0d0
!         source(4) = 0.0d0
!         source(5) = 10.0d0*coef * exp( -( (x-xshift)**2 + (y-yshift)**2 ) / scale )
      RETURN
      END SUBROUTINE compute_source
!
!///////////////////////////////////////////////////////////////////////
!
!
   SUBROUTINE boundary_flux(time,mrtr,bcond,property,num_props,flux)
!
!......................................................................
!     date: 11/17/98
!
!     routines called:
!                      wallbc
!                      inflow
!                      outflow
!
!     applicability: all
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('walladiab')
            CALL wallbc_Inv(mrtr,flux)
         CASE('wallisoth')
            CALL wallbc_isoth(mrtr,flux)
         CASE('inflow')
            CALL inflow_specAll(mrtr,flux)
         CASE('outflow')
            CALL outflow_specAll(mrtr,flux)
         CASE('periodm')
            CALL periodm_specAll(mrtr,flux)
         CASE('periods')
            CALL periods_specAll(mrtr,flux)
      END SELECT
!
      RETURN
      END SUBROUTINE boundary_flux

!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DirichletBC(time,mrtr,bcond,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     routines called:
!                      boundary_values
!                      set_adiabatic wall
!
!     applicability: 3D Navier-Stokes
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('wallisoth')
            CALL set_isothermal_wall(mrtr,Qbound)
         CASE('walladiab')
            CALL set_adiabatic_wall(mrtr,Qbound)
         CASE('inflow')
            CALL boundary_values(mrtr,Qbound)
         CASE('outflow')
            CALL boundary_values(mrtr,Qbound)
         CASE('periodm')
            CALL set_periodm(mrtr,Qbound)
         CASE('periods')
            CALL set_periods(mrtr,Qbound)
      END SELECT
!
      RETURN
      END SUBROUTINE DirichletBC
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Vis_Neumann(time,mrtr,bcond,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     routines called:
!                      boundary_values
!                      set_adiabatic wall
!
!     applicability: 3D Navier-Stokes
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('wallisoth')
            CALL wall_vis_flux(mrtr,Qbound,bcond)
         CASE('walladiab')
            CALL wall_vis_flux(mrtr,Qbound,bcond)
         CASE('inflow')
            CALL boundary_vis_flux(mrtr,Qbound)
         CASE('outflow')
            CALL boundary_vis_flux(mrtr,Qbound)
      END SELECT
!
      RETURN
      END SUBROUTINE Vis_Neumann
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE set_part_par()

!
!    initialization of particle parameters
!
      USE size
      USE constants
      USE part_par
      USE physics
!
!
!     parameters not defined in the euler code
!
      IEVAP  = 0

      SCF    = 0.7d0
      SCO    = 0.7d0
!      DA     = 6.0d0
!      ZE     = 4.0d0
!      CE     = 40.0d0
      RCOEF  = 1.0d0
      YFINIT = 0.0d0
!
!     compute particle parameters (see formulation)
!
      PEF    = SCF*re
      PD0    = dsqrt(18.0d0*TAUP0/(re*RHOP))
      PM0    = pi*RHOP*PD0**3 /6.d0
      PSI    = dble(nprt)*PM0/(2.d0*pi)**3
      CRE    = re*(6.d0/(pi*RHOP))**(1.0d0/3.0d0)
      CTAUP  = (RHOP**(1.0d0/3.0d0) *re/18.0d0)*(6.0d0/pi)**(2.0d0/3.0d0)
      CYFPS  = gamma*A1/((gamma-1.)*TBOIL)
      CF2    = 1.0d0/(3.0d0*pr*A2)
      CF3    = A1/(3.0d0*SCF*A2)
      CF4    = pi*dsqrt(18.0d0/RHOP)/(re**(3.0d0/2.0d0) *SCF)

!   
      nprtmax = 200000

      RETURN
   END SUBROUTINE set_part_par
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE part_initc(drop,nprt)

!
!    initialization of particle parameters
!
      USE size
!
      USE particle_definition
!
      USE constants
      USE part_par
      USE physics
      USE User_Data
!
      TYPE(particle)  :: drop
!
!    initialize particle field
!
      DO j=1,3
         drop%Vp(j)     = 0.0d0
         drop%Vfp(j)    = 0.0d0
      END DO
         drop%Tp       = 14.0d0
         drop%Mp       = PM0
         drop%Tfp      = 14.0d0
         drop%Rhofp    = 1.0d0
         drop%Yffp     = 1.0d0
!
         drop%onoff     = 1
!

      RETURN
   END SUBROUTINE part_initc
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE find_processor(xmid,ymid,zmid,id, nno,ngrid)
!
       USE constants
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
! structured uniform grid in x-direction, and y-direction
! cosine distributiono in z-directions
           CALL read_metis_file(id, nno,ngrid)
      ! ELSE
      !     domlx = 4.0d0
      !     domly = 2.0d0
      !     domlz = 2.0d0
      !     nprocx = 0
      !     nprocy = 1
      !     nprocz = 1
      !     zmid   = 2.0d0*(ACOS(1.0d0-2.0d0*(zmid/2.0d0))/pi)
!
!           n1   = INT(xmid*REAL(nprocx)/domlx)+1
!           n2   = INT(ymid*REAL(nprocy)/domly)+1
!           n3   = INT(zmid*REAL(nprocz)/domlz)+1
!           nno  = (n1-1)*nprocy+n2-1+(n3-1)*nprocx*nprocy  
!
!           IF (ymid < 3.0) THEN
!             nno = 0
!           ELSE
!             nno = 1
!           ENDIF
!     ENDIF

!
      RETURN
   END SUBROUTINE find_processor
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE read_restart_case(rst_unit,time)
!
! read case specific stuff from the restart file
!
       USE User_Data
       USE mpi_par
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!
       INCLUDE "mpif.h"
!
      IF (myid == 0) THEN
         READ(rst_unit) time
         write(*,*) 'time',time
      ENDIF
!
      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
!

       RETURN
      END SUBROUTINE read_restart_case
 
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_case(rst_unit,time)
!
! Write case specific stuff to the restart file
!
       USE User_Data
       USE mpi_par
!       
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!
       INCLUDE "mpif.h"
    
      IF (myid == 0) THEN
         WRITE(rst_unit) time
         write(*,*) 'final time',time
      ENDIF
!
       RETURN
      END SUBROUTINE write_restart_case
!///////////////////////////////////////////////////////////////////

      SUBROUTINE inflow_time_bc(id,d,mrtr,k)

      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: fv(5)
      INTEGER                      :: neighbor


      END SUBROUTINE inflow_time_bc
!/////////////////////////////////////////////////////////////

      SUBROUTINE inflow_cond(x,y,z,rho,rhou,rhov,rhow,rhoe,fv)

      USE constants
      USE User_Data
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION :: fv(5)
      DOUBLE PRECISION :: hlow, hupp, uu, vv, ww, p, rho1, som
      DOUBLE PRECISION :: polyn, hz



      END SUBROUTINE inflow_cond 

!/////////////////////////////////////////////////////////////
      !! Particle inflow subroutines
subroutine partInflow(d,p,dt,ngridl,nprt,myid)
    use domain_definition
    use size 
    use physics
    use particle_definition
    use particleInlets
    implicit none
    ! Domain Types
    type(domain), dimension(ngp)                 :: d            !< Carrier Domain
    ! Particle Types
    type(particle), dimension(npart)             :: p            !< Particle
    
    double precision                             :: uniformRandom
    
    integer, intent(in)                          :: ngridl       !< Local elements
    integer, intent(inout)                       :: nprt
    integer, intent(in)                          :: myid
    double precision, intent(in)                 :: dt
    integer                                      :: i,j,k,id,iface,padd,np
    double precision                             :: u, paddD, yposm, ypos, y1, y2, dy, x1, x2, dx, xpos
    double precision                             :: Fui, Oxi, Nii, Pdi, Fuc, Oxc, Nic, Pdc
    double precision, parameter                  :: Lx1 = 0.38448375335127d0 !< length of the first element
    double precision, parameter                  :: Lx2 = 0.04250d0
    integer, parameter                           :: ideal = 80
    

    
end subroutine partInflow

subroutine partOutflow(d,p,dt,ngridl,nprt,myid)
    use domain_definition
    use size 
    use physics
    use particle_definition
    use particleInlets
    implicit none
    ! Domain Types
    type(domain), dimension(ngp)                 :: d            !< Carrier Domain
    ! Particle Types
    type(particle), dimension(npart)             :: p            !< Particle
    
    integer, intent(in)                          :: ngridl       !< Local elements
    integer, intent(inout)                       :: nprt
    integer, intent(in)                          :: myid
    double precision, intent(in)                 :: dt
    
    integer                                      :: np
    
    
    
end subroutine partOutflow

subroutine initScalars(d,eP,id)
    use size
    use particle_definition
    use domain_definition
    implicit none
    type(domain)                                :: d
    type(ePointers), dimension(ngp,npart)       :: eP
    double precision                            :: x,y,z
    integer                                     :: np, id
    double precision, parameter                 :: pie		= 3.14159265358979323846264338327d0
    double precision                            :: delta
    double precision                            :: T,F,O,N,P,spike
    
    
    do np=1,d%np
        x = eP(id,np)%p%Xp(1)
        y = eP(id,np)%p%Xp(2)
        z = eP(id,np)%p%Xp(3)
        
		!For CV Reactor
        !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        !eP(id,np)%p%scalar(2) = 0.01152d0
        !eP(id,np)%p%scalar(3) = 0.23042d0
        !eP(id,np)%p%scalar(4) = 0.0d0
        !ep(id,np)%p%scalar(5) = 0.75806d0
        
!         !for CV reactor
!         F     = 0.01152d0
!         O     = 0.23042d0
!         N     = 0.75806d0
!         P   = 0.0d0
        !for ramp cavity with injector
        ! F     = 0.000001d0
!         O     = 0.24194d0
!         N     = 0.75806d0
!         P     = 0.0d0

        !for Sod Problem
        if (x > 0.50d0) then
            F     = 0.0115d0
            O     = 0.2304d0
            N     = 0.7581d0
            P   = 0.0d0
        else
            F = 0.0d0
            O = 0.23310d0
            N = 0.76690d0
            P = 0.0d0
        end if

        eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        eP(id,np)%p%scalar(2) = F
        eP(id,np)%p%scalar(3) = O
        eP(id,np)%p%scalar(4) = P
        eP(id,np)%p%scalar(5) = N
        eP(id,np)%p%react = 0.0d0

        ! for ramp cavity
        !T = eP(id,np)%p%Tfp
        !
        !if(x.le.3) then
        !    !T = 1.0d0*(spike-(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N
        !else if(x.ge.6) then
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N
        !!else if(y.le.dsqrt(1.50d0**2.0-(x-4.50d0)**2.0)+0.650d0) then
        !else if(y.le.0.650d0) then
        !    !T = 1.0d0*(spike+(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = 0.0115d0
        !    O = 0.2304d0
        !    N = 0.7581d0
        !    P = 1.0d0-F-O-N
        !else
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N 
        !end if
        !
        !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        !eP(id,np)%p%scalar(2) = F
        !eP(id,np)%p%scalar(3) = O
        !eP(id,np)%p%scalar(4) = P
        !eP(id,np)%p%scalar(5) = N
        ! end of for cavity

        !eP(id,np)%p%scalar(1) = 3.0
        !eP(id,np)%p%scalar(2) = 4.0
        
        !if (y .ge. 0.0d0) then
        !    ep(id,np)%p%scalar(1) = erf(3.14159d0**0.5d0 * y / delta)
        !else
        !    ep(id,np)%p%scalar(1) = 0.0d0
        !end if
        
        !if (y .ge. 0.0d0) then
        !    eP(id,np)%p%scalar(1) = 0.0d0
        !    ep(id,np)%p%scalar(2) = 1.0!erf(sqrt(pie)*y/delta)
        !end if
        !if (y .lt. 0.0d0) then
        !    eP(id,np)%p%scalar(2) = 0.0d0
        !    ep(id,np)%p%scalar(1) = 1.0!- erf(sqrt(pie)*y/delta)
        !end if
        !
        !eP(id,np)%p%scalar(3) = eP(id,np)%p%Rhofp
        
    end do
    
end subroutine initScalars

    