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
!
      SAVE
!
      INTEGER                           :: num_user_keywords = 4
      CHARACTER (LEN=10), DIMENSION(4)  :: user_keyword_list = (/ &
                                                'theta     '    , &
                                                'phi       '    , &
                                                'outlet    '    , &
                                                'slide     '      &
                                           /)
      LOGICAL, DIMENSION(4) :: user_keyword_set = .false.
      END MODULE user_keywords
!
!///////////////////////////////////////////////////////////////////////
!
      MODULE User_Data
         USE size
         SAVE
         DOUBLE PRECISION :: theta = 0.0d0,phi = 0.0d0
         DOUBLE PRECISION :: pout,prestep
         DOUBLE PRECISION :: upper_velocity
         DOUBLE PRECISION :: diss(ngp,nx,ny,nz,2) = 0.0d0 
         DOUBLE PRECISION :: ima((nfour/6+1)**3,ngp,nx) = 0.0d0 !minimum of six domains
         DOUBLE PRECISION :: imb((nfour/6+1)**3,ngp,nx) = 0.0d0
         DOUBLE PRECISION :: imc((nfour/6+1)**3,ngp,nx) = 0.0d0
         DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t        ! viscous
      END MODULE User_Data
!
!///////////////////////////////////////////////////////////////////////
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
         CASE('slide')
            upper_velocity = get_dp_value(input_line)
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
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)
      INTEGER,DIMENSION(1000) :: mortmat
      LOGICAL               :: nodouble,restart

      mortmat(:) = 0
!
!  set periodic boundary conditions
!

!  find the interfaces  at which periodic boundary conditions
!  apply, adjust domain and mortar parameters

IF (.not. restart) THEN

      chick = 1
      IF (chick ==1) THEN
      kkk = 0
      DO idm = 1,ngrid
         IF (d(idm)%bcond(3) == 'inflow') THEN
            d(idm)%bcond(3)  = 'interface'
            d(idm)%ibtype(3) = 1
            mrtr_no          = d(idm)%mortar(1,3)
            d(idm)%which_surface(3) = 0
            DO ids = 1,ngrid
               IF (d(ids)%bcond(5) == 'outflow'              .and.  &
                   d(idm)%corner(1,1) == d(ids)%corner(1,5) .and.  &
                   d(idm)%corner(2,1) == d(ids)%corner(2,5)) THEN

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,5)
                   mortmat(kkk) = mrtr_no_sl

                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids
                   mrtr(mrtr_no)%iface(2) = 5
                   d(ids)%bcond(5)  = 'interface'
                   d(ids)%ibtype(5) = 1
                   d(ids)%which_surface(5) = 0
                   d(ids)%mortar(1,5) = mrtr_no
                   d(ids)%mortar(2,5) = 2
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO idm = 1,ngrid
         IF (d(idm)%bcond(6) == 'inflow') THEN
            d(idm)%bcond(6)  = 'interface'
            d(idm)%ibtype(6) = 1
            mrtr_no          = d(idm)%mortar(1,6)
            d(idm)%which_surface(6) = 0
            DO ids = 1,ngrid
               IF (d(ids)%bcond(4) == 'outflow'              .and.  &
                   d(idm)%corner(2,1) == d(ids)%corner(2,2) .and.  &
                   d(idm)%corner(3,1) == d(ids)%corner(3,2)) THEN

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,4)
                   mortmat(kkk) = mrtr_no_sl

                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids
                   mrtr(mrtr_no)%iface(2) = 4
                   d(ids)%bcond(4)  = 'interface'
                   d(ids)%ibtype(4) = 1
                   d(ids)%which_surface(4) = 0
                   d(ids)%mortar(1,4) = mrtr_no
                   d(ids)%mortar(2,4) = 2
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO idm = 1,ngrid
         IF (d(idm)%bcond(1) == 'inflow') THEN
            d(idm)%bcond(1)  = 'interface'
            d(idm)%ibtype(1) = 1
            mrtr_no          = d(idm)%mortar(1,1)
            d(idm)%which_surface(1) = 0
            DO ids = 1,ngrid
               IF (d(ids)%bcond(2) == 'outflow'              .and.  &
                   d(idm)%corner(1,1) == d(ids)%corner(1,4) .and.  &
                   d(idm)%corner(3,1) == d(ids)%corner(3,4)) THEN

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
      nmortdouble = kkk
!  Clean up mortars
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
      write(*,*) 'chick nmort old and new',nmort,nmorta
      nmort = nmorta
      ENDIF
ENDIF

     DO id=1,ngrid
       CALL  interp_matrix(d(id)%corner,d(id)%cxg(:,1),d(id)%ncg(1),ngrid,id)
     ENDDO
!
!
      RETURN
      END SUBROUTINE user_setup
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE do_user_stuff(time,dt,ngrid,d,mrtr,nmort,stats,stat_unit,m)
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
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
!
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)
!

!
      DOUBLE PRECISION                  :: xwgt(nmax),cxg(nmax),dt
      DOUBLE PRECISION                  :: bmat(nmax,nmax,nmax)
      INTEGER, PARAMETER                :: legendre = 1
      INTEGER                           :: o,m
      INTEGER,DIMENSION(ngp)            :: npoin
     
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
      DOUBLE PRECISION, DIMENSION(9)    :: statall,statalla
      DOUBLE PRECISION                  :: polyn
      REAL                              :: som,som1,som2
      CHARACTER(LEN=32)                 :: fname
      COMPLEX, DIMENSION(nfour,nfour,nfour) :: ufour,vfour,wfour,udum

      INTEGER :: stat_unit

    if(m==1) then
        prestep = 0.0d0
    end if
if( mod(m,1)==0) then

    if (myid==0) write(*,*) 'do user stuff',time
!
!     ------------------------------
!     compute dissipation from velocity derivatives
!     ------------------------------
!
      DO id=1,ngrid
        CALL compute_diss(d(id),id)
      ENDDO
!
!     ------------------------------
!     set up gauss grids and weights
!     ------------------------------
!
      bmat = 0.0d0
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
               sum = sum + 0.1250d00*xwgt(k)*xwgt(m)*xwgt(n)*hi*hj*hl
            END DO
            bmat(i,j,l) = sum
           END DO
         END DO
      END DO
!
!    ----------------------
!    compute averages
!    ----------------------
!
   uav = 0.0d0
   vav = 0.0d0
   wav = 0.0d0
   eav = 0.0d0
   rav = 0.0d0 
   tav = 0.0d0
   epav= 0.0d0
   ediss1= 0.0d0
   ediss2= 0.0d0
!
   DO id =1, ngrid
      DO i = 1,ncg
        DO j = 1,ncg
          DO l = 1,ncg
            u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
            v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
            w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
            rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
            temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
            temp = temp - 0.5d0*(u**2+v**2+w**2)
            temp = temp*gamma*(gamma-1.0d0)
            press = temp*rho/gamma

            uav = uav + u*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
            vav = vav + v*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
            wav = wav + w*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
            eav = eav + d(id)%Q(i,j,l,5)*bmat(i,j,l)*(d(id)%jacob(i,j,l))
            rav = rav + rho*bmat(i,j,l)
            tav = tav + temp*bmat(i,j,l)
            epav= epav + press/(gamma-1.0d0)*bmat(i,j,l)
            ediss1 = ediss1 + diss(id,i,j,l,1)*bmat(i,j,l)
            ediss2 = ediss2 + diss(id,i,j,l,2)*bmat(i,j,l)
          ENDDO
        ENDDO
      ENDDO
   ENDDO
!
   uav = uav/d(1)%jacob(1,1,1)
   vav = vav/d(1)%jacob(1,1,1)
   wav = wav/d(1)%jacob(1,1,1)
   eav = eav/d(1)%jacob(1,1,1)
   rav = rav/d(1)%jacob(1,1,1)
   tav = tav/d(1)%jacob(1,1,1)
   epav= epav/d(1)%jacob(1,1,1)
   ediss1= ediss1/d(1)%jacob(1,1,1)
   ediss2= ediss2/d(1)%jacob(1,1,1)
!
   statall(1) = uav
   statall(2) = vav
   statall(3) = wav
   statall(4) = eav
   statall(5) = rav
   statall(6) = tav
   statall(7) = epav
   statall(8) = ediss1
   statall(9) = ediss2
!
   CALL MPI_ALLREDUCE(statall,statalla,9,MPI_DOUBLE_PRECISION, MPI_SUM, &
                          comm1d, ierr)
!   
   IF (myid==0) THEN
     uav = statalla(1)
     vav = statalla(2)
     wav = statalla(3)
     eav = statalla(4)
     rav = statalla(5)
     tav = statalla(6)
     epav= statalla(7)
     ediss1= statalla(8)
     ediss2= statalla(9)

     uav = uav/(8.0d0*pi**3)
     vav = vav/(8.0d0*pi**3)
     wav = wav/(8.0d0*pi**3)
     eav = eav/(8.0d0*pi**3)
     rav = rav/(8.0d0*pi**3)
     tav = tav/(8.0d0*pi**3)
     epav= epav/(8.0d0*pi**3)
     ediss1= ediss1/(8.0d0*pi**3)
     ediss2= ediss2/(8.0d0*pi**3)

     statall(1) = uav
     statall(2) = vav
     statall(3) = wav
     statall(4) = eav
     statall(5) = rav
     statall(6) = tav
     statall(7) = epav
     statall(8) = ediss1
     statall(9) = ediss2
   ENDIF
!
   CALL MPI_BCAST(statall,9,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
!
!    -------------------------
!    compute root mean squares
!    -------------------------
!
   rav = statall(5)
   tav = statall(6)
   epav= statall(7)
!
   urms = 0.0d0
   vrms = 0.0d0
   wrms = 0.0d0
   rrms = 0.0d0
   trms = 0.0d0
   amarms= 0.0d0
   eprms= 0.0d0
!
   DO id =1, ngrid
      DO i = 1,ncg
        DO j = 1,ncg
          DO l = 1,ncg
            u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
            v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
            w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
            rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
            temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
            temp = temp - 0.5d0*(u**2+v**2+w**2)
            temp = temp*gamma*(gamma-1.0d0)
            press= temp*rho/gamma

            rrms = rrms + ((rho-rav)**2)*bmat(i,j,l)
            trms = trms + ((temp-tav)**2)*bmat(i,j,l)
            urms = urms + (u**2)*bmat(i,j,l)
            vrms = vrms + (v**2)*bmat(i,j,l)
            wrms = wrms + (w**2)*bmat(i,j,l)
            eprms= eprms + ((press/(gamma-1.0d0)-epav)**2)*bmat(i,j,l)
            amarms= amarms + (u**2+v**2+w**2)*bmat(i,j,l)/tav
          ENDDO
        ENDDO
      ENDDO
   ENDDO
!
   rrms = rrms/d(1)%jacob(1,1,1)
   trms = rrms/d(1)%jacob(1,1,1)
   urms = urms/d(1)%jacob(1,1,1)
   vrms = vrms/d(1)%jacob(1,1,1)
   wrms = wrms/d(1)%jacob(1,1,1)
   eprms = eprms/d(1)%jacob(1,1,1)
   amarms = amarms/d(1)%jacob(1,1,1)
!
   statall(1) = urms
   statall(2) = vrms
   statall(3) = wrms
   statall(4) = rrms
   statall(5) = trms
   statall(6) = eprms
   statall(7) = amarms
!
   CALL MPI_ALLREDUCE(statall,statalla,9,MPI_DOUBLE_PRECISION, MPI_SUM, &
                          comm1d, ierr)
!
   IF (myid==0) THEN
     urms = statalla(1)
     vrms = statalla(2)
     wrms = statalla(3)
     rrms = statalla(4)
     trms = statalla(5)
     eprms = statalla(6)
     amarms = statalla(7)

     turbk= 0.5d0*(urms+vrms+wrms)/(8.0d0*pi**3)
     urms = sqrt(urms/(8.0d0*pi**3))
     vrms = sqrt(vrms/(8.0d0*pi**3))
     wrms = sqrt(wrms/(8.0d0*pi**3))
     rrms = sqrt(rrms/(8.0d0*pi**3))
     trms = sqrt(trms/(8.0d0*pi**3))
     eprms = sqrt(eprms/(8.0d0*pi**3))
     amarms = sqrt(amarms/(8.0d0*pi**3))
     dissp = (prestep-turbk)/(dt)
     prestep = turbk

     ! WRITE(20,*) time,uav,vav,wav,tav,rav,eav,epav,ediss1,ediss2
     ! WRITE(22,*) time,urms,vrms,wrms,turbk,trms,rrms,eprms,amarms
      WRITE(22,*) time,turbk
      write(23,*) time,dissp
   ENDIF
end if
!
nspec=0
IF (nspec==1 .and. m==21) THEN
      ufour=(0.,0.)
      vfour=(0.,0.)
      wfour=(0.,0.)
      npoin=0
      DO i=1,nfour
        DO j=1,nfour
          DO k=1,nfour
           xfour =DBLE(i-1)*2.0d0*pi/DBLE(nfour-1)
           yfour =DBLE(j-1)*2.0d0*pi/DBLE(nfour-1)
           zfour =DBLE(k-1)*2.0d0*pi/DBLE(nfour-1)
           DO id=1,ngrid
             IF (xfour>=d(id)%corner(1,1) .and. xfour<=d(id)%corner(1,2) .and. &
                 yfour>=d(id)%corner(2,1) .and. yfour<=d(id)%corner(2,3) .and. &
                 zfour>=d(id)%corner(3,1) .and. zfour<=d(id)%corner(3,5) ) THEN
             npoin(id) = npoin(id) + 1
             npoin1  = npoin(id)
             som =0.0d0
             som1=0.0d0
             som2=0.0d0
             xm=(xfour-d(id)%corner(1,1))/(d(id)%corner(1,2)-d(id)%corner(1,1))
             ym=(yfour-d(id)%corner(2,1))/(d(id)%corner(2,3)-d(id)%corner(2,1))
             zm=(zfour-d(id)%corner(3,1))/(d(id)%corner(3,5)-d(id)%corner(3,1))
             ww=1
             DO m=1,ncg 
               DO n=1,ncg
                  DO o=1,ncg
                    u    = d(id)%Q(m,n,o,2)/d(id)%Q(m,n,o,1)
                    v    = d(id)%Q(m,n,o,3)/d(id)%Q(m,n,o,1)
                    w    = d(id)%Q(m,n,o,4)/d(id)%Q(m,n,o,1)
                     ! ww   = ima(npoin1,id,m)* &
                     !        imb(npoin1,id,n)* &
                     !        imc(npoin1,id,o)
                      som  = som  + u*ww
                      som1 = som1 + v*ww
                      som2 = som2 + w*ww
                      ! som = 1.0d0
                      ! som1 = 2.0d0
                      ! som2 = 3.0d0


                  ENDDO
               ENDDO
             ENDDO
             ufour(i,j,k)=CMPLX(som,0)
             vfour(i,j,k)=CMPLX(som1,0)
             wfour(i,j,k)=CMPLX(som2,0)
           ENDIF
            ENDDO
           ENDDO
         ENDDO
       ENDDO

      CALL MPI_ALLREDUCE(ufour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) ufour=udum
      CALL MPI_ALLREDUCE(vfour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) vfour=udum
      CALL MPI_ALLREDUCE(wfour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) wfour=udum

       IF (myid==0) THEN
        no=INT(time*100.0)
         write(*,*) som,som1,som2
         WRITE(fname, fmt='(a4,i3,a3)') 'fort',23,'.23'
         OPEN(unit=23,file=fname)
            CALL compute_espec(ufour,vfour,wfour,nfour,nfour,nfour,nfour,time)
         CLOSE(23)
        ENDIF
 ENDIF
 
      RETURN
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
      SUBROUTINE initc(id,d,ngrid)
!
!......................................................................
!     Set up the initial conditions as a function of x and y
!     these conditions are for a point source solution
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d(ngp)
      INTEGER, PARAMETER           :: nnc=8
      DOUBLE PRECISION             :: fv(5)
      INTEGER                      :: neighbor,nna(3),nnb(3)
!     DOUBLE PRECISION             :: vel_u1(1:nfour+nnc-2,1:nfour+nnc-2,1:nfour+nnc-2)
!     DOUBLE PRECISION             :: vel_u2(1:nfour+nnc-2,1:nfour+nnc-2,1:nfour+nnc-2)
!     DOUBLE PRECISION             :: vel_u3(1:nfour+nnc-2,1:nfour+nnc-2,1:nfour+nnc-2)
!     DOUBLE PRECISION             :: vel_dens(1:nfour+nnc-2,1:nfour+nnc-2,1:nfour+nnc-2)
!     DOUBLE PRECISION             :: vel_temp(1:nfour+nnc-2,1:nfour+nnc-2,1:nfour+nnc-2)
!     DOUBLE PRECISION             :: xfour(1:nfour+nnc-2)
      DOUBLE PRECISION             :: vel_u1(-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1)
      DOUBLE PRECISION             :: vel_u2(-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1)
      DOUBLE PRECISION             :: vel_u3(-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1)
      DOUBLE PRECISION             :: vel_dens(-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1)
      DOUBLE PRECISION             :: vel_temp(-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1,-nnc/2+2:nfour+nnc/2-1)
      DOUBLE PRECISION             :: xfour(-nnc/2+2:nfour+nnc/2-1)
!     DOUBLE PRECISION,ALLOCATABLE :: uint(:,:,:,:),xint(:,:)
      DOUBLE PRECISION             :: uint(nnc,nnc,nnc,5),xint(nnc,3)
      REAL                         :: som,som1,som2
      DOUBLE PRECISION             :: dumarray1(nfour,nfour,nfour)
      DOUBLE PRECISION             :: dumarray2(nfour,nfour,nfour)
      DOUBLE PRECISION             :: dumarray3(nfour,nfour,nfour)
      DOUBLE PRECISION             :: dumarray4(nfour,nfour,nfour)
      DOUBLE PRECISION             :: dumarray5(nfour,nfour,nfour)
!
!
      vel_u1    =0.0d0
      vel_u2    =0.0d0
      vel_u3    =0.0d0
      vel_dens  =0.0d0
      vel_temp  =0.0d0
      xfour     =0.0d0
!
!     -------------------------------------
!     compute initial velocity distribution
!     uv_init is a single precision routine
!     -------------------------------------
!
!      CALL  uv_init(nfour,vel_u1(1:nfour,1:nfour,1:nfour),vel_u2(1:nfour,1:nfour,1:nfour),  &
!                   vel_u3(1:nfour,1:nfour,1:nfour),vel_dens(1:nfour,1:nfour,1:nfour),    &
!                   vel_temp(1:nfour,1:nfour,1:nfour))
write(*,*) 'before uv_init'
!
      CALL  uv_init(nfour,dumarray1,dumarray2,dumarray3,dumarray4,dumarray5)
write(*,*) 'after uv_init'
      vel_u1(1:nfour,1:nfour,1:nfour) =dumarray1
      vel_u2(1:nfour,1:nfour,1:nfour) =dumarray2
      vel_u3(1:nfour,1:nfour,1:nfour) =dumarray3
      vel_dens(1:nfour,1:nfour,1:nfour) =dumarray4
      vel_temp(1:nfour,1:nfour,1:nfour) =dumarray5
!
      DO i=1,nfour
        xfour(i) =REAL(i-1)*2.0d0*pi/REAL(nfour-1)
      ENDDO
!open(file="N64.dat",unit=23)
!write(23,*) vel_u1
!write(23,*) vel_u2
!write(23,*) vel_u3
!write(23,*) vel_dens
!write(23,*) vel_temp
!close(23)
!
!     -----------------------------------
!     expand domain for interpolation use
!     -----------------------------------
!
      nnca=nnc/2-1
      IF (nnca>0) THEN
! coordinates
        xfour(-nnca+1:0)=xfour(nfour-nnca:nfour-1)-2.0d0*pi
        xfour(nfour+1:nfour+nnca)=xfour(2:nnca+1)+2.0d0*pi
! sides
        vel_u1(-nnca+1:0,:,:) = vel_u1(nfour-nnca:nfour-1,:,:)
        vel_u2(-nnca+1:0,:,:) = vel_u2(nfour-nnca:nfour-1,:,:)
        vel_u3(-nnca+1:0,:,:) = vel_u3(nfour-nnca:nfour-1,:,:)
        vel_dens(-nnca+1:0,:,:) = vel_dens(nfour-nnca:nfour-1,:,:)
        vel_temp(-nnca+1:0,:,:) = vel_temp(nfour-nnca:nfour-1,:,:)

        vel_u1(:,-nnca+1:0,:) = vel_u1(:,nfour-nnca:nfour-1,:)
        vel_u2(:,-nnca+1:0,:) = vel_u2(:,nfour-nnca:nfour-1,:)
        vel_u3(:,-nnca+1:0,:) = vel_u3(:,nfour-nnca:nfour-1,:)
        vel_dens(:,-nnca+1:0,:) = vel_dens(:,nfour-nnca:nfour-1,:)
        vel_temp(:,-nnca+1:0,:) = vel_temp(:,nfour-nnca:nfour-1,:)

        vel_u1(:,:,-nnca+1:0) = vel_u1(:,:,nfour-nnca:nfour-1)
        vel_u2(:,:,-nnca+1:0) = vel_u2(:,:,nfour-nnca:nfour-1)
        vel_u3(:,:,-nnca+1:0) = vel_u3(:,:,nfour-nnca:nfour-1)
        vel_dens(:,:,-nnca+1:0) = vel_dens(:,:,nfour-nnca:nfour-1)
        vel_temp(:,:,-nnca+1:0) = vel_temp(:,:,nfour-nnca:nfour-1)

        vel_u1(nfour+1:nfour+nnca,:,:) = vel_u1(2:nnca+1,:,:)
        vel_u2(nfour+1:nfour+nnca,:,:) = vel_u2(2:nnca+1,:,:)
        vel_u3(nfour+1:nfour+nnca,:,:) = vel_u3(2:nnca+1,:,:)
        vel_dens(nfour+1:nfour+nnca,:,:) = vel_dens(2:nnca+1,:,:)
        vel_temp(nfour+1:nfour+nnca,:,:) = vel_temp(2:nnca+1,:,:)

        vel_u1(:,nfour+1:nfour+nnca,:) = vel_u1(:,2:nnca+1,:)
        vel_u2(:,nfour+1:nfour+nnca,:) = vel_u2(:,2:nnca+1,:)
        vel_u3(:,nfour+1:nfour+nnca,:) = vel_u3(:,2:nnca+1,:)
        vel_dens(:,nfour+1:nfour+nnca,:) = vel_dens(:,2:nnca+1,:)
        vel_temp(:,nfour+1:nfour+nnca,:) = vel_temp(:,2:nnca+1,:)

        vel_u1(:,:,nfour+1:nfour+nnca) = vel_u1(:,:,2:nnca+1)
        vel_u2(:,:,nfour+1:nfour+nnca) = vel_u2(:,:,2:nnca+1)
        vel_u3(:,:,nfour+1:nfour+nnca) = vel_u3(:,:,2:nnca+1)
        vel_dens(:,:,nfour+1:nfour+nnca) = vel_dens(:,:,2:nnca+1)
        vel_temp(:,:,nfour+1:nfour+nnca) = vel_temp(:,:,2:nnca+1)
! corners
        vel_u1(-nnca+1:0,-nnca+1:0,:) = vel_u1(nfour-nnca:nfour-1,nfour-nnca:nfour-1,:)
        vel_u2(-nnca+1:0,-nnca+1:0,:) = vel_u2(nfour-nnca:nfour-1,nfour-nnca:nfour-1,:)
        vel_u3(-nnca+1:0,-nnca+1:0,:) = vel_u3(nfour-nnca:nfour-1,nfour-nnca:nfour-1,:)
        vel_dens(-nnca+1:0,-nnca+1:0,:) = vel_dens(nfour-nnca:nfour-1,nfour-nnca:nfour-1,:)
        vel_temp(-nnca+1:0,-nnca+1:0,:) = vel_temp(nfour-nnca:nfour-1,nfour-nnca:nfour-1,:)

        vel_u1(-nnca+1:0,nfour+1:nfour+nnca,:) = vel_u1(nfour-nnca:nfour-1,2:nnca+1,:)
        vel_u2(-nnca+1:0,nfour+1:nfour+nnca,:) = vel_u2(nfour-nnca:nfour-1,2:nnca+1,:)
        vel_u3(-nnca+1:0,nfour+1:nfour+nnca,:) = vel_u3(nfour-nnca:nfour-1,2:nnca+1,:)
        vel_dens(-nnca+1:0,nfour+1:nfour+nnca,:) = vel_dens(nfour-nnca:nfour-1,2:nnca+1,:)
        vel_temp(-nnca+1:0,nfour+1:nfour+nnca,:) = vel_temp(nfour-nnca:nfour-1,2:nnca+1,:)

        vel_u1(nfour+1:nfour+nnca,-nnca+1:0,:) = vel_u1(2:nnca+1,nfour-nnca:nfour-1,:)
        vel_u2(nfour+1:nfour+nnca,-nnca+1:0,:) = vel_u2(2:nnca+1,nfour-nnca:nfour-1,:)
        vel_u3(nfour+1:nfour+nnca,-nnca+1:0,:) = vel_u3(2:nnca+1,nfour-nnca:nfour-1,:)
        vel_dens(nfour+1:nfour+nnca,-nnca+1:0,:) = vel_dens(2:nnca+1,nfour-nnca:nfour-1,:)
        vel_temp(nfour+1:nfour+nnca,-nnca+1:0,:) = vel_temp(2:nnca+1,nfour-nnca:nfour-1,:)

        vel_u1(nfour+1:nfour+nnca,nfour+1:nfour+nnca,:) = vel_u1(2:nnca+1,2:nnca+1,:)
        vel_u2(nfour+1:nfour+nnca,nfour+1:nfour+nnca,:) = vel_u2(2:nnca+1,2:nnca+1,:)
        vel_u3(nfour+1:nfour+nnca,nfour+1:nfour+nnca,:) = vel_u3(2:nnca+1,2:nnca+1,:)
        vel_dens(nfour+1:nfour+nnca,nfour+1:nfour+nnca,:) = vel_dens(2:nnca+1,2:nnca+1,:)
        vel_temp(nfour+1:nfour+nnca,nfour+1:nfour+nnca,:) = vel_temp(2:nnca+1,2:nnca+1,:)

!x-dir

        vel_u1(:,-nnca+1:0,-nnca+1:0) = vel_u1(:,nfour-nnca:nfour-1,nfour-nnca:nfour-1)
        vel_u2(:,-nnca+1:0,-nnca+1:0) = vel_u2(:,nfour-nnca:nfour-1,nfour-nnca:nfour-1)
        vel_u3(:,-nnca+1:0,-nnca+1:0) = vel_u3(:,nfour-nnca:nfour-1,nfour-nnca:nfour-1)
        vel_dens(:,-nnca+1:0,-nnca+1:0) = vel_dens(:,nfour-nnca:nfour-1,nfour-nnca:nfour-1)
        vel_temp(:,-nnca+1:0,-nnca+1:0) = vel_temp(:,nfour-nnca:nfour-1,nfour-nnca:nfour-1)

        vel_u1(:,-nnca+1:0,nfour+1:nfour+nnca) = vel_u1(:,nfour-nnca:nfour-1,2:nnca+1)
        vel_u2(:,-nnca+1:0,nfour+1:nfour+nnca) = vel_u2(:,nfour-nnca:nfour-1,2:nnca+1)
        vel_u3(:,-nnca+1:0,nfour+1:nfour+nnca) = vel_u3(:,nfour-nnca:nfour-1,2:nnca+1)
        vel_dens(:,-nnca+1:0,nfour+1:nfour+nnca) = vel_dens(:,nfour-nnca:nfour-1,2:nnca+1)
        vel_temp(:,-nnca+1:0,nfour+1:nfour+nnca) = vel_temp(:,nfour-nnca:nfour-1,2:nnca+1)

        vel_u1(:,nfour+1:nfour+nnca,-nnca+1:0) = vel_u1(:,2:nnca+1,nfour-nnca:nfour-1)
        vel_u2(:,nfour+1:nfour+nnca,-nnca+1:0) = vel_u2(:,2:nnca+1,nfour-nnca:nfour-1)
        vel_u3(:,nfour+1:nfour+nnca,-nnca+1:0) = vel_u3(:,2:nnca+1,nfour-nnca:nfour-1)
        vel_dens(:,nfour+1:nfour+nnca,-nnca+1:0) = vel_dens(:,2:nnca+1,nfour-nnca:nfour-1)
        vel_temp(:,nfour+1:nfour+nnca,-nnca+1:0) = vel_temp(:,2:nnca+1,nfour-nnca:nfour-1)

        vel_u1(:,nfour+1:nfour+nnca,nfour+1:nfour+nnca) = vel_u1(:,2:nnca+1,2:nnca+1)
        vel_u2(:,nfour+1:nfour+nnca,nfour+1:nfour+nnca) = vel_u2(:,2:nnca+1,2:nnca+1)
        vel_u3(:,nfour+1:nfour+nnca,nfour+1:nfour+nnca) = vel_u3(:,2:nnca+1,2:nnca+1)
        vel_dens(:,nfour+1:nfour+nnca,nfour+1:nfour+nnca) = vel_dens(:,2:nnca+1,2:nnca+1)
        vel_temp(:,nfour+1:nfour+nnca,nfour+1:nfour+nnca) = vel_temp(:,2:nnca+1,2:nnca+1)

!y-dir
        vel_u1(-nnca+1:0,:,-nnca+1:0) = vel_u1(nfour-nnca:nfour-1,:,nfour-nnca:nfour-1)
        vel_u2(-nnca+1:0,:,-nnca+1:0) = vel_u2(nfour-nnca:nfour-1,:,nfour-nnca:nfour-1)
        vel_u3(-nnca+1:0,:,-nnca+1:0) = vel_u3(nfour-nnca:nfour-1,:,nfour-nnca:nfour-1)
        vel_dens(-nnca+1:0,:,-nnca+1:0) = vel_dens(nfour-nnca:nfour-1,:,nfour-nnca:nfour-1)
        vel_temp(-nnca+1:0,:,-nnca+1:0) = vel_temp(nfour-nnca:nfour-1,:,nfour-nnca:nfour-1)

        vel_u1(-nnca+1:0,:,nfour+1:nfour+nnca) = vel_u1(nfour-nnca:nfour-1,:,2:nnca+1)
        vel_u2(-nnca+1:0,:,nfour+1:nfour+nnca) = vel_u2(nfour-nnca:nfour-1,:,2:nnca+1)
        vel_u3(-nnca+1:0,:,nfour+1:nfour+nnca) = vel_u3(nfour-nnca:nfour-1,:,2:nnca+1)
        vel_dens(-nnca+1:0,:,nfour+1:nfour+nnca) = vel_dens(nfour-nnca:nfour-1,:,2:nnca+1)
        vel_temp(-nnca+1:0,:,nfour+1:nfour+nnca) = vel_temp(nfour-nnca:nfour-1,:,2:nnca+1)

        vel_u1(nfour+1:nfour+nnca,:,-nnca+1:0) = vel_u1(2:nnca+1,:,nfour-nnca:nfour-1)
        vel_u2(nfour+1:nfour+nnca,:,-nnca+1:0) = vel_u2(2:nnca+1,:,nfour-nnca:nfour-1)
        vel_u3(nfour+1:nfour+nnca,:,-nnca+1:0) = vel_u3(2:nnca+1,:,nfour-nnca:nfour-1)
        vel_dens(nfour+1:nfour+nnca,:,-nnca+1:0) = vel_dens(2:nnca+1,:,nfour-nnca:nfour-1)
        vel_temp(nfour+1:nfour+nnca,:,-nnca+1:0) = vel_temp(2:nnca+1,:,nfour-nnca:nfour-1)

        vel_u1(nfour+1:nfour+nnca,:,nfour+1:nfour+nnca) = vel_u1(2:nnca+1,:,2:nnca+1)
        vel_u2(nfour+1:nfour+nnca,:,nfour+1:nfour+nnca) = vel_u2(2:nnca+1,:,2:nnca+1)
        vel_u3(nfour+1:nfour+nnca,:,nfour+1:nfour+nnca) = vel_u3(2:nnca+1,:,2:nnca+1)
        vel_dens(nfour+1:nfour+nnca,:,nfour+1:nfour+nnca) = vel_dens(2:nnca+1,:,2:nnca+1)
        vel_temp(nfour+1:nfour+nnca,:,nfour+1:nfour+nnca) = vel_temp(2:nnca+1,:,2:nnca+1)

      ENDIF
!
!     ----------------------------
!     compute interior point values 
!     with Lagrangian interpolation 
!     of order nnc
!     ----------------------------
!
    DO id=1,ngrid
      DO l = 1,d(id)%ncg(3)
         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               x = d(id)%xg(1,n,m,l)
               y = d(id)%xg(2,n,m,l)
               z = d(id)%xg(3,n,m,l)

               nduma = INT(x/2./pi*REAL(nfour-1))+1 
               ndumb = INT(y/2./pi*REAL(nfour-1))+1 
               ndumc = INT(z/2./pi*REAL(nfour-1))+1 
!
! nnc^th order interp
!
               nna(1)=nduma-(nnc/2-1)
               nna(2)=ndumb-(nnc/2-1)
               nna(3)=ndumc-(nnc/2-1)
               nnb(1)=nna(1)+nnc-1
               nnb(2)=nna(2)+nnc-1
               nnb(3)=nna(3)+nnc-1

               xint(:,1)     = xfour(nna(1):nnb(1))
               xint(:,2)     = xfour(nna(2):nnb(2))
               xint(:,3)     = xfour(nna(3):nnb(3))
               uint(:,:,:,1) = vel_u1(nna(1):nnb(1),nna(2):nnb(2),nna(3):nnb(3))
               uint(:,:,:,2) = vel_u2(nna(1):nnb(1),nna(2):nnb(2),nna(3):nnb(3))
               uint(:,:,:,3) = vel_u3(nna(1):nnb(1),nna(2):nnb(2),nna(3):nnb(3))
               uint(:,:,:,4) = vel_dens(nna(1):nnb(1),nna(2):nnb(2),nna(3):nnb(3))
               uint(:,:,:,5) = vel_temp(nna(1):nnb(1),nna(2):nnb(2),nna(3):nnb(3))

               CALL polin3(xint(:,1),xint(:,2),xint(:,3),uint(:,:,:,1),nnc,nnc,nnc,        &
                            x,y,z,u,dy)
               CALL polin3(xint(:,1),xint(:,2),xint(:,3),uint(:,:,:,2),nnc,nnc,nnc,        &
                            x,y,z,v,dy)
               CALL polin3(xint(:,1),xint(:,2),xint(:,3),uint(:,:,:,3),nnc,nnc,nnc,        &
                            x,y,z,w,dy)
               CALL polin3(xint(:,1),xint(:,2),xint(:,3),uint(:,:,:,4),nnc,nnc,nnc,        &
                            x,y,z,rhofluct,dy)
               CALL polin3(xint(:,1),xint(:,2),xint(:,3),uint(:,:,:,5),nnc,nnc,nnc,        &
                            x,y,z,tempfluct,dy)
               
               rho        = pout*gamma/twall
               rho        = rho + rhofluct
               temp       = twall + tempfluct
               press      = rho*temp/gamma

               d(id)%Q(n,m,l,1) = rho
               d(id)%Q(n,m,l,2) = rho*u
               d(id)%Q(n,m,l,3) = rho*v
               d(id)%Q(n,m,l,4) = rho*w
               d(id)%Q(n,m,l,5) = press/(gamma-1.d0) + 0.5d0*rho*(u**2 + v**2 + w**2) 

            END DO
         END DO
      END DO
    ENDDO
!
!     ----------------------------------------
!     compute the matrix that interpolates
!     the DAK values to a equidistant Fourrier
!     grid
!     ----------------------------------------
!
     DO id=1,ngrid
       CALL  interp_matrix(d(id)%corner,d(id)%cxg(:,1),d(id)%ncg(1),ngrid,id)
     ENDDO
!
      RETURN
      END SUBROUTINE initc
!
!////////////////////////////////////////////////////
!
      SUBROUTINE uv_init(nfour,vel_u1,vel_u2,vel_u3,vel_dens,vel_temp)
!
!
!     The velocity field is initialize using the
!     initialization described in NASA TM 81315 by
!     Rogallo, 'Numerical Experiments in Homogeneous
!     Turbulence'
!     For the dilational component look at Blaisdell,
!     Mansour and Reynolds, Stanfourrd report TF-50,
!     'Numerical simulation of compressible homogeneous
!     turbulence'

!
      USE constants
!
      integer   :: kmax,nx,ny,nz,n1,n2,n3
      parameter(wmol=8.0,wmoh=16.0)
      complex   :: alfa,beta,gam,ran2,ran3,ran4,ran5,ran6
      complex   :: s1(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s2(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s3(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s4(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s5(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: u1(1:nfour,1:nfour,1:nfour)
      complex   :: u2(1:nfour,1:nfour,1:nfour)
      complex   :: u3(1:nfour,1:nfour,1:nfour)
      complex   :: dens(1:nfour,1:nfour,1:nfour)
      complex   :: temp(1:nfour,1:nfour,1:nfour)
      DOUBLE PRECISION   :: vel_u1(nfour,nfour,nfour)
      DOUBLE PRECISION   :: vel_u2(nfour,nfour,nfour)    
      DOUBLE PRECISION   :: vel_u3(nfour,nfour,nfour)                       
      DOUBLE PRECISION   :: vel_dens(nfour,nfour,nfour)                       
      DOUBLE PRECISION   :: vel_temp(nfour,nfour,nfour)                       

!
! input initial paramters
      tav  = 1.0
      chi  = 0.0
      amat = 0.3
      trms = 0.0
      rorms= 0.0
 
! scale the spectra
      ediff= wmoh-wmol
      es   = (1.0-chi)*(amat**2)*tav/ediff
      ed   = chi*es/(1.0-chi)
      et   = (trms**2)/ediff
      er   = (rorms**2)/ediff
       
! start computation initial conditions
      nx = nfour
      ny = nfour
      nz = nfour
      
      kmax=INT(SQRT(2.)*REAL(nx)/3.)                                  

!****Calculate s1 and s2 that is fourier components of u1,u2
      s1=(0.0,0.0)
      s2=(0.0,0.0)
      s3=(0.0,0.0)
      s4=(0.0,0.0) 
      s5=(0.0,0.0)
      DO n3 = -nz/2+1,nz/2
        DO n2 = -ny/2+1,ny/2
         DO n1 = -nx/2+1,nx/2
           IF (n3.eq.0.and.n2.eq.0.and.n1.eq.0) goto 100
           wmag = SQRT(REAL(n1)**2+REAL(n2)**2+REAL(n3)**2)
!*****************Chasnov Spectrum*********
!*****Spectrum****************************
           IF (wmag.ge.wmol.and.wmag.le.wmoh) THEN
             espec=1.0
           ELSE
             espec=0.0
           END IF
!**************************
           CALL RANDOM_NUMBER (dum)
           ran1=2.0*pi*dum
           CALL RANDOM_NUMBER (dum)
           ran2=CMPLX(0.0,dum)
           ran2=2.0*pi*ran2
           CALL RANDOM_NUMBER (dum)
           ran3=CMPLX(0.0,dum)
           ran3=2.0*pi*ran3
           CALL RANDOM_NUMBER (dum)
           ran4=CMPLX(0.0,dum)
           ran4=2.0*pi*ran4
           CALL RANDOM_NUMBER (dum)
           ran5=CMPLX(0.0,dum)
           ran5=2.0*pi*ran5
           CALL RANDOM_NUMBER (dum)
           ran6=CMPLX(0.0,dum)
           ran6=2.0*pi*ran6
!***********************
           alfa = exp(ran2)*cos(ran1)
           alfa = alfa*sqrt(es*espec/(4.0*pi*wmag**2))
           beta = exp(ran3)*sin(ran1)
           beta = beta*sqrt(es*espec/(4.0*pi*wmag**2))
           gam  = sqrt(ed*espec/(4.0*pi*wmag**2))
           gam =  gam*exp(ran4)

           an1  = REAL(n1)
           an2  = REAL(n2)
           an3  = REAL(n3)
           dum  = sqrt(an1**2+an2**2)
       
           IF (n2.eq.0.and.n1.eq.0) THEN
! solenoidal component
             s1(n1,n2,n3)  = alfa
             s2(n1,n2,n3)  = beta
             s3(n1,n2,n3)  = 0.0d0
! dilational component
             s1(n1,n2,n3)  = s1(n1,n2,n3)  + gam*an1/wmag
             s2(n1,n2,n3)  = s2(n1,n2,n3)  + gam*an2/wmag
             s3(n1,n2,n3)  = s3(n1,n2,n3)  + gam*an3/wmag
! density and temperature
             s4(n1,n2,n3)  = (sqrt(er*espec/(4.0*pi*wmag**2)))*exp(ran5)
             s5(n1,n2,n3)  = (sqrt(et*espec/(4.0*pi*wmag**2)))*exp(ran6)
           ELSE
             s1(n1,n2,n3)  = (alfa*wmag*an2+beta*an1*an3)/(wmag*dum)
             s2(n1,n2,n3)  = (beta*an2*an3-alfa*wmag*an1)/(wmag*dum)
             s3(n1,n2,n3)  = -(beta*dum)/wmag
! dilational component
             s1(n1,n2,n3)  = s1(n1,n2,n3)  + gam*an1/wmag
             s2(n1,n2,n3)  = s2(n1,n2,n3)  + gam*an2/wmag
             s3(n1,n2,n3)  = s3(n1,n2,n3)  + gam*an3/wmag
! density and temperature
             s4(n1,n2,n3)  = (sqrt(er*espec/(4.0*pi*wmag**2)))*exp(ran5)
             s5(n1,n2,n3)  = (sqrt(et*espec/(4.0*pi*wmag**2)))*exp(ran6)
           ENDIF
100    ENDDO
      END DO
    ENDDO
    !CALL compute_especa(s1,s2,s3,nx,ny,nz,nfour)

!***************
! reality condition (take the conjugate in an opposing
! quadrant, so that the transform gives no imaginary
! part. This way the energy of the initial spectrum
! is conserved in real space.)
!
    nkut=1
    IF (nkut==1) THEN
    write(*,*) 'chick reality condition'
      DO n3 = 1,nz/2-1
        DO n2 = -ny/2+1,ny/2-1
          DO n1 =-nx/2+1,nx/2-1
            s1(-n1,-n2,-n3) = CMPLX(REAL(s1(n1,n2,n3)),-AIMAG(s1(n1,n2,n3)))
            s1(n1,-n2,-n3)  = CMPLX(REAL(s1(-n1,n2,n3)),-AIMAG(s1(-n1,n2,n3)))
            s1(-n1,n2,-n3)  = CMPLX(REAL(s1(n1,-n2,n3)),-AIMAG(s1(n1,-n2,n3)))
            s1(n1,n2,-n3)   = CMPLX(REAL(s1(-n1,-n2,n3)),-AIMAG(s1(-n1,-n2,n3)))

            s2(-n1,-n2,-n3) = CMPLX(REAL(s2(n1,n2,n3)),-AIMAG(s2(n1,n2,n3)))
            s2(n1,-n2,-n3)  = CMPLX(REAL(s2(-n1,n2,n3)),-AIMAG(s2(-n1,n2,n3)))
            s2(-n1,n2,-n3)  = CMPLX(REAL(s2(n1,-n2,n3)),-AIMAG(s2(n1,-n2,n3)))
            s2(n1,n2,-n3)   = CMPLX(REAL(s2(-n1,-n2,n3)),-AIMAG(s2(-n1,-n2,n3)))

            s3(-n1,-n2,-n3) = CMPLX(REAL(s3(n1,n2,n3)),-AIMAG(s3(n1,n2,n3)))
            s3(n1,-n2,-n3)  = CMPLX(REAL(s3(-n1,n2,n3)),-AIMAG(s3(-n1,n2,n3)))
            s3(-n1,n2,-n3)  = CMPLX(REAL(s3(n1,-n2,n3)),-AIMAG(s3(n1,-n2,n3)))
            s3(n1,n2,-n3)   = CMPLX(REAL(s3(-n1,-n2,n3)),-AIMAG(s3(-n1,-n2,n3)))

            s4(-n1,-n2,-n3) = CMPLX(REAL(s4(n1,n2,n3)),-AIMAG(s4(n1,n2,n3)))
            s4(n1,-n2,-n3)  = CMPLX(REAL(s4(-n1,n2,n3)),-AIMAG(s4(-n1,n2,n3)))
            s4(-n1,n2,-n3)  = CMPLX(REAL(s4(n1,-n2,n3)),-AIMAG(s4(n1,-n2,n3)))
            s4(n1,n2,-n3)   = CMPLX(REAL(s4(-n1,-n2,n3)),-AIMAG(s4(-n1,-n2,n3)))

            s5(-n1,-n2,-n3) = CMPLX(REAL(s5(n1,n2,n3)),-AIMAG(s5(n1,n2,n3)))
            s5(n1,-n2,-n3)  = CMPLX(REAL(s5(-n1,n2,n3)),-AIMAG(s5(-n1,n2,n3)))
            s5(-n1,n2,-n3)  = CMPLX(REAL(s5(n1,-n2,n3)),-AIMAG(s5(n1,-n2,n3)))
            s5(n1,n2,-n3)   = CMPLX(REAL(s5(-n1,-n2,n3)),-AIMAG(s5(-n1,-n2,n3)))
          ENDDO
        ENDDO
      ENDDO
     ENDIF
!****************
! take transform

        u1=(0.0,0.0)
        u2=(0.0,0.0)
        u3=(0.0,0.0)
        dens=(0.0,0.0)
        temp=(0.0,0.0)

        ndat = nx*ny*nz
        isign=1
        call fft(s1,u1,3,ndat,nx,ny,nz,isign)
        call fft(s2,u2,3,ndat,nx,ny,nz,isign)
        call fft(s3,u3,3,ndat,nx,ny,nz,isign)
        call fft(s4,dens,3,ndat,nx,ny,nz,isign)
        call fft(s5,temp,3,ndat,nx,ny,nz,isign)
     
        DO i=1,nx
         DO j=1,ny
           DO k=1,nz
              vel_u1(i,j,k)  =real(u1(i,j,k))
              vel_u2(i,j,k)  =real(u2(i,j,k))
              vel_u3(i,j,k)  =real(u3(i,j,k))
              vel_dens(i,j,k)=real(dens(i,j,k))
              vel_temp(i,j,k)=real(temp(i,j,k))
           ENDDO
         ENDDO
        ENDDO
!        stop

        nspec=1
        IF (nspec==1) THEN
          CALL compute_espec(u1,u2,u3,nx,ny,nz,nx,-0.001)
        ENDIF
!       stop

       RETURN
      END SUBROUTINE uv_init
!
!////////////////////////////////////////////////////
!
      SUBROUTINE compute_espec(u1,u2,u3,nx,ny,nz,nfour,tt)
!
      USE constants
      USE physics
!
      parameter(nknum=256)
      integer   :: nx,ny,nz,nfour,i,j,k,ii,jj,kk,isign,ndat,KMAX
      complex   :: u1(1:nfour,1:nfour,1:nfour)
      complex   :: u2(1:nfour,1:nfour,1:nfour)
      complex   :: u3(1:nfour,1:nfour,1:nfour)
      complex   :: s7(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! dill. Fourier coef
      complex   :: s6(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! vmag  Fourier coef
      complex   :: s3(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! u-vel  Fourier coef
      complex   :: s2(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! v-vel  Fourier coef
      complex   :: s1(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! w-vel  Fourier coef
      real      :: espec(nknum),espec1(nknum),espec2(nknum),espec3(nknum)
      real      :: especd(nknum),dspec(nknum)
      real      :: ur(3),ui(3),an1(3)
      double precision :: tt
! 
!     ----------------------------------------
!     tranfourrm the velocities to Fourier space
!     ----------------------------------------
!
      s6=(0.,0.)
      s3=(0.,0.)
      s2=(0.,0.)
      s1=(0.,0.)

      isign=-1
      ndat=nx*ny*nz
      ndim=3
      call fft(s1,u1,ndim,ndat,nx,ny,nz,isign)
      call fft(s2,u2,ndim,ndat,nx,ny,nz,isign)
      call fft(s3,u3,ndim,ndat,nx,ny,nz,isign)
!

      do k=-nz/2+1,nz/2
        do j=-ny/2+1,ny/2
          do i=-nx/2+1,nx/2
!          write(28,*) u3(i+nx/2,j+nx/2,k+nx/2) 
      an1(1)   = real(i)
      an1(2)   = real(j)
      an1(3)   = real(k)
      ur(1)    = REAL(s1(i,j,k))    
      ur(2)    = REAL(s2(i,j,k))    
      ur(3)    = REAL(s3(i,j,k))      
      ui(1)    = AIMAG(s1(i,j,k))    
      ui(2)    = AIMAG(s2(i,j,k))    
      ui(3)    = AIMAG(s3(i,j,k))    
      dum      = 0.
      wmag     = an1(1)**2 + an1(2)**2 + an1(3)**2
      do ii =1,3
        do jj=1,3
          dum  = dum+an1(ii)*an1(jj)*(ur(ii)*ur(jj)+ui(ii)*ui(jj))
        enddo
      enddo
      dum      = dum/wmag
      s7(i,j,k)=cmplx(dum,0.)
      s1(i,j,k)=cmplx(real(s1(i,j,k))**2 &
                +aimag(s1(i,j,k))**2,0.) 
      s2(i,j,k)=cmplx(real(s2(i,j,k))**2 &
                +aimag(s2(i,j,k))**2,0.)
      s3(i,j,k)=cmplx(real(s3(i,j,k))**2 &
                +aimag(s3(i,j,k))**2,0.)
      s6(i,j,k)=s1(i,j,k)+s2(i,j,k)+s3(i,j,k)
          enddo
        enddo
      enddo
!
!     -------------------------------------
!     compute the spectra in Fourrier space
!     -------------------------------------
!
      espec=0.
      espec1=0.
      espec2=0.
      espec3=0.
      especd=0.
      const=4.*pi/3.
      KMAX =nfour
      KMAX=INT(SQRT(2.)*REAL(nfour)/3.)                                  
      do 40 knum=1,KMAX
      r=float(knum)
      num=0
      do 30 k=-nz/2+1,nz/2
      do 30 j=-ny/2+1,ny/2
      do 30 i=-nx/2+1,nx/2
      w12=real(i)*real(i)
      w22=real(j)*real(j)
      w32=real(k)*real(k)
      wmag=sqrt(w12+w22+w32)
      if((wmag.le.r-.5).or.(wmag.gt.r+.5)) go to 30
20  continue
      num=num+1
      espec1(knum)=espec1(knum)+real(s1(i,j,k))
      espec2(knum)=espec2(knum)+real(s2(i,j,k))
      espec3(knum)=espec3(knum)+real(s3(i,j,k))
      espec(knum)=espec(knum)+real(s6(i,j,k))
      especd(knum)=especd(knum)+real(s7(i,j,k))
30  continue
      if(num.eq.0) go to 40
      corec=float(num)/(const*((r+0.5)**3-(r-0.5)**3))
      if(knum.eq.1) corec=float(num)/(const*(1.5)**3)
      espec(knum) =espec(knum)/corec
      espec1(knum)=espec1(knum)/corec
      espec2(knum)=espec2(knum)/corec
      espec3(knum)=espec3(knum)/corec
      especd(knum)=especd(knum)/corec
40  continue

      do 50 k=1,KMAX
      dspec(k)=2.0*espec(k)*REAL(k)**2 /re
      write(23,*) k,espec(k),espec1(k),espec2(k),espec3(k),especd(k),dspec(k)
50  continue
!
!     -------------------------------------
!     postprocess the spectra to obtain 
!     turbulence statistics
!     -------------------------------------
!
      sum1=0.0
      sum2=0.0
      sum3=0.0
      sum4=0.0
      do 70 k=1,KMAX
      sum1=sum1+espec(k)
      sum2=sum2+dspec(k)
      sum3=sum3+espec(k)/float(k)
      sum4=sum4+especd(k)
70    continue
      q=sqrt(2.0*sum1)
      qk=sum1
      qd=sum4
      uvar=q**2 /3.
      epsi=sum2
      xlength=3.0*pi*sum3/(2.0*q**2)
      xlam=sqrt(15.0*uvar/(re*epsi))
      eta=(1./(epsi*re**3))**(1./4.)
      taue=xlength/q*sqrt(3.)
      tauen=qk/epsi
      relamv=sqrt(uvar)*xlam*re
      ! write(90,*) tt,q,qd,qd/qk,epsi,xlength,xlam,eta,taue,tauen

      RETURN
      END SUBROUTINE compute_espec
!
!////////////////////////////////////////////////////
!
      SUBROUTINE compute_especa(s1,s2,s3,nx,ny,nz,nfour)
!
      USE constants
!
      parameter(nknum=256)
      complex   :: s1(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s2(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s3(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      complex   :: s6(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)
      real      :: espec(nknum),espec1(nknum),espec2(nknum),espec3(nknum)
! 
    ntest=1
    IF(ntest==1)THEN
      do k=-nz/2+1,nz/2
        do j=-ny/2+1,ny/2
          do i=-nx/2+1,nx/2
      s1(i,j,k)=cmplx(real(s1(i,j,k))**2 &
                +aimag(s1(i,j,k))**2,0.) 
      s2(i,j,k)=cmplx(real(s2(i,j,k))**2 &
                +aimag(s2(i,j,k))**2,0.)
      s3(i,j,k)=cmplx(real(s3(i,j,k))**2 &
                +aimag(s3(i,j,k))**2,0.)
      s6(i,j,k)=s1(i,j,k)+s2(i,j,k)+s3(i,j,k)
          enddo
        enddo
      enddo
   ENDIF
    ntest=0
    IF(ntest==1)THEN
      do k=-nz/2+1,nz/2
        do j=-ny/2+1,ny/2
          do i=-nx/2+1,nx/2
      s1(i,j,k)=cmplx(real(s1(i,j,k))**2 &
                ,0.) 
      s2(i,j,k)=cmplx(real(s2(i,j,k))**2 &
                ,0.)
      s3(i,j,k)=cmplx(real(s3(i,j,k))**2 &
                ,0.)
      s6(i,j,k)=s1(i,j,k)+s2(i,j,k)+s3(i,j,k)
          enddo
        enddo
      enddo
   ENDIF


      espec=0.
      espec1=0.
      espec2=0.
      espec3=0.
      const=4.*pi/3.
      KMAX =nfour
      do 40 knum=1,KMAX
      r=float(knum)
      num=0
      do 30 k=-nz/2+1,nz/2-1
      do 30 j=-ny/2+1,ny/2-1
      do 30 i=-nx/2+1,nx/2-1
      w12=real(i)*real(i)
      w22=real(j)*real(j)
      w32=real(k)*real(k)
      wmag=sqrt(w12+w22+w32)
      if((wmag.le.r-.5).or.(wmag.gt.r+.5)) go to 30
20  continue
      num=num+1
      espec1(knum)=espec1(knum)+real(s1(i,j,k))
      espec2(knum)=espec2(knum)+real(s2(i,j,k))
      espec3(knum)=espec3(knum)+real(s3(i,j,k))
      espec(knum)=espec(knum)+real(s6(i,j,k))
30  continue
      if(num.eq.0) go to 40
      corec=2.0*float(num)/(const*((r+0.5)**3-(r-0.5)**3))
      if(knum.eq.1) corec=2.0*float(num)/(const*(1.5)**3)
      espec(knum)=espec(knum)/corec
      espec1(knum)=espec1(knum)/corec
      espec2(knum)=espec2(knum)/corec
      espec3(knum)=espec3(knum)/corec
40  continue

      do 50 k=1,KMAX
      write(24,*) k,espec(k),espec1(k),espec2(k),espec3(k)
50  continue


      RETURN
      END SUBROUTINE compute_especa
!
!////////////////////////////////////////////////////
!
      SUBROUTINE compute_diss(d,id)
!
      USE domain_definition
      USE User_Data
      USE constants
      USE physics
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      TYPE (domain) :: d
!
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,5) :: Qx, Qy, Qz    
      DOUBLE PRECISION :: u,v,w,rho
      DOUBLE PRECISION :: ux,uy,uz,vx,vy,vz,wx,wy,wz
      INTEGER				      :: id
!
! compute derivatives
     DO nv = 1,4
         CALL Gauss_Deriv(d%gmetg,d%Qlgg(:,:,:,nv),d%Qglg(:,:,:,nv), &
                         d%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv), &
                         d%ncg,d%dmx,d%dmy,d%dmz)
     END DO
!
! compute dissipation with Qx for isotropic turbulence
!

!       DO i=1,d%ncg(1)
!         DO j=1,d%ncg(2)
!           DO k=1,d%ncg(3)
!             om1 = Qy(i,j,k,4)-Qz(i,j,k,3)
!             om2 = -Qz(i,j,k,2)+Qx(i,j,k,4)
!             om3 = Qy(i,j,k,2)-Qx(i,j,k,3)
!             om1 = om1**2
!             om2 = om2**2
!             om3 = om3**2
!             diss(id,i,j,k,1)= (om1+om2+om3)/re
!             diss(id,i,j,k,2)= 4.0d0*(Qx(i,j,k,2)+Qy(i,j,k,3)+Qz(i,j,k,4))**2/(3.0d0*re)
!           ENDDO    
!         ENDDO    
!       ENDDO   

      DO k = 1,d%ncg(3)
        DO j = 1,d%ncg(2)
         DO i = 1,d%ncg(1)
                                                                                                               
          rho = d%Q(i,j,k,1)
          u   = d%Q(i,j,k,2)/d%Q(i,j,k,1)
          v   = d%Q(i,j,k,3)/d%Q(i,j,k,1)
          w   = d%Q(i,j,k,4)/d%Q(i,j,k,1)
                                                                                                               
          ux = (Qx(i,j,k,2) - u*Qx(i,j,k,1))/rho
          uy = (Qy(i,j,k,2) - u*Qy(i,j,k,1))/rho
          uz = (Qz(i,j,k,2) - u*Qz(i,j,k,1))/rho
          vx = (Qx(i,j,k,3) - v*Qx(i,j,k,1))/rho
          vy = (Qy(i,j,k,3) - v*Qy(i,j,k,1))/rho
          vz = (Qz(i,j,k,3) - v*Qz(i,j,k,1))/rho
          wx = (Qx(i,j,k,4) - w*Qx(i,j,k,1))/rho
          wy = (Qy(i,j,k,4) - w*Qy(i,j,k,1))/rho
          wz = (Qz(i,j,k,4) - w*Qz(i,j,k,1))/rho

          om1 = wy - vz
          om2 = uz - wx
          om3 = vx - uy
          om1 = om1**2
          om2 = om2**2
          om3 = om3**2

          diss(id,i,j,k,1)= (om1+om2+om3)/re
          diss(id,i,j,k,2)= 4.0d0*(ux+vy+wz)**2/(3.0d0*re)
                                                                                                               
         END DO
       END DO
      END DO
!    
       RETURN
      END SUBROUTINE compute_diss
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE interp_matrix(corner,cxg,ncg,ngrid,id)
!
!......................................................................
!     compute matrix, that interpolates values from the chebyshev
!     element to a global fourrier grid
!......................................................................
!
      USE constants
      USE User_Data
!     
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!    
      INTEGER                         :: ncg,o
      DOUBLE PRECISION,DIMENSION(3,8) :: corner
      DOUBLE PRECISION,DIMENSION(nmax):: cxg
      DOUBLE PRECISION                :: polyn
! 
!
      npoin = 0
      DO i=1,nfour
        DO j=1,nfour
          DO k=1,nfour
           xfour =DBLE(i-1)*2.0d0*pi/DBLE(nfour-1)
           yfour =DBLE(j-1)*2.0d0*pi/DBLE(nfour-1)
           zfour =DBLE(k-1)*2.0d0*pi/DBLE(nfour-1)
          IF (xfour>=corner(1,1) .and. xfour<=corner(1,2) .and. &
              yfour>=corner(2,1) .and. yfour<=corner(2,3) .and. &
              zfour>=corner(3,1) .and. zfour<=corner(3,5) ) THEN
          xm=(xfour-corner(1,1))/(corner(1,2)-corner(1,1))
          ym=(yfour-corner(2,1))/(corner(2,3)-corner(2,1))
          zm=(zfour-corner(3,1))/(corner(3,5)-corner(3,1))
          npoin = npoin + 1
          DO m=1,ncg 
            DO n=1,ncg
               DO o=1,ncg
                 hx=polyn(m,xm,ncg,cxg(:))
                 hy=polyn(n,ym,ncg,cxg(:))
                 hm=polyn(o,zm,ncg,cxg(:))
                 ima(npoin,id,m)=hx
                 imb(npoin,id,n)=hy
                 imc(npoin,id,o)=hm
               ENDDO
            ENDDO
          ENDDO
          ENDIF
        ENDDO
       ENDDO
     ENDDO
!
      RETURN
      END SUBROUTINE interp_matrix
!
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE compute_source(x,y,z,t1,source)
!
!......................................................................
!     physical source term
!......................................................................
!
      USE size
      USE User_Data
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(neq) :: source
!
      
      source = 0.0d0

      RETURN
      END SUBROUTINE compute_source
!
!///////////////////////////////////////////////////////////////////////
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

        !   domlx = 4.0d0
        !   domly = 2.0d0
        !   domlz = 2.0d0
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
       USE mpi
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!
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
       USE mpi
!       
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!

      IF (myid == 0) THEN
         WRITE(rst_unit) time
         write(*,*) 'final time',time
      ENDIF
!
       RETURN
      END SUBROUTINE write_restart_case
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
!       IEVAP  = 0

!       SCF    = 0.7d0
!       SCO    = 0.7d0
!       DA     = 6.0d0
!       ZE     = 4.0d0
!       CE     = 40.0d0
!       RCOEF  = 1.0d0
!       YFINIT = 0.0d0
! !
! !     compute particle parameters (see formulation)
! !
!       PEF    = SCF*re
!       PD0    = dsqrt(18.0d0*TAUP0/(re*RHOP))
!       PM0    = pi*RHOP*PD0**3 /6.d0
!       PSI    = dble(nprt)*PM0/(2.d0*pi)**3
!        CRE    = re*(6.d0/(pi*RHOP))**(1.0d0/3.0d0)
!       CTAUP  = (RHOP**(1.0d0/3.0d0) *re/18.0d0)*(6.0d0/pi)**(2.0d0/3.0d0)
!       CYFPS  = gamma*A1/((gamma-1.)*TBOIL)
!       CF2    = 1.0d0/(3.0d0*pr*A2)
!       CF3    = A1/(3.0d0*SCF*A2)
!       CF4    = pi*dsqrt(18.0d0/RHOP)/(re**(3.0d0/2.0d0) *SCF)

! !   
!       nprtmax = 1102

      RETURN
   END SUBROUTINE set_part_par
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE part_initc(drop,nprt)

!
!    initialization of particle parameters
!
!       USE size
! !
!       USE particle_definition
! !
!       USE constants
!       USE part_par
!       USE physics
!       USE User_Data
! !
!       TYPE(particle)  :: drop
! !!    initialize particle field
! !
!       DO j=1,3
!          drop%Vp(j)    = drop%Vfp(j)
!       END DO
!          drop%Tp       = twall
!          drop%Mp       = PM0
!          drop%Tfp      = twall
!          drop%Rhofp    = 1.0d0
!          drop%Yffp     = 1.0d0
! !
!          drop%onoff     = 1
! !

      RETURN
   END SUBROUTINE part_initc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initScalars(d,eP,id)
    use size
    use particle_definition
    use domain_definition
    implicit none
    type(domain)                                :: d
    type(ePointers), dimension(ngp,npart)       :: eP
    double precision                            :: x,y,z
    integer                                     :: np, id
    double precision, parameter                 :: pie      = 3.14159265358979323846264338327d0
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

