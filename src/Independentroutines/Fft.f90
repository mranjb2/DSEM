!///////////////////////////////////////////////////////////////////////////////////////////////
!////////                                                                               ////////
!////////     Fft.f90                                                                   ////////
!////////                                                                               ////////
!////////     contains:                                                                 ////////
!////////                                                                               ////////
!////////           SUBROUTINE fft()                                                    ////////
!////////           SUBROUTINE fourn()                                                  ////////
!////////                                                                               ////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE fft(s1,u1,NDIM,NDAT,n1,n2,n3,isign)
!
!     computes the fast fourrier transform with single precision.
!     correction factors are used to account for the 
!     changing the range from (-nx/2+1)->(nx/2) to (0)->(nx-1)
!     
!    
!     date: 9/9/02
!     routines called: fourn
!     uses: constants
!     applicability: all
!
      USE constants
!
      INTEGER           :: NDAT,NDIM,n1,n2,n3
      INTEGER           :: i,idum,isign,j,k,l,nn(NDIM)
      REAL      :: data1(2*NDAT),ran1,dum1,dum2,dum3
      complex   :: s1(-n1/2+1:n1/2,-n2/2+1:n2/2,-n3/2+1:n3/2)
      complex   :: u1(1:n1,1:n2,1:n3)
      complex   :: corr1,corr2,corr3,tr1,tr2,tr3
!
      nn(1)=n1
      nn(2)=n2
      nn(3)=n3
!
    ! write(*,*) 'chick isign',isign
    IF (isign==1) THEN  ! discrete transform

      do 14 i=1,nn(3)
        do 13 j=1,nn(2)
          do 12 k=1,nn(1)
            l        = k+(j-1)*nn(1)+(i-1)*nn(2)*nn(1)
            l        = 2*l-1
            ii       = i-nn(3)/2
            jj       = j-nn(2)/2
            kk       = k-nn(1)/2
            data1(l) = REAL(s1(kk,jj,ii)) 
            l        = l+1
            data1(l) = AIMAG(s1(kk,jj,ii)) 
12        continue
13      continue
14    continue
!
      call fourn(data1,nn,NDIM,isign)
!
      do 15  i=1,nn(3)
        dum3 = 2.0d0*pi*REAL((1-n3/2)*(i-1))/REAL(n3)
        corr3 =CMPLX(0.0d0,dum3)
        do 16 j=1,nn(2)
          dum2 = 2.0d0*pi*REAL((1-n2/2)*(j-1))/REAL(n2)
          corr2 =CMPLX(0.0d0,dum2)
          do 17 k=1,nn(1)
            dum1       = 2.0d0*pi*REAL((1-n1/2)*(k-1))/REAL(n1)
            corr1      = CMPLX(0.0d0,dum1)
            l          = k+(j-1)*nn(1)+(i-1)*nn(2)*nn(1)
            l          = 2*l-1
            tr1        = CMPLX(data1(l),data1(l+1))
            data1(l)   = REAL(tr1*exp(corr1+corr2+corr3))
            data1(l+1) = AIMAG(tr1*exp(corr1+corr2+corr3))
            u1(k,j,i)  = CMPLX(data1(l),data1(l+1))
17       continue
16      continue
15    continue

    ELSEIF (isign==-1) THEN  ! inverse transfrom

      do 18 i=1,nn(3)
        dum3 = -2.0d0*pi*REAL((1-n3/2)*(i-1))/REAL(n3)
        corr3 =CMPLX(0.0d0,dum3)
        do 19 j=1,nn(2)
          dum2 = -2.0d0*pi*REAL((1-n2/2)*(j-1))/REAL(n2)
          corr2 =CMPLX(0.0d0,dum2)
          do 20 k=1,nn(1)
            dum1       = -2.0d0*pi*REAL((1-n1/2)*(k-1))/REAL(n1)
            corr1      = CMPLX(0.0d0,dum1)
            l          = k+(j-1)*nn(1)+(i-1)*nn(2)*nn(1)
            l          = 2*l-1
            tr1        = u1(k,j,i)
            data1(l)   = REAL(tr1*exp(corr1+corr2+corr3))
            data1(l+1) = AIMAG(tr1*exp(corr1+corr2+corr3))
20        continue
19      continue
18    continue

      call fourn(data1,nn,NDIM,isign)

      do 21 i=1,nn(3)
        do 22 j=1,nn(2)
          do 23 k=1,nn(1)
            l            = k+(j-1)*nn(1)+(i-1)*nn(2)*nn(1)
            l            = 2*l-1
            ii           = i-nn(3)/2
            jj           = j-nn(2)/2
            kk           = k-nn(1)/2
            data1(l)     = data1(l)
            data1(l+1)   = data1(l+1)
            s1(kk,jj,ii) = CMPLX(data1(l)/NDAT,data1(l+1)/NDAT)
23       continue
22      continue
21    continue
 
    ENDIF
!
      RETURN
      END
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE fourn(data,nn,ndim,isign)
!
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,  &
      k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr

                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
