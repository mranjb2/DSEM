!> @file
!> This file contains  matrix-vector product routines.
!!
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!                                                                       
!                          X-MULTIPLY ROUTINES                          
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!> Matrix multiplication routines                                                                    
      subroutine mxm(a,n1,b,n2,c,n3) 
!---------------------------------------------------------------------- 
!                                                                       
!     matrix-vector product routine.                                    
!     note: use assembly coded routine if available.                    
!                                                                       
!--------------------------------------------------------------------- 
      USE size
!
       implicit double precision (a-h,o-z) 
       dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      if (n2.le.8) then 
         if (n2.eq.1) then 
            call mxm1(a,n1,b,n2,c,n3) 
         elseif (n2.eq.2) then 
            call mxm2(a,n1,b,n2,c,n3) 
         elseif (n2.eq.3) then 
            call mxm3(a,n1,b,n2,c,n3) 
         elseif (n2.eq.4) then 
            call mxm4(a,n1,b,n2,c,n3) 
         elseif (n2.eq.5) then 
            call mxm5(a,n1,b,n2,c,n3) 
         elseif (n2.eq.6) then 
            call mxm6(a,n1,b,n2,c,n3) 
         elseif (n2.eq.7) then 
            call mxm7(a,n1,b,n2,c,n3) 
         else 
            call mxm8(a,n1,b,n2,c,n3) 
         endif 
      elseif (n2.le.16) then 
         if (n2.eq.9) then 
            call mxm9(a,n1,b,n2,c,n3) 
         elseif (n2.eq.10) then 
            call mxm10(a,n1,b,n2,c,n3) 
         elseif (n2.eq.11) then 
            call mxm11(a,n1,b,n2,c,n3) 
         elseif (n2.eq.12) then 
            call mxm12(a,n1,b,n2,c,n3) 
         elseif (n2.eq.13) then 
            call mxm13(a,n1,b,n2,c,n3) 
         elseif (n2.eq.14) then 
            call mxm14(a,n1,b,n2,c,n3) 
         elseif (n2.eq.15) then 
            call mxm15(a,n1,b,n2,c,n3) 
         else 
            call mxm16(a,n1,b,n2,c,n3) 
         endif 
      elseif ( n2 .le. 24 ) then 
          if ( n2 .eq. 17 ) then 
            call mxm17(a,n1,b,n2,c,n3) 
         elseif (n2.eq.18) then 
            call mxm18(a,n1,b,n2,c,n3) 
         elseif (n2.eq.19) then 
            call mxm19(a,n1,b,n2,c,n3) 
         elseif (n2.eq.20) then 
            call mxm20(a,n1,b,n2,c,n3) 
         elseif (n2.eq.21) then 
            call mxm21(a,n1,b,n2,c,n3) 
         elseif (n2.eq.22) then 
            call mxm22(a,n1,b,n2,c,n3) 
         elseif (n2.eq.23) then 
            call mxm23(a,n1,b,n2,c,n3) 
         elseif (n2.eq.24) then 
            call mxm24(a,n1,b,n2,c,n3) 
         endif 
      elseif ( n2 .le. 32 ) then 
         if ( n2 .eq. 25 ) then 
            call mxm25(a,n1,b,n2,c,n3) 
         elseif (n2.eq.26) then 
            call mxm26(a,n1,b,n2,c,n3) 
         elseif (n2.eq.27) then 
            call mxm27(a,n1,b,n2,c,n3) 
         elseif (n2.eq.28) then 
            call mxm28(a,n1,b,n2,c,n3) 
         elseif (n2.eq.29) then 
            call mxm29(a,n1,b,n2,c,n3) 
         elseif (n2.eq.30) then 
            call mxm30(a,n1,b,n2,c,n3) 
         elseif (n2.eq.31) then 
            call mxm31(a,n1,b,n2,c,n3) 
         elseif (n2.eq.32) then 
            call mxm32(a,n1,b,n2,c,n3) 
         endif 
      else 
         n0=n1*n3 
         do 10 i=1,n0 
            c(i,1)=0. 
   10    continue 
         do 100 j=1,n3 
            do 100 k=1,n2 
               bb=b(k,j) 
               do 100 i=1,n1 
                  c(i,j)=c(i,j)+a(i,k)*bb 
  100    continue 
      endif 
                                                                        
      return 
      END                                           
!                                                                       
      subroutine mxm1(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,1),b(1,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j) 
         enddo 
      enddo 
      return 
      END                                           
!                                                                       
      subroutine mxm2(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,2),b(2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      
         enddo 
      enddo 

      return 
      END                                           
      subroutine mxm3(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,3),b(3,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm4(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,4),b(4,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm5(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,5),b(5,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm6(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,6),b(6,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm7(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,7),b(7,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm8(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,8),b(8,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm9(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,9),b(9,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm10(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,10),b(10,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm11(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,11),b(11,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm12(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,12),b(12,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm13(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,13),b(13,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm14(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,14),b(14,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm15(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,15),b(15,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm16(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm17(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm18(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    
         enddo 
      enddo 
      return 
      END                                           
         subroutine mxm19(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    
                enddo 
           enddo 
           return 
      END                                           
      subroutine mxm20(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm21(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm22(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm23(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm24(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm25(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm26(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm27(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm28(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    &
     &             + a(i,28)*b(28,j)                                    
                                                                        
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm29(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    &
     &             + a(i,28)*b(28,j)                                    &
     &             + a(i,29)*b(29,j)                                    
                                                                        
                                                                        
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm30(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    &
     &             + a(i,28)*b(28,j)                                    &
     &             + a(i,29)*b(29,j)                                    &
     &             + a(i,30)*b(30,j)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm31(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    &
     &             + a(i,28)*b(28,j)                                    &
     &             + a(i,29)*b(29,j)                                    &
     &             + a(i,30)*b(30,j)                                    &
     &             + a(i,31)*b(31,j)                                    
          enddo 
      enddo 
!                                                                       
      return 
      END                                           
      subroutine mxm32(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(i,1)*b(1,j)                                      &
     &             + a(i,2)*b(2,j)                                      &
     &             + a(i,3)*b(3,j)                                      &
     &             + a(i,4)*b(4,j)                                      &
     &             + a(i,5)*b(5,j)                                      &
     &             + a(i,6)*b(6,j)                                      &
     &             + a(i,7)*b(7,j)                                      &
     &             + a(i,8)*b(8,j)                                      &
     &             + a(i,9)*b(9,j)                                      &
     &             + a(i,10)*b(10,j)                                    &
     &             + a(i,11)*b(11,j)                                    &
     &             + a(i,12)*b(12,j)                                    &
     &             + a(i,13)*b(13,j)                                    &
     &             + a(i,14)*b(14,j)                                    &
     &             + a(i,15)*b(15,j)                                    &
     &             + a(i,16)*b(16,j)                                    &
     &             + a(i,17)*b(17,j)                                    &
     &             + a(i,18)*b(18,j)                                    &
     &             + a(i,19)*b(19,j)                                    &
     &             + a(i,20)*b(20,j)                                    &
     &             + a(i,21)*b(21,j)                                    &
     &             + a(i,22)*b(22,j)                                    &
     &             + a(i,23)*b(23,j)                                    &
     &             + a(i,24)*b(24,j)                                    &
     &             + a(i,25)*b(25,j)                                    &
     &             + a(i,26)*b(26,j)                                    &
     &             + a(i,27)*b(27,j)                                    &
     &             + a(i,28)*b(28,j)                                    &
     &             + a(i,29)*b(29,j)                                    &
     &             + a(i,30)*b(30,j)                                    &
     &             + a(i,31)*b(31,j)                                    &
     &             + a(i,32)*b(32,j)                                    
                                                                        
          enddo 
      enddo 
!                                                                       
      return 
      END                                           
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!                                                                       
!                MATRIX-VECTOR MULTIPLY ROUTINES                        
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!                                                                       
      subroutine mxv(lda,a,n1,b,n2,c) 
!---------------------------------------------------------------------- 
!                                                                       
!     matrix-vector product routine.                                    
!     note: use assembly coded routine if available.                    
!                                                                       
!---------------------------------------------------------------------  
      USE size
!
      implicit double precision (a-h,o-z) 
       dimension a(lda,*),b(*),c(*) 
!                                                                       
      if (n2.le.8) then 
         if (n2.eq.1) then 
            call mxv1(lda,a,n1,b,c) 
         elseif (n2.eq.2) then 
            call mxv2(lda,a,n1,b,c) 
         elseif (n2.eq.3) then 
            call mxv3(lda,a,n1,b,c) 
         elseif (n2.eq.4) then 
            call mxv4(lda,a,n1,b,c) 
         elseif (n2.eq.5) then 
            call mxv5(lda,a,n1,b,c) 
         elseif (n2.eq.6) then 
            call mxv6(lda,a,n1,b,c) 
         elseif (n2.eq.7) then 
            call mxv7(lda,a,n1,b,c) 
         else 
            call mxv8(lda,a,n1,b,c) 
         endif 
      elseif (n2.le.16) then 
         if (n2.eq.9) then 
            call mxv9(lda,a,n1,b,c) 
         elseif (n2.eq.10) then 
            call mxv10(lda,a,n1,b,c) 
         elseif (n2.eq.11) then 
            call mxv11(lda,a,n1,b,c) 
         elseif (n2.eq.12) then 
            call mxv12(lda,a,n1,b,c) 
         elseif (n2.eq.13) then 
            call mxv13(lda,a,n1,b,c) 
         elseif (n2.eq.14) then 
            call mxv14(lda,a,n1,b,c) 
         elseif (n2.eq.15) then 
            call mxv15(lda,a,n1,b,c) 
         else 
            call mxv16(lda,a,n1,b,c) 
         endif 
      elseif ( n2 .le. 24 ) then 
          if ( n2 .eq. 17 ) then 
            call mxv17(lda,a,n1,b,c) 
         elseif (n2.eq.18) then 
            call mxv18(lda,a,n1,b,c) 
         elseif (n2.eq.19) then 
            call mxv19(lda,a,n1,b,c) 
         elseif (n2.eq.20) then 
            call mxv20(lda,a,n1,b,c) 
         elseif (n2.eq.21) then 
            call mxv21(lda,a,n1,b,c) 
         elseif (n2.eq.22) then 
            call mxv22(lda,a,n1,b,c) 
         elseif (n2.eq.23) then 
            call mxv23(lda,a,n1,b,c) 
         elseif (n2.eq.24) then 
            call mxv24(lda,a,n1,b,c) 
         endif 
      elseif ( n2 .le. 32 ) then 
         if ( n2 .eq. 25 ) then 
            call mxv25(lda,a,n1,b,c) 
         elseif (n2.eq.26) then 
            call mxv26(lda,a,n1,b,c) 
         elseif (n2.eq.27) then 
            call mxv27(lda,a,n1,b,c) 
         elseif (n2.eq.28) then 
            call mxv28(lda,a,n1,b,c) 
         elseif (n2.eq.29) then 
            call mxv29(lda,a,n1,b,c) 
         elseif (n2.eq.30) then 
            call mxv30(lda,a,n1,b,c) 
         elseif (n2.eq.31) then 
            call mxv31(lda,a,n1,b,c) 
         elseif (n2.eq.32) then 
            call mxv32(lda,a,n1,b,c) 
         endif 
      else 
         do 10 i = 1,n1 
            c(i)=  0. 
   10    continue 
         do 100 k=1,n2 
            bb=b(k) 
            do 100 i=1,n1 
               c(i)=c(i)+a(i,k)*bb 
  100    continue 
      endif 
                                                                        
      return 
      END                                           
!                                                                       
      subroutine mxv1(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
         c(i) = a(i,1)*b(1) 
      enddo 
      return 
      END                                           
!                                                                       
      subroutine mxv2(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
         c(i) = a(i,1)*b(1) + a(i,2)*b(2) 
      enddo 
      return 
      END                                           
      subroutine mxv3(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
         c(i) = a(i,1)*b(1)                                             &
     &          + a(i,2)*b(2)                                           &
     &          + a(i,3)*b(3)                                           
      enddo 
      return 
      END                                           
      subroutine mxv4(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
         c(i) = a(i,1)*b(1)                                             &
     &          + a(i,2)*b(2)                                           &
     &          + a(i,3)*b(3)                                           &
     &          + a(i,4)*b(4)                                           
      enddo 
      return 
      END                                           
      subroutine mxv5(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
         c(i) = a(i,1)*b(1)                                             &
     &          + a(i,2)*b(2)                                           &
     &          + a(i,3)*b(3)                                           &
     &          + a(i,4)*b(4)                                           &
     &          + a(i,5)*b(5)                                           
      enddo 
      return 
      END                                           
      subroutine mxv6(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        
      enddo 
      return 
      END                                           
      subroutine mxv7(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        
      enddo 
      return 
      END                                           
      subroutine mxv8(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        
      enddo 
      return 
      END                                           
      subroutine mxv9(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
      do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        
      enddo 
      return 
      END                                           
      subroutine mxv10(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      
         enddo 
      return 
      END                                           
      subroutine mxv11(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      
         enddo 
      return 
      END                                           
      subroutine mxv12(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      
         enddo 
      return 
      END                                           
      subroutine mxv13(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      
         enddo 
      return 
      END                                           
      subroutine mxv14(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      
         enddo 
      return 
      END                                           
      subroutine mxv15(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      
         enddo 
      return 
      END                                           
      subroutine mxv16(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      
         enddo 
      return 
      END                                           
      subroutine mxv17(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      
         enddo 
      return 
      END                                           
      subroutine mxv18(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      
         enddo 
      return 
      END                                           
         subroutine mxv19(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      
                enddo 
      return 
      END                                           
      subroutine mxv20(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      
             enddo 
      return 
      END                                           
      subroutine mxv21(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      
             enddo 
      return 
      END                                           
      subroutine mxv22(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv23(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv24(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv25(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv26(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      
             enddo 
      return 
      END                                           
      subroutine mxv27(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv28(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      &
     &             + a(i,28)*b(28)                                      
                                                                        
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv29(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      &
     &             + a(i,28)*b(28)                                      &
     &             + a(i,29)*b(29)                                      
                                                                        
                                                                        
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv30(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      &
     &             + a(i,28)*b(28)                                      &
     &             + a(i,29)*b(29)                                      &
     &             + a(i,30)*b(30)                                      
                                                                        
             enddo 
      return 
      END                                           
      subroutine mxv31(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      &
     &             + a(i,28)*b(28)                                      &
     &             + a(i,29)*b(29)                                      &
     &             + a(i,30)*b(30)                                      &
     &             + a(i,31)*b(31)                                      
          enddo 
!                                                                       
      return 
      END                                           
      subroutine mxv32(lda,a,n1,b,c) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
      dimension a(lda,*),b(*),c(*) 
!                                                                       
         do i=1,n1 
            c(i) = a(i,1)*b(1)                                          &
     &             + a(i,2)*b(2)                                        &
     &             + a(i,3)*b(3)                                        &
     &             + a(i,4)*b(4)                                        &
     &             + a(i,5)*b(5)                                        &
     &             + a(i,6)*b(6)                                        &
     &             + a(i,7)*b(7)                                        &
     &             + a(i,8)*b(8)                                        &
     &             + a(i,9)*b(9)                                        &
     &             + a(i,10)*b(10)                                      &
     &             + a(i,11)*b(11)                                      &
     &             + a(i,12)*b(12)                                      &
     &             + a(i,13)*b(13)                                      &
     &             + a(i,14)*b(14)                                      &
     &             + a(i,15)*b(15)                                      &
     &             + a(i,16)*b(16)                                      &
     &             + a(i,17)*b(17)                                      &
     &             + a(i,18)*b(18)                                      &
     &             + a(i,19)*b(19)                                      &
     &             + a(i,20)*b(20)                                      &
     &             + a(i,21)*b(21)                                      &
     &             + a(i,22)*b(22)                                      &
     &             + a(i,23)*b(23)                                      &
     &             + a(i,24)*b(24)                                      &
     &             + a(i,25)*b(25)                                      &
     &             + a(i,26)*b(26)                                      &
     &             + a(i,27)*b(27)                                      &
     &             + a(i,28)*b(28)                                      &
     &             + a(i,29)*b(29)                                      &
     &             + a(i,30)*b(30)                                      &
     &             + a(i,31)*b(31)                                      &
     &             + a(i,32)*b(32)                                      
                                                                        
          enddo 
!                                                                       
      return 
      END                                           
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!                                                                       
!                          Y-MULTIPLY ROUTINES                          
!                                                                       
!////////////////////////////////////////////////////////////////////// 
!                                                                       
      subroutine mxmy(a,n1,b,n2,c,n3) 
!---------------------------------------------------------------------- 
!                                                                       
!     matrix-vector product routine.                                    
!     note: use assembly coded routine if available.                    
!                                                                       
!---------------------------------------------------------------------  
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      if (n2.le.8) then 
         if (n2.eq.1) then 
            call mxm1y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.2) then 
            call mxm2y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.3) then 
            call mxm3y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.4) then 
            call mxm4y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.5) then 
            call mxm5y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.6) then 
            call mxm6y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.7) then 
            call mxm7y(a,n1,b,n2,c,n3) 
         else 
            call mxm8y(a,n1,b,n2,c,n3) 
         endif 
      elseif (n2.le.16) then 
         if (n2.eq.9) then 
            call mxm9y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.10) then 
            call mxm10y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.11) then 
            call mxm11y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.12) then 
            call mxm12y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.13) then 
            call mxm13y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.14) then 
            call mxm14y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.15) then 
            call mxm15y(a,n1,b,n2,c,n3) 
         else 
            call mxm16y(a,n1,b,n2,c,n3) 
         endif 
      elseif ( n2 .le. 24 ) then 
          if ( n2 .eq. 17 ) then 
            call mxm17y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.18) then 
            call mxm18y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.19) then 
            call mxm19y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.20) then 
            call mxm20y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.21) then 
            call mxm21y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.22) then 
            call mxm22y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.23) then 
            call mxm23y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.24) then 
            call mxm24y(a,n1,b,n2,c,n3) 
         endif 
      elseif ( n2 .le. 32 ) then 
         if ( n2 .eq. 25 ) then 
            call mxm25y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.26) then 
            call mxm26y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.27) then 
            call mxm27y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.28) then 
            call mxm28y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.29) then 
            call mxm29y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.30) then 
            call mxm30y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.31) then 
            call mxm31y(a,n1,b,n2,c,n3) 
         elseif (n2.eq.32) then 
            call mxm32y(a,n1,b,n2,c,n3) 
         endif 
      else 
         n0=n1*n3 
         do 10 i=1,n0 
            c(i,1)=0. 
   10    continue 
         do 100 j=1,n3 
            do 100 k=1,n2 
               bb=b(k,j) 
               do 100 i=1,n1 
                  c(i,j)=c(i,j)+a(i,k)*bb 
  100    continue 
      endif 
                                                                        
      return 
      END                                           
!                                                                       
      subroutine mxm1y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1) 
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm2y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm3y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm4y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm5y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm6y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm7y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm8y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm9y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm10y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm11y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm12y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm13y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm14y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm15y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    
         enddo 
      enddo 
      return 
      END                                           
      subroutine mxm16y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    
         enddo 
      enddo 
      return 
      END                                           
                                                                        
      subroutine mxm17y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    
         enddo 
      enddo 
      return 
      END                                           
         subroutine mxm18y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    
         enddo 
      enddo 
      return 
      END                                           
         subroutine mxm19y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    
                enddo 
           enddo 
           return 
      END                                           
      subroutine mxm20y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm21y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm22y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm23y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm24y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm25y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm26y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm27y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm28y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    &
     &             + a(j,28)*b(i,28)                                    
                                                                        
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm29y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    &
     &             + a(j,28)*b(i,28)                                    &
     &             + a(j,29)*b(i,29)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm30y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    &
     &             + a(j,28)*b(i,28)                                    &
     &             + a(j,29)*b(i,29)                                    &
     &             + a(j,30)*b(i,30)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm31y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    &
     &             + a(j,28)*b(i,28)                                    &
     &             + a(j,29)*b(i,29)                                    &
     &             + a(j,30)*b(i,30)                                    &
     &             + a(j,31)*b(i,31)                                    
             enddo 
      enddo 
      return 
      END                                           
      subroutine mxm32y(a,n1,b,n2,c,n3) 
!                                                                       
      USE size
!
      implicit double precision (a-h,o-z) 
       double precision  a(n1,n2),b(n2,n3),c(n1,n3) 
                                                                        
!                                                                       
      do j=1,n3 
         do i=1,n1 
            c(i,j) = a(j,1)*b(i,1)                                      &
     &             + a(j,2)*b(i,2)                                      &
     &             + a(j,3)*b(i,3)                                      &
     &             + a(j,4)*b(i,4)                                      &
     &             + a(j,5)*b(i,5)                                      &
     &             + a(j,6)*b(i,6)                                      &
     &             + a(j,7)*b(i,7)                                      &
     &             + a(j,8)*b(i,8)                                      &
     &             + a(j,9)*b(i,9)                                      &
     &             + a(j,10)*b(i,10)                                    &
     &             + a(j,11)*b(i,11)                                    &
     &             + a(j,12)*b(i,12)                                    &
     &             + a(j,13)*b(i,13)                                    &
     &             + a(j,14)*b(i,14)                                    &
     &             + a(j,15)*b(i,15)                                    &
     &             + a(j,16)*b(i,16)                                    &
     &             + a(j,17)*b(i,17)                                    &
     &             + a(j,18)*b(i,18)                                    &
     &             + a(j,19)*b(i,19)                                    &
     &             + a(j,20)*b(i,20)                                    &
     &             + a(j,21)*b(i,21)                                    &
     &             + a(j,22)*b(i,22)                                    &
     &             + a(j,23)*b(i,23)                                    &
     &             + a(j,24)*b(i,24)                                    &
     &             + a(j,25)*b(i,25)                                    &
     &             + a(j,26)*b(i,26)                                    &
     &             + a(j,27)*b(i,27)                                    &
     &             + a(j,28)*b(i,28)                                    &
     &             + a(j,29)*b(i,29)                                    &
     &             + a(j,30)*b(i,30)                                    &
     &             + a(j,31)*b(i,31)                                    &
     &             + a(j,32)*b(i,32)                                    
                                                                        
             enddo 
      enddo 
      return 
      END                                           
