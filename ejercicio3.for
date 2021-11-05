      program integral
      implicit none
      real*8 

      open(50, file='gl10.txt')
      do i=1, 9
        read(50, *) x(i), alfa(i)

      stop
      end program
!-------------------------------------------------------------------
      function f1(x)
      implicit none
      real*8 f1, x

      f1=5/4*(1-x**2/25)

      return
      end function
!-------------------------------------------------------------------
      function f2(x)
      implicit none
      real*8 f2, x, pi

      pi=4.d0*atan(1.d0)
      f2=5.d0/8.d0*(3*(1.d0-x**2/25.d0)+2*(1.d0-x**2/25.d0)**(1.d0/4.d0)
     &+(1-cos(4*pi*x/5)**2))

      return
      end function