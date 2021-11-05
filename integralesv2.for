      program integrales
      implicit none
      
      real*8 pi
      real*8 a, b

C$$$$$$       print*, "Dame el punto inicial: "
C$$$$$$       read*, a
C$$$$$$       print*, "Dame el punto final: "
C$$$$$$       read*, b


      pi=4*atan(1.d0)

      a=0.d0
      b=2*pi

      call punto_medio(pi, a, b)
      call trapecio(pi, a, b)
      call sinson(pi, a, b)
      call newton_cotes(pi, a, b)
      call Gauss_Legendre_2puntos(pi, a, b)
      call Gauss_Legendre_10puntos(pi, a, b)
      call Regla_trapezoidal(pi, a, b)
      stop
      end program
!----------------------------------------------------------------------
      subroutine punto_medio(pi, a, b)
      implicit none
      real*8 pi, h, integral, f, a, b

      h=(b-a)
      integral=h*f(a+h/2)
      print*, "Punto medio: ", integral
      
      return
      end subroutine

      subroutine trapecio(pi, a, b)
      implicit none
      real*8 pi, h, integral, f, a, b, fa

      h=b-a

      if (a.lt.0.001)then
        fa=1
        integral=h/2*(fa+f(b))
      else
        integral=h/2*(f(a)+f(b))
      end if
      print*, "Trapecio: ", integral
      
      return
      end subroutine
!----------------------------------------------------------------------
      subroutine sinson(pi, a, b)
      implicit none
      real*8 pi, h, integral, f, a, b, fa, x1, h2

      h=b-a
      x1=a+h/2
      h2=b-x1

      if(a.lt.0.001)then
        fa=1
        integral=h2/3*(fa+4*f(x1)+f(b))
      else
        integral=h2/3*(f(a)+4*f(x1)+f(b))
      end if
      print*, "Simpson: ", integral
      
      return
      end subroutine
!----------------------------------------------------------------------
      subroutine newton_cotes(pi, a, b)
      implicit none
      real*8 pi, h, integral, f, a, b, fa, x1, h2, a0, a1, a2, a3, a4
      real*8 x0, x1, x2, x3, x4

      h = (b - a)/4
      x0 = a
      x1 = a + h
      x2 = a + 2*h
      x3 = a + 3*h
      x4 = a + 4*h

      a0=7*pi/45
      a1=32*pi/45
      a2=4*pi/15
      a3=32*pi/45
      a4=7*pi/45

      integral=a0*1+a1*f(x1)+a2*f(x2)+a3*f(x3)+a4*f(x4)
      print*, "Newton Cotes: ", integral

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine Gauss_Legendre_2puntos(pi, a, b)
      implicit none
      integer i, n
      real*8 alfa, x, pi, a, b, integral, suma, f, sum_i

      dimension alfa(0:100), x(0:100)
      
      n=1
      open(40, file='gl02.txt')
      do i=0, n
        read(40, *) x(i), alfa(i)
      end do

      suma=0.d0
      
      do i=0, n
        sum_i=alfa(i)*f((b+a+(b-a)*x(i))/2)
        suma=suma+sum_i
      end do

      integral=((b-a)/2)*suma

      print*, "Gauss-Legendre (n=1)", integral

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine Gauss_Legendre_10puntos(pi, a, b)
      implicit none
      integer i, n
      real*8 alfa, x, pi, a, b, integral, suma, f, sum_i

      dimension alfa(0:100), x(0:100)
      
      n=9
      open(41, file='gl10.txt')
      do i=0, n
        read(41, *) x(i), alfa(i)
      end do

      suma=0.d0
      
      do i=0, n
        sum_i=alfa(i)*f((b+a+(b-a)*x(i))/2)
        suma=suma+sum_i
      end do

      integral=((b-a)/2)*suma

      print*, "Gauss-Legendre (n=9)", integral

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine Regla_trapezoidal(pi, a, b)
      implicit none
      integer i, n, nn
      real*8 x, pi, a, b, h, suma, integral, f

      dimension x(0:100)
      
      nn=30.d0

      h=(b-a)/nn
      x(0)=a
      do i=1, nn
        x(i)=x(i-1)+h
      end do
      suma=0.d0
      do i=1, nn-1
        suma=suma+f(x(i))
      end do

      integral=((b-a)/nn)*(((1)+f(x(nn)))/2+suma)

      print*, "Regla Trapezoidal:", integral     


      return
      end subroutine
!----------------------------------------------------------------------


      function f(x)
      implicit none
      real*8 f, x

      f=sin(x)/x

      return
      end function
      