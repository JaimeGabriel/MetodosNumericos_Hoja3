      Program derivadas
      Implicit none

      integer i
      real*8 h, x, f, g
      h=1
      x=10

      print*, "ACELERACIÓN EN X=10"

      print*, "Formula simetrica"
      do h=100, 10, -10
        g=(f(x+h/100.d0)-f(x))/(h/100.d0)
        print*, "Valor de la derivada para un paso de ", h/100, "-->", g
      end do
      print*,"-------------------------------------"

      print*, "Formula no simetrica"
      do h=100, 10, -10
        g=(f(x+h/100.d0)-f(x-h/100.d0))/(2*h/100.d0)
        print*, "Valor de la derivada para un paso de ", h/100, "-->", g
      end do

      print*,"-------------------------------------"

      print*, "Extrapolacion de Richardson"
      do h=100, 10, -10
        g=(-1/(6*(h/100.d0)))*(f(x+h/100.d0)-8*f(x+((h/100.d0)/2.d0))
     &-f(x-h/100.d0)+8*f(x-((h/100.d0)/2.d0)))
        print*, "Valor de la derivada para un paso de ", h/100, "-->", g
      end do


      stop
      end program

      function f(x)
      implicit none

      real*8 x, z0, vt, g, f

      z0=1000
      g=10
      vt=80

      f=-10*(2/dcosh(x/8))**2

      return
      end