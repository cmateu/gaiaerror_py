c----------------------------------------------------------------------
c  This is 'const_ast.h'
c  Astronomical and physical constants (SI units)
c  ep0 is the reference epoch for the astrometric parameters,
c  reckoned from J2000.0
c
      real*8 kt,Gte
      PARAMETER (pc = 3.086d18, !en cm
     +           Gte = 6.670d-8,  ! en  dines cm2/g2 (cm3/(s2 g)) (Bowers)
     +           time=0.10225d0) ! canvi de km/s/kpc a 1/(100 Myr)

      parameter (kt=4.7404705d0)
      parameter (a1=(62.87124882d0*deg))
      parameter (a2=(32.93680516d0*deg))
      parameter (Ag=(282.8594813d0*deg))
      parameter (Dg=(27.12825111d0*deg))



