c-----------------------------------------------------------------------
c  This is 'const_math.h'
c  Mathematical constants and general units
c
      REAL*8 pi, pihalf, twopi, fourpi, deg, arcsec, mas, uas, 
     +       yr, day, sec, km, unit_se
      PARAMETER (pi = 3.14159265358979324d0,
     +           pihalf = pi/2d0,
     +           twopi = 2d0*pi,
     +           fourpi = 4d0*pi,
     +           deg = pi/180d0, 
     +           arcsec = deg/3600d0,
     +           mas = 1d-3*arcsec,
     +           uas = 1d-6*arcsec,
     +           yr = 365.25d0,
     +           day = 1d0,
     +           sec = 1d0/86400d0,
     +           km = 1000d0,
     +           unit_se = mas)
