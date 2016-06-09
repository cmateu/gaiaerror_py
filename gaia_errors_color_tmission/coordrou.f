


      subroutine Carte_to_UVWH(VX,VY,VZ,U,V,W)
C *******************************************************************
C CHANGES FROM (VX,VY,VZ) CARTESIAN GALACTOCENTRIC VELOCITY VECTOR TO
C (U,V,W) HELIOCENTRIC VELOCITY VECTOR 
c units: km/s
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

C ROTATION OF VX,VY OF ANGLE 0. TO PUT X-AXIS ALIGNED WITH SUN-GC 
C LINE
      ANGLE=0.!PI/2.D0
      CALL rot_horari(angle,VX,VY,VX2,VY2)
      VZ2=VZ
C SUBSTRACT THE DIFFERENTIAL ROTATION OF THE SUN
      CALL dV2(R0,0.d0,W)
      U2=VX2
      V2=VY2-W*R0
      W2=VZ2
C SUBSTRACT THE MOTION OF THE SUN WITH RESPECT TO THE LSR
      U=U2-USOL
      V=V2-VSOL
      W=W2-WSOL
      RETURN
      END


C  --------------------------------------------------------------------
c     SUBROUTINE carte_to_equatorial(X,Y,Z,ALPHA,DELTA)
c     Input cartesian coordinates: (X,Y,Z) in kpc
c     Output: equatorial coordinates: (alpha,delta) in radiants
c     
c---------------------------------------------------------------------
      SUBROUTINE carte_to_equatorial(X,Y,Z,ALPHA,DELTA)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

c XYZ, centered on the Sun
       sX=X+R0
       sY=Y!-R0
       sZ=Z

c R,l,b heliocentric
       gr=dsqrt(sx*sx+sy*sy+sz*sz)
       gb=asin(sz/gr)
cc x>0, y<0:
c	if ((sx.gt.0).and.(sy.lt.0)) then
c       gl=atan(dabs(sx/sy))
c       endif
cc x>0, y>0:
c	if ((sx.gt.0).and.(sy.gt.0)) then
c       gl=atan(dabs(sy/sx))+pi/2.d0
c       endif
cc x<0, y>0:
c	if ((sx.lt.0).and.(sy.gt.0)) then
c       gl=atan(dabs(sx/sy))+pi
c       endif
cc x<0, y<0:
c	if ((sx.lt.0).and.(sy.lt.0)) then
c       gl=atan(dabs(sy/sx))+3*pi/2
c       endif
       call carte_to_galactic(X,Y,Z,gl,gb)
C alpha, delta heliocentric
       cb=Dcos(gb)
       sb=Dsin(gb)
       cl=Dcos(gl)
       sl=Dsin(gl)
      
       sd=cb*Dsin(gl-a2)*Dsin(a1)+
     +     sb*Dcos(a1)
       Delta=Dasin(sd)
      
       yk=cb*Dsin(gl-a2)*Dcos(a1)-
     .     sb*Dsin(a1)
       xk=cb*Dcos(gl-a2)
       Alpha=mod(Datan2(yk,xk) +Ag,2.d0*pi)
       RETURN
       END



C  --------------------------------------------------------------------
c     SUBROUTINE carte_to_galactic(X,Y,Z,l,b)
c     Input cartesian coordinates: (X,Y,Z) in kpc
c     Output: galactic coordinates: (l,b) in radiants
c     -Pi<l<Pi; -Pi/2<b<Pi/2
c     checked!! :)
c---------------------------------------------------------------------
      SUBROUTINE carte_to_galactic(X,Y,Z,gl,gb)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

	tol=1.d-3


c XYZ, centered on the Sun
       sX=X+R0
       sY=Y!-R0
       sZ=Z

c R,l,b heliocentric
       gr=dsqrt(sx*sx+sy*sy+sz*sz)
       gb=dasin(sz/gr)
c        write(*,*)sz,gr,gb/deg
c x>0, y<0:
	if ((sx.ge.0.d0).and.(sy.lt.0.d0)) then
       gl=2.d0*pi-datan(dabs(sy/sx))
	if (sx.lt.tol) gl=3.d0*pi/2.d0
       endif
c x>0, y=>0:
	if ((sx.ge.0.d0).and.(sy.ge.0.d0)) then
       gl=datan(dabs(sy/sx))!+pi/2.d0
	if (sx.lt.tol) gl=pi/2.d0
       endif
c x<0, y=>0:
	if ((sx.lt.0.d0).and.(sy.ge.0.d0)) then
       gl=pi-datan(dabs(sy/sx))
	if (sx.gt.-1.d0*tol) gl=pi/2.d0
       endif
c x<0, y<0:
	if ((sx.lt.0.d0).and.(sy.lt.0.d0)) then
       gl=pi+datan(dabs(sy/sx))
	if (sx.gt.-1.d0*tol) gl=3.d0*pi/2.d0
       endif


c sx=0, y=0
	if((dabs(sx).le.tol).and.(dabs(sy).le.tol)) then
	gl=0.d0
	end if
c       write(28,*)sx,sy,gl/deg,gb
       RETURN
       END






      SUBROUTINE UVWH_TO_equatorials(U,V,W,A)
C*************************************************************************
C USING TRIGONOMETRIC MATRIX COMPUTE (VR,MUA,MUD) FROM (U,V,W)HELIOCENTRIC
c input velocities in km/s, parallax, xpi, in mas, angles in radiants
C a(1)=alpha
C a(2)=delta    the first three are input variables
C a(3)=xpi    
C a(4)=mua*cos(delta)
C a(5)=mud      the last three are output variables
C a(6)=Vr
c output proper motions in mas/yr and radial velocity in km/s
C*************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

      call get_fi(a(1),a(2),fi)
      SFI=DSIN(FI)
      CFI=DCOS(FI)
      call equ_to_gal(a(1),a(2),gl,gb)
C APPLY TRANSFORMATION MATRIX TO (U,V,W)INPUT
      cb=Dcos(gb)
      sb=Dsin(gb)
      cl=Dcos(gl)
      sl=Dsin(gl)

      B11=CL*CB 
      B12=CL*SB*SFI-SL*CFI
      B13=-(CL*SB*CFI+SL*SFI)
      B21=SL*CB
      B22=SL*SB*SFI+CL*CFI
      B23=-(SL*SB*CFI-CL*SFI)
      B31=SB
      B32=-CB*SFI
      B33=CB*CFI

      
      a(6)=B11*U+B21*V+B31*W
      a(4)=B12*U+B22*V+B32*W
      a(5)=B13*U+B23*V+B33*W
      gr=1.d3/a(3)
C CORRECT UNITS 
      a(4)=a(4)*1.d3/(kt*gr)
      a(5)=a(5)*1.d3/(kt*gr)


      RETURN
      END


      SUBROUTINE equatorials_TO_UVWH(A,U,V,W)
C*************************************************************************
C USING TROGONOMETRIC MATRIX COMPUTE (U,V,W)HELIOCENTRIC FROM (VR,MUA,MUD) 
c input observed variables in equatorial coordinates
C a(1)=alfa
C a(2)=delta    
C a(3)=xpi    
C a(4)=mua*cos(delta)
C a(5)=mud      
C a(6)=Vr
c output (U,V,W) heliocentric
c units: angles in radians, proper motions in mas/yr, Vr and U,V,W in km/s
C*************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      dimension A(6)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

C COMPUTE PARALLACTIC ANGLE FI
      call get_fi(a(1),a(2),fi)
      SFI=DSIN(FI)
      CFI=DCOS(FI)
      call equ_to_gal(a(1),a(2),gl,gb)
C APPLY TRANSFORMATION MATRIX TO (MUA,MUD,VR) INPUT
      cb=Dcos(gb)
      sb=Dsin(gb)
      cl=Dcos(gl)
      sl=Dsin(gl)

      B11=CL*CB 
      B12=CL*SB*SFI-SL*CFI
      B13=-(CL*SB*CFI+SL*SFI)
      B21=SL*CB
      B22=SL*SB*SFI+CL*CFI
      B23=-(SL*SB*CFI-CL*SFI)
      B31=SB
      B32=-CB*SFI
      B33=CB*CFI

C CORRECT FOR UNITS
      gr=1.d3/a(3)
      XMUAA=A(4)*1.d-3*GR*KT
      XMUDD=A(5)*1.d-3*GR*KT

      U=A(6)*B11+XMUAA*B12+XMUDD*B13
      V=A(6)*B21+XMUAA*B22+XMUDD*B23
      W=A(6)*B31+XMUAA*B32+XMUDD*B33
      RETURN
      END



      SUBROUTINE UVWH_TO_CARTE(UH,VH,WH,VX,VY,VZ)
C********************************************************************
C COMPUTES (VX,VY,VZ) FROM (U,V,W) HELIOCENTRIC
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'



C ADD THE DIFFERENTIAL ROTATION OF THE SUN
      CALL dV2(R0,0.d0,W)
      U1=UH
      V1=VH+W*R0
      W1=WH
C ADD THE MOTION OF THE SUN WITH RESPECT TO THE LSR
      U2=U1+USOL
      V2=V1+VSOL
      W2=W1+WSOL
C ROTATE 0 DEG TO OBTAIN CARTESIAN COORD
      ANGLE=0.!PI/2.D0
      CALL ROT_ANTIHORARI(ANGLE,U2,V2,VX,VY)
      VZ=W2
      RETURN
      END


c **********************************************************************

c subroutine equatorial_to_galactic(alfa,delta,gl,gb)
c input: equatorial coordinates (alfa,delta)
c output: galactic coordinates (xl,xb)
c Units: all angles l,b,alpha,delta in radians
c -pi<l<pi
c -------------------------------------------------------------

      subroutine equ_to_gal(alfa,delta,gl,gb)
      implicit real*8(a-h,o-z)
      real*8 lq,lpl
      include 'const_math.h'
      include 'const_ast.h'
      INCLUDE 'const_pot.h'

      T6=60.d0
      T36=3600.d0
      RD=PI/180.d0
      LQ=122.92851d0*RD
      ALFAG=(12.d0+51.d0/T6+26.2755d0/T36)*15.d0*RD
      DECG=(27.d0+7.d0/T6+41.704d0/T36)*RD
      SDG=dSIN(DECG)
      CDG=dCOS(DECG)
 
      sd=dsin(delta)
      cd=dcos(delta)

      sb=sdg*sd+cdg*cd*dcos(alfa-alfag)
      cb=sqrt(1.d0-sb**2.d0)
      gb=datan2(sb,cb)
      slpl=(cd*dsin(alfa-alfag))/cb
      clpl=(cdg*sd-sdg*cd*dcos(alfa-alfag))/cb
      lpl=datan2(slpl,clpl)
      if(lpl.le.0.d0)lpl=lpl+2.d0*pi
      gl=lq-lpl



      return
      end
*****************************************************
c input: equatorial coordinates (alfa,delta)
c output: galactic coordinates (l,b)

c Units: all angles l,b,alpha,delta in radians
c 0<l<2pi
************************************************
      subroutine lb_cal(al,de,l,b)

      implicit NONE
      real*8 pi, e,fa,fl, e2,fa2,fl2, e1,fa1,fl1,
     .        Equinox,l,b, al,de, w
      parameter ( pi = 3.1415926535897931159979634D+00 )
      parameter ( e = 62.8717d0 /180.d0*pi )
      parameter ( fa = 282.8596d0 /180.d0*pi )
      parameter ( fl = 32.93192d0 /180.d0*pi )

         l = datan2(dcos(de)*dsin(al-fa)*dcos(e)+dsin(de)*dsin(e),
     .              dcos(de)*dcos(al-fa)) + fl
         l = l  - 2*pi*int(l/(2*pi))
         if (l .lt. 0.d0)  l = l + 2*pi
         b = dasin(dsin(de)*dcos(e)-dcos(de)*dsin(al-fa)*dsin(e))


      return
      end




C  --------------------------------------------------------------------
c     SUBROUTINE equatorial_to_carte(gr,ALPHA,DELTA,x,y,z)
c     Input: equatorial coordinates: (alpha,delta) in radians
c     Input: gr:heliocentric distance in kpc.
c     Output cartesian coordinates: (X,Y,Z)
c     
c---------------------------------------------------------------------
      SUBROUTINE equatorial_to_carte(gr,ALPHA,DELTA,x,y,z)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

       
      call equ_to_gal(alpha,delta,gl,gb)
      cl=dcos(gl)
      sl=dsin(gl)
      cb=dcos(gb)
      sb=dsin(gb)
c
      x1=gr*cb*cl
      y1=gr*cb*sl
      z=gr*sb

c      call rot_antihorari(PI/2,x1,y1,x2,y2)

c to galactocentric
      x=x1-R0
      y=y1!+R0
      
      RETURN
      END



c********************************************************************
c********************************************************************
C AUXILIARY SUBROUTINES
C--------------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine dV2(R,zeta,w)
c
c     give(R,z) in cylindrical coordinates, it returns the rotation
c     frequency w at that position.
c     it uses the axisymmetric potential give by Allen & Santillan (1991)
c------------------------------------------------------------------------

      SUBROUTINE dV2(R,zeta,w)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'
      INCLUDE 'const_pot.h'

      zz=zeta**2.d0
      RR=R**2.d0


      q1=(zz+bdbd)**(1.d0/2.d0)
      q5=(ad+q1)**2.d0

c BULGE
      Psi1=rKMb/(RR+zz+abuabu)**(3.d0/2.d0)
      dVBdR=Psi1*R


c DISK
      Psi2=rKMd/(RR+q5)**(3.d0/2.d0)
      dVDdR=Psi2*R

c HALO
      XG=((RR+zz)**(1.d0/2.d0)/ah)**1.02d0
      Psi3=rKMh/ah*1.d0/(RR+zz)*XG/(XG+1.d0)
      dVHdR=Psi3*R
   
c ALL COMPONENTS
      dVdR=dVBdR+dVDdR+dVHdR

      w=dsqrt(1.d0/R*dabs(dVdR))

      RETURN
      END

c----------------------------------------------------------------
c     SUBROUTINE rot_positive(angle,X,Y,X2,Y2)
c
c     it performs a clockwise rotation of angle "angle"
c     input: rotation angle in radians
c            (x,y)
c     output: rotated coordinates (x2,y2)
c----------------------------------------------------------------
      SUBROUTINE rot_positive(angle,X,Y,X2,Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
    
      s=dsin(angle)
      c=dcos(angle)

      X2=X*c+Y*s
      Y2=-X*s+Y*c
 
      RETURN
      end



c----------------------------------------------------------------
c     SUBROUTINE rot_antihorari(angle,X,Y,X2,Y2)
c     transposed transformation matrix between two ref. frames.(?)
c     it performs a counter-clockwise rotation of angle "angle"
c     input: rotation angle in radians
c            (x,y)
c     output: rotated coordinates (x2,y2)
c----------------------------------------------------------------
      SUBROUTINE rot_antihorari(angle,X,Y,X2,Y2)
      IMPLICIT REAL*8 (A-H,O-Z)
    
      s=dsin(angle)
      c=dcos(angle)

      X2=X*c+Y*s
      Y2=-X*s+Y*c
 
      RETURN
      end

c----------------------------------------------------------------
c     SUBROUTINE rot_horari(angle,X,Y,X2,Y2)
c     transformation matrix between two ref. frames.(?)
c     it performs a clockwise rotation of angle "angle"
c     input: rotation angle in radians
c            (x,y)
c     output: rotated coordinates (x2,y2)
c----------------------------------------------------------------
      SUBROUTINE rot_horari(angle,X,Y,X2,Y2)
c rotacion en sentido horario de un angulo de "angle"
      IMPLICIT REAL*8 (A-H,O-Z)
    
      s=dsin(angle)
      c=dcos(angle)

      X2=X*c-Y*s
      Y2=X*s+Y*c
 
      RETURN
      end


c-------------------------------------------------------------------
      subroutine Carte_to_UVW(X,Y,Z,VX,VY,VZ,U,V,W)
C *******************************************************************
C CHANGES FROM (VX,VY,VZ) CARTESIAN GALACTOCENTRIC VELOCITY VECTOR TO
C (U,V,W) VELOCITY VECTOR WITH U POINTING TOWARDS TO GC, V TOWARDS 
C GALACTIC ROTATION AND W PERPENDICULAR TO THE GALACTIC PLANE
C WE NEED (X,Y,Z) TO COMPUTE THE POSITION ANGLE THETA.
c 
c Not used in the main code, but useful.
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)

      PI=4.D0*DATAN(1.D0)

C POSITION ANGLE
      THETA=DATAN2(X,Y)
C ROTATION OF VX,VY OF ANGLE 270-THETA COUNTERCLOCKWISE
c      ANGLE=3.D0*PI/2.D0-THETA
c      CALL rot_antihorari(angle,VX,VY,VX2,VY2)
      ANGLE=PI/2.D0+THETA
      CALL rot_horari(angle,VX,VY,VX2,VY2)
C SUBSTRACT THE DIFFERENTIAL ROTATION OF THE STAR    
      R=DSQRT(X*X+Y*Y)
      CALL dV2(R,Z,W)
      U=VX2
      V=VY2  ! -W*R
      W=VZ
      RETURN
      END


c----------------------------------------------------------
	subroutine gal_mu(dist,gl,gb,UH,VH,WH,xmub,xmuls,vr)
c***********************************************************
c Calculate the galactic proper motions using
cheliocentric distance and heliocentric V,U,W and l and b.
c Inputs: UH,VH,WH(km/s), gl, gb(radians), heliocentric distance(dist) in Kpc
c outputs: xMub=mu_b, xMuls= Mu_l*cosb in arcsec/yr
c***********************************************************
	IMPLICIT REAL*8(A-H,O-Z)
        INCLUDE 'const_math.h'
	INCLUDE 'const_ast.h'

	coef=1.d0/(kt*dist*1000.d0)

	xmub=coef*(-dsin(gb)*dcos(gl)*UH-dsin(gb)*dsin(gl)*VH
     &          +dcos(gb)*WH)
	xmuls=coef*(-dsin(gl)*UH+dcos(gl)*VH)
        vr=dcos(gb)*dcos(gl)*UH+dcos(gb)*dsin(gl)*VH+dsin(gb)
     &        *WH
	return
	end
 
c******************************************************************

c     subroutine get_fi
c input: (alfa,delta) in radians
c output: fi parallactic angle in radians

C  --------------------------------------------------------------------
      subroutine get_fi(alfa,delta,fi)
      implicit real*8(a-h,o-z)
      include 'const_math.h'
      include 'const_ast.h'

      x=Dcos(Dg)*Dsin(Delta)*Dsin(Alfa-Ag)+Dsin(Dg)*Dcos(Delta)
      y=Dcos(Dg)*Dcos(Alfa-Ag)
c Parallactic angle
      Fi=Datan2(y,x)
      return
      end

c******************************************************************
c subrutine pmeq_to_pmgal
c From J. Bovy's GalPy
      subroutine pmeq_to_pmgal(pmra,pmdec,ra,dec,pml,pmb)
      implicit real*8(a-h,o-z)
      include 'const_math.h'
      include 'const_ast.h'

c Bovy's constants
c      theta=122.932*pi/180.
c      dec_ngp=27.12825*pi/180.
c      ra_ngp=192.85948*pi/180.
c Merce's constants
      T6=60.d0
      T36=3600.d0
      RD=PI/180.d0
      theta=122.92851d0*pi/180.
      ra_ngp=(12.d0+51.d0/T6+26.2755d0/T36)*15.d0*RD
      dec_ngp=(27.d0+7.d0/T6+41.704d0/T36)*RD

      if(dec.eq.dec_ngp)dec=dec+1.e-16
      
      sindec_ngp= dsin(dec_ngp)
      cosdec_ngp= dcos(dec_ngp)
      sindec= dsin(dec)
      cosdec= dcos(dec)
      sinrarangp= dsin(ra-ra_ngp)
      cosrarangp= dcos(ra-ra_ngp)

      cosphi= sindec_ngp*cosdec-cosdec_ngp*sindec*cosrarangp
      sinphi= sinrarangp*cosdec_ngp
      xnorm= sqrt(cosphi**2.d0+sinphi**2.d0)
      cosphi= cosphi/xnorm
      sinphi= sinphi/xnorm

      pml= pmra*cosphi + pmdec*sinphi 
      pmb=-pmra*sinphi + pmdec*cosphi

      return 
      end
