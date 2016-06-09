        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(6),AE(6),AO(6),p(2),po(2),pe(2),ap(4),apo(4),ape(4)
	 dimension xab(36),aogal(6),agal(6)
        DIMENSION AGO(6)
	character*150 arxi, arxo
       character afact*1
        INCLUDE 'const_math.h' 
        INCLUDE 'const_ast.h'
        INCLUDE 'const_pop_color.h'
        INCLUDE 'const_pot.h'

c Initializations:
        idum=-15350
  
        WRITE(*,*)'INPUT FILE WITH GALACTOCENTRIC COORDINATES'
        READ(*,*)ARXI
        WRITE(*,*)ARXI

        WRITE(*,*)'OUTPUT FILE TO WRITE OBSERVED GALACTOCENTRIC COORD.'
        READ(*,*)ARXO
        WRITE(*,*)ARXO

        NCI=23
        OPEN(NCI,FILE=ARXI)

        NCO=4
        OPEN(NCO,FILE=ARXO)

	iff=0
	ies=0
	ia=0
	ifa=0

9	continue

c Reading the galactocentric coordinates of the star
        ip=0
10	continue
	read(nci,*,err=10,end=111)X,Y,Z,VX,VY,VZ,aMv,vi_o
c        write(*,*)X,Y,Z,VX,VY,VZ,aMv,vi_o
        ip=ip+1
      
        

c Transform the position of the star from galactocientric to heliocentric 
c coordinates 
         dist=dsqrt((x+R0)*(x+R0)+y*y+z*z)   !in kpc
         distpc=dist*1000.d0                 !in pc
         xpi=1.d0/dist                       !parallax in mas
         CALL carte_to_equatorial(X,Y,Z,ALPHA,DELTA)
         a(1)=alpha
         a(2)=delta
         a(3)=xpi
         
C Transforming velocities from galactocentric to galactic helicocentric 
c coordinates and to equatorial heliocentric coordinates
         call Carte_to_UVWH(VX,VY,VZ,UH,VH,WH)
         CALL UVWH_to_equatorials(UH,VH,WH,A)
         
c galactocentric coordinates to galactic coordinates 
c         call carte_to_galactic(X,Y,Z,gl,gb)
	 call lb_cal(ALPHA,DELTA,gl,gb)
         call gal_mu(dist,gl,gb,UH,VH,WH,xmub,xmuls,vr)
         
c---------------------------------------------------------------
c Drimmel extinction law, input itype=2 (equatorial cood) and 
c afact='y', using scaling factors
         itype=3
         afact='y'
         call  DrimmelAv(dist,gl/deg,gb/deg,itype,afact,Av,ifg,ifac)
c         write(*,*)Av 
c------------------------------------------------------------------
c do not consider stars with Av greater than 40 (fittings of G not considered)
         if(Av.gt.40.d0)goto 10
c  K apparent magnitude of the source 
c         xK=aMk+5.d0*log10(distpc)-5.d0+0.112d0*Av
c  V apparent magnitude of the source 
         V=aMv+5.d0*log10(distpc)-5.d0 + Av
c Compute the observed (V-I) of the source
         vi=factE*Av+vi_o

c Compute the observed (J-K) of the source
c         xJK=factE2*Av+jk_o

c Assign atmospheric parameters to the source
         ap(1)=Teff    
         ap(2)=xlogg  
         ap(3)=FeH    !No radial metallicity gradient is assumed here
         ap(4)=Av     !we hare assuming A0=Av

c Introduce Gaia errors
         jflag=-1                          !weighted errors
c       jflag= 1                          !mean errors
c         write(*,*)V,VI,Av
         call Gaia_errors(jflag,V,VI,a,ao,ae,p,po,pe,ap,apo,
     &                      ape,GRVS,idum)
c         write(*,*)p(1)
c	if (iflag.eq.0) then
c	 iff=iff+1
c	go to 10
c	end if
c---------------------------------------------------------
C transform velocities from equatorial heliocentric to 
c galactocentric 
         CALL equatorials_to_UVWH(ao,Uobs,Vobs,Wobs)
         CALL UVWH_TO_Carte(Uobs,Vobs,Wobs,Vxobs,VYobs,VZobs)  

c transform heliocentric to galactocentric positions 
         distobs=(1.d0/ao(3))
         call equatorial_to_carte(distobs,ao(1),ao(2),xobs,yobs,zobs)
         call lb_cal(ao(1),ao(2),globs,gbobs)
         distobs=1.d0/ao(3) !in kpc

	call gal_mu(distobs,globs,gbobs,Uobs,Vobs,Wobs,xmubo,xmulso
     &              ,vro)
c Computation of (U,V,W) with errors
        call Gal_to_UVW(Xobs,Yobs,Zobs,VXobs,VYobs,VZobs,Uobs,Vobs,Wobs)

c Computation of (U,V,W) without errors
         call Gal_to_UVW(X,Y,Z,VX,VY,VZ,UU1,VV1,WW1)

         Gmag=p(1)
         vrad=a(6)
         vrob=ao(6)
         xpar=a(3)
         xparobs=ao(3)
         relerr=AE(3)/A(3)
         relerr_mu=AE(5)/A(5)
         relerr_vrad=AE(6)/A(6)
        
         if (Gmag.le.20.7) then

         WRITE(NCO,100)Av,V,Gmag,GRVS,relerr,X,Y,Z,VX,VY,VZ,
     &              xpar,gl/deg,gb/deg,dist,xmuls*1000.d0,
     &              xmub*1000.d0,vrad,Xobs,Yobs,Zobs,VXobs,VYobs,VZobs
     &              ,xparobs,globs/deg,gbobs/deg,distobs
     &              ,xmulso*1000.d0,xmubo*1000.d0,vrob
         endif
 100     FORMAT(31(F18.10,1X))
         GOTO 10
 111      CLOSE(NCI)
        CLOSE(nco)

        END




