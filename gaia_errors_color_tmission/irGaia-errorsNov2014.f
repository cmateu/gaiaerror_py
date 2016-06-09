c ---------------------------------------------------------------------
c
c Subroutine that computes and assign Gaia errors to observables
c Errors are computed following Gaia Science Performance web page:
c http://www.cosmos.esa.int/web/gaia/science-performance
c
c Universitat de Barcelona
c Contact: Merce Romero-Gomez, Josep Manel Carrasco 
c
c ---------------------------------------------------------------------

c Updated: October 2014

c Inputs:
c ---------

c   * External tables:
c     ----------------
c	'gfactor-Oct2012.dat': geometrical factors (Table 2 of the web)
c	'gfactor-Jun2013.dat': geometrical factors (Table 6 of the web)
c 	'TableVr-Oct2014.dat': Radial Velocity coefficients

c	Input files for the definition of constants: 
c   		const_math.h : mathematical constants
c   		const_ast.h  : astronomical constants 


c   * Parameters from the star:
c     -------------------------

c     jflag:
c     	The code offers two options depending of the type of errors 
c     	you want to apply (see website): 
c     		jflag=1, errors computed from mean geometrical factors (see website)
c     		jflag=-1,errors computed considering the scanning law of the satellite 
c             		(the geometrical factors are then computed from the ecliptic
c             		latitude and the number transits)

c     V: V apparent magnitude 
c     VI: observed (V-I) colour index (affected by extinction), in magnitudes

c     a(i): Astrometric actual values for the star (not affected by errors)
c         a(1): Equatorial heliocentric right ascension (units: radians)
c         a(2): Equatorial heliocentric declination: delta (units: radians) 
c         a(3): Parallax: pi (units: mas)
c         a(4): Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         a(5): Equatorial proper motions: mud (units: mas/yr)
c         a(6): Radial velocity: Vr (units: km/s) 

c     ap(i): Atmospheric parameters actual values for the star (no affected by errors) 
c         ap(1): Effective temperature (units: kelvin)
c         ap(2): surface gravitye (logg) (units: dex) 
c         ap(3): [Fe/H] (units: dex) 
c         ap(4): Absorption A0 (Av~A0, see Bailer-Jones et al., 2011)



c Output values: 
c -------------
c Parameters of the star as observed by Gaia (affected by errors)
c Gaia errors assigned to each parameter:  

c The output parameters of the star as observed by Gaia are: 

c   Astrometric (+Vr) data affected by errors  
c     ao(i): Astrometric values for the star affected by Gaia errors
c         ao(1): Observed Equatorial heliocentric right ascension (units: radians)
c         ao(2): Observed Equatorial heliocentric declination: delta (units: radians) 
c         ao(3): Observed Parallax: pi (units: mas)
c         ao(4): Observed Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         ao(5): Observed Equatorial proper motions: mud (units: mas/yr)
c         ao(6): Observed Radial velocity: Vr (units: km/s) 
c     ae(i): Gaia standard deviation (1-sigma) 
c         ae(1): Standard deviation in Equatorial heliocentric right ascension in true arcs on the sky:: alpha_{*}=alpha*cos(delta)  (units: mas)
c         ae(2): Standard deviation in Equatorial heliocentric declination: delta (units: mas) 
c         ae(3): Standard deviation in Parallax: pi (units: mas)
c         ae(4): Standard deviation in Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         ae(5): Standard deviation in Equatorial proper motions: mud (units: mas/yr)
c         ae(6): Standard deviation in Radial velocity: Vr (units: km/s) 

c   Photometric data affected by errors (G, G_BP, G_RP, G_RVS)
c     p(i): Photometric actual values for the star (no affected by errors)
c         p(1): G magnitude (units: magnitudes) 
c         p(2): G_BP-G_RP   (units: magnitudes) 
c     po(i): Photometric values of the star affected by Gaia errors
c         po(1): Observed G magnitude (units: magnitudes) 
c         po(2): Observed G_BP-G_RP   (units: magnitudes) 
c     pe(i): Gaia standard deviation (errors, 1-sigma)
c         pe(1): Standard deviation in G magnitude (units: magnitudes) 
c         pe(2): Standard deviation in G_BP-G_RP   (units: magnitudes) 

c   Atmospheric parameters data affected by errors: 
c     apo(i): Atmospheric parameters of the star affected by Gaia errors
c         apo(1): Observed Effective temperature (units: kelvin)
c         apo(2): Observed surface gravity (logg) (units: dex) 
c         apo(3): Observed [Fe/H] (units: dex) 
c         apo(4): Observed Absorption A0 (Av~A0, see Bailer-Jones et al., 2011) (units: magnitudes) 
c     ape(i): Gaia standard deviation (1-sigma)
c         ap(1): Standard deviation in Effective temperature (units: kelvin)
c         ap(2): Standard deviation in surface gravity (logg) (units: dex) 
c         ap(3): Standard deviation in [Fe/H] (units: dex) 
c         ap(4): Standard deviation in Absorption A0 (Av~A0, see Bailer-Jones et al., 2011) (units: magnitudes) 


c **********************************************************************
        subroutine Gaia_errors(jflag,V,VI,a,ao,ae,p,po,pe,
     &             ap,apo,ape,GRVS,Llavor)

	    implicit real*8 (a-h,o-z)
        dimension a(6),ao(6),ae(6),p(2),po(2),pe(2),g(5)
        dimension ap(4),apo(4),ape(4)
        dimension xvi(11),xavr(11),xbvr(11)
        dimension xsb(20),xNOBS(20),xga(20),xgd(20),xgpi(20),xgmua(20),
     &            xgmud(20)
        real*8 gasdev
        character*10 xchar
        logical ifirst
        data ifirst/.true./
c        save xsb,xNOBS,xga,xgd,xgpi,xgmua,xgmud,xnmean,gmeana,gmeand,
c     &       gmeanpi,gmeanmua,gmeanmud
        save xsb,xga,xgd,xgpi,xgmua,xgmud
        INCLUDE 'const_math.h'
        INCLUDE 'const_ast.h'
        COMMON / tableVr / xvi, xavr, xbvr 
        integer Llavor

        if(ifirst)then
c Table 2 from the Gaia Science Performance Webpage. Only geometrical factors
c NOT including the effect of the number of passages (including 6.2% dead time)
c$$$        open(3,file='gfactor-Oct2012.dat')
c$$$            read(3,*)
c$$$            do i=1,20
c$$$              read(3,*)xsb(i),xk,xk,xNOBS(i),xga(i),xgd(i),xgpi(i),
c$$$     &                 xgmua(i),xgmud(i)
c$$$            enddo
c$$$            read(3,*)xchar,xk,xk,xnmean,gmeana,gmeand,gmeanpi,gmeanmua,
c$$$     &               gmeanmud
c$$$            close(3)

c Table 6 from the Gaia Science Performance Webpage. It includes the 
c geometrical factors and the effect of the number of passages
c (including 0% dead time, which affects the number of passages,N_obs, and
c also the mean number. Then it applies the same formula as above)
        open(3,file='gfactor-Jun2013.dat')
            read(3,*)
            do i=1,20
              read(3,*)xsb(i),xga(i),xgd(i),xgpi(i),xgmua(i),xgmud(i)
            enddo
            read(3,*)xchar,gmeana,gmeand,gmeanpi,gmeanmua,gmeanmud
            close(3)

       open(3,file='TableVr-Oct2014.dat')
            read(3,*)
            do i=1,11
              read(3,*)xvi(i),xavr(i),xbvr(i)
            enddo
            close(3)
        ifirst=.false.
        endif


c Get Gaia errors in parallax,magnitudes and Vr 
c        write(*,*)'aki'
        call ErrorsPiMag(v,vi,p,pe,empi,a(3))
        call ErrorsVr(v,vi,ae(6),GRVS)


c       Output errors on parallax and proper motion will be in mas
c        empi = empi /1e3

c       Mean errors (from mean g-geometrical factor)
        if(jflag.gt.0)then
          ae(3)=empi
          ae(1)=gmeana*ae(3)
          ae(2)=gmeand*ae(3)
          ae(4)=gmeanmua*ae(3)
          ae(5)=gmeanmud*ae(3)
        endif

c       Errors depending on the scanning law 
c       (number of passages and the geometrical factor)
        if(jflag.lt.0)then
c Version using Table 2            
c$$$          sbeta=dabs(0.9175*dsin(a(2))-0.3978*dcos(a(2))*dsin(a(1)))
c$$$
c$$$          call lininter(20,xsb,xNobs,sbeta,gns)
c$$$          call lininter(20,xsb,xga,sbeta,g(1))
c$$$          call lininter(20,xsb,xgd,sbeta,g(2))
c$$$          call lininter(20,xsb,xgpi,sbeta,g(3))
c$$$          call lininter(20,xsb,xgmua,sbeta,g(4))
c$$$          call lininter(20,xsb,xgmud,sbeta,g(5)) !!!change here
c$$$
c$$$          xnob=dsqrt(xnmean/gns)
c$$$
c$$$c          write(*,*)sbeta,g(1),xnob*g(1)
c$$$
c$$$          do i=1,5
c$$$          ae(i)=xnob*g(i)*empi
c$$$          enddo 

c Version using Table 6
          sbeta=dabs(0.9175*dsin(a(2))-0.3978*dcos(a(2))*dsin(a(1)))

          call lininter(20,xsb,xga,sbeta,g(1))
          call lininter(20,xsb,xgd,sbeta,g(2))
          call lininter(20,xsb,xgpi,sbeta,g(3))
          call lininter(20,xsb,xgmua,sbeta,g(4))
          call lininter(20,xsb,xgmud,sbeta,g(5)) 

          do i=1,5
          ae(i)=g(i)*empi
          enddo 
        endif

c       Computation of the errors in the atmospheric parameters
c       See Liu et al., 2012
c       A second order polinomial on G has been fitted to the
c       |Teff_eaneas_pq - Teff_real|

        call errorsAp(p(1),ape)
	
	
c       Computation of the observed astrometric quantities 

c       The error in right ascension ae(1) denotes true arc on the sky
c       so the right ascension shall be converted to that before the 
c       random error is assigned 
c       alpha_{*}=alpha*cos(delta) 
 
        a(1)=a(1)*dcos(a(2))

c	The errors in alpha* and  delta are in mas whereas alpha and delta are
c       in radiants, so (alpha*,delta)  are converted from radiants to mas  
c       before the random error is assigned 
 
	a(1)=a(1)/mas
	a(2)=a(2)/mas
	
        do i=1,6
        ao(i)=gasdev(Llavor,a(i),ae(i))
        enddo 

c       The true and observed (alpha*,delta) are converted from mas to radiants 
	a(1)=a(1)*mas
	a(2)=a(2)*mas
	ao(1)=ao(1)*mas
	ao(2)=ao(2)*mas

c       and the alpha_{*} is converted to alpha
        ao(1)=ao(1)/dcos(a(2))
c       we also convert back the input value
        a(1)=a(1)/dcos(a(2))

c       Computation of the observed photometric quantities 
        do i=1,2
        po(i)=gasdev(Llavor,p(i),pe(i))
        enddo

c       Computation of the observed atmospheric parameters
        do i=1,4
        apo(i)=gasdev(Llavor,ap(i),ape(i))
        enddo
        return
	end


c******************************************************************
c     Subroutine ErrorsAp
c
c     A second order polinomial on G has been fitted to the
c     |Teff_eaneas_pq - Teff_real|
c     Commented, there is the option to fitting a second order
c     polinomial to the sigmas of the Aeneas pq-model 
c     (see tableA4 of Liu & Bailer-Jones 2012)
c    
c     input data
c     p(1)=G           Gaia apparent magnitude (G)

c     output data
c     ape(1):          Uncertainty in Teff (Kelvin) 
c     ape(2):          Uncertainty in logg (dex) 
c     ape(3):          Uncertainty in [Fe/H] (dex) 
c     ape(4):          Uncertainty in A0 (mag) 

c------------------------------------------------------------------

      subroutine errorsAp(G,ape)

      implicit real*8 (a-h,o-z)
      dimension ape(4)

c     Fitting |Teff_eaneas_pq - Teff_real|
      ape(1)= 630.3d0   - 107.54d0  *G + 4.9682d0    *G**2
      ape(2)= 0.49865d0 - 0.04432d0 *G + 0.0017055d0 *G**2
      ape(3)= 0.81192d0 - 0.12226d0 *G + 0.0056669d0 *G**2
      ape(4)= 0.56387d0 - 0.093435d0*G + 0.0042024d0 *G**2

c     Fitting sigTeff (Eaneas pq-model)
c      ape(1)=67.289d0   - 7.8492d0  *G + 0.34687d0   *G**2
c      ape(2)= 0.31915d0 - 0.044319d0*G + 0.0020501d0 *G**2
c      ape(3)= 0.71709d0 - 0.115d0   *G + 0.0050649d0 *G**2
c      ape(4)= 0.17787d0 - 0.026371d0*G + 0.0011578d0 *G**2

      return 
      end


c******************************************************************
c     Subroutine ErrorsPiMag
c
c     input data
c     v=V           Johnson V apparent magnitude
c     vi=(V-I)c     Johnson-Cousins V-I colour

c     output data
c     p(1)=G                true Gaia G apparent magnitude
c     p(2)=G_BP-G_RP        true Gaia colour
c     pe(1)=sigG            uncertainty in G (in mag units)
c     pe(2)=sigBPRP         uncertainty in G_BP-G_RP colour (in mag units)
c     ae(3)=empi            uncertainty in parallax (in microarcsec units)

c------------------------------------------------------------------
      subroutine ErrorsPiMag(v,vi,p,pe,empi,xpi)

      implicit real*8 (a-h,o-z)
      dimension p(2),pe(2)

c---- Gaia photometry,
      p(2)=-0.0660d0+1.2061d0*vi-0.0614d0*vi*vi+0.0041d0*vi**3
      p(1)=V-0.0257d0-0.0924d0*vi-0.1623d0*vi*vi+0.0090d0*vi**3
      if(p(1).gt.20.) then
         empi=999999.999
         pe(1)=99.999
         pe(2)=99.999
         return
      endif

c---- Gaia error in parallax (assuming 70 transits) 
c---- A - before commissioning error model
      z=MAX(10**(0.4*(12.0-15.0)), 10**(0.4*(p(1)-15.0)))
c      empi=dsqrt(9.3+658.1*z+4.568*z**2)*(0.986+(1.0-0.986)*vi)

c-----B - after commissioning 
c----- new fitting provided by Rygl,Antoja,DeBruijne et al (October 2014)
c----- it assumes no dependance on the colour
c      empi=dsqrt(-11.5+706.1*z+32.6*z**2)*(0.986+(1.0-0.986)*vi)
c      empi=0.1*xpi
      empi=0.07*xpi
c      empi=0.16*xpi
c      empi=0.20*xpi
c      write(*,*)empi,xpi
      pe(1)=0.001*dsqrt(0.02076*z**2+2.7224*z+0.004352)
      pe(1)=pe(1)/dsqrt(70.d0)

      z=MAX(10**(0.4*(11.0-15.0)), 10**(0.4*(p(1)-15.0)))

      a=-0.003201*vi**3+0.0589*vi**2+0.3353*vi+0.7927
      b=-0.001019*vi**3+0.0244*vi**2+0.1756*vi+1.4684
      c=-0.004093*vi**3+0.0740*vi**2+0.2834*vi-3.4772
      sigBP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigBP=sigBP/dsqrt(70.d0)
      
      a=-0.006560*vi**3+0.1080*vi**2-0.6296*vi+1.4470
      b=-0.003280*vi**3+0.0540*vi**2-0.3148*vi+1.7856
      c=-0.007992*vi**3+0.1482*vi**2-0.7544*vi-3.7232
      sigRP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigRP=sigRP/dsqrt(70.d0)

c     We assume there is no correlation among colours
      pe(2)=dsqrt(sigBP**2+sigRP**2)

      return
      end



c******************************************************************
c     Subroutine ErrorsVr
c
c     input data
c     v=V           Johnson V apparent magnitude
c     vi=(V-I)c     Johnson-Cousins V-I colour

c     output data
c     ae(6)=Sigr            uncertainty in Vr (in km/s units)

c------------------------------------------------------------------
      subroutine ErrorsVr(v,vi,sigVr,GRVS)

      implicit real*8 (a-h,o-z)
      dimension xavr(11),xbvr(11),xvi(11)
      COMMON / tableVr / xvi, xavr, xbvr 
c
c     Vr errors 
c     
c---- Before commissioning errors     
      GRVS=V-0.0119D0-1.2092D0*VI+0.0188D0*VI*VI+0.0005D0*VI*VI*VI
      call lininter(11,xvi,xavr,vi,avr)
      call lininter(11,xvi,xbvr,vi,bvr)
      sigVr=1.d0+bvr*exp(avr*(V-14))   

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



c **********************************************************************
c          subroutine lininter(n,xa,ya,x,y)
c input:
c   xa(n) and ya(n) such that ya(i)=f(xa(i)), and xa ordered in 
c                   increasing sense
c   x, point where you want to know the corresponding y=f(x)
c output:
c   y
C  --------------------------------------------------------------------

           subroutine lininter(n,xa,ya,x,y)
           implicit real*8(a-h,o-z)
           dimension xa(n),ya(n)
          
           if(x.gt.xa(n))then
             y=ya(n)
             return
           else if(x.lt.xa(1))then
             y=ya(1)
             return
           else if(xa(1).le.x.and.x.le.xa(n))then
             i=1
             do while (xa(i).lt.x)
               i=i+1
             enddo
             i=i-1     !!!change here
             x0=xa(i)
             x1=xa(i+1)
             y0=ya(i)
             y1=ya(i+1)
             y=((x-x0)/(x1-x0))*y1+((x1-x)/(x1-x0))*y0
             return
           endif

           end



c **********************************************************************
c Coordinate transformation: from galactic to equatiorial 
c subroutine galactic_to_equatorial
c input: xi(l,b,pi,mul*cosb,mub,vr)
c output: xo(alpha,delta,pi,mua*cos(delta),mud,vr)
c Units:
c angles l,b,alpha,delta in radians
c proper motions in arcsec/yr
c note that pi(xi(3)) and vr(xi(6)) do not change
c------------------------------------------------------------------------

      subroutine galactic_to_equatorial(xi,xo)
      implicit real*8(a-h,o-z)
      dimension xi(6),xo(6)
      include 'const_math.h'
      include 'const_ast.h'

c l,b to alpha,delta
      cb=Dcos(xi(2))
      sb=Dsin(xi(2))
      cl=Dcos(xi(1))
      sl=Dsin(xi(1))

      sd=cb*Dsin(xi(1)-a2)*Dsin(a1)+sb*Dcos(a1)
      xo(2)=Dasin(sd)

      y=cb*Dsin(xi(1)-a2)*Dcos(a1)-sb*Dsin(a1)
      x=cb*Dcos(xi(1)-a2)
      xo(1)=Datan2(y,x) +Ag

c mulcosb,mub to muas,mud
      call get_fi(xo(1),xo(2),fi)
      cfi=dcos(fi)
      sfi=dsin(fi)
      xo(4)=cfi*xi(4)-sfi*xi(5)
      xo(5)=sfi*xi(4)+cfi*xi(5)
c
      xo(3)=xi(3)
      xo(6)=xi(6)
      return
      end

c **********************************************************************
c Coordinate transformation: from equatiorial to galactic 
c subroutine equatorial_to_galactic(xi,xo)
c input: xi(alpha,delta,pi,mua*cos(delta),mud,vr)
c output: xo(l,b,pi,mul*cosb,mub,vr)
c Units:
c angles l,b,alpha,delta in radians
c proper motions in arcsec/yr
c note that pi(xi(3)) and vr(xi(6)) do not change
c -------------------------------------------------------------

      subroutine equatorial_to_galactic(xi,xo)
      implicit real*8(a-h,o-z)
      dimension xi(6),xo(6)
      include 'const_math.h'
      include 'const_ast.h'

c alpha,delta to l,b
      sa=dsin(xi(1)-Ag)
      ca=dcos(xi(1)-Ag)
      sd=dsin(xi(2))
      cd=dcos(xi(2))
      sb=sd*dcos(a1)-cd*dsin(a1)*sa
      xo(2)=dasin(sb)
      xo(1)=datan((sd-dcos(a1)*sb)/(dsin(a1)*cd*ca))+a2

c muastar,mud to mulstar,mub
      call get_fi(xi(1),xi(2),fi)
      cfi=dcos(fi)
      sfi=dsin(fi)
      xo(4)= cfi*xi(4)+sfi*xi(5)
      xo(5)=-sfi*xi(4)+cfi*xi(5)
c
      xo(3)=xi(3)
      xo(6)=xi(6)
      return
      end
c -------------------------------------------------------------
