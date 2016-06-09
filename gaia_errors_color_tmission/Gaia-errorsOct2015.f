c ---------------------------------------------------------------------
c
c Subroutine that computes and assign Gaia errors to observables
c Errors are computed following Gaia Science Performance web page:
c http://www.rssd.esa.int/index.php?project=GAIA&page=Science_Performance
c 
c Note: As discussed in Gaia-JDB-022, an overall end-of-mission scientific 
c contingency margin of 20% is included in the Astrometric and Photometric errors. No 
c correlation is assumed between astrometric errors.
c
c The code allows to compute the astrometric and photometric errors for a fraction of the 
c operation data. The adopted strategy is described in Gaia-C9-TN-UB-RMC-001-1.
c Warning: Only end-of-mission errors for radial velocity and astrophysical parameters are 
c provided.
c
c Universitat de Barcelona
c Contact: Merce Romero-Gomez, Josep Manel Carrasco, Roger Mor
c
c ---------------------------------------------------------------------

c Updated: September 2015

c Inputs:
c ---------

c   * month: Length of operation data released, in months.
c            month=60 means end-of-mission (5 years mission).

c     WARNING: It is important to keep in mind that the astrometric errors computed for a 
c fraction of operation data (month<60) using this strategy will probably be understimated. 
c This is because they are computed assuming that the instruments are going to be perfectly 
c calibrated at time "month". To account for this fact, we introduce:
c
c     CAfactor: Calibration Astrometric factor.
c     CAfactor=1.0 No Calibration Astrometric factor applied.
c     CAfactor=1.2 means Astrometric errors are increased in 20%.      




c   * External tables:
c     ----------------
c	'gfactor-Jun2013.dat': geometrical factors (Table 6 of the web)
c 	'TableVr-Jun2015.dat': Radial Velocity coefficients

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
c         ape(1): Standard deviation in Effective temperature (units: kelvin)
c         ape(2): Standard deviation in surface gravity (logg) (units: dex) 
c         ape(3): Standard deviation in [Fe/H] (units: dex) 
c         ape(4): Standard deviation in Absorption A0 (Av~A0, see Bailer-Jones et al., 2011) (units: magnitudes) 


c **********************************************************************
        subroutine Gaia_errors(month,CAfactor,jflag,V,VI,a,ao,ae,p,po,pe,
     &             ap,apo,ape,GRVS,Llavor)


	    implicit real*8 (a-h,o-z)
        dimension a(6),ao(6),ae(6),p(2),po(2),pe(2),g(5)
        dimension ap(4),apo(4),ape(4)
        dimension xvi(11),xavr(11),xbvr(11)
        dimension xsb(20),xNOBS(20),xga(20),xgd(20),xgpi(20),xgmua(20),
     &            xgmud(20)
        Double precision gasdev, month, factorL, VI, Xfactormua
	Double precision Xfactormud, L, CAfactor
        character*10 xchar
        logical ifirst
        data ifirst/.true./
        save xsb,xga,xgd,xgpi,xgmua,xgmud
        INCLUDE 'const_math.h'
        INCLUDE 'const_ast.h'
        COMMON / tableVr / xvi, xavr, xbvr 
        integer Llavor


	factorL=sqrt(60.d0/month)
	L=month/12.d0
        if(ifirst)then

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


   

c Exponential fit from the Spectroscopic Performace Table provided in the Gaia Science 
c Performance webpage. 

       open(3,file='TableVr-Jun2015.dat')
            read(3,*)
            do i=1,11
              read(3,*)xvi(i),xavr(i),xbvr(i)
            enddo
            close(3)
        ifirst=.false.
        endif


c Get Gaia errors in parallax,magnitudes and Vr 

        call ErrorsPiMag(factorL,CAfactor,v,vi,p,pe,empi)
        call ErrorsVr(v,vi,ae(6),GRVS)


c       Output errors on parallax and proper motion will be in mas
        empi = empi /1e3

c       Mean errors (from mean g-geometrical factor)
        if(jflag.gt.0)then
	  Xfactormua=gmeanmua*5. 
 	  Xfactormud=gmeanmud*5. 
          ae(3)=empi
          ae(1)=gmeana*ae(3)
          ae(2)=gmeand*ae(3)
          ae(4)=(Xfactormua/L)*ae(3)
          ae(5)=(Xfactormud/L)*ae(3)
        endif

c       Errors depending on the scanning law 
c       (number of passages and the geometrical factor)
        if(jflag.lt.0)then


          sbeta=dabs(0.9175*dsin(a(2))-0.3978*dcos(a(2))*dsin(a(1)))

          call lininter(20,xsb,xga,sbeta,g(1))
          call lininter(20,xsb,xgd,sbeta,g(2))
          call lininter(20,xsb,xgpi,sbeta,g(3))
          call lininter(20,xsb,xgmua,sbeta,g(4))
          call lininter(20,xsb,xgmud,sbeta,g(5)) 

          do i=1,3
          ae(i)=g(i)*empi
          enddo

	  Xfactormua=g(4)*5
	  Xfactormud=g(5)*5
          ae(4)=(Xfactormua/L)*empi                
	  ae(5)=(Xfactormud/L)*empi                  
		 
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
c     WARNING: See Bailer-Jones et al. (2013) table 7 for a more updated values.
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
c     factorL=sqrt(60/month)
c     CAfactor: Calibration Astrometric Factor
c     v=V           Johnson V apparent magnitude
c     vi=(V-I)c     Johnson-Cousins V-I colour

c     output data
c     p(1)=G                true Gaia G apparent magnitude
c     p(2)=G_BP-G_RP        true Gaia colour
c     pe(1)=sigG            uncertainty in G (in mag units)
c     pe(2)=sigBPRP         uncertainty in G_BP-G_RP colour (in mag units)
c     ae(3)=empi            uncertainty in parallax (in microarcsec units)

c------------------------------------------------------------------
      subroutine ErrorsPiMag(factorL,CAfactor,v,vi,p,pe,empi)

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


c Astrometric Error
c After Comissioning Nov 2014
      z=MAX(10**(0.4*(12.0-15.0)), 10**(0.4*(p(1)-15.0)))
      empi=dsqrt(-1.631+680.766*z+32.732*z**2)*(0.986+(1.0-0.986)*vi)
      empi=empi*factorL*CAfactor

c Photometric Errors
      pe(1)=0.001*dsqrt(0.04895*z**2+1.8633*z+0.0001985)
      pe(1)=(pe(1)/dsqrt(70.d0))*factorL

      z=MAX(10**(0.4*(11.0-15.0)), 10**(0.4*(p(1)-15.0)))

      a=-0.000562*vi**3+0.044390*vi**2+0.355123*vi+1.043270
      b=-0.000400*vi**3+0.018878*vi**2+0.195768*vi+1.465592
      c= 0.000262*vi**3+0.060769*vi**2-0.205807*vi-1.866968
      sigBP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigBP=sigBP/dsqrt(70.d0) * factorL
      
      a=-0.007597*vi**3+0.114126*vi**2-0.636628*vi+1.615927
      b=-0.003803*vi**3+0.057112*vi**2-0.318499*vi+1.783906
      c=-0.001923*vi**3+0.027352*vi**2-0.091569*vi-3.042268
      sigRP=0.001*dsqrt((10.**a)*(z**2)+(10.**b)*z+10.**c)
      sigRP=sigRP/dsqrt(70.d0)* factorL

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
c     GRVS         Spectrograph magnitude
c------------------------------------------------------------------
      subroutine ErrorsVr(v,vi,sigVr,GRVS)

      implicit real*8 (a-h,o-z)
      dimension xavr(11),xbvr(11),xvi(11)
      COMMON / tableVr / xvi, xavr, xbvr 
c
c     Vr errors 
c     
     
      GRVS=V-0.0119D0-1.2092D0*VI+0.0188D0*VI*VI+0.0005D0*VI*VI*VI
      call lininter(11,xvi,xavr,vi,avr)
      call lininter(11,xvi,xbvr,vi,bvr)

      sigVr=1.d0+bvr*exp(avr*(V-12.7d0))
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
             i=i-1    
             x0=xa(i)
             x1=xa(i+1)
             y0=ya(i)
             y1=ya(i+1)
             y=((x-x0)/(x1-x0))*y1+((x1-x)/(x1-x0))*y0
             return
           endif

           end




