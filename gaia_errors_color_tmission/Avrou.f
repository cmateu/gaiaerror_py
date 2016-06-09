c     ------------------------------------------------------------
c
c     Computation of the extinction: SWG-RD-03 (Drimmel, 2002)
c
c     Fortran version developed by Univ. Barcelona (January, 2003)
c     This version does not include the COBE Fortran routines for
c     computing the pixel number
c
c     The last grid of the Orion Arm has been rebuit at a higher
c     resolution to get rid of the artifacts that were present at 
c     around l=75
c
c MRG, Nov 2012: Modification of the program to a subroutine with input
c     d: heliocentric distance in kpc
c     (l,b): galactic longitud and latitude or (a,d): equatorial coordinates
c     itype: flag to determine type of coordinates
c            2 for equatorial
c            3 for galactic
c            it also has possibility to use ecliptic coord (itype=1)
c     afact: character 'y' or 'n' to determine whether we use scaling 
c            factors or not
c Output:
c     abs: the absorption at the give position
c     ------------------------------------------------------------
      subroutine DrimmelAv(dist,xx,yy,itype,afact,xabs,iflag,ifac)
      implicit double precision (a-h,o-z)

      integer npix
      real*8 l,ll
      parameter(zsun=0.015d0, xsun=-8.5d0, pi=3.141592654d0,
     *          npix=393216, deg_rad=pi/180.d0)
      character afact*1
      dimension a(2,2,2),delta(3),n(3),fac(3)

c     3D Grids from SWG-RD-02
      real*4 avgrid(151,151,51), avloc(101,201,51)

c     3D Grids from SWG-RD-03

c     Large scale 
      real*4 avdisk(151,151,51), avspir(151,151,51),
     *avori(76,151,51)

c     Local scale
      real*4 avori2(101,201,51), avdloc(31,31,51)

c     Rescaling factors
      real*4 glon(npix), glat(npix), rf(npix)
      integer ncomp(npix),pixnum(npix)

c     ncomp=1 ---> disk
c     ncomp=2 ---> spiral
c     ncomp=3 ---> orion


c     Reading the 3D grids:

c      write(*,*) 'Do you want to use rescaling factors (SWG-RD-03)? (Y/N
c     *)'
c      read(*,*) afact
c      write(*,*)itype,afact
      if((afact.eq.'n').or.(afact.eq.'N')) then

	open(1,file='avgrid.dat',status='old',form='unformatted')
	open(2,file='avloc.dat',status='old',form='unformatted')
        read(1) avgrid
        read(2) avloc
	close(1)
	close(2)

      else

        open(11,file='avdisk.dat',status='old',form='unformatted')
        open(12,file='avspir.dat',status='old',form='unformatted')
        open(13,file='avori.dat',status='old',form='unformatted')
        open(14,file='avori2.dat',status='old',form='unformatted')
        open(15,file='avdloc.dat',status='old',form='unformatted')
        read(11) avdisk
        read(12) avspir
        read(13) avori
        read(14) avori2
        read(15) avdloc
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)

        open(16,file='rf_allsky.dat',status='old',form='unformatted')
        read(16) pixnum,ncomp,glon,glat,rf
        close(16)


      endif


c     ------------------------------------------------------------
c
c     Coordinates input and conversion to galactic coordinates 
c
c     ------------------------------------------------------------

c      write(*,*) ''
c      write(*,*) "Coordinates of the object:"
c      write(*,*) '-------------------------'
c      write(*,'(A,$)') 'Enter a distance from the Sun, d (in Kpc): '
c      read (*,*) d

c95    write (*,*) ' Kind of coordinates to be introduced:'
c      write (*,*) '    1 = Ecliptic (lecl,becl)'
c      write (*,*) '    2 = Equatorial (alpha,beta)'
c      write (*,*) '    3 = Galactic (l,b)'
c      write(*,'(A,$)') ''
c      read(*,*) itype
c      write(*,*)''
c      if((itype.ne.1).and.(itype.ne.2).and.(itype.ne.3)) goto 95

        if(itype.eq.1) then
c40      write(*,'(A,$)')'Enter ecliptic longitude (-180 < l_ecl < 180):
c     * '
c        read(*,*) xx
          if(xx.gt.180.d0.or.xx.lt.-180.d0) then
          write(*,*) 'ecliptic longitude out of range'
          stop
c          goto 40
          endif

c41      write(*,'(A,$)') 'Enter ecliptic latitude (-90 < b_ecl < 90): '
c        read(*,*) yy
          if(yy.gt.90.d0.or.yy.lt.-90.d0) then
          write(*,*) 'ecliptic latitude out of range'
          stop
c          goto 41
          endif

        endif


        if(itype.eq.2) then
c30      write(*,'(A,$)') 'Enter alpha (0h < alpha < 24h): '
c        read(*,*) xx
          if(xx.gt.24.d0.or.xx.lt.0.d0) then
          write(*,*) 'alpha out of range'
          stop
c          goto 30
          endif
        xx=xx*360.d0/24.d0   

c71      write(*,'(A,$)') 'Enter delta (-90 < delta < 90): '
c        read(*,*) yy
          if(yy.gt.90.d0.or.yy.lt.-90.d0) then
          write(*,*) 'delta out of range'
          stop
c          goto 71
          endif

        endif

        if(itype.eq.3) then
c60      write(*,'(A,$)') 'Enter galactic longitude (0 < l < 360): '
c        read(*,*) xx
          if(xx.gt.360.d0.or.xx.lt.0.d0) then
          write(*,*) 'l out of range'
          stop
c          goto 60
          endif

c61      write(*,'(A,$)') 'Enter b (-90 < b < 90): '
c        read(*,*) yy
          if(yy.gt.90.d0.or.yy.lt.-90.d0) then
          write(*,*) 'b out of range'
          stop
c          goto 61
          endif

	endif



c     Conversion to galactic coordinates (l,b)

      call coord(itype,xx,yy,l,b) 

      d=dist

c     To avoid sin/cos equal zero in the denominator

      if(b.eq.0.d0.or.b.eq.90.d0.or.b.eq.-90.d0)b=b+1.d-10

      if(l.eq.0.d0.or.l.eq.180.d0.or.
     *   l.eq.90.d0.or.l.eq.270.d0)l=l+1.d-10



c     ------------------------------------------------------------
c
c     If the star is located out of the large-scale grid, its 
c     distance is reduced to the maximum size of the grid 
c
c     ------------------------------------------------------------


      ll=l*deg_rad
      bb=b*deg_rad

c     Greatest distance of the grid in z 
         d1=0.5d0/dabs(dsin(bb))-zsun/dsin(bb)

C     Greatest distance of the grid in x

c to avoid the last column of the input file where the
c extinction is wronge, we changed the coefficient of
c 15.d0 to be 14.8d0 :

         d2=14.8d0/dabs(dcos(ll))-xsun/dcos(ll)
	 d2=d2/dabs(dcos(bb))


c     Greatest distance of the grid in y
         d3=15.d0/dabs(dsin(ll))
	 d3=d3/dabs(dcos(bb))


c     We choose the lower maximum distance
         dmax=min(d1,d2)
         dmax=min(dmax,d3)


c     If the star distance is greater than dmax
c     we decrease this distance to d=dmax
c	write(12,123)d1,d2,d3,dmax,d
c123	format(5f8.3)

         if(d.gt.dmax) then
             d=dmax
	     iflag=1
	    else 
	    iflag=0
         endif



c     ------------------------------------------------------------
c
c     Conversion from (l,b,d) to (x,y,z)   
c
c     ------------------------------------------------------------

       x=d*dcos(bb)*dcos(ll)
       y=d*dcos(bb)*dsin(ll)
       z=d*dsin(bb)+zsun

       x_gc = x + xsun
       x_he = x


c      Defined only to consider the non-centered origin of the 
c      large scale grid for the Orion Arm

       n_cen = 0


c     ------------------------------------------------------------
c
c     Computation of the Av without considering rescaling factors
c     (SWG-RD-02)
c
c     ------------------------------------------------------------

      if((afact.eq.'y').or.(afact.eq.'Y')) go to 998 

c     large scale grid:

      if ((dabs(x).gt.1.d0).or.(dabs(y).gt.2.d0)) then 

	  delta(1)= 0.2d0
	  delta(2)= 0.2d0
	  delta(3)=0.02d0

          n(1)=151
          n(2)=151
          n(3)= 51

          call point (x_gc,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 

          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avgrid(i+ip-1,j+jp-1,k+kp-1)

          end do
          end do
          end do

      else 

c     local scale grid:

	  delta(1)=0.02d0
	  delta(2)=0.02d0
	  delta(3)=0.02d0

          n(1)=101
          n(2)=201
          n(3)= 51

          call point (x_he,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 

          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avloc(i+ip-1,j+jp-1,k+kp-1)

          end do
          end do
          end do

      endif 

      call interp (a,f_i,f_j,f_k,xabs)

      go to 999

998   continue 

c     ------------------------------------------------------------
c
c     Computation of the rescaling factors  (SWG-RD-03)
c
c     ------------------------------------------------------------

      call factors (l,b,pixnum,ncomp,glon,glat,rf,fac,ifac)

c for stars close to l=12.53, b=0.88  the 'rf_allsky.dat' file doesn't have
c value for their 'ncomp and 'rf' so we either can ignore these stars
c or giving them the same values as the ones in the nearst grid. by
c uncommenting the 4 lines below, you are ignoring these stars.

c	if(ifac.eq.0) then
c	return
c	end if
c	write(*,*)'fac',fac(1),fac(2),fac(3)
c     ------------------------------------------------------------
c
c     Computation of the Av considering rescaling factors(SWG-RD-03)
c
c     ------------------------------------------------------------

c     Disk component:

c       Large scale disk component:

        if ((dabs(x).gt.0.75d0).or.(dabs(y).gt.0.75d0)) then 

	  delta(1)= 0.2d0
	  delta(2)= 0.2d0
	  delta(3)=0.02d0

          n(1)=151
          n(2)=151
          n(3)= 51

          call point (x_gc,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 

          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avdisk(i+ip-1,j+jp-1,k+kp-1)

          end do
          end do
          end do

        else 

c       Local scale disk component:

	  delta(1)=0.05d0
	  delta(2)=0.05d0
	  delta(3)=0.02d0

          n(1)=31
          n(2)=31
          n(3)=51

          call point (x_he,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 

          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avdloc(i+ip-1,j+jp-1,k+kp-1)

          end do
          end do
          end do

        endif 

        call interp (a,f_i,f_j,f_k,abs_d)
c	write(*,*)'1',abs_d
c     Spiral component:

	  delta(1)= 0.2d0
	  delta(2)= 0.2d0
	  delta(3)=0.02d0

          n(1)=151
          n(2)=151
          n(3)= 51

          call point (x_gc,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 

          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avspir(i+ip-1,j+jp-1,k+kp-1)

          end do
          end do
          end do 

        call interp (a,f_i,f_j,f_k,abs_s)
c	write(*,*)'2',abs_s
c     Orion arm component:

c       Large scale orion arm component:

       if ((dabs(x).gt.1.d0).or.(dabs(y).gt.2.d0)) then

c     Greatest distance of the grid in z
         dori1=0.5d0/dabs(dsin(bb))-zsun/dsin(bb)

c     Greatest distance of the grid in y
         dori2=3.75d0/dabs(dsin(ll))
         dori2=dori2/dabs(dcos(bb))

c     Greatest distance of the grid in x
         if(cos(ll).gt.0.d0)dori3=2.375/dabs(dcos(ll))
         if(cos(ll).lt.0.d0)dori3=1.375/dabs(dcos(ll))
         dori3=dori3/dabs(dcos(bb))

c     We choose the lower maximum distance
         dmax_ori=min(dori1,dori2)
         dmax_ori=min(dmax_ori,dori3)

c        print*, 'dmax_ori=', dmax_ori

c     If the star distance is greater than dmax_ori
c     we decrease this distance to d=dmax_ori

         if(d.gt.dmax_ori) then
             d_ori=dmax_ori
         else
         d_ori=d
         endif

c     Conversion from (l,b,d) to (x,y,z)

       x_ori=d_ori*dcos(bb)*dcos(ll)
       y_ori=d_ori*dcos(bb)*dsin(ll)
       z_ori=d_ori*dsin(bb)+zsun

       x_gc_ori = x_ori + xsun

          delta(1)= 0.05d0
          delta(2)= 0.05d0
          delta(3)=0.02d0

          n(1)= 76 
          n(2)=151
          n(3)= 51

          n_cen=1


          call point (x_gc_ori,y_ori,z_ori,n,delta,i,j,k,
     *                n_cen,f_i,f_j,f_k)


          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          if(x_gc_ori.gt.0.d0) then
	  a(ip,jp,kp)=0.d0
	  else
          a(ip,jp,kp)=avori(i+ip-1,j+jp-1,k+kp-1)
	  endif

          end do
          end do
          end do

        else 

c       Local scale orion arm component:

	  delta(1)=0.02d0
	  delta(2)=0.02d0
	  delta(3)=0.02d0

          n(1)=101
          n(2)=201
          n(3)=51

          call point (x_he,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 


          do ip = 1,2
          do jp = 1,2
          do kp = 1,2

          a(ip,jp,kp)=avori2(i+ip-1,j+jp-1,k+kp-1)
          end do
          end do
          end do

        endif 

        call interp (a,f_i,f_j,f_k,abs_o)
c	write(*,*)'3',abs_o
c       Computation of the total absorption 

        xabs = fac(1)*abs_d + fac(2)*abs_s + fac(3)*abs_o

c       write(*,*) ''
c       write(*,*) 'Disk=',fac(1),' X ', abs_d,'=',fac(1)*abs_d
c       write(*,*) 'Spiral=',fac(2),' X ', abs_s,'=',fac(2)*abs_s
c       write(*,*) 'Orion=',fac(3),' X ', abs_o,'=',fac(3)*abs_o


999   continue

c        write(*,*) ''
c        write(*,*) 'Av=',xabs
c        write(*,*) ''
c      write(*,*)'aki',abs
      return
      end 

c     ..........................................................
c      
c     Subroutine for the computation of the rescaling factors
c     ..........................................................

      subroutine factors (l,b,pixnum,ncomp,glon,glat,rf,fac,ifac)
      implicit double precision (a-h,o-z)
      parameter(pi=3.141592654d0, deg_rad=pi/180.d0, npix=393216)
      real*4 glon(npix),glat(npix),rf(npix)
      real*8 fac(3),ll,l
      integer ncomp(npix),pixnum(npix)

c          fac(1)    !  disk component rescaling factor
c          fac(2)    !  spiral component rescaling factor
c          fac(3)    !  orion component rescaling factor 

      do i=1,3
      fac(i)=1.d0
      enddo

      distmin=999999.d0
      ll=l*deg_rad
      bb=b*deg_rad

        do i=1,npix

        gl=glon(i)*deg_rad
        gb=glat(i)*deg_rad

	dist=dacos(dsin(ll)*dsin(gl)+dcos(ll)*dcos(gl)*dcos((bb-gb)))

            if(dist.lt.distmin) then
              distmin=dist	
              imin=i
	    endif

	end do

c it seems that for some directions there is no data available
c so we set a flag to be zero for such cases.

       if(imin.lt.npix) then 
        icomponent=ncomp(imin)
        fac(icomponent)=rf(imin)
	ifac=1
      else if(imin.ge.npix) then
	icomponent=2
        fac(icomponent)=1.7584d0
	ifac=0
	write(*,*) l,b
	end if
c     print *, glon(imin), glat(imin),pixnum(imin)
      return 
      end

c     ..........................................................
c     
c     Subroutine for the transformation to (l,b)  coordinates  
c     ..........................................................

      subroutine coord (itype, xx, yy, l, b) 
      implicit double precision (a-h,o-z)
      real*8 l, l_gp
      parameter(ep=23.389291111d0, pi=3.141592654d0,
     *          alphag=192.859458333d0, deltag=27.12825d0, 
     *          l_gp=122.9285d0, deg_rad=pi/180.d0)

c     itype = 1        !Ecliptic Coordinates
c     itype = 2        !Equatorial Coordinates
c     itype = 3        !Galactic Coordinates  

      if(itype.eq.3) then
        l=xx
	b=yy
      else 

c       Conversion from ecliptic to equatorial

        if(itype.eq.1) then

          xx=xx*deg_rad
          yy=yy*deg_rad
          epp=ep*deg_rad

	  delta=(dasin(dsin(yy)*dcos(epp)+
     *           dcos(yy)*dsin(epp)*dsin(xx)))

	  alpha=dacos(dcos(yy)*dcos(xx)/dcos(delta))

        else
          alpha=xx*deg_rad
          delta=yy*deg_rad

        endif

c	Conversion from equatorial to galactic

          a_g=alphag*deg_rad
          d_g=deltag*deg_rad

	b=dasin(dsin(delta)*dsin(d_g)+
     *    dcos(delta)*dcos(d_g)*dcos((a_g-alpha)))/deg_rad
	l=l_gp-dasin(-dcos(delta)*dsin((a_g-alpha))
     *    /dcos(b))/deg_rad

      endif

      return 
      end 

c     ..........................................................

c     Computation of the grid point & interpolation factors 
c     ..........................................................

      subroutine point (x,y,z,n,delta,i,j,k,n_cen,f_i,f_j,f_k) 
      implicit double precision (a-h,o-z)
      dimension delta(3),n(3)


        if(n_cen.ne.1) then 
            exacti=x/delta(1)+(n(1)-1.d0)/2.d0+1.d0  
        else
            exacti=x/delta(1)+2.5d0*(n(1)-1.d0)+1.d0   
        endif 

        exactj=y/delta(2)+(n(2)-1.d0)/2.d0+1.d0   
        exactk=z/delta(3)+(n(3)-1.d0)/2.d0+1.d0

	i=int(exacti)
	j=int(exactj)
	k=int(exactk)


        if(i.eq.n(1))i=i-1
        if(j.eq.n(2))j=j-1
        if(k.eq.n(3))k=k-1

        if(i.eq.0)i=i+1
        if(j.eq.0)j=j+1
        if(k.eq.0)k=k+1

	f_i=exacti-i
	f_j=exactj-j
	f_k=exactk-k

      return 
      end 
c     ..........................................................

c     Linear interpolation inside the grids
c     ..........................................................

      subroutine  interp (a,f_i,f_j,f_k,xabs)
      implicit double precision (a-h,o-z)
      dimension a(2,2,2)

        
        p1=a(1,1,1)+(a(2,1,1)-a(1,1,1))*f_i

        p2=a(1,2,1)+(a(2,2,1)-a(1,2,1))*f_i
        
        p3=a(1,1,2)+(a(2,1,2)-a(1,1,2))*f_i

        p4=a(1,2,2)+(a(2,2,2)-a(1,2,2))*f_i

        p12=p1+(p2-p1)*f_j

        p34=p3+(p4-p3)*f_j

        xabs=p12+(p34-p12)*f_k


      return 
      end 
c     ..........................................................
