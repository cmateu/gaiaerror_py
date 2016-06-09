C  (C) Copr. 1986-92 Numerical Recipes Software 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ran2(idum)
c     random double precision number on [0,1]
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1.0d0/im1,imm1=im1-1,
     *     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *     ntab=32,ndiv=1+imm1/ntab,eps=1.2d-7,rnmx=1.0d0-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum=idum+im1
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function gasdev(idum,mean,var)
c     adapted from numerical recipes function of same name
c     to include arbitrary mean and variance for gaussian dist.
c     price is that we do not bother to save one of the dist. 
c     numbers (gset)
c     generated since it is presumed that the mean and/or 
c     variance will likely change between calls
      integer idum
      double precision gasdev,mean,var,rsq,v1,v2,ran2,var2
      var2 = var*var
 1    v1=var*(2.0d0*ran2(idum)-1.0d0)
      v2=var*(2.0d0*ran2(idum)-1.0d0)
      rsq=v1*v1+v2*v2
      if (rsq.ge.var2.or.rsq.eq.0.0d0) goto 1
      gasdev=v2*sqrt(-2.0d0*log(rsq/var2)/(rsq/var2)) + mean
      return
      end
