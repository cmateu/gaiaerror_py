c-----------------------------------------------------------------------
c  This is 'const_pot.h'
c  Constants related to galactic potential

      REAL*8 R0,w0,Mb,Md,Mh,Mbar,K
      
c valors generals      
      PARAMETER (R0 = 8.5d0, !valors de Palous et al. 1993 en kpc
     +           w0 = 25.8498243d0, !valors de Palous et al. 1993(fins ara hi havia 25.88)
     +           Usol = 10.00d0, ! velocitat solar (Dehen i Binney 1998) (abans (9,12,7))
     +           Vsol = 5.25d0,
     +           Wsol = 7.17d0)
     
     
c valors pel potencial no perturbat
      PARAMETER (ZMsol = 1.989d33,! Msol en gm
     +           abu = 0.387d0,! en Kpc
     +           ad = 5.3178d0,
     +           ah = 12.0d0,
     +           bd = 0.25d0,
     +           Mb = 1.406d10,!en Msol
     +           Md = 8.561d10,!en Msol
     +           Mh = 1.071d11)!en Msol
     
c friccio dinamica
      PARAMETER (fric=0.0d0)!coeficient de friccio dinamica
     
c definició braços espirals:
      PARAMETER (fr0 = 0.05d0,!
     +           Rsag = 8.225d0,!Posicio brac, de Sagitari, suposant una fase de l'estructura espiral de 330 graus i una distancia interbraç de 3.3 kpc
     +           pitch=6.0d0*deg,!pitch angle
     +           zm=2.0d0,!numero braços
     +           Omegaa=30.00d0)!velocitat rotacio      
c defibicio barra
      PARAMETER (Mbar = 1.0d9,
     +           ab=2.381d0,!abar/bbar
     +           ac=3.030d0,!abar/cbar
     +           qb=5.0d0,!qbar
     +           tb0=5.0d0,
     +           angle0=20.0d0*deg,
     +           Omegab=58.55d0,
     +           td0=100.d0)
                
c calculs per fer el programa mes rapid:
       PARAMETER(w0w0=w0*w0,
     +           w02=2.d0*w0,
     &           K=Gte*ZMsol/pc/1.d13,
     +           abuabu=abu*abu,
     +           bdbd=bd*bd,
     +           rKMb=K*Mb,
     +           q9=-3.d0*rKMb,
     +           rKMd=K*Md,
     +           q8=-3.d0*rKMd,
     +           rKMh=K*Mh,
     +           abab=ab*ab)

       PARAMETER( RRww=R0*R0*w0*w0,
     +           RRwwfr0=RRww*fr0)

c definicio braços espirals transitoris:
      PARAMETER(pitcht=6.d0*deg,
     +           zmt=2.d0,
     +           eps=1.d0,
     +           t0=100.d0,
     +           sigts=0.35d0,
     +           sigw=w0/4.d0,
     +           Ns=12)

       PARAMETER(RRwwfr0eps=RRwwfr0*eps)

c definicio braços espirals d estructures tipus 2+4:
       PARAMETER(Nse=2)

c valors per les constants de les espirals 
         tanpitch=dtan(pitch)
         Phi0=Pi+zm/tanpitch*dlog(Rsag/R0)
         zmtanpitch=zm/tanpitch

c valors per les constants de les espirals transitories
         tanpitcht=dtan(pitcht)
         zmttanpitcht=zmt/tanpitcht


