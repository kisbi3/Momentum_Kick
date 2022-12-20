c     This is ppbb.f
c                     that is beased on ryieldbb.f,  
c     integration in b0 (arxiv:0901.0726)
c     include two-D transverse expansion
c     This   ppbb.f  program is for pp collision
c     while  ryieldbb.f,     is for AA collision

c     Note that the output aNRidge is charged particle number per ntrig
c                    aNRidge=(2/3)ANkick, aNkick includes both charged+neutral

c  to compile f95 -o ppbb.x ppbb.f ryieldsub.f?
c                       ryieldsub.f contains basic unchanged subroutines 

cc
CC    SUBROUTINE aaxsec
CC    
CC    
CC    TO CALCULATE Tab(b) and sampling function XI(b)
cc    Tab(b)=INT dba Ta(ba) * Tb(b-ba), just as in Glauber
cc    xi(b)=int(from 0 to b) 2pi*b*db(1-(1-sigma*Tab)**AB) /sigma(AB) 
cc    xi range from 0 to 1
cc       FOR THE COLLISION OF B ON A
CC       B PROJECTILE NUCLEONS TO COLLIDE WITH A TARGET NUCLEONS
CC    BOR(IB) IS IN UNITS OF RADIUS RA WITH STEP DB/RA
CC    B  (IB) IS IN UNITS OF FM WITH STEP DB
CC    WE NEED A TABLE OF B AND TB TO GO TO B ABOUT 1.4R
CC    BNNOR(IBNN) IS THE B FOR NUCLEUS-NUCLEUS SYSTEM WITH STEP DBNN/RA
CC    BNN  (IBNN) IS THE B FOR NUCLEUS-NUCLEUS SYSTEM WITH STEP DBNN
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /SPLINT/  y2(100) 
      COMMON /BANDTB/  BB(100),TB(100),DBB,BMAXB,
     2                 BA(100),TA(100),DBA,BMAXA,IBMAX
      common /extrap/  dtadb(100),dtbdb(100)
      COMMON /NULNUL/  BNN(100),BNNOR(100),TAB(100),XI(100),
     2                 DBNN,DBNNOR,IBNNMX,IBPRNT
      COMMON /XSECNN/  B,A,SIGIN
      common /radii /  RA,RB
      common /bcoord/  bbpx(100),bbpy(100),bboa(100,100),bbob(100,100),
     x                 ttoa(100,100),ttob(100,100)        
      common /paramt/  zeta,dnptondnpar,sigjp,t0   
      common /paramb/  bsep        
      common /exnorm/  ridyldanorm,chgmulanorm
             ! zeta=exp(-zeta*N),dnptondpar=dNparton/dnparpant,sigmjp,t0
             ! bsep=impact parameter
             ! t0  = starttime of r-expansion, measured rel to jet&z-expstart
             ! ridyldanorm -- for (ridge  yield      ) normalization
             ! chgmulanorm -- for (charge multiplcity) normalization
      common /surviv/  totcolf(100),totparf(100),
     x                 dcolda (100),dparda (100)
      common /colpar/  ancolxy(100,100),anparxy(100,100),
     x                 btxy   (100,100),atxy   (100,100),
     x                 kcon   (100,100)      ! kont=0/1 for zero/non-zero cont
      common /impact/  IBNNwhch,ibxyplt,IPRNk,IPLNk
                  ! write out IBNNwhch if ibxyplt=1,write and plot Nk if =1 
      common /IBNNnk/  IBNNkick0,IBNNdelta   ! get IBNNkick0+n*IBNNdelta
      common /Iphis /  phis0deg,dphisdeg,nphismx  ! phis0(deg),dphis,nphismx
      common /Nkdist/  aNktotfis,aNtrgfis,aNcolfis !aNk-all,aNtrg,aNcol[fis]
      common /NRidge/  aNRbfi(100,20),aNtrgbfi(100,20),RAAbfi(100,20),
     >                 aNRb(100),RAAb(100),IBNNkick(100) 
              ! RAAb(B)=Ntrig(B)/tot colision number, gives jet quenching
      DIMENSION dNpsi(100),dNDY(100)
      dimension et(100),dpsidet(100),dDYdet(100)

      DIMENSION BOR(100),TNN(100)
      DIMENSION ENDPTS(2),BSCR(100)
      Character*4 SHAPB,SHAPA
	open ( 5,file='ppbb.05',status='unknown')
        open ( 6,file='ppbb.06',status='unknown') 
        open ( 8,file='ppbb.08',status='unknown')  ! plot(chm(Parpt)#,yield)(B)
        open ( 9,file='ppbb.09',status='unknown')  ! plot(chm(Parpt)#,yie)(B,fi)
        open (10,file='ppbb.10',status='unknown')  ! plot(chm(Parpt)#,RAA)(B)
        open (11,file='ppbb.11',status='unknown')  ! plot(chm(Parpt)#,RAA)(B,fi)

        write (8,701)
        write (9,701)
        write (10,701)
        write (11,701)
 701    format('@TYPE XY')

CC  B IS PROJECTILE, A IS TARGET
CC  IBPRNT=0/1 NO PRINT/PRINT TAB,WNM FOR B INTERVALS
CC  NMAX=MAX NO OF PROJECTILE NUCLEONS IN A TUBE WHICH COLLIDE
CC  MMAX=MAX NO OF TARGET     NUCLEONS IN A TUBE WHICH COLLIDE
CC  IN THIS PROGRAM, WE DIMENSION CUT-OFF NMAX,MMAX .LE. THAN 20
CC  SIGIN=NUCLEON-NUCLEON INELASTIC CROSS-SECTION IN FM**2 (E.G. 2.94)
CC  IF SHAP_='WDSX',A WOOD-SAXON DENSITY IS USED,NEED R0 AND DIF
CC  IF SHAP_='GAUS',A GAUSSIAN DENSITY IS USED,NEED R0RMS(&R0 FOR GRID PTS)
CC  DBOR,DZOR=DB,AND IBMAX USED IN TABULATION OF TA(B) AND TB(B)
CC  DZOR,IZMAX=DZ,IZMAX USED TO GET THE INTEGRAL FOR TA(B) AND TB(B)
CC  IBNNMX=MAX NO OF NUCLEUS-NUCLEUS B CONSIDERED
CC  IBNNPR=IN OUPUT, PRINT EVERY IBNNPR B VALUES OUTPUT
      amp=0.938272d0              ! proton mass
 1    read (5,*,  end=999) zeta,dnptondnpar,sigjp,t0
                   ! zeta=anttenuation exp(-zeta*N)
                   ! dnptondnpar=dnparton/dnparpant e.g. 21
                   ! sigjp=sigma(jet-par in fm**2) 
                   ! t0=starttime of r-hydro-evolution expansion fm/c
                   ! tt=measured wrt jet&zhydrostart
      write(6,603)         zeta,dnptondnpar,sigjp,t0 
 603  format(' program ppbb.f, to calcualte kicked N_parton(phi_s)',
     x   /,' attenuation factor,exp(-zeta*N),zeta=',f10.5,
     x   /,' dnptondnpar=dnparton/dnparpant=',f10.5,
     x   /,' sigjp=sigma(jet-p`arton xsec in fm**2(=10mb))',f10.5,
     x   /,'    and t0=starttime of rhydro evolutn, (fm/c)=',f10.5) 
      write(6,*) ' sig(jet-parton in mb)=',10*sigjp
      READ (5,*) SIGIN    ! in fm**2,NN-sigIN 
      write(6,*) ' SIGIN (input is and kept in fm**2)=',SIGIN 
      WRITE(6,615) SIGIN
 615  format(' NN-SIGIN=', f10.5,' normalized to 1 in pp collisn') 

      READ (5,*)   R0B,R0A
      WRITE(6,632) R0B,R0A
      READ (5,*)   DBOR,IBMAX
      WRITE(6,640) DBOR,IBMAX
      READ (5,*)   DBNNOR,IBNNMX,IBNNPR
      WRITE(6,648) DBNNOR,IBNNMX,IBNNPR

    5 READ (5,*,END=999) B,A,IEOPT,EGEV,IBPRNT ! B,A,IEOPT=0plab/1cms,GEV,ibPr 
      write(6,649)       B,A,IEOPT,EGEV,IBPRNT
 649  format(' B=',f7.2,'  a=',f7.2,/,'   IEOPT=',i3,
     x       ' (0 for plab on fixed targ/ 1 for sqrts(NN))',
     x       /,' Gev=',f10.3,' IBPRNT=',i5)
      if (ieopt .eq. 1)   go to 6 
      plab=EGEV
      s=2.d0*amp**2+2.d0*amp*dsqrt(plab**2+amp**2)
      cms=dsqrt(s)
      go to 8
 6    cms=EGEV
      s  =cms**2
      elab=(cms**2-2.d0*amp**2)/(2.d0*amp)
      plab=dsqrt(elab**2-amp**2)
 8    gamma=0.5d0*cms/amp     ! in colider frame 

      write (6,*) 'plab/GeV, sqrt(s)/GeV=',plab,cms
      anpi=1.5*( 0.88d0+0.44d0*dlog(s)+0.118d0*dlog(s)**2-1.5d0 )  !no. pi/coll
      ang=anpi                     ! number of gluons per NN collision
      d=(4.d0*3.1415926536d0*R0A**3/3.d0)/SIGIN !d=nucleon-nucln sprtn in tube
      dt=2.d0*0.93890595d0*d/dsqrt(s) ! dt=2m(n)d/sqrt(s)  in fm/c

      Acons=0.75+0.38*dlog(cms)
      alcons=3.5+0.7*dlog(cms)
      btcons=0.27+0.037*dlog(cms)
      amt=dsqrt(0.137d0**2+btcons**2)
      ybeam=dlog((plab/amp)+dsqrt((plab/amp)**2+1.d0))
      ycm=ybeam/2.d0
      xp=(amt/amp)*dexp(-ycm)
      dNchdy=Acons*(1-xp)**(2.*alcons)
      dNpidy=1.5d0*dNchdy              ! all pions
      dinitNN=gamma*dNpidy/(d*sigin)

      write (6,*) ' d=',d,' fm, dt=',dt, ' fm/c'
      write (6,*) ' dNchdy(NN) at y=0 dNchdy is',dNchdy
      write (6,*) ' den init (pis/fm3) dinitNN=',dinitNN
      write (6,*) ' transverse mass pi, amt=',amt
      write (6,*) 

      read (5,*)   IBNNwhch,ibxyplt,IPRNk,IPLNk
                 ! IBNNwhch, output npar & ncol, IPR dPdnk or IPL dPdNk if=1
      write(6,*) 'IBNNwhch, for intermediate printout',IBNNwhch
      write(6,*) 'ibxyplt write IBNNwhch case if non-zero',ibxyplt  
                                ! write 2D par,col vs bxy for BNNwhc
      write(6,*) ' IPRNk=',IPRNk,' (printout dPdNk if 1)'
      write(6,*) ' IPLNk=',IPLNk,' (plotout  dPdNk if 1)'
      read (5,*)  IBNNkick0,IBNNdelta  ! get Nkick for IBNNkick0+n*IBNNdelta
      write(6,*)  ' calculate Nkick for IBNNkick0+n*IBNNdelta' 
      write(6,*)  ' IBNNkick0=', IBNNkick0
      write(6,*)  '  and every IBNNdelta=',IBNNdelta
      read (5,*)  phis0deg,dphisdeg,nphismx      ! input in degree  
      write(6,*) ' phis=angle between path and reaction plane (degree)'
      write(6,*) ' phis0  =',phis0deg
      write(6,*) ' dphis  =',dphisdeg
      write(6,*) ' nphismx=',nphismx

      read (5,*)   ridyldanorm,chgmulanorm
      write(6,*) ' ridge -yield-anorm=',ridyldanorm 
      write(6,*) ' chgmul-yield-anorm=',chgmulanorm 
         ! ridyldanorm, so scaledd-ridge-yield = ridyldanorm*(th_yield)
         ! chgmulanorm, so chgmulplicy-yield   = chgmulanorm * totpar (IBNN)
         ! charged multiplicity in CMS cut, e.g. 96, per n-parton-partcipnt 

 703    format('&')

      RA=R0A*A**0.333333D0
      RAG=RA/gamma
      DBA=DBOR*RA
      FACTA=3.D0/(6.283185307D0*RA**3)
      DO 35 IB=1,IBMAX
      BA(IB)=(DFLOAT(IB-1)+0.5d0)*DBA
      if (RA .gt. BA(IB))   TA(IB)=FACTA*dsqrt(RA**2-BA(IB)**2) 
      if (RA .le. BA(IB))   TA(IB)=0.d0         ! sharp cut-off model
   35 CONTINUE
 40   CALL dLINe (BA,TA,IBMAX,dtadb) ! linear extrapolatn 
      BMAXA=BA(IBMAX)
      RB=R0B*B**0.333333D0
      RBG=RB/gamma
      DBB=DBOR*RB
      DBNN=DBNNOR*(RA+RB)
      FACTB=3.D0/(6.283185307D0*RB**3)
      DO 65 IB=1,IBMAX
      BB(IB)=(DFLOAT(IB-1)+0.5d0)*DBB
      if (RB .gt. BB(IB))   TB(IB)=FACTB*dsqrt(RB**2-BB(IB)**2)
      if (RB .le. BB(IB))   TB(IB)=0.d0         ! sharp cut-off model
   65 CONTINUE
 70   CALL dLINe (BB,TB,IBMAX,dtbdb) ! linear extrapolatn 
      BMAXB=BB(IBMAX)

      DO 140 IBNN=1,IBNNMX
      BSEP       =(DFLOAT(IBNN-1)+0.5D0)*DBNN
      BNN  (IBNN)=BSEP
      BNNOR(IBNN)=(DFLOAT(IBNN-1)+0.5D0)*DBNNOR
 140  continue

      IF (IBPRNT .EQ. 0)  GO TO 76
      WRITE(6,*) ' (B, T(B) )  FOR BEAM   NUCLEUS '
      DO 72 IB=1,IBMAX,4
      IBJMX=MIN(IB+3,IBMAX)
      do IBJ=IB,IBJMX
      WRITE(6,656)   BB(IBJ),TB(IBJ)
      enddo
   72 CONTINUE
      WRITE(6,*) ' (B, T(B) )  FOR TARGET NUCLEUS '
      DO 74 IB=1,IBMAX,4
      IBJMX=MIN(IB+3,IBMAX)
      do IBJ=IB,IBJMX
      WRITE(6,656)   BA(IBJ),TA(IBJ)
      enddo
   74 CONTINUE

 76   sigmaAB=3.14159d0*(RA+RB)**2

      WRITE(6,658)   sigmaAB

 85   call parcol        ! parcol gives the participant & collision numbers 

      write (6,*) ' IB   b   part(f)  ovlpA  par/A  part/ovlpA',
     x        ' col(f) col/A col/ovlpA'
      Rab=0.d0
      avpart=0.d0
      do 86 IBNN=1,IBNNMX
      if (totcolf(IBNN) .le. 0.d0)  go to 86      
      costa=(BNN(IBNN)**2+RA**2-RB**2)/(2.d0*BNN(IBNN)*RA)
      thetaA=dacos(costa)
      costb=(BNN(IBNN)**2+RB**2-RA**2)/(2.d0*BNN(IBNN)*RB)
      thetaB=dacos(costb)
      overla=thetaA*RA**2-RA**2*costa*dsqrt(1.d0-costa**2)
     x      +thetaB*RB**2-RB**2*costb*dsqrt(1.d0-costb**2)
         write (6, 666) IBNN,BNN(IBNN),totparf(IBNN),
     x                         overla,dparda(IBNN),totparf(IBNN)/overla,
     x                         totcolf(IBNN),dcolda(IBNN)

      avpart = avpart + totparf(IBNN) * BNN(IBNN) 
      Rab=Rab + BNN(IBNN)*totcolf(IBNN)/SIGIN
 86   continue
      
      avpart=avpart*6.2831853d0*dbNN/sigmaAB ! average number of participants 
      Rab=Rab*6.2831853d0*dbNN/A*B      ! Rab should be A*B if no absorption
 666  format(I2,f7.3,2f8.3,f6.3,f7.3,2f9.3,2f6.3)
       WRITE (6,678) Rab,Rab/(A*B)
 678   FORMAT(/,' Rab=',E12.5,' Rab/A*B=',f10.5,' (if no absrptn,~1)')
       write (6,*) ' #partcp aver over b=',avpart
      go to 1
  999 STOP
CC
  608 FORMAT(2I7)
  620 FORMAT(3X,A4,3X,A4)
  632 FORMAT(' FOR BEAM   B,(FM) R0=',F7.3,
     2 /,    ' FOR TARGET A,(FM) R0=',F7.3)
  640 FORMAT('   DB/RA     =',F10.3,', IBMAX=',I7)
  648 FORMAT(' DBNN/(RA+RB)=',F10.3,', IBNNMX=',I7,' IBNNPR=',I4)
  656 FORMAT(1X,4('(',F5.2,1X,E10.4,') '))
  658 FORMAT(/,' sigma_AB=',f10.3,'  fm**2')
  664 FORMAT(' (',F8.3,2E10.3,')','(',F8.3,2E10.3,')')
      END
CC

      SUBROUTINE parcol
cc    this gives the particpant number, colision number in tube of sigin 
cc    this give particpant number/A and collision numbers/A
cc    this gives aNRidge(B,phis), <average NK>/ntrig(B) averaged over phis 
CC    copy from and modify from INTDB_TAB 
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SPLINT/  y2(100) 
      COMMON /BANDTB/  BB(100),TB(100),DBB,BMAXB,
     2                 BA(100),TA(100),DBA,BMAXA,IBMAX
      common /extrap/  dtadb(100),dtbdb(100)
      COMMON /NULNUL/  BNN(100),BNNOR(100),TAB(100),XI(100),
     2                 DBNN,DBNNOR,IBNNMX,IBPRNT
      COMMON /XSECNN/  B,A,SIGIN
      common /radii /  RA,RB
      common /bcoord/  bbpx(100),bbpy(100),bboa(100,100),bbob(100,100),
     x                 ttoa(100,100),ttob(100,100)        
      common /paramt/  zeta,dnptondnpar,sigjp,t0   
      common /paramb/  bsep        
             ! zeta=exp(-zeta*N),dnptondpar=dNparton/dnparpant,sigjp,t0 
             ! bsep=impact parameter
             ! t0  =hydro start time, before that, all static
      common /surviv/  totcolf(100),totparf(100),
     x                 dcolda (100),dparda (100)
      common /colpar/  ancolxy(100,100),anparxy(100,100),
     x                 btxy   (100,100),atxy   (100,100),
     x                 kcon   (100,100)      ! kont=0/1 for zero/non-zero cont
                                ! We use bb=bb'+b/2, bb' is used as coord.
      common /impact/  IBNNwhch,ibxyplt,IPRNk,IPLNk
      common /IBNNnk/  IBNNkick0,IBNNdelta   ! get IBNNkick0+n*IBNNdelta
      common /Iphis /  phis0deg,dphisdeg,nphismx  ! phis0(deg),dphis,nphismx
      common /Nkdist/  aNktotfis,aNtrgfis,aNcolfis !aNk-all,aNtrg,aNcol[fis]
      common /NRidge/  aNRbfi(100,20),aNtrgbfi(100,20),RAAbfi(100,20),
     >                 aNRb(100),RAAb(100),IBNNkick(100) 
      common /exnorm/  ridyldanorm,chgmulanorm

      do IBBPX =1,IBMAX
         bbpx(ibbpx)=bb(ibbpx)      ! this is just preparing the mesh, bb(fm)
      enddo
      do IBBPY =1,IBMAX
         bbpy(ibbpy)=bb(ibbpy)      ! this is just preparing the mesh, bb(fm)
      enddo  
           ! the coordinates of the source point is bb(ibbpx), bb(ibbpy)
           ! the coords bbpx, bbpy are set up around O, mid-pt bet'n two nuclei
           ! where the pt relative to OA and OB need recalculation as below 
      Rab=0.d0
      iBkik=0   ! start the ankick(iBkik,nphis) for NK(B,phis) array counter
  
      DO 140 IBNN=1,IBNNMX
      BSEP = BNN (IBNN)     ! now get it back
      DO 6 IBBPX =1,IBMAX
        do 5 IBBPY =1,IBMAX
         bbx=bbpx(ibbpx)+0.5d0*bsep         ! bbx=bbpx+bsep/2  in fm
         bby=bbpy(ibbpy)                    ! bby=bbpy
                                            ! bbpx,bbpy relative to O
                                            ! bbx, bby relative to OB (nucl B)
         bob=dsqrt( bbx**2 + bby**2 )       ! Bmag(bx,by) relative to OB  
         bbob(ibbpx,ibbpy)=bob              ! B(bx,by) relative to OB
c
c     In order for interpolation to work,
c     we need to choose bb(ibmax), and ba(ibmax) covers R+6a distances
c      
 1       call dLINt(BB, TB, IBMAX, dtbdb, BOB, TOB)     ! linear extrapolotn
 101     ttob(ibbpx,ibbpy)=tob
         BOA=DSQRT(BOB**2+BSEP**2-2.D0*BSEP*bbx)
         bboa(ibbpx,ibbpy)=boa         
 2       call dLINt(BA, TA, IBMAX, dtadb, BOA, TOA) 
 3       ttoa(ibbpx,ibbpy)=toa
 5      continue
 6    continue

      sumcolf2=0.d0    ! total number in floating pt
      sumparf2=0.d0    ! total number of partcpant (float)
      sumparcolf2=0.d0
      do 100 ibbpx=1,ibmax
      sumcolf1=0.d0
      sumparf1=0.d0
      sumparcolf1=0.d0
      do  70 ibbpy=1,ibmax
      BOB=bbob(ibbpx,ibbpy) 
      TOB=ttob(ibbpx,ibbpy)
      BOA=bboa(ibbpx,ibbpy)
      TOA=ttoa(ibbpx,ibbpy)
      bt=b*tob*sigin
      at=a*toa*sigin
      ancolxy(ibbpx,ibbpy)=0.
      anparxy(ibbpx,ibbpy)=0.
      btxy(ibbpx,ibbpy)=0.
      atxy(ibbpx,ibbpy)=0.
      kcon(ibbpx,ibbpy)=0    ! kcon=0, if content is zero
      IF (TOB .LE. 0.D0)  go to 70
      IF (TOA .LE. 0.D0)  go to 70
 50   ancol=bt*at            ! number of collisions  (floating number)
      anpar=bt+at
      ancolxy(ibbpx,ibbpy)=ancol
      anparxy(ibbpx,ibbpy)=anpar
      kcon   (ibbpx,ibbpy)=1    ! kcon=1, if content is non-zero,
      btxy(ibbpx,ibbpy)=bt
      atxy(ibbpx,ibbpy)=at
      sumcolf1=sumcolf1+ancol
      sumparf1=sumparf1+anpar
      sumparcolf1=sumparcolf1+anpar*ancol  ! flotaing point
   70 CONTINUE
      sumcolf2=sumcolf2+sumcolf1*dbb
      sumparf2=sumparf2+sumparf1*dbb
      sumparcolf2=sumparcolf2+sumparcolf1*dbb   ! floating point   
  100 CONTINUE
      totcolf(IBNN)=sumcolf2*4.D0*DBB/SIGIN ! total number in floating pt
      totparf(IBNN)=sumparf2*4.D0*DBB/SIGIN ! total number of partcpnt (float)
      dcolda (IBNN)=totcolf(IBNN)/sigin     ! average Ncol/unit_area
      if (sumcolf2 .gt. 0.d0) 
     >   dparda (IBNN)=sumparcolf2/(sumcolf2*SIGIN) ! average Npar/unit_area
      if (sumcolf2 .le. 0.d0) dparda(IBNN)=0.d0      

      if ( mod(IBNN,IBNNdelta) .ne. IBNNkick0)  go to 125 
      iBkik=iBkik+1      ! icounter for the ankick(ikcik,nphis) quantities 
cc    get Nkick here

       sumaNR=0.d0
       sumRAA=0.d0
       do 110 nphis=1,nphismx
       phisdeg=phis0deg+(nphis-1)*dphisdeg    ! in degree
       phis=phisdeg*3.14159d0/180.d0          ! in radian

       call getNkick (phis)  ! get aver # kicked partons along phis 
                             ! get aNktotfis=Nk(phis),total # of assoc chrgd Nk
                             ! get aNtrgfis =N_trg(phis),total# of trg

           aNkavg=aNktotfis/aNtrgfis ! average N_kick(charge+neutral)
           aNRavg=aNkavg*2.d0/3.d0       ! average aNRidge, charged only
           aNRbfi  (iBkik,nphis)=aNRavg  
          ! aNRbfi(iBkik,nphis)=(ridge number/trig)(B,phis),aver ovr plane-src
           aNtrgbfi(iBkik,nphis)=aNtrgfis  
          ! aNtrgbfi(B,phis)=number of trig particles(B,phi),aver ovr plane-src
           Raabfi  (iBkik,nphis)=aNtrgbfi(iBkik,nphis)/totcolf(IBNN)

           IBNNkick(iBkik)      =IBNN

           if (IPRNk .ne. 1)  go to 105
           write (6,*)
           chgmul=chgmulanorm*totparf(IBNN)
           write (6,669) BNN(IBNN),chgmul,phisdeg,aNRbfi(iBkik,nphis)
 669       format('  B(fm)=',f8.3,' chgmul=',f10.3,' fi(o)=',f6.1,
     >               'chNR/trg=',f12.5) 
                                    ! quanty(b,fis,(2/3)Nk/trg(b,fis)
           write (6,*) '  trg(B,phis)=',aNtrgbfi(iBkik,nphis) !   ntrig (b,fis)
           write (6,*) '  Raa(B,phis)=',Raabfi  (iBkik,nphis) !   Raa   (b,fis)
           write (6,*) ' totcolf(IBNN),aNcolfis=',totcolf(IBNN),aNcolfis

 105       sumaNR=sumaNR+aNRbfi(iBkik,nphis)
           sumRAA=sumRaa+RAAbfi(iBkik,nphis) 
 110    enddo

       write (9,718) BNN(IBNN)/(RA+RB)
 718   format('# (phis(deg),th-yield#/trig(B,fi)) B/(RA+RB)(B)=', f10.3)
       do nphis=1,nphismx
          phisdeg=phis0deg+(nphis-1)*dphisdeg    ! in degree
          write (9,*) phisdeg,aNRbfi(iBkik,nphis)
       enddo
       write (9,703)
       write (9,7718) BNN(IBNN)/(RA+RB)
 7718  format('# (phis(deg),scl-yield#/trig(B,fi)) B/(RA+RB)(B)=',f10.3)
       do nphis=1,nphismx
          phisdeg=phis0deg+(nphis-1)*dphisdeg    ! in degree
          write (9,*) phisdeg,ryldanorm*aNRbfi(iBkik,nphis)
       enddo
       write (9,703)

       write (11,719) totparf(IBNN)
 719   format('# (phis(deg), RAA(B,fi)) for partpant #(B)=', f10.3)
       do nphis=1,nphismx
          phisdeg=phis0deg+(nphis-1)*dphisdeg    ! in degree
          write (11,*) phisdeg,RAAbfi(iBkik,nphis)
       enddo
       write (11,703)

           aNRb(iBkik)=sumanR/dfloat(nphismx) ! anRidge(B),avg over phis
           RAAb(iBkik)=sumRAA/dfloat(nphismx)

cc    plot & write out B profile here for IBNNwhch
 125  if (IBNN    .ne. IBNNwhch)  go to 135
      if (ibxyplt .ne. 1)         go to 135
      write (6,*) ' DBNN=',DBNN,' fm, IBMAX=',IBMAX
      write (6,*) ' IBNNwhch=',IBNNwhch,' BNNwhchb/(RA+RB)=',BNNwhch
      write (6,*) '   or BNN=',IBNNwhch*DBNN,' fm'
      write (6,*) ' beampar(bx,by) writeout every 5 pts or',5*dbb,' fm'

      do ibbpy=3,ibmax,5     ! to start at 3 to symmetrize output every 5
            write (6,666) (btxy(ibbpx,ibbpy),ibbpx=3,ibmax,5)
      enddo
      write (6,*) ' targetpar(bx,by) writeout evy 5 pts or',5*dbb,' fm'
      do ibbpy=3,ibmax,5
            write (6,666) (atxy(ibbpx,ibbpy),ibbpx=3,ibmax,5)
      enddo
      write (6,*) ' npar(bx,by) writeout every 5 points or',5*dbb,' fm'
      do ibbpy=3,ibmax,5
            write (6,666) (anparxy(ibbpx,ibbpy),ibbpx=3,ibmax,5)
      enddo
      write (6,*) ' ncol(bx,by) writingout every 5 points'
      do ibbpy=3,ibmax,5
            write (6,666) (ancolxy(ibbpx,ibbpy),ibbpx=3,ibmax,5)
      enddo
 666  format(10f7.2)

 135  continue
 140  CONTINUE

        if (IPLNk .ne. 1 )  go to 170
        write (8,728) 
 728    format('# (participant#(B), theory-yield=chrgN_ridge/ntrig(B))')
      iBkikmx=iBkik
      do iBkik=1,iBkikmx
         IBNNindx=IBNNkick(iBkik)
         write (6,*) 
         chgmul = chgmulanorm * totparf(IBNNindx)
         ridyld = ridyldanorm * aNRb(iBkik)
         write (6,687) BNN(IBNNindx),totparf(IBNNindx),chgmul
 687     format(' B(fm)=',f8.3,' Npart=',f10.3,' chgmul=',f12.4)
         write (6,688) aNRb(iBkik),ridyld
 688     format('       NRb/Ntrg=',f12.4,' scl-ridydl/Ntrg=',f12.4)
         write (6,*) ' Ra(B)=',RAAb(iBkik)
      enddo
       write ( 8,703)
        write (8,7728) 
 7728 format('# (chgmul#(B),scl-yield=yldanrom*chgN_rdge/ntrg(B))')
      do iBkik=1,iBkikmx
         IBNNindx=IBNNkick(iBkik)
         chgmul=chgmulanorm*totparf(IBNNindx)  ! chgmul is scled-chged multplcy
         ridyld=ridyldanorm*aNRb(iBkik)        ! ridyld is scled ridge yield
         write(8,*) chgmul,ridyld
      enddo
       write ( 8,703)
        write (8,7730) 
 7730   format('# (chgmul#(B),(chgN+neut)_rdge/ntrg(B))')
      do iBkik=1,iBkikmx
         IBNNindx=IBNNkick(iBkik)
         chgmul=chgmulanorm*totparf(IBNNindx)  ! chgmul is scled-chged multplcy
         write(8,*) chgmul,(3./2.)*aNRb(iBkik) ! mult vs #kicked partn(chg+neu)
      enddo
       write ( 8,703)

       write (10,729) 
       write ( 6,729) 
 729  format('# (centrality(%)(B), RAA(B)=Ntrig/tot coll number(B))')
      do iBkik=1,iBkikmx
         IBNNindx=IBNNkick(iBkik)
         centrality=100.d0*BNN(IBNNindx)**2/(RA+RB)**2    ! centrality in %
         write ( 6,*) centrality,RAAb(iBkik)
         write (10,*) centrality,RAAb(iBkik)
      enddo
       write (10,703)
       write (10,730) 
 730   format('# (participant#(B), RAA(B)=Ntrig/tot coll number(B))')
      do iBkik=1,iBkikmx
         IBNNindx=IBNNkick(iBkik)
         chgmul=chgmulanorm*totparf(IBNNindx)
         write (10,*) chgmul,RAAb(iBkik)   ! chgmul is scled-chged multplcy
      enddo
       write (10,703)

 703   format('&')

 170   RETURN                  
      END


      subroutine getNkick (phis)     
cc    get dP/dN_kick as a function of phi_s (angle wrt reaction plane)
cc    get average Nkick
cc    output is in dPdNk vs aNk
                               ! get aver # kicked partons along phis 
      IMPLICIT REAL*8(A-H,O-Z)
      common /Nkdist/  aNktotfis,aNtrgfis,aNcolfis !aNk-all,aNtrg,aNcol[fis]
      COMMON /XSECNN/  B,A,SIGIN
      COMMON /BANDTB/  BB(100),TB(100),DBB,BMAXB,
     2                 BA(100),TA(100),DBA,BMAXA,IBMAX
      common /radii /  RA,RB
      common /bcoord/  bbpx(100),bbpy(100),bboa(100,100),bbob(100,100),
     x                 ttoa(100,100),ttob(100,100)        
      common /paramt/  zeta,dnptondnpar,sigjp,t0   
      common /paramb/  bsep        
             ! zeta=exp(-zeta*N),dnptondpar=dNparton/dnparpant,sigmjp,t0
             ! bsep=impact parameter
      common /surviv/  totcolf(100),totparf(100),
     x                 dcolda (100),dparda (100)
      common /colpar/  ancolxy(100,100),anparxy(100,100),
     x                 btxy   (100,100),atxy   (100,100),
     x                 kcon   (100,100)      ! kont=0/1 for zero/non-zero cont
      cs=1.d0/dsqrt(3.d0)       ! speed of sound cs = 1/sqrt(3)
      cssq=0.33333333d0

      dL=DBB          ! choose to increment dL same as DBB (in fm)
      sumcol  =0.d0
      sumtrg  =0.d0
      sumNktot=0.d0
      do 70 ibbpx= ibmax,-ibmax,-1
         if (ibbpx .eq. 0)  go to 70    ! the index is 1 to ibmax only 
         iaapx=iabs(ibbpx)
         bsrcx=bbpx(iaapx)      ! |x| position of source of jet
      do 60 ibbpy=-ibmax,ibmax
         if (ibbpy .eq. 0)  go to 60    ! the index is 1 to ibmax only 
         iaapy=iabs(ibbpy)
         bsrcy=bbpy(iaapy)      ! |y| position of source of jet
         if (kcon(iaapx,iaapy) .eq. 0)  go to 60   ! no contribution if zero

         call interp(bbpx,dbb,ibmax,bbpy,dbb,ibmax,ancolxy,
     >               bsrcx,bsrcy,ancolq)

         if (ancolq .le. 0.d0)  go to 60   ! skip if no jet source

         sum=0.d0
         aL=t0+0.5d0*dL         ! only after medium is produced can jet be absr
 10      qlx=isign(1,ibbpx)*bsrcx+aL*dcos(phis) 
                                 !    x-position of Q(vec)=src(vec)+al(vec)
         qly=isign(1,ibbpy)*bsrcy+aL*dsin(phis) 
                                 !    y-position of Q(vec)=src(vec)+al(vec)    
        bp=dsqrt(qlx**2+qly**2)
        phip=datan(qly/qlx)
        tt=aL            ! tt=time variable,trjtry,measured after jet&zexpstart
                         ! t0=zr-hydrostart=const, measured after jet&zexpstart
        tR=tt-t0         ! tR=TRiemann=varbl=tt-t0 measured aftr  r-hydrostart
 20      bbqx=dabs(qlx)         ! find the x of P=|QL| to get parcpant#
         bbqy=dabs(qly)          ! find the y of P=|QL| to get parcpant#
         fbt=(t0/tt)**(3.d0*cssq)
 40      call interp(bbpx,dbb,ibmax,bbpy,dbb,ibmax,anparxy,
     >               bbqx,bbqy,anparq)
         if (anparq .le.  1.d-10)  go to 55
         sum=sum+fbt*anparq/sigin
         aL=aL+dL
         go to 10

 50      continue
 55      ankxy=sum*dL*dnptondnpar*sigjp/(2.d0*t0) ! sigjp is kept in fm**2
                   ! ankxy=N_kickpartonnumber(bbpx,bbpy; phis)
                   ! nb: factor of two in deonminator, left and right movement
                   ! dnptondnpar=dnparton/dnpparpnt
                   ! sigjp      =sigma(jet-parton collision) in fm**2
                   ! r-hydrostar initial time t0, mesured wrt z-hydro time tt 
                   ! now get dncoldb at source point

         dncoldb=ancolq/sigin    ! per unit area of fm**2

         weight=dncoldb*dbb*dbb                 ! dbxy=dbnn*dbnn area
         ! put ankxy in bin with weight                     
         sumcol  =sumcol  +weight               ! binary coliision number
         sumtrg  =sumtrg  +weight*      dexp(-zeta*ankxy)
         sumNktot=sumNktot+weight*ankxy*dexp(-zeta*ankxy)
 60   continue 
 70   continue
      if (sumtrg .gt. 0.d0)  go to 80
      aNcolfis =0.d0
      aNtrgfis =0.d0
      aNktotfis=0.d0
      return
 80   aNcolfis =sumcol        ! Ncollision(B,fis), B is understood
      aNtrgfis =sumtrg        ! Ntrigger  (B,fis)
      aNktotfis=sumNktot      ! Nktot     (B,fis)= toal(chrg+neutral)(B,fis)
      return
      end


      subroutine dLINe (BA, TA, IBMAX, dtadb)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension   BA(100),TA(100),dtadb(100)

      do i=1,ibmax-1
         dtadb(i)=(TA(i+1)-TA(i))/(BA(i+1)-BA(i))
      enddo
      return
      end

      subroutine dLINt (BA, TA, IBMAX, dtadb, boa, toa)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension   BA(100),TA(100),dtadb(100)
      if (boa .lt. BA(IBMAX))  go to 20
      toa=0.d0
      return
 20   if (boa .ge. BA(1)    )  go to 30
      toa=TA(1)
      return
 30   iL=1
      do i=2,ibmax
         if (BA(i) .ge. boa) go to 40
         iL=i
      enddo
 40   iR=iL+1
      toa=ta(iL)+(boa-ba(iL))*dtadb(iL)
      return
      end


      SUBROUTINE interp(xpts,dx,ixmax,ypts,dy,iymax,tabxy,x,y,f)
cc    from Eq. 25.2.56 of NBS table Abramowitz and Stegun
cc    symmetrical,  F(x,y)=F(|x|,|y|)  in our case

      IMPLICIT REAL*8(A-H,O-Z)
      dimension xpts(100),ypts(100),tabxy(100,100)

      ix=dabs(x)/dx+0.50001d0
      if (ix .ge. ixmax)  go to 10
      iy=dabs(y)/dy+0.50001d0
      if (iy .ge. iymax)  go to 10
      if (ix .le. 1)  ix=1
      if (iy .le. 1)  iy=1
      p=(x-xpts(ix))/dx
      q=(y-ypts(iy))/dy
      f=(1.d0-p)*(1.d0-q)*tabxy(ix,iy  )+p*(1.d0-q)*tabxy(ix+1,iy)
     >        +q*(1.d0-p)*tabxy(ix,iy+1)+p*q*tabxy(ix+1,iy+1)
      return
 10   f=0.d0
      return
      end

