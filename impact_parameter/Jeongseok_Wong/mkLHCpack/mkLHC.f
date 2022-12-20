cc    mkLHC.f   based on momkick02.f, with no major modifications
cc    momkick02.f   version 01,  July, 2009
cc        
cc    to calculate the distribtuion of ridge particles
cc    based on (1) ridge are medium partons kicked by the jet
cc             (2) ridge particle distribution depends on: q,fRNk
cc             (3) total chargd associated=(chrgd ridge yield)+fJ*(ppjet yield)
cc                  q   =momentum kick
cc                  fRNk=attenuatn fR*(# of total kicked med partons) in AA
cc                  fJ  =attenuation factor for ppjet in AA collision
cc   initial parton distribution specified by parameters: am,aa,Ti,amd
cc   initial ppjet yield specified by parameters:anjet,Tjet,sigfi0,sigeta0,ama

cc    when using this program, please refer to 
cc             C. Y. Wong, Phys. Rev. C 78, 064905 (2008),
cc        and  C. Y. Wong, Phys. Rev. C 80, 034908 (2009).
cc    fortran program written by 
cc             Dr. Cheuk-Yin Wong 
cc             Oak Ridge National Laboratory
cc             Oak Ridge, TN 37831-5373
cc             Email:  wongc@ornl.gov
cc   To compile:   f95 -o mkLHC.x mkLHC.f 

      IMPLICIT REAL*8 (A-H,O-Z)

      dimension Peta (100,100,100),Petaj(100,100,100)   !Peta=ridge,Petaj=ppjet
         ! Peta (Deta,Dfi,pt)=dF/dDeta dDfi pt dpt= normalize ridge distrbtn
         ! Peta is normalized to 1; int dDeta dDfi ptdpt  Peta(...) = 1
         ! Petaj(Deta,Dfi,pt) is ppjet yield = dNch_ppjet/dDeta dDfi pt dpt
         ! Petaj is normalized to anjet = num of ppjet ass charged particles 
      dimension fiv  (100),zetav(100),ptv(100)          !fi,zeta,pt of ass part
      dimension wetaj(100),wzeta(100),facc(100),wfi(100) !facc=f_accptnc-crrtn 
      dimension dNptdpt(100),dNptdptj(100),dNptdpttot(100)
      dimension dNdpt  (100),dNdptj  (100),dNdpttot  (100)
      dimension dNdfi  (100),dNdfij  (100),dNdfitot  (100)
      dimension dNdeta (100),dNdetaj (100),dNdetatot (100)
      dimension dNdetadfi(100,100)    ! integrated over pt*dpt

      open ( 5,file='mkLHC.05',status='unknown')
      open ( 6,file='mkLHC.06',status='unknown') 
      open ( 7,file='mkLHC.07',status='unknown') !plt dN/ptdpt,dN/dpt,etc
      open ( 8,file='mkLHC.08',status='unknown') !plt dN/dfi
      open ( 9,file='mkLHC.09',status='unknown') !plt dN/deta(eta),R,J,J+R
      open (10,file='mkLHC.10',status='unknown') !plt dNtot(dfi,deta) 

 1     read (5,*,end=999)  izetaCase,iaccCor  
            ! the meaning of the zeta variable depends on izetaCase
            ! If izetaCase=0, zeta =etaAss,iaccCor is not used(e.g.PHENIX-Jia)
            ! If izetaCase=1, zeta =etaAss-etajet=Delta_eta(e.g.STAR)
            ! If izetaCase=1 and iaccCor=1,use tringlr Delta_eta acptnce crrctn
            ! If izetaCase=1,and iaccCor=0,not use tringl Delta_eta acpt crrctn
       write (6,*) 
       write (6,*) ' EXPERIMENTAL SETUP & ETA ACCEPTANCE INPUT'
       write (6,*) ' Input choice: izetaCase =',izetaCase
       write (6,*) ' Input choice: iaccCor   =',iaccCor
       write (6,*) ' The meaning of zeta depends on izetaCase'
       write (6,*) ' Cases are:'
       write (6,*) ' If izetaCase=0,zeta=etaAss & iaccCor is not used '
       write (6,*) ' If izetaCase=1,zeta=etaAss-etajet=Delta_eta(STAR)'
       write (6,*) ' If izetaCase=1 &iaccCor=0, no tnglr accpt corrctns'
       write (6,*) ' If izetaCase=1 &iaccCor=1,use tnglr accpt corrctns'
       read (5,*)  etajmn,etajmx,netaj,nsdetaj
            ! etajetmn= eta jet min; etajmx=eta jet max;
            ! netaj   = number of etajet
            ! nsdetaj = # of sides of the etaj (e.g. 1 side or 2 sides) 
       write(6,*)  '(Jet) etajmn,etajmx,netaj=',etajmn,etajmx,netaj
       write(6,*)  '# of above sides of (Jet),nsdetaj=',nsdetaj

       read (5,*)  etaAssmn, etaAssmx
       write(6,*)  '(Associated) etaAssmn,etaAssmx=',etaAssmn,etaAssmx

       read (5,*)  zetamn,zetamx,nzeta,nsdzeta !zeta=what? depnds on izetaCase:
       write(6,*)  '(zetaRange)zetamn,zetamx,nzeta=',zetamn,zetamx,nzeta 
       write(6,*)   '# of sides of above zeta set  =',nsdzeta

      read (5,*,end=999)  q, fRNk, fJ  ! momentum kick condition parameters
cc     q   = magnitude of momentum kick per jet-medium-parton kick in GeV
cc     fRNk=fR*<N_k>=ridge_attenuation fR * number_of kicked partons
cc     fJ  =attenuation factor,to multiply ppjet to get AA jetcomponent/trigger
cc       implicitly assume:  (total ridge particle) = number_of kicked partons
cc       also assume:observd ridge part=charged= (2/3)*(total ridge particles)
       write(6,*) 
       write(6,*) 'INPUT MOMENTUM KICK PARAMETERS:'
       write(6,*) ' q=',q,'GeV = momentum kick in jet directn(GeV)'
       write(6,*) ' fRNk=',fRNk,'= fRNk = fR * <total Nkick >'
       frNkch=0.6666667d0*fRNk    
       write(6,*) ' fRNkch=',fRNkch,'=(2/3)*fR*<total Nk>'
       write(6,*) '                =totl chrgd ridge particles'
       write(6,*) ' fJ=',fJ,'=survval factor of pp jet in AA collisns'

       read (5,*)  am,aa,Ti,amd    ! normalized initial distribution parameters
cc     nromalized initial momentum distribution is
cc      dF/dyptdpt = A_norm (1-x)^aa * exp(-amT/Ti)/amdT; 
cc                 amT =sqrt(am^2 +pt^2)
       write(6,*) 
       write(6,*) 'INPUT INITIAL PARTON MOMENTUM DISTRBTN PARAMETERS:'
       write(6,*) 'am=',am,'GeV '
       write(6,*) ' initial y-dis=(1-x)**aa' 
       write(6,*) 'aa=',aa,'= expon index of (1-x)**aa init distrbtn'
       write(6,*) 'Ti=',Ti,' GeV = initial Tempertaure'    
       write(6,*) 'amd=',amd,' GeV = mass in denominator of distrbtn'


       read (5,*)  sqrts, PTtrig          ! sqrt(s_NN), average of PTtrig
       write(6,*) 
       write(6,*) 'EXPERIMENTAL sqrt(s) and PT_trig'
       write(6,*)  ' sqrt(s_NN) (GeV)=', sqrts
       emrat=sqrts/0.93890595d0
       yb=dlog(emrat+dsqrt(emrat**2-1.d0))
       write(6,*)  ' ybeam=-ytarget=',yb
       write(6,*)  ' PT_trig (GeV)=', PTtrig 

       read (5,*)  aNjet0, dNjetdPT, Tjet0, dTjetdPT
       read (5,*)  sigfi0, sigeta0, ama  
            ! ppjet distribtn parameters
            ! aNjet,Tjet =# and T of jet partons,(integrated ovr eta, fi)
            ! aNjet = aNjet0 + dNjetdPT * PTtrig
            ! Tjet  = Tjet0  + dTjetdPT * PTtrig
            ! sigfi0 = fi  wdth of jet partons,(intg ovr eta, pt sectn)
            ! sigeta0= eta wdth of jet partons,(intg ovr eta, pt sectn)
            ! ama for sigfi=sigfi0*ama/sqrt(ama**2+pt**2), same with sigeta
       write(6,*) 
       write(6,*) 'INPUT ppjet MOMENTUM DISTRIBUTION PARAMETERS'
       write(6,*) 'aNjet0=',aNjet0
       write(6,*) 'Tjet0 =',Tjet0
       write(6,*) 'dNjetdPT=',dNjetdPT
       write(6,*) 'dTjetdPT=',dTjetdPT
       aNJet = aNjet0 + dNjetdPT * PTtrig
       TJet  = Tjet0  + dTjetdPT * PTtrig
       write(6,*) 
       write(6,*) 'aNjet = aNjet0 + dNjetdPT * PTtrig'
       write(6,*) 'TJet  = Tjet0  + dTjetdPT * PTtrig'
       write(6,*) 'aNjet=',aNjet,',totl chrged partl asscatd with ppjet' 
       write(6,*) 'Tjet =',Tjet, ' GeV= Tempertaure of jet component'
       write(6,*) 'sigfi =sigfi0 * ama /sqrt(ama^2+pt**2)'
            !        sigfi =sigfi0 * ama /sqrt(ama**2 +pt**2)
       write(6,*) 'sigeta =sigeta0 * ama /sqrt(ama^2+pt**2)'
            !        sigeta =siget0 * ama /sqrt(ama**2 +pt**2)
       write(6,*) ' sigfi0=',sigfi0,'  sigeta0=',sigeta0
       write(6,*) 'ama=',ama,' GeV= mass ama in above formula'

       write(6,*) '-------------------------------------------------'
       write(6,*) ' BRIEF SUMMARY OF FORMULAS'
       write(6,*) ' dNtotal=dNridge_component + dNAAjet_component'
       write(6,*) ' dNridge/ptdptdetadphi=(2/3)*fRNk*dF/ptdptdetadphi'
       write(6,*) ' dNAAjet/ptdptdetadphi=   fJ * dNppjet/ptdptdetadphi'
       write(6,*) ' fRNk  =fR * < Nk >=ridgeatten*tot num kiked parton' 
       write(6,*) ' fRNkch=(2/3)fR*<Nk>=ridatten*totchrgnkiked parton' 
       write(6,*) '   observed(charged) ridge particle = fRNkch' 
       write(6,*) 
       write(6,*) ' Normlz init (chrgd and neutrl) momentum dstrbtn is'
       write(6,*) ' dF/dyptdpt = A_norm (1-x)^aa * exp(-amT/Ti)/amdT'
       write(6,*) '                 amT =sqrt(am^2 +pt^2)'
       write(6,*) '                 amdT=sqrt(amd^2+pt^2)'
       write(6,*) '             x=amT*exp(|y|-yb)/amb; yb=ybeam, amb=am'
       write(6,*) 
       write(6,*) ' pp jt associated (charged) yield is ' 
       write(6,*) 'dNpp/dDeta dDfi ptdpt=anjet*ANormj*exp(-amT/Tjet))'
       write(6,*) '    *exp{-Delta_fi^2 /[2*sigfi^2 ]}'
       write(6,*) '    *exp{-Delta_eta^2/[2*sigeta^2]}'
       write(6,*) '    /[2*pi*sigfi*sigeta]'
       write(6,*) '    where sigfi =sigfi0 *ama/sqrt(ama**2+pt**2)'    
       write(6,*) '          sigeta=sigeta0*ama/sqrt(ama**2+pt**2)'    
       write(6,*) 'dNchtot/...= fRNkch*dF/... + fJ*dNpp/...'
       write(6,*) '--------------------------------------------------'
       write(6,*) 

       write(6,*) 'EXPERIMENTAL CONDITIONS AND ACCEPTANCE'

       if (netaj .gt. 1)   go to  15
       detaj=0.d0
       wetaj(1)=1.d0*nsdetaj
       aNumetaj=1.d0                    !aNumetaj=numof etaj jets tobesumedover
       go to 20
 15    detaj=(etajmx-etajmn)/dfloat(netaj-1)
          wetaj    (1)=0.5d0*nsdetaj           ! weights for etaj integration
       do izetaj=2,netaj-1
          wetaj(izetaj)= 1.d0*nsdetaj
       enddo
          wetaj(netaj)=0.5d0*nsdetaj
       aNumetaj=dfloat(netaj-1)*nsdetaj !aNumetaj=numof etaj jets tobesumedover

 20    write(6,*)  '(Jet)               detaj=',detaj


         ! if izetaCase=1, eta=Deta=etaAss-etajet or etaAss=eta+etajet
         ! if izetaCase=0, eta=etaAss             or etaAss=eta
         ! nsdzeta = 1 only 1 zeta(positive) side included, =2 both side incld

       if (nzeta .gt. 1)   go to  25
       dzeta=0.d0
       zetav(1)=zetamn
       wzeta(1)=1.d0*nsdzeta
       facc(1)=1.d0            ! normally facc is one unless we correct for it
       go to 30

 25    dzeta=(zetamx-zetamn)/dfloat(nzeta-1)
          do izeta=1,nzeta
            zeta=zetamn+(izeta-1)*dzeta
            zetav(izeta)=zeta
            wzeta(izeta)=1.0d0*nsdzeta
            facc(izeta)=1.d0    ! normally facc is one unless we correct for it
          enddo
       wzeta   (1)=0.5d0*nsdzeta          ! weights for zetaj integration
       wzeta(nzeta)=0.5d0*nsdzeta

       if (iacCcor  .eq. 0)  go to 30  ! when izetaCase=1,choice to crrct/ornot 
c      calcualte acceptance correcetion to be used when (izetaCase=1 and iacccor=1)
       x1mn=etaAssmn
       x1mx=etaAssmx
cc       x2mn=etajmn
cc       x2mx=etajmx
       x2mn=x1mn             ! we should choose etaAssmx greater than etajmx
       m2mx=x1mx             ! we should choose|etaAssmn| greater than |etajmn|
       deltaxmn=x1mn-x2mx    ! deltax=x1-x2
       deltaxmx=x1mx-x2mn
          do 28 izeta=1,nzeta
             deltax=zetav(izeta)
            if ((deltax.le.deltaxmn).or.(deltax.ge.deltaxmx))  go to 998
            if (deltax .gt. (x1mn-x2mn) ) go to 26
               facc(izeta)=(x2mx-x2mn)/(deltax-(x1mn-x2mx))
            go to 28
 26         if (deltax .gt. (x1mx-x2mx) ) go to 27
               facc(izeta)=1.d0
            go to 28
 27            facc(izeta)=(x2mx-x2mn)/((x1mx-x2mn)-deltax)
 28      continue
 30    write (6,*)  '(Output zeta)        dzeta=', dzeta

       read (5,*)  fimin, dfi, nfi, nsdfi  ! final particle phi
       write (6,*)  'fimin, dfi, nfi =', fimin, dfi, nfi    
       write(6,*)   '# of above zeta sides  nsdfi=', nsdfi

       if (nfi .gt. 1)   go to  35
       wfi(1)=1.d0*nsdfi
       go to 40
 35    wfi  (1)=0.5d0*nsdfi           ! weights for fij integration
       do ifi=2,nfi-1
          wfi(ifi)=1.0d0*nsdfi
       enddo
          wfi(nfi)=0.5d0*nsdfi

 40    read (5,*)   ptmin, dpt, npt ! final particle pt 
       write(6,*)  'ptmin, dpt, npt =', ptmin, dpt, npt
       write(6,*)
       read (5,*)   pttrmin, dpttr, npttr ! final trigger particle pt 
       write(6,*)  'trgger: ptmin, dpt, npt =', pttrmin, dpttr, npttr

       read (5,*)  iptplt,ifiplt,ietaplt,ifietaplt
                  ! if any of the above=0, no plot
                  ! if 1 for iptplt,ifiplt,ietaplt, then plot integrated plot
                  ! if ifiplt=1, plot dNch/ptdpt (ppjet,AAridge,AAjet+ridge) 
                  ! if ifiplt=2, plot dNch/dpt   (ppjet,AAridge,AAjet+ridge) 
                  ! if ifietaplt=0 no plot of yield(phi,eta)
                  ! if ifietaplt=1 plot dn/dphideta(phi,zeta), fi=x, zeta=y
                  ! if ifietaplt=2 plot dn/dphidets(phi,zeta), zeta=x, fi=y
       write (6,*)  'iptplt=   ',iptplt    ! if=1 plot yield vs pt, 
       write (6,*)  'ifiplt=   ',ifiplt    ! if=1 plot yield vs phi,
       write (6,*)  'ietaplt=  ',ietaplt   ! if=1 plot yield vs zeta,
       write (6,*)  'ifietaplt=',ifietaplt ! if=1 get yield vs fi-eta,@e.g.pt=2GeV?
       write (7,701)
       write (8,701)
       write (9,701)
       write (10,701)
 701   format('@type xy')
 703   format('&') 
       dyy=yb/100.d0
       sum=0.d0
       do ipt=1,400
          pt=dfloat(ipt-1)*0.01
          amti=dsqrt(am**2+pt**2)
          do 55 iy=1,100
          yy=0.d0+(dfloat(iy)-0.5)*dyy
          xx=(amti/am)*dexp(yy-yb)
          if (xx .gt. 1.d0)  go to 55
          term=(1.d0-xx)**aa*dexp(-amti/Ti)/dsqrt(amd**2+pt**2)
          sum=sum+term*pt
 55    continue
       enddo
       sum=6.2831853*sum*0.01d0*dyy*2.d0     ! factor of 2 for both sides

       aNorm=1.d0/sum     ! this is A_norm, normalizatn for initial mom dstrbtn

       aNormj=1.d0/( Tjet * (am+Tjet) * dexp(-am/Tjet) )  ! normlizn con ppjet

        do 100 izeta=1,nzeta
            zeta=zetav(izeta)

            do 90  ifi=1,nfi
               fi=fimin+(ifi-1)*dfi
               fiv(ifi)=fi

               do 80  ipt=1,npt
                  pt=ptmin+(ipt-1)*dpt
                  ptv(ipt)=pt
                                    ! these are final p
                                    ! to beam3-z,1&2 are Trnsv-perp
                  sigfi =sigfi0 *ama /dsqrt(ama **2+pt**2)
                  sigeta=sigeta0*ama /dsqrt(ama **2+pt**2)
                  amtj=dsqrt(am**2+pt**2)

                  sumdndeta =0.d0
                  sumdndetaj=0.d0

               do 70 ietaj=1,netaj   ! average over all etajet distributn
                  etaj=etajmn+(ietaj-1)*detaj    ! eta of jet trigger
                  if (izetaCase .eq. 1)  go to 62 
                  etaAss=zeta        ! now izetaCase .eq. 0 then  zeta=output=etaAss
                  DelEta=etaass-etaj ! 
                  go to 65
 62               DelEta=zeta        ! now izetaCase .eq. 1 then  zeta=output=etaAss-etaj
                  etaAss=DelEta+etaj !                            etaAss=DelEta+etaj
 65               if (etaAss .lt. etaAssmn ) go to 70
                  if (etaAss .gt. etaAssmx ) go to 70
                                     ! etaAss is within etaAss acceptance
                  p3=pt*dsinh(etaAss)  ! p in 3 (z-direction, beam) uses etaAss
                  p2=pt*dsin(fi)       ! p in 2 (y-direction, beamXjet)
                  p1=pt*dcos(fi)       ! p in 1 (x-direction, jet)
                  amt=dsqrt(am**2+p1**2+p2**2)
                  Ef =dsqrt(am**2+p1**2+p2**2+p3**2)
                  y  =0.5d0*dlog((Ef+p3)/(Ef-p3))
                              ! p(final)=p=p(init)+q, p1i i for initial
                p1i=p1-q/dcosh(etaj)   ! etaj is eta of jet
                p2i=p2
                p3i=p3-q*dsinh(etaj)/dcosh(etaj)
                amti=dsqrt(am**2+p1i**2+p2i**2)  ! to beam3,1&2 are Trnsv-perp
                amtiden=dsqrt(amd**2+p1i**2+p2i**2)     ! denominator amd 
                Ei  =dsqrt(am**2+p1i**2+p2i**2+p3i**2)
                yi  =0.5d0*dlog((Ei+p3i)/(Ei-p3i)) 
                xx=(amti/am)*dexp(dabs(yi)-yb)
                if (xx .ge. 1.d0)  go to 70
                termy=aNorm*(1.d0-xx)**aa      ! normalization const aNorm used
     >               *(dexp(-amti/Ti)/amtiden)
     >               *Ef/Ei
                termeta =termy*dsqrt(1.d0-am**2/(amt**2*dcosh(y)**2))
                sumdndeta =sumdndeta +wetaj(ietaj)*termeta  ! multplied by weight

cc       evaluate ppjet distribution below
                termetaj=aNjet * aNormj * dexp(-amtj/Tjet)
     >               *dexp(-fi **2/(2.d0*sigfi **2 ))
     >               *dexp(-DelEta**2/(2.d0*sigeta**2 ))
     >               /(6.2831853d0*sigfi*sigeta)
                sumdndetaj=sumdndetaj+wetaj(ietaj)*termetaj  ! multplied by weight

 70             continue

                sumdndeta =sumdndeta /aNumetaj  ! assume etajet=uniform
                sumdndetaj=sumdndetaj/aNumetaj  ! assume etajet=uniform
                                        !aNumetaj=numof etaj jets tobesumedover
         Peta (izeta,ifi,ipt)=sumdndeta   ! normlized ridge dF/deta dphi pt dpt
                                          ! averaged over incident eta_jet
         Petaj(izeta,ifi,ipt)=sumdndetaj  ! jet dN/deta ddfi pt dpt
         if (izetaCase .eq. 0)  go to 80
         if (iaccCor   .eq. 0)  go to 80
         Peta (izeta,ifi,ipt)=Peta (izeta,ifi,ipt)*facc(izeta)
         Petaj(izeta,ifi,ipt)=Petaj(izeta,ifi,ipt)*facc(izeta)
 80      continue
 90      continue
         write (6,*) 
 100     continue

cc  ----------------------------------------------------
cc   following is for plotting dN/ptdpt or dN/dpt 

       if (iptplt  .eq. 0)   go to 105    ! no pt plotting if iptplt is 0
       do ipt=1,npt
          pt=ptmin+(ipt-1)*dpt
          sumdNptdpt =0.d0
          sumdNptdptj=0.d0
          do izeta=1,nzeta
            do ifi=1,nfi
      sumdNptdpt =sumdNptdpt +wzeta(izeta)*wfi(ifi)*Peta (izeta,ifi,ipt)
      sumdNptdptj=sumdNptdptj+wzeta(izeta)*wfi(ifi)*Petaj(izeta,ifi,ipt) 
            enddo
          enddo
          sumtm  =sumdNptdpt *dzeta*dfi ! normalized      dN/pt dpt
          sumtmj =sumdNptdptj*dzeta*dfi ! jet      dN/pt dpt
          dNptdpt   (ipt)=fRNkch*sumtm    ! observed AA ridge dNch/pt dpt 
          dNptdptj  (ipt)=       sumtmj    ! charged ppjet dNchpp/pt dpt 
          dNptdpttot(ipt)=dNptdpt(ipt)+fJ*dNptdptj(ipt) 
                                          ! observed AA ridge dNch/pt dpt 
          dNdpt   (ipt)=pt*dNptdpt(ipt) ! observed AA ridge dNch/dpt 
          dNdptj  (ipt)=pt*dNdpt  (ipt) ! observed ppjet dNch/dpt 
          dNdpttot(ipt)=pt*dNptdpttot(ipt) ! observed AA ridge dNch/dpt 
       enddo

       if (iptplt  .eq. 2)   go to 103
       write (7,710)
 710   format('#  ppjet:    [pt, dNch/ptdpt(pt)] ')
       do ipt=1,npt
          write (7,*) ptv(ipt),dNptdptj(ipt)  ! plot out jet pp dN/ptdpt
       enddo
       write (7,703)
       write (7,711)
 711   format('#  AA-ridge: [pt, dNch/ptdpt(pt)] ' ) 
       do ipt=1,npt
         write (7,*) ptv(ipt),dNptdpt(ipt) ! plot out integrated yields dN/ptdpt
       enddo
       write (7,703)
       write (7,712)
 712   format('#  AA jet+ridge: [pt, dNch/ptdpt(pt)]' )
       do ipt=1,npt
          write (7,*) ptv(ipt),dNptdpttot(ipt)
       enddo
       write (7,703)
       go to 105

 103   write (7,713)
 713   format('#  pppjet:    [pt, dNch/  dpt(pt)] ')
       do ipt=1,npt
          write (7,*) ptv(ipt),dNdptj(ipt)  ! plot out jet pp dN/ptdpt
       enddo
       write (7,703)
       write (7,714)
 714   format('#  AA-ridge: [pt, dNch/  dpt(pt)] ' ) 
       do ipt=1,npt
        write (7,*) ptv(ipt),dNdpt(ipt) ! plot out integrated yields dN/ptdpt
       enddo
       write (7,703)
       write (7,715)
 715   format('#  AA jet+ridge: [pt, dNch/  dpt(pt)]' )
       do ipt=1,npt
          write (7,*) ptv(ipt),dNdpttot(ipt)
       enddo
       write (7,703)

cc  ----------------------------------------------------
cc   following is for plotting dNch/dDelta_fi

 105   if (ifiplt  .eq. 0)   go to 130   ! no plot
       if (ifiplt  .eq. 1)   go to 110   ! plot intgd zeta, pt; J+R

 110  do ifi=1,nfi
          fi=fimin+(ifi-1)*dfi
          sumdNdfi =0.d0
          sumdNdfij=0.d0
          do izeta=1,nzeta
           do ipt=1,npt
            pt=ptmin+(ipt-1)*dpt
            term =wzeta(izeta)*pt*Peta (izeta,ifi,ipt)
            termj=wzeta(izeta)*pt*Petaj(izeta,ifi,ipt)
            if (ipt .ne. 1)   go to 111  
            term =0.5d0*term     ! the weight for ipt=1 is 0.5
            termj=0.5d0*termj     ! the weight for ipt=1 is 0.5
            go to 112
 111        if (ipt .ne. npt) go to 112
            term =0.5d0*term     ! the weight for ipt=1 is 0.5
            termj=0.5d0*termj     ! the weight for ipt=1 is 0.5
 112        sumdNdfi =sumdNdfi +term
            sumdNdfij=sumdNdfij+termj
           enddo
          enddo
          sumtm     =sumdNdfi *dzeta*dpt ! dN/dfi
          sumtmj    =sumdNdfij*dzeta*dpt ! dN/dfi
          dNdfi (ifi)=fRNKch*sumtm      ! AA charged ridge dNdfi
          dNdfij(ifi)=       sumtmj     ! ppjet dNdfi 
          dNdfitot(ifi)=dNdfi(ifi)+fJ*dNdfij(ifi) ! AA jet+ridge
       enddo
       write (8,720)   zetamn, zetamx 
 720   format('# ppjet: [fi, dNch/dfi], zetamn=',F10.3,' zetamx=',f10.3) 
       if (nsdfi.eq.1)  go to 115
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfij(ifi)  !plot negative part of phi
       enddo
 115   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfij(ifi)  !plot positive part of phi
       enddo
       write (8,703)

       write (8,721) 
 721   format('# AA ridge: [fi,  dNch/dfi]') 
       if (nsdfi.eq.1)  go to 117
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfi(ifi)  !plot negative part of phi
       enddo
 117   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfi(ifi)  !plot positive part of phi
       enddo
       write (8,703)

       write (8,722)
 722   format('# AA jet+ridge: [fi,  dNch/dfi]')
       if (nsdfi.eq.1)  go to 119
       do ifi=nfi,2,-1
          write (8,*) -fiv(ifi),dNdfitot(ifi)  !plot negative part of phi
       enddo
 119   do ifi=1,nfi
          write (8,*)  fiv(ifi),dNdfitot(ifi)  !plot positive part of phi
       enddo
       write (8,703)

cc  ----------------------------------------------------
cc   following is for plotting dNch/dDelta_eta

 130   if (ietaplt .eq. 0)  go to 160  !no dn/deta plot
       if (ietaplt .eq. 1)  go to 140  !plot dn/deta, intgtd fi, pt; J+R 

cc     plot dN/dzeta here, integrated over fi, pt 
 140   do izeta=1,nzeta
          sumdNdeta =0.d0
          sumdNdetaj=0.d0
          do ifi=1,nfi
           do ipt=1,npt
            pt=ptmin+(ipt-1)*dpt
            term =wfi(ifi)*pt*Peta (izeta,ifi,ipt)
            termj=wfi(ifi)*pt*Petaj(izeta,ifi,ipt)
            if (ipt .ne. 1)   go to 141  
            term =0.5d0*term     ! the weight for ipt=1 is 0.5
            termj=0.5d0*termj     ! the weight for ipt=1 is 0.5
            go to 142
 141        if (ipt .ne. npt) go to 142
            term =0.5d0*term     ! the weight for ipt=1 is 0.5
            termj=0.5d0*termj     ! the weight for ipt=1 is 0.5
 142        sumdNdeta =sumdNdeta +term
            sumdNdetaj=sumdNdetaj+termj
           enddo
          enddo
          sumtm     =sumdNdeta *dfi*dpt ! dN/dfi
          sumtmj    =sumdNdetaj*dfi*dpt ! dN/dfi
 144      dNdeta (izeta)=fRNkch *sumtm  ! AA ridge (charged) 
          dNdetaj(izeta)=        sumtmj ! pp jet charged 
          dNdetatot(izeta)=dNdeta(izeta)+fJ*dNdetaj(izeta) !add the two incl fJ 
       enddo

       write (9,740)  iaccCor 
       write (6,740)  iaccCor 
 740   format('#  ppjet:[Deta,dNch/dDeta],iaccCor(1yes/0no)=',I2)
       if (nsdzeta.eq.1)  go to 146
       do izeta=nzeta,2,-1
          write (9,*) -zetav(izeta),dNdetaj(izeta)  !plot negative part of eta
          write (6,*) -zetav(izeta),dNdetaj(izeta)  !plot negative part of eta
       enddo
 146   do izeta=1,nzeta
          write (9,*)  zetav(izeta),dNdetaj(izeta)  !plot positive part of zeta
          write (6,*)  zetav(izeta),dNdetaj(izeta)  !plot positive part of zeta
       enddo
       write (9,703)
       write (6,*) 

       write (9,742)  iaccCor 
       write (6,742)  iaccCor 
 742   format('# AAridge:[Deta,dNch/dDeta],iaccCor(1yes/0no)=',I2)
       if (nsdzeta.eq.1)  go to 145
       do izeta=nzeta,2,-1
          write (9,*) -zetav(izeta),dNdeta(izeta)  !plot negative part of zeta
          write (6,*) -zetav(izeta),dNdeta(izeta)  !plot negative part of zeta
          sum=sum+dNdeta(izeta)*dzeta
       enddo
 145   sum=0.d0
       do izeta=1,nzeta
          write (9,*)  zetav(izeta),dNdeta(izeta)  !plot positive part of zeta
          write (6,*)  zetav(izeta),dNdeta(izeta)  !plot positive part of zeta
          sum=sum+dNdeta(izeta)*dzeta
       enddo
       write (9,703)
       write (6,*) ' 1-side chg-rid-yld from zetamn->zetamx=', sum
       write (6,*) ' 2-side chg-rid-yld from zetamn->zetamx=', 2.d0*sum
       write (6,*) 

       write (9,744)  iaccCor 
       write (6,744)  iaccCor 
 744   format('# AAjet+ridge:[Deta,dNch/dDeta],iaccCor(1yes/0no)=',I2)
       if (nsdzeta.eq.1)  go to 147
       do izeta=nzeta,2,-1
          write (9,*) -zetav(izeta),dNdetatot(izeta) !plot negative part of zeta
          write (6,*) -zetav(izeta),dNdetatot(izeta) !plot negative part of zeta
       enddo
 147   do izeta=1,nzeta
          write (9,*)  zetav(izeta),dNdetatot(izeta) !plot positive part of eta
          write (6,*)  zetav(izeta),dNdetatot(izeta) !plot positive part of eta
       enddo
       write (9,703)

cc  ------------------------------------------------------------
cc   following is for plotting 3D plot dNch/dDelta_fi dDelta_eta


 160   if (ifietaplt  .eq. 0)   go to 170  


cc     plot  fi=x, zeta=y, yield=z

 165   do izeta=1,nzeta
            do ifi=1,nfi
                  sum=0.d0
               do ipt=1,npt 
                  term=fRNkch * Peta (izeta,ifi,ipt)*ptv(ipt)*dpt
     >                +fJ     * Petaj(izeta,ifi,ipt)*ptv(ipt)*dpt 
                  sum=sum+term
               enddo
            dNdetadfi(izeta,ifi)=sum
            enddo
                 ! sum = int_ptmn^ptmx (pt dpt) dNtot/[d(Deta)d(Dphi) pt dpt]
       enddo

       do izeta=nzeta,2,-1
          do ifi=nfi,2,-1
          sum=dNdetadfi(izeta,ifi)
          if (ifietaplt.eq.1) write (10,*) -fiv(ifi),-zetav(izeta),sum
          if (ifietaplt.eq.2) write (10,*) -zetav(izeta),-fiv(ifi),sum
          enddo
          do ifi=1,nfi
          sum=dNdetadfi(izeta,ifi)
          if (ifietaplt.eq.1) write (10,*)  fiv(ifi),-zetav(izeta),sum
          if (ifietaplt.eq.2) write (10,*) -zetav(izeta), fiv(ifi),sum
          enddo
          write (10,*)    ! blank line important to separate lines
       enddo

       do izeta=1,nzeta
          do ifi=nfi,2,-1
          sum=dNdetadfi(izeta,ifi)
          if (ifietaplt.eq.1) write (10,*) -fiv(ifi),zetav(izeta),sum
          if (ifietaplt.eq.2) write (10,*) zetav(izeta),-fiv(ifi),sum
          enddo
          do ifi=1,nfi
          sum=dNdetadfi(izeta,ifi)
          if (ifietaplt.eq.1) write (10,*)  fiv(ifi),zetav(izeta),sum
          if (ifietaplt.eq.2) write (10,*) zetav(izeta), fiv(ifi),sum
          enddo
          write (10,*)    ! blank line important to separate lines
       enddo

 170   continue
       go to 1      ! read another set of data cards   
 998   write (6,*) ' Delta-eta=etaAss-etajet is out of range, stop'
 999   stop
       end

