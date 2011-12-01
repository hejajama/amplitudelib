      subroutine dlib(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi) 

************************************************************************
************************************************************************
* Stefan Kretzer (kretzer@pa.msu.edu) Nov 2001:                        *
*                                                                      *
* Library of FF parametrizations updated to be posted on               *
* http://www.pv.infn.it/~radici/FFdatabase/                            *
* site maintained by Marco Radici (radici@pv.infn.it)                  *
* and Rainer Jakob (rainer@theorie.physik.uni-wuppertal.de)            *
* within the EU network:                                               *
* "Hadronic Physics with High Energy Electromagnetic Probes"           *
*                                                                      *
* Also posted on that site are the individual FF sets included in this *
* library as distributed by the corresponding authors.                 *
************************************************************************
* K, KKP and BFGW fragmentation functions for                          *
* iparton -> ihadron                                                   *
*                                                                      *
* Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        *
* "dff(iparton)" of dimension dff(-5:5)                                *
*                                                                      *
* iparton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar      *
*                                                                      *
* [Q2] = GeV^2                                                         *
*                                                                      *
* ffset = 1,2,3 means K, KKP, BFGW                                     *
* Details are to be found in the corresponding references:             *
* S.Kretzer, Phys.Rev.D62, 054001 (2000)                               *
* B.A.Kniehl, G.Kramer, B.Potter, Nucl.Phys.B582, 514 (2000)           *
* L.Bourhis, M.Fontannaz, J.P.Guillet, M.Werlen,Eur.Phys.J.C19,89(2001)*
*                                                                      *
* fforder=0,1 is LO, NLO(MSbar)                                        *
* no LO set for BFGW                                                   *
*                                                                      *
* ihadron=1,2,3,4,5 is pi,K,h,p,n                                      *
* no pi,K sets for BFGW                                                *
* no p,n sets for K, BFGW                                              *
*                                                                      *
* icharge=0,1,2,3 is 0,+,-,+&-                                         *
* note: 3 = +&- = charge sum (NOT average) in this library             *
* for the neutral (icharge=0) particles K^0, n:                        *
* icp = 1,2,3 chooses between particle, anti-particle or sum of both   *
* icp is inactive for pi,h,p and charged Kaons                         *
*                                                                      *
* Following KKP, FFs into {pi^0;(anti-)K^0} and into {pi^+/-;K^+/-}    *
* are related by isospin.                                              *
*                                                                      *
* ipi = 1,2,3 is a flag for BFGW; inactive for K, KKP:                 *
* ipi = 1: best fit     (formula (8))  of BFGW                         *
* ipi = 2: large Ng set (formula (9))  of BFGW                         *
* ipi = 3: low Ng set   (formula (10)) of BFGW                         *
************************************************************************
* If options (ihadron, icharge etc.) are chosen which do not exist for *
* a given parametrization (ffset) or if a flag is chosen outside its   *
* range as defined above then the code returns a corresponding warning *
* and stops.                                                           *
************************************************************************
************************************************************************

      implicit double precision ( a - z )
      dimension dff(-5:5)  
      dimension uqff(2), dqff(2), sqff(2), cqff(2), bqff(2)
      dimension dh(0:10)
      integer ffset, fforder, iparton, ipi, ipisk, ifini, old1, 
     >        iset, ihadron, icharge, old2, old3, ih, isu2, icp
      character*12 readme
      COMMON / FRAGINI / IFINI
      common / order / old1
      common / old / old2, old3
      data readme/ '; see MANUAL' /
      ipisk = ipi  
      
      if ( .not. ( (fforder .eq. 0).or.(fforder .eq. 1) ) ) then
      write(6,*) 'choose fforder=0,1 for LO,NLO'//readme
      stop
      endif      

      if ( (ihadron .gt. 3 ) .and. (.not. (ffset .eq. 2)) ) then
      write(6,*) 'no p,n FFs (ihadron=4,5) for K, BFGW (ffset=1,3)'
     > //readme
      stop
      endif

      if ( ihadron .lt. 1 .or. ihadron .gt. 5 ) then 
      write(6,*)
     > 'set ihadron=1,2,3,4,5 for pi,K,h,p,n'//readme
      stop
      endif    

      if ( (icharge .eq. 0 ) .and. 
     >   (  ihadron .eq. 3 .or. ihadron .eq. 4 ) ) then
       write(6,*) 'incl hadrons (ihadron=3) or protons (ihadron=4)'
     > //' must be charged (icharge=1,2,3)'//readme
      stop
      endif 

      if ( icharge .lt. 0 .or. icharge .gt. 3 ) then
      write(6,*) 'set icharge=0,1,2,3 for 0,+,-,+&-'//readme
      stop
      endif    

      if (
     > ( (ihadron .eq. 2 .or. ihadron .eq. 5) .and. icharge .eq. 0 )
     >  .and. ( icp .lt. 1 .or. icp .gt. 3 )  ) then
      write(6,*) 'for neutral (icharge=0) K,n (ihadron=2,5) choose'
     >//' particle, anti-particle or sum (icp=1,2,3)'//readme 
      stop
      endif

      if ( ffset .eq. 1 ) then

      norm = 1.d0
      isu2 = 0
 
      if ( .not. ( fforder  .eq. old1 ) .or.
     >     .not. ( ihadron  .eq. old2 ) .or.
     >     .not. ( icharge  .eq. old3 )
     >   ) ifini = 0
      old1 = fforder
      old2 = ihadron
      old3 = icharge

      iset = 2*ihadron + fforder -1
      if ( ihadron .eq. 1 .and. icharge .eq. 0 ) then
      norm = 0.5d0
      icharge = 3
      endif
      if ( ihadron .eq. 2 .and. icharge .eq. 0 ) then
      icharge = icp
      isu2 = 1
      endif
 
      call PKHFF(ISET,ICHARGE,Z,Q2,uqff,dqff,sqff,cqff,bqff,gff)

      dff(-5     ) = bqff(2) * norm 
      dff(-4     ) = cqff(2) * norm
      dff(-3     ) = sqff(2) * norm
      dff(-2+isu2) = dqff(2) * norm
      dff(-1-isu2) = uqff(2) * norm
      dff( 0     ) = gff     * norm
      dff( 1+isu2) = uqff(1) * norm
      dff( 2-isu2) = dqff(1) * norm
      dff( 3     ) = sqff(1) * norm
      dff( 4     ) = cqff(1) * norm
      dff( 5     ) = bqff(1) * norm
      return

      elseif ( ffset .eq. 2 ) then

      norm = 2.d0

      qs = dsqrt(Q2)
      if ( icharge .eq. 3 ) then
      if ( ihadron .eq. 5 ) then
      write(6,*)
     > 'set icharge=0 for neutrons (ihadron=5)'//readme
      stop
      endif
      ih = ihadron
      if ( ihadron .eq. 3) then
      ih = 7
      norm = 1.d0
      endif
      elseif ( icharge .eq. 0 ) then
       if ( ihadron .eq. 1 ) then
       ih = 5
       norm = 1.d0
       elseif ( ihadron .eq. 2 ) then
       if ( .not. icp .eq. 3 ) then 
       write(6,*)
     > 'only sum of particle & anti-particle (icp=3)'
     > //'for KKP K^0 set'//readme
       stop
       endif 
       ih = 3
       elseif ( ihadron .eq. 5 ) then
       if ( .not. icp .eq. 3 ) then
       write(6,*)
     > 'only sum of particle & anti-particle (icp=3)'
     > //'for KKP neutron set'//readme
       stop
       endif 
       ih = 6
       else
       write(6,*) 
     >'only neutral (icharge=0) pions (ihadron=1), Kaons (ihadron=2)'
     > //' and neutrons (ihadron=5)'//readme
       stop
       endif
      else
      write(6,*) 
     >'icharge must be 0 (neutral) or 3 (charge sum) for KKP (ffset=2)'
     >//readme
      stop
      endif

      call kkp(ih,fforder,z,qs,dh)
      dff(-5) = dh(10) * norm
      dff(-4) = dh( 8) * norm
      dff(-3) = dh( 6) * norm
      dff(-2) = dh( 4) * norm
      dff(-1) = dh( 2) * norm
      dff( 0) = dh( 0) * norm
      dff( 1) = dh( 1) * norm
      dff( 2) = dh( 3) * norm
      dff( 3) = dh( 5) * norm
      dff( 4) = dh( 7) * norm
      dff( 5) = dh( 9) * norm
      return

      elseif ( ffset .eq. 3) then

      if (fforder .eq. 0) then
      write(6,*) 'only NLO (fforder=1) set for BFGW (ffset=3)'//readme
      stop
      endif

      if (.not. ihadron .eq. 3) then
      write(6,*) 'only incl hadron set (ihadron=3) for BFGW (ffset=3)'
     > //readme
      stop
      endif      

      if (.not. icharge .eq. 3) then
      write(6,*) 'only charge sum (icharge=3) for BFGW (ffset=3)'
     > //readme
      stop
      endif

      if ( ipisk .lt. 1 .or. ipisk .gt. 3 ) then
      write(6,*)
     > 'set ipi=1,2,3 for central,large,small Ng set of BFGW (ffset=3)'
     > //readme
      stop
      endif

      call fonfra(z,ipisk,Q2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     > ,xdbp,xdbbp,xdgp)
      dff(-5) = xdbbp / z
      dff(-4) = xdcp  / z
      dff(-3) = xdsp  / z
      dff(-2) = xddbp / z
      dff(-1) = xdubp / z
      dff( 0) = xdgp  / z
      dff( 1) = xdup  / z
      dff( 2) = xddp  / z
      dff( 3) = xdsp  / z
      dff( 4) = xdcp  / z
      dff( 5) = xdbp  / z
      return

      else
 
      write(6,*) 'choose ffset=1,2,3 for K,KKP,BFGW'//readme
      stop      

      endif

      return

      end

