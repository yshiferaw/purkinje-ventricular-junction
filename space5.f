c     ================================================================
c     CARDIAC ELECTROPHYSIOLOGY SIMULATION
c     double precision code - organized for readability

c     ifx space5.f jun3.f avcurrents7.f acur6.f
c     ================================================================

      implicit double precision (a-h,o-z)

c     ================================================================
c     PARAMETER DEFINITIONS
c     ================================================================
      parameter (Lx=100,Ly=6)      ! horizontal long strip
      parameter (Lx1=100,Ly1=52)    ! 2D sheet
      parameter(nstim=2)

c     ================================================================
c     SIMULATION CONTROL ARRAYS
c     ================================================================
      double precision tstim(0:nstim), tstimx(0:nstim)
c      double precision mask(110,60)
c      maskca(400,400)
c      double precision zmask(1:Lx+Lx1),zmask2(1:Lx+Lx1),zmask3(1:Lx+Lx1)
     
       double precision zmask1(1:Lx+Lx1), zmask2(1:Lx+Lx1)
       double precision zmask3(1:Lx+Lx1), zmask4(1:Lx+Lx1)
       double precision zmask5(1:Lx+Lx1), zmask6(1:Lx+Lx1)
       double precision zmask7(1:Lx+Lx1)
      
      character*15 filenm1(500)
      character*15 filenm2(500)
      character*15 filenm3(500)
      character*15 filenm4(500)
      character*15 filenm5(500)
      character*15 filenm6(500)	
      character*15 filenm7(500)	

c     ================================================================
c     TISSUE 1 VARIABLE DECLARATIONS
c     ================================================================
      double precision v(0:Lx+1,0:Ly+1), vnew(0:Lx+1,0:Ly+1)
      double precision dv(0:Lx+1,0:Ly+1),vold(0:Lx+1,0:Ly+1)

      double precision xm(Lx,Ly), xh(Lx,Ly), xj(Lx,Ly), xr(Lx,Ly)
      double precision xs1(Lx,Ly), qks(Lx,Ly), xkur(Lx,Ly), ykur(Lx,Ly)
      double precision xtof(Lx,Ly), ytof(Lx,Ly)

      double precision cb(Lx,Ly),ci(Lx,Ly),csrb(Lx,Ly),csri(Lx,Ly)
      double precision cnsr(Lx,Ly), po(Lx,Ly), c1(Lx,Ly), c2(Lx,Ly)
      double precision xi1(Lx,Ly), xi2(Lx,Ly), xi2s(Lx,Ly),c1s(Lx,Ly)
      double precision c2s(Lx,Ly), xi1s(Lx,Ly), cit(Lx,Ly), cbt(Lx,Ly)
      double precision pi(Lx,Ly), pb(Lx,Ly), pox(Lx,Ly), pos(Lx,Ly)
      double precision ra(Lx,Ly)
      double precision poxx(Lx,Ly), pixx(Lx,Ly)
      
      integer nsb(Lx,Ly), nsi(Lx,Ly)
     
c     ================================================================
c     TISSUE 2 VARIABLE DECLARATIONS
c     ================================================================
      double precision vz(0:Lx1+1,0:Ly1+1), vnewz(0:Lx1+1,0:Ly1+1)
      double precision dvz(0:Lx1,0:Ly1),voldz(0:Lx1,0:Ly1)

      double precision xmz(Lx1,Ly1), xhz(Lx1,Ly1), xjz(Lx1,Ly1)
      double precision xs1z(Lx1,Ly1), qksz(Lx1,Ly1), xkurz(Lx1,Ly1)
      double precision ykurz(Lx1,Ly1), xrz(Lx1,Ly1)
      double precision xtofz(Lx1,Ly1), ytofz(Lx1,Ly1)
  
      double precision cbz(Lx1,Ly1),ciz(Lx1,Ly1),csrbz(Lx1,Ly1)
      double precision csriz(Lx1,Ly1), cnsrz(Lx1,Ly1)
      double precision poz(Lx1,Ly1), c1z(Lx1,Ly1), c2z(Lx1,Ly1)
      double precision xi1z(Lx1,Ly1), xi2z(Lx1,Ly1), xi2sz(Lx1,Ly1)
      double precision c1sz(Lx1,Ly1)
      double precision c2sz(Lx1,Ly1), xi1sz(Lx1,Ly1), citz(Lx1,Ly1)
      double precision cbtz(Lx1,Ly1)
      double precision piz(Lx1,Ly1), pbz(Lx1,Ly1), poxz(Lx1,Ly1)
      double precision posz(Lx1,Ly1)
      double precision raz(Lx1,Ly1)
      integer nsbz(Lx1,Ly1), nsiz(Lx1,Ly1)
     
c     ! ***************************************************
     
      double precision wna(Lx,Ly),wica(Lx,Ly)
      double precision wnaz(Lx1,Ly1),wicaz(Lx1,Ly1)
      double precision wika(Lx,Ly), wikaz(Lx1,Ly1)
      
      double precision poa(Lx,Ly), poa_inf(Lx,Ly)
      double precision poaz(Lx1,Ly1), poaz_inf(Lx1,Ly1)
     
c     ================================================================
c     OPEN OUTPUT FILES
c     ================================================================

     
        open(unit=20,file='vs10.dat',status='unknown')
        open(unit=21,file='vs100.dat',status='unknown')
        open(unit=22,file='vs190.dat',status='unknown')
        
        open(unit=30,file='icas10.dat',status='unknown')
  
 
       call generate_filenamesv(filenm1)
       call generate_filenamesk(filenm2)
       call generate_filenamesca(filenm3)
       call generate_filenamesna(filenm4)
       call generate_filenamespo(filenm5)
       call generate_filenamespof(filenm6)
       call generate_filenamespofx(filenm7)
 	
 
 
c     ================================================================
c     STIMULATION PROTOCOL
c     ================================================================


      bcl=640.0

      duration =dfloat(nstim)*bcl
            
      tstim(0)=0.0
      tstimx(0)=1.0
      
      do i=1,nstim
      tstim(i)=tstim(i-1)+bcl
      tstimx(i)=tstim(i)+1.0d0
      enddo 

      dt=0.1d0/1.50
      nstep = duration/dt
      
      t=0.0
     
      kku=1
     
      nbt=1000  ! total number of channels in cells in the strip
      nit=6000

      nbtz=1000 ! total number of channels in cells in the square tissue
      nitz=6000
      
c     ================================================================
c     INITIALIZE ARRAYS
c     ================================================================
c     tissue 1
      call initCa1(Lx,Ly,cb,ci,csrb,csri,cnsr,po,c1,c2,xi1,
     #    xi2,ra,pos,c1s,c2s,xi1s,xi2s,cit,cbt,nsb,nsi)

     
c     tissue 2
      call initCa1(Lx1,Ly1,cbz,ciz,csrbz,csriz, cnsrz,
     #  poz,c1z,c2z,xi1z,xi2z,raz,posz,c1sz,c2sz,xi1sz,
     #  xi2sz,citz,cbtz,nsbz,nsiz)

   
c    tissue 1 and 2 share the same voltage initialization 

	 call initV(Lx,Ly,v,xm,xh,xj,xr,xs1,qks,xkur,
     #  ykur,xtof,ytof,vold)

        call initV(Lx1,Ly1,vz,xmz,xhz,xjz,xrz,xs1z,qksz,
     #  xkurz,ykurz,xtofz,ytofz,voldz)

c     ================================================================
c     BIOPHYSICAL PARAMETERS
c     ================================================================

         xnao=136.0d0           ! external Na (mM)
         xki=140.0d0            ! internal K (mM)
         xko=5.40d0             ! external K (mM)
         cao=1.8d0              ! external Ca (mM)

         temp=308.0d0           ! temperature (K)
         xxr=8.314d0            ! gas constant
         xf=96.485d0            ! Faraday's constant
         frt=xf/(xxr*temp)

         xmx=-2.0/250.
         xnai=xmx*bcl+16.0

	 aleak=0.001d0 
 
 	 gica=1.500
 	 gnaca=0.70
 	 pbx=0.3
 
	 gicaz=1.8
	 gnacaz=0.7d0
	 pbxz=0.3
          
c *****************************************************************        
         
        ! volumes for tissue 1  
        vi=1.0  
	vb=0.3
	vbi=vb/vi

	vq=30.0
	visr=vq
	vsrin=(1.0/visr)*vi

	vbsr=vq
	vsrbound=(1.0/vbsr)*vb

	vbisr=(vsrbound/vsrin)

	vnsr=vq
	vnsrin=(1.0/vnsr)*vi	
	
c  --------------------------------------------------------------
	  ! volumes for tissue 1  & 2
        viz=1.0  
	vbz=0.3
	vbiz=vbz/viz

	vqz=30.0
	visrz=vqz
	vsrinz=(1.0/visrz)*viz

	vbsrz=vqz
	vsrboundz=(1.0/vbsrz)*vbz

	vbisrz=(vsrboundz/vsrinz)

	vnsrz=vqz
	vnsrinz=(1.0/vnsrz)*viz	
	
            
c     ================================================================
c     MAIN INTEGRATION LOOP
c     ================================================================
       do 1 ncount = 0, nstep
       
         time = dfloat(ncount)*dt

	    mv=mod(time,bcl)
            if (mv.eq.0) then
            rux=t/bcl
            iz=int(rux)
            endif 

c        ============================================================
c        TISSUE 1 PROCESSING (horizontal strip)
c        ============================================================
           do  iy=1,Ly
            do  ix=1,Lx

	! add heterogeneity here

c              ================================================
c              CALCIUM HANDLING
c              ================================================
c              convert total Ca to free
               call xfree(cit(ix,iy),ci(ix,iy))
               call xfree(cbt(ix,iy),cb(ix,iy))

c              fraction of clusters with sparks
               pi(ix,iy)=dfloat(nsi(ix,iy))/dfloat(nit)
               pb(ix,iy)=dfloat(nsb(ix,iy))/dfloat(nbt)

c              uptake currents
	
	       vupb=0.1d0
	       vupi=0.1d0
               call uptake(cb(ix,iy),vupb,xupb)
               call uptake(ci(ix,iy),vupi,xupi)

c              ================================================
c              Na-Ca EXCHANGE AND L-TYPE Ca CURRENT
c              ================================================
               call inaca(v(ix,iy),frt,xnai,xnao,cao,cb(ix,iy),xinacaq)
              
               xinacaq=gnaca*xinacaq

               pox(ix,iy)=po(ix,iy)+pos(ix,iy)
               call ica(v(ix,iy),frt,cao,cb(ix,iy),pox(ix,iy),rca,xicaq)
               xicaq=gica*140.0*xicaq

c              ================================================
c              CALCIUM SPARKS AT BOUNDARY
c              ================================================
              
                ab=50.0*0.3
		csrx=850.0d0
		phisr=1.0/(1.0+(csrx/csrb(ix,iy))**10) ! cutoff rate bellow 500 muM

		alphab=ab*dabs(rca)*po(ix,iy)*phisr ! spark on rate due to LCC 
		bts=1.0/30.0   ! spark off rate 

c              update LCC Markov states
               call markov(dt,v(ix,iy),cb(ix,iy),c1(ix,iy),
     +          c2(ix,iy),xi1(ix,iy),xi2(ix,iy),po(ix,iy),
     +          c1s(ix,iy),c2s(ix,iy),
     +          xi1s(ix,iy), xi2s(ix,iy),pos(ix,iy),alphab,bts)

c	      poxx(ix,iy)=po(ix,iy)+pos(ix,iy)+c1(ix,iy)+c1s(ix,iy)
c              xlv=xi1(ix,iy)+xi1s(ix,iy)
c              +xi1s(ix,iy)+xi2s(ix,iy))
c              pixx(ix,iy)=xlv
              
c              boundary RyR current
               gsrb=0.01d0/2.5  
               xryrb=gsrb*csrb(ix,iy)*pb(ix,iy)

c              ================================================
c              INTERIOR CALCIUM SPARKS
c              ================================================
               
               	csrxx=900.0d0
	      	chx=850.0d0
		phi=(csri(ix,iy)/chx)/(1.0+(csrxx/csri(ix,iy))**10) ! cutoff rate bellow 1000 muM
	
	! activation of sparks

		gi=0.25
		pra=0.05  ! increase threshold here
		xra=(pra/pi(ix,iy))**6
		ra_inf=phi/(1.0+xra)

		dra=(ra_inf-ra(ix,iy))/60.0
	
		rr=0.02*3.0*1.5
	
		pbinf=1.0/(1.0+(pbx/pb(ix,iy))**10)

		ar1=(aleak/500.0)*(csri(ix,iy)/1000.0)
		ar2=rr*pbinf*phi
		ar3=gi*ra(ix,iy)

		alphai=ar1+ar2+ar3

		gryri=0.015
		xryri=gryri*pi(ix,iy)*csri(ix,iy)

		btsi=1.0/70.0
          
c              ================================================
c              CALCIUM CONCENTRATION UPDATES
c              ================================================
         	
	xsarc=-xicaq+xinacaq
	
	tausri=50.0
	dff=(cnsr(ix,iy)-csri(ix,iy))/tausri

	tau1=10.0
	tau2=10.0

	dfbi=(cb(ix,iy)-ci(ix,iy))/tau1
	dfbisr=(csrb(ix,iy)-cnsr(ix,iy))/tau2

	dcbt=xryrb-xupb+xsarc-dfbi
	dcsrb=vbsr*(-xryrb+xupb)-dfbisr   

	dcit=xryri-xupi+vbi*dfbi               	
	dcsri=visr*(-xryri)+dff
	
	dnsr=vnsr*(xupi)-dff+vbisr*dfbisr

	cbt(ix,iy)=cbt(ix,iy)+dcbt*dt
	cit(ix,iy)=cit(ix,iy)+dcit*dt
	
	csrb(ix,iy)=csrb(ix,iy)+dcsrb*dt	
	csri(ix,iy)=csri(ix,iy)+dcsri*dt

	cnsr(ix,iy)=cnsr(ix,iy)+dnsr*dt
	
	ra(ix,iy)=ra(ix,iy)+dra*dt

            
c              ================================================
c              STOCHASTIC SPARK EVOLUTION
c              ================================================
               nsbx=nsb(ix,iy)
               call binevol(nbt,nsbx,alphab,bts,dt,ndeltapx,
     &              ndeltamx)

               if(ndeltamx.gt.nsbx.or.ndeltapx.gt.nbt) then
                 nsb(ix,iy)=0
               else
                 nsb(ix,iy)=nsb(ix,iy)+ndeltapx-ndeltamx
               endif

               nsix=nsi(ix,iy)
               call binevol(nit,nsix,alphai,btsi,dt,ndeltapy,
     &              ndeltamy)
               
               if(ndeltamy.gt.nsix.or.ndeltapy.gt.nit) then
                 nsi(ix,iy)=0
               else
                 nsi(ix,iy)=nsi(ix,iy)+ndeltapy-ndeltamy
               endif

c              ================================================
c              VOLTAGE DYNAMICS WITH ADAPTIVE TIME STEP
c              ================================================
               wca=12.0
               xinaca=wca*xinacaq
               xica=2.0d0*wca*xicaq

c              time step adjustment
               adq=dabs(dv(ix,iy))
               if(adq.gt.25.0d0) then
                 mstp=10
               else
                 mstp=1
               endif
               hode=dt/dfloat(mstp)

               do iii=1 , mstp

c                voltage dependent currents
                 call ina(hode,v(ix,iy),frt,xh(ix,iy),xj(ix,iy),
     &            xm(ix,iy),xnai,xnao,xina,pon,pon_inf)

                 call ikr(hode,v(ix,iy),frt,xko,xki,
     &            xr(ix,iy),xikr)

                 call iks(hode,v(ix,iy),frt,xnao,xnai,xko,xki,
     &            xs1(ix,iy),qks(ix,iy),xiks)

                 call ik1(hode,v(ix,iy),frt,xki,xko,xik1)

                 call ikur(hode,v(ix,iy),frt,xki,xko,xkur(ix,iy),
     &            ykur(ix,iy),xikur)

                 call ito(hode,v(ix,iy),frt,xki,xko,xtof(ix,iy),
     &            ytof(ix,iy),xito)

                 call inak(v(ix,iy),frt,xko,xnao,xnai,xinak)
  
  
  		    stim=0.0d0 
		 if(time.ge.tstim(iz).and.time.le.tstimx(iz)) then
		 
		 if(ix.le.10) then 
                       stim=80.0d0      
                     else
                       stim=0.0
                 endif         
                 
                      endif

	  if(iz.eq.2) then
                stim=0.0
                endif 
                
	 xkzz=xik1+xikr+xiks+xito+xikur+xinak
                
          wna(ix,iy)=-xina*1000.0  ! record the sodium current here
	  wica(ix,iy)=-xica 
	  wika(ix,iy)=xkzz     
	  
	 poa(ix,iy)=pon
	 poa_inf(ix,iy)=pon_inf
 	
            dvh=-(xina+xkzz+xinaca+xica)+stim
     
                 v(ix,iy)=v(ix,iy)+dvh*hode

                enddo
                

             end do
         end do
         
         
         ! do updating of tissue 2
         
       
           do  iy=1, Ly1
            do  ix=1, Lx1

c              ================================================
c              CALCIUM HANDLING
c              ================================================
c              convert total Ca to free
               call xfree(citz(ix,iy),ciz(ix,iy))
               call xfree(cbtz(ix,iy),cbz(ix,iy))

c              fraction of clusters with sparks
               piz(ix,iy)=dfloat(nsiz(ix,iy))/dfloat(nitz)
               pbz(ix,iy)=dfloat(nsbz(ix,iy))/dfloat(nbtz)

c              uptake currents
	
	       vupb=0.1d0
	       vupi=0.1d0
               call uptake(cbz(ix,iy),vupb,xupb)
               call uptake(ciz(ix,iy),vupi,xupi)

c              ================================================
c              Na-Ca EXCHANGE AND L-TYPE Ca CURRENT
c              ================================================
          call inaca(vz(ix,iy),frt,xnai,xnao,cao,cbz(ix,iy),xinacaq)
              
               xinacaq=gnacaz*xinacaq

               poxz(ix,iy)=poz(ix,iy)+posz(ix,iy)
          call ica(vz(ix,iy),frt,cao,cbz(ix,iy),poxz(ix,iy),rca,xicaq)
               xicaq=gicaz*140.0*xicaq

c              ================================================
c              CALCIUM SPARKS AT BOUNDARY
c              ================================================
              
                ab=50.0*0.3
		csrx=850.0d0
		phisr=1.0/(1.0+(csrx/csrbz(ix,iy))**10) ! cutoff rate bellow 500 muM

		alphab=ab*dabs(rca)*poz(ix,iy)*phisr ! spark on rate due to LCC 
		bts=1.0/30.0   ! spark off rate 

c              update LCC Markov states
               call markov(dt,vz(ix,iy),cbz(ix,iy),c1z(ix,iy),
     +          c2z(ix,iy),xi1z(ix,iy),xi2z(ix,iy),poz(ix,iy),
     +          c1sz(ix,iy),c2sz(ix,iy),
     +          xi1sz(ix,iy), xi2sz(ix,iy),posz(ix,iy),alphab,bts)
     
     		xrr=5.0*(c1z(ix,iy)+c1sz(ix,iy))

c              boundary RyR current
               gsrb=0.01d0/2.5  
               xryrb=gsrb*csrbz(ix,iy)*pbz(ix,iy)

c              ================================================
c              INTERIOR CALCIUM SPARKS
c              ================================================
               
               	csrxx=900.0d0
	      	chx=850.0d0
		phi=(csriz(ix,iy)/chx)/(1.0+(csrxx/csriz(ix,iy))**10) ! cutoff rate bellow 1000 muM
	
	! activation of sparks

		gi=0.25
		pra=0.05  ! increase threshold here
		xra=(pra/piz(ix,iy))**6
		ra_inf=phi/(1.0+xra)

		dra=(ra_inf-raz(ix,iy))/100.0
	
		rr=0.02*3.0*1.5
	
		pbinf=1.0/(1.0+(pbxz/pbz(ix,iy))**10)

		ar1=(aleak/500.0)*(csriz(ix,iy)/1000.0)
		ar2=rr*pbinf*phi
		ar3=gi*raz(ix,iy)

		alphai=ar1+ar2+ar3

		gryri=0.015
		xryri=gryri*piz(ix,iy)*csriz(ix,iy)

		btsi=1.0/70.0
          
c              ================================================
c              CALCIUM CONCENTRATION UPDATES
c              ================================================
         	
 
	xsarc=-xicaq+xinacaq
	
	tausri=50.0
	dff=(cnsrz(ix,iy)-csriz(ix,iy))/tausri

	tau1=10.0
	tau2=10.0

	dfbi=(cbz(ix,iy)-ciz(ix,iy))/tau1
	dfbisr=(csrbz(ix,iy)-cnsrz(ix,iy))/tau2

	dcbt=xryrb-xupb+xsarc-dfbi
	dcsrb=vbsrz*(-xryrb+xupb)-dfbisr   

	dcit=xryri-xupi+vbiz*dfbi               	
	dcsri=visrz*(-xryri)+dff
	
	dnsr=vnsrz*(xupi)-dff+vbisrz*dfbisr

	cbtz(ix,iy)=cbtz(ix,iy)+dcbt*dt
	citz(ix,iy)=citz(ix,iy)+dcit*dt
	
	csrbz(ix,iy)=csrbz(ix,iy)+dcsrb*dt	
	csriz(ix,iy)=csriz(ix,iy)+dcsri*dt

	cnsrz(ix,iy)=cnsrz(ix,iy)+dnsr*dt
	
	raz(ix,iy)=raz(ix,iy)+dra*dt

            
c              ================================================
c              STOCHASTIC SPARK EVOLUTION
c              ================================================
               nsbx=nsbz(ix,iy)
               call binevol(nbtz,nsbx,alphab,bts,dt,ndeltapx,
     &              ndeltamx)

               if(ndeltamx.gt.nsbx.or.ndeltapx.gt.nbtz) then
                 nsbz(ix,iy)=0
               else
                 nsbz(ix,iy)=nsbz(ix,iy)+ndeltapx-ndeltamx
               endif

               nsix=nsiz(ix,iy)
               call binevol(nitz,nsix,alphai,btsi,dt,ndeltapy,
     &              ndeltamy)
               
               if(ndeltamy.gt.nsix.or.ndeltapy.gt.nitz) then
                 nsiz(ix,iy)=0

               else
                 nsiz(ix,iy)=nsiz(ix,iy)+ndeltapy-ndeltamy
               endif
               

c              ================================================
c              VOLTAGE DYNAMICS WITH ADAPTIVE TIME STEP
c              ================================================
               wca=12.0d0
               xinaca=wca*xinacaq
               xica=2.0d0*wca*xicaq

               adq=dabs(dvz(ix,iy))
               if(adq.gt.25.0d0) then
                 mstp=10
               else
                 mstp=1
               endif
               hode=dt/dfloat(mstp)

               do iii=1 , mstp
               
              
                 call ina(hode,vz(ix,iy),frt,xhz(ix,iy),xjz(ix,iy),
     &            xmz(ix,iy),xnai,xnao,xina,pon,pon_inf)

                 call ikr(hode,vz(ix,iy),frt,xko,xki,
     &            xrz(ix,iy),xikr)

                 call iks(hode,vz(ix,iy),frt,xnao,xnai,
     &            xko,xki,xs1z(ix,iy),qksz(ix,iy),xiks)

                 call ik1z(hode,vz(ix,iy),frt,xki,xko,xik1)

                 call ikur(hode,vz(ix,iy),frt,xki,xko,xkurz(ix,iy),
     &            ykurz(ix,iy),xikur)

                 call ito(hode,vz(ix,iy),frt,xki,xko,xtofz(ix,iy),
     &            ytofz(ix,iy),xito)

                 call inak(vz(ix,iy),frt,xko,xnao,xnai,xinak)

c  	stimulation current
        
        	stimx=0.0

	xkzz=xik1+xikr+xiks+xito+xikur+xinak
        
          wnaz(ix,iy)=-xina*1000.0          
          wicaz(ix,iy)=-xica
          wikaz(ix,iy)=xkzz
          
           poaz(ix,iy)=pon
	 poaz_inf(ix,iy)=pon_inf
                              
            dvh=-(xina+xkzz+xinaca+xica)+stimx
          
         
                 vz(ix,iy)=vz(ix,iy)+dvh*hode

               enddo

             end do
         end do       
               
c        ============================================================
c        BOUNDARY CONDITIONS & COUPLING
c        ============================================================
        call Euler_Forward(Lx,Ly,Lx1,Ly1,
     &  v,vnew,vz,vnewz,dt)

        call Euler_Forward(Lx,Ly,Lx1,Ly1,
     &  v,vnew,vz,vnewz,dt)

c        ============================================================
c        UPDATE VOLTAGE DERIVATIVES
c        ============================================================

        do iy=1,Ly
          do ix=1,Lx
           dv(ix,iy)=(vnew(ix,iy)-vold(ix,iy))/dt
           vold(ix,iy)=vnew(ix,iy)
          enddo
        enddo

        do iy=1,Ly1
          do ix=1,Lx1
           dvz(ix,iy)=(vnewz(ix,iy)-voldz(ix,iy))/dt
           voldz(ix,iy)=vnewz(ix,iy)
          enddo
        enddo

c        ============================================================
c        DATA OUTPUT
c        ============================================================
        if(mod(ncount,5).eq.0.and.t.gt.0.) then

          write(20,*)  t, v(10,2)
          write(21,*)  t, v(100,2)
          write(22,*)  t, vz(80,2)
          
          write(30,*)  t, wica(20,2)
        
        endif

c        ============================================================
c        VISUALIZATION OUTPUT
c        ============================================================

	ttz=0.0*bcl
	
c        if(mod(time,15.).eq.0.and.time.gt.ttz) then

	ty1=640.0
	ty2=ty1+640.0

	if(mod(time,2.).eq.0.and.time.gt.ty1.and.time.lt.1400.) then

        open(unit=1,file=filenm1(kku),status='unknown')
        open(unit=2,file=filenm2(kku),status='unknown')
	open(unit=3,file=filenm3(kku),status='unknown')
	open(unit=4,file=filenm4(kku),status='unknown')
	open(unit=5,file=filenm5(kku),status='unknown')
	open(unit=6,file=filenm6(kku),status='unknown')
	open(unit=7,file=filenm7(kku),status='unknown')

	do ixx=1,Lx
	zmask1(ixx)=v(ixx,3)
	zmask2(ixx)=wika(ixx,3)
	zmask3(ixx)=wica(ixx,3)
	zmask4(ixx)=wna(ixx,3)
	zmask5(ixx)=poa(ixx,3)
	zmask6(ixx)=poa_inf(ixx,3)
	zmask7(ixx)=poa_inf(ixx,3)*poa(ixx,3)
	
	enddo
	
	iz3=Ly1/2

	
	do ixx=1,Lx1

	zmask1(Lx+ixx)=vz(ixx,iz3)
	zmask2(Lx+ixx)=wikaz(ixx,iz3)
	zmask3(Lx+ixx)=wicaz(ixx,iz3)
	zmask4(Lx+ixx)=wnaz(ixx,iz3)
	zmask5(Lx+ixx)=poaz(ixx,iz3)
	zmask6(Lx+ixx)=poaz_inf(ixx,iz3)
	zmask7(Lx+ixx)=poaz_inf(ixx,iz3)*poaz(ixx,iz3)
	
	
	enddo 
	
	do ixx=1,Lx+Lx1

	write(1,*) ixx, zmask1(ixx)
	write(2,*) ixx, zmask2(ixx)
	write(3,*) ixx, zmask3(ixx)
	write(4,*) ixx, zmask4(ixx)
	write(5,*) ixx, zmask5(ixx)
	write(6,*) ixx, zmask6(ixx)
	write(7,*) ixx, zmask7(ixx)
	
	enddo 

          kku=kku+1

            endif


         t=t+dt

 1     continue

 500  stop

      end 
        
        
        
      
c *****************************************************************************************************

	subroutine generate_filenamesv(filenm1)
        implicit none

       ! Output array of filenames
       character*15 filenm1(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm1(ip)='vb'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm1(ip)='vb'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm1(ip)='vb'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine      
c **************************************************************************
      
      subroutine generate_filenamesk(filenm2)
        implicit none

       ! Output array of filenames
       character*15 filenm2(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm2(ip)='ik'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm2(ip)='ik'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm2(ip)='ik'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine      

c **************************************************************************      
 
      subroutine generate_filenamesca(filenm3)
        implicit none

       ! Output array of filenames
       character*15 filenm3(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm3(ip)='icaz'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm3(ip)='icaz'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm3(ip)='icaz'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine            
      
        subroutine generate_filenamesna(filenm4)
        implicit none

       ! Output array of filenames
       character*15 filenm4(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm4(ip)='inaz'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm4(ip)='inaz'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm4(ip)='inaz'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine            
      
       subroutine generate_filenamespo(filenm5)
        implicit none

       ! Output array of filenames
       character*15 filenm5(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm5(ip)='po'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm5(ip)='po'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm5(ip)='po'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine      
      	      
        subroutine generate_filenamespof(filenm6)
        implicit none

       ! Output array of filenames
       character*15 filenm6(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm6(ip)='poinf'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm6(ip)='poinf'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm6(ip)='poinf'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine            
      
       subroutine generate_filenamespofx(filenm7)
        implicit none

       ! Output array of filenames
       character*15 filenm7(500)
      
       ! Local variables
         character*1 cam(0:9)
         integer ip, ifil1, ifil2, ifil3
      
      ! Initialize digit characters
      	cam(0) = '0'
      	cam(1) = '1'
      	cam(2) = '2'
      	cam(3) = '3'
      	cam(4) = '4'
      	cam(5) = '5'
      	cam(6) = '6'
      	cam(7) = '7'
      	cam(8) = '8'
      	cam(9) = '9'

      ! Generate filenames for indices 1-999
      	do ip = 1,500
          ! Extract individual digits using modulo
          ifil1 = mod(ip,10)        ! ones digit
          ifil2 = mod(ip/10,10)     ! tens digit
          ifil3 = mod(ip/100,10)    ! hundreds digit

          ! Format filename based on number of digits
          if(ip.lt.10) then
              ! Single digit (1-9)
              filenm7(ip)='pox'//cam(ifil1)//'.dat'
          else if(ip.lt.100) then
              ! Two digits (10-99)
              filenm7(ip)='pox'//cam(ifil2)//cam(ifil1)//'.dat'
          else
              ! Three digits (100-999)
              filenm7(ip)='pox'//cam(ifil3)//cam(ifil2)//
     &                    cam(ifil1)//'.dat'
          endif
      	enddo

      	return
      	end subroutine            
      
      
