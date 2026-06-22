
      subroutine ica_luorudy_simple(dt,v,d_gate,f_gate)
      implicit double precision (a-h,o-z)

      
c     ================================================================
c     ACTIVATION GATE (d)
c     ================================================================

c      vrx=10.0

      vrx=10.0
      arx=7.0
      alpha=2.0	
	
      d_inf = 1.0d0 / (1.0d0 + dexp(-(v + vrx)/arx))
      
      
      tau_d = 1.0d0 / (1.0d0 + dexp(-(v + 10.0d0)/6.24d0))
      tau_d = tau_d * (1.0d0 - dexp(-(v + 10.0d0)/6.24d0)) / 
     &        (0.035d0 * (v + 10.0d0))
      
      
      
      if (dabs(v + 10.0d0) .lt. 1.0d-6) then
         tau_d = 2.5d0  ! avoid division by zero
      endif
    
      if (tau_d .lt. 0.5d0) tau_d = 0.5d0


	tau_d=tau_d*alpha
	
c     ================================================================
c     INACTIVATION GATE (f)
c     ================================================================

c      vrf=28.0d0
c      vrf=27.0
c      arf=5.
 
 	vrf=27.0
 	arf=6.1
 	beta=0.23
 	
      
      f_inf = 1.0d0 / (1.0d0 + dexp((v + vrf)/arf))
      
      
      tau_f = 9.0d0 / (0.0197d0 * dexp(-0.0337d0**2 * 
     &        (v + 10.0d0)**2) + 0.02d0)

	tau_f=tau_f*beta

c     ================================================================
c     UPDATE GATES
c     ================================================================
      d_gate = d_inf - (d_inf - d_gate) * dexp(-dt/tau_d)
      f_gate = f_inf - (f_inf - f_gate) * dexp(-dt/tau_f)

      
      return
      end	


	
c **************************************************************************************
	
        subroutine ina(hode,v,frt,xh,xj,xm,xnai,xnao,xina,po,po_inf)
       	implicit double precision (a-h,o-z)

       	gna=12.0d0
       	gna=16.0d0
       	
       	XKMCAM=0.3d0
       	deltax=-0.18d0
       
       	ena = (1.0d0/frt)*dlog(xnao/xnai) ! sodium reversal potential

       	am = 0.32d0*(v+47.13d0)/(1.0d0-dexp(-0.1d0*(v+47.13d0)))
       	bm = 0.08d0*dexp(-v/11.0d0)

c	camfact= 1.0d0/(1.0d0+(XKMCAM/caM)**4)
c	vshift = -3.25*camfact

	camfact=0.	
	vshift=0.0d0
	
	vx=v-vshift		
	 	
        if(vx.lt.(-40.0d0)) then
         ah=0.135*dexp((80.0+vx)/(-6.8d0))
         bh=3.56*dexp(0.079*vx)+310000.0d0*dexp(0.35d0*vx)

	 aj1a=-127140.0*dexp(0.2444*vx)
	 aj1b=0.00003474d0*dexp(-0.04391d0*vx)
	 aj1c=(vx+37.78)/(1.0d0+dexp(0.311*(vx+79.23)))

	aj=(1.0d0+camfact*deltax)*(aj1a-aj1b)*aj1c
  	bj=(0.1212*dexp(-0.01052*vx))/(1.0+dexp(-0.1378d0*(vx+40.14d0)))

        else

         ah=0.0d0
         bh=1.0d0/(0.13d0*(1.0d0+dexp((vx+10.66)/(-11.1d0))))
         aj=0.0d0

	 bja=0.3*dexp(-0.0000002535d0*vx)
	 bjb=1.0+dexp(-0.1d0*(vx+32.0d0))

	bj=bja/bjb
        endif          
	
	tauhx=1.0d0/(ah+bh)  
	
	tauh=tauhx
	
      	taujx=1.0d0/(aj+bj)
      	
      	tauj=taujx
      	
      	taum=1.0d0/(am+bm)
 
	xina=gna*xh*xj*xm*xm*xm*(v-ena)


c	xh_inf=ah/(ah+bh)
c	xh_inf=am/(am+bm)
c	xh_inf=aj/(aj+bj)
	

	xh_inf=ah/(ah+bh)
	xj_inf=aj/(aj+bj)
	xm_inf=am/(am+bm)


c        po=xm**3*xh
c	po_inf=(xm_inf**3)*xh_inf

	po=xm**3
	po_inf=xh
c  xm_inf


        xh = ah/(ah+bh)-((ah/(ah+bh))-xh)*dexp(-hode/tauh)
        xj = aj/(aj+bj)-((aj/(aj+bj))-xj)*dexp(-hode/tauj)
        xm = am/(am+bm)-((am/(am+bm))-xm)*dexp(-hode/taum)

	return
	end 


c *****************************************************************************

	subroutine ikr(hode,v,frt,xko,xki,xr,xikr)
	implicit double precision (a-h,o-z)


      	ek = (1.0d0/frt)*dlog(xko/xki) ! K reversal potential
    
      	gss=dsqrt(xko/5.40)
      	xkrv1=0.00138d0*(v+7.0d0)/( 1.0-dexp(-0.123*(v+7.0d0))  )
      	xkrv2=0.00061d0*(v+10.0d0)/(dexp( 0.145d0*(v+10.0d0))-1.0d0)
      	taukr=1.0d0/(xkrv1+xkrv2)

      	xkrinf=1.0d0/(1.0d0+dexp(-(v+50.0d0)/7.5d0))

      	rg=1.0d0/(1.0d0+dexp((v+33.0d0)/22.4d0))
              
      	gkr=0.007836d0*1.0  ! Ikr conductance 
c	gkr=gkr*1.6

      	xikr=gkr*gss*xr*rg*(v-ek)

      	xr=xkrinf-(xkrinf-xr)*dexp(-hode/taukr)       
    

	return
	end 

c *****************************************************************************

	subroutine iks(hode,v,frt,xnao,xnai,xko,xki,xs1,qks,xiks)
	implicit double precision (a-h,o-z)

       prnak=0.018330d0
      
       qks_inf=0.6d0*0.3*1.6
       tauqks=1000.0d0
      
      eks=(1.0d0/frt)*dlog((xko+prnak*xnao)/(xki+prnak*xnai))
      xs1ss=1.0/(1.0+dexp(-(v-1.50d0)/16.70d0))
      xs2ss=xs1ss

      tauxs=1.0d0/(0.0000719*(v+30.0d0)/(1.0d0-dexp(
     +-0.148d0*(v+30.0)))+0.000131d0
     + *(v+30.0d0)/(dexp(0.0687d0*(v+30.0d0))-1.0d0))
 
c	gksx=1.50d0     
        gksx=0.200d0 ! Iks conductance 
       xiks=gksx*qks*xs1**2*(v-eks) 

       xs1=xs1ss-(xs1ss-xs1)*dexp(-hode/tauxs)
       qks=qks+hode*(qks_inf-qks)/tauqks

	return
	end 

c *****************************************************************************

	subroutine ik1(hode,v,frt,xki,xko,xik1)
	implicit double precision (a-h,o-z)

	ek = (1.0d0/frt)*dlog(xko/xki)     

      gkix=0.50d0*0.4
      gki=gkix*(dsqrt(xko/5.4))
      aki=1.02/(1.0+dexp(0.2385*(v-ek-59.215)))
      bki=(0.49124*dexp(0.08032*(v-ek+5.476))+dexp(0.061750
     +*(v-ek-594.31)))/(1.0+dexp(-0.5143*(v-ek+4.753)))
      xkin=aki/(aki+bki)
      xik1=gki*xkin*(v-ek)

	return
	end 

c *****************************************************************************

	subroutine ik1z(hode,v,frt,xki,xko,xik1)
	implicit double precision (a-h,o-z)

	ek = (1.0d0/frt)*dlog(xko/xki)     

c      gkix=0.50d0*0.2 
       gkix=0.5d0*0.4

      gki=gkix*(dsqrt(xko/5.4))
      aki=1.02/(1.0+dexp(0.2385*(v-ek-59.215)))
      bki=(0.49124*dexp(0.08032*(v-ek+5.476))+dexp(0.061750
     +*(v-ek-594.31)))/(1.0+dexp(-0.5143*(v-ek+4.753)))
      xkin=aki/(aki+bki)
      xik1=gki*xkin*(v-ek)

	return
	end 


c *****************************************************************************

	subroutine ikur(hode,v,frt,xki,xko,xkur,ykur,xikur)
	implicit double precision (a-h,o-z)
	! Grandi removes Ito slow in the Atria
	! Here, we replace it with Ikur


	Gkur=0.05*0.7
c	Gkur=Gkur*1.2  ! Final value: 0.042
	Gkur=0.05*0.7*0.9

	ek = (1.0d0/frt)*dlog(xko/xki)

	rh1=(v+6.)/(-8.6)
	xkurss=1./(1.+dexp(rh1))
        tauxkur=9.0d0/(1.+dexp((v+5.0d0)/12.0))+0.5

	ykurss  =1.0d0/(1.0d0+dexp((v+7.5)/10.0d0))
	tauykur = 590.0d0/(1.+dexp((v+60.)/10.0))+3050.0d0

	xkur = xkurss-(xkurss-xkur)*dexp(-hode/tauxkur)
	ykur = ykurss-(ykurss-ykur)*dexp(-hode/tauykur)

	xikur = Gkur*xkur*ykur*(v-ek)

	return
	end
	
c *****************************************************************************

       	subroutine ito(hode,v,frt,xki,xko,xtof,ytof,xito)
	implicit double precision (a-h,o-z)
	

	ek = (1.0d0/frt)*dlog(xko/xki)

	gtof=0.03d0*1.5*0.8
	gtof=gtof  ! Final value: 0.0324

	rg1=-(v+1.0)/11.0
	xtof_inf=1./(1.+dexp(rg1))
	
	rg3=-(v/30.)**2
	txf=3.5*dexp(rg3)+1.5
	

	rg2=(v+40.5)/11.5
	ytof_inf=1.0/(1.+dexp(rg2))

	rg4=-((v+52.45)/15.8827)**2
	tyf=25.635*dexp(rg4)+24.14

	xitof=gtof*xtof*ytof*(v-ek)
	
        xtof = xtof_inf-(xtof_inf-xtof)*dexp(-hode/txf)
        ytof = ytof_inf-(ytof_inf-ytof)*dexp(-hode/tyf)

	xito=xitof ! total ito current (itos removed in Grandi)


	return
	end 
	
c *****************************************************************************

       	subroutine itoz(hode,v,frt,xki,xko,xtof,ytof,xito)
	implicit double precision (a-h,o-z)
	

	ek = (1.0d0/frt)*dlog(xko/xki)

	gtof=0.03d0*1.5*0.8
	gtof=gtof

	rg1=-(v+1.0)/11.0
	xtof_inf=1./(1.+dexp(rg1))
	
	rg3=-(v/30.)**2
	txf=3.5*dexp(rg3)+1.5

	rg2=(v+40.5)/11.5
	ytof_inf=1.0/(1.+dexp(rg2))

	rg4=-((v+52.45)/15.8827)**2
	tyf=25.635*dexp(rg4)+24.14


	xitof=gtof*xtof*ytof*(v-ek)
	
        xtof = xtof_inf-(xtof_inf-xtof)*dexp(-hode/txf)
        ytof = ytof_inf-(ytof_inf-ytof)*dexp(-hode/tyf)

	xito=xitof ! total ito current (itos removed in Grandi)


	return
	end 	
	

c *****************************************************************************


	subroutine inak(v,frt,xko,xnao,xnai,xinak)
	implicit double precision (a-h,o-z)

c -------  Inak (sodium-potassium exchanger) following Shannon --------------

	 xibarnak=1.50d0

         xkmko=1.5d0 ! these are Inak constants adjusted to fit
c                  ! the experimentally measured dynamic restitution curve
         xkmnai=12.0d0
         hh=1.0d0  ! Na dependence exponent

         sigma = (dexp(xnao/67.3d0)-1.0d0)/7.0d0  
      fnak = 1.0d0/(1+0.1245*dexp(-0.1*v*frt)
     +           +0.0365*sigma*dexp(-v*frt)) 
      xinak = xibarnak*fnak*(1./(1.+(xkmnai/xnai)**hh))*xko/(xko+xkmko)

	return
	end 


	


