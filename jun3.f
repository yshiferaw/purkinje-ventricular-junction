

      subroutine initCa1(Lx,Ly,cb,ci,csrb,csri,cnsr,po,c1,c2,xi1,
     #    xi2,ra,pos,c1s,c2s,xi1s,xi2s,cit,cbt,nsb,nsi)

      implicit double precision (a-h,o-z)
      
c     Calcium concentration arrays
      double precision cb(Lx,Ly),ci(Lx,Ly),csrb(Lx,Ly),csri(Lx,Ly)
      double precision cnsr(Lx,Ly), cit(Lx,Ly), cbt(Lx,Ly)
      
c     L-type calcium channel states
      double precision po(Lx,Ly), c1(Lx,Ly), c2(Lx,Ly)
      double precision xi1(Lx,Ly), xi2(Lx,Ly), xi2s(Lx,Ly)
      double precision c1s(Lx,Ly), c2s(Lx,Ly), xi1s(Lx,Ly)
      double precision pos(Lx,Ly)
      
c     Spark-related variables
      double precision ra(Lx,Ly)
      integer nsb(Lx,Ly), nsi(Lx,Ly)

c     Initial concentration values
      cxinit = 800.0d0    ! Initial SR calcium (μM)
      cixx = 0.3d0       ! Initial cytosolic calcium (μM)

      do ix = 1, Lx
        do iy = 1, Ly
            
c         ================================================
c         CALCIUM CONCENTRATIONS
c         ================================================
          cb(ix,iy) = cixx      ! boundary cytosolic Ca
          ci(ix,iy) = cixx      ! interior cytosolic Ca
          csrb(ix,iy) = cxinit  ! boundary SR Ca
          csri(ix,iy) = cxinit  ! interior SR Ca
          cnsr(ix,iy) = cxinit  ! NSR Ca
 
c         ================================================
c         L-TYPE Ca CHANNEL STATES (facing cytosol)
c         ================================================
          po(ix,iy) = 0.0d0     ! open probability
          c1(ix,iy) = 0.0d0     ! closed state 1
          c2(ix,iy) = 1.0d0     ! closed state 2 (initial)
          xi1(ix,iy) = 0.0d0    ! inactivated state 1
          xi2(ix,iy) = 0.0d0    ! inactivated state 2

c         ================================================
c         L-TYPE Ca CHANNEL STATES (facing spark)
c         ================================================
          pos(ix,iy) = 0.0d0    ! open probability
          c1s(ix,iy) = 0.0d0    ! closed state 1
          c2s(ix,iy) = 0.0d0    ! closed state 2
          xi1s(ix,iy) = 0.0d0   ! inactivated state 1
          xi2s(ix,iy) = 0.0d0   ! inactivated state 2

c         ================================================
c         SPARK VARIABLES
c         ================================================
          ra(ix,iy) = 0.0d0     ! spark activation variable
          nsb(ix,iy) = 5        ! number of boundary sparks
          nsi(ix,iy) = 5        ! number of interior sparks

c         ================================================
c         CONVERT FREE TO TOTAL CALCIUM
c         ================================================
          call total(ci(ix,iy), cit(ix,iy))
          call total(cb(ix,iy), cbt(ix,iy))

        enddo
      enddo

      return
      end

	 subroutine initCa2(Lx,Ly,cb,ci,csrb,csri,cit,cbt,d_gate,f_gate,nsb)

      implicit double precision (a-h,o-z)
      
c     Calcium concentration arrays
      double precision cb(Lx,Ly),ci(Lx,Ly),csrb(Lx,Ly),csri(Lx,Ly)
      double precision cit(Lx,Ly), cbt(Lx,Ly)
      double precision d_gate(Lx,Ly), f_gate(Lx,Ly)
      

c     Spark-related variables
      double precision ra(Lx,Ly)
      integer nsb(Lx,Ly)

c     Initial concentration values
      cxinit = 800.0d0    ! Initial SR calcium (μM)
      cixx = 0.3d0       ! Initial cytosolic calcium (μM)

      do ix = 1, Lx
        do iy = 1, Ly
            
c         ================================================
c         CALCIUM CONCENTRATIONS
c         ================================================
          cb(ix,iy) = cixx      ! boundary cytosolic Ca
          ci(ix,iy) = cixx      ! interior cytosolic Ca
          csrb(ix,iy) = cxinit  ! boundary SR Ca
          csri(ix,iy) = cxinit  ! interior SR Ca
   
c         ================================================
c         L-TYPE Ca CHANNEL STATES (facing cytosol)
c         ================================================
       
         d_gate(ix,iy)=0.01
         f_gate(ix,iy)=0.9


c         ================================================
c         SPARK VARIABLES
c         ================================================
          
          nsb(ix,iy) = 5        ! number of boundary sparks
      
c         ================================================
c         CONVERT FREE TO TOTAL CALCIUM
c         ================================================
          call total(ci(ix,iy), cit(ix,iy))
          call total(cb(ix,iy), cbt(ix,iy))

        enddo
      enddo

      return
      end



c     ================================================================
c     SUBROUTINE: INITIALIZE VOLTAGE AND GATING VARIABLES
c     ================================================================

      subroutine initV(Lx,Ly,v,xm,xh,xj,xr,xs1,qks,xkur,
     #    ykur,xtof,ytof,vold)

      implicit double precision (a-h,o-z)

c     Voltage arrays
      double precision v(0:Lx+1,0:Ly+1), vold(0:Lx,0:Ly)
      
c     Sodium channel gating variables
      double precision xm(Lx,Ly), xh(Lx,Ly), xj(Lx,Ly)
      
c     Potassium channel gating variables
      double precision xr(Lx,Ly), xs1(Lx,Ly), qks(Lx,Ly)
      double precision xkur(Lx,Ly), ykur(Lx,Ly)
      double precision xtof(Lx,Ly), ytof(Lx,Ly)

      do ix = 1, Lx
        do iy = 1, Ly
            
c         ================================================
c         MEMBRANE VOLTAGE
c         ================================================
          v(ix,iy) = -90.0d0        ! resting potential (mV)
          vold(ix,iy) = v(ix,iy)    ! copy for time stepping

c         ================================================
c         SODIUM CHANNEL GATING
c         ================================================
          xm(ix,iy) = 0.001d0       ! activation gate (closed)
          xh(ix,iy) = 1.0d0         ! fast inactivation (open)
          xj(ix,iy) = 1.0d0         ! slow inactivation (open)

c         ================================================
c         POTASSIUM CHANNEL GATING
c         ================================================
c         IKr (rapid delayed rectifier)
          xr(ix,iy) = 0.0d0         ! activation gate (closed)

c         IKs (slow delayed rectifier)
          xs1(ix,iy) = 0.3d0        ! activation gate
          qks(ix,iy) = 0.2d0        ! calcium modulation

c         IKur (ultra-rapid delayed rectifier)
          xkur(ix,iy) = 0.01d0      ! activation gate
          ykur(ix,iy) = 1.0d0       ! inactivation gate

c         Ito (transient outward)
          xtof(ix,iy) = 0.02d0      ! fast activation
          ytof(ix,iy) = 0.8d0       ! slow inactivation

        enddo
      enddo

      return
      end


c     ================================================================
c     SUBROUTINE: EULER FORWARD DIFFUSION (TWO TISSUES ONLY)
c     ================================================================
      subroutine Euler_forward(LLx,LLy,LLx1,LLy1,
     #    v,vnew,vz,vnewz,dt)
      implicit double precision (a-h,o-z)
        
c     Voltage arrays for two tissues
      double precision v(0:LLx+1,0:LLy+1), vnew(0:LLx+1,0:LLy+1)
      double precision vz(0:LLx1+1,0:LLy1+1), vnewz(0:LLx1+1,0:LLy1+1)

c     ================================================
c     DIFFUSION PARAMETERS
c     ================================================
      dx = 0.015d0      ! spatial step in x direction (cm)
      dy = 0.015d0      ! spatial step in y direction (cm)
      Dfu1 =0.0005
      
c      0.001  ! tissue 1 diffusion
      Dfu2 =0.0005  ! tissue 2 diffusion      
c      0.0005d0/2.0d0    ! diffusion coefficient

      sx1 = Dfu1*dt/4.0d0/dx/dx
      sy1 = Dfu1*dt/4.0d0/dy/dy

      sx2 = Dfu2*dt/4.0d0/dx/dx
      sy2 = Dfu2*dt/4.0d0/dy/dy
      
      sxd = Dfu1*dt/4.0d0/dx/dx  ! couples the two pieces of tissue
c	sxd=0.0d0

c     ================================================
c     FIRST HALF TIME STEP (dt/4) - TISSUE 1
c     ================================================

c     Set corner boundary conditions
      v(0,0) = v(2,2)
      v(0,LLy+1) = v(2,LLy-1)
      v(LLx+1,0) = v(LLx-1,2)
      v(LLx+1,LLy+1) = v(LLx-1,LLy-1)

c     Set edge boundary conditions (no flux)
      do ix = 1, LLx
        v(ix,0) = v(ix,2)
        v(ix,LLy+1) = v(ix,LLy-1)
      enddo
  
      do iy = 1, LLy
        v(0,iy) = v(2,iy)
        v(LLx+1,iy) = v(LLx-1,iy)
      enddo

c     Interior points diffusion

	nyz=LLy1/2-LLy/2

      do ix = 1, LLx
      
      sxp=0.0d0
      if(ix.eq.LLx) then
      sxp=sxd
      endif 
          
        do iy = 1, LLy
        
        iz1=nyz+iy
        
          vnew(ix,iy) = v(ix,iy) + sx1 * (v(ix+1,iy) + v(ix-1,iy)
     #        - 2.0d0 * v(ix,iy)) + sy1 * (v(ix,iy+1) + v(ix,iy-1)
     #        - 2.0d0 * v(ix,iy))+sxp*(vz(1,iz1)-v(LLx,iy))
     
     
        enddo
      enddo


c     ================================================
c     FIRST HALF TIME STEP (dt/4) - TISSUE 2
c     ================================================
c     Set corner boundary conditions

      vz(0,0) = vz(2,2)
      vz(0,LLy1+1) = vz(2,LLy1-1)
      vz(LLx1+1,0) = vz(LLx1-1,2)
      vz(LLx1+1,LLy1+1) = vz(LLx1-1,LLy1-1)

c     Set edge boundary conditions (no flux)
      do ix = 1, LLx1
        vz(ix,0) = vz(ix,2)
        vz(ix,LLy1+1) = vz(ix,LLy1-1)
      enddo
  
      do iy = 1, LLy1
        vz(0,iy) = vz(2,iy)
        vz(LLx1+1,iy) = vz(LLx1-1,iy)
      enddo

	iu1=nyz+1
	iu2=nyz+LLy
	
c     Interior points diffusion
      do ix = 1, LLx1 
        do iy = 1, LLy1
        
      sxp=0.0d0
      if(ix.eq.1.and.iy.le.iu2.and.iy.ge.iu1) then
      sxp=sxd
      endif 
       
        iu3=iy-nyz
        
          vnewz(ix,iy) = vz(ix,iy) + sx2*(vz(ix+1,iy) + vz(ix-1,iy)
     #        - 2.0d0 * vz(ix,iy)) + sy2*(vz(ix,iy+1) + vz(ix,iy-1)
     #        - 2.0d0 * vz(ix,iy)) + sxp*(v(LLx,iu3)-vz(1,iy))
     
 
     
        enddo
      enddo



c     ================================================
c     SECOND HALF TIME STEP (dt/4) - TISSUE 1
c     ================================================
c     Set corner boundary conditions
      vnew(0,0) = vnew(2,2)
      vnew(0,LLy+1) = vnew(2,LLy-1)
      vnew(LLx+1,0) = vnew(LLx-1,2)
      vnew(LLx+1,LLy+1) = vnew(LLx-1,LLy-1)

c     Set edge boundary conditions
      do ix = 1, LLx
        vnew(ix,0) = vnew(ix,2)
        vnew(ix,LLy+1) = vnew(ix,LLy-1)
      enddo

      do iy = 1, LLy
        vnew(0,iy) = vnew(2,iy)
        vnew(LLx+1,iy) = vnew(LLx-1,iy)
      enddo

	

      do ix = 1, LLx
      
      sxp=0.0d0
      if(ix.eq.LLx) then
      sxp=sxd
      endif 
          
        do iy = 1, LLy
        
      	 iz1=nyz+iy
        
        ru1=sx1*(vnew(ix+1,iy) + vnew(ix-1,iy)- 2.0d0 * vnew(ix,iy))
        ru2=sy1*(vnew(ix,iy+1) + vnew(ix,iy-1)- 2.0d0 * vnew(ix,iy))
        ru3=sxp*(vnewz(1,iz1)-vnew(LLx,iy))
        
        v(ix,iy) = vnew(ix,iy)+ru1+ru2+ru3
           
        enddo
      enddo



c     ================================================
c     SECOND HALF TIME STEP (dt/4) - TISSUE 2
c     ================================================
c     Set corner boundary conditions
      vnewz(0,0) = vnewz(2,2)
      vnewz(0,LLy1+1) = vnewz(2,LLy1-1)
      vnewz(LLx1+1,0) = vnewz(LLx1-1,2)
      vnewz(LLx1+1,LLy1+1) = vnewz(LLx1-1,LLy1-1)

c     Set edge boundary conditions
      do ix = 1, LLx1
        vnewz(ix,0) = vnewz(ix,2)
        vnewz(ix,LLy1+1) = vnewz(ix,LLy1-1)
      enddo

      do iy = 1, LLy1
        vnewz(0,iy) = vnewz(2,iy)
        vnewz(LLx1+1,iy) = vnewz(LLx1-1,iy)
      enddo

	iu1=nyz+1
	iu2=nyz+LLy

c     Interior points diffusion
      do ix = 1, LLx1
        do iy = 1, LLy1
        
      sxp=0.0d0
      if(ix.eq.1.and.iy.le.iu2.and.iy.ge.iu1) then
      sxp=sxd
      endif 
        
        iu3=iy-nyz
        
	ru1=sx2*(vnewz(ix+1,iy)+vnewz(ix-1,iy)-2.0d0*vnewz(ix,iy))
        ru2=sy2*(vnewz(ix,iy+1)+vnewz(ix,iy-1)-2.0d0*vnewz(ix,iy))
	ru3=sxp*(vnew(LLx,iu3)-vnewz(ix,iy))

     	vz(ix,iy)=vnewz(ix,iy)+ru1+ru2+ru3
     
        enddo
      enddo


      return
      end
