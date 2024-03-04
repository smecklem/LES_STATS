#include "us3d_header.h"

c  ***********************************************************************************
c     This user subroutine illustrates using an external function to read additional
c     solution variables and to include them in the output of a full volume solution.
c
c     Written by Heath B. Johnson                           Last modified: 2016-03-16
c  ***********************************************************************************
      module my_functions
      use POB_MODULE
      use POST_MODULE
      use US3D_EXTRAS

#include "interface_exv_fun.f"

      contains

!      ! --------------------------------------------------------------
!      ! User-defined function to read and store stats
      function exv_rates(pob,np,ige,s,nexv,exv_names) result (rates)
      implicit none
      type(pob_t), intent(IN), target :: pob
      integer, intent(IN) :: np,nexv
      integer, dimension(np), intent(IN) :: ige
      real*8, dimension(:,:), target, intent(IN):: s
      character(LEN=*), dimension(nexv), intent(IN) :: exv_names
      real*8, dimension(nexv,np):: stats
      real*8, dimension(nexv,np) :: rates
!
!      integer :: ns,ier
!      real*8, dimension(:,:), allocatable :: sm,sp
!
!      ns= pob%ns
!
!      ! -- Allocate storage for reading mean and stats
!      allocate(sm(10,np),sp(12,np))
!
!      ! -- Read mean and stats
!      call us3d_data_r(pob%rid,'stats_mean',pob%nel,np,sm,ige,0,ier)
!      if (ier/=0) stop
!
!      call us3d_data_r(pob%rid,'stats_stat',pob%nel,np,sp,ige,0,ier)
!      if (ier/=0) stop
!
!      stats(1:3,:)= sm(ns+1:ns+3,:)
!      stats(4:6,:)= sp(ns+1:ns+3,:)
!
!      ! -- Cleanup
!      deallocate(sm)
!      deallocate(sp)
!
      Integer :: ier,n,i1,i2,i3,i4,j1,j2,j3,j4,i,m,ns,icsd(pob%gas%nr,7)
      Real(8),dimension(:),allocatable :: t,p,g,a,hf,cv,tv
      Real(8),dimension(:,:), allocatable :: rhos
      Real(8) :: tt,ttv,ti,tvi,tab,tabi,tabl,taf,tafi,tafl,rvib
      Real(8) :: gg,ccv,pp,rcvi,xf,xb,rki,rkeq,drk,rkf,rkfp,rkb,rkinf
      Real(8) :: rkfs,rkbs,rkps,rkfm,rkbm,rkpm,rate,dkft,dkbt,drt
      Real(8) :: tr,tri,ws(pob%ns),rs(pob%ns),xs(pob%ns)
      Real(8) :: aeta(pob%gas%nr),athet(pob%gas%nr)

      real*8 rkb0,rkf0,rkbinf,tmp1,ifac,mfac,bfac
      real*8 ea_0,ea_inf,cfm_inf
      Real*8 athetu,athetu0,cfmu,cfmu0,aetau,aetau0,tbeu(30)
      real*8 cfm_0(pob%ns)
      real*8 rkfs0,rkbs0,rkfm0,rkbm0,pr,prb
      Real*8 atroe,t1,t2,t3,Fcent,f1,cc,nn,ftp
      Real(8) :: cfm(pob%gas%nr,pob%ns)
      Character(100) :: word
      Type(gas_t) :: gas

      Allocate(t(np),p(np),g(np),a(np),hf(np),cv(np),tv(np))
      Allocate(rhos(pob%ns,np))

      ! Setup some stuff that I couldn't import because models and
      ! simvars and chemistry aren't initialised in postpar
      Gas  = pob%gas
      Ns   = pob%ns
      Aeta = gas%aeta
      Athet= gas%athet
      Cfm  = gas%cfm
      Icsd = gas%icsd
      Rates(:,:) = 0.0d0

      Select Case(pob%ivib)
        Case Default 
          rvib = 0.0d0
        Case(1,2)   
          rvib = 1.0d0
      End Select

      Rhos(:,:) = s(1:pob%ns,:)
      T(:)   = s(pob%vid_t,:)
      Tv(:)  = s(pob%vid_tv,:)
      Call gas_therm2(gas,np,rhos,t,p,g,a,hf,cv)


      Do i = 1,np   ! Loop through cells in this chunk (probably 200,000 of 'em)
        tr  = t(i)
        tri = 1.0d0/tr
        tt  = t(i) !max(t(i), chem_reac_tmin)
        ttv = tv(i)!max(tv(i),chem_reac_tmin)
        ttv = ttv*rvib + tt*(1.0d0 - rvib)
        ti  = 1.0d0/tt
        tvi = 1.0d0/ttv

        gg  = g(i)
        ccv = cv(i)
        pp  = p(i)
        rs(:) = rhos(:,i)
        rcvi = 1.0d0/(sum(rs)*cv(i))
        ws = 0.0d0

        do m = 1,ns
          xs(m) = rs(m)*gas%smwi(m) ! Mole Densities
        enddo

        ! Compute reaction rates and source terms
        Do n = 1,gas%nr  ! Loop over each reaction

          if(gas%irxon(n) .eq. 0) cycle        ! Skip reaction if user requests

          taf  = (tt**gas%txf(n))*(ttv**(1.0d0-gas%txf(n)))
          tafi = 1.0d0/taf
          tafl = log(taf)
          tab  = (tt**gas%txb(n))*(ttv**(1.0d0-gas%txb(n)))
          tabi = 1.0d0/tab
          tabl = log(tab)

        Select Case(gas%irxtyp(n))
          Case Default                         ! Reaction undefined
            word = 'UNKNOWN'
            write(*,*)  gas%irxtyp(n)
            goto 901

          Case(1)                              ! Dissociation reactions
              !write(*,*) i1,j1,j2

            i1 = icsd(n,1)
            j1 = icsd(n,2)
            j2 = icsd(n,3)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq
            !rkb  = exp(aeta(n)*tabl - athet(n)*tabi)*rkeq ! Save an exponential

            xf = xs(i1)
            xb = xs(j1)*xs(j2)

            rkfs = 0.0d0
            rkps = 0.0d0
            rkbs = 0.0d0
            do m = 1,ns
              rkfm = cfm(n,m)*rkf
              rkpm = cfm(n,m)*rkfp
              rkbm = cfm(n,m)*rkb
              rkfs = rkfs + rkfm*xs(m)
              rkps = rkps + rkpm*xs(m)
              rkbs = rkbs + rkbm*xs(m)
            enddo

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate

          Case(2)                              ! Single partner dissociation reactions

            i1 = icsd(n,1)                     ! Reactant
            i2 = icsd(n,2)                     ! Collision partner
            j1 = icsd(n,3)                     ! Product 1
            j2 = icsd(n,4)                     ! Product 2

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)        ! Forward rate coefficient @ Ta
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)        ! Forward rate coefficient @ T
            rkb  = rkfp*rkeq                                    ! Backward rate coefficient

            xf = xs(i1)
            xb = xs(j1)*xs(j2)

            rkfm = cfm(n,1)*rkf
            rkbm = cfm(n,1)*rkb
            rkfs = rkfm         *xs(i2)
            rkps = cfm(n,1)*rkfp*xs(i2)
            rkbs = rkbm         *xs(i2)

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate

          Case(3)                              ! Recombination reactions

            i1 = icsd(n,1)                     ! Reactant 1
            i2 = icsd(n,2)                     ! Reactant 2
            j1 = icsd(n,3)                     ! Product

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)            ! Forward rate coefficient @ Taf
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)            ! Forward rate coefficient @ Tab
            rkb  = rkfp*rkeq                                    ! Backward rate coefficient

            xf = xs(i1)*xs(i2)
            xb = xs(j1)

            rkfs = 0.0d0
            rkps = 0.0d0
            rkbs = 0.0d0
            do m = 1,ns
              rkfm = cfm(n,m)*rkf
              rkpm = cfm(n,m)*rkfp
              rkbm = cfm(n,m)*rkb
              rkfs = rkfs + rkfm*xs(m)
              rkps = rkps + rkpm*xs(m)
              rkbs = rkbs + rkbm*xs(m)

            enddo

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate


          Case(4)                              ! Single partner recombination reactions

            i1 = icsd(n,1)                     ! Reactant 1
            i2 = icsd(n,2)                     ! Reactant 2
            i3 = icsd(n,3)                     ! Collision partner
            j1 = icsd(n,4)                     ! Product

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)            ! Forward rate coefficient @ T
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)            ! Forward rate coefficient @ Ta
            rkb  = rkfp*rkeq                                    ! Backward rate coefficient

            xf = xs(i1)*xs(i2)
            xb = xs(j1)

            rkfm = cfm(n,1)*rkf
            rkbm = cfm(n,1)*rkb
            rkfs = rkfm         *xs(i3)
            rkps = cfm(n,1)*rkfp*xs(i3)
            rkbs = rkbm         *xs(i3)

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate

          Case(5)                              ! Exchange reactions

            i1 = icsd(n,1)
            i2 = icsd(n,2)
            j1 = icsd(n,3)
            j2 = icsd(n,4)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq

            xf = xs(i1)*xs(i2)
            xb = xs(j1)*xs(j2)

            rkfs = cfm(n,1)*rkf
            rkps = cfm(n,1)*rkfp
            rkbs = cfm(n,1)*rkb

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate

          Case(6)                              ! Electron impact ionization reactions

            i1 = icsd(n,1)
            i2 = icsd(n,2)                     ! should be electrons
            j1 = icsd(n,3)
            j2 = icsd(n,4)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq

            xf = xs(i1)
            xb = xs(j1)*xs(j2)

            rkfs = cfm(n,1)*rkf
            rkps = cfm(n,1)*rkfp
            rkbs = cfm(n,1)*rkb

            rate = (rkbs*xb - rkfs*xf)*xs(i2)

            ws(i1) = ws(i1) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate

          Case(12)                              ! 3 product, 2 reactant reactions

            i1 = icsd(n,1)
            i2 = icsd(n,2)
            j1 = icsd(n,3)
            j2 = icsd(n,4)
            j3 = icsd(n,5)
            j4 = icsd(n,6)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq

            xf = xs(i1)*xs(i2)
            xb = xs(j1)*xs(j2)*xs(j3)

            rkfs = cfm(n,1)*rkf
            rkps = cfm(n,1)*rkfp
            rkbs = cfm(n,1)*rkb

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate
            ws(j3) = ws(j3) - rate
            ws(j4) = ws(j4) - rate

          Case(10)                              ! 3 product, 2 reactant reactions

            i1 = icsd(n,1)
            i2 = icsd(n,2)
            j1 = icsd(n,3)
            j2 = icsd(n,4)
            j3 = icsd(n,5)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq

            xf = xs(i1)*xs(i2)
            xb = xs(j1)*xs(j2)*xs(j3)

            rkfs = cfm(n,1)*rkf
            rkps = cfm(n,1)*rkfp
            rkbs = cfm(n,1)*rkb

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate
            ws(j3) = ws(j3) - rate

          Case(11)                              ! 2 product, 3 reactant reactions

            i1 = icsd(n,1)
            i2 = icsd(n,2)
            i3 = icsd(n,3)
            j1 = icsd(n,4)
            j2 = icsd(n,5)

            call getkeq_lewis_rates(gas,n,tab,rkeq)

            rkf  = exp(aeta(n)*tafl - athet(n)*tafi)
            rkfp = exp(aeta(n)*tabl - athet(n)*tabi)
            rkb  = rkfp*rkeq

            xf = xs(i1)*xs(i2)*xs(i3)
            xb = xs(j1)*xs(j2)

            rkfs = cfm(n,1)*rkf
            rkps = cfm(n,1)*rkfp
            rkbs = cfm(n,1)*rkb

            rate = rkbs*xb - rkfs*xf

            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(i3) = ws(i3) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate

          Case(14)
            !write(*,*) icsd(n,1)
            i1 = 6  !icsd(ir,1)
            j1 = 7  !icsd(ir,2)
            j2 = 18 !icsd(ir,3)
            tafi = ti
            tafl = log(tt)

            xf = xs(6)
            xb = xs(7)*xs(18)
      
            cfm_inf = 630000000000000.0
            ea_inf  = 52334.8
      
            do m=1,ns
               cfm_0(m) = 1.0d+14
            enddo
            cfm_0(2)  = 1.0d+14*3.0
            cfm_0(6)  = 1.0d+14*6.5
            cfm_0(13) = 1.0d+14*0.75
            cfm_0(14) = 1.0d+14*1.5
            cfm_0(21) = 1.0d+14*1.0
            cfm_0(17) = 1.0d+14*6.5
            cfm_0(16) = 1.0d+14*0.4
            cfm_0(15) = 1.0d+14*0.4
      
            ea_0  = 4.32769d+04
      
            rkinf = cfm_inf*exp(-ea_inf*tafi) !beta = 0
            rkf   = exp(-ea_0*tafi)
            
            !write(*,*) 'rkinf, rkf'
            !write(*,*) rkinf, rkf
            rkfs = 0.0d0
            do m = 1,ns
              rkfm = cfm_0(m)*rkf
              rkfs = rkfs + rkfm*xs(m)
            enddo
            
            ifac = rkinf**2/(rkinf+rkfs)**2
           
            mfac = rkfs**2/(rkinf+rkfs)**2
            bfac = rkinf*rkfs/(rkinf+rkfs)
            rate = -bfac*xf
      
            !write(*,*) xf, xb
            !write(*,*) 'ifac, bfac, mfac, rate'
            !write(*,*) ifac, bfac, mfac, rate
            
            ws(i1) = ws(i1) + rate
            ws(j1) = ws(j1) - rate
            ws(j2) = ws(j2) - rate
            
          
            !write(*,*) ' 1   drt, aa(i1,i1), smwi(i1)'
            !write(*,*) drt, aa(i1,i1), smwi(i1)
           !write(*,*) ' 1   drt, aa(i1,i1), smwi(i1)'
            !write(*,*) drt, aa(i1,i1), smwi(i1)
      
            !write(*,*) 'ws(i1),ws(j1),ws(j2)    '
            !write(*,*) ws(i1),ws(j1),ws(j2)    
      
            !write(*,*) 'dkfu,dfkm'
            !write(*,*) dkfu,dkfm
            !!RECOM!IF ( ir == 40 ) THEN
            !write(*,*) 'solving reaction 40'
      
            i1 = 7  !icsd(ir,1)
            i2 = 18 !icsd(ir,2)
            j1 = 6  !icsd(ir,3)
            tafi = ti
            tafl = log(tt)
            !write(*,*) i1,i2,j2
      
            xf = xs(7)*xs(18)
            xb = xs(6)
      
            cfm_inf = 5200000000.0
            ea_inf  = -659.2
      
            do m=1,ns
               cfm_0(m) = 8.2500d+08
            enddo
            cfm_0(2)  = 8.2500d+08*3.0
            cfm_0(6)  = 8.2500d+08*6.5
            cfm_0(13) = 8.2500d+08*0.75
            cfm_0(14) = 8.2500d+08*1.5
            cfm_0(21) = 8.2500d+08*1.0
            cfm_0(17) = 8.2500d+08*6.5
            cfm_0(16) = 8.2500d+08*0.4
            cfm_0(15) = 8.2500d+08*0.4
      
            ea_0  = -9.71717d+03
      
            rkinf = cfm_inf*exp(-ea_inf*tafi) !beta = 0
            rkf   = exp(-ea_0*tafi)
            
            rkfs = 0.0d0
            do m = 1,ns
              rkfm = cfm_0(m)*rkf
              rkfs = rkfs + rkfm*xs(m)
            enddo
            
            ifac = rkinf**2/(rkinf+rkfs)**2
           
            mfac = rkfs**2/(rkinf+rkfs)**2
            bfac = rkinf*rkfs/(rkinf+rkfs)
            rate = -bfac*xf
            
            !write(*,*) xf, xb
            !write(*,*) 'ifac, bfac, mfac, rate'
            !write(*,*) ifac, bfac, mfac, rate
            ws(i1) = ws(i1) + rate
            ws(i2) = ws(i2) + rate
            ws(j1) = ws(j1) - rate
            !write(*,*) rate
          Case(9)                              ! General, unbalanced reactions - RESERVED
            word = 'General'
            goto 901

        End Select

        !rates(n,i) = rate ! For rate of each reaction (mole/m3/sec?)
      
        Enddo ! End Reaction Loop

        do m = 1,ns ! Uncomment for species formation rates
            rates(m,i) = gas%smw(m)*ws(m)*gas%hs(m)! Units are confusing (kg/m3/sec?)-nowj/m3/s
            rates(ns+1,i) = rates(ns+1,i) + gas%smw(m)*ws(m)*gas%hs(m)
        enddo

      Enddo ! End cell loop
      write(*,*) 'made it to here'

      ! -- Cleanup
      !Deallocate(t,p,g,a,hf,cv,tv)
      !Deallocate(rhos)
      Return

 901  Write(*,fmt='(1x,3a)') '***',word,' reacts not supported.'
      Write(*,*) '*** Error: Problem in subroutine chem. Exiting...'
      Stop




      end function exv_rates
c***********************************************************************
c **** Copy of subroutine to compute Keq from Lewis fits for postpar'
c      Adapted from US3D v0.5.2mdb
c      
c      Written by:    Michael Barnhardt  (michael.d.barnhardt@nasa.gov)
c      Created   :    09.02.11
c      Modified  :    Graham Candler 07.02.12: changed lewis to CEA form
c      Modified  :    Nick Gibbons   18.05.15: new inputs + F90 reformat
c***********************************************************************
      Subroutine getkeq_lewis_rates(gas,jr,t,rkeq)

      Use GASPROPS,only: gas_ugcon
      Use GAS_T_MODULE, only : gas_t

      Use models, only: rkeq_ex_limiter
      Use chemistry
      !Use simvars, only: rkeq_ex_limiter
      Implicit none

      Real(8), intent(in) :: t
      Integer, intent(in) :: jr
      Type(gas_t), intent(in) :: gas
      Real(8), intent(out) :: rkeq

      Integer :: n,nn,jj,idir
      Real(8) :: tt,ti,ti2,ti3,tln,t2,t3,t4,rr,thd
      Real(8) :: l1,l2,l3,dl2,dl3
      Real(8) :: a(9),exx
      
      Thd= 1.0d0/3.0d0

      ! This routine originally did some weird stuff with the input
      ! tt that caused some undefined behaviour with the GNU
      ! optimising compiler. I've changed the formatting to be a
      ! bit more sensible and to match the F90 feel of the newer
      ! parts of the code. NNG
      Tt  = t                     
      If(tt .gt. 2.0d4) tt = 2.0d4  ! limit T used in Lewis fits

      T2  = tt*tt
      T3  = t2*tt
      T4  = t3*tt
      Ti  = 1.0d0/tt
      Ti2 = ti*ti
      Ti3 = ti2*ti
      Tln = log(tt)

      If(tt .le. 1.0d3) then
        jj = 1
      Else if(tt .gt. 1.0d3 .AND. tt .le. 6.0d3) then
        jj = 2
      Else
        jj = 3
      Endif

      Exx = 0.0d0
      Do n = 1,gas%ns
        if (.not.(gas%lcrs(jr,1,n).or.gas%lcrs(jr,2,n))) cycle
        a(:)  = gas%lewis(n,jj,:)

        l2  =      -a(1)*ti2 +     a(2)*ti*tln +        a(3)    +
     &        0.5d0*a(4)*tt  + thd*a(5)*t2     + 0.25d0*a(6)*t3 +
     &        0.2d0*a(7)*t4  +     a(8)*ti
        l3  = -0.5d0*a(1)*ti2 -       a(2)*ti +     a(3)*tln +
     &               a(4)*tt  + 0.5d0*a(5)*t2 + thd*a(6)*t3  +
     &        0.25d0*a(7)*t4  +       a(9)

        exx = exx + (gas%scrs(jr,1,n) + gas%scrs(jr,2,n))*( l2  - l3)
      Enddo

      Rr   = (1.0d5*ti/gas_ugcon)**(-gas%snu(jr))
      Exx  = max(-rkeq_ex_limiter,min(exx,rkeq_ex_limiter))
      Rkeq = rr*exp(exx)

      Return
      End Subroutine getkeq_lewis_rates
!
c  ***********************************************************************************
c     Subroutine to build user-defined postprocessing functions
c
c     Input
c     -----
c     pob       - A postprocessor object
c
c     Output
c     ------
c     ier    - Returns nonzeo on an error
c
c  ***********************************************************************************

      Subroutine my_user_main_postpar(pob,ier)
      use POB_MODULE
      use POST_MODULE
      !use chemrates
      !use models, only : gas
      Implicit none
      type(pob_t), intent(INOUT) :: pob
      integer, intent(OUT) :: ier

      logical, dimension(pob%nsv)        :: mask
      integer :: nexv,ns,n,nr
      character(LEN=10) :: numstr 
      character(LEN=10),allocatable,dimension(:) :: exv_names

      ier= 0

      ! Also change the chemrates.f behaviour on line 356
      ns= pob%ns   ! Comment in this for species formation
      nexv = ns + 1
      !nr = pob%gas%nr ! or this for reaction rates 
      !nexv = nr

      allocate(exv_names(nexv))

      ! -- Set TRUE for all solution variables you want in the output
      mask(:)= .false.                       ! Everything excluded by default

      ! -- Set names of additional variables to add to output
      exv_names(:)(:)= ' '

      do n = 1,ns
         exv_names(n) = 'w'//trim(pob%gas%sp_list(n))
      enddo
      
      exv_names(ns+1) = 'H'

      !do n = 1,nr
      !   if(n.lt.10) then
      !      write(unit=numstr,fmt='(i1)') n
      !   else
      !      write(unit=numstr,fmt='(i2)') n
      !   endif
      !   exv_names(n) = 'r'//trim(numstr)
      !enddo

      !call pob_write_fvs_ex(pob,'rates.vtk',nexv,exv_names,exv_rates,ier,mask=mask)      
      call pob_write_fvs_ex(pob,'rates.plt',nexv,exv_names,exv_rates,ier,mask=mask)      

      deallocate(exv_names)
      ier = 0
      return
      !Implicit none
      !type(pob_t), intent(INOUT) :: pob
      !integer, intent(OUT)       :: ier

      !integer, parameter                 :: nexv=6
      !logical, dimension(pob%nsv)        :: mask
      !character(LEN=10), dimension(nexv) :: exv_names

      !! -- Set TRUE for all solution variables you want in the output
      !mask(:)= .false.                       ! Everything excluded by default

      !!mask(pob%vid_r)= .true.               ! Include density
      !!mask(pob%vid_u:pob%vid_u+2)= .true.   ! Include velocity components
      !!mask(pob%vid_t)= .true.               ! Include temperature

      !! -- Set names of additional variables to add to output
      !exv_names(:)(:)= ' '

      !exv_names(1)= 'um'
      !exv_names(2)= 'vm'
      !exv_names(3)= 'wm'

      !exv_names(4)= 'up'
      !exv_names(5)= 'vp'
      !exv_names(6)= 'wp'

      !call pob_write_fvs_ex(pob,'stats.plt',nexv,exv_names,exv_stats,ier,mask=mask)

      !return
      end subroutine my_user_main_postpar

      end module my_functions

c  ***********************************************************************************
c     User initialization routine
c  ***********************************************************************************
      subroutine user_initialize(component,debug,ier)
      use POST_MODULE, only : user_main_postpar
      use my_functions, only : my_user_main_postpar
      Implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier

      ier= 0

      user_main_postpar => my_user_main_postpar

      return
      end subroutine user_initialize

