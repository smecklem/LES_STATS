cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This file provides a means of allowing the users of US3D to compile
c their own routines and to access the data-structure that is already
c in place, as well as to create their own. In principle, this part of
c the code allows for the operation of an entirely different algorithm
c that is in sync with the solver of US3D and is MPI-aware.
c
c The infrastructure consists of several routines and a data module.
c There is a routine that is called after US3D has fully initialized,
c one that is executed before anything is done by the solver at every
c time-step, one that is executed after all operations of the solver
c within a time-step have been completed, one that is executed every
c "nplot" steps if nplot>0 and, finally, one that is called
c before the MPI is shut down by the solver and the code exits.
c
c The routines may access the available modules by uncommenting lines
c and the user is responsible for the smooth operation of all routines.
c within this file. If the need to use a file arises, a large handle
c number should be used to avoid conflicts, or use the IO_FUNCS module
c and make a call to io_find_unused_lun(lun).
c
c The default version of this file is used during compilation and it
c is compiled to a dynamically-loaded library that gets installed into
c the $US3D_HOME directory.  The user may build their own copy of this
c library by running the command us3d-build-user in the current working
c directory.  If a copy of the library is found in the working directory
c at runtime, it will be used instead of the default version.
c
c     US3D modules that you might want to use
c     ---------------------------------------
c     Use MPI
c     Use HDF5
c     Use IO_FUNCS
c     Use TXT_FUNCS
c     Use US3D_EXTRAS
c     Use POB_MODULE
c     Use POST_MODULE
c
c     Use mpivars
c     Use switches
c     Use sizing
c     Use constants
c     Use models
c     Use connect
c     Use geometry
c     Use simvars
c     Use chemistry
c     Use turbulence
c     Use userdata
c
c Ioannis Nompelis <nompelis@arm.umn.edu>                    04.09.2006
c Modified by Heath B. Johnson                Last modified: 2013-10-30
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include "us3d_header.h"

! -- Uncomment this if you want to use the MoVing Grids and Dynamics
! -- module and user subroutines.   You can find this file in the
! -- $US3D_HOME/src directory.  You can either #include this file from
! -- your working directory, or copy and paste the subroutines.
#include "us3d_user_mvg.f"

! -- Uncomment this if you want to use additional special-purpose
! -- user subroutines.  You can find this file in the $US3D_HOME/src
! -- directory.  You can either #include this file from your working
! -- directory, or copy and paste the subroutines you want to use.
#include "us3d_user_special.f"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Module used for data structures required by user-define calls
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       MODULE userdata

       implicit none
       save

       END MODULE userdata

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to perform some startup functions for anything that
c **** might be needed in user subroutines.  This routine is called
c **** after precomp.  It is called before init and user_init.  This
c **** would be a good place to initialize user-defined data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine user_startup(ier)

      Implicit none
      integer, intent(OUT) :: ier

      ier= 0

c *** execution may go here, as well as calls to other routines
c *** this is a typical place for allocating user-defined data storage

      return
      end subroutine user_startup

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** This routine is called at the end of startup, after the main
c **** init routine.  Here you can initialize the flow and boundary
c **** conditions for example.  This routine is also called after
c **** grid tailoring in case you want to reinitialize the flow in
c **** some special way.
c
c     Input
c     -----
c     cmod  - Calling module
c     iop   - Operation, passed from calling module
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_init(cmod,iop)

      Implicit none

      Character(LEN=*), intent(IN) :: cmod
      Integer, intent(IN)          :: iop

      select case (cmod)

      ! --- If called from the solver initialization.  iop= iconr
      case ('startup')
         !select case (iop)
         !end select

      ! --- If called from the tailoring routine.  iop= igti(ipass)
      case ('tailor')
         !select case (iop)
         !end select

      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed before anything is done by the solver
c **** at every timestep.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_pre

      Implicit none

c *** execution may go here, as well as calls to other routines

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed after everything is done by the
c **** solver within a timestep.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_post

      Implicit none

c *** execution may go here, as well as calls to other routines


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every "nplot" steps if nplot>0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_inter

      Implicit none

c *** execution may go here, as well as calls to other routines

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every time step before solution is
c **** updated.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_preupdate

      Implicit none


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every time step just after solution
c **** is updated.
c ****
c **** Output
c **** ------
c **** new   - Whether new output was generated by this routine and
c ****         therefore whether new convergence header lines should
c ****         be written to the screen.
c ****
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_postupdate(new)

      Implicit none

      logical, intent(OUT) :: new

      new= .false.

      return
      end Subroutine user_postupdate
      
c  ***********************************************************************************
c     Subroutine that is executed after calcrate.  Here all gradients have been
c     computed, and exchanges made.  This would be a good place to add source terms.
c
c     Input
c     -----
c     iform - Whether to store the RHS.  This is the variable iform that is passed
c             into calcrate.
c
c     Output
c     ------
c
c  ***********************************************************************************
      subroutine user_postcalcrate(iform)

      implicit none
      
      integer, intent(IN) :: iform
      
      return
      end subroutine user_postcalcrate

c  ***********************************************************************************
c     Subroutine to explicitly set the inviscid wall flux for the US3D boundary
c     condition tap US3D_EXPINVFLUX = 2.  This routine is called from inside the
c     US3D calcrate routine.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     ne    - The number of equations, and the size of the first index of fxi
c     j1,j2 - The face range for this group
c
c     Output
c     ------
c     fxi   - Explicitly-set inviscid wall fluxes, fxi(1:ne,j1:j2)
c
c  ***********************************************************************************
      subroutine user_expinvwallflux(ibc,ibt,ne,j1,j2,fxi)

      implicit none

      integer, intent(IN) :: ibc,ibt,ne,j1,j2
      real(8), dimension(ne,j1:j2), intent(OUT) :: fxi

      ! -- Explicitly set inviscid boundary face fluxes here.  Loop over face
      ! -- range j=j1:j2, and store in flux indices fxi(:,j1:j2).

      return
      end subroutine user_expinvwallflux

c  ***********************************************************************************
c     Subroutine to explicitly set the inviscid wall flux for the US3D boundary
c     condition tap US3D_EXPINVFLUX = 2.  This routine is called from inside the
c     US3D calcrate routine.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c     igrad - The value of igrad that was passed into subroutine gradints
c
c     Output
c     ------
c     grad  - Store values in module variable grad(:,:,:)
c
c  ***********************************************************************************
      subroutine user_gradbc(ibc,ibt,j1,j2,igrad)

      use sizing, only : ngr
      use simvars, only : grad

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2,igrad

      ! -- Explicitly set ghost cell gradients here.  Loop over face
      ! -- range j1:j2, and store in grad indices ii: grad(1:3,1:ngr,ii).

      return
      end subroutine user_gradbc

c  ***********************************************************************************
c     Subroutine to perform any required boundary condition initialization for
c     user-defined boundary conditions.  This routine should have a select case on
c     the variable ibt, and then loops over faces j1:j2.  If a problem is starting
c     fresh, then this routine might set wall solution variables in bvar.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c
c  ***********************************************************************************
      subroutine user_initbc(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Perform a select case on ibt, then loop over face range j1:j2, set values
      ! -- of all variables bvar on first initialization.  This might be skipped on
      ! -- problem restarts.

      return
      end subroutine user_initbc

c  ***********************************************************************************
c     Subroutine to set general user-specified boundary conditions.  This routine is
c     called from inside the US3D expinvbc routine for any values of ibt not built-in.
c
c     When setting the inviscid BC's, you should set all the variables in the ghost
c     cells such that averaging to the boundary gives the right inviscid condition
c     at that boundary.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c
c  ***********************************************************************************
      subroutine user_expinvbc(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Loop over face range j1:j2, set values of all variables in ghost cell ii

      return
      end subroutine user_expinvbc

c  ***********************************************************************************
c     Subroutine to set general user-specified boundary conditions.  This routine is
c     called from inside the US3D visflux routine for any values of ibt not built-in.
c
c     When setting viscous BC's, you should set the variables in the ghost cells
c     such that averaging to the boundary faces gives the correct viscous condition
c     at the face.  You also need to set the values of rmu, rkap, and rkav to the
c     values that give the right average value at the boundary if it was not already
c     handles in user_walltrans.  In addition to setting the solution, you should
c     set dw2 if needed, and the solution variables in bvar (1:ne) should be set to
c     the primitive variable solution at the face.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c
c  ***********************************************************************************
      subroutine user_visfluxbc(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Loop over face range j1:j2, set values of all variables in ghost cell ii
      ! -- if they are different from the values set for inviscid flux calculation,
      ! -- also set values in bvar

      return
      end subroutine user_visfluxbc

c  ***********************************************************************************
c     Subroutine to set general user-specified boundary conditions.  This routine is
c     called from inside the US3D implicitbc routine for any values of ibt not built-
c     in.
c
c     In this routine, you should set the terms in the at matrices relating the plus
c     and minus fluxes at face j with the solutions at i and ii.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c
c  ***********************************************************************************
      subroutine user_implicitbc(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Loop over face range j1:j2 and modify at(:,:,i) and at(:,:,ii)

      return
      end subroutine user_implicitbc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every time that dataio is called,
c **** so that the user can also perform any required data IO.  Inside
c **** US3D, HDF5 is opened before calling this routine and closed
c **** afterwards, so it does not need to be opened and closed here.
c **** If you want to add anything into the solver restart file, you can
c **** use this routine.  Uncomment the lines below to access the file.
c ****
c **** Input
c **** -----
c **** cio   - Set to either 'r' or 'w' for reading or writing,
c ****         depending on how dataio was called.
c **** fname - The US3D solution filename
c **** irnum - The current run number used for reading or writing
c ****
c **** Output
c **** ------
c **** ier   - Return nonzero on an error
c ****
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_dataio(cio,fname,irnum,ier)

      use mpivars

      Implicit none
      Integer, intent(IN)          :: irnum
      Integer, intent(OUT)         :: ier
      character(LEN=1), intent(IN) :: cio
      character(LEN=*), intent(IN) :: fname

      ier= 0

      return
 999  if (id==0) write(6,*) '*** Error in subroutine user_dataio'
      ier= 1
      return
      end

c  ***********************************************************************************
c     Subroutine to be executed after the CFD solver update to write convergence data
c
c     Output
c     ------
c     ncv   - (Optional) Number of convergence variables
c     cvars - (Optional) Convergence variable names
c     cvals - (Optional) Values of convergence variables
c
c  ***********************************************************************************
      subroutine user_conv(ncv,cvars,cvals)

      implicit none

      integer, intent(OUT), optional :: ncv
      real*8, intent(OUT), dimension(:), optional :: cvals
      character(LEN=*), dimension(:), intent(OUT), optional :: cvars

      ! --- Default action.  No additional variables
      if (present(ncv)) then
         ncv= 0
      endif

      return
      end subroutine user_conv

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed from the US3D finalize routine.  Use
c **** the parameter iop to determine the stop condition.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_finalize(iop)

      Implicit none

      Integer, intent(IN) :: iop

      select case (iop)
      case (0)   ! Normal completion

      case (1)   ! NaN condition

      case (2)   ! Abort condition, code will call STOP

      end select

      return
      end

c **********************************************************************
c     Subroutine to compute mixture viscosity, thermal conductivity,
c     vibrational energy thermal conductivity at a number of points.
c     This routine is called from the gas property routines for
c     unknown values of icase.
c
c     Input
c     -----
c     gas        - A gas_t object
c     icase      - Method to use.  Comes from input file ivmod.
c     rhos(1:ns) - Species densities
c     t          - translational-rotational temperature
c     tv         - Vibrational temperature
c     ib         - Beginning index
c     ie         - Ending index
c
c     Output
c     ------
c     rmu    - Mixture viscosity
c     rkap   - Mixture translational-rotational thermal conductivity
c     rkav   - Mixture vibrational energy thermal condictivity
c
c     To use this routine, make a loop over cell indices ib:ie and
c     use the solution variables (rhos,t,tv) and the gas properties
c     object to compute the values of rmu, rkap, and rkav that you
c     want to model.  If you wish, use the input variable icase to
c     select between different models, or program just one.
c
c **********************************************************************
      Subroutine user_trans(gas,rhos,t,tv,ib,ie,rmu,rkap,rkav)
      use GAS_T_MODULE
      Implicit none
      type(gas_t), intent(IN) :: gas
      integer, intent(IN) ::ib,ie
      real*8, dimension(:), intent(IN) :: t,tv
      real*8, dimension(:,:), intent(IN) :: rhos
      real*8, dimension(:), intent(OUT) :: rmu,rkap,rkav

c --- Program the gas mixture transport properties that you want here
      ! Diffusivity
      Select Case(gas%idmod)

      Case default
         !!! - Unknown transport property method
         write(6,*) '*** Unknown user diffusivity method= ',gas%idmod
         stop
      End Select

      ! Viscosity
      Select Case(gas%ivmod)

      Case default
         !!! - Unknown transport property method
         write(6,*) '*** Unknown user viscosity method= ',gas%ivmod
         stop
      End Select

      return
      end subroutine user_trans

c  ***********************************************************************************
c     Subroutine that is executed in visflux() to set the conditions at the wall
c     (inside ghosts) for transport properties calculations.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c     iwall - The walue of iwall for this BC
c
c  ***********************************************************************************
      Subroutine user_walltrans(ibc,ibt,j1,j2,iwall)

      use mpivars

      Implicit none
      Integer, intent(IN) :: ibc,ibt,j1,j2,iwall

c
c *** treatment of the boundary variables for properties calculations
c
      Select Case(iwall)

      Case Default
         if(id.eq.0) print*,'Problematic assignment of wall conditions, iwall= ',iwall
         if(id.eq.0) print*,'Stopping in subroutine user_walltrans'
         call finalize(2)
      End Select

      return
      end

c  ***********************************************************************************
c     Subroutine that is executed in visflux() to treat the wall fluxes
c     in a non-trivial manner.  This routine is called for wall boundary
c     conditions where iwall is not handled internally.
c
c     Input
c     -----
c     ibc   - The boundary condition number in the bcs structure
c     ibt   - The boundary type from the zone definition: zdefs(6,iz)
c     j1,j2 - The face range for this group
c     iwall - The walue of iwall for this BC.
c
c  ***********************************************************************************

      Subroutine user_viswallbc(ibc,ibt,j1,j2,iwall)

      Use mpivars
      Use models

      Implicit none
      integer, intent(IN) :: ibc,ibt,j1,j2,iwall

c
c *** treatment of the wall boundary conditions
c
      Select Case(iwall)

      Case Default
         if(id.eq.0) print*,'Problematic assignment of wall conditions, iwall= ',iwall
         if(id.eq.0) print*,'Stopping in subroutine user_viswallbc'
         call finalize(2)
      End Select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Subroutine that is called in calcrate() to set the switch for
c     whether or not to add dissipative fluxes when ibase>0.  This
c     routine is called if an unknown value of idiss_g is specified
c     in the solver input file.  The switch needs to be set for each
c     equation, 1:ne.
c
c     Input
c     -----
c     k    - Current index in jmaster
c     j    - Current face
c     i    - Left neighbor cell
c     ii   - Right neighbor cell
c     ml   - Mach number to the left
c     mr   - Mach number to the right
c     divl - Divergence to the left
c     divr - Divergence to the right
c
c     Output
c     ------
c     dfac - The switch.  Valid values are 0.0-1.0, (0.0=off, 1.0=on)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_dswitch(k,j,i,ii,ml,mr,divl,divr,dfac)

      Use sizing, only   : ne
      Use mpivars, only  : id
      Use switches, only : ibase,idiss_glob_on

      Implicit none

      Integer, intent(IN)  :: k,j,i,ii
      Real(8), intent(IN)  :: ml,mr,divl,divr
      Real(8), dimension(ne), intent(OUT) :: dfac

      if (id==0) then
         write(6,*) '*** user_dswitch called for ibase,idiss_g= ',ibase,idiss_glob_on
         write(6,*) '*** without this user routine being built'
      endif
      stop

      ! -- Set dfac here (1.0 = on)
      dfac(:)= 1.0d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Subroutine that is called in calcrate() to set the inviscid fluxes
c     when ibase=4.
c
c     Input
c     -----
c     k    - Current index in jmaster
c     j    - Current face
c     i    - Left neighbor cell
c     ii   - Right neighbor cell
c
c     Output
c     ------
c     fxi - The fluxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine user_flux(k,j,i,ii,
     &                     ns, nv, nt, ne, ngr, rsl, rsr, rl, rr,
     &                     ul, ur, vl, vr, wl, wr, tl, tr, tvl, tvr, evl, evr,
     &                     pl, pr, rntl, rntr, grdxl, grdxr, grdxl2, grdxr2, 
     &                     gsp, cvs, hs, sss, ssx, fxi)

      Use mpivars, only  : id

      Implicit none

      Integer, intent(IN)  :: k,j,i,ii
      Integer, intent(IN)  :: ns,nv,nt,ne,ngr
      Real(8), intent(IN)  :: rl,rr,ul,ur,vl,vr,wl,wr,pl,pr
      Real(8), intent(IN)  :: tl,tr,tvl,tvr,evl,evr,rntl,rntr,sss
      Real(8), intent(IN)  :: grdxl(ngr),grdxr(ngr),grdxl2(ngr),grdxr2(ngr)
      Real(8), dimension(ns), intent(IN)  :: rsl,rsr,gsp,hs,cvs
      Real(8), dimension(3), intent(IN)   :: ssx

      Real(8), dimension(ne), intent(OUT) :: fxi

      if (id==0) then
         write(6,*) '*** user_flux called without this user routine being built'
      endif
      stop

      ! -- Build fxi here
      fxi(:)= 0.0d0

      return
      end

c  ***********************************************************************************
c     Use this routine if you want to change the blending function that is used to
c     modify Jacobians for Crank-Nicolson time integration
c  ***********************************************************************************
      subroutine user_cn_blend()
      implicit none

      write(6,*) '*** user_cn_blend called without this user routine being built'
      stop

      return
      end subroutine user_cn_blend

c  ***********************************************************************************
c     Subroutine to code in any custom management routines for US3D.  When triggered,
c     the value of lun is set to the current logical unit reading the management
c     section of the input file.  You can use this to read any additional parameters
c     that you want to act on, similar to the CFL list.  Only processor zero should
c     do the reading.
c
c     Input
c     -----
c     lun - On processor 0, the LUN to the management file
c
c     Output
c     ------
c     return_paused - Whether to return paused or to continue the simulation
c     ier  - Returns nonzeo on an error
c
c  ***********************************************************************************

      Subroutine user_manage(lun,return_paused,ier)
      Implicit none
      integer, intent(IN)  :: lun
      integer, intent(OUT) :: ier
      logical, intent(OUT) :: return_paused

      ier = 0
      return_paused= .false.

      write(6,*) '-- No user management routines defined'

      return
      end subroutine user_manage

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

      Subroutine user_postpar(pob,ier)
      use POB_MODULE
      use POST_MODULE
      Implicit none
      type(pob_t), intent(INOUT) :: pob
      integer, intent(OUT) :: ier

      ier= 0

      write(6,*) '-- No user postpar routines defined'

      return
      end subroutine user_postpar

c  ***********************************************************************************
c     Subroutine to give a starting point for all kinds of US3D hacking.
c
c     Input
c     -----
c     opts   - Options string passed in from the command line
c     debug  - Logical flag whether debugging is on
c
c     Output
c     ------
c     ier    - Returns nonzeo on an error
c
c  ***********************************************************************************

      Subroutine user_hacks(opts,debug,ier)
      Implicit none
      logical, intent(IN) :: debug
      character(LEN=*), intent(IN) :: opts
      integer, intent(OUT) :: ier

      ier= 0

      write(6,*) '-- No user hacking routines defined'

      return
      end subroutine user_hacks

c  ***********************************************************************************
c     Subroutine to return the version of US3D with which this version of the user
c     library was compiled.  This is just a check to make sure that you are not
c     using an old version of libus3d_user.so when you try to run.
c
c     Input
c     -----
c
c     Output
c     ------
c     cvers - Returns the version of US3D at compilation time
c
c  ***********************************************************************************

      Subroutine user_vcheck(cvers)
      use TXT_FUNCS
      Implicit none
      character(LEN=*), intent(OUT) :: cvers

      call txt_str_copy(US3D_MAIN_VERSION,cvers)
      
      return
      end subroutine user_vcheck
