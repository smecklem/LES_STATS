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
!#include "us3d_user_mvg.f"

! -- Uncomment this if you want to use additional special-purpose
! -- user subroutines.  You can find this file in the $US3D_HOME/src
! -- directory.  You can either #include this file from your working
! -- directory, or copy and paste the subroutines you want to use.
#include "us3d_user_special.f"
#include "stats_module.f"
#include "combined-user.f"

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Module used for data structures required by user-define calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       MODULE userdata_stats

       implicit none
       integer,parameter :: istats=1 ! flag for stats
       save

       END MODULE userdata_stats

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Module used for all the LES stats stuff
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       MODULE userdata_stats_routines

      implicit none
      save

      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to perform some startup functions for anything that
c **** might be needed in user subroutines.  This routine is called
c **** after precomp.  It is called before init and user_init.  This
c **** would be a good place to initialize user-defined data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine my_user_startup(ier)
      use userdata_stats
      use stats
      Implicit none
      integer, intent(OUT) :: ier

      ier= 0
      if (istats.gt.0) call init_stats()

      return
      end subroutine my_user_startup

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
      use MPI
      Implicit none

      Character(LEN=*), intent(IN) :: cmod
      Integer, intent(IN)          :: iop
      Integer :: ier
      select case (cmod)

      ! --- If called from the solver initialization.  iop= iconr
      case ('startup')
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

      Subroutine user_main_pre
      Implicit none
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed after everything is done by the
c **** solver within a timestep.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_main_post
      implicit none
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every "nplot" steps if nplot>0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_main_inter

      Implicit none

c *** execution may go here, as well as calls to other routines

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every time step before solution is
c **** updated.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine user_update_pre

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

      Subroutine my_user_update_post(new)
      use stats
      use userdata_stats
      Implicit none
      
      logical, intent(OUT) :: new

      if(istats.gt.0) call update_stats(istats)   ! *** for statistics

      new= .false.

      return
      end Subroutine my_user_update_post
      
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
      subroutine user_calcrate_post(iform)
      implicit none
      
      integer, intent(IN) :: iform

      return
      end subroutine user_calcrate_post

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
      subroutine user_bc_expinv_wallflux(ibc,ibt,ne,j1,j2,fxi)

      implicit none

      integer, intent(IN) :: ibc,ibt,ne,j1,j2
      real(8), dimension(ne,j1:j2), intent(OUT) :: fxi

      ! -- Explicitly set inviscid boundary face fluxes here.  Loop over face
      ! -- range j=j1:j2, and store in flux indices fxi(:,j1:j2).

      return
      end subroutine user_bc_expinv_wallflux

c  ***********************************************************************************
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
      subroutine user_bc_grad(ibc,ibt,j1,j2,igrad)

      use sizing, only : ngr
      use simvars, only : grad

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2,igrad

      ! -- Explicitly set ghost cell gradients here.  Loop over face
      ! -- range j1:j2, and store in grad indices ii: grad(1:3,1:ngr,ii).

      return
      end subroutine user_bc_grad

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
      subroutine user_init_bc(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Perform a select case on ibt, then loop over face range j1:j2, set values
      ! -- of all variables bvar on first initialization.  This might be skipped on
      ! -- problem restarts.

      return
      end subroutine user_init_bc

c  *************************************************************************************
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
      subroutine user_bc_visflux(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Loop over face range j1:j2, set values of all variables in ghost cell ii
      ! -- if they are different from the values set for inviscid flux calculation,
      ! -- also set values in bvar

      return
      end subroutine user_bc_visflux

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
      subroutine user_bc_implicit(ibc,ibt,j1,j2)

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2

      ! -- Loop over face range j1:j2 and modify at(:,:,i) and at(:,:,ii)

      return
      end subroutine user_bc_implicit

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed every time that dataio is called,
c **** so that the user can also perform any required data IO.  HDF5 is
c **** opened before calling this routine and closed afterwards.  If
c **** you want to add anything into the solver restart file, you can
c **** do so here.  Uncomment the lines below to access the file.
c ****
c **** Input
c **** -----
c **** cio   - Set to either 'r' or 'w' for reading or writing,
c ****         depending on how dataio was clled.
c **** fname - The US3D solution filename
c **** irnum - The current run number used for reading or writing
c ****
c **** Output
c **** ------
c **** ier   - Return nonzero on an error
c ****
c **** 13/10/2023 - This is an accumulation of stats_dataio, stats_read, 
c ****              and stats_write bc hdf5 files weren't being written 
c **** Sarah M      properly. It may be destilled into appropriate 
c ****              subfunctions at a later data. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine my_user_dataio(cio,fname,irnum,ier)

      use MPI
      use HDF5
      use H5_EXTRAS
      use US3D_EXTRAS

      use us3d_global
      use us3d_varkinds
      use userdata_stats
      use stats

      use mpivars
      use switches
      use restartio, only : read_us3d_qvals,write_us3d_qvals
      use connect, only : ugrid
      use sizing, only : nel,neg

      Implicit none
      integer, intent(IN)          :: irnum
      integer, intent(OUT)         :: ier
      character(LEN=1), intent(IN) :: cio
      character(LEN=*), intent(IN) :: fname

      integer :: i
      integer :: hdf_err,ihave,irnumg,istat
      integer(HID_T) :: rid
      logical :: new
      type(us3d_sfile) :: sfile

      real*8 :: testval
      real*8, pointer, dimension(:,:) :: dump

      integer :: epg,global_epg
      real*8, allocatable, dimension(:,:) :: myqdum
      integer, allocatable, dimension(:) :: ige

      integer :: have_mydata
      integer, parameter :: nuv=1  ! CHANGE NUMBER OF USER VARIABLES HERE
      character(LEN=200) :: rpath,dpath
      character(40) :: dname
      integer :: olun= 6

      ! -- Set the name of your dataset here
      dname= 'stats-mean'

      ier= 0
      dpath(:)= ' '
      
      if (id==0) write(6,*) '-- Inside user_dataio, case: '//trim(cio)

      do while (id==0)

         ! -- Initialize solution file
         new= .false.
         if (us3d_debug) write(6,*) '-- Initializing solution file: "'//trim(fname)//'"'
         call us3d_sfile_init(fname,new,sfile,ier)
         if (ier/=0) exit

         if (irnum>sfile%nruns) then
            write(6,*) '*** Run ',irnum,' does not yet exist in this file'
            ier= 1
            exit
         endif

         ! ---Open current run group
         irnumg= irnum
         call us3d_sfile_orun(sfile,irnumg,rid,ier,rpath=rpath)
         if (ier/=0) exit

         call h5ex_att_getd(rid,'have_'//trim(dname),0,ihave,ier)
         if (ier/=0) exit

         write(6,*) '-- Have my data? ',ihave
         if (ihave==1) then
            new= .false.
         else
            new= .true.
         endif

         if (cio=='w') then
            ihave= 1
            call h5ex_att_add(rid,'have_'//trim(dname),ihave,ier)
            call h5ex_att_add(rid,'stats_time',stat_time,ier)
            call h5ex_att_add(rid,'unst_time',unst_time,ier)
            call h5ex_att_add(rid,'stats_mean_vars',nmvar,ier)
            if (ier/=0) exit
         endif

         if (cio=='r'.and.ihave==1) then
            call h5ex_att_get(rid,'stats_time',stat_time,ier)
            call h5ex_att_get(rid,'unst_time',unst_time,ier)
            call h5ex_att_get(rid,'stats_mean_vars',nmvar,ier)
            if (ier/=0) exit
         endif

         ! --- Close run group and solution file
         call h5gclose_f(rid, hdf_err)
         call us3d_sfile_cnd(sfile,ier)

         exit
      enddo

      call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 999

      call MPI_BCAST(new,1,MPI_LOGICAL,0,icomw,ier)
      call MPI_BCAST(ihave,1,MPI_INTEGER,0,icomw,ier)
      call MPI_BCAST(rpath,len(rpath),MPI_CHARACTER,0,icomw,ier)

      dpath= trim(rpath)//'/'//trim(dname)

      write(6,*) '-- dpath= "'//trim(dpath)//'"'

      ! --- Add any required IO here
      Select Case(cio)
      Case('w')

         ! -- Write cell-based data.  Assume sized mydata(my_nv,epg) on each processor
         !!epg= nel + neg                   ! This processor number of interior plus shared/ghost
         !!global_epg= ugrid%nc + ugrid%ng  ! number of interior + ghost cells in global grid
         allocate(myqdum(nmvar,nel),ige(nel),STAT=istat)

         if (istat/=0) ier= 1
         call us3d_check_int(icomw,MPI_MAX,ier)
         if (ier/=0) goto 901

         do i= 1,nel
            ige(i)= ugrid%ige(i)

            myqdum(1:nmvar,i) = mean_var(1:nmvar,i)
            !rdum(1:nsvar,i) = stat_var(1:nsvar,i)
         enddo

         ! -- Called with rpath, we expect to open and close the file in us3d_h5data_pw
         call write_us3d_qvals(dpath,myqdum,.false.,new,.false.,ier)
         if (ier/=0) goto 999

         if (allocated(myqdum)) deallocate(myqdum)
         if (allocated(ige)) deallocate(ige)
         if (id==0) write(6,*) '== Successfully wrote "'//trim(dname)//'"'

      Case('r')

         if (ihave==1) then

            !epg= nel + neg                   ! This processor number of interior plus shared/ghost
            allocate(myqdum(nmvar,nel),ige(nel),STAT=istat)

            write(6,*) 'size(qva)= ',size(myqdum)
            write(6,*) 'size(qva)= ',size(ige)

            if (us3d_debug.and.id==0) write(olun,*) 'Reading solution variables in parallel'
            !*** Read solution data globally according to the global map
            ige(1:nel)= ugrid%ige(1:nel)

            dump => read_us3d_qvals(dpath,nuv,.false.,ier,qva=myqdum)
            !dump => read_us3d_qvals(dpath,nuv,.false.,ier,qva=ige)

            if (ier/=0) goto 999

            ! -- Do something with the data you read.  For now we just test that we
            ! -- read it properly.
            do i=1,nel
               testval = myqdum(1,i) - dble(ugrid%ige(i))
               if (abs(testval).gt.1.0d-20) then
                  write(6,*) '*** Failed to properly load data: ',
     &               myqdum(1,i),dble(ugrid%ige(i)),abs(testval)
                  ier= 1
               endif
            enddo
            write(6,*) id,' read all data properly'

            if (allocated(myqdum)) deallocate(myqdum)
            if (allocated(ige)) deallocate(ige)
            if (id==0) write(6,*) '== Successfully read "'//trim(dname)//'"'

         else
            if (id==0) write(6,*) '-- No "'//trim(dname)//'" present for reading'
         endif

      end select

      return
 901  if (id==0) write(olun,*) '*** Error allocating memory in subroutine stats_read'
 999  if (id==0) write(6,*) '*** Error in subroutine my_user_dataio'
      ier= 1
      return
      end subroutine my_user_dataio

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

      integer :: i,n

      ! --- Set the number of user-specified convergence variables
      if (present(ncv)) then
          ncv = 0
      endif

      ! -- Set the names of the variables
      if (present(cvars)) then
      endif

      if (present(cvals)) then
      endif ! end if (present(cvals))

      return
      end subroutine user_conv

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine that is executed from the US3D finalize routine.  Use
c **** the parameter iop to determine the stop condition.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine my_user_finalize(iop)
      use userdata_stats
      use stats
      use mpivars, only : id
      Implicit none

      Integer, intent(IN) :: iop
      integer :: ier

      if (istats.gt.0) call shutdown_stats(ier) 
      select case (iop)
      case (0)   ! Normal completion

      case (1)   ! NaN condition

      case (2)   ! Abort condition, code will call STOP

      end select

      return
      end subroutine my_user_finalize

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
      Subroutine my_user_transport(gas,rhos,t,tv,ib,ie,rmu,rkap,rkav)
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
      end subroutine my_user_transport


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
      Subroutine user_bc_walltrans(ibc,ibt,j1,j2,iwall)

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

      Use sizing, only   : ne,nel,nv
      Use mpivars, only  : id
      Use switches, only : ibase,idiss_glob_on, ncr
      Use simvars, only  : cuv,t,tv
      Use geometry, only : dw, xcf
      Use models, only   : ns

      Implicit none

      Integer, intent(IN)  :: k,j,i,ii
      Real(8), intent(IN)  :: ml,mr,divl,divr
      Real(8), dimension(ne), intent(OUT) :: dfac


      Real*8,parameter     :: eps_ducros=1.0d-12, thresh_ducros=0.75d0
      Real*8,parameter     :: delta_dr = 0.5d0, kappa_dr = 10.0d0

      REAL*8  :: dmach,sw1,sw2,ddfac,ddfac2,rdum,hack_x,hack_y,hack_z,hack_r
      INTEGER :: iii,n

      ! *** Ducros switch
      sw1  = divl**2/(divl**2 + cuv(i)**2  + eps_ducros)
      sw2  = divr**2/(divr**2 + cuv(ii)**2 + eps_ducros)

      sw1  = min(dabs(sw1/thresh_ducros), 1.0d0)
      sw2  = min(dabs(sw2/thresh_ducros), 1.0d0)
      ddfac = max(sw1, sw2)

      ! *** Mach switch
      dmach = dabs(mr-ml)
      if (dmach .le. delta_dr) then
         ddfac2 = (dmach**2 + delta_dr**2)*(1.0d0 -
     &        exp(-kappa_dr*0.5d0*(ml+mr)))
     &        /(2.0d0*delta_dr)
      else if(dmach.gt.delta_dr.and.dmach.le.1.0d0) then
         ddfac2 = dmach
      else
         ddfac2 = 1.0d0
      endif

      if((ml.le.1.0d0.and.mr.ge.1.0d0).or.(mr.le.1.0d0.and.ml.ge.1.0d0)) then
         ddfac = ddfac2
      end if


      dfac(:) = ddfac 
     
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
c     @Curator: Nick N. Gibbons (n.gibbons@uq.edu.au)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine user_flux(k,j,i,ii,
     &                     ns, nv, nt, ne, ngr, rsl, rsr, rl, rr,
     &                     ul, ur, vl, vr, wl, wr, tl, tr, tvl, tvr, evl, evr,
     &                     pl, pr, rntl, rntr, grdxl, grdxr, grdxl2, grdxr2, 
     &                     gsp, cvs, hs, sss, ssx, fxi)

      Use mpivars, only  : id
      Use models, only   : gas
      Use simvars, only  : faci,t,tv
      Use switches, only : ivib, ncr
      Use GASPROPS, only : gas_vibels1
      Use sizing, only   : nel
      Use geometry, only : xcf

      Implicit none
      Integer, intent(IN)  :: k,j,i,ii
      Integer, intent(IN)  :: ns,nv,nt,ne,ngr
      Real(8), intent(IN)  :: rl,rr,ul,ur,vl,vr,wl,wr,pl,pr
      Real(8), intent(IN)  :: tl,tr,tvl,tvr,evl,evr,rntl(2),rntr(2),sss
      Real(8), intent(IN)  :: grdxl(ngr),grdxr(ngr),grdxl2(ngr),grdxr2(ngr)
      Real(8), dimension(ns), intent(IN)  :: rsl,rsr,gsp,hs,cvs
      Real(8), dimension(3), intent(IN)   :: ssx

      Real(8), dimension(ne), intent(OUT) :: fxi

      Real(8), parameter :: beta = 2.0d0/3.0d0
      Real(8), parameter :: alpha = 1.0d0, delta = -2.0d0/45.0d0, gamm = 16.0d0/15.0d0
      Integer,parameter :: order = 6
      Real(8) :: upl,upr,upav,rav,rupa,pav,kav,el,er,gtr,gtl
      Real(8) :: rhosfl(ns),csfl(ns),rfl,ufl,vfl,wfl,pfl,tfl,evfl,tvfl,rntfl(2)
      Real(8) :: rhosfr(ns),csfr(ns),rfr,ufr,vfr,wfr,pfr,tfr,evfr,tvfr,rntfr(2)
      Real(8) :: rflag,AA,BB,AH,BH,rpsi,csmn,csmx,ht,rfs,evn(ns),cvv(ns)
      Real(8) :: hack_x,hack_y,hack_z,hack_r
      Integer :: n,m,nn,idx
      Logical :: hackflag
  
       fxi(:) = 0.0d0  ! VERY IMPORTANT TO INITIALISE THIS (maybe)
       hackflag = .false.

       rfl = 0.0d0
       rfr = 0.0d0
       evfl= 0.0d0
       evfr= 0.0d0
       rntfl= 0.0d0
       rntfr= 0.0d0
       
       select case(order)
       case(4)
           do n=1,ns
             rhosfl(n) = rsl(n) + grdxl(n)*beta
             rhosfr(n) = rsr(n) + grdxr(n)*beta
             if(rhosfl(n) .lt. 0.0d0) rhosfl(n) = rsl(n)
             if(rhosfr(n) .lt. 0.0d0) rhosfr(n) = rsr(n)
             rfl       = rfl + rhosfl(n)
             rfr       = rfr + rhosfr(n)
           enddo
           csfl(:) = rhosfl(:)/rfl
           csfr(:) = rhosfr(:)/rfr

           ufl = ul + grdxl(ns+1)*beta
           ufr = ur + grdxr(ns+1)*beta
           vfl = vl + grdxl(ns+2)*beta
           vfr = vr + grdxr(ns+2)*beta
           wfl = wl + grdxl(ns+3)*beta
           wfr = wr + grdxr(ns+3)*beta

           tfl = tl + grdxl(ns+4)*beta
           tfr = tr + grdxr(ns+4)*beta
           if(tfl .lt. 0.0d0) tfl = tl
           if(tfr .lt. 0.0d0) tfr = tr

           do n=1,nv
               evfl= evl + grdxl(ns+5+n)*beta
               evfr= evr + grdxr(ns+5+n)*beta
               if(evfl .lt. 0.0d0) evfl = evl
               if(evfr .lt. 0.0d0) evfr = evr
           enddo

           do n=1,nt
               nn = ns+5+nv+nv+n 
               rntfl(n) = rntl(n) + grdxl(nn)*beta
               rntfr(n) = rntr(n) + grdxr(nn)*beta
               if(rntfl(n) .lt. 0.0d0) rntfl(n) = rntl(n)
               if(rntfr(n) .lt. 0.0d0) rntfr(n) = rntr(n)
           enddo
       case(6)
           do n=1,ns
             rhosfl(n) = rsl(n)*alpha + grdxl(n)*gamm + grdxl2(n)*delta
             rhosfr(n) = rsr(n)*alpha + grdxr(n)*gamm + grdxr2(n)*delta
             if(rhosfl(n) .le. 0.0d0) rhosfl(n) = rsl(n)
             if(rhosfr(n) .le. 0.0d0) rhosfr(n) = rsr(n)
             rfl = rfl + rhosfl(n)
             rfr = rfr + rhosfr(n)
           enddo
           csfl(:) = rhosfl(:)/rfl
           csfr(:) = rhosfr(:)/rfr

           ufl = ul*alpha + grdxl(ns+1)*gamm + grdxl2(ns+1)*delta
           ufr = ur*alpha + grdxr(ns+1)*gamm + grdxr2(ns+1)*delta
           vfl = vl*alpha + grdxl(ns+2)*gamm + grdxl2(ns+2)*delta
           vfr = vr*alpha + grdxr(ns+2)*gamm + grdxr2(ns+2)*delta
           wfl = wl*alpha + grdxl(ns+3)*gamm + grdxl2(ns+3)*delta
           wfr = wr*alpha + grdxr(ns+3)*gamm + grdxr2(ns+3)*delta

           tfl = tl*alpha + grdxl(ns+4)*gamm + grdxl2(ns+4)*delta
           tfr = tr*alpha + grdxr(ns+4)*gamm + grdxr2(ns+4)*delta
           if(tfl .le. 0.0d0) tfl = tl
           if(tfr .le. 0.0d0) tfr = tr

           do n=1,nv
               evfl= evl*alpha + grdxl(ns+5+n)*gamm + grdxl2(ns+5+n)*delta
               evfr= evr*alpha + grdxr(ns+5+n)*gamm + grdxr2(ns+5+n)*delta
               if(evfl .le. 0.0d0) evfl = evl
               if(evfr .le. 0.0d0) evfr = evr
           enddo

           do n=1,nt
               nn = ns+5+nv+nv+n
               rntfl(n)= rntl(n)*alpha + grdxl(nn)*gamm + grdxl2(nn)*delta
               rntfr(n)= rntr(n)*alpha + grdxr(nn)*gamm + grdxr2(nn)*delta
               if(rntfl(n) .le. 0.0d0) rntfl(n) = rntl(n)
               if(rntfr(n) .le. 0.0d0) rntfr(n) = rntr(n)
           enddo

       case default 
           write(*,*) "Error in userflux, check order=4,6"
           stop
       end select

       upl     = ufl*ssx(1) + vfl*ssx(2) + wfl*ssx(3)
       upr     = ufr*ssx(1) + vfr*ssx(2) + wfr*ssx(3)
       upav    = 0.5d0*(upl + upr)


       ! *** mass generation fix (band-aid)
       rflag = 0.0d0
       rupa  = 0.0d0
       rav   = 0.0d0

       do n=1,ns
          rpsi = 0.5d0
          csmn = min(csfl(n),csfr(n))
          csmx = max(csfl(n),csfr(n))
         
          if ((csfl(n)+csfr(n)).gt.1.0d-16) then
            
             AA = csfl(n)
             BB = csfr(n)
             AH = AA/(AA+BB)
             BH = BB/(AA+BB)
            
             rflag = sign(1.0d0,(csfr(n) - csfl(n))*upav)
             if(rflag.gt.0.0d0) then
               
                if( (csmx-csmn) .gt. 1.0d-10) then
                   rpsi = 2.0d0*AH*BH*(AA+BB)
                   rpsi = (rpsi-csmn)/(csmx-csmn)
                else
                   rpsi = 0.5d0 
                endif
                rpsi = max(rpsi,0.0d0)
                rpsi = min(rpsi,1.0d0)               
             else
                !rfs = (1.0d0-2.0d0*AH*BH)*(AA+BB)
             endif
          endif
         
          if(csfl(n) .eq. csmn) then
             rhosfl(n) = 2.0d0*(1.0d0-rpsi)*rhosfl(n)
             rhosfr(n) = 2.0d0*rpsi*rhosfr(n)
          else
             rhosfl(n) = 2.0d0*rpsi*rhosfl(n)
             rhosfr(n) = 2.0d0*(1.0d0-rpsi)*rhosfr(n)
          endif
         
          rfs    = 0.5d0*(rhosfl(n) + rhosfr(n))
          rav    = rav  + rfs
          fxi(n) = rfs*upav  ! Set rhos flux here
          rupa   = rupa + fxi(n)
       enddo ! End Mass Gen Fix Section

       evn  = 0.0d0 
       if (ivib==3) call gas_vibels1(gas,0.5*(tfl+tfr),evn,cvv)

       ht   = 0.0d0
       el   = 0.0d0
       er   = 0.0d0
       gtl  = 0.0d0
       gtr  = 0.0d0

       ! This bit edited to use adjusted rhos instead of cs
       do n=1,ns
          ht   = ht + 0.5d0*(rhosfl(n) + rhosfr(n))*(hs(n)+evn(n))
          el   = el + rhosfl(n)*cvs(n)*tfl
          er   = er + rhosfr(n)*cvs(n)*tfr
          gtl  = gtl + rhosfl(n)*gsp(n)
          gtr  = gtr + rhosfr(n)*gsp(n)
       enddo
       pfl = tfl*gtl
       pfr = tfr*gtr
       ht = ht/rav     ! Needs to be adjusted because not cs
       el = el/rav     ! Needs to be adjusted because not cs
       er = er/rav     ! Needs to be adjusted because not cs

       !rav     = 0.5d0*(rfl + rfr)   ! Mass Gen Handles this
       !rupa    = upav*rav            ! Mass Gen Handles this
       pav     = 0.5d0*(pfl + pfr)
       kav     = 0.25d0*(ufl*ufl + vfl*vfl + wfl*wfl
     >                 + ufr*ufr + vfr*vfr + wfr*wfr)

       ! Mass Gen Handles this too
       !do n=1,ns
       !   fxi(n) = 0.5d0*(csfl(n) + csfr(n))*rupa
       !enddo

       fxi(ns+1) = rupa*0.5d0*(ufl+ufr) + pav*ssx(1)
       fxi(ns+2) = rupa*0.5d0*(vfl+vfr) + pav*ssx(2)
       fxi(ns+3) = rupa*0.5d0*(wfl+wfr) + pav*ssx(3)
       fxi(ns+4) = rupa*0.5d0*(ufl*ufr + vfl*vfr + wfl*wfr
     >                         + el + er + evfl + evfr)
     >             +rupa*ht + pav*upav
       do n=1,nv
           fxi(ns+4+n) = 0.5d0*(evfl + evfr)*rupa
       enddo

       !do n=1,nv ! Use MSW for Tv (NNG)
       !  do m=1,ne
       !     fxi(ns+4+n) = fxi(ns+4+n) + aap(ns+4+n,m)*fup(m) + aam(ns+4+n,m)*fum(m)
       !  enddo
       !  fxi(ns+4+n) = fxi(ns+4+n)*faci
       !enddo

       do n=1,nt
           fxi(ns+4+nv+n) = rupa*0.5d0*(rntfl(n) + rntfr(n))
       enddo

      ! MSW already multiplied by face area, so do everything else
      !fxi(1:ns+4) = fxi(1:ns+4)*sss
      !fxi(ns+5+nv:) = fxi(ns+5+nv:)*sss
      fxi(:) = fxi(:)*sss

      return
      end subroutine user_flux

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
c     ier    - Returns nonzeo on an error'
c
c  ***********************************************************************************

      Subroutine my_user_main_postpar(pob,ier)
      use POB_MODULE
      use POST_MODULE
      use stats
      Implicit none
      type(pob_t), intent(INOUT) :: pob
      integer, intent(OUT) :: ier

      logical, dimension(pob%nsv)        :: mask
      integer :: n,nexv,ns
      character(LEN=10) :: numstr
      character(LEN=10),allocatable,dimension(:) :: exv_names

      ier= 0
      ns= pob%ns

      nexv = ns + 13
      allocate(exv_names(nexv))

      ! -- Set TRUE for all solution variables you want in the output
      mask(:)= .false.                       ! Everything excluded by default

      ! -- Set names of additional variables to add to output
      exv_names(:)(:)= ' '

      do n = 1,ns
         if(n.lt.10) then
            write(unit=numstr,fmt='(i1)') n
         else
            write(unit=numstr,fmt='(i2)') n
         endif

         exv_names(n) = 'r'//trim(pob%gas%sp_list(n))//'m'

      enddo

      exv_names(ns+1) = 'um'
      exv_names(ns+2) = 'vm'
      exv_names(ns+3) = 'wm'
      exv_names(ns+4) = 'Mm'
      exv_names(ns+5) = 'ksgsm'
      exv_names(ns+6) = 'tm'
      exv_names(ns+7) = 'localdtm'
      exv_names(ns+8) = 'pm'
      exv_names(ns+9) = 'krm'
      exv_names(ns+10) = 'uf'
      exv_names(ns+11) = 'vf'
      exv_names(ns+12) = 'wf'
      exv_names(ns+13) = 'pf'

      call pob_write_fvs_ex(pob,'stats.vtk',nexv,exv_names,exv_stats,ier,mask=mask)      

      ier = 0
      return
      end subroutine my_user_main_postpar

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

      END MODULE userdata_stats_routines

!  ***********************************************************************************
!     User initialization routine
!  ***********************************************************************************
      subroutine user_initialize(component,debug,ier)
      !!use us3d_user_module, only :user_startup,user_chem_reaction,user_bc_init,user_dataio_post, user_update_post
      use us3d_user_module, only : user_chem_reaction
      use us3d_user_module, only : user_bc_init
      use us3d_user_module, only : user_startup
      use us3d_user_module, only : user_main_postupdate
      use us3d_user_module, only : user_dataio_post
      use us3d_user_module, only : user_finalize
      use POST_MODULE, only : user_main_postpar 

!      use GASPROPS, only : user_transport
      use userdata      ! chemistry
      use userdata_stats    ! turns stats on and off
      use userbcdata    ! boundary conditions
      use geo2D         ! calculate area of polygon for flux matching
      use userdata_stats_routines ! stats integration to runtime
      
      implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier

      ier= 0
      user_chem_reaction => my_user_chem
      user_bc_init => my_user_bc_init
      user_startup => my_user_startup    
      user_main_postupdate =>  my_user_update_post 

      ! Trying this to collect stats
      if(istats.gt.0) then
            user_dataio_post => my_user_dataio 
      endif
      
      user_finalize => my_user_finalize
      user_main_postpar => my_user_main_postpar

      return
      end subroutine user_initialize
