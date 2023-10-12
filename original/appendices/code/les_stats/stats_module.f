cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     STATS MODULE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      MODULE stats
      implicit none
      save

      integer :: nmvar

      real(8) :: stat_time,unst_time
      real(8), allocatable,dimension(:,:) :: mean_var
      real(8), allocatable,dimension(:) :: delta

      integer, private :: olun= 6
      logical, private :: debug= .false.

      CONTAINS
c *****************************************************
      subroutine init_stats()
      use models, only: ns
      use mpivars, only : id
      Use sizing, only : nel
      implicit none

      nmvar = ns + 13

      if(id.eq.0) write(*,*) 'init stats ns = ', ns

      allocate(mean_var(nmvar,nel))
      allocate(delta(nel))
      call filterwidth()

      stat_time = 0.0d0
      unst_time = 0.0d0
      mean_var(:,:) = 0.0d0
      return
      end subroutine init_stats

c *****************************************************
c I made this with the best of intentions, but it segfaults????
c *****************************************************
      subroutine shutdown_stats(ier)
      use mpivars, only : id
      implicit none
      integer, intent(out) :: ier

      ier = -99
      
      !if (.not.allocated(mean_var)) then
      !    write(*,*) "id: ", id,"mean_var free fuckup @: ",loc(mean_var)
      !    stop
      !else
      !    deallocate(mean_var,stat=ier)
      !endif

      !if (.not.allocated(delta)) then
      !    write(*,*) "id: ", id,"delta free fuckup @: ",loc(delta)
      !    stop
      !else
      !     deallocate(delta, stat=ier)
      !endif

      return
      end subroutine shutdown_stats
 
c *****************************************************
      subroutine update_stats(istats)
      Use simvars
      use sizing, only : nel
      use mpivars, only : id
      Use turbulence, only : tcv1
      use  models, only : ns
      implicit none

      integer, intent(IN) :: istats

      integer :: i,n
      real(8) :: tmkgv, ksgs, eds, kr, Ret, Dat, tturb
      real(8) :: uf, vf, wf, k, nu, vel, pf, tkmgv
      
      stat_time = stat_time + dt
      if (istats.eq.2) unst_time = unst_time + dt

      do i = 1,nel
        do n = 1,ns
           mean_var(n,i)  = mean_var(n,i) + rhos(n,i)*dt ! Do we really need these?
        enddo

        call subscales(r(i), rmu(i), rnu(1,i), delta(i), tkmgv, ksgs, eds)

        vel = sqrt(u(i)**2 + v(i)**2 + w(i)**2)
        mean_var(ns+1,i)  = mean_var(ns+1,i) + u(i)*dt
        mean_var(ns+2,i)  = mean_var(ns+2,i) + v(i)*dt
        mean_var(ns+3,i)  = mean_var(ns+3,i) + w(i)*dt
        mean_var(ns+4,i)  = mean_var(ns+4,i) + vel/c(i)*dt
        mean_var(ns+5,i)  = mean_var(ns+5,i) + ksgs*dt 
        mean_var(ns+6,i)  = mean_var(ns+6,i) + t(i)*dt 
        mean_var(ns+7,i)  = mean_var(ns+7,i) + localdt(i)*dt
        mean_var(ns+8,i)  = mean_var(ns+8,i) + p(i)*dt

        if (istats.eq.2) then
            uf = u(i) - mean_var(ns+1,i)/stat_time
            vf = v(i) - mean_var(ns+2,i)/stat_time
            wf = w(i) - mean_var(ns+3,i)/stat_time
            pf = p(i) - mean_var(ns+8,i)/stat_time
            kr = 0.5*(uf**2 + vf**2 + wf**2)

            mean_var(ns+9,i)  = mean_var(ns+9,i) + kr*dt ! Resolved TKE
            mean_var(ns+10,i)  = mean_var(ns+10,i) + uf*dt
            mean_var(ns+11,i)  = mean_var(ns+11,i) + vf*dt
            mean_var(ns+12,i)  = mean_var(ns+12,i) + wf*dt
            mean_var(ns+13,i)  = mean_var(ns+13,i) + pf*dt
        endif

      enddo

      return
      end subroutine update_stats


c***********************************************************************
c **** Compute the turbulent timescales needed for TCI formulae
c      
c      Written by:    Nick Gibbons (CfH, 15.11.03 n.gibbons@uq.edu.au)
c      Modified  :    
c***********************************************************************

       Subroutine subscales(ra, rma, rnui, deltai, tkmgv, ksgs, eds)
       Use turbulence, only : tcv1
       Implicit none

       real(8), intent(in) :: ra, rma, rnui, deltai 
       real(8), intent(out) :: tkmgv, ksgs, eds

       real(8) :: ri,rmat,chi,chi3,fv1,nut,nu,rnat
       real(8) :: tsgs,ukmgv,usgs

       ! Calculate kinematic turbulent viscosity from SA formulas
       ri   = 1.0d0/ra                  ! inverse density
       !rma  = rmu(i)                   ! dynamic viscosity (mu)
       rnat = rnui                      ! nu_hat (SAvar, rmuti is nu_hat*r)
       rnat = max(rnat,0.0d0)           ! clip (not for SAneg itrb==2)
       rmat = rnat*ra                   ! nu_hat*r = mu_hat
       chi  = rmat/rma                  ! nu_hat*r/(nu*r)
       chi3 = chi*chi*chi
       fv1  = chi3/(chi3 + tcv1**3)
       nu   = rma*ri                    ! kinematic viscosity
       nut  = rmat*fv1*ri               ! kinematic turbulent viscosity

       ! Calculate subgrid tke from HIT approximation
       ksgs = nut*nut/(0.07d0*deltai)**2            ! subgrid tke
       eds = max(0.931*ksgs**(1.5d0)/deltai,1d-32)  ! subgrid dissipation
       tkmgv= sqrt(nu/eds)                          ! kolmogorov timescale
       !ukmgv= sqrt(sqrt(nu*eds))                   ! kolmogorov velocity
       !usgs = max(ukmgv,sqrt(2*ksgs/3.0d0))        ! subgrid velocity
       !tsgs = max(tkmgv,delta/usgs)                ! subgrid time scale
       !tmix = sqrt(tkmgv*tsgs)                     ! mixing time, or something

       end Subroutine subscales

c***********************************************************************
c **** Compute the LES filter width based on gehre_phd's method
c      
c      Written by:    Nick Gibbons (CfH, 15.11.03 n.gibbons@uq.edu.au)
c      Modified  :    
c***********************************************************************
      Subroutine filterwidth()
      Use geometry, only : xcn
      use connect, only : ien
      use sizing, only : nel
      Implicit none


      real(8) :: vec(3),maxlength
      integer :: ier,i,j,cell

      maxlength = -1.0d0

      do cell = 1,nel ! Loop over cells
          do i = 1,8 ! Loop over nodes in this cell
              do j = 1,8
                  vec(:) = xcn(:,ien(i,cell)) - xcn(:,ien(j,cell))
                  maxlength = max(maxlength, sqrt(vec(1)**2+vec(2)**2+vec(3)**2))
              enddo
          enddo
          if (maxlength.lt.1d-32) then
              write(*,*) "Filterwidth too small!", maxlength
          endif
          delta(cell) = maxlength ! Find delta in thismodule's scope
          maxlength = -1.0d0
      enddo

      end Subroutine filterwidth

c *****************************************************
c  Subroutine to co-ordinate IO for stats data
c  @curator: Nick Gibbons
c  Notes: I inherited this from Rolf, but I beleive it
c  to have originateed from someone at UoM
c *****************************************************
      subroutine stats_dataio(solf,cio,inum,h5oc,ier)
      use MPI
      use HDF5
      use us3d_logging
      use US3D_EXTRAS
      implicit none

      Integer, intent(IN) :: inum
      Integer, intent(OUT) :: ier
      character(LEN=*), intent(IN) :: solf
      character(LEN=1), intent(IN) :: cio
      logical, intent(IN) :: h5oc

      integer :: istat,hdf_err
      real(8) :: rtime

      ier= 0

      if (h5oc) then
         call h5open_f(hdf_err)
         if (hdf_err<0) stop
      endif

c
c *** Read or write the solution
c
      rtime = MPI_WTIME()

      Select Case(cio)

      Case('w')
         call us3d_wlog('== Writing stats HDF5 solution to: '//trim(solf),
     &      l=LOG_STATE,i=0,aei=0,e=olun)
         call stats_write(solf,inum,ier)
         if (ier/=0) goto 999

      Case('r')
         call us3d_wlog('== Reading stats HDF5 solution from: '//trim(solf),
     &      l=LOG_STATE,i=0,aei=0,e=olun)
         call stats_read(solf,inum,ier)
         if (ier/=0) goto 999

      End Select

      rtime = MPI_WTIME() - rtime
      if (us3d_cklog(i=0,l=LOG_TIMING)) then
         write(loglun,fmt='(1x,a,g15.5,1x,a)') '-- Data IO took ',rtime,'seconds'
      endif

      if (h5oc) then
         call h5close_f(hdf_err)
      endif

      call flush(olun)

      return
 999  if (id==0) write(olun,*) '*** Error in subroutine stats_dataio'
      ier= 1
      return
      end Subroutine stats_dataio
c *****************************************************
c     Routine to read the stats solution file in HDF5
c     format.  Note, only interior cells are read.
c
c     Input
c     -----
c     fname   - The HDF5 filename to read
c     inum    - (Optional) solution number.  Defaults to 1.
c
c     Output
c     ------
c     ier     - Returns nonzero on an error
c *****************************************************
      subroutine stats_read(fname,inum,ier)
      use MPI
      use HDF5
      use connect, only: ugrid
      use sizing, only : nel
      use us3d_extras
      use mpivars, only: icomw, id
      use switches, only : us3d_debug

      implicit none
      character(LEN=*) :: fname
      integer, intent(OUT) :: ier
      integer, intent(IN) :: inum

      integer :: n,e,s,hdf_err,istat,ihave_stats
      integer :: global_nbf,ier2
      integer :: old_nmvar,old_nsvar,iw
      integer(HID_T) :: rid
      logical :: new,attr_exists
      integer :: mpistat(MPI_STATUS_SIZE)
      integer, allocatable, dimension(:) :: ige

      type(us3d_sfile) :: sfile

      ier= 0

c     *** Each node stores their cell mapping
      allocate(ige(nel),STAT=istat)        ! Local to global cell mapping
      if (istat/=0) ier= 1

      call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 901

c     *** Set the solution number integer and string
      if (us3d_debug.and.id==0) write(olun,*) '-- Will read solution from run: ',inum
c'
c     *** Initializes file and groups
      new= .false.
      if (us3d_debug.and.id==0) write(olun,*) '-- Opening solution file'
      call us3d_sfile_init(fname,new,sfile,ier)
      if (ier/=0) goto 999

c     *** Open or create group for this run
      call us3d_sfile_orun(sfile,inum,rid,ier)
      if (ier/=0) goto 999

      if (sfile%nruns<inum) then
         if (id==0) then
            write(olun,*) '*** Unable to restart from run ',inum
            write(olun,*) '*** when solution file says it has ',sfile%nruns
         endif
         goto 999
      endif

      ihave_stats= 0
      call h5aexists_by_name_f(rid,'.','ihave_stats',attr_exists,hdf_err)
      if (attr_exists) then
         call h5ex_att_get(rid,'ihave_stats',ihave_stats,ier)
      endif

      if (ihave_stats/=1) then
         if (id==0) then
            write(olun,*) ' -- Run not listed as having a stats solution to read'
            write(olun,*) ' -- Will wait till next data_write to create h5 tree'
         endif
         return
      endif

c     *** Read solution attributes
      call h5ex_att_get(rid,'stats_time',stat_time,ier)
      call h5ex_att_get(rid,'unst_time',unst_time,ier)
      call h5ex_att_get(rid,'stats_mean_vars',old_nmvar,ier)
      !call h5ex_att_get(rid,'stats_stat_vars',old_nsvar,ier)

      if (old_nmvar/=nmvar ) then !.or. old_nsvar/=nsvar) then
         if (id==0) then
            write(olun,*) '*** Mismatch in number of variables'
            write(olun,*) '*** Old solution had nmvar= ',old_nmvar
            write(olun,*) '*** New solution has nmvar= ',nmvar
            !write(olun,*) '*** Old solution had nsvar= ',old_nsvar
            !write(olun,*) '*** New solution has nsvar= ',nsvar
            write(olun,*) '*** Unable to restart!'
         endif
         goto 999
      endif

c     *** Read solution data globally according to the global map
      ige(1:nel)= ugrid%ige(1:nel)

      if (us3d_debug.and.id==0) write(olun,*) 'Reading solution variables in parallel'
c'
      call us3d_data_r(rid,'stats-mean',ugrid%nc,nel,mean_var,ige,0,ier)
      if (ier/=0) goto 999
      !call us3d_data_r(rid,'stats-stat',ugrid%nc,nel,stat_var,ige,0,ier)
      !if (ier/=0) goto 999

      call MPI_BARRIER(icomw,ier2)

      if (us3d_debug.and.id==0) write(olun,*) '-- Finished reading stats'
c'
c     *** Check for any errors in above section
 802  call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 999

c     *** Close run group
      call h5gclose_f(rid, hdf_err)
      if (hdf_err<0) goto 999

c     *** Close file
      if (id==0.and.us3d_debug) write(olun,*) '-- Closing solution file'
      call us3d_sfile_cnd(sfile,ier)

c     *** Cleanup
      if (allocated(ige)) deallocate(ige)

      return
 901  if (id==0) write(olun,*) '*** Error allocating memory in subroutine stats_read'
c'
      call finalize(2)
 999  if (id==0) write(olun,*) '*** Error in subroutine stats_read'
      ier= 1
      return
      end subroutine stats_read
c *****************************************************
c     Routine to write the stats solution file in HDF5
c     format.  Note, only interior cells are written.
c
c     Input
c     -----
c     fname   - The HDF5 filename to write
c     inum    - (Optional) solution number.  Defaults to 1.
c
c     Output
c     ------
c     ier     - Returns nonzero on an error
c *****************************************************
      subroutine stats_write(fname,inum,ier)
      use MPI
      use HDF5
      use connect, only: ugrid
      use sizing, only : nel
      use us3d_extras
      use mpivars, only: icomw, id
      use switches, only : us3d_debug
      implicit none
      integer, intent(IN) :: inum
      character(LEN=*), intent(IN) :: fname
      integer, intent(OUT) :: ier

      integer :: i,n,e,s,is,hdf_err,nel_max,istat,iw,ie1,ie2,ict,ist,iemax
      integer :: ier2
      logical :: new,attr_exists,iexist
      integer :: ihave_stats,mpistat(MPI_STATUS_SIZE)
      integer, allocatable, dimension(:) :: ige
      integer(HID_T) :: rid
      real(8), allocatable, dimension(:,:) :: qdum!, rdum

      type(us3d_sfile) :: sfile

      ier= 0
      hdf_err= 0

c     *** Node zero initializes file and groups
      if (id==0) then

c        *** Set the solution number integer and string
         if (us3d_debug) write(olun,*) '-- Will write solution to run: ',inum

c        *** Decide whether to create a new solution file
         new= .true.
         inquire(file=fname, exist=iexist)
         if (iexist) new= .false.

         if (us3d_debug) write(olun,*) '-- Initializing solution file: "'//trim(fname)//'"'
         call us3d_sfile_init(fname,new,sfile,ier)
         if (ier/=0) goto 201

c        *** Open or create group for this run '
         call us3d_sfile_orun(sfile,inum,rid,ier)
         if (ier/=0) goto 201

c        *** Determine whether this run has a stats solution in it
         ihave_stats= 0
         new= .true.
         call h5aexists_by_name_f(rid,'.','ihave_stats',attr_exists,hdf_err)
         if (attr_exists) then
            call h5ex_att_get(rid,'ihave_stats',ihave_stats,ier)
         endif
         
         if (ihave_stats==1) then
            new= .false.
         else
            new= .true.
         endif
         
c        *** Store solution attributes
         if (us3d_debug) write(olun,*) '-- Storing attributes'

c        *** Attributes stored for this run
         call h5ex_att_add(rid,'ihave_stats',1,ier)
         call h5ex_att_add(rid,'stats_time',stat_time,ier)
         call h5ex_att_add(rid,'unst_time',unst_time,ier)
         call h5ex_att_add(rid,'stats_mean_vars',nmvar,ier)
         !call h5ex_att_add(rid,'stats_stat_vars',nsvar,ier)

      endif

 201  if (hdf_err<0) ier= 1
      call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 999

      !allocate(qdum(nmvar,nel),rdum(nsvar,nel),ige(nel),STAT=istat)
      allocate(qdum(nmvar,nel),ige(nel),STAT=istat)
      if (istat/=0) ier= 1
      call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 901

      do i= 1,nel
         ige(i)= ugrid%ige(i)

         qdum(1:nmvar,i) = mean_var(1:nmvar,i)
         !rdum(1:nsvar,i) = stat_var(1:nsvar,i)
      enddo

      call us3d_h5data_pw('stats-mean',ugrid%nc,qdum,ige,nel,.false.,icomw,new,ier,hid=rid)
      if (ier/=0) goto 999

      if (allocated(qdum)) deallocate(qdum)
      if (allocated(ige)) deallocate(ige)

      if (us3d_debug.and.id==0) write(olun,*) '-- Finished writing stats'
c'
c     *** Check for any errors in above section
 801  call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 999

c     *** Cleanup
      if (id==0) then

c        *** Close run group
         call h5gclose_f(rid, hdf_err)
         if (hdf_err<0) goto 999

c        *** Close file
         if (us3d_debug) write(olun,*) '-- Closing solution file'
         call us3d_sfile_cnd(sfile,ier)

      endif

      return
 901  if (id==0) write(olun,*) '*** Error allocating memory in subroutine stats_write'
c'
      call finalize(2)
 999  if (id==0) write(olun,*) '*** Error in subroutine stats_write'
      ier= 1
      return
      end subroutine stats_write
c *****************************************************
      function exv_stats(pob,np,ige,s,nexv,exv_names) result (stats)
      use POB_MODULE
      use POST_MODULE
      use mpivars, only : id
      implicit none

      type(pob_t), intent(IN) :: pob
      integer, intent(IN) :: np,nexv
      integer, dimension(np), intent(IN) :: ige
      real(8), dimension(pob%nsv,np), target, intent(IN) :: s
      character(10), dimension(nexv), intent(IN) :: exv_names
      real(8), dimension(nexv,np) :: stats

      integer :: ns,ier,n,i,j,cell
      real(8), dimension(:,:), allocatable :: sm,xcn
      integer, dimension(:,:), allocatable :: ien
      integer, dimension(:), allocatable   :: ign,iglob
      real(8) :: stat_time,unst_time,vec(3),maxlength

      ns= pob%ns
      write(*,*) "Beginning exv_stats"

      ! -- Allocate storage for reading mean and stats
      allocate(sm((ns+13),np))

      ! -- Read mean and stats
      call h5ex_att_get(pob%rid,'stats_time',stat_time,ier)
      call h5ex_att_get(pob%rid,'unst_time',unst_time,ier)
      if (ier/=0) stop
      !     write(*,*) 'Time of statistics taken is:',stat_time,'seconds'

      call us3d_data_r(pob%rid,'stats-mean',pob%nel,np,sm,ige,0,ier)
      if (ier/=0) stop

      !call us3d_data_r(pob%rid,'stats-stat',pob%nel,np,sp,ige,0,ier)
      !if (ier/=0) stop

      ! -- Store the mean
      stats(1:(ns+8),:)= sm(1:(ns+8),:)/stat_time
      if (unst_time.eq.0.d0) then
          stats((ns+9):(ns+13),:) = 0.0d0
      else
          stats((ns+9):(ns+13),:) = sm((ns+9):(ns+13),:)/unst_time
      endif

      deallocate(sm)

      !! -- Find Filter Widths !Maybe you can.
      !if (id/=0) then
      !    write(*,*) "Node co-ords are tricky in parralel."
      !    write(*,*) "Only run this routine with one processor!"
      !    stop
      !endif
      !allocate(ien(0:8,pob%nel))
      !allocate(iglob(pob%nel))
      !allocate(xcn(3,pob%nno))
      !allocate(ign(pob%nno))
      !do n=1,pob%nno
      !    ign(n) = n
      !enddo
      !do n=1,pob%nel
      !    iglob(n) = n
      !enddo
      !call pob_get_ien2(pob,iglob,ien,ier,pob%nel)
      !call pob_get_xcn2(pob,ign,xcn,ier,pob%nno)

      !maxlength = -1.0d0

      !do n = 1,np ! Loop over cells
      !    cell = ige(n) ! cell is the global cell index
      !    do i = 1,8 ! Loop over nodes in this cell
      !        do j = 1,8
      !            vec(:) = xcn(:,ien(i,cell)) - xcn(:,ien(j,cell))
      !            maxlength = max(maxlength, sqrt(vec(1)**2+vec(2)**2+vec(3)**2))
      !        enddo
      !    enddo
      !    stats(ns+12,n) = maxlength
      !    maxlength = -1.0d0
      !enddo

      !! -- Cleanup
      !deallocate(ien)
      !deallocate(xcn)
      !deallocate(ign)
      !deallocate(iglob)
      end function exv_stats

      END MODULE stats
