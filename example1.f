#include "us3d_header.h"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Module used for data structures required by user-define calls
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      MODULE userdata

      implicit none
      save

      real(8) :: global_xmin,global_xmax

      contains

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This example of user_dataio illustrates how you can add storage of
c  your own variables into the US3D solver restart file.  This example
c  illustrates saving your variables on writing, and loading your
c  variables on restarting, including using the current run number.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
c Written by Heath B. Johnson                 Last modified: 2016-03-16
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "us3d_header.h"

      Subroutine my_user_dataio(cio,fname,irnum,ier)

      use MPI
      use HDF5
      use H5_EXTRAS
      use US3D_EXTRAS

      use us3d_global
      use us3d_varkinds

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
      integer :: hdf_err,ihave,irnumg
      integer(HID_T) :: rid
      logical :: new
      type(us3d_sfile) :: sfile

      real*8 :: testval
      real*8, pointer, dimension(:,:) :: dump

      integer :: epg,global_epg
      real*8, allocatable, dimension(:,:) :: myqdum

      integer :: have_mydata
      integer, parameter :: nuv=1  ! CHANGE NUMBER OF USER VARIABLES HERE
      character(LEN=200) :: rpath,dpath
      character(40) :: dname

      ! -- Set the name of your dataset here
      dname= 'mydata'

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
            if (ier/=0) exit
         endif

         ! --- Close run group and solution file
         call h5gclose_f(rid, hdf_err)
         call us3d_sfile_cnd(sfile,ier)

         exit
      enddo

      call us3d_check_int(icomw,MPI_MAX,ier)
      if (ier/=0) goto 999

!      call MPI_BCAST(new,1,MPI_LOGICAL,0,icomw,ier)
!      call MPI_BCAST(ihave,1,MPI_INTEGER,0,icomw,ier)
!      call MPI_BCAST(rpath,len(rpath),MPI_CHARACTER,0,icomw,ier)

      dpath= trim(rpath)//'/'//trim(dname)

      write(6,*) '-- dpath= "'//trim(dpath)//'"'

      ! --- Add any required IO here
      Select Case(cio)
      Case('w')

         ! -- Write cell-based data.  Assume sized mydata(my_nv,epg) on each processor
         epg= nel + neg                   ! This processor number of interior plus shared/ghost
         global_epg= ugrid%nc + ugrid%ng  ! number of interior + ghost cells in global grid
         allocate(myqdum(nuv,epg))

         ! -- Generate some fake data
         do i=1,epg
            myqdum(1,i) = dble(ugrid%ige(i))
         enddo

         ! -- Called with rpath, we expect to open and close the file in us3d_h5data_pw
         call write_us3d_qvals(dpath,myqdum,.false.,new,.false.,ier)
         if (ier/=0) goto 999

         deallocate(myqdum)
         if (id==0) write(6,*) '== Successfully wrote "'//trim(dname)//'"'

      Case('r')

         if (ihave==1) then

            epg= nel + neg                   ! This processor number of interior plus shared/ghost
            allocate(myqdum(nuv,epg))

            write(6,*) 'epg, size(qva)= ',epg,size(myqdum)

            dump => read_us3d_qvals(dpath,nuv,.false.,ier,qva=myqdum)
            if (ier/=0) goto 999

            ! -- Do something with the data you read.  For now we just test that we
            ! -- read it properly.
            do i=1,epg
               testval = myqdum(1,i) - dble(ugrid%ige(i))
               if (abs(testval).gt.1.0d-20) then
                  write(6,*) '*** Failed to properly load data: ',
     &               myqdum(1,i),dble(ugrid%ige(i)),abs(testval)
                  ier= 1
               endif
            enddo
            write(6,*) id,' read all data properly'

            deallocate(myqdum)
            if (id==0) write(6,*) '== Successfully read "'//trim(dname)//'"'

         else
            if (id==0) write(6,*) '-- No "'//trim(dname)//'" present for reading'
         endif

      end select

      return
 999  if (id==0) write(6,*) '*** Error in subroutine my_user_dataio'
      ier= 1
      return
      end subroutine my_user_dataio
      END MODULE userdata

c  ***********************************************************************************
c     User initialization routine
c  ***********************************************************************************
      Subroutine user_initialize(component,debug,ier)
      use us3d_user_module, only : user_dataio_post
      use userdata, only : my_user_dataio
      Implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier

      ier= 0

      user_dataio_post => my_user_dataio

      return
      end subroutine user_initialize

