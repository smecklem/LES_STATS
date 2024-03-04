!  ***********************************************************************************
!     Partition a grid for US3D into nparts.  You can use this routine to test other
!     partitioning methods.  The us3d_part_t object has access to the casefile_t
!     object used in preprocessing and all the variables contained in it.
!
!     Upon successful completion of this routine, the following should be defined:
!
!     up%nparts          - The number of partitions
!     up%part(1,1:my_nc) - The partition number for each cell (0:nparts-1)
!
!     Input
!     -----
!     up     - us3d_part_t derived type
!     nparts - Desired number of partitions
!     imeth  - Partitioning method
!
!     Output
!     ------
!     ier    - Returns nonzeo on an error
!
!  ***********************************************************************************
      subroutine user_partition(up,nparts,imeth,ier)
      use US3D_PART_TYPE
      !use US3D_PARTITIONING
      implicit none
      
      integer, intent(IN) :: nparts,imeth
      integer, intent(OUT) :: ier
   
      type(us3d_part_t) :: up
      
      integer :: id,nc,my_nc,nproc,istat
      
      ier  = 0
      nproc= up%nproc
      
      call MPI_COMM_RANK(up%iwcom,id,ier)
      if (ier/=0) goto 999
      
      ! -- User routines here -------------
      if (nparts==1) then

         if (id==0) write(6,*) '-- Only one partition requested'
         allocate(up%part(1,my_nc))

         up%nparts= 1
         up%part(1,:)= 0
         return
         
      else
      
         if (id==0) write(6,*) '*** No user partitioning routine defined!'
         goto 999

      endif
      
      return
 999  if (id==0) write(6,*) '*** Error in subroutine user_partition'
      ier= 1
      return
      end subroutine user_partition

!  ***********************************************************************************
!     Subroutine that is executed inside calcrate after exchange is called and before
!     the call to compute gradients.  To activate this routine, add into user_startup:
!
!     use switches, only : do_user_pregrad
!     do_user_pregrad= .true.
!
!     Input
!     -----
!     iform - Whether to store the RHS.  This is the variable iform that is passed
!             into calcrate.
!
!     Output
!     ------
!
!  ***********************************************************************************
      subroutine user_pregrad(iform)

      implicit none
      
      integer, intent(IN) :: iform
      
      return
      end subroutine user_pregrad
