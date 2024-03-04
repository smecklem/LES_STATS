cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  User routines specific to MoVing Grids and Dynamics
c
c  Currently the dynamics routines will replace whatever
c  you do to the grid, if they are used.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ***********************************************************************************
c     Custom routine for MoVing the Grid
c
c  ***********************************************************************************
c  To get to this routine, set iMVG = 1 (or anything other than zero) 
c     (iMVG is in the module switches)
c
c  Use this routine to move the grid points to where evere you want them
c  Other routines will take care of recomputing the metrics, and adjusting the fluxes.
c  
c  Potentialy useful variables in module mv_grid
c
c  xcn_orig(1:3,1:nno)            -- original grid nodes
c  xce_b_orig(1:3,nel+1:net)      -- original ghost positions
c
c  By now, the time step dt should be computed.
c
c  ***********************************************************************************
      subroutine user_main_MVG()
        use sizing, only : nno
        use simvars, only : time,dt
        use geometry, only : xcn
        use mv_grid, only : xcn_orig

      implicit none
      Integer i
      Real*8   scale_time


c
c  --- Silly example, wobble the whole grid in the x-direction
c      by +/- 10%, at ~1000 Hz
c
c     scale_time = 2*3.14159/0.001
c
c     do i = 1,nno
c      xcn(1,i) = (1.0d0 + 0.1d0*DSIN(scale_time*(time+dt)))*xcn_orig(1,i)
c      xcn(2:3,i) = xcn_orig(2:3,i)
c     enddo
c 
c
c
c

      return
      end subroutine user_main_MVG

c  ***********************************************************************************
c  ***********************************************************************************
c     Custom routine for calculating force and moments in the Dynamics Routine
c  ***********************************************************************************
      subroutine user_dyn_forces(iuser)
        use mv_dyn
      implicit none
      Integer iuser

      iuser = 0  ! set to non-zero to skip regular routine

      return
      end subroutine user_dyn_forces
c  ***********************************************************************************
c     Custom routine for Linear part of the Dynamics Routine
c  ***********************************************************************************
      subroutine user_dyn_linear(iuser)
      use mv_dyn
      implicit none
      Integer iuser

      iuser = 0  ! set to non-zero to skip regular routine

      return
      end subroutine user_dyn_linear
c  ***********************************************************************************
c     Custom routine for Rotational part of the Dynamics Routine
c  ***********************************************************************************
      subroutine user_dyn_rotation(iuser)
      use mv_dyn
      implicit none
      Integer iuser

      iuser = 0  ! set to non-zero to skip regular routine

      return
      end subroutine user_dyn_rotation
c  ***********************************************************************************
c      Custom routine for rotating sections of the grid
c  ***********************************************************************************
      subroutine user_grid_rotation(iuser)
      use mv_dyn
      implicit none
      Integer iuser

      iuser = 0  ! set to non-zero to skip regular routine

      return
      end subroutine user_grid_rotation

c  ***********************************************************************************
c  ---------------------------------

c  ---------------------------------



