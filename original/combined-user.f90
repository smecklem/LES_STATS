!  ***********************************************************************************
! This file provides a means of allowing the users of US3D to compile
! their own routines and to access the data-structure that is already
! in place, as well as to create their own. In principle, this part of
! the code allows for the operation of an entirely different algorithm
! that is in sync with the solver of US3D and is MPI-aware.
!
! The default version of this file is used during compilation and it
! is compiled to a dynamically-loaded library that gets installed into
! the $US3D_HOME directory.  The user may build their own copy of this
! library by running the command us3d-build-user in the current working
! directory.  If a copy of the library is found in the working directory
! at runtime, it will be used instead of the default version.  You can
! verify this by running the command:
!
! us3d --show-libs
!
! If the need to use a file arises in your user subroutines, use the
! IO_FUNCS module and make a call to io_find_unused_lun(lun).
!
!     Common modules that you might want to use
!     -----------------------------------------
!     Use MPI
!     Use HDF5
!     Use IO_FUNCS
!     Use TXT_FUNCS
!     Use US3D_EXTRAS
!     Use POB_MODULE
!     Use POST_MODULE
!
!     US3D modules that you might want to use
!     ---------------------------------------
!     Use mpivars
!     Use switches
!     Use sizing
!     Use constants
!     Use models
!     Use connect
!     Use geometry
!     Use simvars
!     Use chemistry
!     Use turbulence
!     Use userdata
!
! Look in $US3D_HOME/include for interface files that you can include when building
! your user subroutines, and check the documentation and examples for usage tips.
!  ***********************************************************************************

#include "us3d_header.h"

!  ***********************************************************************************

!  ***********************************************************************************
!c***********************************************************************
!c **** Subroutine to compute Keq and derivative from Lewis fits
!c      Adapted from US3D v0.5.2mdb
!c      
!c      Written by:    Michael Barnhardt  (michael.d.barnhardt@nasa.gov)
!c      Created   :    09.02.11
!c      Modified  :    Graham Candler 07.02.12: changed lewis to CEA form
!c***********************************************************************
      Subroutine getkeq_lewis(jr,tt,rkeq,drk)

      Use GASPROPS,only: gas_ugcon
      Use models
      Use chemistry
      Use simvars

      Implicit none
      Real*8 tt,rkeq,drk

      Integer n,nn,jj,jr,idir
      Real*8 ts,ti,ti2,ti3,tln,t2,t3,t4,rr,thd
      Real*8 l1,l2,l3,dl2,dl3
      Real*8 a(9),exx,dex
      
      thd= 1.0d0/3.0d0

      ts  = tt                                  ! store T
      if(tt .gt. 2.0d4) tt = 2.0d4              ! limit T used in Lewis fits

      t2  = tt*tt
      t3  = t2*tt
      t4  = t3*tt
      ti  = 1.0d0/tt
      ti2 = ti*ti
      ti3 = ti2*ti
      tln = log(tt)

      if(tt .le. 1.0d3) then
        jj = 1
      else if(tt .gt. 1.0d3 .AND. tt .le. 6.0d3) then
        jj = 2
      else
        jj = 3
      endif

      exx = 0.0d0
      dex = 0.0d0
      do n = 1,ns
        if (.not.(gas%lcrs(jr,1,n).or.gas%lcrs(jr,2,n))) cycle
        a(:)  = gas%lewis(n,jj,:)

!        ! NOT CURRENTLY USED
!c       l1  = a(1)*ti2 + a(2)*ti + a(3)     +
!c    &        a(4)*tt  + a(5)*t2 + a(6)*t3  +
!c    &        a(7)*t4
        l2  =      -a(1)*ti2 +     a(2)*ti*tln +a(3)+ &
              0.5d0*a(4)*tt  + thd*a(5)*t2     + 0.25d0*a(6)*t3 +&
              0.2d0*a(7)*t4  +     a(8)*ti
        dl2 = 2.0d0*a(1)*ti3     + a(2)*ti2*(1.0d0 - tln) + 0.5d0*a(4)+&
              2.0d0*thd*a(5)*tt  + 0.75d0*a(6)*t2 + 0.8d0*a(7)*t3 - a(8)*ti2

        l3  = -0.5d0*a(1)*ti2 -       a(2)*ti +     a(3)*tln +&
                     a(4)*tt  + 0.5d0*a(5)*t2 + thd*a(6)*t3  +&
              0.25d0*a(7)*t4  +       a(9)
        dl3 = a(1)*ti3 + a(2)*ti2 + a(3)*ti + a(4) + a(5)*tt +&
              a(6)*t2  + a(7)*t3

        exx = exx + (gas%scrs(jr,1,n) + gas%scrs(jr,2,n))*( l2  - l3)
        dex = dex + (gas%scrs(jr,1,n) + gas%scrs(jr,2,n))*(dl2 - dl3)
      enddo

      rr   = (1.0d5*ti/gas_ugcon)**(-gas%snu(jr))
      exx  = max(-rkeq_ex_limiter,min(exx,rkeq_ex_limiter))
      rkeq = rr*exp(exx)
      drk  = rkeq*(gas%snu(jr)*ti + dex)

      tt   = ts                                           ! restore T

      return
      End Subroutine getkeq_lewis
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Module used for data structures required by userdata calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       MODULE geo2D

       implicit none
       save

      contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Sort a list of 2D points into a nice clockwise polygon
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: points  - List of points in 2D to be sorted
!              npoints - number of points in npoints
!
! **** Outputs: polygon - Clockwise arrangement of points
!      
! **** Notes: This routine also glues a copy of the first point to the
!      end of the shape. MAKE SURE len(polygon) can accomodate this.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine makepolygon(points,npoints,polygon)
       implicit none
 
       real(8), intent(in),dimension(:,:) :: points
       integer, intent(in) :: npoints
       real(8), intent(out),dimension(:,:):: polygon

       real(8) :: c(2), angle, tempangle, angles(24)
       integer :: i, worki, tempi, j, maxlocisweird(1)
       real(8),parameter :: pi=3.14159265359d0

       ! Compute Centroid
       c(:) = 0.0d0
       angles(:) = -1000.0
       do i= 1,npoints
           c(:) = c(:) + points(:,i)
       enddo
       c(:) = c(:)/npoints

       ! Calculate the angle from the centroid of each point
       do i = 1,npoints
           polygon(:,i) = points(:,i)-c(:)
           angle = atan2(polygon(2,i), polygon(1,i))
           if (angle .lt. 0.0) angle = 2*pi + angle
           angles(i) = angle
       enddo
       
       ! Explicit Sort because no std function ***grumbles***
       do i = 1, npoints
           maxlocisweird = maxloc(angles) ! What???! Why is this a one element array
           worki = maxlocisweird(1)
           polygon(:,i) = points(:,worki)
           angles(worki) = -1000.0
       enddo
           
       ! Add a copy of the last point for neat edges
       polygon(:,npoints+1) = polygon(:,1)

       end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Calculate the area of an arbitrary clockwise polygon
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: poly  - Points forming shape verticies (m,m)
!              npoly - Number of points forming verticies (3+)
!
! **** Ouputs: shapearea  - Area of Polygon (m2)
!      
! **** Notes: Assumes that shape has a copy of the start at the end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine shapearea(np, poly, area)
       implicit none
 
       real(8), intent(in),dimension(:,:) :: poly
       integer, intent(in) :: np
       real(8), intent(out) :: area
 
       real(8) :: c(2),x1,x2,y1,y2!,shapearea
       integer :: i
 
       c(:) = 0.0d0 ! Centroid
       area = 0.0d0
       do i = 1, np
           c(:) = c(:) + poly(:,i)
       enddo
       c(:) = c(:)/np
       
       ! Compute the area of each point-pair-centroid triangle
       do i=2,np+1
           x1 = poly(1,i-1) - c(1)
           y1 = poly(2,i-1) - c(2)
           x2 = poly(1,i) - c(1)
           y2 = poly(2,i) - c(2)
           area = area + 0.5*abs(x1*y2 - x2*y1)
       enddo

       end subroutine

       end module geo2D

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Module used for data structures required by user-define calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       MODULE userdata

       implicit none
       save

      contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Compute intersections of a plane with a hexahedral solid
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: hexcell  - Coordinates of each corner (m,m,m)
!              normal   - Hessian Normal form Normal (m,m,m)
!              d        - Hessian Normal form distance to origin (m)
!
! **** Outputs: interpoly - List of intersection points (m,m,m)
!               npoly     - Number of points found (0-6)
!      
! **** Notes: This routine embodies some assumptions about the order of            
!             points that aren't true for unstructured cells
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine hexintersect(hexcell, normal, d, interpoly, npoly)
       implicit none
    
       real(8), intent(in) :: hexcell(3,8), normal(3), d
       real(8), intent(out) :: interpoly(3,6)
       integer, intent(out) :: npoly
    
       integer :: i, srt(12), fin(12), j, idx
       real(8) :: e(3),denom,L,xp,yp,zp,xh,yh,zh,p(3)
       logical :: dupcheck(8)
    
       npoly = 0
    
       ! Correct nodes to walk the edges???
       srt = (/ 1,2,3,4,5,6,7,8,1,2,3,4 /)
       fin = (/ 2,3,4,1,6,7,8,5,5,6,7,8 /)
       dupcheck(:) = .false.
    
       do i = 1,12
           e = hexcell(:,fin(i)) - hexcell(:,srt(i))
    
           ! Check for parallelism
           denom = dot_product(normal,e)
           if (abs(denom).lt.1d-16) cycle 
    
           ! Compute intercept
           L = (d - dot_product(normal,hexcell(:,srt(i))))/denom
    
           if ((L.le.1.0d0).and.(L.ge.0.0d0)) then ! Easy Case
               if (L.eq.0.0d0) then
                   if (dupcheck(srt(i))) cycle 
                   p(:) = hexcell(:,srt(i))
                   dupcheck(srt(i)) = .true.
    
               else if (L.eq.1.0d0) then
                   if (dupcheck(fin(i))) cycle 
                   p(:) = hexcell(:,fin(i))
                   dupcheck(fin(i)) = .true.
               else
                   p(:) = hexcell(:,srt(i)) + L*e(:)
               endif
    
               npoly = npoly + 1
               interpoly(:,npoly) = p(:)
           endif
       enddo
       if (npoly.gt.6) then
           write(*,*) "ERROR detected!"
           write(*,*) "npoly .gt. 6" , npoly
       endif
    
       end subroutine
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Compute plane parameters from a cluster of points. 
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: plane  - Points making up boundary face (m,m,m)
!
! **** Outputs: normal  - Vector Normal to the plane
!               d       - Distance normal to plane from origin
!      
! **** Notes: Note that the fourth point isn't actually used.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine hessiannormal(plane, normal, d)
       implicit none
 
       real(8), intent(in) :: plane(3,4)
       real(8), intent(out) :: normal(3), d
 
       real(8) :: a(3), b(3), norm

       a(:) = plane(:,2) - plane(:,1)
       b(:) = plane(:,3) - plane(:,1)

       normal(1) = a(2)*b(3) - b(2)*a(3)
       normal(2) = b(1)*a(3) - a(1)*b(3)
       normal(3) = a(1)*b(2) - b(1)*a(2)

       norm = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
       if (norm==0.0) then
           write(*,*) "Nonreal normal vector detected"
           write(*,*) plane(:,1)
           write(*,*) plane(:,2)
           write(*,*) plane(:,3)
           write(*,*) plane(:,4)
           stop
       endif
       normal(:) = normal(:)/norm

       d = dot_product(normal(:),plane(:,1))
 
       end subroutine
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Top level control routine for finding intersection areas in 3D
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: plane  - Points making up boundary face (m,m,m)
!              normal - Normal Vector to boundary face (m,m,m)
!              poly   - Point cloud of 3D intersections (m,m,m)
!              npoly  - Number of points in poly        (3-6)
!
! **** Outputs: area  - 2D intersection area for cell and plane face (m2)
!      
! **** Notes: This is a horrible routine...
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine intersectarea(plane,normal,poly,npoly,area)
       use geo2D
       implicit none
 
       real(8), intent(in) :: plane(3,4), normal(3), poly(3,6)
       integer, intent(in) :: npoly
       real(8), intent(out):: area
 
       real(8) :: shape1(2,6), reshape1(2,7), shape2(2,6), reshape2(2,7)
       real(8) :: sixplane(3,6), ishape(2,24), ipoints(2,24), ireshape(2,25)
       integer :: i,ni,np,nf
       !logical :: isinside
 
       ! Flatten both planes into the same co-ordinate system
       sixplane(:,1:4) = plane(:,:) ! Needs same length as poly for later
       call flatten(sixplane, 4, poly, npoly, normal, shape1, shape2)

       ! Sort the points in the interpoly array into a polygon
       call makepolygon(shape2, npoly, reshape2)
       call makepolygon(shape1, 4    , reshape1)
 
       ! Compute inside points for plane
       ni = 0
       do i = 1,4
           if (isinside(reshape1(:,i), reshape2,npoly)) then
              ni = ni + 1
              ishape(:,ni) = reshape1(:,i)
           endif
       enddo

       ! Compute inside points for hexplane
       do i = 1,npoly
           if (isinside(reshape2(:,i), reshape1,4)) then
              ni = ni + 1
              ishape(:,ni) = reshape2(:,i)
           endif
       enddo

       ! Compute intersections
       call intersections(reshape1,4,reshape2,npoly,ipoints,np)
       nf = ni + np
       
       if (nf.gt.24) then
           write(*,*) "Too Many Points: segfault iminent"
           write(*,*) "Points found: ni", nf
       endif

       if (nf.eq.0) then
           area = 0.0
           return
       endif

       ! Sort point cloud into clockwise polygon
       ishape(:,ni+1:ni+np) = ipoints(:,1:np) ! Append
       call makepolygon(ishape,nf,ireshape)
       call shapearea(nf,ireshape,area)
 
       end subroutine
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Take two 3D point clouds in the same plane and rotate into 2D
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: shape1  - First Point Cloud (m,m,m)
!              n1      - Number of points in shape1 (Should be 4)
!              shape2  - Second point Cloud  (m,m,m
!              n2      - Number of Points in shape2 (3-6)
!
! **** Outputs: out1   - Flattened result if shape1 (m,m)
!               out2   - Flattened result of shape2 (m,m)
!      
! **** Notes: I called them shapes but they don't have to be sorted.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine flatten(shape1, n1, shape2, n2, normal, out1, out2)
       implicit none
 
       real(8), intent(in) :: shape1(3,6), shape2(3,6), normal(3)
       integer, intent(in) :: n1, n2
       real(8), intent(out) :: out1(2,6), out2(2,6)
 
       integer :: i
       real(8),parameter :: pi=3.14159265359d0
       real(8) :: angle, a(3,4), b(3,6), t(3,4), p(3,6),O(3), c, s
 
 
       O(:) = shape1(:,1) ! Pick an arbitrary point for the origin 
       do i =1,n1
           a(:,i) = shape1(:,i) - O(:)
       enddo

       do i =1,n2
           b(:,i) = shape2(:,i) - O(:)
       enddo

       angle = atan2(normal(2), normal(1))
       c = cos(-angle)
       s = sin(-angle)

       ! Rotate about Z direction so that normal is in z-x plane
       do i = 1,n1
           t(1,i) = c*a(1,i) - s*a(2,i)
           t(2,i) = s*a(1,i) + c*a(2,i)
           t(3,i) =   a(3,i) 
       enddo
       do i = 1,n2
           p(1,i) = c*b(1,i) - s*b(2,i)
           p(2,i) = s*b(1,i) + c*b(2,i)
           p(3,i) =   b(3,i) 
       enddo

       ! Rotate about y direction so that normal is along the z axis
       angle = atan2(normal(3), sqrt(normal(1)**2+normal(2)**2)) - pi/2.0
       c = cos(angle)
       s = sin(angle)
       
       do i = 1,n1
           out1(1,i) = c*t(1,i) + s*t(3,i)
           out1(2,i) =   t(2,i)
           !out1(3,i) =-s*t(1,i) + c*t(3,i)
       enddo
       do i = 1,n2
           out2(1,i) = c*p(1,i) + s*p(3,i)
           out2(2,i) =   p(2,i)
           !out2(3,i) =-s*p(1,i) + c*p(3,i)
       enddo
 
       end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Check if point is inside an arbitrary polygon. (2D)
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: point  - Check if this point is inside (m,m)
!              poly   - Collection of clockwise points (m,m)
!              npoly  - Number of points in poly (3-6)
!
! **** Ouputs: Logical Flag if point is definitely inside
!      
! **** Notes: Assumes that each shape has a copy of the start at the end
!      This routine deliberately returns false for cases where the point
!      is on the edge of the shape. They are caught by intersections.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       function isinside(point, poly, npoly)
       implicit none
 
       real(8), intent(in) :: point(2), poly(2,7)
       integer, intent(in) :: npoly
       logical :: isinside
 
       integer :: i
       real(8) :: angle, edge(2), conn(2), edge_angle, conn_angle
       real(8),parameter :: pi=3.14159265359d0

       ! Assuming that the points are arranged rotationally
       do i = 2,npoly+1
         edge(:) = poly(:,i)-poly(:,i-1)
         conn(:) = point(:) - poly(:,i-1)
         edge_angle = atan2(edge(2), edge(1))
         conn_angle = atan2(conn(2), conn(1))
         angle = conn_angle - edge_angle
         
         if (angle.lt.0.0)  angle = 2*pi + angle
         if (angle.gt.2*pi) angle = angle - 2*pi
         if (angle.lt.pi) then
             isinside = .false.
             return
         endif
       enddo
 
       isinside = .true.
       return
       end function
       
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Compute the points formed by intersecting the lines of two shapes
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: shape1  - First shape (m,m)
!              n1      - Number of points in first shape (Should be 4)
!              shape2  - Second shape (m,m)
!              n2      - Number of points in second shape (3-6)
!
! **** Outputs: interpoints  - Point Cloud of intersecting points
!               npoints      - Numbe of intersections found
!      
! **** Notes: The returned points are not sorted. (call makepolygon)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine intersections(shape1,n1,shape2,n2,interpoints,npoints)
       implicit none
 
       real(8), intent(in) :: shape1(2,7), shape2(2,7)
       integer, intent(in) :: n1,n2
       real(8), intent(out):: interpoints(2,24)
       integer, intent(out):: npoints
  
       integer :: i,j
       real(8) :: x1,x2,x3,x4,y1,y2,y3,y4,d1,d2,m1,m2,m3,xi,yi

       npoints = 0
       do i = 2,n1+1        ! Weird Loop because of the extra point at the end
           do j = 2,n2+1
               x1 = shape1(1,i)
               y1 = shape1(2,i)
               x2 = shape1(1,i-1)
               y2 = shape1(2,i-1)
               x3 = shape2(1,j)
               y3 = shape2(2,j)
               x4 = shape2(1,j-1)
               y4 = shape2(2,j-1)

               ! Do duplicate testing
               if ((x2==x4).and.(y2==y4)) then 
                   interpoints(1,npoints+1) = x2
                   interpoints(2,npoints+1) = y2
                   npoints = npoints + 1
                   cycle
               endif
               if ((x1==x3).and.(y1==y3)) cycle 
               if ((x1==x4).and.(y1==y4)) cycle
               if ((x2==x3).and.(y2==y3)) cycle
               
               ! Calculate Intersection of Infinite lines
               m1 = x1*y2 - y1*x2
               m2 = x3*y4 - y3*x4
               m3 = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
               if (m3==0) cycle  ! Parallel Lines
               xi = (m1*(x3-x4) - (x1-x2)*m2)/m3
               yi = (m1*(y3-y4) - (y1-y2)*m2)/m3
               
               ! Check to see if in bounds
               d1 = sqrt((x2-x1)**2 + (y2-y1)**2)
               d2 = sqrt((x4-x3)**2 + (y4-y3)**2)
               m1 = sqrt((xi-(x2+x1)/2.0)**2 + (yi-(y2+y1)/2.0)**2)
               m2 = sqrt((xi-(x4+x3)/2.0)**2 + (yi-(y4+y3)/2.0)**2)
               if ((m1>d1/2.0).or.(m2>d2/2.0)) cycle
               if (npoints.eq.23) then
                   write(*,*) "Too manyintersecting points!"
                   stop
               endif
               interpoints(1,npoints+1) = xi
               interpoints(2,npoints+1) = yi
               npoints = npoints + 1
           enddo
       enddo
 
       end subroutine
 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Calculate the flux vector normal to a surface
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: rhos  - Partial Densities of Other Cell (kg/m3)
!              vel   - Velocity vector of Other Cell   (m/s,m/s,m/s)
!              t     - Temperature inside Other Cell   (K)
!              normal- Normal Vector of surface        (m,m,m)
!
! **** Ouputs: fluxnormal - Vector of ne length containing fluxes
!      
! **** Notes: The rotation that happens here is needed because later we
!      want to reverse the process: i.e. Calculating some primitives to
!      match a flux vector. This is very hard in 3D, but easy if the 
!      flux vector is aligned with the surface. (Probably)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       function fluxnormal(rhos,vel,t,tv,trb,normal)
       use models, only : ns,gas
       use switches, only: ivib,itrb
       use gasprops
       implicit none
 
       real(8), intent(in) :: rhos(ns), vel(3), t, tv, trb(2), normal(3)
       real(8) :: fluxnormal(ns+7),v1(3),u,v,w,angle,c,s
       real(8) :: p,g,a,hf,cv,r,re,ke,ev
 
       ! Rotate about Z direction so that normal is in z-x plane
       angle = atan2(normal(2), normal(1))
       c = cos(-angle)
       s = sin(-angle)

       v1(1) = c*vel(1) - s*vel(2)
       v1(2) = s*vel(1) + c*vel(2)
       v1(3) =   vel(3) 

       ! Rotate about y direction so that normal is along the x axis
       angle = atan2(normal(3), sqrt(normal(1)**2+normal(2)**2))
       c = cos(angle) ! Really subtle issue here. The 'angle' above
       s = sin(angle) ! is in a weird lefthanded reference frame.
       
       u = c*v1(1) + s*v1(3)
       v =   v1(2)
       w =-s*v1(1) + c*v1(3)

       ! Reconstruct Variables as needed
       call gas_therm1(gas,rhos,t,p,g,a,hf,cv)
       r = sum(rhos)
       select case(ivib)
       case(0)
           ev=0.0d0
       case(1)
           call gas_vibe1(gas, rhos, tv, ev)
       case(2)
           call gas_vibel1(gas,rhos, tv, ev)
       case(3)
           call gas_vibel1(gas,rhos, tv, ev)
           ! Note gas therm includes ev in it's cv calculation
           cv = dot_product(rhos,gas%cvs)/r ! Transrot only
       end select

       ke = 0.5d0*(u**2 + v**2 + w**2) 
       re = r*(hf + ke + cv*t + ev) + trb(2) ! Include TKE
       
       ! Construct the normal Flux
       fluxnormal(1:ns) = rhos(:)*u
       fluxnormal(ns+1) = r*u**2 + p
       fluxnormal(ns+2) = r*u*v
       fluxnormal(ns+3) = r*u*w
       fluxnormal(ns+4) = u*(re + p)
       fluxnormal(ns+5) = r*u*ev
       fluxnormal(ns+6) = u*trb(1)
       fluxnormal(ns+7) = u*trb(2)

       end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Take a set of fluxes and return the primitive variables 
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: flux  - Flux vector (rhos..., momx,momy,momz,E)
!
! **** Ouputs: rhos  - Partial densities vector      (kg/m3...)
!              vel   - Velocity vector in grid frame (m/s,m/s,m/s)
!              ret   - Total energy                  (J/kg)
!      
! **** Notes: Put in complex thermo/turbulence stuff later
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine inverseflux(flux,normal,rhos,vel,ret,t,tv,trb,rav)
       use models, only : ns, gas
       use switches, only: ivib,itrb
       implicit none
 
       real(8), intent(in) :: flux(ns+7), normal(3),rav
       real(8), intent(out) :: rhos(ns), vel(3), ret, t, tv, trb(2)

       real(8) :: knum,kdem,hform,fns,k,a,b,c,u1,u2,u,s,angle
       real(8) :: cv,v,w,p,ev,tke,r1,r2,rsve(ns,1),tve(1),eve(1)
       integer :: i

       knum = 0.0
       kdem = 0.0
       hform= 0.0
       Fns = 0.0
       cv  = 0.0

       do i = 1,ns
           knum = knum + flux(i)*gas%gsp(i)
           kdem = kdem + flux(i)*gas%cvs(i)
           hform= hform + flux(i)*gas%hs(i)
           Fns = Fns + flux(i)
           cv = cv + flux(i)*gas%cvs(i)
       enddo

       k = knum/kdem + 1.0d0
       cv = cv/Fns

       ! Energy sum equation is quadratic in u
       a = Fns*(0.5d0-(k/(k-1)))
       b = flux(ns+1)*k/(k-1)
       c = (flux(ns+2)**2 + flux(ns+3)**2)/(2*Fns) - flux(ns+4) + hform
       c = c + flux(ns+5) + flux(ns+7) ! Add vibe and TKE (might be zero)

       u1= (-b + sqrt(b**2 - 4*a*c))/2.0d0/a 
       u2= (-b - sqrt(b**2 - 4*a*c))/2.0d0/a
       
       ! Choose the velocity that is gives the closest density match
       r1 = Fns/u1
       r2 = Fns/u2

       if (abs(r1-rav).lt.abs(r2-rav)) then
           u = u1
       else
           u = u2
       endif

       ! Calculate primitives
       rhos(:) = flux(1:ns)/u
       p = flux(ns+1) - Fns*u
       ret = flux(ns+4)/u - p
       ev = flux(ns+5)/Fns
       t = p/(sum(rhos)*cv*(k-1))
       trb(1) = flux(ns+6)/u
       trb(2) = flux(ns+7)/u
       select case(ivib)
       case(0)
           tv = t
       case(1)
           tv = gas_tv(rhos,ev,t)
       case(2) ! uncommented next line and commented remainder for setonix int.
           tv = gas_tv(rhos,ev,t)
           !rsve(:,1) = rhos(:)
           !eve(1) = ev
           !call gas_tvel2(gas, 1, rsve, eve, tve)
           !tv = tve(1)
       case(3)
           tv = t !gas_tvel(rhos,ev,t) These can be slightly different??
       end select

       vel(1) = u
       vel(2) = flux(ns+2)/Fns
       vel(3) = flux(ns+3)/Fns

       ! Rotate about y direction so that normal is not along the x axis
       angle = atan2(normal(3), sqrt(normal(1)**2+normal(2)**2))
       c = cos(-angle) ! Weird axes thing here again
       s = sin(-angle)
       u = c*vel(1) + s*vel(3)
       v =   vel(2)
       w =-s*vel(1) + c*vel(3)

       ! Rotate about Z direction so that normal is no longer in z-x plane
       angle = atan2(normal(2), normal(1))
       c = cos(angle)
       s = sin(angle)
       vel(1) = c*u - s*v
       vel(2) = s*u + c*v
       vel(3) =   w

       end subroutine
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Compute the maximum, minimum extent of this face
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: j1,j2  - Start and ending face indecies (integer)
!
! **** Ouputs: x1    - Lowest X,Y,Z of any node      (m,m,m)
!              x2    - Highest X,Y,Z of any node      (m,m,m)
!      
! **** Notes: This is needed for optimisation purposes only
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       subroutine bbx(j1, j2, x1, x2) 
       use connect, only : ifn
       use geometry, only : xcn
       implicit none
 
       integer, intent(in) :: j1, j2
       real(8), intent(out)  :: x1(3), x2(3)

       integer :: j,jj,n
       real(8) :: maxx,minx,maxy,miny,maxz,minz,point(3)

       maxx = -1d16
       maxy = -1d16
       maxz = -1d16
       minx =  1d16
       miny =  1d16
       minz =  1d16

       do j=j1,j2
           do jj=1,4
               n = ifn(jj,j)
               point(:) = xcn(:,n)
               maxx = max(maxx,point(1))
               minx = min(minx,point(1))
               maxy = max(maxy,point(2))
               miny = min(miny,point(2))
               maxz = max(maxz,point(3))
               minz = min(minz,point(3))
           enddo
       enddo

       x1(1) = minx
       x1(2) = miny
       x1(3) = minz
       x2(1) = maxx
       x2(2) = maxy
       x2(3) = maxz

       end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Find the SHO temperature from the vibe energy
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: gass  - Gas Object 
!              rhos  - Species partial densities       (kg/m3,...)
!              t     - Transrotational Temperature     (K)
!              ev    - Vibrational Energy              (J/m3)
!
! **** Ouputs: gas_tv - Vibrational Temperature        (K)
!      
! **** Notes: Code copied from time_integration.f. Uses t as a guess.
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       function gas_tv(rhos, ev, t)
       use models, only : ns, gas
       !use gasprops, only gas_t
       implicit none

       real(8), intent(IN) :: rhos(ns),t,ev

       integer :: it,m,n
       real(8) :: gas_tv,rt,ri,q1,q2,t1,t2,t3,t4

       rt = t ! This may require some damping

       do it = 1,20 
          ri = 1.0d0/rt
          q1 = sum(rhos)*ev
          q2 = 0.0d0
          do m = 1,gas%nvm
             n  = gas%idvs(m)
             t1 = exp(gas%thetv(m)*ri)
             t2 = gas%thetv(m)/(t1 - 1.0d0)
             t3 = rhos(n)*gas%gsp(n)
             t4 = t2*ri
             q1 = q1 - t3*t2             ! Ev - Ev(k)
             q2 = q2 + t3*t1*t4*t4       ! dEv/dTv(k)
          enddo
          if (abs(q1/q2)/abs(rt).lt.1d-6) exit
          rt = rt + q1/q2
          if (it.eq.20) then
              write(*,*) "Convergence of Tv failed in inverseflux"
              write(*,*) "Temp: ", rt
              stop
          endif
       enddo

       gas_tv = rt + q1/q2


       end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Find the CEA vibe-electronic temperature from the energy
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: gass  - Gas Object 
!              rhos  - Species partial densities       (kg/m3,...)
!              t     - Transrotational Temperature     (K)
!              ev    - Vibrational Energy              (J/m3)
!
! **** Ouputs: gas_tvel - Vibronic Temperature        (K)
!      
! **** Notes: Code copied from gasprops with additions for robustness
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!       function gas_tvel(rhos, ev, t)
!       use models, only : ns, gas
!       use gas_constants, only : gas_ugcon
!       implicit none
!
!       real(8), intent(IN) :: rhos(ns),t,ev
!       integer :: i,jj,n,it
!       real(8) :: q1,q2,tt,t1,t2,t3,t4,t5,ti,ti2,tln,a(9),gas_tvel
!       real(8) :: evels,cvels
!
!       tt = t
! 
!       do it = 1,20
!          tt = max(tt,1.0d0) ! Clip to prevent NaN during convergence
!          jj  = 2
!          if (tt.le.1.0d03) jj = 1
!          if (tt.gt.6.0d03) jj = 3
!          q1  = sum(rhos)*ev
!          q2  = 0.0d0
!          t2  = tt*tt
!          t3  = t2*tt
!          t4  = t3*tt
!          t5  = t4*tt
!          ti  = 1.0d0/tt
!          ti2 = ti*ti
!          tln = log(tt)
! 
!          do n = 1,gas%ns
!              a(:)  = gas%lewis(n,jj,:)
!              evels =(- a(1)*ti + a(2)*tln + a(3)*tt + 0.5d0*a(4)*t2 + a(5)*t3/3.0d0
!     .               + 0.25d0*a(6)*t4 + 0.2d0*a(7)*t5 + a(8))*gas_ugcon*gas%smwi(n)
!     .               - gas%cps(n)*tt + gas%dh298(n,jj) - gas%hs(n)
!              cvels =(  a(1)*ti2 + a(2)*ti + a(3) + a(4)*tt + a(5)*t2
!     .               + a(6)*t3 + a(7)*t4)*gas_ugcon*gas%smwi(n) - gas%cps(n)
!              q1    = q1 - rhos(n)*evels
!              q2    = q2 + rhos(n)*cvels
!          enddo
! 
!          if (abs(q1/q2)/abs(tt).lt.1d-6) exit
!          tt = tt + q1/q2 ! possibly consider damping this correction
!          if (it.eq.20) then
!              write(*,*) "Convergence of Tv-el failed in inverseflux"
!              write(*,*) "Temp: ", tt
!              tt = 0.0d0
!              q1 = -1.0d0
!              q2 = 1.0d0
!          endif
!       enddo
! 
!       gas_tvel = tt + q1/q2
!
!       end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Check that the old solution file matches the new one
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: gass  - Gas Object 
!              rhos  - Species partial densities       (kg/m3,...)
!              t     - Transrotational Temperature     (K)
!              ev    - Vibrational Energy              (J/m3)
!
! **** Ouputs: gas_tvel - Vibronic Temperature        (K)
!      
! **** Notes: Code copied from gasprops with additions for robustness
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       function check_consistent(othersol)
       use models, only : ns
       use sizing, only : nv,nt
       use US3D_EXTRAS, only : us3d_sfile
       use h5_extras, only : h5ex_att_get
       implicit none

       type(us3d_sfile),intent(IN)  :: othersol
       logical :: check_consistent
       integer :: nto,nvo,nso,ier

       check_consistent = .true.
       nso = -999
       nto = -999
       nvo = -999

       call h5ex_att_get(othersol%siid,'ns',nso,ier)
       if (ns.ne.nso) then 
           write(*,*) "Number of species does not match new grid."
           check_consistent=.false.
           return
       endif

       call h5ex_att_get(othersol%siid,'nveq',nvo,ier)
       if (nv.ne.nvo) then
          write(*,*) "Number of vibration terms does not match new grid"
          check_consistent=.false.
          return
       endif

       !call h5ex_att_get(othersol%siid,'nteq',nto,ier)
       !if (nt.ne.nto) then 
       !  write(*,*) "Number of turbulence terms does not match new grid"
       !  check_consistent=.false.
       !  return
       !endif

       return
       end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! **** Check that a cell is inside or crossing the bbx
! **** @author: Nick Gibbons, CfH 15/01/28 (n.gibbons@uq.edu.au)
!
! **** Inputs: - x1  Lower x,y,z of bounding box
!              - x2  Upper x,y,z of bounding box
!              - xcc  Nodes of this cell
!
! **** Ouputs: - logical flag is the cell overlaps the 
!      
! **** Notes: 
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       function cell_bbx_overlap(x1,x2,xcc)
       implicit none

       real(8),intent(IN)  :: x1(3), x2(3), xcc(3,8)
       logical :: cell_bbx_overlap
       integer :: i
       real(8) :: maxi,mini

       do i=1,3
           ! test if all Xi is completely outside
           maxi = maxval(xcc(i,:))
           mini = minval(xcc(i,:))
           if ((maxi.lt.x1(i)).or.(mini.gt.x2(i))) then
               cell_bbx_overlap= .false.
               return
           endif
       enddo

       ! If we reached here it means that the cloud of points must
       ! overlap the bbx in all three dimensions
       cell_bbx_overlap = .true.
       return
       end function
      subroutine my_user_chem(ns,nq,i,ir,irm,irxt,xs,ws,aa,rate_data)
!      use us3d_varkinds
!     Use mpivars
     Use switches
     Use sizing, only :  nt,nv,ne
     Use constants
     Use models, only : smw,hs,smwi,gas,cvs
     !Use us3dh_chem_reaction
!     Use connect
     Use geometry, only : voli
     Use simvars, only : u,v,w,t,p,cv,gcon,r
     Use chemistry
!     Use turbulence
     
      integer, intent(IN)  :: ns,nq,i,ir,irm,irxt
      real*8 tab,tabi,tabl,taf,tafi,tafl,drk,rkeq,dkfm,dkfu
      real*8 atroe,t1,t2,t3,Fcent,f1,cc,nn,ftp,tlog
      real*8 rs(ns),cvv(ns),ti,tt,drt,rkinf,gg,cfm_inf,ea_0,ea_inf
      real*8 rkfs,rkbs,rkfm,rkbm,rkb0,rkf0,rkbinf,tmp1,ifac,mfac,bfac
      real*8 rkfs0,rkbs0,rkfm0,rkbm0,rate,pr,prb
      Real*8 ccv,pp,vi,uu,vv,ww,absu,revst,rcvi,rcvv,rcvt,xf,xb
      Real*8 athetu,athetu0,cfmu,cfmu0,aetau,aetau0,rkf,rkb,tbeu(30)
      real*8 dtdrs(nq),cfm_0(ns),dtf(nq)
      real(rk), intent(IN)  :: xs(ns)
      real(rk), intent(INOUT) :: ws(nq),aa(nq,nq)
      real(rk), intent(OUT), dimension(:,:) :: rate_data

      integer m,zz,n,nruu !nru = number of reactions user
      integer i1,i2,i3,j1,j2,j3
   
! *** translational temperature and TTv derivatives wrt U
!
! add in after reaction rates 
      ! have to consider the new big U .... 
      !cvs = gas%cvs
      gg  = gcon(i)
      ccv = cv(i)
      uu  = u(i)
      vv  = v(i)
      ww  = w(i)
      pp  = p(i)
      !write(*,*) 'pp,gg,ccv,uu,vv,ww'
      !write(*,*) pp,gg,ccv,uu,vv,ww
      vi  = 1.0d0/voli(i)
      !rs(:)  = rhos(:,i)
      absu = 0.5d0*(uu*uu + vv*vv + ww*ww)
      !write(*,*) 'absu'
      !write(*,*) absu

      rcvi = 1.0d0/(r(i)*cv(i))

      do n = 1,ns
        dtdrs(n) = (absu - hs(n) - cvs(n)*tt)*rcvi
      enddo
      dtdrs(ns+1) = -uu*rcvi
      dtdrs(ns+2) = -vv*rcvi
      dtdrs(ns+3) = -ww*rcvi
      dtdrs(ns+4) =     rcvi 
      do n=1,nv
        dtdrs(ns+4+n) =  -rcvi 
      enddo
      do n=1,nt
        dtdrs(ns+4+nv+n) = 0.0d0
      enddo
      !write(*,*) dtdrs

      !write(*,*) dtdrs

      tt=t(i)
      ti = 1.0d0/tt
      tlog = log(tt)

      !write(*,*) dtdu 
      dtf(:) = dtdrs(:)

      !call getkeq_lewis(ir,t(i),rkeq,drk)
      !write(*,*) icsd(ir,:)
      !IF ( ir == 39 ) THEN

      !write(*,*) 'solving reaction 39'
      ! A + M  <=> B + C + M
      i1 = 6  !icsd(ir,1)
      j1 = 7  !icsd(ir,2)
      j2 = 18 !icsd(ir,3)
      tafi = ti
      tafl = tlog
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
      do m = 1,ns
        rkfm = cfm_0(m)*rkf
        t1 = ifac*(-rkfm*xf)*smwi(m)
        aa(i1,m) = aa(i1,m) + t1
        aa(j1,m) = aa(j1,m) - t1
        aa(j2,m) = aa(j2,m) - t1

      enddo
      
      mfac = rkfs**2/(rkinf+rkfs)**2
      bfac = rkinf*rkfs/(rkinf+rkfs)
      rate = -bfac*xf

      !write(*,*) xf, xb
      !write(*,*) 'ifac, bfac, mfac, rate'
      !write(*,*) ifac, bfac, mfac, rate
      
      ws(i1) = ws(i1) + rate
      ws(j1) = ws(j1) - rate
      ws(j2) = ws(j2) - rate
      
      !dkfu = ifac*rkfs*tafi*(aeta(n) + tafi*athet(n))*xf
      !dkfm = mfac*rkinf*tafi*(aeta(n) + tafi*gas%txb(n))*xf
      dkfu = ifac*rkfs*tafi*(tafi*ea_0)*xf
      dkfm = mfac*rkinf*tafi*(tafi*ea_inf)*xf

      
      do m = 1,nq
         t1 = -(dkfu+dkfm)*dtf(m)
         aa(i1,m) = aa(i1,m) + t1
         aa(j1,m) = aa(j1,m) - t1
         aa(j2,m) = aa(j2,m) - t1
      enddo
      
      !write(*,*) ' 1   drt, aa(i1,i1), smwi(i1)'
      !write(*,*) drt, aa(i1,i1), smwi(i1)
      drt = bfac*smwi(i1)
      aa(i1,i1) = aa(i1,i1) - drt
      aa(j1,i1) = aa(j1,i1) + drt
      aa(j2,i1) = aa(j2,i1) + drt
      
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
      tafl = tlog
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
      do m = 1,ns
        !write(*,*) i1,i2,j1,m
        rkfm = cfm_0(m)*rkf
        t1 = ifac*(-rkfm*xf)*smwi(m)
        aa(i1,m) = aa(i1,m) + t1
        aa(i2,m) = aa(i2,m) + t1
        aa(j1,m) = aa(j1,m) - t1

      enddo
      
      mfac = rkfs**2/(rkinf+rkfs)**2
      bfac = rkinf*rkfs/(rkinf+rkfs)
      rate = -bfac*xf
      
      !write(*,*) xf, xb
      !write(*,*) 'ifac, bfac, mfac, rate'
      !write(*,*) ifac, bfac, mfac, rate
      ws(i1) = ws(i1) + rate
      ws(i2) = ws(i2) + rate
      ws(j1) = ws(j1) - rate
      
      !dkfu = ifac*rkfs*tafi*(aeta(n) + tafi*athet(n))*xf
      !dkfm = mfac*rkinf*tafi*(aeta(n) + tafi*gas%txb(n))*xf

      dkfu = ifac*rkfs*tafi*(tafi*ea_0)*xf
      dkfm = mfac*rkinf*tafi*(tafi*ea_inf)*xf

      !write(*,*) 'dkfu,dfkm'
      !write(*,*) dkfu,dkfm
      
      do m = 1,nq
         t1 = -(dkfu+dkfm)*dtf(m)
         aa(i1,m) = aa(i1,m) + t1
         aa(i2,m) = aa(i2,m) + t1
         aa(j1,m) = aa(j1,m) - t1
      enddo
      
      !write(*,*) ' 2   drt, aa(i1,i2), smwi(i2)'
      !write(*,*) drt, aa(i1,i2), smwi(i2)
      drt = bfac*smwi(i1)*xs(i2)
      aa(i1,i1) = aa(i1,i1) - drt
      aa(i2,i1) = aa(i2,i1) - drt
      aa(j1,i1) = aa(j1,i1) + drt

      drt = bfac*smwi(i2)*xs(i1)
      aa(i1,i2) = aa(i1,i2) - drt
      aa(i2,i2) = aa(i2,i2) - drt
      aa(j1,i2) = aa(j1,i2) + drt 

      !write(*,*) ' 2   drt, aa(i1,i2), smwi(i2)'
      !write(*,*) drt, aa(i1,i2), smwi(i2)
                           

      !write(*,*) aeta(38),gas%txf(38)
      !rkf  = exp(aeta(ir)*tafl - athet(ir)*tafi)
      !rkinf = gas%txf(ir)*exp(aeta(ir)*tafl - gas%txb(ir)*tafi) ! highjacked txf=A_inf and txb=Einf/Ru


      !ELSE
      !write(*,*) 'Not solving this reaction'
      !ENDIF

      !write(*,*) 'reached end of reactions'
      
      end subroutine my_user_chem
 

       END MODULE userdata

       MODULE userbcdata 

       implicit none
       save

      contains

!  ***********************************************************************************
!     Subroutine to perform any required boundary condition initialization for
!     user-defined boundary conditions.  This routine should have a select case on
!     the variable ibt, and then loops over faces j1:j2.  If a problem is starting
!     fresh, then this routine might set wall solution variables in bvar.
!
!     Input
!     -----
!     ibc   - The boundary condition number in the bcs structure
!     ibt   - The boundary type from the zone definition: zdefs(6,iz)
!     j1,j2 - The face range for this group
!
!  ***********************************************************************************
      subroutine my_user_bc_init(ibc,ibt,j1,j2,ier)
      use US3D_EXTRAS
      use sizing
      use constants
      use connect
      use geometry
      use simvars
      use models
      use switches
      use HDF5
      use userdata

      implicit none

      integer, intent(IN) :: ibc,ibt,j1,j2
      integer, intent(OUT) :: ier

      real(8) :: inflowplane(3,4), hexcell(3,8),interpoly(3,6),normal(3)
      real(8) :: d,area,areatot,rs(ns),vel(3),flux(ns+7),ti,ret,ke,cvt
      real(8) :: fluxtot(ns+7),hfg,x1(3),x2(3),px(3),tvi,trb(2),evi,rav
      integer :: nn_o,nc_o,nf_o,ng_o,itype,ii,nsv_o,i,j,jj,n,npoly
      integer :: ngl,iif,iloc

      real(8), allocatable :: xcn_other(:,:),data_other(:,:)
      integer, allocatable :: ige(:),ien_other(:,:),igl(:)

      type(gridh5_file) :: otherfile
      type(us3d_sfile)  :: othersol
      integer(HID_T)    :: rid

      ! -- Perform a select case on ibt, then loop over face range j1:j2, set values
      ! -- of all variables bvar on first initialization.  This might be skipped on
      ! -- problem restarts.
      select case (ibt)
      case(-12) ! 2D interpolated Inflow (@author: Nick Gibbons)
        !if (iconr/=0) return
        ! Find the bounding box of this face (For optimisation)
        call bbx(j1, j2, x1, x2)

        ! Open the donor grid file and read solution information
        call h5open_f(ier)
        if (ier/=0) goto 998
        call us3d_gridh5f_init(bcs%bcv(ibc)%cvals(1),.false.,otherfile,ier)
        if (ier/=0) goto 998
        call us3d_file_gir(otherfile%giid,nn_o,nc_o,nf_o,ng_o,ier)

        ! Do allocation of node arrays
        allocate(xcn_other(3,nn_o),ien_other(0:8,nc_o),ige(nn_o),igl(nc_o),STAT=ier)
        if (ier/=0) goto 998

        do i = 1,nn_o 
            ige(i) = i   ! Pointless array made because of reasons
        enddo

        ! Read Node positions from other file (Can we do this without a barrier? Yes!)
        call us3d_data_r(otherfile%fid,'xcn',nn_o,nn_o,xcn_other,ige,0,ier)
        call us3d_gridh5f_cnd(otherfile,ier) ! Close grid.h5

        ! Now read ien from conn file (Should be fine to reuse ige)
        call us3d_gridh5f_init(bcs%bcv(ibc)%cvals(2),.false.,otherfile,ier)
        if (ier/=0) goto 998
        call us3d_data_r(otherfile%fid,'ien',nc_o,nc_o,ien_other,ige,0,ier)
        call us3d_gridh5f_cnd(otherfile,ier) ! Close conn.h5

        ! Determine the cells which are relevant to our face.
        ngl = 0
        do i=1,nc_o
            ! Assemble cell node co-ords
            do n=1,8
                hexcell(:,n) = xcn_other(:,ien_other(n,i))
            enddo 
            if (cell_bbx_overlap(x1,x2,hexcell)) then
               igl(ngl+1) = i
               ngl = ngl + 1
            endif
        enddo    ! Possibly resize igl
        if (ngl.eq.0) then
           write(*,*) "Didn't find any cells inside the boundary box..."
           goto 998
        endif

        ! Now import solution data from the otherfile 
        call us3d_sfile_init(bcs%bcv(ibc)%cvals(3),.false.,othersol,ier)
        if (ier/=0) goto 998
        if (.not.check_consistent(othersol)) goto 998
        call h5ex_att_get(othersol%siid,'nsv',nsv_o,ier)
        if (ier/=0) goto 998
        allocate(data_other(nsv_o,nc_o),STAT=ier)
        if (ier/=0) goto 998

        ! We only want to load cells in the igl array
        call us3d_sfile_orun(othersol,1,rid,ier)
        call us3d_data_r(rid,'interior',nc_o,ngl,data_other,igl,0,ier)
        if (ier/=0) goto 998

        call us3d_sfile_cnd(othersol,ier)
        call h5close_f(ier) ! Finished data readin

        ! Loop through inflow plane points
        do j = j1,j2
            jj = j-nif
            i = ife(j,2)
            areatot = 0.0d0
            rav = 0.0d0
            fluxtot(:) = 0.0d0
            iif = 0

            ! Assemble boundary plane nodes
            do n = 1,4
               inflowplane(:,n) = xcn(:,ifn(n,j))
            enddo
            call hessiannormal(inflowplane, normal, d)

            ! Loop through other grid cells, note that the solution data
            ! uses a different indexing system to the other grid geometry arrays
            ! such as xcn and ifn. This is messy, but prevents memory
            ! issues when loading a large other grid on many processors. The 
            ! scheme is also a drastic improvement in speed because of the 
            ! reduced number of cells that need to be searched.
            do iloc=1,ngl       ! Local cell id for reading data_other
                ii = igl(iloc)  ! Global cell id for reading everything else 
           
                ! Assemble a hexcell node coords
                 do n = 1,8
                    hexcell(:,n) = xcn_other(:,ien_other(n,ii))
                 enddo

                 ! Determine if cell intersects infinite plane
                 call hexintersect(hexcell,normal,d,interpoly,npoly)
                 if (npoly.lt.3) cycle
                 if (npoly.gt.6) then
                     write(*,*) "Impossible number of points found."
                     goto 998
                 endif

                 ! The tricky bit: find intersectarea, which might be zero!
                 call intersectarea(inflowplane,normal,interpoly,npoly,area)
                 if (area.lt.1d-32) cycle

                 ! If we've reached here cell ii,iloc intersects this face
                 rs(:)  = data_other(1:ns,iloc)       ! kg/m3
                 vel(:) = data_other(ns+1:ns+3,iloc)  ! m/s
                 ti     = data_other(ns+4,iloc)       ! K
                 select case(ivib)
                     case(0)
                         tvi = ti
                     case(1:2)
                         tvi = data_other(ns+5,iloc)
                     case(3)
                         tvi = ti
                 end select
                 select case(itrb)
                     case (0)
                       trb(:)=0.0d0
                     case (1:5)
                       trb(1) = data_other(ns+4+nv+1, iloc)         ! sa-mut=rho*nu_hat
                       trb(2) = 0.0d0 ! No kinetic energy with SA
                     case (10)
                       trb(1) = data_other(ns+4+nv+1, iloc)         ! sst-rhow
                       trb(2) = data_other(ns+4+nv+2, iloc)         ! sst-rhok
                     case default
                      write(*,*) "Invalid Turbulence Model for ibt=-12"
                      goto 998
                 end select

                 flux = fluxnormal(rs,vel,ti,tvi,trb,normal)
                 fluxtot = fluxtot + flux*area
                 rav = rav + sum(rs)*area ! Use this to decide between solutions later
                 areatot = areatot + area

                 iif = iif + 1
            enddo

            if (iif.eq.0) then
                jj = j-1-nif
                i = ife(j-1,2)
                areatot = 0.0d0
                rav = 0.0d0
                fluxtot(:) = 0.0d0
                iif = 0
    
                ! Assemble boundary plane nodes
                do n = 1,4
                   inflowplane(:,n) = xcn(:,ifn(n,j-1))
                enddo
                call hessiannormal(inflowplane, normal, d)
    
                ! Loop through other grid cells, note that the solution data
                ! uses a different indexing system to the other grid geometry arrays
                ! such as xcn and ifn. This is messy, but prevents memory
                ! issues when loading a large other grid on many processors. The 
                ! scheme is also a drastic improvement in speed because of the 
                ! reduced number of cells that need to be searched.
                do iloc=1,ngl       ! Local cell id for reading data_other
                    ii = igl(iloc)  ! Global cell id for reading everything else 
               
                    ! Assemble a hexcell node coords
                     do n = 1,8
                        hexcell(:,n) = xcn_other(:,ien_other(n,ii))
                     enddo
    
                     ! Determine if cell intersects infinite plane
                     call hexintersect(hexcell,normal,d,interpoly,npoly)
                     if (npoly.lt.3) cycle
                     if (npoly.gt.6) then
                         write(*,*) "Impossible number of points found."
                         goto 998
                     endif
    
                     ! The tricky bit: find intersectarea, which might be zero!
                     call intersectarea(inflowplane,normal,interpoly,npoly,area)
                     if (area.lt.1d-32) cycle
    
                     ! If we've reached here cell ii,iloc intersects this face
                     rs(:)  = data_other(1:ns,iloc)       ! kg/m3
                     vel(:) = data_other(ns+1:ns+3,iloc)  ! m/s
                     ti     = data_other(ns+4,iloc)       ! K
                     select case(ivib)
                         case(0)
                             tvi = ti
                         case(1:2)
                             tvi = data_other(ns+5,iloc)
                         case(3)
                             tvi = ti
                     end select
                     select case(itrb)
                         case (0)
                           trb(:)=0.0d0
                         case (1:5)
                           trb(1) = data_other(ns+4+nv+1, iloc)         ! sa-mut=rho*nu_hat
                           trb(2) = 0.0d0 ! No kinetic energy with SA
                         case (10)
                           trb(1) = data_other(ns+4+nv+1, iloc)         ! sst-rhow
                           trb(2) = data_other(ns+4+nv+2, iloc)         ! sst-rhok
                         case default
                          write(*,*) "Invalid TurbModel for ibt=-12"
                          goto 998
                     end select
    
                     flux = fluxnormal(rs,vel,ti,tvi,trb,normal)
                     fluxtot = fluxtot + flux*area
                     rav = rav + sum(rs)*area ! Use this to decide between solutions later
                     areatot = areatot + area
    
                     iif = iif + 1
                enddo
    
    
                write(*,*) "Error:Faceisoutsidedonorgrid.CAUTIONHACK"
                !goto 998
            endif
            flux = fluxtot/areatot
            rav = rav/areatot
            call inverseflux(flux,normal,rs,vel,ret,ti,tvi,trb,rav)

           if (sum(rs).lt.0.0d0) then
               write(*,*) "Weird result in face: ", j
               write(*,*) "Rhos: ", rs(:)
               write(*,*) "T:", ti
               write(*,*) "Tv:", tvi
               write(*,*) "vel:", vel
               write(*,*) "re:", ret
               write(*,*) "flux:", flux(:)
               write(*,*) "Cells in calculation: ", iif
               stop
           endif

            select case (itrb)
            !case(0) 
            !   rmut(:,i) = 1.0d-16
            case (1:5)
               rmut(1,i) = trb(1)
            case (10)
               rmut(1,i) = trb(1)
               rmut(2,i) = trb(2)
            end select


            r(i) = 0.0 ! Looks good. Copy results into the ghost cells
            do n = 1,ns
               r(i) = r(i) + rs(n)
            enddo
            cs(:,i) = rs(:)/r(i)
            
            u(i)   = vel(1)
            v(i)   = vel(2)
            w(i)   = vel(3)
            t(i)   = ti
            tv(i)  = tvi
            ! With rhos, uvw, T and Tv, all other variables are calculated here
            ! Put this into a subroutine and pretty it.
            call gas_therm1r(gas,rs(1:ns),ti,p(i),gcon(i),c(i),hfg,cv(i))

            select case (ivib)
                case(0) 
                   ev(i) = 0.0d0
                case(1)
                   call gas_vibe1(gas,rs(1:ns),tv(i),ev(i))  
                case(2)
                   call gas_vibel1(gas,rs(1:ns),tv(i),ev(i))
                case(3)
                   call gas_vibel1(gas,rs(1:ns),tv(i),ev(i))
                   cvt = cv(i) ! Gas therm sets cv and something weird happens
                   cv(i) = dot_product(rs,gas%cvs)/r(i)
            end select 

            ke  = 0.5d0*(vel(1)**2 + vel(2)**2 + vel(3)**2)
            re(i) = r(i)*(ke + cv(i)*t(i) + ev(i) + hfg) + trb(2)
            rnu(:,i) = rmut(:,i)/r(i)

            if (ivib.eq.3) then
                cv(i) = cvt
                ev(i) = 0.0d0
            endif

            bvar(1:ns,jj) = rs(1:ns)
            bvar(ns+1,jj) = u(i)
            bvar(ns+2,jj) = v(i)
            bvar(ns+3,jj) = w(i)
            bvar(ns+4,jj) = t(i)
            do n= 1,nv   ! FIX THIS (Why?)
               bvar(ns+4+n,jj) = tv(i)
            enddo
            do n= 1,nt 
               bvar(ns+4+nv+n,jj) = rmut(n,i)
            enddo

        enddo ! End Face Loop

        if (allocated(xcn_other))  deallocate(xcn_other)
        if (allocated(ien_other))  deallocate(ien_other)
        if (allocated(ige))        deallocate(ige)
        if (allocated(igl))        deallocate(igl)
        if (allocated(data_other)) deallocate(data_other)
      end select ! End of User BC Initialization
      return
 998  write(6,*) " *** Error detected in userbc subroutine.***"
      write(*,*) " Check filenames and grid compatibility "
      write(*,*) " IBC is: ", ibc
      stop
      return
      end subroutine my_user_bc_init !"


       END MODULE userbcdata





!  ***********************************************************************************
!     User initialization routine
!  ***********************************************************************************

      subroutine user_initialize(component,debug,ier)
      use us3d_user_module, only : user_chem_reaction, user_bc_init
      use userdata
      use userbcdata !, only : my_user_bc_init
      use geo2D
      implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier

      ier= 0
      user_chem_reaction => my_user_chem
      user_bc_init => my_user_bc_init

      return
      end subroutine user_initialize


