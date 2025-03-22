
program xpk3

   !! Driver program for circle fitting problem using
   !! a robust estimator based on the Hough transform
   !! (Sect. 5.5)

   use pikaia_module, only: pikaia_class
   use,intrinsic :: iso_fortran_env, only: wp => real64

   implicit none

   integer,parameter :: n = 2

   integer :: i, status
   real(wp) :: x(n), f
   type(pikaia_class) :: solver
   real(wp) :: sigma
   integer :: ndata
   real(wp),dimension(:),allocatable :: xdat, ydat !! from the data file

   real(wp),dimension(n),parameter :: xl = 0.0_wp
   real(wp),dimension(n),parameter :: xu = 1.0_wp

! initialize
   call solver%init(n,xl,xu,fit3,status,&
      ivrb  = 1, &
      iseed = 123456 )
! Read in synthetic data
   call finit()

! Now call pikaia
   call solver%solve(x,f,status)

! Print the results
   write(*,*) ' status: ',status
   write(*,*) '      x: ',x
   write(*,*) '      f: ',f

contains

   subroutine finit()
      !! Reads in synthetic dataset
      use json_module

      type(json_file) :: json

      ! Same error estimate in x and y for all data points
      sigma=0.05_wp

      call json%load(filename='test/syndat3.json')
      call json%get('Ndata',  Ndata)
      call json%get('xdat',   xdat)
      call json%get('ydat',   ydat)
      call json%destroy()

   end subroutine finit

   subroutine fit3(me,x,f)

      !! Fitness function for circle fitting problem using
      !! a robust estimator based on the Hough transform
      !! (Sect. 5.5)

      class(pikaia_class),intent(inout) :: me
      real(wp),dimension(:),intent(in)  :: x
      real(wp),intent(out)              :: f

      INTEGER :: i
      REAL(wp) :: x0 , y0 , distdat , r1 , r2 , sum

      !---------- 1. rescale input variables:
      !           0 <= x(1,2) <= 3
      x0 = X(1)*3.0_wp
      y0 = X(2)*3.0_wp
      !---------- 2. compute merit function
      r1 = 1.0_wp - Sigma/2.0_wp
      r2 = 1.0_wp + Sigma/2.0_wp
      sum = 0.
      DO i = 1 , Ndata
         !     (a) compute distance d_j from center to data point
         distdat = sqrt((x0-Xdat(i))**2+(y0-Ydat(i))**2)
         !     (b) increment Hough merit function
         IF ( distdat>=r1 .AND. distdat<=r2 ) sum = sum + 1.0_wp
      ENDDO
      !---------- 3. equate fitness to Hough merit function
      f = sum

   END subroutine fit3

end program xpk3
!***************************************************************

