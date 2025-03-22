
program xpk2

   !! Driver program for ellipse fitting problem (Sect. 5.4)

   use pikaia_module, only: pikaia_class
   use,intrinsic :: iso_fortran_env, only: wp => real64

implicit none

integer,parameter :: n = 5

integer ::  i, status
real(wp) :: x(n), f
type(pikaia_class) :: solver
integer :: ndata
real(wp),dimension(:),allocatable :: xdat, ydat !! from the data file

real(wp),dimension(n),parameter :: xl = 0.0_wp
real(wp),dimension(n),parameter :: xu = 1.0_wp
real(wp),parameter :: pi = acos(-1.0_wp)

! Initializations
!     Set control variables (50 individuals for 200 generations
!     under Select-random-delete-worst reproduction plan)
call solver%init(n,xl,xu,fit2,status,&
                  ngen  = 200,&
                  np    = 50,&
                  iseed = 654321, &
                  irep  = 3 )
                  ! ctrl(9)=0.0   !? what was this ?
                  ! ctrl(10)=3    !? is this irep ?

call finit()

! Now call pikaia
call solver%solve(x,f,status)

! Print the results
write(*,*) ' status: ',status
write(*,*) '      x: ',x
write(*,*) '      f: ',f

contains

subroutine finit()

   !! Reads in synthetic dataset (see Figure 5.9)

   use json_module

   type(json_file) :: json

   call json%load(filename='test/syndat2.json')
   call json%get('Ndata',  Ndata)
   call json%get('xdat',   xdat)
   call json%get('ydat',   ydat)
   call json%destroy()

end subroutine finit

pure function fatan(yy,xx)
   !! Returns arctangent in full circle

   real(wp),intent(in) :: yy,xx
   real(wp) :: fatan

   real(wp) :: a1

   a1=atan(yy/xx)
   if(xx < 0.0_wp) a1 = a1+pi
   if(a1 < 0.0_wp) a1 = 2.0_wp*pi+a1
   fatan=a1

   end function fatan

   subroutine fit2(me,x,f)

      !! Fitness function for ellipse fitting problem (Sect. 5.4)
      !!
      !! Ellipse parameters are:
      !! x(1)=x_0, x(2)=y_0, x(3)=a, x(4)=b, x(5)=theta_0

      class(pikaia_class),intent(inout) :: me
      real(wp),dimension(:),intent(in)  :: x
      real(wp),intent(out)              :: f

         INTEGER :: i
         REAL(wp) :: x0 , y0 , a2 , b2 , theta0 , distdat , angdat , rthj , sum

      !---------- 1. rescale input variables:
      !           0 <= x(1...4) <= 2, 0 <= x(5) <= pi
         x0 = X(1)*2.0_wp
         y0 = X(2)*2.0_wp
         a2 = (X(3)*2.0_wp)**2
         b2 = (X(4)*2.0_wp)**2
         theta0 = X(5)*PI
      !---------- 2. compute merit function
         sum = 0.0_wp
         DO i = 1 , Ndata
      !     (a) compute distance d_j from center to data point
            distdat = sqrt((x0-Xdat(i))**2+(y0-Ydat(i))**2)
      !     (b) compute angle theta_j of segment center---data point
            angdat = fatan((Ydat(i)-y0),(Xdat(i)-x0)) - X(5)*PI
      !     (c) compute radius r(\theta_j) of ellipse at that angle
            rthj = sqrt(a2*b2/(a2*sin(angdat)**2+b2*cos(angdat)**2))
      !     (d) increment DR merit function
            sum = sum + (distdat-rthj)**2
         ENDDO
      !---------- 3. equate fitness to inverse of DR merit function
         f = 1.0_wp/sum

      END subroutine fit2

end program xpk2