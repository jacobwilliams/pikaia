program xpk1b

   !! Driver program for non-linear least-squares problem
   !! (Sect. 5.3)

   use pikaia_module, only: pikaia_class
   use,intrinsic :: iso_fortran_env, only: wp => real64

implicit  none

integer,parameter :: n = 17

integer :: seed, i, status
real(wp) :: x(n), f
type(pikaia_class) :: solver
real :: f_data(200), t(200) , sigma ! single precision from file
integer :: ndata, m

real(wp),dimension(n),parameter :: xl = 0.0_wp
real(wp),dimension(n),parameter :: xu = 1.0_wp

! Initializations

! Set control variables (use defaults except for population size)
   call solver%init(n,xl,xu,fit1b,status,&
                    np    = 50,&
                    iseed = 13579)
   call finit()

! Now call pikaia
call solver%solve(x,f,status)

! Print the results
write(*,*) ' status: ',status
write(*,*) '      x: ',x
write(*,*) '      f: ',f

contains

subroutine finit()

   !! Reads in synthetic dataset (see Figure 5.4)

real :: f0(200),vdum(11),delt
integer :: i, iunit

open(newunit=iunit,file='test/syndat1.i3e',form='unformatted',status='OLD',convert='big_endian')
read(iunit) ndata
read(iunit) (vdum(i),i=1,11)
read(iunit) (t(i),i=1,ndata)
read(iunit) (f0(i),i=1,ndata)
read(iunit) (f_data(i),i=1,ndata)
close(iunit)

! Use 5 Fourier modes for the fit
m=5

! same error bar for all point
sigma=5.
delt=t(3)-t(1)

end subroutine finit

subroutine fit1b(me,x,f)

   !! Fitness function for non-linear least squares problem
   !! (Sect. 5.3)

   class(pikaia_class),intent(inout) :: me
   real(wp),dimension(:),intent(in)  :: x
   real(wp),intent(out)              :: f

   integer,parameter :: MMAX=10
   real(wp),parameter :: pi = acos(-1.0_wp)

   INTEGER :: N , i , j
   REAL(wp) :: amp(MMAX) , per(MMAX) , a , b , sum , sum2 , nyp

   !---------- 1. rescale input variables:
      a = X(1)*10.0_wp
      b = X(2)*100.0_wp
      nyp = 2.0_wp*(T(2)-T(1))
      DO j = 1 , M
         amp(j) = X(3*j)*100.0_wp
         per(j) = X(3*j+1)*(50.0_wp-nyp) + nyp
      ENDDO
   !---------- 2. compute chi**2
      sum = 0.0_wp
      DO i = 1 , Ndata
         sum2 = 0.0_wp
         DO j = 1 , M
            sum2 = sum2 + amp(j)*sin(2.*PI*(T(i)/per(j)+X(3*j+2)))
         ENDDO
         sum = sum + ((a*T(i)+b+sum2-F_data(i))/Sigma)**2
      ENDDO
   !---------- 3. define fitness
      f = 1.0_wp/sum

   END subroutine fit1b

end program xpk1b