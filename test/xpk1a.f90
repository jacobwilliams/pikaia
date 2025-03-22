PROGRAM xpk1a
!================================================================
!     Driver program for linear least-squares problem
!     (Sect. 5.2)
!================================================================
   use pikaia_module, only: pikaia_class
   use,intrinsic :: iso_fortran_env, only: wp => real64

   IMPLICIT NONE

   INTEGER N , seed , i , status
   PARAMETER (N=2)

   REAL(wp) :: x(N) , f
   type(pikaia_class) :: solver
   real :: f_data(200), t(200) , sigma ! single precision from file
   integer :: ndata

   real(wp),dimension(n),parameter :: xl = 0.0_wp
   real(wp),dimension(n),parameter :: xu = 1.0_wp

! Initializations

! Set control variables (evolve 50 individuals over 100
! generations, use defaults values for other input parameters)
   call solver%init(n,xl,xu,fit1a,status,&
                    ngen  = 100,&
                    np    = 50,&
                    iseed = 123456)
   CALL finit()

! Now call pikaia
   call solver%solve(x,f,status)

! Print the results
   WRITE (*,*) ' status: ' , status
   WRITE (*,*) '      x: ' , x
   WRITE (*,*) '      f: ' , f

contains

SUBROUTINE finit()

   !! Reads in synthetic dataset (see Figure 5.4)

   real :: f0(200) , vdum(11)
   integer :: i
   integer :: iunit

   ! todo: convert these data files to JSON...
   OPEN (newunit=iunit,FILE='test/syndat1.i3e',FORM='unformatted',STATUS='old',convert='big_endian')

   READ (iunit) Ndata
   READ (iunit) (vdum(i),i=1,11)
   READ (iunit) (T(i),i=1,Ndata)
   READ (iunit) (f0(i),i=1,Ndata)
   READ (iunit) (f_data(i),i=1,Ndata)
   Sigma = 5.0

   close(iunit)

END SUBROUTINE finit

subroutine fit1a(me,x,f)

   !! Fitness function for linear least squares problem
   !! (Sect. 5.2)

   class(pikaia_class),intent(inout) :: me
   real(wp),dimension(:),intent(in)  :: x
   real(wp),intent(out)              :: f

      INTEGER :: i
      REAL(wp) :: a , b , sum

   ! 1. rescale input variables:
      a = X(1)*10.0_wp
      b = X(2)*100.0_wp
   ! 2. compute chi**2
      sum = 0.0_wp
      DO i = 1 , Ndata
         sum = sum + ((a*T(i)+b-f_data(i))/Sigma)**2
      ENDDO
   ! 3. define fitness
      f = 1.0_wp/sum

   END subroutine fit1a

END PROGRAM xpk1a
