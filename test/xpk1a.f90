PROGRAM xpk1a
   !! Driver program for linear least-squares problem
   !! (Sect. 5.2)

   use pikaia_module, only: pikaia_class
   use,intrinsic :: iso_fortran_env, only: wp => real64

   IMPLICIT NONE

   INTEGER :: i , status
   integer, PARAMETER :: N = 2

   REAL(wp) :: x(N) , f
   type(pikaia_class) :: solver
   real(wp) :: sigma
   integer :: ndata
   real(wp),dimension(:),allocatable :: f_data, t !! from the data file
   real(wp),dimension(n),parameter :: xl = 0.0_wp
   real(wp),dimension(n),parameter :: xu = 1.0_wp

! Initializations

! Set control variables (evolve 50 individuals over 100
! generations, use defaults values for other input parameters)
   call solver%init(n,xl,xu,fit1a,status,&
      ngen  = 100,&
      np    = 50,&
      ivrb  = 1, &
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

      use json_module

      type(json_file) :: json

      Sigma = 5.0_wp

      call json%load(filename='test/syndat1.json')
      call json%get('Ndata',  Ndata)
      call json%get('t',      t)
      call json%get('f',      f_data)
      call json%destroy()

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
