
      function fit1b(n,x)
c==========================================================
c     Fitness function for non-linear least squares problem
c     (Sect. 5.3)
c==========================================================
      implicit none
      common/data/ f(200),t(200),sigma,ndata,m
      integer      n,m,mmax,ndata,i,j
      real         x(n),fit1b,f,t,sigma,amp,per,a,b,
     +             sum,sum2,nyp,pi
      parameter    (mmax=10,pi=3.1516926536)
      dimension    amp(mmax),per(mmax)
c---------- 1. rescale input variables:
      a=x(1)*10.
      b=x(2)*100.
      nyp=2.*(t(2)-t(1))
      do 10 j=1,m
         amp(j)=x(3*j)*100.
         per(j)=x(3*j+1)*(50.-nyp)+nyp
  10  continue
c---------- 2. compute chi**2
      sum=0.
      do 1 i=1,ndata
         sum2=0.
         do 2 j=1,m
            sum2=sum2+amp(j)*sin(2.*pi*(t(i)/per(j)+x(3*j+2)))
    2    continue
         sum=sum+ ( (a*t(i)+b+sum2-f(i))/sigma )**2
    1 continue
c---------- 3. define fitness
      fit1b=1./sum

      return
      end

