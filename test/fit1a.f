
      function fit1a(n,x)
c========================================================
c     Fitness function for linear least squares problem
c     (Sect. 5.2)
c========================================================
      implicit none
      common/data/ f(200),t(200),sigma,ndata
      integer      n,i,ndata
      real         x(n),fit1a,f,t,sigma,a,b,sum
c---------- 1. rescale input variables:
      a=x(1)*10.
      b=x(2)*100.
c---------- 2. compute chi**2
      sum=0.
      do 1 i=1,ndata
         sum=sum+ ( (a*t(i)+b-f(i))/sigma )**2
    1 continue
c---------- 3. define fitness
      fit1a=1./sum

      return
      end


