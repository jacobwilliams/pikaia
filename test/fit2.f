
      function fit2(n,x)
c=============================================================
c     Fitness function for ellipse fitting problem (Sect. 5.4)
c
c     Ellipse parameters are:
c     x(1)=x_0, x(2)=y_0, x(3)=a, x(4)=b, x(5)=theta_0
c=============================================================
      implicit none
      common/data/ xdat(200),ydat(200),ndata
      integer      n,ndata,i
      real         x(n),fit2,xdat,ydat,fit2,fatan,x0,y0,a2,b2,
     +             theta0,distdat,angdat,rthj,sum,pi
      parameter    (pi=3.1415926536)
c---------- 1. rescale input variables:
c           0 <= x(1...4) <= 2, 0 <= x(5) <= pi
      x0    = x(1)*2.
      y0    = x(2)*2.
      a2    =(x(3)*2.)**2
      b2    =(x(4)*2.)**2
      theta0= x(5)*pi
c---------- 2. compute merit function
      sum=0.
      do 1 i=1,ndata
c     (a) compute distance d_j from center to data point
         distdat=sqrt((x0-xdat(i))**2+(y0-ydat(i))**2)
c     (b) compute angle theta_j of segment center---data point
         angdat=fatan((ydat(i)-y0),(xdat(i)-x0))-x(5)*pi
c     (c) compute radius r(\theta_j) of ellipse at that angle
         rthj=sqrt(a2*b2/(a2*sin(angdat)**2+b2*cos(angdat)**2))
c     (d) increment DR merit function
         sum=sum+(distdat-rthj)**2
    1 continue
c---------- 3. equate fitness to inverse of DR merit function
      fit2=1./sum

      return
      end

