
      program xpk3
c============================================================
c     Driver program for circle fitting problem using
c     a robust estimator based on the Hough transform
c     (Sect. 5.5)
c============================================================
      implicit  none
      integer   n, seed, i, status
      parameter (n=2)
      real      ctrl(12), x(n), f, fit3
      external  fit3
c
c     First, initialize the random-number generator
c
      seed = 123456
      call rninit(seed)
c
c     Read in synthetic data
c
      call finit
c
c     Set control variables (use 
c     
      do 10 i=1,12
         ctrl(i) = -1
   10 continue

c     Now call pikaia
      call pikaia(fit3,n,ctrl,x,f,status)
c
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',f
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)

      end
c***************************************************************
      subroutine finit
c===============================================================
c     Reads in synthetic dataset
c===============================================================
      implicit     none
      common/data/ xdat(200),ydat(200),sigma,ndata
      real         xdat,ydat,sigma
      integer      ndata,i

      open(1,file='syndat3.i3e',form='unformatted')

      read(1) ndata
      read(1) (xdat(i),i=1,ndata)
      read(1) (ydat(i),i=1,ndata)

c     Same error estimate in x and y for all data points
      sigma=0.05

      return
      end

      include "fit3.f"

      include "pikaia.f"
