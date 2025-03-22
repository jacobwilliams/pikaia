
      program xpk2
c==================================================================
c     Driver program for ellipse fitting problem (Sect. 5.4)
c==================================================================
      implicit  none
      integer   n, seed, i, status
      parameter (n=5)
      real      ctrl(12), x(n), f, fit2
      external  fit2
c
c     First, initialize the random-number generator
c
      seed=654321
      call rninit(seed)
c
c     Initializations
c
      call finit
c
c     Set control variables (50 individuals for 200 generations
c     under Select-random-delete-worst reproduction plan)
c
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=50
      ctrl(2)=200
      ctrl(9)=0.0
      ctrl(10)=3

c     Now call pikaia
      call pikaia(fit2,n,ctrl,x,f,status)
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
c     Reads in synthetic dataset (see Figure 5.9)
c===============================================================
      implicit     none
      common/data/ xdat(200),ydat(200),ndata
      real         xdat,ydat
      integer      ndata,i

      open(1,file='syndat2.i3e',form='unformatted')

      read(1) ndata
      read(1) (xdat(i),i=1,ndata)
      read(1) (ydat(i),i=1,ndata)

      return
      end

c***********************************************************
      function fatan(yy,xx)
c===========================================================
c     Returns arctangent in full circle
c===========================================================
      real yy,xx,fatan,a1,pi
      data pi/3.1415926536/

      a1=atan(yy/xx)
      if(xx.lt.0.) a1=a1+pi
      if(a1.lt.0.) a1=2.*pi+a1
      fatan=a1

      return
      end

      include "fit2.f"

      include "pikaia.f"


