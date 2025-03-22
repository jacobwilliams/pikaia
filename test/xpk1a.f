      program xpk1a
c================================================================
c     Driver program for linear least-squares problem
c     (Sect. 5.2)
c================================================================
      implicit  none
      integer   n, seed, i, status
      parameter (n=2)
      real      ctrl(12), x(n), f, fit1a
      external  fit1a
c
c     First, initialize the random-number generator
c
      seed = 123456
      call rninit(seed)
c
c     Initializations
c
c     open(1,file='fake.i3e',form='unformatted')
      call finit
c
c     Set control variables (evolve 50 individuals over 100
c     generations, use defaults values for other input parameters)
c
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=50
      ctrl(2)=100
c
c     Now call pikaia
c
      call pikaia(fit1a,n,ctrl,x,f,status)
c
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',f
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)

      end
c*********************************************************************
      subroutine finit
c
c     Reads in synthetic dataset (see Figure 5.4)
c
      implicit none
      common/data/ f(200),t(200),sigma,ndata
      dimension    f0(200),vdum(11)
      real         f,t,sigma,f0,vdum
      integer      ndata,i

      open(1,file='syndat1.i3e',form='unformatted')

      read(1) ndata
      read(1) (vdum(i),i=1,11)
      read(1) (t(i),i=1,ndata)
      read(1) (f0(i),i=1,ndata)
      read(1) (f(i),i=1,ndata)
      sigma=5.

      return
      end

      include "fit1a.f"

      include "pikaia.f"
