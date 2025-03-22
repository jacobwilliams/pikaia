      program xpk1b
c==================================================================
c     Driver program for non-linear least-squares problem
c     (Sect. 5.3)
c==================================================================
      implicit  none
      integer   n, seed, i, status
      parameter (n=17)
      real      ctrl(12), x(n), f, fit1b
      external  fit1b
c
c     First, initialize the random-number generator
c
      seed=13579
      call rninit(seed)
c
c     Initializations
c
      call finit
c
c     Set control variables (use defaults except for population size)
c      
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=50

c     Now call pikaia
      call pikaia(fit1b,n,ctrl,x,f,status)
c
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',f
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)
c
      end
c***************************************************************
      subroutine finit
c===============================================================
c     Reads in synthetic dataset (see Figure 5.4)
c===============================================================
      implicit     none
      common/data/ f(200),t(200),sigma,ndata,m
      dimension    f0(200),vdum(11)
      real         f,t,f0,vdum,sigma,delt
      integer      ndata,m,i

      open(1,file='syndat1.i3e',form='unformatted')

      read(1) ndata
      read(1) (vdum(i),i=1,11)
      read(1) (t(i),i=1,ndata)
      read(1) (f0(i),i=1,ndata)
      read(1) (f(i),i=1,ndata)

c     Use 5 Fourier modes for the fit
      m=5

c     same error bar for all point
      sigma=5.
      delt=t(3)-t(1)

      return
      end

      include "fit1b.f"

      include "pikaia.f"
