      program reweightfit
*      
* Plan
*
* 1. Read in unweighted data sampled from alternative model (P2) using 
*    histogram
*
* 2. Chi-squared fit to the unweighted data based on reweighting MC events 
*    from model generated with model P1. 
*
      implicit none
*
* Graham W. Wilson,     December 9th, 2021
*

      integer nwpawc
      parameter (nwpawc=1000000)
      real hmemor
      common/pawc/hmemor(nwpawc)
      
      include 'exptdefn.f'
      
      integer n
      integer iev,rtype,ibin1,ibin2
      double precision x1,x2,y1,y2,z1,z2,weight,weightp,rwt
      
      integer i

      external fcn

      call hlimit(nwpawc)
      call hidopt(0,'stat')
      
* Read data from histogram to be fitted 
* and store in expt common block (see exptdefn.f include file) for use in fcn
      call hrget(0,
     +  '/home/graham/BeamSpectrumGenerator/100k/testbc-2-2.hbook',' ')
      call hprint(107)
      call hunpak(107,cont1,'hist',1)
      call hprint(108)
      call hunpak(108,cont2,'hist',1)                
      print *,'Data to fit                        x1          x2 '
      do i=1,nbins
         print *,'Bin ',i,' entries ',nint(cont1(i)),nint(cont2(i))
      enddo
      
* Next read the pre-generated MC events and store relevant info

      open(unit=21,file='1M/testbc-2-2.dat',status='old')
      n=0
   10 continue
      read(21,*,end=999)iev,rtype,x1,x2,y1,y2,z1,z2,
     +                   weight,weightp,rwt,ibin1,ibin2
      n = n+1
      goto 10
  999 continue
      close(21)
      
      print *,'Number of events read ',n  

      open(unit=8,file='reweightfit.dat',status='old')

      call mintio(8,6,7)
  
      call minuit(fcn)

      end

      subroutine fcn(npar,grad,fval,x,iflag)
* Calculate the chi-squared      
      implicit none
      integer npar,iflag
      double precision grad,x,fval
      dimension grad(*),x(*)
      double precision ppeak,pbody,alphab,betab,alphaa,betaa
      double precision bbody(4),barms(4)
      double precision normb,norma
      integer i
      include 'exptdefn.f'

      fval = 0.0d0

* Copy fit parameters into more amenable variables
      ppeak = x(1)
      pbody = x(2)
      alphab = x(3)
      betab = x(4)
      alphaa = x(5)
      betaa = x(6)
      call cpbeta(bbody,alphab,betab,normb)
      call cpbeta(barms,alphaa,betaa,norma)
      
      if(iflag.eq.1)then
         print *,'Minimization Initialization '
         print *,x(1)
         print *,'fval0 = ',fval
         do i=1,nbins
            print *,'Bin ',i,' entries ',nint(cont1(i)),nint(cont2(i))
         enddo
      endif

      if(iflag.eq.3)then
         print *,'Minimization Finished '
         print *,x(1)
         print *,'fval = ',fval
      endif
      
      end
      
      include 'cpbeta.f'
      
