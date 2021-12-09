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
      
      integer i

      external fcn

      call hlimit(nwpawc)
      call hidopt(0,'stat')
      
* Read data from histogram to be fitted 
* and store in expt common block (see exptdefn.f include file) for use in fcn
      call hrget(0,
     +  '/home/graham/BeamSpectrumGenerator/100k/testbc-2-2.hbook',' ')
      call hprint(107)
      call hunpak(107,cont,'hist',1)      
      print *,'Data to fit '
      do i=1,nbins
         print *,'Bin ',i,' entries ',nint(cont(i))
      enddo
      
* Next get the pre-generated MC events.      

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
            print *,'Bin ',i,' entries ',nint(cont(i))            
         enddo
      endif

      if(iflag.eq.3)then
         print *,'Minimization Finished '
         print *,x(1)
         print *,'fval = ',fval
      endif
      
      end
      
      include 'cpbeta.f'
      
