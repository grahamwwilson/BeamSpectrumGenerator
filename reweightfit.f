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
      include 'mcdefn.f'

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
      open(unit=22,file="1m/testbc-1-1.dat",status='old')
      call readheader(22)
      n=0
   10 continue
      read(22,*,end=999)iev,rtype,x1,x2,y1,y2,z1,z2,
     +                   weight,weightp,rwt,ibin1,ibin2
*      print *,'iev,ibin1,ibin2 ',iev,ibin1,ibin2
      n = n+1
* Store info
      iv(1,n) = ibin1
      iv(2,n) = ibin2
      iv(3,n) = rtype
      v(1,n) = x1
      v(2,n) = y1
      v(3,n) = z1
      v(4,n) = x2
      v(5,n) = y2
      v(6,n) = z2
      v(7,n) = weight
      goto 10
  999 continue
      close(22)
      
      print *,'Number of events read ',n  

      open(unit=8,file='reweightfit.dat',status='old')

      call mintio(8,6,7)
  
      call minuit(fcn)

      end

      subroutine readheader(lun)
      implicit none
      integer lun
      integer nheader,version,nparameters,nevs,seedm,seedl
      double precision mu1,mu2,s1,s2,pnorm1,pnorm2,alphab,betab
      double precision alphaa,betaa
      
      read(lun,*)nheader
      read(lun,*)version
      read(lun,*)nparameters
      read(lun,*)nevs
      read(lun,*)seedm
      read(lun,*)seedl
      read(lun,*)mu1
      read(lun,*)mu2
      read(lun,*)s1
      read(lun,*)s2
      read(lun,*)pnorm1
      read(lun,*)pnorm2
      read(lun,*)alphab
      read(lun,*)betab
      read(lun,*)alphaa
      read(lun,*)betaa
      
      end
      
      subroutine fcn(npar,grad,fval,x,iflag)
* Calculate the chi-squared      
      implicit none
      integer npar,iflag
      double precision grad,x,fval
      dimension grad(*),x(*)
      double precision ppeak,pbody,alphab,betab,alphaa,betaa
      double precision parm
      double precision pregion
      double precision bbody(4),barms(4)
      double precision normb,norma
      integer ia,ib,ja,jb
      parameter (ia=1, ib=2,ja=3,jb=4)      
      double precision wtp,py1,py2
      integer i,j,ibin
      include 'exptdefn.f'
      include 'mcdefn.f'
      include 'fcnlocal.f'
      double precision ff
      double precision w1sum(nbins),w2sum(nbins)
      double precision rwt
      double precision chisq1,chisq2
      double precision ndata1,nexp1,var1,chi1
      double precision ndata2,nexp2,var2,chi2
      integer myflag

* Copy fit parameters into more amenable variables
      ppeak = x(1)
      pbody = x(2)
      alphab = x(3)
      betab = x(4)
      alphaa = x(5)
      betaa = x(6)
      ff = x(7)
      parm = 0.5d0*(1.0d0-ppeak-pbody)
      call cpbeta(bbody,alphab,betab,normb)
      call cpbeta(barms,alphaa,betaa,norma)
      print *,'Pars: ',x(1),x(2),x(3),x(4),x(5),x(6),x(7)
      
      fval = 0.0d0
      do j=1,nbins
         w1sum(j) = 0.0d0
         w2sum(j) = 0.0d0
      enddo
      
* Loop over MC events
      do i=1,nmc
         myflag = 0
* Define variables
         ibin1 = iv(1,i)
         ibin2 = iv(2,i)
         rtype = iv(3,i)
         x1 = v(1,i)
         y1 = v(2,i)
         z1 = v(3,i)
         x2 = v(4,i)
         y2 = v(5,i)
         z2 = v(6,i)
         wt = v(7,i)
* Calculate appropriate weight factor for the new parameters
         if(rtype.eq.1)then
* peak
            wtp = ppeak
            pregion = ppeak
         elseif(rtype.eq.2)then
* body
            py1=((1d0-y1)**bbody(ja))*(y1**bbody(jb))/normb
            py2=((1d0-y2)**bbody(ja))*(y2**bbody(jb))/normb      
            wtp = pbody*py1*py2
            pregion = pbody
         elseif(rtype.eq.3)then
* arms 
            py1=((1d0-y1)**barms(ja))*(y1**barms(jb))/norma
            wtp = parm*py1
            pregion = parm
         else
            py2=((1d0-y2)**barms(ja))*(y2**barms(jb))/norma
            wtp = parm*py2
            pregion = parm         
         endif
* Check for infinity issue
         if(wtp-1d0.eq.wtp)then
*            print *,'Infinity issue. Event: ',i,' region: ',rtype
            myflag = 1
         endif
         
         rwt =wtp/wt
* Sum the weights for the x1 and x2 distributions
         if(myflag.eq.0)then
            w1sum(ibin1) = w1sum(ibin1) + rwt
            w2sum(ibin2) = w2sum(ibin2) + rwt
         endif
      enddo

* Now construct the chi-squared 
      chisq1 = 0.0d0
      chisq2 = 0.0d0
      do j=1,nbins
         ndata1 = dble(cont1(j))
         nexp1  = ff*w1sum(j)
         var1 = nexp1
         chi1 = (ndata1 - nexp1)/sqrt(var1)
         chisq1 = chisq1 + chi1*chi1
         
         ndata2 = dble(cont2(j))
         nexp2  = ff*w2sum(j)
         var2 = nexp2
         chi2 = (ndata2 - nexp2)/sqrt(var2)         
         chisq2 = chisq2 + chi2*chi2
      enddo
     
      print *,'Chisq values ',chisq1,chisq2,chisq1+chisq2
      fval = chisq1 + chisq2
      
      if(iflag.eq.1)then
         print *,'Minimization Initialization '
         print *,'fval0 = ',fval
         do ibin=1,nbins
            print *,'Bin ',ibin,' entries ',
     +              nint(cont1(ibin)),nint(cont2(ibin))
         enddo
      endif

      if(iflag.eq.3)then
         print *,'Minimization Finished '
         print *,'fval = ',fval
      endif
      
      end
      
      include 'cpbeta.f'
      
