      program rwbinl
*      
* Plan
*
* 1. Read in unweighted using histogram
*
* 2. Fit to the unweighted data based on reweighting MC events 
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
      
      include 'exptdefn2.f'
      include 'mcdefn.f'

* Toggle ichoice to fit different GP or pseudo-data data-sets
      integer ichoice
      include 'ichoice.f'
      
* Toggle jchoice to use different reweighting MC data-sets
      integer jchoice
      include 'jchoice.f'

      integer iev,rtype,ibin1,ibin2,icomb
      double precision x1,x2,y1,y2,z1,z2,weight,weightp,rwt
      
      integer lun
      parameter (lun = 22)
      integer ntoread
      
      integer i

      external fcn

      call hlimit(nwpawc)
      call hidopt(0,'stat')
      
* Read data from histogram to be fitted 
* and store in expt common block (see exptdefn.f include file) for use in fcn

      print *,'ichoice set to ',ichoice
      print *,'jchoice set to ',jchoice

      if(ichoice.eq.-999)then
         call hrget(0,
     +  '/home/graham/BeamSpectrumGenerator/1m/testbc-2-2.hbook',' ')
      elseif(ichoice.eq.0)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run5.hbook',' ')
      elseif(ichoice.eq.1)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run5-minuszv.hbook',' ')
      elseif(ichoice.eq.2)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run5-pluszv.hbook',' ')
      elseif(ichoice.eq.3)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run5-smallzv.hbook',' ')
      elseif(ichoice.eq.4)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run5-bigzv.hbook',' ')
      elseif(ichoice.eq.5)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run7.hbook',' ')
      elseif(ichoice.eq.6)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run7-minuszv.hbook',' ')
      elseif(ichoice.eq.7)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run7-pluszv.hbook',' ')
      elseif(ichoice.eq.8)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run7-smallzv.hbook',' ')
      elseif(ichoice.eq.9)then
         call hrget(0,
     +  '/home/graham/gpDigest/gplumi-run7-bigzv.hbook',' ')     
      else
         print *,'Looks like input dataset is not found!'
      endif

      call hunpak(111,cont,'hist',1)
      ntotdata = 0         
      do i=1,nbins
         ntotdata = ntotdata+ cont(i)
      enddo
      print *,'Number of data events being fitted ',ntotdata
      
* Next read the pre-generated MC events and store relevant info
      print *,'MC array size is: ',nmcarray
      if(jchoice.eq.0)then
         open(unit=lun,file="10m/testbc-1-0.dat",status='old')
      elseif(jchoice.eq.1)then
         open(unit=lun,file="10m/testbc-1-1.dat",status='old')
      elseif(jchoice.eq.2)then
         open(unit=lun,file="1m/testbc-1-1.dat",status='old')         
      endif
      call readheader(lun,ntoread)
      if(ntoread.gt.nmcarray)then
         print *,'MC array not big enough ',nmcarray
         stop
      endif
      nmc=0
   10 continue
      read(lun,*,end=999)iev,rtype,x1,x2,y1,y2,z1,z2,
     +                   weight,weightp,rwt,ibin1,ibin2,icomb
*      print *,'iev,ibin1,ibin2 ',iev,ibin1,ibin2
      nmc = nmc+1
* Store info
      iv(1,nmc) = ibin1
      iv(2,nmc) = ibin2
      iv(3,nmc) = rtype
      iv(4,nmc) = icomb
      v(1,nmc) = x1
      v(2,nmc) = y1
      v(3,nmc) = z1
      v(4,nmc) = x2
      v(5,nmc) = y2
      v(6,nmc) = z2
      v(7,nmc) = weight
      goto 10
  999 continue
      close(lun)
      
      print *,'Number of MC events read ',nmc
      if(nmc.ne.ntoread)then
         print *,'Mis-match in events read!'
      else
         print *,'Expected number of MC events has been read'
      endif

      open(unit=8,file='reweightlfit.dat',status='old')

      call mintio(8,6,7)
  
      call minuit(fcn)

      end

      subroutine readheader(lun,ntoread)
      implicit none
      integer lun
      integer nheader,version,nparameters,nevs,seedm,seedl
      integer ntoread
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
      
      ntoread = nevs
      print *,'Number of MC reweighting events being read = ',ntoread
      
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
      include 'exptdefn2.f'
      include 'mcdefn.f'
      include 'fcnlocal.f'
      integer icomb
* fscale: defined by equation 26 in Poss-Sailer      
      double precision fscale
      double precision wsum(nbins),wsumsq(nbins),wsumtot
      integer nentries(nbins)
      double precision rwt
      double precision chisq
      double precision ndata,nmodel,var,chi
      double precision chisqi
      integer nbig
      integer myflag

* Copy fit parameters into more amenable variables
      ppeak = x(1)
      pbody = x(2)
      alphab = x(3)
      betab = x(4)
      alphaa = x(5)
      betaa = x(6)
      parm = 0.5d0*(1.0d0-ppeak-pbody)
      call cpbeta(bbody,alphab,betab,normb)
      call cpbeta(barms,alphaa,betaa,norma)

      fval = 0.0d0
      do j=1,nbins
         wsum(j) = 0.0d0
         wsumsq(j) = 0.0d0
         nentries(j) = 0
      enddo
      
* Loop over MC events
      do i=1,nmc
         myflag = 0
* Define variables
         ibin1 = iv(1,i)
         ibin2 = iv(2,i)
         rtype = iv(3,i)
         icomb = iv(4,i)
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
* Sum the reweighting weights for the (x1,x2) distribution bin
         if(myflag.eq.0)then
            wsum(icomb) = wsum(icomb) + rwt
            wsumsq(icomb) = wsumsq(icomb) + rwt*rwt
            nentries(icomb) = nentries(icomb) + 1
         endif
      enddo

* Sum ALL the weights
      wsumtot = 0.0d0
      do i=1,nbins
         wsumtot = wsumtot + wsum(i)
      enddo
      
      fscale = dble(ntotdata)/wsumtot

* Now construct the chi-squared
* FIXME. Should also include statistical uncertainty on model prediction

      chisq = 0.0d0
      do j=1,nbins
         ndata = dble(cont(j))
         nmodel  = fscale*wsum(j)
* Uncertainties
* 1. statistical uncertainty on finite data sample. 
*    Use model value for data statistics. => variance = nmodel
* 2. statistical uncertainty on the prediction
*    Use variance = fscale**2*wsumsq(j)
         var = nmodel + fscale*fscale*wsumsq(j)
         chi = (ndata - nmodel)/sqrt(var)
         chisq = chisq + chi*chi
      enddo
      print *,'f= ',chisq,
     +        ' Pars: ',x(1),x(2),x(3),x(4),x(5),x(6),' scale ',fscale
      fval = chisq
      
      if(iflag.eq.1)then
         print *,'Minimization Initialization '
         print *,'fval0 = ',fval
      endif

      if(iflag.eq.3)then
         print *,'Minimization Finished '
         print *,'fval = ',fval
         nbig = 0
         do j=1,nbins
            ndata = dble(cont(j))
            nmodel  = fscale*wsum(j)
            var = nmodel + fscale*fscale*wsumsq(j)
            chi = (ndata - nmodel)/sqrt(var)
            chisqi = chi*chi
            if(chisqi.gt.9.0d0)then
               print *,'Big chisq ',j,ndata,nmodel,chi,chisqi
               nbig = nbig + 1
            endif
         enddo
         print *,'No. of bins exceeding 3 sigma (expect 0.33%)',nbig
      endif
      
      end
      
      include 'cpbeta.f'
      
