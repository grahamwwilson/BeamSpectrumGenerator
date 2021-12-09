* Use two separate functions.
* Currently hard-coded to 100 bins for each.

      integer function findbinx1(x)
* Find which bin x belongs to based on x1bins.dat bin boundaries      
      implicit none
      double precision x
      integer ivalue
      integer nbins
      parameter (nbins=100)
      integer i,ibin
      double precision xval(nbins)
      data ivalue/-1/
      save ivalue
      save xval
      
* Initialize xval array by reading file
      if(ivalue.eq.-1)then
         open(unit=21,file='x1bins-100-1M-1-0.dat',status='old')
         do i=1,nbins
            read(21,*)ibin,xval(i)
         enddo
         close(21)
         ivalue = 1
      endif
      
      findbinx1 = 0
      do i=1,nbins-1
         if(x.le.xval(i))then
            findbinx1 = i
            return
         endif
      enddo
      
      if(x.gt.xval(nbins-1))then
         findbinx1 = nbins
         return
      endif
      
      end
      
      integer function findbinx2(x)
* Find which bin x belongs to based on x2bins.dat bin boundaries      
      implicit none
      double precision x
      integer ivalue
      integer nbins
      parameter (nbins=100)
      integer i,ibin
      double precision xval(nbins)
      data ivalue/-1/
      save ivalue
      save xval
      
* Initialize xval array by reading file
      if(ivalue.eq.-1)then
         open(unit=21,file='x2bins-100-1M-1-0.dat',status='old')
         do i=1,nbins
            read(21,*)ibin,xval(i)
         enddo
         close(21)
         ivalue = 1
      endif
      
      findbinx2 = 0
      do i=1,nbins-1
         if(x.le.xval(i))then
            findbinx2 = i
            return
         endif
      enddo
      
      if(x.gt.xval(nbins-1))then
         findbinx2 = nbins
         return
      endif
      
      end
