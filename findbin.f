      integer function findbin(x)
* Find which bin x belongs to based on x1bins.dat bin boundaries      
      implicit none
      double precision x
      integer ivalue
      integer nbins
      parameter (nbins=100)
      integer i,ibin
      double precision x1val(nbins)
      data ivalue/-1/
      save ivalue
      save x1val
      
* Initialize x1val array by reading file
      if(ivalue.eq.-1)then
         open(unit=21,file='x1bins-1M.dat',status='old')
         do i=1,nbins
            read(21,*)ibin,x1val(i)
         enddo
         close(21)
         ivalue = 1
      endif
      
      findbin = 0
      do i=1,nbins-1
         if(x.le.x1val(i))then
            findbin = i
            return
         endif
      enddo
      
      if(x.gt.x1val(nbins-1))then
         findbin = nbins
         return
      endif
      
      end
