      program makebins
*
* Read in files with 1M events, sorted either by x1 or x2.
* Save every 10,000th such value to a file.
* Thus these values can be used to define 100 bins in x1 
* and 100 bins in x2 with approximately equal probability content 
* under the hypothesis used to generate the initial file.
*
* The beam energy spread is usually different between x1 and x2 - so the 
* bin edges should be systematically different.
*
      implicit none

      open(unit=31,file='testbc-1-0-1M.x1sorted',status='old')
      call doit(31,51,1)
      
      open(unit=32,file='testbc-1-0-1M.x2sorted',status='old')      
      call doit(32,52,2)
      
      end 
      
      subroutine doit(lin,lout,iwhich)
      implicit none
      integer lin,lout,iwhich
      
      integer iev,rtype
      integer nread
      integer ibin
      double precision x1,x2,y1,y2,z1,z2,weight,weightp,rwt      
      
      nread = 0
      ibin = 0
      
 10   continue      
      read(lin,*,end=999)iev,rtype,x1,x2,y1,y2,z1,z2,weight,weightp,rwt
      nread = nread + 1
      if(mod(nread,10000).eq.0)then
         ibin = ibin + 1
         if(iwhich.eq.1)then
            write(lout,*)ibin,x1
         elseif(iwhich.eq.2)then
            write(lout,*)ibin,x2         
         endif
      endif
      goto 10
      
  999 continue
      close(lin)
      
      end
      

      
      
