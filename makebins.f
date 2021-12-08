      program makebins
*
* Read in file sorted by x1 of the electron beam 
* with 1,000,000 entries, and save every 10,000th such x1 value to a file.
* Thus these values can be used to define 100 bins in x1 with 
* approximately equal probability content under the 
* hypothesis used to generate the initial file
*
      implicit none
      
      integer iev,rtype
      integer nread
      integer ibin
      double precision y1,y2,x1,x2,z1,z2
      
* Read in sorted file      
      open(unit=31,file='testbc-2-x1ordered-1M.dat',status='old')
      
      nread = 0
      ibin = 0
      
 10   continue      
      read(31,*,end=999)iev,rtype,y1,y2,x1,x2,z1,z2
      nread = nread + 1
      if(mod(nread,10000).eq.0)then
         ibin = ibin + 1
         write(32,*)ibin,y1
      endif
      goto 10
      
  999 continue
  
      close(31)
      
      end
