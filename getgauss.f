      subroutine getgauss(z1,z2)
* RNORMX based implementation using RANLUX
      implicit none
      integer len 
      parameter (len=2000)
      real rvec(len)
      double precision z1, z2
      integer ivalue
      include 'seeds.f'
      external ranlux
      data ivalue/-1/
      save ivalue
      save rvec
      
* RNG initialization
      if(ivalue.eq.-1)then
         call rluxgo(3,seedl,0,0)
      endif
      
* RNG array filling/replenishment and array index reset      
      if(ivalue.eq.-1.or.ivalue.gt.len)then
         call rnormx(rvec,len,ranlux)
         ivalue = 1
      endif
      
      z1 = dble(rvec(ivalue))
      z2 = dble(rvec(ivalue+1))
      
* Set to next Gaussian random number pair in array
      ivalue = ivalue + 2
      
      end
