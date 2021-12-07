      subroutine rngm (r)
* RM48 based implementation (CERNLIB V116)
      implicit none
      integer len 
      parameter (len=10000)
      double precision rvec(len)
      double precision r
      include 'seeds.f'      
      integer ivalue
      data ivalue/-1/
      save ivalue
      save rvec
      
*      print *,'ivalue = ',ivalue
      
* RNG initialization
      if(ivalue.eq.-1)then
         call rm48in(seedm,0,0)
      endif
      
* RNG array filling/replenishment and array index reset      
      if(ivalue.eq.-1.or.ivalue.gt.len)then
         call rm48(rvec,len)
         ivalue = 1
      endif
      
      r = rvec(ivalue)     
*      print *,'r = ',r
      
* Set to next random number in array
      ivalue = ivalue + 1
      
      end
