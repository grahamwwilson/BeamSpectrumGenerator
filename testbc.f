      program testbc
* V1. test
* V2. testb.  Add BES.
* V3. testbc. Separate into peak, arm, and body with specified probabilities.
* V4. testbc. But add separate slope parameters for arm and body.
      implicit none
      double precision x1,x2
      external rng
      real rtype
      double precision a1(0:7)
      double precision u,girceb
      double precision x1m,x2m
      double precision rg1,rg2
      double precision z1,z2
      double precision y1,y2
      double precision s1,s2
      double precision pnorm(2)
      parameter (s1=0.190d-2, s2=0.152d-2)
      logical lhbook
      parameter (lhbook=.true.)
      integer iev
      include 'hinit.f'
* peak, body
      data pnorm/0.26307d0,0.28151d0/
      data a1 /
     $   0.49808d+00,  0.54613d+00,  0.12287d+02, -0.62756d+00, 
     $   0.42817d+00, -0.69120d+00,  0.17067d+02,  0.51143d+00 /       
      
      print *,'Beam energy spread (BES) parameters'
      print *,'s1 (electron) = ',s1
      print *,'s2 (positron) = ',s2
      
      print *,'Partitioning of events'
      print *,'p(peak) = ',pnorm(1)
      print *,'p(body) = ',pnorm(2)
      print *,'p(arms) = ',1.0d0-pnorm(1)-pnorm(2),
     +                     0.5d0*(1.0d0-pnorm(1)-pnorm(2))
     
      print *,'Slope parameters'
      print *,a1(2),a1(3)
      
      if(lhbook)then
         call hlimit(nwpawc)
         call bookh
      endif
      
      x1m=0d0
      x2m=0d0
      
* Generate two random numbers according to the CIRCE functions 
* for scaled beam energy after potential beamstrahlung.
* Currently only need the three parameters a1(0), a1(2), a1(3) 
* in this minimalist implementation   
      
      do iev=1,10000000
      
* Get two normally distributed (standardized random numbers)
         call getgauss(rg1,rg2)
* Scaled beam energy after uncorrelated Gaussian BES
         z1 = 1.0d0 + s1*rg1
         z2 = 1.0d0 + s2*rg2
      
* Beamstrahlung part      
         call rng(u)
         if (u .le. pnorm(1)) then 
* peak
            x1 = 1d0
            x2 = 1d0
            rtype = 1.0                                            
         elseif (u .le. pnorm(1)+pnorm(2))then                                                        
* body         
            x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
            x2 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
            rtype = 2.0
         elseif (u. le. 0.5d0*(1d0+pnorm(1)+pnorm(2)))then
* arm1      
            x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)            
            x2 = 1d0
            rtype = 3.0
         else
* arm2
            x1 = 1d0
            x2 = 1d0 - girceb (0d0, 1d0-x2m, a1(3)+1d0, a1(2)+1d0, rng)            
            rtype = 4.0
         endif
     
* Scaled energy after BES and beamstrahlung
         y1 = z1*x1
         y2 = z2*x2

         if(lhbook)then
            call hfill(101,real(y1),0.0,1.0)
            call hfill(102,real(y2),0.0,1.0)
            call hfill(103,real(sqrt(y1*y2)),0.0,1.0)
            call hfill(104,real(abs(y1-y2)),0.0,1.0)
            call hfill(105,rtype,0.0,1.0)
            call hfill(106,real(y1-y2),0.0,1.0)            
         endif
 
      enddo
      
      if(lhbook)then
         call hrput(0,'testbc.hbook','NT')
      endif
      
      end
      subroutine bookh
      implicit none
      call hidopt(0,'stat')
      call hbook1(101,'x1 ',1101,-0.0005,1.1005,0.0)
      call hbook1(102,'x2 ',1101,-0.0005,1.1005,0.0)
      call hbook1(103,'sqrt(x1*x2)',1101,-0.0005,1.1005,0.0)      
      call hbook1(104,'|x1-x2|',1100,0.000,1.100,0.0)
      call hbook1(105,'Type ',4,0.5,4.5,0.0)
      call hbook1(106,'x1-x2',2200,-1.100,1.100,0.0)      

      end
      
      include 'rng.f'
      include 'girceb.f'
      include 'getgauss.f'
