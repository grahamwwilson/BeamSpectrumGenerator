      program test
      implicit none
      double precision x1,x2
      external rng
      double precision a1(0:7)
      double precision u,girceb
      double precision x1m,x2m
      integer iev
      include 'hinit.f'
      data a1 /
     $   0.49808e+00,  0.54613e+00,  0.12287e+02, -0.62756e+00, 
     $   0.42817e+00, -0.69120e+00,  0.17067e+02,  0.51143e+00 /       
      
      call hlimit(nwpawc)
      
      call bookh
      
      x1m=0d0
      x2m=0d0
      
* Generate two random numbers according to the CIRCE functions 
* Currently only need the three parameters a1(0), a1(2), a1(3) 
* in this minimalist implementation   
      
      do iev=1,1000000
      
         call rng(u)
         if (u .le. a1(0)) then 
            x1 = 1d0                                                    
         else                                                        
            x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
         endif
         
         call rng(u)
         if (u .le. a1(0)) then 
            x2 = 1d0                                                    
         else                                                        
            x2 = 1d0 - girceb (0d0, 1d0-x2m, a1(3)+1d0, a1(2)+1d0, rng)
         endif
         
*         print *,iev,' x1 = ',x1,' x2 = ',x2 

         call hfill(101,real(x1),0.0,1.0)
         call hfill(102,real(x2),0.0,1.0)
         call hfill(103,real(sqrt(x1*x2)),0.0,1.0)
         call hfill(104,real(abs(x1-x2)),0.0,1.0)
 
      enddo
      
      call hrput(0,'test.hbook','NT')
      
      end
      subroutine bookh
      implicit none
      call hidopt(0,'stat')
      call hbook1(101,'x1 ',1001,-0.0005,1.0005,0.0)
      call hbook1(102,'x2 ',1001,-0.0005,1.0005,0.0)
      call hbook1(103,'sqrt(x1*x2)',1001,-0.0005,1.0005,0.0)      
      call hbook1(104,'|x1-x2|',1001,-0.0005,1.0005,0.0)

      end
      
      include 'rng.f'
      include 'girceb.f'
