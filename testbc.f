      program testbc
* V1. test
* V2. testb.  Add BES.
* V3. testbc. Separate into peak, arm, and body with specified probabilities.
* V4. testbc. Add separate slope parameters for arm and body.
* V5. testbc  Add output file with header portion.
      implicit none
      double precision x1,x2
      external rng
      integer nevs
      parameter (nevs=1000000)
*      parameter (nevs=187100)
      integer nheader
      parameter (nheader=16)
      integer version
      parameter (version=7)
      integer nparameters
      parameter (nparameters=10)
      integer rtype
      integer findbin
      double precision betabody(4)
      double precision betaarms(4)
      double precision betabodyp(4)
      double precision betaarmsp(4)      
      double precision u,girceb
      double precision y1m,y2m
      double precision rg1,rg2
      double precision z1,z2
      double precision y1,y2
      double precision mu1,mu2
      parameter (mu1=1.0d0, mu2=1.0d0)
      double precision s1,s2
      double precision bnormbody,bnormarms
      double precision pnorm(2)
      
      double precision bnormbodyp,bnormarmsp
      double precision pnormp(2)      
      
      double precision pregion,py1,py2
      double precision weight,weight0
      double precision pregionp,py1p,py2p
      double precision weightp,weightb
      double precision rwt
      double precision dy1,dy2
      integer ibin   
      integer flag   
      external double precision dgamma
      parameter (s1=0.190d-2, s2=0.152d-2)
      integer ia,ib,ja,jb
      parameter (ia=1, ib=2,ja=3,jb=4)
      logical lhbook
      parameter (lhbook=.true.)
      include 'seeds.f'
      integer iev
      include 'hinit.f'
* Note CIRCE Beta distribution is  a1*x^a2*(1-x)^a3 (no longer used)
*      double precision a1(0:7)
*      data a1 /
*     $   0.49808d+00,  0.54613d+00,  0.17161d+02, -0.65863d+00, 
*     $   0.42817d+00, -0.69120d+00,  0.13707d+02, -0.63926d+00 /         

*     Beamstrahlung parameters from include file
      include 'betapars.f'
      include 'betaalt.f'

      print *,'Beam energy spread (BES) parameters'
      print *,'s1 (electron) = ',s1
      print *,'s2 (positron) = ',s2
      
      print *,'Partitioning of events'
      print *,'p(peak) = ',pnorm(1)
      print *,'p(body) = ',pnorm(2)
      print *,'p(arms) = ',1.0d0-pnorm(1)-pnorm(2),
     +                     0.5d0*(1.0d0-pnorm(1)-pnorm(2))
     
      print *,'Slope parameters'
*      print *,'Arms: ',a1(2),a1(3)
*      print *,'Body: ',a1(6),a1(7)
      
      print *,'Body (alpha,beta): ',betabody(ia),betabody(ib)      
      print *,'Arms (alpha,beta): ',betaarms(ia),betaarms(ib)

      bnormbody = dgamma(betabody(ia))*dgamma(betabody(ib))/
     +            dgamma(betabody(ia)+betabody(ib))
      bnormarms = dgamma(betaarms(ia))*dgamma(betaarms(ib))/
     +            dgamma(betaarms(ia)+betaarms(ib))
      print *,'B(alp,beta) = ',bnormbody,bnormarms    

* Alternative functions
      print *,'Bodyp (alpha,beta): ',betabodyp(ia),betabodyp(ib)      
      print *,'Armsp (alpha,beta): ',betaarmsp(ia),betaarmsp(ib)

      bnormbodyp = dgamma(betabodyp(ia))*dgamma(betabodyp(ib))/
     +            dgamma(betabodyp(ia)+betabodyp(ib))
      bnormarmsp = dgamma(betaarmsp(ia))*dgamma(betaarmsp(ib))/
     +            dgamma(betaarmsp(ia)+betaarmsp(ib))
      print *,'B(alp,beta) p= ',bnormbodyp,bnormarmsp    

      if(lhbook)then
         call hlimit(nwpawc)
         call bookh
      endif
      
      y1m=0d0
      y2m=0d0
      
* Generate two random variates (x1,x2) for Ei/Enominal, following a 
* probability density function p(x1,x2).
* The pdf is made up of four regions (the region choice is 
* governed by the multinomial distribution), 
* and in each the random variate x = y*z
* where z reflects the beam energy spread (present in each case), and 
* where y reflects the potential beamstrahlung.
* z ~ Ga(mean,sigma)
* y ~ 1 - B(alpha,beta)
* Use the two slope parameters (alpha, beta) for the arms beta distribution 
* and the two slope parameters (alpha, beta) for the body beta distribution
      
      write(45,*)nheader
      write(45,*)version
      write(45,*)nparameters
      write(45,*)nevs
      write(45,*)seedm
      write(45,*)seedl
      write(45,*)mu1
      write(45,*)mu2
      write(45,*)s1
      write(45,*)s2
      write(45,*)pnorm(1)
      write(45,*)pnorm(2)
      write(45,*)betabody(ia)
      write(45,*)betabody(ib)
      write(45,*)betaarms(ia)
      write(45,*)betaarms(ib)
      
      do iev=1,nevs
      
* Get two normally distributed (standardized random numbers)
         call getgauss(rg1,rg2)
* Scaled beam energy after uncorrelated Gaussian BES
         z1 = mu1 + s1*rg1
         z2 = mu2 + s2*rg2
      
* Beamstrahlung part      
         call rng(u)
         if (u .le. pnorm(1)) then 
* peak
            rtype = 1
            pregion = pnorm(1)
            y1 = 1d0
            y2 = 1d0
            weight0 = 1d0
            
            pregionp = pnormp(1)
            weightb = 1d0              
         elseif(u .le.pnorm(1)+pnorm(2))then                                                       
* body         
            rtype = 2
            pregion = pnorm(2)
            y1 = 1d0-girceb(0d0,1d0-y1m,betabody(ia),betabody(ib),rng)
            py1=((1d0-y1)**betabody(ja))*(y1**betabody(jb))/bnormbody
            y2 = 1d0-girceb(0d0,1d0-y1m,betabody(ia),betabody(ib),rng)
            py2=((1d0-y2)**betabody(ja))*(y2**betabody(jb))/bnormbody
            weight0 = py1*py2
            
            pregionp = pnormp(2)
        py1p=((1d0-y1)**betabodyp(ja))*(y1**betabodyp(jb))/bnormbodyp
        py2p=((1d0-y2)**betabodyp(ja))*(y2**betabodyp(jb))/bnormbodyp                   
            weightb = py1p*py2p
            
         elseif (u. le. 0.5d0*(1d0+pnorm(1)+pnorm(2)))then
* arm1      
            rtype = 3
            pregion = 0.5d0*(1.0d0-pnorm(1)-pnorm(2))
            y1 = 1d0-girceb(0d0,1d0-y1m,betaarms(ia),betaarms(ib),rng)
            py1=((1d0-y1)**betaarms(ja))*(y1**betaarms(jb))/bnormarms             
            y2 = 1d0
            weight0 = py1
            
            pregionp = 0.5d0*(1.0d0-pnormp(1)-pnormp(2))
        py1p=((1d0-y1)**betaarmsp(ja))*(y1**betaarmsp(jb))/bnormarmsp
            weightb = py1p           
         else
* arm2
            rtype = 4
            pregion = 0.5d0*(1.0d0-pnorm(1)-pnorm(2))            
            y1 = 1d0
            y2 = 1d0-girceb(0d0,1d0-y2m,betaarms(ia),betaarms(ib),rng)
            py2=((1d0-y2)**betaarms(ja))*(y2**betaarms(jb))/bnormarms
            weight0 = py2
            
            pregionp = 0.5d0*(1.0d0-pnormp(1)-pnormp(2))
        py2p=((1d0-y2)**betaarmsp(ja))*(y2**betaarmsp(jb))/bnormarmsp
            weightb = py2p                        
         endif
         
* Weight factor for initial parametrization
         weight = pregion*weight0
* Weight factor for new parametrization
         weightp = pregionp*weightb
         
         flag = 0
         
* Check for infinities
         if(weight-1d0.eq.weight)then
            print *,'Infinity for weight  in event ',iev
            flag = flag + 1
         endif
         if(weightp-1d0.eq.weightp)then
            print *,'Infinity for weightp in event ',iev
            flag = flag + 2
         endif                     
         
         if(flag.eq.0)then
            rwt = weightp/weight
         else
            rwt = pregionp/pregion
            print *,'Setting rwt to ',rwt,             
     +              ' in event ',iev,' region ',rtype   
         endif
         
* Scaled energy after BES and beamstrahlung
         x1 = y1*z1
         x2 = y2*z2
         
* Save information for each event to file
         write(45,*)iev,rtype,x1,x2,y1,y2,z1,z2,weight,weightp,rwt

         if(lhbook)then
         
            ibin = findbin(x1)
            
            call hfill(101,real(x1),0.0,1.0)
            call hfill(102,real(x2),0.0,1.0)
            call hfill(103,real(sqrt(x1*x2)),0.0,1.0)
            call hfill(104,real(abs(x1-x2)),0.0,1.0)
            call hfill(105,real(rtype),0.0,1.0)
            call hfill(106,real(x1-x2),0.0,1.0)
            call hfill(107,real(ibin),0.0,1.0)
            call hfill(200+10*rtype+1,real(x1),0.0,1.0)
            call hfill(200+10*rtype+2,real(x2),0.0,1.0)  
            
            call hfill(1101,real(x1),0.0,real(rwt))
            call hfill(1102,real(x2),0.0,real(rwt))
            call hfill(1103,real(sqrt(x1*x2)),0.0,real(rwt))
            call hfill(1104,real(abs(x1-x2)),0.0,real(rwt))
            call hfill(1105,real(rtype),0.0,real(rwt))
            call hfill(1106,real(x1-x2),0.0,real(rwt))
            call hfill(1107,real(ibin),0.0,real(rwt))
            call hfill(1200+10*rtype+1,real(x1),0.0,real(rwt))
            call hfill(1200+10*rtype+2,real(x2),0.0,real(rwt))            
                      
         endif
 
      enddo
      
      close(45)
      
      if(lhbook)then
*         call histdo
         call hrput(0,'testbc.hbook','NT')
         call hprint(107)
         call hprint(1107)
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
      call hbook1(107,'Equiprobability x1 bin',100,0.5,100.5,0.0)
      call hbook1(211,'x1 R1',1101,-0.0005,1.1005,0.0)
      call hbook1(212,'x2 R1',1101,-0.0005,1.1005,0.0)
      call hbook1(221,'x1 R2',1101,-0.0005,1.1005,0.0)
      call hbook1(222,'x2 R2',1101,-0.0005,1.1005,0.0)
      call hbook1(231,'x1 R3',1101,-0.0005,1.1005,0.0)
      call hbook1(232,'x2 R3',1101,-0.0005,1.1005,0.0)
      call hbook1(241,'x1 R4',1101,-0.0005,1.1005,0.0)
      call hbook1(242,'x2 R4',1101,-0.0005,1.1005,0.0)
      
* Weighted distributions
      call hbook1(1101,'x1 ',1101,-0.0005,1.1005,0.0)
      call hbook1(1102,'x2 ',1101,-0.0005,1.1005,0.0)
      call hbook1(1103,'sqrt(x1*x2)',1101,-0.0005,1.1005,0.0)      
      call hbook1(1104,'|x1-x2|',1100,0.000,1.100,0.0)
      call hbook1(1105,'Type ',4,0.5,4.5,0.0)
      call hbook1(1106,'x1-x2',2200,-1.100,1.100,0.0)
      call hbook1(1107,'Equiprobability x1 bin',100,0.5,100.5,0.0)
      call hbook1(1211,'x1 R1',1101,-0.0005,1.1005,0.0)
      call hbook1(1212,'x2 R1',1101,-0.0005,1.1005,0.0)
      call hbook1(1221,'x1 R2',1101,-0.0005,1.1005,0.0)
      call hbook1(1222,'x2 R2',1101,-0.0005,1.1005,0.0)
      call hbook1(1231,'x1 R3',1101,-0.0005,1.1005,0.0)
      call hbook1(1232,'x2 R3',1101,-0.0005,1.1005,0.0)
      call hbook1(1241,'x1 R4',1101,-0.0005,1.1005,0.0)
      call hbook1(1242,'x2 R4',1101,-0.0005,1.1005,0.0)      
           
      end
      
      include 'rng.f'
      include 'girceb.f'
      include 'getgauss.f'
      include 'findbin.f'
