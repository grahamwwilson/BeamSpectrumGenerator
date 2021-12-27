      program reframe
* Translate fit parameters into (alpha,beta) parameters used in generator
      implicit none
      
      double precision bpars(4)
      double precision mean,rms,norm
      
      print *,' ' 
      print *,'Specify beta mean in double precision'
      read(5,*)mean
*      print *,'Read mean value of ',mean

      print *,' ' 
      print *,'Specify beta rms in double precision'
      read(5,*)rms
*      print *,'Read rms value of ',rms
      
      call cpbetap(bpars,mean,rms,norm)
      
      print *,'Input values ',mean,rms
      print *,'beta ',bpars(1),bpars(2),bpars(3),bpars(4)
      
      end
      
      include 'cpbetap.f'
