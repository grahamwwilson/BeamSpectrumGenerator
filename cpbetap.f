      subroutine cpbetap(bpars,mean,rms,norm)
* We are fitting here using the mean 
* and rms parametrization instead of alpha and beta.
* See https://en.wikipedia.org/wiki/Beta_distribution#Alternative_parameterizations
*
* But programmatically we just need to calculate alpha and beta from the 
* mean and rms quantities.
* 
      implicit none
      double precision bpars(4)
      double precision mean,rms,var,nu
      double precision alpha,beta
      double precision norm
      external double precision dgamma
      
      nu = (mean*(1.0d0-mean)/(rms*rms)) - 1.0d0
      alpha = mean*nu
      beta = nu-alpha
      
      bpars(1) = alpha
      bpars(2) = beta
      bpars(3) = alpha - 1.0d0
      bpars(4) = beta - 1.0d0
      
      norm = dgamma(alpha)*dgamma(beta)/dgamma(alpha+beta)
      
      end
