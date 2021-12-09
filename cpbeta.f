      subroutine cpbeta(bpars,alpha,beta,norm)
      implicit none
      double precision bpars(4)
      double precision alpha,beta
      double precision norm
      external double precision dgamma
      
      bpars(1) = alpha
      bpars(2) = beta
      bpars(3) = alpha - 1.0d0
      bpars(4) = beta - 1.0d0
      
      norm = dgamma(alpha)*dgamma(beta)/dgamma(alpha+beta)
      
      end
