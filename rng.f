      subroutine rng (r)
      implicit none
      double precision r

*      call rngl(r)
      call rngm(r)
     
      end
      
      include 'rngl.f'
      include 'rngm.f'
