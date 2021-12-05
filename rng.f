      subroutine rng (r)
      implicit none
      double precision r

      call rngl(r)
     
      end
      
      include 'rngl.f'
