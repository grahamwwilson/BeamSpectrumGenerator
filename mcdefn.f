      integer nmc
      parameter (nmc=1000000)
      integer nvars,mvars
      parameter (nvars = 7, mvars = 3)
      double precision v(nvars,nmc)
      integer iv(mvars,nmc)
      common/mc/v,iv
