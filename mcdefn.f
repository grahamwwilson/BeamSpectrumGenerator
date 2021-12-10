      integer nmc
      parameter (nmc=10000000)
      integer nvars,mvars
      parameter (nvars = 7, mvars = 4)
      double precision v(nvars,nmc)
      integer iv(mvars,nmc)
      common/mc/v,iv
