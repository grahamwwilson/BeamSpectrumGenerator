* Note. Likely not a problem if the array is much bigger than needed.
      integer nmcarray
      parameter (nmcarray=10000000)
      integer nmc
      integer nvars,mvars
      parameter (nvars = 7, mvars = 4)
      double precision v(nvars,nmcarray)
      integer iv(mvars,nmcarray)
      common/mc/v,iv,nmc
