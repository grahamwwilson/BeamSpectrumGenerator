      double precision function girceb (xmin, xmax, a, b, rng)          ADFI0846
      implicit none                                                     ADFI0847
      double precision xmin, xmax, a, b                                 ADFI0848
      external rng                                                      ADFI0849
      double precision t, p, u, umin, umax, x, w                        ADFI0850
      if ((a .gt. 1d0) .or. (b .lt. 1d0)) then                          ADFI0851
         girceb = -1d0                                                  ADFI0852
         print *,'beta-distribution expects a<=1<=b'
         return                                                         ADFI0854
      endif                                                             ADFI0855
      t = (1d0 - a) / (b + 1d0 - a)                                     ADFI0856
      p = b*t / (b*t + a * (1d0 - t)**b)                                ADFI0857
      if (xmin .le. 0d0) then                                           ADFI0858
         umin = 0d0                                                     ADFI0859
      elseif (xmin .lt. t) then                                         ADFI0860
         umin = p * (xmin/t)**a                                         ADFI0861
      elseif (xmin .eq. t) then                                         ADFI0862
         umin = p                                                       ADFI0863
      elseif (xmin .lt. 1d0) then                                       ADFI0864
         umin = 1d0 - (1d0 - p) * ((1d0 - xmin)/(1d0 - t))**b           ADFI0865
      else                                                              ADFI0866
         umin = 1d0                                                     ADFI0867
      endif                                                             ADFI0868
      if (xmax .ge. 1d0) then                                           ADFI0869
         umax = 1d0                                                     ADFI0870
      elseif (xmax .gt. t) then                                         ADFI0871
         umax = 1d0 - (1d0 - p) * ((1d0 - xmax)/(1d0 - t))**b           ADFI0872
      elseif (xmax .eq. t) then                                         ADFI0873
         umax = p                                                       ADFI0874
      elseif (xmax .gt. 0d0) then                                       ADFI0875
         umax = p * (xmax/t)**a                                         ADFI0876
      else                                                              ADFI0877
         umax = 0d0                                                     ADFI0878
      endif                                                             ADFI0879
      if (umax .lt. umin) then                                          ADFI0880
         girceb = -1d0                                                  ADFI0881
         return                                                         ADFI0882
      endif                                                             ADFI0883
 10   continue                                                          ADFI0884
      call rng (u)                                                      ADFI0885
      u = umin + (umax - umin) * u                                      ADFI0886
      if (u .le. p) then                                                ADFI0887
         x = t * (u/p)**(1d0/a)                                         ADFI0888
         w = (1d0 - x)**(b-1d0)                                         ADFI0889
      else                                                              ADFI0890
         x = 1d0 - (1d0 - t) * ((1d0 - u)/(1d0 - p))**(1d0/b)           ADFI0891
         w = (x/t)**(a-1d0)                                             ADFI0892
      endif                                                             ADFI0893
         call rng (u)                                                   ADFI0894
      if (w .le. u) goto 10                                             ADFI0895
      girceb = x                                                        ADFI0896
      end                                                               ADFI0897
