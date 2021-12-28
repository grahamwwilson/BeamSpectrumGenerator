      subroutine readheader(lun,ntoread)
      implicit none
      integer lun
      integer nheader,version,nparameters,nevs,seedm,seedl
      integer ntoread
      double precision mu1,mu2,s1,s2,pnorm1,pnorm2,alphab,betab
      double precision alphaa,betaa
      
      read(lun,*)nheader
      read(lun,*)version
      read(lun,*)nparameters
      read(lun,*)nevs
      read(lun,*)seedm
      read(lun,*)seedl
      read(lun,*)mu1
      read(lun,*)mu2
      read(lun,*)s1
      read(lun,*)s2
      read(lun,*)pnorm1
      read(lun,*)pnorm2
      read(lun,*)alphab
      read(lun,*)betab
      read(lun,*)alphaa
      read(lun,*)betaa
      
      ntoread = nevs
      print *,'Number of MC reweighting events being read = ',ntoread
      
      end
