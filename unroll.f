      subroutine unroll(id1,id2)
      implicit none
      integer id1,id2
      real contents(400),errors(400)
      
* Unpack 2d histogram      
      call hunpak(id1,contents,'hist',1)
      call hunpke(id1,errors,'hist',1)      
      
* Fill into 1d histogram
      call hpak(id2,contents)
*      call hpake(id2,errors)
      
      end
