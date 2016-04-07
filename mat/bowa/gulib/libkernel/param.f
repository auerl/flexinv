        subroutine param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,
     &compgrid,nlatzomax,iadapt)
c
c this version is parameterized exactly as julias version
c and allows to distinguish between adaptive and regular
c parameterization, L.A. Aug. 2012
c
	implicit real*8 (a-h,o-z)
        integer nlatzones0
        integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)
        integer nsqrs0(nlatzomax),nsqtot0(nlatzomax+1)
	parameter(pi=3.1415926536)

c
c this is ludwigs preferred regular parameterization
c
        if(iadapt.eq.0)then
           print*,"came to ludwigs prefered parameterization"
           numto=0
           colat=-eq_incr/2.
           do k=1,nlatzones
              colat=colat+eq_incr
              theta=(colat/180.)*pi
              ! for this latitudinal zone, compute number of blocks (nsqrs)
              deltalon=eq_incr/(sin(theta))
              nsqrs(k)=(360./deltalon)+1
              ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
              if(iswit.eq.1)then
                 if(360./nsqrs(k).ge.compgrid)then
 100             if((mod(360./nsqrs(k),compgrid).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)+1
                    goto 100
                 else
                 endif
              elseif(360./nsqrs(k).lt.compgrid)then
 101             if((mod(compgrid,360./nsqrs(k)).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)-1
                    goto 101
                 else
                 endif
              endif
           else ! modified by lud
              if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
	   endif
         if(MOD(NSQRS(K),2).ne.0)then
              stop "nsqrs has to be even"
         endif
         ! print regular grid parameterization
         print*,k,nsqrs(k)
         nsqtot(k)=numto
         numto=numto+nsqrs(k)
         enddo
	 nsqtot(nlatzones+1)=numto         
c
c this is julias adaptive grid, which is slighl different from mine
c
         elseif(iadapt.eq.1)then
           print*,"came to julias prefered parameterization"
            n1layer0=0
            nlatzones0=180./compgrid
            colat=-compgrid/2.

            do k=1,nlatzones0
               colat=colat+compgrid
               theta=(colat/180.)*pi
               ! for this latitudinal zone, compute number of blocks (nsqrs)
               deltalon=compgrid/(sin(theta))
               nsqrs0(k)=(360./deltalon)+1
               if(mod(nsqrs0(k),2).ne.0)nsqrs0(k)=nsqrs0(k)-1
               ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
               if(iswit.eq.1)then
                  if(360./nsqrs0(k).ge.compgrid)then
 102              if(mod(360./nsqrs0(k),compgrid).ne.0)then
                     nsqrs0(k)=nsqrs0(k)+1
                     goto 102
                  else
                  endif
               elseif(360./nsqrs0(k).lt.compgrid)then
 103              if(mod(compgrid,360./nsqrs0(k)).ne.0)then
                     nsqrs0(k)=nsqrs0(k)-1
                     goto 103
                  else
                  endif
               endif
	       endif
               
               if(mod(nsqrs0(k),2).ne.0)stop "nsqrs has to be even"
               ! print rough parameterization acc. julia
               print*,k,nsqrs0(k)
               nsqtot0(k)=n1layer0
               n1layer0=n1layer0+nsqrs0(k)
            enddo

            ! define fine grid
            fact=compgrid/eq_incr
            n1layer=0
            print*,"ratio of coarser to finest grid is",fact
            do k=1,nlatzones
               k0=((k-1)/fact)+1
               nsqrs(k)=nsqrs0(k0)*fact
               nsqtot(k)=n1layer
               n1layer=n1layer+nsqrs(k)
               print*, "k=",k,"k0=",k0," nsqrs=", nsqrs(k),nsqrs0(k0)," lon=", 360./nsqrs(k)
               print*, "latitudinal zone =",k, "lon=", 360./nsqrs(k), "#pixels=", nsqrs(k) 
            enddo
            nsqtot(nlatzones+1)=n1layer
            print*,'Number of pixel with finest parameterization: n1layer=',n1layer
            numto=n1layer
            
         endif
       
         return
         end
