      include 'gemodes.h'
      character*80 getunx,structurename
      parameter (twopi=6.2831853072)
      character*80 infile
c
      call chekcl('|     :r:1:Earth model|')
      structurename=getunx('',1,nbyts)
      infile=structurename(1:lnblnk(structurename))//'.bin'
      call openfl(1,infile,1,0,0,ierr,-1)
c
    1 continue
      write(6,"('type depth')") 
      read(5,*) depth
      call rewfl(1,istat)
      imodel=1
      call getgemodes(1,imodel,.true.,.false.,.false.,depth)
c
      nn=0
      ii=2
      fac0=1.-depth/6371.
      do ii=2,3
        if(ii.eq.2) open(1,file=structurename(1:lnblnk(structurename))//'.l')
        if(ii.eq.3) open(1,file=structurename(1:lnblnk(structurename))//'.r')
        do ll=0,1000
          call findmode(imodel,nn,ii,ll,imod)
          if(imod.gt.0) then
            pvel0=premgephvelo(ii-1,omegamod(imod,imodel))
            delc=100.*(pvelmod(imod,imodel)-pvel0)/pvel0
            write(1,"(i2,i4,6e12.4)") ii,ll,umod(imod,imodel),upmod(imod,imodel),uppmod(imod,imodel),
     #      vmod(imod,imodel),vpmod(imod,imodel),vppmod(imod,imodel)
            if(ll.eq.200) then
              if(ii.eq.2) then
              fac=haccmod(imod,imodel)
            write(6,"(i2,i4,6e12.4)") ii,ll,fac*umod(imod,imodel)/fac0,fac*upmod(imod,imodel),fac*uppmod(imod,imodel)*fac0,
     #      fac*vmod(imod,imodel)/fac0,fac*vpmod(imod,imodel),fac*vppmod(imod,imodel)*fac0
              else
              fac=1.
            write(6,"(i2,i4,6e12.4)") ii,ll,fac*umod(imod,imodel)/fac0,fac*upmod(imod,imodel),fac*uppmod(imod,imodel)*fac0,
     #      fac*vmod(imod,imodel)/fac0,fac*vpmod(imod,imodel),fac*vppmod(imod,imodel)*fac0
              endif
            endif
          endif
        enddo
        close(1)
      enddo
      go to 1
      end
