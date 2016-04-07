      include 'gemodes.h'
      character*80 getunx,structurename
      parameter (twopi=6.2831853072)
      character*80 infile
      common /plevel/ iprtlv
c
      call chekcl('|     :r:1:Earth model '//
     #            '|   -v:o:1:Verbosity level [0]|')
      read(getunx('-v',1,nbyts),*) iprtlv
      write(6,"('type name of file to read:')") 
      structurename=getunx('',1,nbyts)
      infile=structurename(1:lnblnk(structurename))//'.bin'
c      read(5,"(a)") infile
      call openfl(1,infile,1,0,0,ierr,-1)
      write(6,"('type depth')") 
ccc      read(5,*) depth
      depth=1.
      call rewfl(1,istat)
      imodel=1
      call getgemodes(1,imodel,.false.,.false.,.false.,depth)
      write(6,"('type 1 for phase velocity')")
      write(6,"('type 2 for group velocity')")
      write(6,"('type 3 for phase vel. difference from PREM')")
ccc      read(5,*) isel
      isel=3
c
      nn=0
      ii=2
      do ii=2,3
        if(ii.eq.2) open(1,file=structurename(1:lnblnk(structurename))//'.love')
        if(ii.eq.3) open(1,file=structurename(1:lnblnk(structurename))//'.rayl')
        do ll=0,1000
          call findmode(imodel,nn,ii,ll,imod)
          if(imod.gt.0) then
            pvel0=premgephvelo(ii-1,omegamod(imod,imodel))
            delc=100.*(pvelmod(imod,imodel)-pvel0)/pvel0
            if(isel.eq.3) then
              write(1,"(2f8.3,3i5,4f15.8)") 1000.*omegamod(imod,imodel)/twopi,delc,
     #           nn,ii,ll,omegamod(imod,imodel),pvelmod(imod,imodel),gvelmod(imod,imodel),smallqmod(imod,imodel)
            else if(isel.eq.2) then
              write(1,"(2f8.3,3i5,2f15.8)") 1000.*omegamod(imod,imodel)/twopi,
     #           gvelmod(imod,imodel),nn,ii,ll,omegamod(imod,imodel),pvelmod(imod,imodel)
            else if(isel.eq.1) then
              write(1,"(2f8.3,3i5,2f15.8)") 1000.*omegamod(imod,imodel)/twopi,
     #           pvelmod(imod,imodel),nn,ii,ll,omegamod(imod,imodel),pvelmod(imod,imodel)
            endif
          endif
        enddo
        close(1)
      enddo
      end
