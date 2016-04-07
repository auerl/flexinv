      subroutine initmode(lu,imodel,imod)
      common /plevel/ iprtlv
        common /sizeint/ intsize
      include 'gemodes.h'
      include 'gemodl.h'
c
      imodf=imodefil(imod,imodel)
      ipos=ibytfrst(imodel)+(imodf-1)*nbytesmo(imodel)
c
      if(nlevmodel(imodel).ne.nlevoutmodel(imodel)) then
        write(6,"('cannot get eigenf -- not all levels present')")
        write(6,"(2i10)") nlevmodel(imodel),nlevoutmodel(imodel)
        stop
      endif
      if(iprtlv.gt.4) then
        write(6,"('about to read eigenfunction info at:',2i10)") lu,ipos
      endif
c
c---- put the eigenfunction information into the common block
c
      ishif=33+intsize*3!lb twb 2009
      do ilin=1,nlevmodel(imodel)
c        ibyt=ipos+45+(ilin-1)*24
        ibyt=ipos+ishif+(ilin-1)*24
c        write(6,"('ilin, ibyt:',2i10)") ilin,ibyt
        call bffi(lu,1,f1,4,iostat,nr,ibyt)
        call bffi(lu,1,f2,4,iostat,nr,ibyt+4)
        call bffi(lu,1,f3,4,iostat,nr,ibyt+8)
        call bffi(lu,1,f4,4,iostat,nr,ibyt+12)
        call bffi(lu,1,f5,4,iostat,nr,ibyt+16)
        call bffi(lu,1,f6,4,iostat,nr,ibyt+20)
        u(ilin)=f1
        up(ilin)=f2
        v(ilin)=f3
        vp(ilin)=f4
        p(ilin)=f5
        pp(ilin)=f6
      enddo
      nord=  iovermod(imod,imodel)
      jcom=  itypemod(imod,imodel)
      lord=  ilordmod(imod,imodel)
      wcom=  omegamod(imod,imodel)
      qbar=  smallqmod(imod,imodel)
      cgp=   gvelmod(imod,imodel)
      avert= vaccmod(imod,imodel)
      ahor=  haccmod(imod,imodel)
      phis=  potmod(imod,imodel)
      if(iprtlv.gt.2) then
        write(6,"(3i6,6e12.3)") nord,jcom,lord,wcom,qbar,cgp,avert,ahor,phis
      endif
      if(nord.eq.0.and.lord.gt.3.and.lord.lt.7) then
        if(iprtlv.gt.2) write(6,"(i4,6f10.5)")
     #         jcom,u(moho),up(moho),v(moho),vp(moho),p(moho),pp(moho)
      endif
      return
      end
      
      
