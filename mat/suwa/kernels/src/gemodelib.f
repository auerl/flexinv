      subroutine getgemodes(lu,imodel,fundonly,overonly,acconly,depth)
      include 'gemodes.h'
c        parameter(intsize=8)! twb lb 2009. set to 4 on a 32bit machine
        common /sizeint/ intsize
      common /plevel/ iprtlv
c
      logical fundonly
      logical overonly
      logical acconly
      character*80 modname
      character*80 title
      integer ifanis
      integer ideck
      real tref
      integer nlev
      integer nlevout
      integer nic
      integer noc
      real vpv(mxlevels),vph(mxlevels),vsv(mxlevels),vsh(mxlevels),eta(mxlevels)
      real rad(mxlevels),rho(mxlevels),qmu(mxlevels),qka(mxlevels)
      integer ilevflag(mxlevels)
      dimension modetype(mxmodes)
      dimension lord(mxmodes)
      dimension nord(mxmodes)
      integer ifmodel(mxmodels)
      data ifmodel /mxmodels*0/

        intsize=8
cTEST
c        print*,"in getgemodes"
c        print*,"depth=",depth
      if(imodel.gt.mxmodels) then
        stop 'trying to read too many models'
      else
        ifmodel(imodel)=1
      endif
      nmodels=0
c
c---- read in all header stuff
c
      call bffi(lu,1,modname,80,ierr,nr,0)
      call bffi(lu,1,depmin,4,ierr,nr,0)
      call bffi(lu,1,depmax,4,ierr,nr,0)
      call bffi(lu,1,title,80,ierr,nr,0)
      call bffi(lu,1,ifanis,4,ierr,nr,0)
      call bffi(lu,1,tref,4,ierr,nr,0)
      call bffi(lu,1,ideck,4,ierr,nr,0)
      call bffi(lu,1,nlev,4,ierr,nr,0)
      call bffi(lu,1,nic,4,ierr,nr,0)
      call bffi(lu,1,noc,4,ierr,nr,0)
      call bffi(lu,1,nmoho,4,ierr,nr,0)
cTEST
c        print*,"model :"
c        print*,modname,depmin,depmax,nlev
c        pause

      do ilev=1,nlev
        call bffi(lu,1,rad(ilev),4,iostat,nr,0)
        call bffi(lu,1,rho(ilev),4,iostat,nr,0)
        call bffi(lu,1,vpv(ilev),4,iostat,nr,0)
        call bffi(lu,1,vph(ilev),4,iostat,nr,0)
        call bffi(lu,1,vsv(ilev),4,iostat,nr,0)
        call bffi(lu,1,vsh(ilev),4,iostat,nr,0)
        call bffi(lu,1,eta(ilev),4,iostat,nr,0)
        call bffi(lu,1,qmu(ilev),4,iostat,nr,0)
        call bffi(lu,1,qka(ilev),4,iostat,nr,0)
        call bffi(lu,1,ilevflag(ilev),4,iostat,nr,0)
      enddo
      call bffi(lu,1,nlevout,4,ierr,nr,0)
      call bffi(lu,1,lmin,4,ierr,nr,0)
      call bffi(lu,1,lmax,4,ierr,nr,0)
      call bffi(lu,1,nbran,4,ierr,nr,0)
      call bffi(lu,1,wmin,4,ierr,nr,0)
      call bffi(lu,1,wmax,4,ierr,nr,0)
      call bffi(lu,1,eps,4,ierr,nr,0)
      call bffi(lu,1,wgrav,4,ierr,nr,0)
cTEST
c        print*,"checking  model parameters"
c        print*,nlevout,lmin,lmax,nbran,wmin,wmax,eps,wgrav
c        pause
c
c---- put it in the first common block
c
      do i=1,mxmodels
        nmodels=nmodels+ifmodel(i)
      enddo
      nlevmodel(imodel)=nlev
      nlevoutmodel(imodel)=nlevout
      ninnercore(imodel)=nic
      noutercore(imodel)=noc
      nmohorovic(imodel)=nmoho
      ifanismodel(imodel)=ifanis
      ifdeckmodel(imodel)=ideck
      trefmodel(imodel)=tref
      lminmodel(imodel)=lmin
      lmaxmodel(imodel)=lmax
      nbranmodel(imodel)=nbran
      fminmodel(imodel)=wmin
      fmaxmodel(imodel)=wmax
      epsmodel(imodel)=eps
      wgravmodel(imodel)=wgrav
      do ilev=1,nlev
        vpvlayer(ilev,imodel)=vpv(ilev)
        vphlayer(ilev,imodel)=vph(ilev)
        vsvlayer(ilev,imodel)=vsv(ilev)
        vshlayer(ilev,imodel)=vsh(ilev)
        etalayer(ilev,imodel)=eta(ilev)
        radlayer(ilev,imodel)=rad(ilev)
        rholayer(ilev,imodel)=rho(ilev)
        qmulayer(ilev,imodel)=qmu(ilev)
        qkalayer(ilev,imodel)=qka(ilev)
        ioutlayer(ilev,imodel)=ilevflag(ilev)
      enddo
      modelname(imodel)=modname
      if(iprtlv.gt.1) then
        write(6,"(a80)") modname
        write(6,"(a80)") title
        write(6,"('ifanis ideck  nlev nlevo   nic   noc nmoho      
     &tref')")
        write(6,"(7i6,f10.2)") ifanis,ideck,nlev,nlevout,nic,noc,nmoho,tref
        write(6,"(' lmin lmax nbra           wmin           wmax       
     &     eps          wgrav')")
        write(6,"(3i5,4e15.5)") lmin,lmax,nbran,wmin,wmax,eps,wgrav
      endif
c
      do ilev=1,nlev
        if(ilevflag(ilev).eq.1) then
          if(iprtlv.gt.2) then
            write(6,"(i4,i2,8f10.3)") ilev,ilevflag(ilev),rad(ilev),
     #      vpv(ilev),vph(ilev),vsv(ilev),vsh(ilev),eta(ilev),qmu(ilev),
     #      qka(ilev)
          endif
        endif
      enddo
c
c---- bracket the source
c
cTEST
        print*,"nlev=",nlev,"rad(nlev)=",rad(nlev)
      radius=rad(nlev)-depth
      ir1=0
      ir2=0
      do i=1,nlev
        if(rad(i).ge.radius.and.ir2.eq.0) then
          ir2=i
          ir1=i-1
        endif
      enddo
      if(ir2.gt.0) then
        hn=(rad(ir2)-rad(ir1))/6371.
        dr=(radius-rad(ir1))/6371.
        hn2=1./(hn**2)
        hn3=1./(hn**3)
        if(iprtlv.gt.1) then
          write(6,"('radii -- below, source, above')")
          write(6,"(3f10.3)") rad(ir1),radius,rad(ir2)
          write(6,"('  nlev nlevo   ir1   ir2')")
          write(6,"(4i6)") nlev,nlevout,ir1,ir2
          write(6,"('hn,dr,hn2,hn3',4e12.3)") hn,dr,hn2,hn3
        endif
      else
        print*,"radius=",radius
        print*,"depth=",depth
        write(6,"('source not bracketed:',g15.5)")radius
        stop
      endif
      if(ilevflag(ir1).gt.0.and.ilevflag(ir2).gt.0) then
          ilin=0
          do i=1,ir1
            if(ilevflag(i).gt.0) ilin=ilin+1
          enddo
          radiusmod(imodel)=1.0-depth/6371.
          depthmod(imodel)=depth
          if(iprtlv.gt.1) then
            write(6,"('level below is ilin:',i4)") ilin
            write(6,"('level above is ilin+1:',i4)") ilin+1
            write(6,"('depth and normalized radius:',f8.2,f8.5)") 
     #               depthmod(imodel),radiusmod(imodel)
          endif
      else
        write(6,"('source level not in catalog:',3g15.5)") depth,depmin,depmax
        stop
      endif
c
c---- scan the modes
c
      call getflpos(lu,ipos0,istat)
      nbytpmod=(11+6*nlevout)*4
c      nbytpmod=(11+6*nlevout)*intsize !lb twb 2009
      nbytesmo(imodel)=nbytpmod
      ibytfrst(imodel)=ipos0
cTEST
      if(iprtlv.gt.1) then
        write(6,"('byte position before modes:',i5)") ipos0
        write(6,"('bytes per mode:',i5)") nbytpmod
      endif
      iostat=2
      nmod=0
      ipos=ipos0
      nnmax=0
      llmax=0
      nmodrad=0
      nmodsph=0
      nmodtor=0
      do while(iostat.eq.2)
        call bffi(lu,1,nn,4,iostat,nr,ipos+1)
        if(iostat.ne.3)then
          nmod=nmod+1
          if(nmod.gt.mxmodes) then
            stop 'too many modes'
          endif
          call bffi(lu,1,ityp,4,iostat,nr,ipos+intsize+1)!lb twb 2009
          call bffi(lu,1,ll,4,iostat,nr,ipos+intsize*2+1)!lb twb 2009
          modetype(nmod)=ityp
          lord(nmod)=ll
          nord(nmod)=nn
          ipos=ipos+nbytpmod
          if(ll.eq.0.and.ityp.eq.1)then
             nmodrad=nmodrad+1
          elseif(ll.gt.0.and.ityp.eq.3)then
             nmodsph=nmodsph+1
          elseif(ll.gt.0.and.ityp.eq.2)then
             nmodtor=nmodtor+1
          else
             print*,"mode should be spheroidal toroidal or radial"
             stop
          endif
          nnmax=max0(nn,nnmax)
          llmax=max0(ll,llmax)
cTEST
c        write(*,"(a5,1x,6(i6,2x))")"mode:",ityp,nn,ll,nmodrad,nmodsph,nmodtor
        endif
      enddo
c------------end of loop while iostat.eq.2
      if(iprtlv.gt.1) then
        write(6,"('total number of modes:',i6)") nmod
        write(6,"('total number of radial modes:',i6)") nmodrad
        write(6,"('total number of spheroidal modes:',i6)") nmodsph
        write(6,"('total number of toroidal modes:',i6)") nmodtor
        write(6,"('maximum n:',i4)") nnmax
        write(6,"('maximum l:',i4)") llmax
      endif
c
c---- loop on modes and fill buffers in common
c
      nmoddone=0
      imodbuffer=0
c
c---- write out the radial modes
c
      do imod=1,nmod
        if(modetype(imod).eq.1.and.lord(imod).eq.0) then
          nmoddone=nmoddone+1
          if(.not.fundonly) then
            imodbuffer=imodbuffer+1
            call putmodebuffer(lu,imod,imodbuffer,imodel,
     #            hn,hn2,hn3,dr,ilin,ipos0,nbytpmod)
          endif
        endif
      enddo
c
c---- write out the spheroidal modes
c
      nn=0
      do while(nn.le.nnmax)
        do imod=1,nmod
          if(modetype(imod).eq.3.and.lord(imod).gt.0.and.nord(imod).eq.nn) then
            nmoddone=nmoddone+1
            if(.not.overonly) then
              imodbuffer=imodbuffer+1
              call putmodebuffer(lu,imod,imodbuffer,imodel,
     #            hn,hn2,hn3,dr,ilin,ipos0,nbytpmod)
            endif
          endif
        enddo
        nn=nn+1
      enddo
c
c---- write out the toroidal modes
c
      nn=0
      do while(nn.le.nnmax)
        do imod=1,nmod
          if(modetype(imod).eq.2.and.lord(imod).gt.0.and.nord(imod).eq.nn) then
            nmoddone=nmoddone+1
            if(.not.overonly) then
              imodbuffer=imodbuffer+1
              call putmodebuffer(lu,imod,imodbuffer,imodel,
     #            hn,hn2,hn3,dr,ilin,ipos0,nbytpmod)
            endif
          endif
        enddo
        nn=nn+1
      enddo
c
      if(nmod.ne.nmoddone) then
        print*,"nmod different from nmoddone"
        write(6,"(3i6)") nmod,nmoddone
        stop
      endif
      nmodemod(imodel)=imodbuffer
c
      return
      end
c-----------------------------------------------------------------------
      subroutine putmodebuffer(lu,imod,imodb,imodel,
     #         hn,hn2,hn3,dr,ilin,ipos0,nbytpmod)
      include 'gemodes.h'
        common /sizeint/ intsize
      dimension buff(4)
      dimension buff2(4)
c
c      write(6,"(3i5,4g15.5,3i5)") imod,imodb,imodel,
c     #         hn,hn2,hn3,dr,ilin,ipos0,nbytpmod
      ipos=ipos0+(imod-1)*nbytpmod
      call bffi(lu,1,nn,4,iostat,nr,ipos+1)
      if(iostat.ne.3) then
        call bffi(lu,1,ityp,4,iostat,nr,ipos+1+intsize)
        call bffi(lu,1,ll,4,iostat,nr,ipos+1+intsize*2)
        call bffi(lu,1,f1,4,iostat,nr,ipos+1+intsize*3)
        isize=4
        call bffi(lu,1,f2,4,iostat,nr,ipos+1+intsize*3+isize)
        call bffi(lu,1,f3,4,iostat,nr,ipos+1+intsize*3+isize*2)
        call bffi(lu,1,f4,4,iostat,nr,ipos+1+intsize*3+isize*3)
        call bffi(lu,1,f5,4,iostat,nr,ipos+1+intsize*3+isize*4)
        call bffi(lu,1,f6,4,iostat,nr,ipos+1+intsize*3+isize*5)
        call bffi(lu,1,f7,4,iostat,nr,ipos+1+intsize*3+isize*6)
        call bffi(lu,1,f8,4,iostat,nr,ipos+1+intsize*3+isize*7)
cTEST
c        write(*,"(2(i4,1x),8(e12.5,2x)))")ityp,ll,f1,f2,f3,f4,f5,f6,f7,f8
c        pause
c
c---- eigenfunction info
c
        ishif=33+intsize*3
        call bffi(lu,1,buff, 16,iostat,nr,ipos+ishif+(ilin-1)*24)
        call bffi(lu,1,buff2,16,iostat,nr,ipos+ishif+ ilin*24)

        call splinege(buff(1),buff(2),buff2(1),buff2(2),dr,hn,hn2,hn3,u,up,upp)
        call splinege(buff(3),buff(4),buff2(3),buff2(4),dr,hn,hn2,hn3,v,vp,vpp)
cc        write(6,"(3i4,6g15.5)") nn,ityp,ll,f1,f2,f3,f4,f5,f6
c        write(6,"('above',4g15.5)") buff
c        write(6,"('at it',6g15.5)") u,up,v,vp,upp,vpp
c        write(6,"('below',4g15.5)") buff2
c        pause
      else
        stop 'error in putmodebuffer'
      endif
      imodefil(imodb,imodel)=imod
      iovermod(imodb,imodel)=nn
      itypemod(imodb,imodel)=ityp
      ilordmod(imodb,imodel)=ll
      omegamod(imodb,imodel)=f1
      smallqmod(imodb,imodel)=f2
      gvelmod(imodb,imodel)=f3
      pvelmod(imodb,imodel)=6371.*f1/(float(ll)+0.5)
      vaccmod(imodb,imodel)=f4
      haccmod(imodb,imodel)=f5
      vdismod(imodb,imodel)=f6
      hdismod(imodb,imodel)=f7
      potmod(imodb,imodel)=f8
      umod(imodb,imodel)=u
      upmod(imodb,imodel)=up
      uppmod(imodb,imodel)=upp
      vmod(imodb,imodel)=v
      vpmod(imodb,imodel)=vp
      vppmod(imodb,imodel)=vpp
      return
      end
c
      subroutine splinege(f1,fp1,f2,fp2,dr,hn,hn2,hn3,v,vp,vpp)
c
c     f1 is f(x1)
c     fp1 is df/dx(x1)
c     f2 is f(x2)
c     fp2 is df/dx(x2)
c     dr is (xx-x1) --- where xx is the desired x value
c     hn is x2-x1
c     hn2 is 1/(x2-x1)**2
c     hn3 is 1/(x2-x1)**3
c
c     v is f(xx)
c     vp is df/dx(xx)
c     vp2 is d2f/dx2(xx)
c
      a=hn3*(hn*(fp1+fp2)+2.*(f1-f2))
      b=hn2*(3.*(f2-f1)-hn*(fp2+2.*fp1))
      v=f1+dr*(fp1+dr*(b+dr*a))
      vp=fp1+dr*(2.*b+3.*dr*a)
      vpp=2.*b+6.*dr*a
      return
      end
c-----------------------------------------------------------------------
      subroutine findmode(imodel,nn,ii,ll,ifound)
      include 'gemodes.h'
      ifound=-1
      do imod=1,nmodemod(imodel)
        if(iovermod(imod,imodel).eq.nn) then
          if(itypemod(imod,imodel).eq.ii) then
            if(ilordmod(imod,imodel).eq.ll) then
	      ifound=imod
            endif
          endif
        endif
      enddo
      return
      end
