c
c ***
c ***
c

      implicit double precision (a-h,o-z)
      character*80 name
      common/earmod/ rnorm,coef(4,8,20),xb(20),xt(20),numlyr,nic,
     # noc,moho,nl(20)
      dimension xtnew(20)
      dimension xbnew(20)
      dimension nlnew(20)
c
      real*4 vpshallow(8),vsshallow(8),rhshallow(8),thshallow(7)
      real*4 tpmo
      character*2 type
      character*60 desc
      character*5 crmod
      character*80 polyfile
      character*80 getunx
      character*80 crustfoldname
      character*128 inflb
c
      character*80 string
      character*80 dbsdir
      logical abspath
      logical exists
c
      call chekcl('|   -R:r:1:Reference Earth model file'//
     #            '|   -a:o:0:Use absolute path for Ref. Earth'//
     #            '|   -C:r:1:Code for type of crust AA_AA'//
     #            '|   -M:r:1:Crust 2.0 model file'// 
     #            '|   -v:o:1:Verbosity level[0]|')

      crustfoldname=getunx('-M',1,nbyts)
      string=getunx('-a',0,nbyts)

      ldbs=lnblnk(crustfoldname)
      inflb=crustfoldname(1:ldbs)//"/CNtype2_key.txt"

      if(nbyts.eq.0) then 
        abspath=.true.
      endif

      if(abspath) then
        polyfile=getunx('-R',1,nbyts)

c---------polyfile will be PREM-type reference model (e.g. GREM_0)
      else
        call getenv('HRVDBS',dbsdir)
        ldbs=lnblnk(dbsdir)
        string=getunx('-R',1,nbyts)
        lstr=lnblnk(string)
        polyfile=dbsdir(1:ldbs)//'/REMS_1D/'//string(1:lstr)
      endif
      lpoly=lnblnk(polyfile)
c
      inquire(file=polyfile,exist=exists)
      if(.not.exists) then
        write(6,"(a)") polyfile(1:lpoly)
        stop 'no such 1-D reference model'
      endif

      open(13,file=polyfile)
c
      call rdearmod(13,icm,ici,nlyr,rnorm)
      write(6,"('number of layers defined:',i3)") numlyr
      do ilyr=1,numlyr
        write(6,"('sublayers between radii:',2i4,2f9.0)") 
     #       ilyr,nl(ilyr),1000.*xb(ilyr),1000.*xt(ilyr)
      enddo
c
c---- moho depth should be top of 10th layer
c
      write(6,"('moho depth?',f10.2)") rnorm-xt(10)
c
c---- find parameters for shallow (crustal) structure
c
      crmod=getunx('-C',1,nbyts)
      if(nbyts.ne.5) stop 'wrong crustal type'
c      
      type=crmod(1:2)
      ibath=0
      ielev=0
      if(crmod(4:4).eq.'B') then
          read(crmod(5:5),"(i1)") ibath
      else if(crmod(4:4).eq.'E') then
          read(crmod(5:5),"(i1)") ielev
      endif
      write(6,"('type:',a2)") type

      call cstruc_20(13,type,desc,vpshallow,vsshallow,rhshallow,thshallow,tpmo,inflb,ierr)

      if(ierr.ne.0) then
        write(6,"('error from cstruc_20')")
        stop
      endif
      write(6,"('crust:',a)") desc
      write(6,"(8f6.2)") vpshallow
      write(6,"(7f6.2)") thshallow

c
c---- remove ice on top of water
c
      if(thshallow(1).gt.0..and.thshallow(2).gt.0.) thshallow(1)=0.

c
c---- change water thickness
c
      thshallow(2)=float(ibath)
c
      thick=0.d0
      do i=1,7
        if(thshallow(i).gt.0.5) thick=thick+dble(thshallow(i))
      enddo
      write(6,"('thickness of new shallow structure:',f10.2)") thick
c
      nlytot=17
c
c---- fill in the earth model
c
      xbnew(1)=0.
      do ilyr=1,nlytot
        nlnew(ilyr)=0
        if(ilyr.lt.10) then
            xtnew(ilyr)=xt(ilyr)
            if(ilyr.lt.3) then
              nlnew(ilyr)=nl(ilyr)
            else if(ilyr.lt.8) then
              nlnew(ilyr)=nl(ilyr)*2
            else 
              nlnew(ilyr)=nl(ilyr)*4
            endif
        else if(ilyr.eq.10) then
            xtnew(ilyr)=rnorm-thick+dfloat(ielev)
            if(xtnew(ilyr).ne.xbnew(ilyr)) nlnew(ilyr)=10
        else 
            if(thshallow(18-ilyr).gt.0.5) then
              xtnew(ilyr)=xtnew(ilyr-1)+thshallow(18-ilyr)
            else
              xtnew(ilyr)=xtnew(ilyr-1)
            endif
            if(xtnew(ilyr).ne.xbnew(ilyr)) nlnew(ilyr)=8
        endif
        xbnew(ilyr+1)=xtnew(ilyr)
      enddo

c
c---- now write out the model
c
      ntot=0
      do ilyr=1,nlytot
          write(6,"('sublayers between radii:',2i4,2f9.0)") 
     #       ilyr,nlnew(ilyr),1000.*xbnew(ilyr),1000.*xtnew(ilyr)
          ntot=ntot+nlnew(ilyr)
      enddo
c
        nmoho=0
        do i=1,10
          nmoho=nmoho+nlnew(i)
        enddo
c
        name=crmod(1:5)

        write(6,"(a)") name(1:lnblnk(name))
        write(6,"(a)") crmod(1:lnblnk(crmod))

        open(1,file=name)

        write(6,"('Didnt get here')") 

        write(1,"(a)") name(1:lnblnk(name))
        write(1,"('1 1. 1 1')")
        write(1,"(7i5)") ntot,nlnew(1),nlnew(1)+nlnew(2),nmoho
c        
        do ilyr=1,nlytot
          nnn=nlnew(ilyr)
          if(nnn.gt.1) then
            rstep=(xtnew(ilyr)-xbnew(ilyr))/dfloat(nnn-1)
            do j=1,nnn
              ras=xbnew(ilyr)+dfloat(j-1)*rstep
	      ras2=ras
              if(j.eq.1) ras2=ras+1.d-5
              if(j.eq.nnn) ras2=ras-1.d-5
              if(ilyr.le.10) then

c
c---- the following is a hack to extrapolate gradients. Should only happen
c---- for layer 10 (the lithosphere)
c

                if(ras.gt.xt(ilyr)) then
                  if(ilyr.ne.10) then
                    write(6,"('something strange:',i10)") ilyr
                    stop
                  endif
                  ras3=2.*xt(ilyr)-ras2
                  ras2=xt(ilyr)-1.d-5
cc                  write(6,"('above below',3f12.2)") ras,ras2,ras3
                  call evanisomod(ras2,rho,vpv,vsv,qka,qmu,vph,vsh,eta)
                  call evanisomod(ras3,rho2,vpv2,vsv2,qka2,qmu2,vph2,vsh2,eta2)
cc                  write(6,"('at, below, above',3f12.2)") rho,rho2,2.*rho-rho2
                  rho=2.*rho-rho2
                  vpv=2.*vpv-vpv2
                  vsv=2.*vsv-vsv2
                  qka=2.*qka-qka2
                  qmu=2.*qmu-qmu2
                  vph=2.*vph-vph2
                  vsh=2.*vsh-vsh2
                  eta=2.*eta-eta2
                else
                  call evanisomod(ras2,rho,vpv,vsv,qka,qmu,vph,vsh,eta)
                endif
                if(qmu.gt.1.d-6) qmu=1./qmu
                qka=1./qka
c
                ianisotropic=1
                if(ianisotropic.ne.1) then
                  aaaa=rho*vph**2
                  cccc=rho*vpv**2
                  xnnn=rho*vsh**2
                  xlll=rho*vsv**2
                  ffff=eta*(aaaa-2.d0*xlll)
                  xkap=(1.d0/9.d0)*(4.d0*aaaa+cccc+4.d0*ffff-4.d0*xnnn)
                  xmu=(1.d0/15.d0)*(aaaa+cccc-2.d0*ffff+5.d0*xnnn+6.d0*xlll)
                  if(ilyr.eq.9.or.ilyr.eq.10) then
                    eta=1.d0
                  else
                    eta=1.d0
                  endif
                  vpv=dsqrt((xkap+(4.d0/3.d0)*xmu)/rho)
                  vph=vpv
                  vsv=dsqrt(xmu/rho)
                  vsh=vsv
                endif
c
              else
                rho=rhshallow(18-ilyr)
                vpv=vpshallow(18-ilyr)
                vsv=vsshallow(18-ilyr)
                vph=vpv
                vsh=vsv
c
                eta=1.d0
                qmu=600.d0
                if(18-ilyr.eq.2) qmu=0.
                qka=57822.5
              endif
              write(1,"(f8.0,3f9.2,2f9.1,2f9.2,f9.5)")
     #          1000.d0*ras,1000.d0*rho,1000.d0*vpv,1000.d0*vsv,
     #          qka,qmu,1000.d0*vph,1000.d0*vsh,eta
            enddo
          endif
        enddo
        close(1)
      end
