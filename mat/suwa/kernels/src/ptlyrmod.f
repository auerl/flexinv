c--uses the routine radial_basis
      implicit double precision (a-h,o-z)
      character*80 infile
      character*80 outfile
      character*80 model
      character*80 pertmodel
      character*80 pertfile
      character*80 layerfile
      character*80 layerfilename
      character*80 line2
      parameter (mxlyr=1000)
      dimension ras(mxlyr)
      dimension rho(mxlyr)
      dimension vpv(mxlyr)
      dimension vph(mxlyr)
      dimension vsv(mxlyr)
      dimension vsh(mxlyr)
      dimension qka(mxlyr)
      dimension qmu(mxlyr)
      dimension eta(mxlyr)
      character*1 type
      real*4 radius4
      real*4 value4
c
      parameter (maxker=30)
      parameter (maxpar=8)
      dimension parm(maxker,maxpar)
      dimension layers(30)
c
      character*80 getunx
c
      call chekcl('|   -I:r:1:Input Earth model'//
     #            '|   -P:r:1:Mantle perturbation file'//
     #            '|   -L:r:1:Layer file'// 
     #            '|   -O:r:1:Output Earth model'//
     #            '|   -v:o:1:Verbosity level[0]|')
c
      infile=getunx('-I',1,nbyts)
      pertfile=getunx('-P',1,nbyts)
      outfile=getunx('-O',1,nbyts)
      layerfilename=getunx('-L',1,nbyts)
      layerfile=layerfilename(1:lnblnk(layerfilename))

      open(333,file=layerfile,status='old')
      read(333,*),numker
      print*,'Number of layers ',numker
      layers(1)=0.
      do il=2,numker+1
         read(333,*),layer
         layers(il)=layer
      enddo
      close(333)
      
      open(1,file=infile)
c
      read(1,"(a)") model
      lmod=lnblnk(model)
      imod1=1
      do while(model(imod1:imod1).eq.' '.and.imod1.le.lmod)
        imod1=imod1+1
      enddo
      read(1,"(a)") line2
      read(1,*) ntot,nic,noc,moho
      if(moho.eq.0) then
        stop 'moho level has to be defined'
      endif
      do il=1,ntot
              read(1,"(f8.0,3f9.2,2f9.1,2f9.2,f9.5)")
     #          ras(il),rho(il),vpv(il),vsv(il),
     #          qka(il),qmu(il),vph(il),vsh(il),eta(il)
      enddo
      close(1)
c
c---- check nic, noc, moho
c
      write(6,"('ic boundary:',f12.3)") ras(nic)
      write(6,"('            ',f12.3)") ras(nic+1)
      write(6,"('oc boundary:',f12.3)") ras(noc)
      write(6,"('            ',f12.3)") ras(noc+1)
      write(6,"('moho       :',f12.3)") ras(moho)
      write(6,"('            ',f12.3)") ras(moho+1)
c
c---- choose perturbation
c
      open(1,file=pertfile)
      read(1,*) type,nparm,nkern
      write(6,*) type,nparm,nkern
      do i=1,nkern
        read(1,"(8e12.4)") (parm(i,j),j=1,nparm)
      enddo
      close(1)
c
      do il=noc+1,moho
        radius4=sngl(ras(il)*0.001d0)
        do ikern=1,nkern
 	  istatus=0	! Lapo 20Feb04
          call radial_basis(radius4,type,ikern,value4,istatus,layers)
          if(istatus.ne.0) then
            write(6,"('radial_basis returns',i10)") istatus
            stop
          endif
          do iparm=1,nparm
            if(iparm.eq.0) then
            else if(iparm.eq.1) then
              rho(il)=rho(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.2) then
              vph(il)=vph(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.3) then
              vpv(il)=vpv(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.4) then
              vsh(il)=vsh(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.5) then
              vsv(il)=vsv(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.6) then
              eta(il)=eta(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.7) then
              qka(il)=qka(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else if(iparm.eq.8) then
              qmu(il)=qmu(il)*(1.d0+parm(ikern,iparm)*dble(value4)*0.01d0)
            else
              stop
            endif
          enddo
        enddo
      enddo
c
c---- write out the new perturbed model
c
      write(pertmodel,"(a,'_',a1)") 
     #      model(imod1:lmod),type
      open(2,file=outfile)

      write(2,"(a)") pertmodel(1:lnblnk(pertmodel))
      write(2,"(a)") line2(1:lnblnk(line2))
      write(2,"(4i5)") ntot,nic,noc,moho
      do il=1,ntot
              write(2,"(f8.0,3f9.2,2f9.1,2f9.2,f9.5)")
     #          ras(il),rho(il),vpv(il),vsv(il),
     #          qka(il),qmu(il),vph(il),vsh(il),eta(il)
        
      enddo
      close(2)
      end
