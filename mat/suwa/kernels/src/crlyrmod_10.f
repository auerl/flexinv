* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                             *
*     This set of algorithm creates 1D reference model        *
*     profiles files by combining crust1 with anisotropic     *
*     prem. The output files are simple ascii files ready     *
*     to be processed with our mode based kernel code         *
*                                                             *
*     (c) Ludwig Auer & Lapo Boschi 2014                      *
*                                                             *
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      program crlyrmod

        implicit double precision (a-h,o-z)
        
        common/earmod/ rnorm,coef(4,8,20),xb(20),xt(20)
        common/earmod/ numlyr,nic,noc,moho,nl(20)

        parameter(np=9,nlo=360,nla=180)      

        character*120 fname
        character*120 cdirname
        character*120 refmfile
        character*120 odirname
        character*80 getunx
        character*3 miter ! iteration index or index of starting model

        dimension vp(np,nla,nlo)
        dimension vs(np,nla,nlo)
        dimension bnd(np,nla,nlo)
        dimension rhoin(np,nla,nlo)

        real*4 vpout(8)
        real*4 vsout(8)
        real*4 rhout(8)
        real*4 thout(8)
        real*4 topo

c
c---- read the reference model once
c

        call chekcl('|   -R:r:1:Reference Earth model file path'//
     #              '|   -C:r:1:Absolute CRUST1 model dir path'// 
     #              '|   -O:r:1:Kernel output directory path'// 
     #              '|   -v:o:1:Verbosity level[0]|')

        cdirname=getunx('-C',1,nbyts)
        refmfile=getunx('-R',1,nbyts)
        odirname=getunx('-O',1,nbyts)

        open(13,file=trim(refmfile))
        
        call rdearmod(13,icm,ici,nlyr,rnorm)
        write(6,"('Number of layers in original PREM:',i3)") numlyr
        do ilyr=1,numlyr
           write(6,"('sublayers between radii:',2i4,2f9.0)") 
     &     ilyr,nl(ilyr),1000.*xb(ilyr),1000.*xt(ilyr)
        enddo
        ! moho depth should be top of 10th layer
        write(6,"('moho depth?',f10.2)") rnorm-xt(10)

c
c---- read the crustal model CRUST1 once
c

        print*,' .... reading all maps ... ' 
        open(51,file=trim(cdirname)//'/crust1.vp')
        open(52,file=trim(cdirname)//'/crust1.vs')
        open(53,file=trim(cdirname)//'/crust1.rho')
        open(54,file=trim(cdirname)//'/crust1.bnds')
        do j=1,nla
           do i=1,nlo
              read(51,*)(vp(k,j,i),k=1,np)
              read(52,*)(vs(k,j,i),k=1,np)
              read(53,*)(rhoin(k,j,i),k=1,np)
              read(54,*)(bnd(k,j,i),k=1,np)
           enddo
        enddo
        close(51)
        close(52)
        close(53)
        close(54)

        ! Iteration index
        itern=0
        write(miter,'(i3.3)')itern

        ! Looop over all 64800 locations
        nker=0
        do ilat=1,nla ! 1 to 180         
           do ilon=1,nlo ! 1 to 360

              ilon2=ilon+180
              if (ilon2>360) then ! shift grid by 180
                 ilon2=ilon2-360
              end if
              topo=bnd(1,ilat,ilon2)
              do i=1,np-1
                 thout(i)=bnd(i,ilat,ilon2)-bnd(i+1,ilat,ilon2) ! thickness
                 vpout(i)=vp(i,ilat,ilon2) ! vp
                 vsout(i)=vs(i,ilat,ilon2) ! vs
                 rhout(i)=rhoin(i,ilat,ilon2) ! rho
              end do
              nker=nker+1
              print*,"# of kernel",nker
              write(fname,"('/L_',i3.3,i3.3,'_')")ilat,ilon
              call extrmod(vpout,vsout,thout,rhout,topo,fname,odirname)
           end do
        end do

      end program crlyrmod
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine extrmod(vpshallow,vsshallow,thshallow,
     &                   rhshallow,elevation,fname,odirname)

        implicit double precision (a-h,o-z)

        common/earmod/ rnorm,coef(4,8,20),xb(20),xt(20)
        common/earmod/ numlyr,nic,noc,moho,nl(20)

        parameter(ianisotropic=1) ! our model is anisotropic

        real*4 vpshallow(9)
        real*4 vsshallow(9)
        real*4 rhshallow(9)
        real*4 thshallow(8)
        real*4 elevation

        dimension xtnew(20)
        dimension xbnew(20)
        dimension nlnew(20)

        character*80 fname
        character*120 odirname
c
c---- flip ice and water
c            
        tmpvs=vsshallow(2)
        tmpvp=vpshallow(2)
        tmprh=rhshallow(2)
        tmpth=thshallow(2)

        vpshallow(2)=vpshallow(1)
        vsshallow(2)=vsshallow(1)
        rhshallow(2)=rhshallow(1)
        thshallow(2)=thshallow(1)

        vpshallow(1)=tmpvp
        vsshallow(1)=tmpvs
        rhshallow(1)=tmprh
        thshallow(1)=tmpth

c
c---- remove ice on top of water
c

        if(thshallow(1).gt.0..and.thshallow(2).gt.0.) thshallow(1)=0.

c
c---- only consider layers thicker than 500 meters
c

        thick=0.d0
        do i=1,8 ! run over all layers in crust1, used to be 7 in crust2
           if(thshallow(i).gt.0.5) thick=thick+dble(thshallow(i))
        enddo
        write(6,"('thickness of new shallow structure:',f10.2)") thick

c
c--- nlaytot is essentially the 10 layers below moho in prem plus
c--- the 8 layers contained in crust1
c
        nlytot = 18

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
                 nlnew(ilyr)=nl(ilyr)*2 ! this just increases the nr of sublayers
              else 
                 nlnew(ilyr)=nl(ilyr)*4 ! this just increases the nr of sublayers
              endif
           else if(ilyr.eq.10) then

              ! top of 10th layer is the moho depth in prem
              ! in case new moho depth (depends on "thick")
              ! is different, the new number of sublayers is
              ! 10. This is basically always the case

              xtnew(ilyr)=rnorm-thick+elevation

              if(xtnew(ilyr).ne.xbnew(ilyr)) nlnew(ilyr)=10

           else 

              ! only consider layers which are thicker than 500m
              if(thshallow(nlytot+1-ilyr).gt.0.5) then
                 xtnew(ilyr)=xtnew(ilyr-1)+thshallow(nlytot+1-ilyr)
              else
                 xtnew(ilyr)=xtnew(ilyr-1)
              endif

              if(xtnew(ilyr).ne.xbnew(ilyr)) nlnew(ilyr)=8 ! number of sublayers=8

           endif
           xbnew(ilyr+1)=xtnew(ilyr)
        enddo
c
c---- now write out the model
c
        ntot=0
        do ilyr=1,nlytot
           write(6,"('sublayers between radii:',2i4,2f9.0)") 
     &     ilyr,nlnew(ilyr),1000.*xbnew(ilyr),1000.*xtnew(ilyr)
           ntot=ntot+nlnew(ilyr)
        enddo

c
c---- determine all moho layers
c        
        nmoho=0
        do i=1,10
          nmoho=nmoho+nlnew(i)
        enddo

        open(1,file=trim(odirname)//trim(fname))
        write(1,"(a)") fname(1:lnblnk(fname))
        write(1,"('1 1. 1 1')")
        write(1,"(7i5)") ntot,nlnew(1),nlnew(1)+nlnew(2),nmoho

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
                  call evanisomod(ras2,rho,vpv,vsv,qka,qmu,vph,vsh,eta)
                  call evanisomod(ras3,rho2,vpv2,vsv2,qka2,qmu2,vph2,vsh2,eta2)
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
              else
                rho=rhshallow(nlytot+1-ilyr)
                vpv=vpshallow(nlytot+1-ilyr)
                vsv=vsshallow(nlytot+1-ilyr)
                vph=vpv
                vsh=vsv
                eta=1.d0
                qmu=600.d0
                if(nlytot+1-ilyr.eq.2) qmu=0.
                qka=57822.5
              endif
              write(1,"(f8.0,3f9.2,2f9.1,2f9.2,f9.5)")
     #          1000.d0*ras,1000.d0*rho,1000.d0*vpv,1000.d0*vsv,
     #          qka,qmu,1000.d0*vph,1000.d0*vsh,eta
            enddo
          endif
        enddo
        close(1)

      end subroutine extrmod
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






