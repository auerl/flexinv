c
c ***
c *** 
c ***
c
      character*1 chtesttype(3)

      parameter (maxker=100)
      parameter (maxpar=8)
      common/mdpert/omegaref,gvelref,numker,numpar,parm(maxker,maxpar)
c
      character*80 getunx,structurename,layerfilename,xivoigt
      parameter (twopi=6.2831853072)
      character*80 infile,layerfile
      character*80 outfile
      common /plevel/ iprtlv

      dimension layers(30)

      chtesttype(2)="T"
      chtesttype(3)="S"

      call chekcl('|   -I:r:1:Earth model '//
     #            '|   -L:r:1:Layer file '// 
     #            '|   -P:r:1:Xi-vS yes [y] or no [n]'// 
     #            '|   -v:o:1:Verbosity level [0]|')

      structurename=getunx('-I',1,nbyts)
      infile=structurename(1:nbyts)//'.bin'

      layerfilename=getunx('-L',1,nbyts)
      layerfile=layerfilename(1:nbyts)

      xivoigt=getunx('-P',1,nbyts)
      print*,"Xi-Voigt average parameterization? ",xivoigt
      
      open(333,file=layerfile,status='old')
      read(333,*),numker
      print*,'Number of layers ',numker
      layers(1)=0.
      do il=2,numker+1
         read(333,*),layer
         layers(il)=layer
      enddo
      close(333)

      call openfl(1,infile,1,0,0,ierr,-1)
      outfile=structurename(1:lnblnk(structurename))//'_kernels'
      open(12,file=outfile)
c
c---- set imodel to 1 and depth to 2 km (arbitrary, but required)
c---- getgemodes fills the common blocks modelpara and modepara
c
      lu1=1
      imodel=1
      depth=2.
      call getgemodes(lu1,imodel,.false.,.false.,.false.,depth)
c
c---- call the program which splines the earth model and
c---- integrates the ellipticity of figure (similar to modl.f)
c
      call initmodel(imodel)
c
c---- loop on modes
c
c      nn=0 !overtone number

      do nn=0,6
      do ii=2,3 ! spheroidal, toroidal (or vice-versa?)
        do ll=0,1000 ! angular degree

c---- maps overtone number, mode type, angular degree to a single index imod
          call findmode(imodel,nn,ii,ll,imod)

c---- now integrate eigenfunctions of mode imod to find sensitivity kernels
          if(imod.gt.0) then !negative imod means we dont use that mode but interpolate others

c---- call program that calculates eigenfunctions and puts them in common/mode/
            call initmode(lu1,imodel,imod)

c---- calculate kernel and store result
            call integrate(layers,xivoigt)

            if(ii.eq.2) then 
              write(12,"(1x,i2,1x,a1,i4,2e15.7,i3,i3)") nn,'T',ll,omegaref,gvelref,numker,numpar
            else
              write(12,"(1x,i2,1x,a1,i4,2e15.7,i3,i3)") nn,'S',ll,omegaref,gvelref,numker,numpar
            endif
            do i=1,numker
              write(12,"(i3,8e12.4)") i,(parm(i,j),j=1,numpar)
            enddo
c
          endif
c
        enddo
      enddo
      enddo!loop over overtone number
      call closfl(lu1,ierror)
      close(12)
      end
