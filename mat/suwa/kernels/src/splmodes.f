      include 'gemodes.h'
      character*80 structure1
      character*80 structure2
      parameter (twopi=6.2831853072)
      character*80 infile
c
      parameter (maxspl=40)
      parameter (maxefs=12)
      dimension splef1(maxspl,maxefs,2)
      dimension splef2(maxspl,maxefs,2)
c
      parameter (maxpts=100)
      dimension omarr(maxpts)
      dimension farr(maxpts)
      dimension var1(maxpts)
      dimension var2(maxpts)
c
      character*8 formx,formyl,formyr
      data formx,formyl,formyr /'(i3)  ','(f8.5)  ','(i3)  '/
      data ix1,ix2,iy1,iy2 /700,3600,400,2600/
      logical more
c
      write(6,"('type name of first structure file')")
      read(5,"(a)") structure1
      write(6,"('type name of second structure file')")
      read(5,"(a)") structure2
c
      call openfl(1,structure1,1,0,0,ierr,-1)
      call openfl(2,structure2,1,0,0,ierr,-1)
c
    1 continue
      write(6,"('type first depth')") 
      read(5,*) depth1
      write(6,"('type second depth')") 
      read(5,*) depth2
      call rewfl(1)
      call rewfl(2)
      call getgemodes(1,1,.true.,.false.,.false.,depth1)
      call getgemodes(2,2,.true.,.false.,.false.,depth2)
c
c      omsp0=twopi/333.
      omsp0=twopi/500.
      domsp=twopi/750.
c      ommin=twopi/333.
      ommin=twopi/500.
c      ommax=twopi/40.
      ommax=twopi/40.
c
      write(6,"('type wave type 1-Love, 2-Rayleigh')")
      read(5,*) lora
      write(6,"('type 1 to plot ratios')")
      read(5,*) ifratio
c
c---- call the function which calculates b-splines of eigenfunctions
c
      call spleigen(1,lora,ommin,ommax,omsp0,domsp,
     #              nspl,splef1,maxspl,maxefs)
      call spleigen(2,lora,ommin,ommax,omsp0,domsp,
     #              nspl,splef2,maxspl,maxefs)
c
      more=.true.
      do while(more)
c
c---- decide what to plot
c
      
        write(6,"('plot vertical acceleration at the surface ---------------- 1')")
        write(6,"('plot horizontal acceleration at the surface -------------- 2')")
        write(6,"('plot U at the source ------------------------------------- 3')")
        write(6,"('plot dU/dr at the source --------------------------------- 4')")
        write(6,"('plot V (W) at the source --------------------------------- 5')")
        write(6,"('plot dV/dr (dW/dr) at the source ------------------------- 6')")
        read(5,"(i3)") iselect
c
        npts=50
        dom=(ommax-ommin)/float(npts-1)
        do ipt=1,npts
          omarr(ipt)=ommin+float(ipt-1)*dom
          farr(ipt)=1000.*omarr(ipt)/twopi
        enddo
c
        if(iselect.eq.1) ivar=9
        if(iselect.eq.2) ivar=10
        if(iselect.eq.3) ivar=5
        if(iselect.eq.4) ivar=6
        if(iselect.eq.5) ivar=7
        if(iselect.eq.6) ivar=8
        do ipt=1,npts
          omega=omarr(ipt)
            var1(ipt)=ecbspl(omega,omsp0,domsp,nspl,splef1(1,ivar,lora))
            var2(ipt)=ecbspl(omega,omsp0,domsp,nspl,splef2(1,ivar,lora))
            if(ifratio.eq.1) then
              var1(ipt)=var1(ipt)/var2(ipt)
              var2(ipt)=var1(ipt)
            endif
        enddo
        ampmax=-1.e20
        ampmin=1.e20
        do ipt=1,npts
          ampmax=amax1(ampmax,var1(ipt))
          ampmin=amin1(ampmin,var1(ipt))
          ampmax=amax1(ampmax,var2(ipt))
          ampmin=amin1(ampmin,var2(ipt))
        enddo
c
        write(6,"('max and min',2e12.3)") ampmin,ampmax
        if(ifratio.eq.1) then
          if(ampmin.gt.1.) ampmin=0.95
          if(ampmax.lt.1.) ampmax=1.05
        else
          if(ampmin.lt.0..and.ampmax.lt.0.) ampmax=0.
          if(ampmin.gt.0..and.ampmax.gt.0.) ampmin=0.
        endif
c  
        xmin=farr(1)
        xmax=farr(npts)
        ymin=ampmin
        ymax=ampmax
c
        tfx=5.
        sx=5.
        sy=(ymax-ymin)/5.
        tfy=ymin
        call chterm(5,6)
        call getcol('splmodes_colors',13)
        call lincol(1)
        call filcol(0)
        call thicktext(1.2)
        call twindo(ix1,ix2,iy1,iy2)
        call dwindo(xmin,xmax,ymin,ymax)
        call linlinboxl(xmin,xmax,ymin,ymax,tfx,sx,
     #       tfy,sy,formx,formyl,70)
c
c---- compare the eigenfunctions with what is obtained for each mode
c
        do ivar=1,2
	  if(ivar.eq.1) call lincol(2)
	  if(ivar.eq.2) call lincol(3)
          do ipt=1,npts
            xx=farr(ipt)
            if(ivar.eq.1) then
              yy=var1(ipt)
            else if(ivar.eq.2) then
              yy=var2(ipt)
            endif
            if(ipt.eq.1) then
              call movea(xx,yy)
            else
              call drawa(xx,yy)
            endif
          enddo
          call tsend
        enddo
      enddo
      end
c
c
      subroutine linlinboxl(xmin,xmax,ymin,ymax,tfx,sx,
     1tfy,sy,formx,formy,iheight)
      character*(*) formx,formy
      character*80 string
      call movea(xmin,ymin)
      call drawa(xmax,ymin)
      call drawa(xmax,ymax)
      call drawa(xmin,ymax)
      call drawa(xmin,ymin)
      iterm=4
      x=tfx
      length=llen(formx)
      intgr=0
      do 11 i=1,length
        if(formx(i:i).eq.'i'.or.formx(i:i).eq.'I') then
          intgr=1
          go to 1
        endif
   11 continue
    1 continue
      if(ibetween(xmin,xmax,x).eq.1) then
        call movea(x,ymin)
        call drwrel(0,15*iterm)
        if(formx(1:1).ne.' ') then
          if(intgr.eq.0) write(string,formx)x
          if(intgr.eq.1) write(string,formx)nint(x)
          call movea(x,ymin)
          call censtr(ileft,string,iheight)
          idummy=nint(float(iheight)*1.5)
          if(iheight.ne.0)call textr(-ileft,-idummy,iheight,string,0.)
        endif
        x=x+sx
        goto 1
      endif
      y=tfy
    2 continue
c      write(6,"('ymin,ymax,y',3f10.3)")ymin,ymax,y
      if(ibetween(ymin,ymax,y).eq.1) then
        call movea(xmax,y)
        call drwrel(-15*iterm,0)
        y=y+sy
        go to 2
      endif
      x=tfx
    3 continue
      if(ibetween(xmin,xmax,x).eq.1) then
        call movea(x,ymax)
        call drwrel(0,-15*iterm)
        x=x+sx
        go to 3
      endif
      y=tfy
      length=llen(formy)
      intgr=0
      do 14 i=1,length
        if(formy(i:i).eq.'i'.or.formy(i:i).eq.'I') then
          intgr=1
          go to 4
        endif
   14 continue
    4 continue
      if(ibetween(ymin,ymax,y).eq.1) then
        call movea(xmin,y)
        call drwrel(15*iterm,0)
        call movea(xmin,y)
        call tsend
        if(formy(1:1).ne.' ') then
          if(intgr.eq.0) write(string,formy)y
          if(intgr.eq.1) write(string,formy)nint(y)
          call rjuststr(ileft,string,iheight)
          idummy=nint(float(iheight)*0.25)
          if(iheight.ne.0)
     #      call textr(-ileft-iheight/4,-idummy,iheight,string,0.)
        endif
        y=y+sy
        go to 4
      endif
      return
      end
      subroutine linlinboxr(xmin,xmax,ymin,ymax,tfx,sx,
     1tfy,sy,formx,formy,iheight)
      character*(*) formx,formy
      character*80 string
      call movea(xmin,ymin)
      call drawa(xmax,ymin)
      call drawa(xmax,ymax)
      call drawa(xmin,ymax)
      call drawa(xmin,ymin)
      iterm=4
      x=tfx
      length=llen(formx)
      intgr=0
      do 11 i=1,length
        if(formx(i:i).eq.'i'.or.formx(i:i).eq.'I') then
          intgr=1
          go to 1
        endif
   11 continue
    1 continue
      if(ibetween(xmin,xmax,x).eq.1) then
        call movea(x,ymin)
c        call drwrel(0,15*iterm)
        if(formx(1:1).ne.' ') then
          if(intgr.eq.0) write(string,formx)x
          if(intgr.eq.1) write(string,formx)nint(x)
          call movea(x,ymin)
          call censtr(ileft,string,iheight)
          idummy=nint(float(iheight)*1.5)
c          if(iheight.ne.0)call textr(-ileft,-idummy,iheight,string,0.)
        endif
        x=x+sx
        goto 1
      endif
      y=tfy
    2 continue
c      write(6,"('ymin,ymax,y',3f10.3)")ymin,ymax,y
      if(ibetween(ymin,ymax,y).eq.1) then
        call movea(xmax,y)
        call drwrel(-15*iterm,0)
        y=y+sy
        go to 2
      endif
      x=tfx
    3 continue
      if(ibetween(xmin,xmax,x).eq.1) then
        call movea(x,ymax)
        call drwrel(0,-15*iterm)
        x=x+sx
        go to 3
      endif
      y=tfy
      length=llen(formy)
      intgr=0
      do 14 i=1,length
        if(formy(i:i).eq.'i'.or.formy(i:i).eq.'I') then
          intgr=1
          go to 4
        endif
   14 continue
    4 continue
      if(ibetween(ymin,ymax,y).eq.1) then
        call movea(xmax,y)
        call drwrel(-15*iterm,0)
        call movea(xmax,y)
        call tsend
        if(formy(1:1).ne.' ') then
          if(intgr.eq.0) write(string,formy)y
          if(intgr.eq.1) write(string,formy)nint(y)
c          call rjuststr(ileft,string,iheight)
          idummy=nint(float(iheight)*0.25)
          if(iheight.ne.0)
     #      call textr(iheight/4,-idummy,iheight,string,0.)
c     #      call textr(-ileft-iheight/4,-idummy,iheight,string,0.)
        endif
        y=y+sy
        go to 4
      endif
      return
      end
