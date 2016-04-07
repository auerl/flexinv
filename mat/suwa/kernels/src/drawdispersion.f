      character*80 xtext,ytextl,ytextr
      character*80 filein
      logical exists,condition
      dimension polygon(2,100)
      parameter (maxdim=1000)
      dimension c(maxdim,2,100)
      dimension np(100)
      dimension cmax(2),cmin(2)
      dimension mnemx(5),mnemy(5)
      common/parx/ xmin,xmax,xfir,xinc
      common/paryl/ yminl,ymaxl,yfirl,yincl
      common/paryr/ yminr,ymaxr,yfirr,yincr
      data mnemx/4,'xmin','xmax','xfir','xinc'/
      data mnemy/4,'ymin','ymax','yfir','yinc'/
      character*8 formx,formyl,formyr
      logical more
      data formx,formyl,formyr /'(i3)  ','(i3)  ','(i3)  '/
      data ix1,ix2,iy1,iy2 /700,3600,400,2600/
      character*1 cdummy
      character*80 string
c
      character*80 earth
c
      parameter (maxspl=25)
      dimension omarr(maxspl)
      dimension delcgc(maxspl)
      dimension delc(maxspl)
      dimension delcgcspl(3,maxspl)
      dimension delcspl(3,maxspl)
      parameter (twopi=6.2831853)
c
      nfiles=0
      more=.true.
      do while(more)
        write(6,"('type name of next input file')")
        read(5,"(a)") filein
        if(filein(1:1).eq.' ') then
          more=.false.
        else
          inquire(file=filein,exist=exists)
          if(.not.exists) then
            write(6,"('file does not exist')")
          else
            nfiles=nfiles+1
            open(1,file=filein)
            condition=.true.
            i=0
            do while (condition)
              read(1,"(a)",iostat=ios) string
c              write(6,"(a)") string
              if(ios.eq.0) then
                iskip=0
                do j=1,lnblnk(string)
                  if(string(j:j).eq.':') iskip=j
                enddo
                if(string(1:1).ne.'#') then
                  read(string(iskip+1:lnblnk(string)),*,iostat=ios)
     #            c(i+1,1,nfiles),c(i+1,2,nfiles)
                endif
              endif
              if(ios.ne.0) then
                condition=.false.
              else
                if(string(1:1).ne.'#') i=i+1
              endif
            enddo
            close(1)
            np(nfiles)=i
            write(6,"(i5,' points')") np(nfiles)
          endif
        endif
      enddo
c
      do ic=1,2
        cmax(ic)=-1.e30
        cmin(ic)=1.e30
        do ifile=1,nfiles
          do j=1,np(ifile)
            cmax(ic)=amax1(c(j,ic,ifile),cmax(ic))
            cmin(ic)=amin1(c(j,ic,ifile),cmin(ic))
          enddo
        enddo
        write(6,"('column ',i1,':',g15.5,' -- ',g15.5)") ic,cmin(ic),cmax(ic)
      enddo
c
      write(6,"('type 1 for column 1 to be x axis')")
      iord=1
c      read(5,"(i1)") iord
      if(iord.eq.1) then
        ixc=1
        iyc=2
      else
        ixc=2
        iyc=1
      endif
      write(6,"('parameters for x axis:')")
      call setupg('drawdispersion_x1',5,6,13,mnemx,xmin)
      write(6,"('type label for x-axis')")
      xtext='frequency (mHz)'
c      read(5,"(a)") xtext
      write(6,"('type format for x-labels')")
      formx='(i3)'
c     read(5,"(a)") formx
   
      write(6,"('parameters for left y axis:')")
      call setupg('drawdispersion_y1',5,6,13,mnemy,yminl)
      write(6,"('type label for left y-axis')")
      ytextl='\\delta\\c\/c (%)'
c      read(5,"(a)") ytext
      write(6,"('type format for y-labels')")
      formyl='(f5.1)'
c      read(5,"(a)") formy
c
c      write(6,"('parameters for right y axis:')")
c      call setupg('drawgtr1_y2',5,6,13,mnemy,yminr)
c      write(6,"('type label for right y-axis')")
c      ytextr='CMT residual misfit'
c      read(5,"(a)") ytext
c      write(6,"('type format for y-labels')")
c      formyr='(f5.3)'
c      read(5,"(a)") formy
c
      call chterm(5,6)
      call getcol('drawdispersion_colors',13)
      call lincol(1)
      call filcol(0)
      call thicktext(1.2)
      call twindo(ix1,ix2,iy1,iy2)
      call dwindo(xmin,xmax,yminl,ymaxl)
c
c      call boxfill(6,r1,r2,r3,r4,150,150)
      polygon(1,1)=xmin
      polygon(2,1)=yminl
      polygon(1,2)=xmin
      polygon(2,2)=ymaxl
      polygon(1,3)=xmax
      polygon(2,3)=ymaxl
      polygon(1,4)=xmax
      polygon(2,4)=yminl
      call filcol(0)
      call polyfil(0,polygon(1,1),polygon(2,1),2,4)
c
      call linwdt(6)
      call texth(1,1,10,'\\duplex\\',0.)
      call linlinbox(xmin,xmax,yminl,ymaxl,xfir,xinc,yfirl,yincl,formx,formyl,70)
      call censtr(ileft,xtext,90)
      ixmid=(ix2+ix1)/2
      call texth(ixmid-ileft,iy1-220,90,xtext,0.)
      call censtr(ileft,ytextl,90)
      iymid=(iy2+iy1)/2
      call texth(ix1-300,iymid-ileft,90,ytextl,90.)
c
c      call twindo(ix1,ix2,iy1,iy2)
c      call dwindo(xmin,xmax,yminr,ymaxr)
c      call linwdt(6)
c      call linlinboxr(xmin,xmax,yminl,ymaxl,xfir,xinc,yfirl,yincl,formx,formyl,70)
c      call linlinboxr(xmin,xmax,yminr,ymaxr,xfir,xinc,yfirr,yincr,formx,formyr,70)
c      call censtr(ileft,ytextr,90)
c      iymid=(iy2+iy1)/2
c      call texth(ix2+350,iymid+ileft,90,ytextr,-90.)
      
      call linwdt(3)
      call lincol(3)
      call filcol(0)
      do ifile=1,nfiles
        if(ifile.eq.1) then
          call twindo(ix1,ix2,iy1,iy2)
          call dwindo(xmin,xmax,yminl,ymaxl)
c        else if(ifile.eq.2) then
c          call twindo(ix1,ix2,iy1,iy2)
c          call dwindo(xmin,xmax,yminr,ymaxr)
        endif
        write(6,"('type symbol size,type,fill color for file',i2)") ifile
        isiz=20
        ityp=ifile
c        if(ifile.eq.1) then
c          ityp=3
c          icol=31
c        else if(ifile.eq.2) then
c          ityp=2
c          icol=33
c        else if(ifile.eq.3) then
c          ityp=1
c          icol=35
c        endif
c        read(5,*) isiz,ityp,icol
c        call filcol(icol)
      do ip=1,np(ifile)
        if(ip.eq.1) then
          call movea(c(ip,ixc,ifile),c(ip,iyc,ifile))
        else
          call drawa(c(ip,ixc,ifile),c(ip,iyc,ifile))
        endif
      enddo
        do ip=1,np(ifile)
c          call symbfil(c(ip,ixc,ifile),c(ip,iyc,ifile),isiz,ityp)
        enddo
        call tsend
      enddo
c
      earth='/home/columbus/ekstrom/DISPERSION'
      lu=13
      call getemodl(lu,earth,ierror)
    1 continue
      write(6,"('type 1 to draw local dispersion from observed maps')")
      read(5,"(i1)") ifso
      if(ifso.eq.1) then
      write(6,"('type 1 for love waves')")
      read(5,"(i1)") lora
      if(lora.ne.1) lora=2
      write(6,"('type geocentric lat lon')")
      read(5,*) elat,elon
      call getdelcpntspl(elat,elon,lora,omarr,delc,delcspl,maxspl,nspl,ierror)
      if(ierror.eq.0) then
        do i=1,30
          freq=float(i)*0.001
          omega=twopi*freq
          pert=rsple(1,nspl,omarr,delc,delcspl,omega)
          call symbfil(float(i),pert,20,1)
        enddo
      endif
      endif
      go to 1
      call tsend
      call addtoimage
      call adjcol('drawdispersion_colors',5,6,13)
      call saveterm
      call quitterm
      end
 
c=====================================================
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
