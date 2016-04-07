      character*80 name,name1,name2,string,getunx
      integer*4 iret
      logical exists
      integer unlink
      real*8 dffirst,dflast,dprecis,dfgrav
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common /geout/ depmin,depmax,ioutflag(1000),nlevout
c
      call chekcl('|     :r:1:Earth model|')
      iin=7
      write(6,"('enter spherical-earth model file ')")
      name=getunx('',1,nbyts)
c     read(5,"(a)") name
      open(iin,file=name,status='old',form='formatted',iostat=iret)
      if(iret.ne.0) call error_$print(iret)
c
      iout=8
      write(6,"('type base name for output file')")
ccc      read(5,"(a)") string
      string=name
      name1=string(1:lnblnk(string))//'.out'
      name2=string(1:lnblnk(string))//'.bin'
      open(iout,file=name1,form='formatted',iostat=iret)
      if(iret.ne.0) call error_$print(iret)
c
      write(6,"('type minimum and maximum depth for output file')")
cccc      read(5,*) depmin,depmax
      depmin=0.
      depmax=3.
c
      ioeig=9
c
      inquire(file=name2,exist=exists)
      if(exists) idummy=unlink(name2)
      open(ioeig,file=name2)
      close(ioeig)
      call openfl(ioeig,name2,4,0,0,istat,-1)
c
      call bffo(ioeig,1,name,80,istat,0)
      call bffo(ioeig,1,depmin,4,istat,0)
      call bffo(ioeig,1,depmax,4,istat,0)
c
      call model(iin,iout,ioeig)
      close(iin)
      write(6,"('maximum frequency in millihertz')")
cccc      read(5,*) dflast
      dflast=70.d0
      dffirst=0.1d0
      dfgrav=100.
      lfirst=2
      llast=2000
      dprecis=1.d-7
      ngren=1
c
      do jcom=2,3
        call wtable(iout,ioeig,jcom,
     #      dffirst,dflast,lfirst,llast,dprecis,dfgrav,ngren)
      enddo
      close(iout)
      close(ioeig)
      stop
      end
c
c---------------------------------------------
c
      subroutine error_$print(iret)
      print *,'read error',iret
      stop
      end
c
c---------------------------------------------
c
      subroutine wtable(iout,ioeig,jjcom,
     #       dffirst,dflast,lfirst,llast,dprecis,dfgrav,ngren)
c*** makes up table of frequencies ***
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      common/mtab/we(430),de(430),ke(430),wtry(250),bm(250),um(250)
      dimension wt(215)
      data nmx/215/,inss/5/
      data ifirst/0/
      jcom=jjcom
      cmhz=pi/500.d0
      stepf=1.d0
c      print *,'enter eps and wgrav'
c      read(*,*) eps,wgrav
      eps=dprecis
      wgrav=dfgrav
c---
      eps1=eps
      eps2=eps
      wgrav=wgrav*cmhz
      if(ifirst.eq.0) then
        write(iout,100) eps,eps1,wgrav
        ifirst=1
      endif
  100 format(/,'integration precision =',g12.4,'  root precision =',
     +   g12.4,'  gravity cut off =',g12.4,' rad/s',///,6x,'mode',
     +   8x,'w(rad/s)',7x,'w(mhz)',10x,'t(secs)',6x,'grp vel(km/s)',
     +   8x,'q',13x,'raylquo',/)
      call steps(eps)
c      print *,'enter lmin,lmax,wmin,wmax,nbran'
c      read(*,*) lmin,lmax,wmin,wmax,nbran
      lmin=lfirst
      lmax=llast
      wmin=dffirst
      wmax=dflast
      nbran=ngren
c---
      wmin=wmin*cmhz
      wmax=wmax*cmhz
c
      call geout1(ioeig,eps,wgrav,lmin,lmax,wmin,wmax,nbran)
c
      if(nbran.le.0.or.nbran.gt.nmx) nbran=nmx
      do 1 i=1,250
      bm(i)=0.d0
    1 wtry(i)=0.d0
      if(lmin.le.0) lmin=1
      nev=2
      wt(1)=wmin
      wt(2)=wmax
      if(jcom.ne.1) goto 6
      lmin=0
      lmax=0
    6 l=lmin-1
   10 l=l+1
      if(l.gt.2.and.l.lt.10) l=l+1
      if(l.ge.11.and.l.lt.50) l=l+4
      if(l.ge.51.and.l.lt.300) l=l+9
      if(l.ge.301.and.l.lt.450) l=l+14
      if(l.gt.451) l=l+19
      if(l.gt.lmax.or.wt(1).ge.wmax) return
      knsw=1
      maxo=inss
      fl=l
      fl1=fl+1.d0
      fl2=fl+fl1
      fl3=fl*fl1
      sfl3=dsqrt(fl3)
c*** determine mode count ***
      we(1)=wt(1)
      call detqn(we(1),ke(1),de(1),0)
      imax=2*nmx
      do 15 i=2,imax
      we(i)=wmax
   15 ke(i)=-10
      do 20 i=2,nev
      call entry(wt(i),imax,kei)
      nmode=kei-ke(1)
      kkk=i
      if(nmode.ge.nbran) goto 25
   20 continue
   25 continue
c      print 900,wt(1),wt(kkk),nmode,nbran
  900 format(' count between ',g14.6,' and',g14.6,' =',i8,' : max =',i4)
      if(nmode.le.0) return
      if(nmode.gt.nbran) nmode=nbran
      imax=2*nmode
c*** fill up table using bisection ***
      indx=2
      ichk=ke(1)+1
   35 if(ke(indx).ne.ichk) goto 40
      indx=indx+2
      if(indx.gt.imax) goto 60
      ichk=ichk+1
      goto 35
   40 i1=indx-1
   45 indx=indx+2
      if(indx.ge.imax) goto 46
      ichk=ichk+1
      if(ke(indx).ne.ichk) goto 45
   46 i2=min0(indx,imax)
      wtst=0.5d0*(we(i2)+we(i1))
      if((we(i2)-we(i1))/wtst.lt.eps2) goto 50
      call entry(wtst,imax,ktst)
      indx=i1+1
      ichk=ke(i1)+1
      goto 35
   50 print 901,we(i1),ke(i1),we(i2),ke(i2)
  901 format('problem in table : ',2(g16.8,i5))
      j1=i1+1
      j2=i2-1
      do 55 i=j1,j2
   55 de(i)=1.d0
      indx=i2+2
      if(indx.ge.imax) goto 60
      ichk=ke(i2)+1
      goto 35
   60 nev=imax
c      print 902,(i,we(i),ke(i),de(i),i=1,nev)
  902 format(i5,1pd17.7,i6,1pd17.7)
c*** find roots ***
      knsw=0
      maxo=8
      call rotspl(nev,eps1,wmin,wmax,wt,iout,ioeig)
      goto 10
      end
c
c---------------------------------------------
c
      subroutine entry(w,imax,kei)
      implicit real*8(a-h,o-z)
      common/mtab/we(430),de(430),ke(430),wtry(250),bm(250),um(250)
      call detqn(w,kei,dei,0)
      indx=min0(max0(2*(kei-ke(1)),1),imax)
      if(indx.eq.1.and.we(1).lt.w) goto 10
      if(indx.eq.imax.and.we(imax).gt.w) goto 10
      if(kei.ne.ke(indx)) goto 5
      if(we(indx).gt.w) goto 10
      indx=indx+1
      if(we(indx).lt.w) goto 10
      return
    5 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      indx=indx+1
   10 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      return
      end
c
c---------------------------------------------
c
      subroutine rotspl(nev,eps1,wmin,wmax,wt,iout,ioeig)
c*** find roots by spline interpolation ***
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/mtab/we(430),de(430),ke(430),wtry(250),bm(250),um(250)
      dimension x(20),det(20),qx(3,20),wrk(60),wt(1),ichar(4)
      data tol/1.d-9/,itmax/15/,ichar/' s',' t',' s',' c'/
      knev=1
      nmo=nev/2
      do 100 i=1,nmo
      k1=2*i
      k=k1-1
      if(de(k)*de(k1).gt.0.d0) goto 100
      nord=ke(k1)
      if(l.eq.1) nord=nord+1
      nind=nord+1
      det(1)=de(k)
      det(2)=de(k1)
      x(1)=we(k)
      x(2)=we(k1)
      if(wtry(nind).le.x(1)) goto 5
      if(wtry(nind).ge.x(2)) goto 5
      call detqn(wtry(nind),knt,ftry,0)
      if(ftry*det(1).lt.0.d0) goto 5
      x(1)=wtry(nind)
      det(1)=ftry
    5 c=x(1)
      if(dabs(det(1)).gt.dabs(det(2))) c=x(2)
c      print 910,x(1),det(1),x(2),det(2)
  910 format(4g18.10)
   10 j=1
      m=2
      ntry=2
      b=0.5d0*(x(1)+x(2))
   15 t=dabs(b*eps1)
      if(dabs(b-c).lt.t) goto 65
      call detqn(b,knt,fb,0)
c      print 900,b,fb
  900 format(2g20.12)
      ind=1
      do 20 m=2,ntry
      ind=ind+1
   20 if(b.lt.x(m)) goto 25
   25 ntry=ntry+1
      j2=ntry
   30 j1=j2-1
      x(j2)=x(j1)
      det(j2)=det(j1)
      if(j1.eq.ind) goto 35
      j2=j2-1
      goto 30
   35 x(ind)=b
      det(ind)=fb
      idn=0
      do 40 m=2,ntry
      idn=idn+1
      iup=idn+1
   40 if(det(idn)*det(iup).le.0.d0) goto 45
   45 ind=iup
      if(dabs(det(idn)).lt.dabs(det(iup))) ind=idn
      c=x(ind)
      if(ntry.ge.itmax) goto 60
      call dsplin(ntry,x,det,qx,wrk)
      del=-det(ind)/qx(1,ind)
   50 delx=-det(ind)/(qx(1,ind)+del*qx(2,ind))
      if(dabs(delx-del).lt.tol) goto 55
      if(del*delx.lt.0.d0) goto 60
      del=delx
      goto 50
   55 b=c+delx
      if(b.ge.x(idn).and.b.le.x(iup)) goto 15
   60 x(1)=x(idn)
      x(2)=x(iup)
      det(1)=det(idn)
      det(2)=det(iup)
      goto 10
c*** write out frequencies ***
   65 call detqn(b,knt,fb,1)
c      print 900,b,fb
      wdiff=(b-wray*wn)/b
      tcom=2.d0*pi/b
      wmhz=1000.d0/tcom
      gcom=vn*cg/1000.d0
      qmod=0.d0
      if(qinv.gt.0.d0) qmod=1.d0/qinv
      print 200,nord,ichar(jcom),l,b,wmhz,tcom,gcom,qmod,wdiff
      write(iout,200) nord,ichar(jcom),l,b,wmhz,tcom,gcom,qmod,wdiff
  200 format(i5,a2,i5,6g16.7)
      call modout(b,qmod,gcom,ioeig)
      if(nord.eq.0) wmin=b
      if(bm(nind).le.0.d0) goto 70
      blin=b+cg*wn
      bp2=5.d0*bm(nind)-4.d0*b+2.d0*(um(nind)+2.d0*cg*wn)
      bdiff=dmax1(dabs(blin-bp2),eps1*bp2)
      bup=bp2+bdiff
      bdwn=bp2-bdiff
      wtry(nind)=bdwn
      goto 71
   70 bup=b+2.d0*cg*wn
      wtry(nind)=b
   71 bm(nind)=b
      um(nind)=cg*wn
      if(bup.ge.wmax) goto 100
      knev=knev+1
      wt(knev)=bup
  100 continue
      nev=knev+1
      wt(1)=wmin
      wt(nev)=wmax
      return
      end
c
c---------------------------------------------------------------
c
      subroutine modout(wcom,qmod,gcom,ioeig)
c        output modified to Harvard(JW) catalog convention
      implicit real*8(a-h,o-z)
      real*4 abuf,buf,ww,qq,gc,au,av,usur,vsur,ps
      dimension abuf(6009)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/eifx/a(14,1000),idum(1000)
      common/buf$/nn,jj,ll,ww,qq,gc,au,av,usur,vsur,ps,buf(6000)
      data bigg,tau,rhobar/6.6723e-11,1000.0,5515.0/
      data fot/1.33333333333333d0/
      equivalence (nn,abuf)
c
      common /geout/ depmin,depmax,ioutflag(1000),nlevout
      real*4 depmin,depmax
c
      dsfl3=1.0
      if (jcom.ne.1) dsfl3=1.0/sfl3

      if (jcom.eq.2.or.jcom.eq.4) then
        do i= 1,n
          a(3,i)=a(1,i)
          a(4,i)=a(2,i)
          a(1,i)=0.
          a(2,i)=0.
        enddo
      elseif (jcom.eq.1) then
        do i= 1,n
          a(3,i)=0.
          a(4,i)=0.
        enddo
      endif
      
c             calculating acceleration at surface
      esur=a(3,nsl)*dsfl3
      au= -((wsq+2.*fot)*a(1,nsl) + fl1*a(5,nsl) )*(wn*wn)
      av= -(wsq*esur-a(1,nsl)*fot-a(5,nsl))*(wn*wn)
c             gravitational potential at surface
      usur=sngl(a(1,nsl))
      vsur=sngl(a(3,nsl)*dsfl3)
      ps=a(5,nsl)
      
      nn=nord
      ll=l
      ww=wcom
      qq=1.0/qmod
      gc=gcom
      jj=jcom
      if(jcom.ne.2) goto 5
      do 1 i=1,noc
      a(3,i)=0.d0
    1 a(4,i)=0.d0
    5 continue
      icount=0
      j=1
      do 10 i=1,n
        if(ioutflag(i).eq.1) then
          buf(j)=a(1,i)
          buf(j+1)=a(2,i)
          buf(j+2)=a(3,i)*dsfl3
          buf(j+3)=a(4,i)*dsfl3
          buf(j+4)=a(5,i)
          buf(j+5)=a(6,i)
          j=j+6
          icount=icount+1
        endif
   10 continue
      if(icount.ne.nlevout) write(6,"('error in modout',2i10)") n,icount
      if(j.ne.6*nlevout+1) write(6,"('error in modout',2i10)") j,icount
      call bffo(ioeig,1,nn,(11+6*nlevout)*4,istat,0)
      return
      end
c
c-----------------------------------------------
c
      subroutine geout1(ioeig,deps,dwgrav,lmin,lmax,dwmin,dwmax,nbran)
      real*8 deps,dwgrav,dwmin,dwmax
      real*4  eps, wgrav, wmin, wmax
      integer ioeig
      integer lmin
      integer lmax
      integer nbran
      integer ifirst
      data ifirst/0/
      save ifirst
c
      if(ifirst.eq.0) then
        eps=sngl(deps)
        wgrav=sngl(dwgrav)
        wmin=sngl(dwmin)
        wmax=sngl(dwmax)
        call bffo(ioeig,1,lmin,4,istat,0)
        call bffo(ioeig,1,lmax,4,istat,0)
        call bffo(ioeig,1,nbran,4,istat,0)
        call bffo(ioeig,1,wmin,4,istat,0)
        call bffo(ioeig,1,wmax,4,istat,0)
        call bffo(ioeig,1,eps,4,istat,0)
        call bffo(ioeig,1,wgrav,4,istat,0)
        ifirst=1
      endif
      return
      end
