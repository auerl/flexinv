      subroutine calc_tessa(nf, arco, arlo, db)
c
c Computes tesslations on the earth with given frequency
c The equation for the number of nodes is  freq*freq*10+2
c edited from an earlier program called tessa and is used for
c spline inversions and stacking of SS precursors, J.G., 2000.
c Input:  
c	 nf  ---  frequency number (total = nf*nf*10+2)
c
c output:
c     arco  ---   latitude of the node
c     arlo  ---   longitude of the node
c
      
      implicit double precision (a-h,o-z)

c *** basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    


      integer nf
      dimension co1(50,50),xlo1(50,50),co67(50,50),xlo67(50,50)
      dimension co16(50,50),xlo16(50,50),colat(150,250),xlong(150,250)
      dimension arco(30000),arlo(30000),i_near(30,30000),k_near(30000)
      dimension x_near(20,30000)
      character*256 outfile

c
c check dimension
c
      if(maxfreq.lt.nf) stop 'too many horizontal nodes, quit...'
      pi=3.141592654
      rad=pi/180.
      c72=cos(rad*72.)
      base=acos(c72/(1.-c72))
      db=base*1.09/float(nf)
      call triangle1(nf,co1,xlo1)
      call triangle67(nf,co67,xlo67)
      call triangle16(nf,co16,xlo16)
      colat(1,1)=0.
      xlong(1,1)=0.
      i=1
      j=1
      nknt=1
 1011 format(2f10.3)
      arco(nknt)=colat(1,1)
      arlo(nknt)=xlong(1,1)
      do 1 i=2,nf+1
      j1=1
      j2=i
      knt=0
      do 1 n=1,5
      if(n.gt.1) j1=2
      do 1 j=j1,i
      if(n.eq.5.and.j.eq.i) go to 1
      knt=knt+1
      nknt=nknt+1
      colat(i,knt)=co1(i,j)
      xlong(i,knt)=xlo1(i,j)+float(n-1)*rad*72.
      if(xlong(i,knt).gt.pi) xlong(i,knt)=xlong(i,knt)-2.*pi
      arco(nknt)=colat(i,knt)
      arlo(nknt)=xlong(i,knt)
 1001 format(3i6,2f10.3)
    1 continue
      i1=nf+1
      do 2 i=2,nf+1
      knt=0
      do 2 n=1,5
      j1=1
      j2=nf+1
      if(n.gt.1) j1=2
      do 2 j=j1,j2
      if(n.eq.5.and.j.eq.nf+1) go to 2
      knt=knt+1
      nknt=nknt+1
      colat(nf+i,knt)=co67(i,j)
      xlong(nf+i,knt)=xlo67(i,j)+float(n-1)*rad*72.
      if(xlong(nf+i,knt).gt.pi) xlong(nf+i,knt)=xlong(nf+i,knt)-2.*pi
c      write(1,1001) i+nf,knt,nknt,colat(nf+i,knt)/rad,
c     1xlong(nf+i,knt)/rad
       arco(nknt)=colat(nf+i,knt)             
       arlo(nknt)=xlong(nf+i,knt)
    2 continue
      do 3 i=2,nf+1
      knt=0
      do 3 n=1,5
      j1=1
      if(n.gt.1) j1=2
      j2=nf+2-i
      do 3 j=j1,j2
      if(n.eq.5.and.j.eq.j2) go to 3
      knt=knt+1
      nknt=nknt+1
      colat(2*nf+i,knt)=co16(i,j) 
      xl=xlo16(i,j)+float(n-1)*rad*72.
      if(xl.gt.pi) xl=xl-2.*pi
      xlong(2*nf+i,knt)=xl
c      write(1,1001) i+2*nf,knt,nknt,colat(2*nf+i,knt)/rad,
c     1xlong(2*nf+i,knt)/rad  
      arco(nknt)=colat(2*nf+i,knt)
      arlo(nknt)=xlong(2*nf+i,knt)      
    3 continue
      avg=0.
      knt=0
      do 5 n=1,nknt
      k_near(n)=0
      c1=cos(arco(n))
      s1=sin(arco(n))
      cc1=cos(arlo(n))
      ss1=sin(arlo(n))
      do 6 i=1,nknt
      if(i.eq.n) go to 6
      c2=cos(arco(i))
      s2=sin(arco(i))
      cc2=cos(arlo(i))
      ss2=sin(arlo(i))
      cdelta=c1*c2+s1*s2*(cc1*cc2+ss1*ss2)
      if(abs(cdelta).gt.1)  then
      cdelta=sign(1.,cdelta)
      endif
      delta=acos(cdelta)
      if(delta.gt.1.5*db) go to 6
      avg=avg+delta/db
      knt=knt+1
      k_near(n)=k_near(n)+1
      k=k_near(n)
      i_near(k,n)=i
      x_near(k,n)=delta/db
    6 continue
      k=k_near(n)
      arco(n)=(0.5*pi-arco(n))/rad
      arlo(n)=arlo(n)/rad
 1003 format(i6,2f10.3,i6,15i6)  
    5 continue
      avg=avg/float(knt)
 1010 format(i10,f10.4)
      return
      end

      subroutine triangle1(nf,colat,xlong)
      implicit double precision (a-h,o-z)
      dimension colat(50,50),xlong(50,50)
      pi=3.141592654
      rad=pi/180.
      c72=cos(rad*72.)
      base=acos(c72/(1.-c72))
      db=base/float(nf)
      colat(1,1)=0.
      xlong(1,1)=0.
      do 1 i=2,nf+1
      xlong(i,1)=0.
      y=db*float(i-1)
      colat(i,1)=y
      x=(cos(y))**2+c72*(sin(y))**2
      xx=acos(x)
      xs=xx/float(i-1)
      sang=sin(y)*sin(72*rad)/sin(xx)
      ang=asin(sang)
      do 1 j=2,i
      xm=xs*float(j-1)
      cl=cos(y)*cos(xm)+sin(y)*sin(xm)*cos(ang)
      acl=acos(cl)
      colat(i,j)=acl
      slong=sang*sin(xm)/sin(acl)
      xlong(i,j)=asin(slong)
    1 continue
      do 2 i=1,nf+1
c      write(*,1000) i-1,(colat(i,j)/rad,xlong(i,j)/rad,j=1,i)
 1000 format(i4,14f7.2/(4x,14f7.2))
    2 continue
      return
      end
      subroutine triangle67(nf,colat,xlong)
      implicit double precision (a-h,o-z)
      dimension colat(50,50),xlong(50,50)
      pi=3.141592654
      rad=pi/180.
      c72=cos(rad*72.)
      s72=sin(rad*72.)
      c144=cos(rad*144.)
      s144=sin(rad*144.)
      base=acos(c72/(1.-c72))
c      write(*,*) ' type n-frequency'
c      read(*,*) nf
      db=base/float(nf)
      colat(1,1)=base
      xlong(1,1)=0.
      do 10 i=2,nf+1
      xlong(i,1)=0.
      x=db*float(i-1)
      ccl=cos(base)*cos(x)+sin(base)*sin(x)*c72
      alat=acos(ccl)
      colat(1,i)=alat
      xlong(1,i)=asin(s72*sin(x)/sin(alat))
c      write(*,*) alat/rad, xlong(1,i)/rad
   10 continue
      do 11 i=2,nf+1
      x=db*float(i-1)
      ccl=cos(base)*cos(x)+sin(base)*sin(x)*c144
      accl=acos(ccl)
      colat(i,1)=accl
      sgamma=sin(x)*s144/sin(accl)
      gamma=asin(sgamma)
c      write(*,*) accl/rad,sgamma,gamma/rad
      xlong(i,1)=gamma
      if(i.ge.nf+1) go to 11
      top=72.*rad-2.*gamma
      y=(cos(accl))**2+cos(top)*(sin(accl))**2
      ay=acos(y)
      cang=(ccl-ccl*y)/(sin(ay)*sin(accl))
      sang=sin(top)*sin(accl)/sin(ay)
      ang=acos(cang)
c      write(*,*) i, ang/rad, ay/rad, top/rad
      day=ay/float(nf+1-i)
      do 12 j=2,nf+2-i
      yy=day*float(j-1)
      cclat=cos(colat(i,1))*cos(yy)+sin(colat(i,1))*sin(yy)*cos(ang)
      aclat=acos(cclat)
      colat(i,j)=aclat
      sfi=asin(sin(yy)*sang/sin(aclat))
      xlong(i,j)=xlong(i,1)+sfi
   12 continue
   11 continue
      do 14 i=1,nf
      ii=nf+2-i
      do 14 j=1,nf+1-i
      jj=nf+2-j
      colat(ii,jj)=pi-colat(i,j)
      xlong(ii,jj)=108.*rad-xlong(i,j)
   14 continue
      do 13 i=1,nf+1
c      write(*,1000) i-1,(colat(i,j)/rad,xlong(i,j)/rad,j=1,nf+1)
 1000 format(i4,14f7.2/(4x,14f7.2))
   13 continue
      return
      end
      subroutine triangle16(nf,colat1,xlong1)
      implicit double precision (a-h,o-z)
      dimension colat(50,50),xlong(50,50),colat1(50,50),xlong1(50,50)
      pi=3.141592654
      rad=pi/180.
      xl36=36.*rad
      c72=cos(rad*72.)
      base=acos(c72/(1.-c72))
c      write(*,*) ' type n-frequency'
c      read(*,*) nf
      db=base/float(nf)
      colat(1,1)=0.
      xlong(1,1)=0.
      do 1 i=2,nf+1
      xlong(i,1)=0.
      y=db*float(i-1)
      colat(i,1)=y
      x=(cos(y))**2+c72*(sin(y))**2
      xx=acos(x)
      xs=xx/float(i-1)
      sang=sin(y)*sin(72*rad)/sin(xx)
      ang=asin(sang)
      do 1 j=2,i
      xm=xs*float(j-1)
      cl=cos(y)*cos(xm)+sin(y)*sin(xm)*cos(ang)
      acl=acos(cl)
      colat(i,j)=acl
      slong=sang*sin(xm)/sin(acl)
      xlong(i,j)=asin(slong)
    1 continue
      do 3 i=1,nf+1
      ii=nf+2-i
      do 3 j=1,i
      colat1(ii,j)=pi-colat(i,j)
      xlong1(ii,j)=xl36+xlong(i,j)
    3 continue
      do 2 i=1,nf+1
c      write(*,1000) i-1,(colat1(i,j)/rad,xlong1(i,j)/rad,j=1,nf+2-i)
 1000 format(i4,14f7.2/(4x,14f7.2))
    2 continue
      return
      end
