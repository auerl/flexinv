      subroutine initmodel(imodel)
      common /plevel/ iprtlv
c
c---- some constants: big G, Earth rotation, Pi, Average density
c
      parameter (bigg=6.6723e-11)
      parameter (capom=7.292115e-05)
      parameter (pi=3.1415926)
      parameter (rhon=5514.3)
c
      include 'gemodes.h'
      include 'gemodl.h'
c
c---- some work arrays
c
      dimension wk(3,nknts)
      dimension coef(5)
      dimension y(3),y1(3),y2(3),y3(3),ydot(3)
c
c---- copy parameters from gemodes common block
c
      n=nlevmodel(imodel)
      nic=ninnercore(imodel)
      noc=noutercore(imodel)
      moho=nmohorovic(imodel)
      do i=1,n
        r(i)=radlayer(i,imodel)*1000.
        rho(i)=   rholayer(i,imodel)
        acon(i)=  vphlayer(i,imodel)
        ccon(i)=  vpvlayer(i,imodel)
        lcon(i)=  vsvlayer(i,imodel)
        ncon(i)=  vshlayer(i,imodel)
        fcon(i)=  etalayer(i,imodel)
        if(qmulayer(i,imodel).ne.0.) then
          qshear(i)=1./qmulayer(i,imodel)
        endif
        if(qkalayer(i,imodel).ne.0.) then
          qkappa(i)=1./qkalayer(i,imodel)
        endif
      enddo
      rn=r(n)
      ifanis=ifanismodel(imodel)
c
      if(ifanis.eq.0) then
        do i=1,n
          acon(i)=ccon(i)
          ncon(i)=lcon(i)
          fcon(i)=1.
        enddo
      endif
c
c---- find top of topmost solid layer
c
      nsl=n
      do while(lcon(nsl).eq.0.)
        nsl=nsl-1
      enddo
c
c---- normalize radii and density
c
      do i=1,n
        r(i)=r(i)/rn
        rho(i)=rho(i)/rhon
      enddo
c
c---- spline the density
c
      call rspln(1,n,r,rho,qrho,wk)
c
c---- integrate to get ellipticity and gravity
c
c
      y(1) = 1.
      y(2) = 0.
      y(3) = rho(1)
      rr   = 1.e-5
c
      do 100 i=2,n
      rstep=r(i)-r(i-1)
      if(abs(rstep).gt.1.e-5) goto 110
      ell(i)=ell(i-1)
      eta(i)=eta(i-1)
      g(i)=g(i-1)
      goto 100
c
  110 do 101 ii=1,3
  101 y1(ii)=y(ii)
      rr1=rr
      if(rr+rstep.gt.r(i)) rstep=r(i)-rr
      iback=1
      goto 1100
c
 1201 do 201 ii=1,3
  201 y2(ii)=y(ii)
      rr2=rr
  901 rstep=rstep*.5
      do 321 ii=1,3
  321 y(ii)=y1(ii)
      rr=rr1
      iback=2
      goto 1100
c
 1202 do 401 ii=1,3
  401 y3(ii)=y(ii)
      rr3=rr
      iback=3
      goto 1100
c
 1203 do 501 ii=1,3
      if(abs(y(ii)-y2(ii)).gt..5e-5) goto 601
  501 continue
      goto 701
c
  601 do 801 ii=1,3
  801 y2(ii)=y3(ii)
      rr2=rr3
      goto 901
  701 if(abs(rr-r(i)).lt.1.e-5) goto 1300
      rstep=4.*rstep
      goto 110
c
 1300 ell(i)=y(1)
      eta(i)=rr*y(2)/y(1)
      g(i)=4.*rr*y(3)/3.
      goto 100
c
 1100 k=krunge(3,y,ydot,rr,rstep)
      if(k.ne.1) goto (1201,1202,1203),iback
c
      ydot(1)=y(2)
      t=rr-r(i-1)
      im1=i-1
      rhot=rho(im1)+t*(qrho(1,im1)+t*(qrho(2,im1)+t*qrho(3,im1)))
      ydot(3)=3.*(rhot-y(3))/rr
      ydot(2)=-2.*(ydot(3)*y(1)+3.*rhot*y(2))/(y(3)*rr)
      goto 1100
c
  100 continue
c
c
 1400 factr=.75*g(n)/r(n)
      rhobar=rhon*factr
      wn=(pi*bigg*rhobar)**.5
      vn=rn*wn
      gn=rn*wn**2
      fac=2.5*(capom/wn)**2/(4.*ell(n)*(eta(n)+2.)/3.)
      g(1)=0.
      ell(1)=1.
      do 1401 i=1,n
      rho(i)=rho(i)/factr
      g(i)=g(i)/factr
      ell(i)=ell(i)*fac
      do 1401 j=1,3
 1401 qrho(j,i)=qrho(j,i)/factr
c
c
c
      do 7 i=1,n
      ccon(i)=rho(i)*(ccon(i)/vn)**2
      lcon(i)=rho(i)*(lcon(i)/vn)**2
      acon(i)=rho(i)*(acon(i)/vn)**2
      ncon(i)=rho(i)*(ncon(i)/vn)**2
    7 fcon(i)=fcon(i)*(acon(i)-2.*lcon(i))
c
c
      call rspln(1,n,r,ccon,qccon,wk)
      call rspln(1,n,r,lcon,qlcon,wk)
      call rspln(1,n,r,acon,qacon,wk)
      call rspln(1,n,r,ncon,qncon,wk)
      call rspln(1,n,r,fcon,qfcon,wk)
c
      ndisc=0.
      do 56 i=2,n
      if(abs(r(i)-r(i-1)).gt.1.e-5) goto 56
      ndisc=1+ndisc
      ndsc(ndisc)=i-1
   56 continue
      ndisc=1+ndisc
      ndsc(ndisc)=n
c
      write(6,909) rhobar,1./ell(n),nic,noc,moho,nsl,n
  909 format(///' mean density =',f10.3,'     surface ellipticity ='
     #   ,' 1/',f7.3/' nic =',i4,'   noc =',i4,'   moho =',i4
     #   ,'   nsl =',i4,'   n =',i4///)
      i=0
      rlas=-10000.
      lines=9
  903 continue
      if(iprtlv.gt.3) then
        write(6,902)
      endif
  902 format(1x,'level',
     1 3x,'radius',7x,'rho',8x,'vpv',8x,'vph',7x,'vsv',
     2 7x,'vsh',6x,'eta',7x,'g',7x,'1/qm',5x,'1/qk'
     3  ,'      1/ell    etaell'/)
  906 i=1+i
      if(i.gt.n) goto 99
      rr=r(i)*rn
      rrho=rho(i)*rhobar
      vpv=sqrt(ccon(i)/rho(i))*vn
      vph=sqrt(acon(i)/rho(i))*vn
      vsv=sqrt(lcon(i)/rho(i))*vn
      vsh=sqrt(ncon(i)/rho(i))*vn
      eeta=fcon(i)/(acon(i)-2.*lcon(i))
      gg=g(i)*gn
      if(abs(rlas-rr).lt.500.) then
        if(iprtlv.gt.3) then 
          write(6,912)
        endif
      endif
  912 format(1x,132('-'))
      if(rr.eq.rlas) lines=1+lines
      rlas=rr
      if(iprtlv.gt.3) then
        write(6,9021) i,rr,rrho,vpv,vph,vsv,vsh,eeta,gg
     #   ,qshear(i),qkappa(i),ell(i),eta(i)
      endif
 9021 format(2x,i3,f11.1,f10.3,2f11.3,2f10.3,4f9.5,f11.7,f9.5)
      lines=1+lines
      if(lines.lt.57) goto 906
      lines=1
      goto 903
   99 continue
      return
      end
c
c
c
      function krunge(n,y,f,x,h)
c ??? some sort of integration or interpolation? x is 
c incremented on second and fourth calls. Resets itself after 
c 5'th call.
c
c input: (these are guesses)
c  n   = number of points
c  f() = function evaluated at each point
c  x   = independent variable
c  h   = step size
c output:
c  y() = 
c   
c correction: put savey in save statement, b.w. 8/1/90
c
      dimension phi(6),savey(6),y(6),f(6)
      data m/0/
      save savey
c
      m = m + 1
      goto(1,2,3,4,5), m
c***
1     krunge=1
      return
c***
2     do 22 j=1,n
      savey(j) = y(j)
      phi(j)   = f(j)
22    y(j) = savey(j)+0.5d0*h*f(j)
      x = x + 0.5d0*h
      krunge = 1
      return
c***
3     do 33 j=1,n
      phi(j) = phi(j) + 2.d0*f(j)
33    y(j)   = savey(j)+0.5d0*h*f(j)
      krunge = 1
      return
c***
4     do 44 j=1,n
      phi(j) = phi(j)+2.0d0*f(j)
44    y(j)   = savey(j)+h*f(j)
      x      = x + 0.5d0*h
      krunge = 1
      return
c***
5     do 55 j=1,n
55    y(j) = savey(j)+(phi(j)+f(j))*h/6.d0
      m    = 0
      krunge = 0
      return
      end
