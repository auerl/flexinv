      subroutine integr(n,qints,iflag,jpath,i12)
c   subroutine to integrate contributions to delta,time, and
c   time derivatives (returned in qints) for layer n.
c     qints(1)=d
c     qints(2)=ddp
c     qints(3)=t
c     qints(4)=tstar
c     qints(5...)=tt model derivatives
c     iflag  =0(if ray traversed entire layer)
c            =1(if ray was reflected at the top of the layer)
c            =2(if ray turned in this layer)
c
      implicit double precision (a-h,o-z)
      dimension qints(*),ind(3)
      common/path$/delt(2,50000,2),nfin(2),kdep ! la 1600->50000
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/xnorm/xnorm
      common/dx/dxtdp,fxtpp,ifl
c---------------- for anisotropic inversion, Yu--------------------------
      common/anisopath/theta_path(50000),rad_path(50000),ntheta !la 3200->50000
c------------------------------------------------------------------------
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      external fqanis
      data ind/1,2,4/
c
c--num =1(only deltas); =2(delta and dX/dp); =4 (+time and tstar)
c--if ider is not 0, then num = 24 (4-24 are for anisotropic kernels)
      num=ind(itim)
c      print*, "num in integr.f", num

      if(ider.ne.0) num=4+nplay
      do i=1,num
        qints(i)=0.d0
      enddo
      if(rdep.eq.6371.d0) kdep=1
      ifl=-1
      iflag=0
      call qtau(xt(n),n,q1)
      if(q1.lt.0.d0)then 
c*** oops, ray already turned
        iflag=1
        return
      endif
      call qtau(xb(n),n,q1)
c*** find turning point -- assumes 1 per layer!
c ----dxt returns the derivative (dxtdp) of the turning radius w.r.t.
c ----parameter p. for fxtpp, see woodhouse & girnius, iii.14 
      if(q1.le.0.d0) then
        call zero1(xnorm,xt(n),xb(n),n,ierr)
        ifl=0
        iflag=2
        if(itim.ge.2) call dxt(xnorm,n,dxtdp,fxtpp)
      else
        xnorm=xb(n)
      end if
      nlev=nl(n)-1
      xdec=(xt(n)-xb(n))/nlev
cTEST
c	write(45,*)xdec
      xdecs=xdec+1.
      x1=xt(n)-xnorm
      sumd=0.
  10    x2=x1
c*** need fine spacing near center for deep turning rays
        if(xnorm.le.10.d0.and.x2.lt.xdecs) xdec=xnorm
        x1=x1-xdec
        if(x2.le.0.d0) go to 20
        if(x1.le.0.d0) x1=0.d0
c----------------------------------------------------
c**** temporary code to write kernels to unit 99
        if(ider.ne.0) then
          y=0.5d0*(x1+x2)+xnorm
          call qtau(y,n,q1)
          qq=sqrt(q1)
          call tder(y,qq)
c
c ---the following was used to output derivative in cubic polynomial
c         write(*,*) y,(tdif(k),k=1,5)
cc  900     format(6g14.6)
c----------------------------------------------------
        end if
        xx1=dsqrt(x1)
        xx2=dsqrt(x2)
cTEST
c	write(45,*)jpath+1,xx1,xx2
	
c ----gauslv is the routine which performs the gauss-legendre integration.
        call gauslv(xx1,xx2,idleg,fqanis,n,qints,num)
        if(xdec.eq.xnorm.and.x1.ne.0.d0) goto 10
c*** save increments in delta as a function of radius
        jpath=jpath+1
        delt(1,jpath,i12)=x1+xnorm
        delt(2,jpath,i12)=qints(1)-sumd
cTEST
c	write(43,*)jpath,qints(1)
c	write(44,"(a6,1x,i6.6,1x,2(f15.7,1x))")"delta:",jpath,delt(1,jpath,i12),delt(2,jpath,i12)
        sumd=qints(1)
c*** The following is a problematic part which produce wrong first layer,
c    YG, 1998 
c        if(x1+xnorm.gt.rdep) kdep=jpath
        if(x2+xnorm.ge.rdep) kdep=jpath
        goto 10
   20   continue
c ----------------------------------------------------------------------
      if(itim.eq.1) return
      if(xnorm.eq.0.d0) return
      if(ifl.ge.0) qints(2)=qints(2)-dxtdp*fxtpp/(dsqrt(xt(n)-xnorm))
      return
      end
