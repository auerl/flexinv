      subroutine prray_kernel(ray,pp,qvec,qray,xtu,iopt)
c*** prints out results for this ray type and writes to unit 9
c*** qvec(1)=delta, qvec(2)=dX/dp, qvec(3)=time, qvec(4)=tstar
      implicit real*8(a-h,o-z)
      character*(*) ray
      dimension qray(4,5,20),qvec(*)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      data if1/1/
      if(iopt.eq.1) then
        print 930
        write(9,930)
  930   format(3x,'Ray',8x,'Distance',5x,'Time',10x,'t*',12x,'dXdp',9x,
     +   'p',5x,'Bot. rad')
        if1=0
      end if
      qdel=qvec(1)
      qddp=qvec(2)
      qtim=qvec(3)
      qtst=qvec(4)
      if(iopt.eq.1) then
	  write(9,920)ray(1:12),qdel,qtim,qtst,qddp,pp,xtu
      endif
c      write(*,921)ray(1:12),qdel,qtim,qddp,pp,xtu
      if(iopt.eq.1) then 
      	if(ider.ne.0)then
       	 	do 170 l=1,numlyr
          		write(9,900)l 
          		do 170 k=1,4
  170         			write(9,910) k,(qray(k,j,l),j=1,5)
      	end if
      endif
  900 format(3x,28("-"),"   l a y e r   ",i2,3x,28("-")
     +  /,3x,"k",8x,"vph",12x,"vpv",12x,"vsh",12x,"vsv",12x,"eta")
  910 format(2x,i2,5(3x,e12.6))
  920 format(1x,a12,3(f8.3,4x),f12.4,4x,f7.4,f10.3)
  921 format(1x,a12,2(f8.3,4x),f12.4,4x,f7.4,f10.3)
      return
      end
