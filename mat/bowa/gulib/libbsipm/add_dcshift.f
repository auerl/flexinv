      subroutine add_dcshift(a,inda,n,num3d,ngpt,dc400,dc670,dotp,iop)
c
c... adds dc shifts of average discontinuity depth to residuals
c
      dimension a(1),coef(1)
      integer inda(1)
c
      dotp=0.0
      itopo=ngpt+num3d
      do i=1, n
	 ind=inda(i)
	 if(ind.gt.num3d.and.ind.le.itopo) then
		if(iop.eq.1) then
			dotp=dotp+a(i)*dc400
		else
			dotp=dotp-a(i)*dc400
		endif
	 endif
	 if(ind.gt.itopo) then
		if(iop.eq.1) then	
			dotp=dotp+a(i)*dc670
		else
			dotp=dotp-a(i)*dc670
		endif
	 endif
      enddo
      return
      end
