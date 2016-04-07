       subroutine handle_2mod(coef1,coef2,ngpt,nrad,iop,coefo1,coefo2)
c  manipulate two models and the operations depends on 
c  input option "iopt":
c  1. add two B-spline models
c  2. subtract two B-spline models
c  3. For a voigt_avg & Vsv-Vsh, recovers Vsh & Vsv
c  4. For a simple_avg & Vsv-Vsh, recovers Vsh & Vsv
c  for options 1 and 2, only coefo1 is outputed coefo2 is dummy.
c  for other two, coefo1=Vsh and coefo2=Vsv
c	
      dimension coef1(1), coef2(1), coefo1(1), coefo2(1)
      integer  ngpt, nrad, iop
      
      ncoef=ngpt*nrad
      print*, 'total number of coefficients= ', ncoef 
      if(iop.eq.1) then
c add models
	do i=1, ncoef
		coefo1(i)=coef1(i)+coef2(i)
	enddo
      endif
      if(iop.eq.2) then
c subtract 2 from 1
	do i=1, ncoef
		coefo1(i)=coef1(i)-coef2(i)
	enddo
      endif
      if(iop.eq.3) then
c solve for Vsv and Vsh from given Viso (coef1) and Vsv-Vsh (coef2)
	do i=1, ncoef
		coefo1(i)=coef1(i)-0.666666666*coef2(i)
		coefo2(i)=coef1(i)+0.333333333*coef2(i)
	enddo
      endif
      if(iop.eq.4) then
c solve for Vsv and Vsh from given Viso (coef1) and Vsv-Vsh (coef2)
c simple average case
	do i=1, ncoef
		coefo1(i)=coef1(i)-0.5*coef2(i)
		coefo2(i)=coef1(i)+0.5*coef2(i)
	enddo
      endif
      if(iop.gt.4.or.iop.lt.1) stop 'option not found!'
      return
      end 
