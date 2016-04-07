      subroutine addarow(res,premtt,res2,nonzero,indarow,arow,
     +	 dc400,dc670,weight,num3d,ngpt,ata,atd,ifdcshift)
c
c... This computes ATA from row of A-matrix
c... Input:	res    ---  residual
c... 		premtt ---  predicted PREM time
c... 		res2   ---  predicted time residual from starting model
c... 		nonzero  ---  number of nonzero A elements
c... 		indarow  ---  array of index of nonzero A elements
c... 		arow  ---  a row of A matrix
c...          The factors below are average coefficients for a spline model
c...	      of a simple DC shift
c... 		dc400 ---  shift factor to the average depth of 400 (e.g., 10 km).
c... 		dc670 ---  shift factor to the average depth of 670 (e.g., -20 km)
c... 		weight ---  weighting factor of this A matrix
c... 		num3d  ---  number of 3D elements
c... 		ngpt   ---  number of horizontal splines
c... 		ifdcshift ---  if correct for average shift in discontinuity depth
c...Output:
c... 		ata, atd    ---  ATA and ATD for inversions,lower triangle
c
      integer indarow(1)
      real arow(1),ata(1),atd(1)
c
c  the following removes perturbation to the average of discontinuities
c  at 400 and 670 km
c

      residual=res+res2
      corr400=0.0
      corr670=0.0
      if(ifdcshift.ne.0) then
	 itopo1=num3d+ngpt
      	 do i=1, nonzero
	    ia=indarow(i)
	    if(ia.gt.num3d.and.ia.le.itopo1) then
		corr400=corr400+arow(i)*dc400
	    endif
	    if(ia.gt.itopo1) then
		corr670=corr670+arow(i)*dc670
	    endif
	 enddo
      endif
      residual=residual-corr400-corr670
      do i=1,nonzero
	 ia=indarow(i)
	 kk=ia*(ia-1)/2+1
	 arowi=arow(i)
	 atd(ia)=atd(ia)+arowi*residual
	 do j=1,nonzero
	 	if(j.le.i) then
			iata=kk+indarow(j)-1
			ata(iata)=ata(iata)+arowi*arow(j)
		endif
	 enddo
      enddo
      return
      end
