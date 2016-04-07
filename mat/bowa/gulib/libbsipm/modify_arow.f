      subroutine modify_arow(arow, indarow, nonzero, natd1, iopt)
      include "../aniso_spline_parallel.h"
c
c  This sub can compute a new A matrix from a given one.
c  purpose is to construct a new set of inversion parameters 
c  rather than the original.  Currently this implement a couple
c  of options only.
c
c  Input:
c	arow     ---- a row of original A matrix
c	indarow  ---- index of non-zero elements
c	nonzero  ---- number of non-zero elements
c	natd1     ---- expected number of element in each parameter
c	iopt     ---- input that controls the operation
c		      1= doe nothing,  2=invert for V_voigt+(Vsv-Vsh),  
c		      3= invert for V_iso+(Vsv-Vsh) 
c
c  Output:
c	arow  ---- the new row will be saved constructed and saved
c		   in the same arrow if iopt is non 1
c	indarow  ---- index for the new non-zero elements
c	nonzero  ---- save the new number of nonzero elements
c 
c
c does nothing
	dimension arow(1), indarow(1)
	integer   nonzero, natd1, iopt
	dimension temp1(mxparm), temp2(mxparm), temp3(mxparm)
	dimension indtopo(mxparm)
	integer   ntoponz
	if(iopt.eq.1) goto 100
	do i=1, mxparm
		temp1(i)=0.0
		temp2(i)=0.0
		temp3(i)=0.0
		indtopo(mxparm)=0
	enddo
	natd2=natd1*2
	ntoponz=0
	do i=1, nonzero
		j=indarow(i)
		if(j.le.natd1) then
c first parameter from original A mat
			temp1(j)=arow(i)
		else if(j.le.natd2) then
c second parameter from original A mat
			temp2(j-natd1)=arow(i)
		else
c save topography portion
			ntoponz=ntoponz+1
			temp3(j)=arow(i)
			indtopo(ntoponz)=j
		endif
	enddo
	if(iopt.eq.2) then
c
c  assume SH is the first element!
c  construct  Viso and Vdiff
c
		do i=1, natd1
			i2=i+natd1
			if(abs(temp1(i)).gt.0.0000001) then
				temp3(i)=temp3(i)+temp1(i)
				temp3(i2)=temp3(i2)-temp1(i)*0.66667
			endif
			if(abs(temp2(i)).gt.0.0000001) then
				temp3(i)=temp3(i)+temp2(i)
				temp3(i2)=temp3(i2)+temp2(i)*0.33333
			endif
		enddo
      	else if(iopt.eq.3) then
c
c regular average and Vdiff
c
		do i=1, natd1
			i2=i+natd1
			if(abs(temp1(i)).gt.0.0000001) then
				temp3(i)=temp3(i)+temp1(i)
				temp3(i2)=temp3(i2)-temp1(i)*0.5
			endif
			if(abs(temp2(i)).gt.0.0000001) then
				temp3(i)=temp3(i)+temp2(i)
				temp3(i2)=temp3(i2)+temp2(i)*0.5
			endif
		enddo
	endif
	j=0
	do i=1, natd2
		if(abs(temp3(i)).gt.0.0000001) then
			j=j+1
			arow(j)=temp3(i)
			indarow(j)=i
		endif
	enddo
c
c put back the topography portion
c
	do i=1, ntoponz
		j=j+1		! j continues from the above
		k=indtopo(i)
		arow(j)=temp3(k)
		indarow(j)=k
	enddo
	nonzero=j
100     return
      end
