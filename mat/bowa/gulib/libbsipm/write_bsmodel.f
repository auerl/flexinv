	subroutine write_bsmodel(io,filenm,npt,nknot,xlat,xlon,nnb,near,xknot,
     +				istart,coef)
c   this subroutine writes out a B-spline model.  It automatically detects if
c   the input radial parameterization is continuous or split model and treats
c   each accordingly, JG2002.
c
	include "../aniso_spline_parallel.h"
	dimension xlat(1),xlon(1),xknot(1),coef(1)
	dimension nnb(1), near(mxleny,1)
c	integer isplayer(mxstrp),isplayer
	integer isplayer(mxstrp)
	character*(*) filenm
	logical   exists
       
c	inquire(file=bsfile(1:lnblnk(bsfile)),exist=exists)
	open(io, file=filenm,status='unknown')
c
c... The following finds the spliting depths
c... Assumption is that there is no "Hole" between spliting depths
	ifreq=sqrt(real((npt-2.0)/10.0))+0.45
	if(nknot.eq.1) then
c...  topography map
		isplit=0
		write(*,"('output topography map...')")		
	else
		call find_rsplit(nknot,xknot,isplayer,isplit)		
	endif
c----------------------------------------------------------------------
	if(isplit.eq.0) then
c... continuous model parameterization or topography model
		write(io, "(i6, i4)") npt, nknot		
	else if (isplit.eq.1) then
		write(io, "('UPPERLOWERM', i6, i4, i4)") npt, isplayer(1), nknot-isplayer(1)
	else 
		write(*, "('Multiple-split model not supported yet...')")
	endif
	do i=1, npt
		write(io,"(i6,2f12.4,15i6)") i,xlat(i),xlon(i),nnb(i),
     #		(near(i,j), j=1,ifreq)
	enddo
	kk=istart
	nlines=npt/5
c write model coefficients
	iknot=1
	do ii=1, nknot
		if(nknot.ne.1) write(io,*) iknot, xknot(ii) ! for non-topography models
		do jj=1, nlines
			write(io,"(5f12.7)") (coef(l), l=kk,kk+4)
			kk=kk+5
		enddo
		write(io,"(2f12.7)") coef(kk), coef(kk+1)
		kk=kk+2
		iknot=iknot+1
	enddo
        close(io)
	return
	end

