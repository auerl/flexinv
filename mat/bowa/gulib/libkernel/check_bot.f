	subroutine check_bot(k,rad,iflag)
c this will check if a point is the bottoming point of the ray.
c it uses the assumption that   rayparm = r/v
c INPUT:	p ----- ray parameter
c		k ----- ray type: 1 = SH;  2 = P;  3 = SV
c	      rad ----- radius
c OUTPUT:
c	    iflag ----- 1 = yes;  0 = no
        implicit real*8(a-h,o-z)
	common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
	common/layr/nl(20),xb(20),xt(20),ifanis,nplay
	dimension val(5)
	iflag = 0
	if(rad.gt.3479.9.and.rad.lt.3480.1) rad = 3480.05
	call xmodelvel(rad, val)
	if(k.eq.1) then
		v = val(3)
	else if(k.eq.3) then
		v = val(4)
	else
		v= val(2)
	endif
	rv = rad/v
c use a cursion of 2.0 for inaccuracy in the calculation,
c subjected to change, YG, 1998.
	if(rv.gt.p-2.d0.and.rv.lt.p+2.d0) then
		iflag = 1
	endif
	return
	end		
