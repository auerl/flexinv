	function isqre_orig(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
c----finds the index of the square where (lat,lon) is
	implicit real*8 (a-h,o-z)
	real*8 lat,lon,loc_incr
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	lazone=(90.-lat)/eq_incr+1
	if((90.-lat).gt.180.)lazone=nlatzones
	if((90.-lat).gt.181.)stop "problems in function isqre"
	if(lazone.gt.nlatzones)then
	   print*,"problems in function isqre, latitude",lazone,lat
	   stop
	endif
	if(lon.lt.0.)lon=360.+lon
	loc_incr=360./float(nsqrs(lazone))
	isqre=(lon/loc_incr)+1
	isqre=isqre+nsqtot(lazone)
	if(isqre.gt.n)then
	   print*,"problems in function isqre, longitude",n
	   stop
	endif
	return
	end
c 
c finds the number of the square where (xlat,xlon) is *************
c
	function isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n,eq_incr)
	implicit real*8 (a-h,o-z)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
      lazone=(90.-xlat)/eq_incr+1
	if(lazone.gt.nlatzones)lazone=nlatzones
	isqre=(xlon/360.)*nsqrs(lazone)+1
	isqre=isqre+nsqtot(lazone)
	return
	end ! end of subroutine isqre
