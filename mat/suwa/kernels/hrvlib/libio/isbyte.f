	subroutine isbyte(k,ibuff,j)
c		needed in subroutine trnslt.
	dimension ibuff(*),mask(4)
 	data mask/z'00ffffff',z'ff00ffff',z'ffff00ff',z'ffffff00'/
	iw=1+j/4	
	ib=1+j-4*(iw-1)
   	ii=and(ibuff(iw),mask(ib))
 	ibuff(iw)=or(ii,ishft(k,8*(4-ib)))
 	return
	end
