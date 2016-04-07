	call rsoinc(raydel,npt,idx)
	do i=1,npt
	   ztmp(i)=z(idx(i))
	enddo
	do i=1,npt
	   z(i)=ztmp(i)
	enddo	

c--find_path routine occasionally repeats points
	k=1
	raydtmp(1)=raydel(1)
	ztmp(1)=z(1)
	do i=2,npt
c	   if(raydel(i).ne.raydel(i-1))then
	   if(abs(raydel(i)-raydel(i-1)).gt.toler)then
	      k=k+1
	      raydtmp(k)=raydel(i)
	      ztmp(k)=z(i)
	   endif
	enddo
	npt=k
	do i=1,npt
	   raydel(i)=raydtmp(i)
	   z(i)=ztmp(i)
	enddo

	call gceq(eplon,eplat,stlon,stlat,raydel,npt,x,y,ndim)

	if(raydel(1).ne.0.)then
	   x(npt+2)=stlon
	   y(npt+2)=stlat
	   z(npt+2)=rearth
	   raydel(npt+2)=ddel
	   do i=npt+1,2,-1
	      y(i)=y(i-1)
	      x(i)=x(i-1)
	      z(i)=z(i-1)
	      raydel(i)=raydel(i-1)
	   enddo
	   x(1)=eplon
	   y(1)=eplat
	   z(1)=rearth-ddep
	   raydel(1)=0.d0
	   npt=npt+2
	else
	   x(npt+1)=stlon
	   y(npt+1)=stlat
	   z(npt+1)=rearth
	   raydel(npt+1)=ddel
	   x(1)=eplon
	   y(1)=eplat
	   z(1)=rearth-ddep
	   raydel(1)=0.d0
	   npt=npt+1
	endif

	do k=1,nint
	   call linint(y,yi,raydel,delint,ndim,npt)
	   call linint_x(x,xi,raydel,delint,ndim,npt)
	   call linint(z,zi,raydel,delint,ndim,npt)
	   do i=1,npt*2-1
	      if(xi(i).lt.0.)xi(i)=360.+xi(i)
	      if(xi(i).gt.360.)xi(i)=xi(i)-360.
	      x(i)=xi(i)
	      y(i)=yi(i)
	      z(i)=zi(i)
	      raydel(i)=delint(i)
	   enddo
	   npt=2*npt-1
	enddo
	
c---------------------project ray path onto grid
	do ila=1,nlay
	   if((z(1).le.rbnd(ila-1)).and.(z(1).ge.rbnd(ila)))then
	   iv0=ila
	   endif
        enddo
	ih0=isqre(y(1),x(1),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
	ind0=(iv0-1)*n1layer+ih0
	call span(ih0,ymi0,yma0,xmi0,xma0,nsqrs,nsqtot,nlatzones,eq_incr)
	x0=x(1)
	y0=y(1)
	z0=z(1)
	d0=0.

	do i=1,npt !loop over all points in ray path
	if(xma0.eq.360.)xma0=0.
	if(xmi0.eq.0.)xmi0=360.
c--determine index of voxel for i-th point on ray path
	   do ila=1,nlay
              if((z(i).le.rbnd(ila-1)).and.(z(i).ge.rbnd(ila)))then
                 iv=ila
              endif
           enddo
	   if(z(i).lt.rbnd(nlay))iv=nlay+1
	   ih=isqre(y(i),x(i),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
	   call span(ih,ymi,yma,xmi,xma,nsqrs,nsqtot,nlatzones,eq_incr)
	   ind=(iv-1)*n1layer+ih
	   if(ind.ne.ind0)then !crossed over to another voxel
	      if(iv.ne.iv0)then !vertical intersection
	            if(z(i).eq.z0)then
	               print*,j,"exception vert"
	               ind=ind0
	               goto11
	            endif
	         ivint=min(iv,iv0)
	         zint=rbnd(ivint)
		 xint=x0+(x(i)-x0)*(zint-z0)/(z(i)-z0)
		 yint=y0+(y(i)-y0)*(zint-z0)/(z(i)-z0)
		 dint=d0+(raydel(i)-d0)*(zint-z0)/(z(i)-z0)
c		 write(*,"(a8,1x,3(f12.6,1x))")"vertical",zint,xint,yint
	      else !horizontal intersection
		     if(dabs(ymi-yma0).lt.toler)then
	            if(y(i).eq.y0)then
	               print*,j,"exception s to n"
	               ind=ind0
	               goto11
	            endif
		    yint=ymi
		    xint=x0+(x(i)-x0)*(yint-y0)/(y(i)-y0)
		    zint=z0+(z(i)-z0)*(yint-y0)/(y(i)-y0)
		    dint=d0+(raydel(i)-d0)*(yint-y0)/(y(i)-y0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"n to s",zint,xint,yint!test
		 elseif(dabs(yma-ymi0).lt.toler)then
	            if(y(i).eq.y0)then
	               print*,j,"exception n to s"
	               ind=ind0
	               goto11
	            endif
		    yint=yma
		    xint=x0+(x(i)-x0)*(yint-y0)/(y(i)-y0)
		    zint=z0+(z(i)-z0)*(yint-y0)/(y(i)-y0)
		    dint=d0+(raydel(i)-d0)*(yint-y0)/(y(i)-y0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"s to n",zint,xint,yint!test
		 elseif(dabs(xmi-xma0).lt.toler)then
	            if(x(i).eq.x0)then
	               print*,j,"exception e to w"
	               ind=ind0
	               goto11
	            endif
	            xint=xmi
		    yint=y0+(y(i)-y0)*(xint-x0)/(x(i)-x0)
		    zint=z0+(z(i)-z0)*(xint-x0)/(x(i)-x0)
		    dint=d0+(raydel(i)-d0)*(xint-x0)/(x(i)-x0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"e to w",zint,xint,yint!test
		 elseif(dabs(xma-xmi0).lt.toler)then
	            if(x(i).eq.x0)then
	               print*,j,"exception w to e"
	               ind=ind0
	               goto11
	            endif
	            xint=xma
		    yint=y0+(y(i)-y0)*(xint-x0)/(x(i)-x0)
		    zint=z0+(z(i)-z0)*(xint-x0)/(x(i)-x0)
		    dint=d0+(raydel(i)-d0)*(xint-x0)/(x(i)-x0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"w to e",zint,xint,yint!test
		 else
		    print*,"problem",j
		    print*,xmi0,xma0,ymi0,yma0
		    print*,xmi,xma,ymi,yma
		    stop "points on ray are too far"
		 endif
	      endif

	      zav=(zint+z0)/2.
	      ds=((dint-d0)/radian)*zav

c------------------------------------------this needs to be redone
	      dqdvshout=dqdvshint-dqdvsh0
	      dqdvsvout=dqdvsvint-dqdvsv0
	      dqdvphout=dqdvphint-dqdvph0
	      dqdvpvout=dqdvpvint-dqdvpv0

	      if(ds.ne.0..and.zint.ge.rcmb.and.z0.ge.rcmb)then
	         do l=1,nstartmod-1
	            if(radmod(l).ge.zav.and.zav.ge.radmod(l+1))then
	               a=(vstart(l+1)-vstart(l))/(radmod(l+1)-radmod(l))
	               b= vstart(l)-a*radmod(l)
		       vref=a*zav+b
		       goto33
		    endif
	         enddo
33	         continue

c------------------------------------increment matrix (entries and indices)
c	         irec=irec+1	
c	         write(111,rec=irec)sngl(-ds/vref)
c	         write(112,rec=irec)ind0
	         irec=irec+1
	         write(111,rec=irec)sngl(dqdvshout)
	         write(112,rec=irec)ind0
	         irec=irec+1
	         write(111,rec=irec)sngl(dqdvsvout)
	         write(112,rec=irec)ind0+nvx
c------------------------------------here easy to add VPH, VPV part!
cTEST	
c	write(95,*)sngl(zav),sngl(vref)
c	write(99,"(5(f12.6,1x))")sngl(zav),sngl(-ds/vref),sngl(dqdvshout),sngl(dqdvsvout),sngl(dqdvshout+dqdvsvout)
c	write(49,"(5(f12.5,1x))")sngl(zav),sngl(-ds/vref),sngl(dqdvphout),sngl(dqdvpvout),sngl(dqdvphout+dqdvpvout)
	
	      endif

	      x0=xint !update x0,y0,z0 
	      y0=yint
	      z0=zint
	      d0=dint
	      iv0=iv !update other variables
	      ih0=ih !not strictly needed
	      ind0=ind 
	      xmi0=xmi
	      xma0=xma
	      ymi0=ymi
	      yma0=yma
	   endif !executed if intersection
11	   continue
	enddo !end of loop over ray path

c------------------------------------receiver voxel
	if(raydel(npt).ne.d0)then
	   zav=(z(npt)+z0)/2.
	   ds=((raydel(npt)-d0)/radian)*zav


c------------------------------------------this needs to be redone
	   dqdvshout=dqdvsh(npt)-dqdvsh0
	   dqdvsvout=dqdvsv(npt)-dqdvsv0
	   dqdvphout=dqdvph(npt)-dqdvph0
	   dqdvpvout=dqdvpv(npt)-dqdvpv0

c	   irec=irec+1 ! why was this up here?
	   if(ds.ne.0.)then ! it has to be nonzero anyway at this point
	      do l=1,nstartmod-1
	         if(radmod(l).ge.zav.and.zav.ge.radmod(l+1))then
	            a=(vstart(l+1)-vstart(l))/(radmod(l+1)-radmod(l))
	            b= vstart(l)-a*radmod(l)
		    vref=a*zav+b
		    goto34
		 endif
	      enddo
34	      continue
c	      write(111,rec=irec)sngl(-ds/vref)
c	      write(112,rec=irec)ind
	      irec=irec+1
	      write(111,rec=irec)sngl(dqdvshout)
	      write(112,rec=irec)ind
	      irec=irec+1
	      write(111,rec=irec)sngl(dqdvsvout)
	      write(112,rec=irec)ind+nvx
cTEST
c	write(99,"(5(f12.6,1x))")sngl(zav),sngl(-ds/vref),sngl(dqdvshout),sngl(dqdvsvout),sngl(dqdvshout+dqdvsvout)
c	write(49,"(5(f12.5,1x))")sngl(zav),sngl(-ds/vref),sngl(dqdvphout),sngl(dqdvpvout),sngl(dqdvphout+dqdvpvout)
	   endif
	endif
