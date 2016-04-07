      double precision function surfreflker(rparm,rdis,rname)
      implicit double precision (a-h,o-z)
c..  dt due to dr Morrelli and Dziewonski for reflected waves at a bondary
c..        dr=dt*r/(2*(eta**2-(rparm*57.3)**2)**(0.5))
c..  note I took the negative of the above due to depth, not radius
c

c *** basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    


      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/coeff/coef(4,8,20)
      common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso
      common/secphases/isectype,rsec,indexsec,lsec
      character*(*) rname
c
c... if not precursor, do nothing.
      if(isectype.ne.1) return
      r=rdis-0.1	! underside reflection, use underside velocity
      y=r/rnorm
      do j=1, numlyr
	      if(r.gt.xb(j).and.r.le.xt(j)) then
		      iq = j
		      ierror = 0
		      goto 999
	      endif
      enddo
999   continue
      call getmod_iso(y,iq,vpeq,vseq)
      j=ichdec(rname(1:1),k)
	
      if(k.eq.1) then
c... PP precursors
	 v=vpeq
      else
c... SS precursors
	 v=vseq
      endif
      eta=rdis/v
      fact=2.0/rdis
      surfreflker=-fact*sqrt(eta*eta-rparm*rparm)
      return
      end
