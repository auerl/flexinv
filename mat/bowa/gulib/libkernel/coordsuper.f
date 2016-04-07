c 
c given cell index on the surface, finds lon and lat of its center
c
	subroutine coordsuper(nbloc,blocla,bloclo,nsqrs,nlatzones,eq_incr)
        parameter(nlatzomax=180.)
	dimension nsqrs(nlatzomax)
	ntot=0
        ! Loop(s) over all the blocks
	do 500 ila=1,nlatzones
           ! Increment latitude
	   rlati=90.-(eq_incr*(ila-1))
           ! Calculate increment in longitude for this band
	   rinlo=(360./nsqrs(ila))
	   do 400 isq=1,nsqrs(ila)
	      rlong=(360./nsqrs(ila))*(isq-1)
	      ntot=ntot+1
	      if(ntot.eq.nbloc)then
	         bloclo=rlong+(rinlo/2.)
	         blocla=rlati-(eq_incr/2.)
	         goto 600
	      endif
400	   continue
500	continue
600	return
	end ! end of subroutine coordsuper
