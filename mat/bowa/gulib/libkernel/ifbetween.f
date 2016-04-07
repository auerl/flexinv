      integer function ifbetween(a, x, b)
      implicit double precision (a-h, o-z)
      ifbetween=0
      if((a.le.x).and.(b.gt.x)) then
	ifbetween=1
      endif
      if((b.lt.x).and.(a.ge.x)) then
	ifbetween=1
      endif
      return
      end
