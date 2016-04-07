      subroutine zero(z,a1,b1,f,ftarg)
c  finds a root (a=z) of  f(a)-ftarg=0 with z
c  between a1,b1 by a combination of bisection and lin interp.
      implicit double precision(a-h,o-z)
      data re/1.d-10/
      a=a1
      fa=f(a)-ftarg
      if(fa.eq.0.d0) then
         z=a
         return
      end if
      b=b1
      fb=f(b)-ftarg
      if(fb.eq.0.d0) then
         z=b
         return
      end if
c*** return if no zero (or not monotonic??)
      if(fa*fb.ge.0.d0) then
        z=-1.d0
        return
      end if
      c=a
      fc=fa
      s=c
      fs=fc
c ==============================================================
   10 h=0.5d0*(b+c)
      t=dabs(h*re)
      if(dabs(h-b).le.t) then
        z=h
        return
      end if 
      if(dabs(fb).gt.dabs(fc)) then
        y=b
        fy=fb
        g=b
        fg=fb
        s=c
        fs=fc
      else
        y=s
        fy=fs
        g=c
        fg=fc
        s=b
        fs=fb
      end if
      if(fy.eq.fs) then
        b=h
      else
        e=(s*fy-y*fs)/(fy-fs)
        if(dabs(e-s).le.t) e=s+dsign(t,g-s)
        if((e-h)*(s-e).lt.0.d0) then
          b=h
        else
          b=e
        end if
      end if
      fb=f(b)-ftarg
      if(fg*fb.ge.0.d0) then
        c=s
        fc=fs
      else
        c=g
        fc=fg
      end if
      goto 10
c =============================================================
      end

