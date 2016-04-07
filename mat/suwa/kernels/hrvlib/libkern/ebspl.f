      function ebspl(x,xarr,np,splcof)
      dimension xarr(np),splcof(np)
c
c---- iflag=1 ==>> second derivative is 0 at end points
c---- iflag=0 ==>> first derivative is 0 at end points
c
      iflag=1
c
c---- first, find out within which interval x falls
c
      interval=0
      ik=1
      do while(interval.eq.0.and.ik.lt.np)
        ik=ik+1
        if(x.ge.xarr(ik-1).and.x.le.xarr(ik)) interval=ik-1
      enddo
      if(x.gt.xarr(np)) interval=np
      if(interval.eq.0) then
cc        write(6,"('low value:',2f10.3)") x,xarr(1)
      else if(interval.gt.0.and.interval.lt.np) then
cc        write(6,"('bracket:',i5,3f10.3)") interval,xarr(interval),x,
cc     &xarr(interval+1)
      else
cc        write(6,"('high value:',2f10.3)") xarr(np),x
      endif
      value=0.
      do ib=1,np
        val=0.
        if(ib.eq.1) then
          r1=(x-xarr(1))/(xarr(2)-xarr(1))
          r2=(xarr(3)-x)/(xarr(3)-xarr(1))
          r4=(xarr(2)-x)/(xarr(2)-xarr(1))
          r5=(x-xarr(1))/(xarr(2)-xarr(1))
          r6=(xarr(3)-x)/(xarr(3)-xarr(1))
         r10=(xarr(2)-x)/(xarr(2)-xarr(1))
         r11=(x-xarr(1))  /(xarr(2)-xarr(1))
         r12=(xarr(3)-x)/(xarr(3)-xarr(2))
         r13=(xarr(2)-x)/(xarr(2)-xarr(1))
          if(interval.eq.ib.or.interval.eq.0) then
               if(iflag.eq.0) val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11 
     &+r13**3
               if(iflag.eq.1) val=0.6667*(r1*r4*r10 + r2*r5*r10 + 
     &r2*r6*r11 + 1.5*r13**3)
          else if(interval.eq.ib+1) then
               if(iflag.eq.0) val=r2*r6*r12
               if(iflag.eq.1) val=0.6667*r2*r6*r12
          else
            val=0.
          endif
        else if(ib.eq.2) then
          rr1=(x-xarr(1))/(xarr(2)-xarr(1))
          rr2=(xarr(3)-x)/(xarr(3)-xarr(1))
          rr4=(xarr(2)-x)/(xarr(2)-xarr(1))
          rr5=(x-xarr(1))/(xarr(2)-xarr(1))
          rr6=(xarr(3)-x)/(xarr(3)-xarr(1))
         rr10=(xarr(2)-x)/(xarr(2)-xarr(1))
         rr11=(x-xarr(1))  /(xarr(2)-xarr(1))
         rr12=(xarr(3)-x)/(xarr(3)-xarr(2))
          r1=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
          r3=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
         r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
          if(interval.eq.ib-1.or.interval.eq.0) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
               if(iflag.eq.1) val=val+0.3333*(rr1*rr4*rr10 + 
     &rr2*rr5*rr10 + rr2*rr6*rr11)
          else if(interval.eq.ib) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
               if(iflag.eq.1) val=val+0.3333*rr2*rr6*rr12
          else if(interval.eq.ib+1) then
               val=r2*r6*r12
          else
               val=0.
          endif
        else if(ib.eq.np-1) then
          rr1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          rr2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          rr3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          rr4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          rr5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          rr7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
          rr8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
          rr9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
c
          r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
          r2=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
          r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
          if(interval.eq.ib-2) then
               val=r1*r3*r7
          else if(interval.eq.ib-1) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
               if(iflag.eq.1) val=val+0.3333*rr1*rr3*rr7
          else if(interval.eq.ib.or.interval.eq.np) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
               if(iflag.eq.1) val=val+0.3333*(rr1*rr3*rr8 + 
     &rr1*rr4*rr9 + rr2*rr5*rr9)
          else
            val=0.
          endif
        else if(ib.eq.np) then
          r1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          r2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          r3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          r4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          r5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          r7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
          r8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
          r9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          r13=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          if(interval.eq.np-2) then
               if(iflag.eq.0) val=r1*r3*r7
               if(iflag.eq.1) val=0.6667*r1*r3*r7
          else if(interval.eq.np-1.or.interval.eq.np) then
               if(iflag.eq.0) val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
     &+ r13**3
               if(iflag.eq.1) val=0.6667*(r1*r3*r8 + r1*r4*r9 + 
     &r2*r5*r9 + 1.5*r13**3)
          else
            val=0.
          endif
        else
          r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
          r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
          r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
          r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
         r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
          if(interval.eq.ib-2) then
               val=r1*r3*r7
          else if(interval.eq.ib-1) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
          else if(interval.eq.ib) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
          else if(interval.eq.ib+1) then
               val=r2*r6*r12
          else
            val=0.
          endif
        endif
        value=value+val*splcof(ib)
      enddo
      ebspl=value
      return
      end

