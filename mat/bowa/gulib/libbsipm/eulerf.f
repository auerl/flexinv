      subroutine neweuler(sth,sph,rth,rph,alpha,beta,gamma,del,pth,pph)
c     Finds the Euler angles alpha,beta,gamma that rotate the coordinate axes
c     so that the z-axis is at the pole of the source-receiver great circle
c     (s x r), and the x-axis is at the source. See Edmonds' Angular Momentum
c     in Quantum Mechanics, page 7 for the angle conventions.
c     input: sth,sph = source coordinates in radians
c     rth,rph = receiver coordinates in radians
c     output: alpha,beta,gamma = euler angles in radians which rotate the
c     original coordinate system to the one with the source-receiver
c     great circle on the equator, the source at (PI/2,0). The minor
c     arc to the receiver is in the positive phi direction.
c     del = source-receiver separation in radians.
c     pth,pph = source-receiver great circle pole location.
c     
      data pi/3.14159265358979/, hpi/1.5707963705063/
c     Get cartesian coordinates for source and receiver
      call cart(sth,sph,sx,sy,sz)
      call cart(rth,rph,rx,ry,rz)
c     { epicentral dist
      del = acos(sx*rx + sy*ry + sz*rz)
      call cross(sx,sy,sz,rx,ry,rz,px,py,pz)
      pth = atan2(sqrt(px*px+py*py),pz)
      if(px.eq.0. .and. py.eq.0.) then
c     special case of pole at z or -z
         pph = 0.
      else
         pph = atan2(py,px)
      endif
      alpha = pph
      beta = pth
c     the x'' axis (call it t) is at pth+pi/2,pph
      ttheta = pth + hpi
      call cart(ttheta,pph,tx,ty,tz)
c     the third Euler angle, gamma, rotates x'' to the source s.
      gamma = acos(sx*tx + sy*ty + sz*tz)
c     form q = x'' x s to check the sign of gamma (q/|q| = +-p/|p|)
      call cross(tx,ty,tz,sx,sy,sz,qx,qy,qz)
      sgn = px*qx + py*qy + pz*qz
      if(sgn .lt. 0.) gamma = -gamma
      return
      end
c     
c-----------------------------------------------------------------------
c     
      subroutine cart(thet,phi,x,y,z)
      
      s = sin(thet)
      x = s*cos(phi)
      y = s*sin(phi)
      z = cos(thet)
      
      return
      end
c     
c-----------------------------------------------------------------------
c     
      subroutine cross(sx,sy,sz,rx,ry,rz,px,py,pz)
      px = sy*rz - sz*ry
      py = sz*rx - sx*rz
      pz = sx*ry - sy*rx
      return
      end
