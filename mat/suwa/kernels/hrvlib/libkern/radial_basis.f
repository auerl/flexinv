      subroutine radial_basis(radius4,type,ikern,value4,istatus)
      character*1 type
      dimension splpts(30)
      dimension splcof(30)
c
      if(type.eq.'A') then
        nspl=20
        splpts(1)=0.
        splpts(2)=50.
        splpts(3)=100.
        splpts(4)=150.
        splpts(5)=200.
        splpts(6)=250.
        splpts(7)=300.
        splpts(8)=400.
        splpts(9)=500.
        splpts(10)=600.
        splpts(11)=700.
        splpts(12)=800.
        splpts(13)=1000.
        splpts(14)=1250.
        splpts(15)=1500.
        splpts(16)=1750.
        splpts(17)=2000.
        splpts(18)=2250.
        splpts(19)=2500.
        splpts(20)=2891.
c
        depth=6371.-radius4
c
        do ik=1,nspl
          splcof(ik)=0.
          if(ik.eq.ikern) splcof(ik)=1.
        enddo
        value4=ebspl(depth,splpts,nspl,splcof)
      else
        istatus=999
      endif
      return
      end
