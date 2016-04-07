      subroutine radial_basis(radius4,type,ikern,value4,istatus,layers)
      character*1 type
      dimension splpts(50),varlr(50),splcof(50)
      dimension laypts(30)
      dimension layers(30)

c *** 50 equally thick layers
      if(type.eq.'C') then ! constant thickness layers
        dr=96. !corresponding to 50 layers
        rtop=6371.-(float(ikern-1)*dr)
        rbot=6371.-(float(ikern)*dr)
    	    if(radius4.le.rtop.and.radius4.gt.rbot)then
	      value4=1.
	    else
	      value4=0.
	    endif
c *** Read vertical parametrization from layers.in
      elseif(type.eq.'F') then
        ! Top and bottom of current layer
        rtop=6371.-dble(layers(ikern))
        rbot=6371.-dble(layers(ikern+1))
        ! Step basis functions
        if (radius4.le.rtop.and.radius4.gt.rbot) then
           value4=1.
        else
           value4=0.
        endif
c *** 26 variable thick layers
      elseif(type.eq.'E') then
        ! Layer definition, 27 laypts = 26 layers
        laypts(1)=0.
        laypts(2)=50.
        laypts(3)=100.
        laypts(4)=150.
        laypts(5)=225.
        laypts(6)=300.
        laypts(7)=375.
        laypts(8)=450.
        laypts(9)=550.
        laypts(10)=650.
        laypts(11)=750.
        laypts(12)=850.
        laypts(13)=975.
        laypts(14)=1100.
        laypts(15)=1225.
        laypts(16)=1350.
        laypts(17)=1500.
        laypts(18)=1650.
        laypts(19)=1800.
        laypts(20)=1950.
        laypts(21)=2100.
        laypts(22)=2250.
        laypts(23)=2375.
        laypts(24)=2500.
        laypts(25)=2625.
        laypts(26)=2750.
        laypts(27)=2891.        
        ! Top and bottom of current layer
        rtop=6371.-dble(laypts(ikern))
        rbot=6371.-dble(laypts(ikern+1))
        ! Step basis functions
        if (radius4.le.rtop.and.radius4.gt.rbot) then
           value4=1.
        else
           value4=0.
        endif
c *** 19 variable thick layers, automatically increasing d
      elseif(type.eq.'D') then ! variable thickness layers
        nlatop=10 ! number of thin layers at top of mantle
        drtop=25. ! their thickness
        nlabum=5  ! number of thicker layers in mid upper mantle
        drum=50.  ! their thickness
        ttop=50.  ! thickness of top layer
        nspl=nlatop+nlabum+4
        rearth=6371.
        varlr(1)=rearth
        varlr(2)=rearth-ttop ! top layer masked by crust
        do ik=3,nlatop+2
           varlr(ik)=varlr(ik-1)-drtop
        enddo
        do ik=nlatop+3,nlatop+nlabum+2
           varlr(ik)=varlr(ik-1)-drum
        enddo
        ik=nlatop+nlabum+2
        varlr(ik+1)=rearth-660.
        if(varlr(ik+1).ge.varlr(nlatop+nlabum+2))then
           print*,"too many layers in um",ik,varlr(ik+1),varlr(ik)
           stop
        endif
        varlr(ik+2)=rearth-1300.
        varlr(ik+3)=rearth-3000.
        ! find in which layer we are and assign basis function value
        rtop=varlr(ikern)
        rbot=varlr(ikern+1)
        do ik=1,nspl
	     if(radius4.le.rtop.and.radius4.gt.rbot)then
	       value4=1.
	     else
	       value4=0.
	     endif
        enddo
c *** 19 spline knots, untested and probably obsolete
      elseif(type.eq.'A') then
        nspl=19
        splpts(1)=0.
        splpts(2)=60.
        splpts(3)=120.
        splpts(4)=180.
        splpts(5)=240.
        splpts(6)=300.
        splpts(7)=450.
        splpts(8)=600.
        splpts(9)=750.
        splpts(10)=900.
        splpts(11)=1050.
        splpts(12)=1200.
        splpts(13)=1350.
        splpts(14)=1500.
        splpts(15)=1750.
        splpts(16)=2000.
        splpts(17)=2250.
        splpts(18)=2500.
        splpts(19)=2891.
        depth=6371.-radius4
        do ik=1,nspl
          splcof(ik)=0.
          if(ik.eq.ikern) splcof(ik)=1.
        enddo
        value4=ebspl(depth,splpts,nspl,splcof)
c *** 20 spline knots, untested and probably obsolete
      elseif(type.eq.'B') then
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
        depth=6371.-radius4
        do ik=1,nspl
          splcof(ik)=0.
          if(ik.eq.ikern) splcof(ik)=1.
        enddo
        value4=ebspl(depth,splpts,nspl,splcof)
c *** otherwise: send error
      else
        istatus=999
      endif
      return
      end
