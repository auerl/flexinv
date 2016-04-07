c----------------------------------------------------------------------
c  This subroutine gets the crustal type (2 characters) corresponding
c  to a given latitude and longitude, based on the model CRUST 2.0
c  by Laske et al. (GE 2001-06-20)
c----------------------------------------------------------------------
c
      subroutine gettype_20 (lu,xlat,xlon,type)
      character*1000 strtwo
      character*80 dbsdir
      character*128 infla
      character*2 type
      character*2 abbv(90,180)
      data irdtp/0/
      save
      
      if (irdtp.eq.0) then
        call getenv('HRVDBS',dbsdir)
        ldbs=lnblnk(dbsdir)
        infla=dbsdir(1:ldbs)//'/CRUST_2.0/CNtype2.txt'
        open (lu,file=infla)
        ios=0
        iline=1
        iuse=0
        read(lu,"(a)") strtwo
        lstr=lnblnk(strtwo)
c        write(6,"('length of line:',i6)") lstr
        do ila=1,90
          read(lu,"(a)") strtwo
          iline=iline+1
        lstr=lnblnk(strtwo)
c        write(6,"('length of line:',2i6)") lstr,iline
          read(strtwo,"(i4,3x,180(a2,3x))") ii,(abbv(ila,i),i=1,180)
          if(ii.ne.(90-(ila-1)*2)) then
            write(6,"('incorrect latitude:',2i6)") ii,90-(ila-1)*2
            stop
          endif
        enddo
        irdtp=1
        close(lu)
      endif
      
      lat=int(1+((90.-xlat)/2.))
      if (lat.eq.91) lat=90
      lon=int(1+((180.+xlon)/2.))
      if (lon.eq.181) lon=1
      type=abbv(lat,lon)
            
      if (lat.gt.91.or.lat.lt.1.or.xlat.lt.-90.) then
        write (6,"('invalid latitude ',f8.2,i4)") xlat,lat
        type='--'
        stop
      endif
      if (lon.lt.1.or.lon.gt.181) then
        write (6,"('invalid longitude ',f8.2,i4)") xlon,lon
        write (6,"('longitude must be between -180 and 180 deg.E')")
        type='--'
        stop
      endif
c      write(6,"('returning from gettype_20:',a2)") type
      return
      end
c
c----------------------------------------------------------------------
c  This subroutine gets the topography and bathymetry at a given 
c  latitude and longitude, based on the 2x2 degree model CRUST 2.0
c  by Laske et al. (GE 2001-06-20)
c----------------------------------------------------------------------
c
      subroutine getbael_20 (lu,xlat,xlon,bael)
      character*2000 strtwo
      character*80 dbsdir
      character*128 infla
      dimension baelarr(90,180)
      data irdtp/0/
      save
c
      if (irdtp.eq.0) then
        call getenv('HRVDBS',dbsdir)
        ldbs=lnblnk(dbsdir)
        infla=dbsdir(1:ldbs)//'/CRUST_2.0/CNelevatio2.txt'
        open (lu,file=infla)
        ios=0
        read(lu,"(a)") strtwo
        lstr=lnblnk(strtwo)
c        write(6,"('length of line:',i6)") lstr
        do ila=1,90
          read(lu,"(a)") strtwo
          lstr=lnblnk(strtwo)
c          write(6,"('length of line:',i6)") lstr
          read(strtwo,*) ii,(baelarr(ila,i),i=1,180)
          if(ii.ne.(90-(ila-1)*2)) then
            write(6,"('incorrect latitude:',2i6)") ii,90-(ila-1)*2
            stop
          endif
        enddo
        irdtp=1
        close(lu)
      endif
      
      lat=int(1+((90.-xlat)/2.))
      if (lat.eq.91) lat=90
      lon=int(1+((180.+xlon)/2.))
      if (lon.eq.181) lon=1
      bael=baelarr(lat,lon)
c      write(6,"('returning from getbael:',f10.2)") bael
      return
      end
c
c----------------------------------------------------------------------
c  This subroutine gets the structure (thicknesses and velocities)
c  for a given crustal type, based on the model CRUST 2.0
c  by Laske et al. (GE 2001-06-20)
c----------------------------------------------------------------------
c
      subroutine cstruc_20(lu,type,desc,vptyp,vstyp,rhtyp,thtp,tpmo,inflb,ierr)
      character*2 type
      character*80 dbsdir
      character*100 junk
      character*128 inflb
      character*4 trash
      character*60 descrip(500)
      character*60 desc
      character*80 crustfoldname
      character*2 code(500)
      dimension veloc1(500,8)
      dimension veloc2(500,8)
      dimension veloc3(500,8)
      dimension thlr(500,7)
      dimension dmoho(500)
      dimension vptyp(8)
      dimension vstyp(8)
      dimension rhtyp(8)
      dimension thtp(7)
      data irdst/0/
      save
      
      if (irdst.eq.0) then
        call getenv('HRVDBS',dbsdir)
c        ldbs=lnblnk(crustfoldname)
c        inflb=dbsdir(1:ldbs)//'/CRUST_2.0/CNtype2_key.txt'
c        inflb='/work/dev/comb_mat_inv/MAT/SUWA/KERNELS/CRUST20/CNtype2_key.txt'  ! Lapo 10 Feb 04
c        inflb=crustfoldname(1:ldbs)//"/CNtype2_key.txt"

        write(6,"(a)"),inflb
        open(lu,file=inflb)

cTEST
c        print*,inflb
c        pause

        ios=0
        iline=0
        iuse=0
        read(lu,"(a)") junk
        read(lu,"(a)") junk
        read(lu,"(a)") junk
        read(lu,"(a)") junk
        read(lu,"(a)") junk
        do while (ios.eq.0)
          iline=iline+1
          read (lu,"(a2,a)",iostat=ios) code(iline),descrip(iline)
          read (lu,*,iostat=ios) (veloc1(iline,i),i=1,8)
          read (lu,*,iostat=ios) (veloc2(iline,i),i=1,8)
          read (lu,*,iostat=ios) (veloc3(iline,i),i=1,8)
          read (lu,*,iostat=ios) (thlr(iline,j),j=1,7),trash,
     #       dmoho(iline)
        enddo
        close(lu)
        write (6,"('read ',i4,' lines - cstruc5')") iline
        irdst=1
      endif
      
      m=1
      do k=1,iline
        if (code(k).eq.type) then
          desc=descrip(k)          
          do i=1,8
            vptyp(i)=veloc1(k,i)
            vstyp(i)=veloc2(k,i)
            rhtyp(i)=veloc3(k,i)
          enddo
          do i=1,7
            thtp(i)=thlr(k,i)
          enddo
          tpmo=dmoho(k)
          m=0
        endif
        if (m.ne.0) m=1
      enddo

      ierr=0 ! hook, la
      return
      end
