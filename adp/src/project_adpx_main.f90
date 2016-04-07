program project_adpx_main

  use project_adpx_module
  implicit none

  real               :: t1,t2,t3,t4,t5,t6
  integer            :: iplot
  character(len=120) :: check 
  character(len=3)   :: project, inversion

! -------------------------------------------------------------
! welcome message

  print*,""
  print*,""
  print*,""
  print*,""  
  print*," _|_|_|    _|_|_|      _|_|          _|  _|_|_|_|    _|_|_|  _|_|_|_|_|"
  print*," _|    _|  _|    _|  _|    _|        _|  _|        _|            _|    "  
  print*," _|_|_|    _|_|_|    _|    _|        _|  _|_|_|    _|            _|    "  
  print*," _|        _|    _|  _|    _|  _|    _|  _|        _|            _|    "  
  print*," _|        _|    _|    _|_|      _|_|    _|_|_|_|    _|_|_|      _|    "
  print*,""
  print*,""
  print*,""
  print*,"Welcome to Project version alpha 0.1 (c) L.Auer, L.Boschi, J.Schaefer 2013"

! -------------------------------------------------------------
! read parameters from stdin

  call cpu_time(t1)

  print*,""
  print*,""
  print*,"--------------------- READ PARAMETERS FROM STDIN ------------------------"
  print*,""
  print*,""

  write(*,'(A,$)'),"Adaptive grid that varies with dept? (0=yes,>0=no)"
  read*,adaptmode
  write(*,*),adaptmode

  write(*,'(A,$)'),"Enter dimension of sparse matrix arrays?"
  read*,nonz
  write(*,*),nonz

  write(*,'(A,$)'),"Enter equatorial cell size in degree:"
  read*,eq_incr
  write(*,'(3f8.3)'),eq_incr

  write(*,'(A,$)'),"Compatible with crustal model of same gridsize (0=no,1=yes): "
  read*,iswit
  write(*,*),iswit

  write(*,'(A,$)'),"Restrict variable grid to larger Europe (1) or USA (2)? (0=no)?"
  read*,outsideREG
  write(*,*),outsideREG

  ! Some switches and definitions and allocations
  refgrid   = 5.0 ! reference grid is always 5.0
  fact      = refgrid/eq_incr ! relation between reference and finest grid
  nlev      = log(refgrid/eq_incr)/log(2.)+1	! number of different pixel sizes
  nlatzomax =180/eq_incr            ! maximal number of latitudinal zones in finest grid 
  write(*,*), "Nr. of pixel sizes in variable grid: ",nlev

  ! Allocate memory
  allocate(ithres(nlev),nsqrs(nlatzomax,nlev),nsqtot(nlatzomax+1,nlev),n1layer(nlev))

  write(*,'(A,$)'),"Enter threshold values for hitcount:"
  read(*,*)(ithres(lev),lev=1,nlev-1)
  write(*,*),(ithres(lev),lev=1,nlev-1)

  write(*,'(A,$)'),"Enter minimum value of hits per cell:"
  read(*,*)ithres(nlev)
  write(*,*),ithres(1:nlev-1)

  write(*,'(A,$)'),"Enter number of layers: "
  read*,nlayi
  write(*,*),nlayi

  write(*,'(A,$)'),"Enter number of layers: "
  read*,writegrid
  write(*,*),writegrid
 
  print*,"Path to VAL array of sub-matrix: (reg/var)"
  read*,namexxx
  read*,namexxxad
  print*, trim(namexxx)
  print*, trim(namexxxad)

  print*,"Path to IND array of sub-matrix: (reg/var)"
  read*,nameind
  read*,nameindad  
  print*, trim(nameind)
  print*, trim(nameindad)

  print*,"Path to PNT vector of sub-matrix (reg/var):"
  read*,namepoi
  read*,namepoiad
  print*, trim(namepoi)
  print*, trim(namepoiad)

  print*,"Path to RHS vector of sub-matrix (reg/var): "
  read*,namerhs
  read*,namerhsad
  print*, trim(namerhs)
  print*, trim(namerhsad)

  if (writeata.eq.1) then
     print*,"Paths to ATA and ATD output file: "
     read*,nameata
     read*,nameatd  
     print*,trim(nameata)
     print*,trim(nameatd)
  end if

  print*,"Path to hitcount/dws file: "
  read*,namehits
  print*,trim(namehits)

  print*,"Files to store grid information: "
  read*,namegrid
  read*,nameadpx
  print*,trim(namegrid)
  print*,trim(nameadpx)

  print*,"Files to store htc/dws after projection: "
  read*,namehtcadgrid
  read*,namenumberadpx
  print*,trim(namehtcadgrid)
  print*,trim(namenumberadpx)
  
  write(*,'(A,$)'),"Enter dimension of right-hand-side vector?"
  read*,m
  write(*,*),m
  r1=1.
  r2=m
  
  nlay=nlayi+nlayouc
  if (nlay.gt.nlaym) then
     print*,'ERROR: number of layers out of bounds, stopping ...',nlaym
     stop
  end if 

  call cpu_time(t2)
  write(*,*), "Time to read from stdin: ",t2-t1



! -------------------------------------------------------------
! Define variable grid parameterization

  print*,""
  print*,""
  print*,"-------------- DEFINING VARIABLE GRID PARAMETERIZATION ------------------"
  print*,""
  print*,""

  info=1
  nlatzones=180/eq_incr ! not sure why this is here (ludwig)
  call param_adpx ! subroutine defines parameterization
  if (info.ne.0) then
     write(*,*)'ERROR: inversion: param_adpx failed, stopping ...'
     goto 9999
  else
     call cpu_time(t3)
     write(*,*), "Variable grid parameterization defined!"
     write(*,*), "Time to define parameterization: ",t3-t2
  endif



! -------------------------------------------------------------
! Project matrix to variable grid

  print*,""
  print*,""
  print*,"---------------- PROJECTING MATRIX TO VARIABLE GRID ---------------------"
  print*,""
  print*,""

  info=1
  call project_matrix ! subroutine projects from regular to adaptive grid
  if (info.ne.0) then 
     write (*,*)'ERROR: project_matrix failed, stopping ...'
     goto 9999
  else
     print*, "Matrix successfully projected!"
  end if
  call cpu_time(t4)
  print*, "Time to project matrix: ",t4-t3
  


! -------------------------------------------------------------
! Finished!


  goto 11 
9999 continue ! jump here upon error
  if(info.ne.0)write(*,*)'... bye!'
  stop

11 continue

  print*,""
  print*,""
  print*," ____ ____ ____ ____ "
  print*,"||D |||o |||n |||e ||"
  print*,"||__|||__|||__|||__||"
  print*,"|/__\|/__\|/__\|/__\|"
  print*,""
  print*,""
  print*,"... check your projected matrices!"

end program project_adpx_main
