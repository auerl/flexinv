module project_adpx_module

  implicit none

  integer               :: info                 ! to check if routines work properly
  character(len=1)   	:: packed

  ! some constants
  real, parameter       :: cmbradius=3480.      !radius CMB
  double precision,parameter :: pi=3.1415926536

  ! anisotropy
  integer,parameter     :: npar=4               !number of anisotropy parameters, npar=1: invert for isotropic Vs

  !  reference grid
  real                  :: refgrid              !pixelsize in reference grid
  integer               :: iswit                ! to choose if grid is compatible with reference grid (0=no,1=yes)
  integer               :: writegrid            ! only write grid information in first proj job

  !  defining the grid, number of parameters
  integer,parameter     :: writeata=0
  integer,parameter	:: layernumber=2	! parameterisation is defined according to hitcount in this layer
 						! the first one should not be chosen since layer 
  						! is in some areas thinner than the crust and than 
 						! there is per definition no hitcount since we do not
  						! invert for the crust
  real                  :: eq_incr              ! size of smallest pixel 
  integer               :: nlatzones            ! number of latitudinal zones
  integer               :: nlatzomax            ! maximal number of latitudinal zones in finest grid
  integer,dimension(:,:),allocatable	:: nsqrs    ! number of pixel in each latitudinal zone dimension is (nlatzomax,nlev) 
  integer,dimension(:,:),allocatable	:: nsqtot   ! total number of pixel before actual colatitudinal zone dimension is (nlatzomax+1,nlev)
  integer,dimension(:)  ,allocatable	:: n1layer  ! total number of pixel in 1 layer at different pixel size levels
  integer,parameter     :: n0max=130000 !32000      ! maximum number of pixel per layer
  integer               :: nlayi                    ! number of layers in mantle, read in

  ! outer core 
  integer,parameter     :: iflaouc=0            ! iflaouc=1/0 include/not outer core structure in the inversion
  integer,parameter     :: nlayouc0=10          ! number of layers in core
  integer,parameter     :: nlayouc=nlayouc0*iflaouc         ! actual number of layers in outer core (takes into account if outer core structure is included in inversion)
  integer,parameter     :: nparoucmax=iflaouc*n0max*nlayouc ! maximum number of parameters to invert for in outer core ????

  ! CMB structure 
  integer,parameter     :: iflacmb=0            ! iflaouc=1/0 include/not CMB structure in the inversion
  integer               :: ngridcmb             ! actual number of parameters to invert for at CMB
  integer,parameter     :: ngridcmbmax=n0max*iflacmb ! maximum number of parameters to invert for at CMB

  ! total 
  integer               :: nlay                 ! total number of layers to invert for
  integer,parameter     :: nlaym=30             ! maximum possible number of layers 
  integer,parameter     :: n1max=n0max*nlaym    ! maximum number of pixel  
  integer,parameter     :: nmax=npar*n1max+ngridcmbmax+nparoucmax ! maximum number of parameters 

  ! for adaptive parameterization 
  integer               :: outsideREG           ! if outside a certain region (eg Europe) biggest pixel should be assumemd set to 1 
  integer               :: fact                 ! relation between reference and finest grid
  integer               :: nlev,lev             ! max. number of different pixel sizes, level index for ithres

  integer,dimension(:),allocatable  :: ithres ! threshold values for hitcount at different levels of parameterization. Dimension is nlev-1
  integer,dimension(:,:),allocatable:: htcad  ! array with hitcount values for the adaptive grid
  integer                           :: n1layer_adpx_total      ! total number of pixel in 1 adaptive parameterized layer 
  integer,dimension(:)              :: n1layer_adpx(nlaym)
  integer                           :: npx_adpx                ! total number of pixel for adaptive parameterized model
  integer                           :: npar_adpx,npar_reg      ! total number of parameter for adaptive parameterized model: npx_adpx*npar

  integer,parameter                 :: nadmax=300000                 ! maximum number of parameters to invert for
  integer                           :: adaptmode

  ! matrix in adaptive grid
  integer,dimension(nmax)           :: inew,inew0      ! new index in adaptive grid for all of the smallest pixel
  real,dimension(:),allocatable     :: ata,atd
  real,dimension(:,:),allocatable   :: ata2


  ! matrices, vectors, hitcounts, ...
  integer                           :: m,nonz   ! maximum dimension of right-hand-side vector and sparse matrix array
  integer,dimension(6)              :: ipo
  real,parameter                    :: relwei=1.0   ! weight for the data
  integer                           :: ndata    ! number of observations in this subset
  integer(8)                        :: jj       ! jj: index matrix-element
  integer,dimension(:,:),allocatable:: htc      ! hitcount
  integer                           :: r1,r2
  integer                           :: i,j

  ! characters
  character(len=200)           :: namexxx,namepoi,nameind,namerdamp,namexxxad,nameindad,namepoiad,namerhsad
  character(len=200)           :: namehits,namegrid,nameadpx,namexsadgrid,namenumberadpx,namehtcadgrid 
  character(len=200)           :: nameata,nameatd
  character(len=200)           :: namerhs,name


end module project_adpx_module
