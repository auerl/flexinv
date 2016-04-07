!     -------------------------------------------------------------
! this contains global parameters for flexinv90

    module flexinv_module
        
        implicit none

        integer*8, dimension(:), allocatable :: pnt0
        integer*4, dimension(:), allocatable :: ind0
        real, dimension(:), allocatable :: rhs0
        real, dimension(:), allocatable :: rhs0b
        real*4, dimension(:), allocatable :: val0

        real, dimension(:), allocatable :: rhs
        integer*8, dimension(:), allocatable :: ind,indb
        integer*8, dimension(:), allocatable :: pnt
        real, dimension(:), allocatable :: val,valb
        real*8, dimension(:), allocatable :: ata
        real*8, dimension(:), allocatable :: atd

        integer*8 n,m,nonz,n0,n1,nlay,irec

	integer, parameter :: nlaym=50 ! maximum number of layers
	integer, parameter :: n0max=50000 ! maximum number of cells in each layer
	integer, parameter :: n1max=n0max*nlaym ! maximum number of voxels per parameter
	integer, parameter :: nparmax=4 ! number of inversion parameters, 2 = anisotropic
	integer, parameter :: nlatzomax=180 ! maximum number of latitudinal zones
        integer, parameter :: ntemp=40000**2
        integer*8 :: npar ! Nr of physical inversion parameters

        ! parameterization
        integer, dimension(nlatzomax)   :: nsqrs
        integer, dimension(nlatzomax+1) :: nsqtot


        ! adaptive parameterization
        integer*8, dimension(:) :: adpx1layer(nlaym)
        integer*8 :: adpxtotal,natamax
        ! variance reduction

	integer, parameter :: nmax=nparmax*n1max ! maximum total number of parameters       

        integer*8 nelm,nelm0,nrhs,nrhs0
        integer*4 iswit,iata,iadapt,nlatzones


!	real, dimension(nonz) :: val

!     -------------------------------------------------------------
!       switches and variables associated with damping

        integer idamp ! damp roughness/model gradient?
        integer inormdmp ! additional mid mantle norm damping?
        integer iadaptnormdmp ! adaptive grid norm damping?



        ! weighting factors
	real, dimension(nparmax) :: wgrad ! horizontal rdamp weight       
	real, dimension(nlaym) :: wgradv ! vertical rdamp factor
	real, dimension(nlaym,nparmax) :: wgradh ! vertical rdamp factor

        real*8, dimension(4) :: area ! area of cells for adaptive norm damping
        real*4, dimension(nmax) :: ndampvec1
        real*4, dimension(nlaym,nparmax) :: wddampvec
        real, dimension(nparmax) :: ndampvec0 ! norm damp weight factor


        real*4 wrdamp,wndamp,wddamp,dampsize
        real*4 af1,af2,af3,af4,ad1,ad2,afv,preani
        real*4 eqincr,refgrid,tinitial,adtime
        real*4 cutoff,cutoffn,relwei,wei2

!     -------------------------------------------------------------
!       some directory and file names

	character*140 namexxx,namepoi,nameind,namerhs ! submatrices
        character*140 outfile,wdir,nomfil,gridinfo,dmpsched_r,dmpsched_d
        character*140 preani_model

    end module flexinv_module
