!------------------------------------------------------------------
!Compile with:
!/scratch/cantal/opt/hdf5-1.8.17/hdf5/bin/h5fc -I /home/astraw/codes/lib -L /home/astraw/codes/lib -lStARTLib -lcfitsio -ffree-line-length-none -static-libgfortran -c AnalysisTool.f90

!/scratch/cantal/opt/hdf5-1.8.17/hdf5/bin/h5fc -I /home/astraw/codes/lib -L /home/astraw/codes/lib -lStARTLib -lcfitsio -ffree-line-length-none -static-libgfortran AnalysisTool.o -o /home/astraw/bin/Analyze
!
!-------------------------------------------------------------------

PROGRAM main

  USE StARTLib
  IMPLICIT NONE
  CHARACTER(len=250) :: inpfile, outfile, bovfile
  CHARACTER(len=40)  :: var, ftype
  INTEGER      :: action, lmax, proj, TStep, imsize(3), ii, i,j,k, im_pos(3), x, y, z, lc(3), rc(3), lmin, xx, yy, Coords(3), nbins, bin, orig_lev, &
       ix, iy, iz, PSF_r, smooth_r
  REAL(kind=4)        :: xc(3),rad, dist, rmin, rmax, rmin2, rmax2, den, temp, dmin, dmax, tmin, tmax, dd, dt, inv_dd, inv_dt, maxv(3), box_size(3), &
       box_orig(3), sigma, stars2SB, ThisGas, minden, xs(3), OpAngle, dist2, pos_(3), oadir(3), mdist
  REAL(kind=8) :: box_le(3), box_re(3), pos(3), ref, cell_le(3), cell_re(3), inv_imsize(3), boxfact(3), avg, em_fact, NtrFracSource
  !REAL(kind=8) :: pos(3)
  REAL(kind=4), ALLOCATABLE :: im(:,:,:), weight(:),image(:,:), weight_SB(:), avrg(:)
  REAL(kind=4)    :: ClumpingFact, redshift_, SSThr, Gamma12, SFThr, nHI, nHII, ThisSB, ne, SaturationDensity, SFT
  REAL(kind=4)    :: GaussianKernel(7)=[0.006,0.061,0.242,0.383,0.242,0.061,0.006], ForceRedshift, ThisStars
  REAL(kind=4), PARAMETER :: fourpi_in_arcsec2=5.3407075e11
  REAL(kind=4), PARAMETER :: arcsec2_to_4pi=1./fourpi_in_arcsec2
  REAL(kind=4), PARAMETER :: LyaPh_to_erg=10.2*1.602176e-12
  REAL(kind=4), ALLOCATABLE :: nc(:)
  LOGICAL                 :: ReadExtra, IncludeCollEx, IncludeRec
  INTEGER ::  status,unit,blocksize,bitpix,naxis,fpixel,nelements,group, rwstatus
  LOGICAL :: simple, extend, ex
  INTEGER, ALLOCATABLE :: naxes(:)
  REAL(kind=8)          :: R(2), v1, v2, rsq

  CALL read_param
  print *, 'reading ', TRIM(inpfile)
  CALL ReadStARTGrids(inpfile, TStep, ReadExtra)

!  DO ii=1,NGrids
!    DO k=1,Grid(ii)%Dims(3)
!      DO j=1,Grid(ii)%Dims(2)
!        DO i=1,Grid(ii)%Dims(1)
!          Grid(ii)%Cell(i,j,k)%Density=MIN(Grid(ii)%Cell(i,j,k)%Density,SFThr)
!        END DO
!      END DO
!    END DO
!  END DO

CONTAINS

SUBROUTINE read_param

    IMPLICIT NONE
    INTEGER :: n,i,iargc
    CHARACTER(len=250) :: opt, arg
    LOGICAL :: bad, ok

    !..set defaults
    action=1
    box_le=0.
    box_re=1.
    lmax=-1 ! it will be fixed later on
    proj=1
    var="f_HI"
    lmin=0   
    xc=0.5
    rmax=1.
    rmin=0.
    SaturationDensity=1.e19
    nbins=100
    IncludeCollEx=.true.
    IncludeRec=.true.
    ForceRedshift=-1
    SFT=1.e19
    ftype="bov"
    sigma=0.
    PSF_r=0
    smooth_r=0
    stars2SB=1
    minden=0
    !nvar=1
    !nfiles=1
    NtrFracSource=0.9d0
    xs=-1
    OpAngle=90.0
    oadir=[0,0,1]
    dmax=1.e9

    n=iargc()
    IF(n<4) THEN
       print *, 'usage: PlotStART  -inp StART_file -out filename'
       print *, 'options: [-act] 1=projected map (default); 2=T-rho plot; 3=spherical avg'
       print *, '         [-box_le] cutting box [default=0,0,0]'
       print *, '         [-box_re] cutting box [default=1,1,1]'
       print *, '         [-lmax]   max level'
       print *, '         [-lmin]   min_level'
       print *, '         [-proj]   projection dim [default=1]'
       print *, '         [-var]    variable [default=f_HI]'
       print *, '         [-cen]    region centre (box units)'
       print *, '         [-rma]    max distance from centre'
       print *, '         [-rmi]    min distance from centre'
       print *, '         [-sdn]    SaturationDensity [def=1.e19]'
       print *, '         [-sft]    star formation (phys.) density thr [default=1.e19]'
       print *, '         [-minden] minimum (phys.) density for estimating Mgas [default=0.]'
       print *, '         [-nbins]  number of bins for spherical avg [default=100]'
       print *, '         [-cex]    include collisional excitations  [default=.true.]'
       print *, '         [-rec]    include recombinations           [default=.true.]'
       print *, '         [-fr]     force redshift to this value [default from input file]'
       print *, '         [-ftype]  output file type (bov/fits) [default=bov]'
       print *, '         [-addnoise]  add noise with this sigma (per res. elem.) [def=0]'
       print *, '         [-psfrad] add gaussian psf (rad in pixel) [default=0]'
       print *, '         [-smrad]  smooth (gauss) over smrad pixels [default=0]'
       print *, '         [-stars2SB]  stellar mass to ly-alpha SB conv fact [def=1]'
       print *, '         [-fHIs]  fHI for a source cell [def=0.9d0]'
       print *, '         [-xs]    source position for opening angle [def= no source]'
       print *, '         [-oa]    source half-cone opening angle in dir "oadir" in degrees [def= no source]'
       print *, '         [-oadir] source opening angle versor dir as defined in Radamesh-1.3 [def=z-dir]'
       print *, '         [-mdist]  max distance from source in box units [def=inf]'
       STOP
    END IF



    DO i=1,n,2
       CALL GETARG(i,opt)
       CALL GETARG(i+1,arg)
       print *, TRIM(opt), " ", TRIM(arg)

       SELECT CASE(TRIM(opt))
         CASE('-inp')
            READ(arg,'(a)') inpfile
         CASE('-out')
            READ(arg,'(a)') outfile
         CASE('-act')
            READ(arg,*) action
         CASE('-box_le')
            READ(arg,*) box_le
         CASE('-box_re')
            READ(arg,*) box_re
         CASE('-lmax')
            READ(arg,*) lmax
         CASE('-proj')
            READ(arg,*) proj
         CASE('-var')
            READ(arg,*) var
         CASE('-lmin')
            READ(arg,*) lmin
         CASE('-cen')
            READ(arg,*) xc(:)
         CASE('-rma')
            READ(arg,*) rmax
         CASE('-rmi')
            READ(arg,*) rmin
         CASE('-sdn')
            READ(arg,*) SaturationDensity
         CASE('-nbins')
            READ(arg,*) nbins
         CASE('-cex')
            READ(arg,*) IncludeCollEx
         CASE('-rec')
            READ(arg,*) IncludeRec
         CASE('-fr')
            READ(arg,*) ForceRedshift
         CASE('-sft')
            READ(arg,*) SFT
         CASE('-ftype')
            READ(arg,*) ftype
         CASE('-addnoise')
            READ(arg,*) sigma
         CASE('-psfrad')
            READ(arg,*) PSF_r
         CASE('-smrad')
            READ(arg,*) smooth_r
         CASE('-stars2SB')
            READ(arg,*) stars2SB
         CASE('-minden')
            READ(arg,*) minden
         CASE('-fHIs')
            READ(arg,*) NtrFracSource
         CASE('-xs')
            READ(arg,*) xs
         CASE('-oa')
            READ(arg,*) OpAngle
         CASE('-oadir')
            READ(arg,*) oadir
         CASE('-mdist')
            READ(arg,*) mdist
         CASE default
            print '("unknown option ",a2," ignored")', TRIM(opt)
       END SELECT
    END DO

    OpAngle=cos(OpAngle)
     
    IF(action==2) THEN
       proj=3
       var='plot'
    END IF

  END SUBROUTINE read_param

END PROGRAM main
