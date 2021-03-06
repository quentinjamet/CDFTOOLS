PROGRAM cdfvFWov
  !!-------------------------------------------------------------------
  !!======================================================================
  !!                     ***  PROGRAM  cdfvFWov  ***
  !!=====================================================================
  !!  ** Purpose : from a section calculate net freshwater transport and its
  !!               overturning component
  !!               section is assumed to be 2 j lines (j and j+1)
  !!
  !!  ** Method  : compute salinity at v point
  !!               compute zonal mean of salinity-vpt and v velocity
  !!               compute total freshwater transport and overturning component
  !!
  !! History : 2.1  : 12/2011  : J. Deshayes  : Original code
  !!           3.0  : 12/2011  : J.M. Molines : Port to 3.0
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class transport
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jk, jt          ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc, ijarg  ! arguments on command line
  INTEGER(KIND=4)                            :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                            :: ierr                ! error status of I/O
  INTEGER(KIND=4)                            :: ncout               ! ncid of output file
  INTEGER(KIND=4)                            :: ikx=1, iky=1, ikz=1 ! dims of netcdf output file
  INTEGER(KIND=4)                            :: nboutput = 3        ! number of output variables
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout      ! output variables properties

  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: rmaskn, rmasks      ! S-mask North and South of V-point
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: rmaskv              ! mask at V point
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zsaln, zsals        ! salinity North and South of v point
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zsalv, zvitv        ! salinity at V point, velocity
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: zFWv                ! work array
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: rdumlon, rdumlat    ! dummy arrays for I/O

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: de3v                ! vertical metrics
  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: de1v                ! horizontal metrics
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE :: dnetvFW, dtotvFW    ! transport array
  REAL(KIND=8), DIMENSION (:),   ALLOCATABLE :: dovFW, dtim         ! overturning array, time
  REAL(KIND=8), DIMENSION (1)                :: dnetv, dnetFW       ! total transport accross section (mass and FW)
  REAL(KIND=8), DIMENSION (1)                :: darea, dareak       ! work space
  REAL(KIND=8), DIMENSION (1)                :: dzonalv, dzonalFW
  REAL(KIND=8)                               :: dztrp, dcellarea
  REAL(KIND=8), PARAMETER                    :: dp_rsal=35.d0       ! reference salinity for freshwater calculation

  CHARACTER(LEN=256)                         :: cf_vfil             ! input V file
  CHARACTER(LEN=256)                         :: cf_sfil             ! input S file
  CHARACTER(LEN=256)                         :: cf_zgr              ! zgr file
  CHARACTER(LEN=256)                         :: cf_hgr              ! hgr file
  CHARACTER(LEN=256)                         :: cf_mask             ! mask file
  CHARACTER(LEN=256)                         :: cf_out    = 'vFWov.nc'  ! output file
  CHARACTER(LEN=256)                         :: cv_netvFW = 'netvFW'    ! output variable 1
  CHARACTER(LEN=256)                         :: cv_totvFW = 'totvFW'    ! output variable 2
  CHARACTER(LEN=256)                         :: cv_ovFW   = 'ovFW'      ! output variable 3
  CHARACTER(LEN=256)                         :: cglobal             ! Global attribute for output file
  CHARACTER(LEN=256)                         :: cldum               ! working char variable

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: stypvar             ! I/O data structure

  LOGICAL                                    :: lchk = .FALSE.      ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line 
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvFWov -v V-secfile -s S-secfile -zgr ZGR-secfile -hgr HGR-secfile '
     PRINT *,'                  ... -msk MSK-secfile [-o OUT-file] [-vvl] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the fresh water transport and its overturning component through'
     PRINT *,'        a section specified by the input files (data and metrics).'
     PRINT *,'      '
     PRINT *,'        All input files corresponds to a zonal extraction holding 2 rows of'
     PRINT *,'        data. In fact the 2 rows are needed for the salinity to be interpolated'
     PRINT *,'        on the V-points.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        All arguments are ''section files'', which are assumed to be files with'
     PRINT *,'        2 zonal lines of data ( j and j+1 ): '
     PRINT *,'       -v V_secfile   : meridional velocity section file.'
     PRINT *,'       -s S_secfile   : salinity section file.'
     PRINT *,'       -zgr ZGR_secfile : mesh_zgr section file '
     PRINT *,'       -hgr HGR_secfile : mesh_hgr section file '
     PRINT *,'       -msk MSK_secfile : mask section file '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-vvl] : use time-varying vertical metrics.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless option -o is used.'
     PRINT *,'       variables : ',TRIM(cv_netvFW),', ',TRIM(cv_totvFW),', ',TRIM(cv_ovFW)
     PRINT *,'       Output file only has time relevant dimension. Other dims are set to 1.'
     PRINT *,'       Degenerated dimensions can be removed with :'
     PRINT *,'           ncwga -a x,y,depthw ',TRIM(cf_out), ' -o out.nc'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg) 
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-v'   ) ; CALL getarg( ijarg, cf_vfil ) ; ijarg=ijarg+1
     CASE ( '-s'   ) ; CALL getarg( ijarg, cf_sfil ) ; ijarg=ijarg+1
     CASE ( '-zgr' ) ; CALL getarg( ijarg, cf_zgr  ) ; ijarg=ijarg+1
     CASE ( '-hgr' ) ; CALL getarg( ijarg, cf_hgr  ) ; ijarg=ijarg+1
     CASE ( '-msk' ) ; CALL getarg( ijarg, cf_mask ) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'   ) ; CALL getarg( ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-vvl' ) ; lg_vvl = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = lchk .OR. chkfile ( cf_vfil )
  lchk = lchk .OR. chkfile ( cf_sfil )
  lchk = lchk .OR. chkfile ( cf_zgr  )
  lchk = lchk .OR. chkfile ( cf_hgr  )
  lchk = lchk .OR. chkfile ( cf_mask )

  IF ( lchk ) STOP 99 ! missing files
  IF ( lg_vvl ) THEN
     cn_fe3v = cf_vfil
     cn_ve3v = cn_ve3vvvl
  ENDIF

  !! get dimensions of input file containing data
  npiglo = getdim(cf_vfil, cn_x)
  npjglo = getdim(cf_vfil, cn_y)
  npk    = getdim(cf_vfil, cn_z)
  npt    = getdim(cf_vfil, cn_t)  

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  IF ( npjglo /= 2 ) THEN
     PRINT *,' ERROR : This program works with section files.'
     PRINT *,'       all data should have j dimension equal to 2 '
     STOP 99
  ENDIF

  ALLOCATE ( dnetvFW(npt), dtotvFW(npt), dovFW(npt), dtim(npt) )
  ALLOCATE ( de1v(npiglo,npjglo), de3v(npiglo,npk))
  ALLOCATE ( zFWv(npiglo,npk), zvitv(npiglo,npk) )
  ALLOCATE ( zsals(npiglo,npk), zsaln(npiglo,npk), zsalv(npiglo,npk) )
  ALLOCATE ( rmasks(npiglo,npk), rmaskn(npiglo,npk), rmaskv(npiglo,npk))
  ALLOCATE ( stypvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE ( rdumlon(1,1), rdumlat(1,1) )

  rdumlon(:,:) = 0.e0  ! dummy longitude
  rdumlat(:,:) = 0.e0  ! dummy latitude

  CALL CreateOutput

  !! load scale factors  as vertical slice for 3D fields
  rmasks(:,:) = getvarxz(cf_mask, cn_tmask, kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  rmaskn(:,:) = getvarxz(cf_mask, cn_tmask, kj=2,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  rmaskv(:,:) = getvarxz(cf_mask, cn_vmask, kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)
  de1v(:,:)   = getvar  (cf_hgr,  cn_ve1v,  klev=1, kpi=npiglo, kpj=npjglo, kimin=1, kjmin=1)
  de3v(:,:)   = getvarxz(cn_fe3v, cn_ve3v,  kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1)

  WHERE ( rmasks /= 0. ) rmasks = 1.
  WHERE ( rmaskn /= 0. ) rmaskn = 1.
  WHERE ( rmaskv /= 0. ) rmaskv = 1.

  !! do calculation for each time step
  DO jt=1,npt
     IF ( lg_vvl ) THEN 
        ! read e3v at each jt
        de3v(:,:) = getvarxz(cn_fe3v, cn_ve3v,  kj=1,   kpi=npiglo, kpz=npk,    kimin=1, kkmin=1, ktime=jt )
     ENDIF
     PRINT *,'jt =',jt
     ! reset cumulative arrays to 0
     zFWv(:,:)   = 0.e0
     dnetv       = 0.d0  ;  dnetFW = 0.d0 
     darea       = 0.d0
     dnetvFW(jt) = 0.d0  ;  dtotvFW(jt) = 0.d0  ; dovFW(jt) = 0.d0

     zvitv(:,:)= getvarxz(cf_vfil, cn_vomecrty, kj=1, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt) 
     zsals(:,:)= getvarxz(cf_sfil, cn_vosaline, kj=1, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt)
     zsaln(:,:)= getvarxz(cf_sfil, cn_vosaline, kj=2, kpi=npiglo, kpz=npk, kimin=1, kkmin=1, ktime=jt)

     DO jk = 1, npk
        DO ji = 1, npiglo
           IF ( rmasks(ji,jk) + rmaskn(ji,jk) /=  0 )  THEN
              zFWv(ji,jk) = ( dp_rsal - ( zsals(ji,jk) * rmasks(ji,jk) + zsaln(ji,jk) * rmaskn(ji,jk) ) &
                   &           / ( rmasks(ji,jk) + rmaskn(ji,jk) ) ) / dp_rsal ! freshwater at Vpoint
           ENDIF
           dcellarea   = de1v(ji,1) * de3v(ji,jk) * rmaskv(ji,jk)
           dztrp       = zvitv(ji,jk) * dcellarea
           dnetvFW(jt) = dnetvFW(jt) + zFWv(ji,jk) * dztrp /1.d6         ! net freshwater transport in Sv
           dnetv       = dnetv       +               dztrp
           dnetFW      = dnetFW      + zFWv(ji,jk) * dcellarea
           darea       = darea       + dcellarea
        END DO !ji
     END DO !jk

     PRINT *,'total mass transport across section =', dnetv / 1.d6,' Sv'

     dnetv  = dnetv  / darea   ! mean velocity across section
     dnetFW = dnetFW / darea   ! mean freshwater along section
     PRINT *,'mean salinity along section =', dp_rsal - dnetFW * dp_rsal,' psu'

     DO jk = 1, npk
        dzonalv  = 0.d0
        dzonalFW = 0.d0
        dareak   = 0.d0 
        DO ji = 1, npiglo
           dcellarea = de1v(ji,1) * de3v(ji,jk) * rmaskv(ji,jk)
           dzonalv   = dzonalv  + ( zvitv(ji,jk) - dnetv  ) * dcellarea
           dzonalFW  = dzonalFW + ( zFWv(ji,jk)  - dnetFW ) * dcellarea
           dareak    = dareak   +                             dcellarea

           dtotvFW(jt) = dtotvFW(jt) + ( zvitv(ji,jk) - dnetv(1) ) * zFWv(ji,jk) * dcellarea /1.d6 
        END DO !ji

        IF ( dareak(1) > 0 ) THEN
           ! overturning freshwater transport in Sv
           dovFW(jt) = dovFW(jt) + dzonalv(1) * dzonalFW(1) / dareak(1) /1.d6  
        ENDIF
     END DO !jk

     PRINT *,'netvFW = ', dnetvFW(jt), ' Sv'
     PRINT *,'totvFW = ', dtotvFW(jt), ' Sv'
     PRINT *,'ovFW   = ', dovFW(jt),   ' Sv'

     ierr = putvar0d( ncout, id_varout(1), REAL(dnetvFW(jt)), ktime = jt )
     ierr = putvar0d( ncout, id_varout(2), REAL(dtotvFW(jt)), ktime = jt )
     ierr = putvar0d( ncout, id_varout(3), REAL(dovFW(jt))  , ktime = jt )

  END DO  !jt

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    !! define output variables
    ipk(:) = 1
    stypvar%rmissing_value    = 0.
    stypvar%valid_min         = -1000.
    stypvar%valid_max         = 1000.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'T'
    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.

    stypvar(1)%cname          = TRIM(cv_netvFW )
    stypvar(2)%cname          = TRIM(cv_totvFW )
    stypvar(3)%cname          = TRIM(cv_ovFW   )

    stypvar(1)%cunits         ='Sv'
    stypvar(2)%cunits         ='Sv'
    stypvar(3)%cunits         ='Sv'

    stypvar(1)%clong_name     = 'Net transport of freshwater across section'
    stypvar(2)%clong_name     = 'Transport of freshwater across section when net mass transport equals 0'
    stypvar(3)%clong_name     = 'Overturning component of freshwater transport across section'

    stypvar(1)%cshort_name    = TRIM(cv_netvFW )
    stypvar(2)%cshort_name    = TRIM(cv_totvFW )
    stypvar(3)%cshort_name    = TRIM(cv_ovFW   )

    WRITE(cglobal,'(a,f5.1,a)' ) 'comment : Reference salinity ', dp_rsal,' transport is positive northward'

    !! prepare output file
    ncout = create      (cf_out, 'none',  ikx,      iky, ikz,       cdep='depthw'                    )
    ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout, cdglobal=TRIM(cglobal)           )
    ierr  = putheadervar(ncout,  cf_vfil, ikx,      iky, ikz,       pnavlon=rdumlon, pnavlat=rdumlat )

    dtim  = getvar1d(cf_vfil, cn_vtimec, npt)
    ierr  = putvar1d(ncout,   dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfvFWov
