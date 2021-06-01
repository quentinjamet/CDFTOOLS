PROGRAM cdf_dynadv_ubs_eddy_mean 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynadv_ubs_eddy_mean  ***
  !!=====================================================================
  !!  ** Purpose : Compute eddy/mean momentum advection trends following UBS advection scheme.
  !!
  !!  ** Method  : Adapt NEMO dynadv_ubs.F90 to CDFTOOLS (cf below for further details)
  !!               following the parameter used in the configuration eNATL60.
  !!
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet & J.M. Molines : Original code
  !!                : 05/2021  : Q. Jamet & J.M. Molines : Turn the computation layer per layer
  !!                                                       to avoid memory issues.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                   :: wp=4
  INTEGER(KIND=4), PARAMETER                   :: pnvarout = 4             ! number of output variables
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index for computation of grad/div
  INTEGER(KIND=4)                              :: jkkm1=1, jkk=2, jkkp1=3  ! for swapping the levels
  INTEGER(KIND=4)                              :: jpkk=3                   ! number of level to load (jkkm1, jkk, jkkp1)
  INTEGER(KIND=4)                              :: ncout_u, ncout_v         ! ncid of output file
  INTEGER(KIND=4)                              :: ncout_ke1, ncout_ke2     ! ncid of output file
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke1            ! id of output variables (ke-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke2            ! id of output variables (ke-comp)

  REAL(wp), PARAMETER                          :: gamma1 = 1._wp/3._wp     ! =1/4 quick      ; =1/3  3rd order UBS
  REAL(wp), PARAMETER                          :: gamma2 = 1._wp/32._wp    ! =0   2nd order  ; =1/32 4th order centred
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t, u,v)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshnm                    ! now sea surface height - mean
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: tmask, fmask             ! mesh masks
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: umask, vmask             ! mesh masks
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wn                       ! 3D vert. velocity (now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wnm                      ! 3D vert. velocity (now) - mean
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE     :: un, vn               ! 3D hz. velocity (now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE     :: tmpu, tmpv               ! tmeporary
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: unm, vnm               ! 3D hz. velocity (now) - mean
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_u, adv_z_u         ! hor. and vert. advection of u-mom. (outputs)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_v, adv_z_v         ! hor. and vert. advection of v-mom. (outputs)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_ke, adv_z_ke       ! hor. and vert. advection of KE     (outputs)

  CHARACTER(LEN=256)                           :: cf_tt                    ! temperature netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_uu                    ! zonal vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_um                    ! MEAN zonal vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_vv                    ! merd vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_vm                    ! MEAN merd vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_ww                    ! vert. vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_wm                    ! MEAN vert. vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_ssh                   ! Sea surface height
  CHARACTER(LEN=256)                           :: cf_sshm                  ! MEAN sea surface height
  CHARACTER(LEN=255)                           :: cf_mh                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask       netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u ='adv_u.nc'     ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v ='adv_v.nc'     ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke1='adv_ke1.nc'  ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke2='adv_ke2.nc'  ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar1                 ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes (ke1-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar4                 ! structure for attibutes (ke2-comp)

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynadv_ubs -t T-file -u U-file -um Um-file -v V-file -vm Vm-file ...'
     PRINT *,'          -w W-file -wm Wm-file -ssh SSH-file -sshm SSHm-file...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-u -o_v OUT-file-v -o_ke1 OUT-file-ke1 -o_ke2 OUT-file-ke2'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the eddy/mean decomposition of momentum/KE advection trend following UBS advection scheme.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for temperature (for mesh only)'
     PRINT *,'       -u U-file          : netcdf file for      zonal velocity'
     PRINT *,'       -um Um-file        : netcdf file for MEAN zonal velocity'
     PRINT *,'       -v V-file          : netcdf file for      meridional velocity'
     PRINT *,'       -vm Vm-file        : netcdf file for MEAN meridional velocity'
     PRINT *,'       -w W-file          : netcdf file for      vertical velocity'
     PRINT *,'       -wm Wm-file        : netcdf file for MEAN vertical velocity'
     PRINT *,'       -ssh SSH-file      : netcdf file for      SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -sshm SSHm-file    : netcdf file for MEAN SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       -o_u OUT-file      : netcdf file for advection term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for advection term for v-momentum'
     PRINT *,'       -o_ke1 OUT-file    : netcdf file for advection term for KE1'
     PRINT *,'       -o_ke2 OUT-file    : netcdf file for advection term for KE2'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      '
     STOP
  ENDIF

  cglobal = 'Partial step computation'

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tt ) ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_uu ) ; ijarg=ijarg+1
     CASE ('-um'       ) ; CALL getarg( ijarg, cf_um ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vv ) ; ijarg=ijarg+1
     CASE ('-vm'       ) ; CALL getarg( ijarg, cf_vm ) ; ijarg=ijarg+1
     CASE ('-w'        ) ; CALL getarg( ijarg, cf_ww ) ; ijarg=ijarg+1
     CASE ('-wm'       ) ; CALL getarg( ijarg, cf_wm ) ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_ssh) ; ijarg=ijarg+1
     CASE ('-sshm'     ) ; CALL getarg( ijarg, cf_sshm); ijarg=ijarg+1
        !
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh   ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz   ) ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask ) ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy ) ; ijarg=ijarg+1
        ! options
     CASE ( '-full' ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke1'  ) ; CALL getarg(ijarg, cf_out_ke1); ijarg = ijarg + 1
     CASE ( '-o_ke2'  ) ; CALL getarg(ijarg, cf_out_ke2); ijarg = ijarg + 1
     CASE ( '-nc4'  ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  !-- get dimensions (all files must have the same dimension that U-file) --
  jpiglo = getdim (cf_uu, cn_x)
  jpjglo = getdim (cf_uu, cn_y)
  jpk    = getdim (cf_uu, cn_z)
  jpt    = getdim (cf_uu, cn_t)
  jpim1 = jpiglo-1
  jpjm1 = jpjglo-1
  jpkm1 = jpk-1
  
  !-- summary --
  PRINT *, 'jpiglo =', jpiglo
  PRINT *, 'jpjglo =', jpjglo
  PRINT *, 'jpk    =', jpk
  PRINT *, 'jpt    =', jpt

  !-- Allocate --
  ! mesh
  ALLOCATE( deptht(jpk)                   , depthu(jpk)                   , depthv(jpk)                    )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)     , nav_lat_t(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)     , nav_lat_u(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)     , nav_lat_v(jpiglo, jpjglo)      )
  ALLOCATE( ht_0(jpiglo, jpjglo)                                           )
  ALLOCATE( e1t(jpiglo, jpjglo)           , e2t(jpiglo, jpjglo)            )
  ALLOCATE( e1u(jpiglo, jpjglo)           , e2u(jpiglo, jpjglo)            )
  ALLOCATE( e1v(jpiglo, jpjglo)           , e2v(jpiglo, jpjglo)            )
  ALLOCATE( e12t(jpiglo, jpjglo)                                           )
  ALLOCATE( r1_e12u(jpiglo, jpjglo)       , r1_e12v(jpiglo, jpjglo)        )
  ALLOCATE( e3t_0(jpiglo, jpjglo)    , e3t(jpiglo, jpjglo)          )
  ALLOCATE( e3u_0( jpiglo, jpjglo)   , e3v_0(jpiglo, jpjglo)     )
  ALLOCATE( e3u(jpiglo, jpjglo)      , e3v(jpiglo, jpjglo)       )
  ALLOCATE( tmask(jpiglo, jpjglo)    , fmask(jpiglo, jpjglo)     )
  ALLOCATE( umask(jpiglo, jpjglo)    , vmask(jpiglo, jpjglo)     )
  !! variables
  ALLOCATE( sshn(jpiglo, jpjglo)                                           )
  ALLOCATE( un(jpiglo, jpjglo, jpkk)       , vn(jpiglo, jpjglo, jpkk)        )
  ALLOCATE( tmpu(jpiglo, jpjglo, jpkk)     , tmpv(jpiglo, jpjglo, jpkk)    )
  ALLOCATE( wn(jpiglo, jpjglo, jpkk)                                        )
  ALLOCATE( sshnm(jpiglo, jpjglo)                                           )
  ALLOCATE( unm(jpiglo, jpjglo, jpkk)       , vnm(jpiglo, jpjglo, jpkk)        )
  ALLOCATE( wnm(jpiglo, jpjglo, jpkk)                                        )
  !
  ALLOCATE( adv_h_u(jpiglo, jpjglo)  , adv_z_u(jpiglo, jpjglo)   )
  ALLOCATE( adv_h_v(jpiglo, jpjglo)  , adv_z_v(jpiglo, jpjglo)   )
  ALLOCATE( adv_h_ke(jpiglo, jpjglo) , adv_z_ke(jpiglo, jpjglo)  )


  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  nav_lon_t    = getvar(cf_tt, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tt, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tt, cn_vdeptht , jpk)
  nav_lon_u    = getvar(cf_uu, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u    = getvar(cf_uu, 'nav_lat', 1, jpiglo, jpjglo)
  depthu       = getvar1d(cf_uu, cn_vdepthu , jpk)
  nav_lon_v    = getvar(cf_vv, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v    = getvar(cf_vv, 'nav_lat', 1, jpiglo, jpjglo)
  depthv       = getvar1d(cf_vv, cn_vdepthv , jpk)
  ht_0(:,:)    = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  !- hz. mesh -
  e1t(:,:)     = getvar(cf_mh  , 'e1t'  , 1, jpiglo, jpjglo)
  e2t(:,:)     = getvar(cf_mh  , 'e2t'  , 1, jpiglo, jpjglo)
  e1u(:,:)     = getvar(cf_mh  , 'e1u'  , 1, jpiglo, jpjglo)
  e2u(:,:)     = getvar(cf_mh  , 'e2u'  , 1, jpiglo, jpjglo)
  e1v(:,:)     = getvar(cf_mh  , 'e1v'  , 1, jpiglo, jpjglo)
  e2v(:,:)     = getvar(cf_mh  , 'e2v'  , 1, jpiglo, jpjglo)
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))

  !-- Creat output netcdf files to fill in --
  PRINT *, '-- Creat output --'
  CALL CreateOutput

  DO jk = 1, jpkm1
   PRINT *, '-- klayer: ', jk

   !-- load vert. mesh (at rest) and masks (dommsk.f90) --
   e3t_0(:,:) = getvar(cf_mz  , 'e3t_0' , jk, jpiglo, jpjglo)
   e3u_0(:,:) = e3t_0(:,:)
   e3v_0(:,:) = e3t_0(:,:)
   !DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
      DO jj = 1, jpjm1
         DO ji = 1, jpim1   ! vector opt.
            e3u_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji+1,jj) )
            e3v_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji,jj+1) )
         END DO
      END DO
   !END DO
   tmask(:,:) = getvar(cf_mask, 'tmask' , jk, jpiglo, jpjglo)
   umask(:,:) = getvar(cf_mask, 'umask' , jk, jpiglo, jpjglo )
   vmask(:,:) = getvar(cf_mask, 'vmask' , jk, jpiglo, jpjglo )
   !! fmask = 2 on lateral boundaries for no-slip bdy conditions on vorticity !!
   fmask(:,:) = getvar(cf_mask, 'fmask' , jk, jpiglo, jpjglo )

   DO jt = 1,1
     sshn(:,:)     = getvar(cf_ssh  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
     e3t(:,:) = e3t_0(:,:) * (1 + sshn/ht_0)
     !- at u- and v- pts (domvvl.F90) -
     e3u(:,:) = e3u_0(:,:)
     e3v(:,:) = e3v_0(:,:)
     DO jj = 1, jpjm1
        DO ji = 1, jpim1   ! vector opt.
           e3u(ji,jj) = e3u_0(ji,jj) + 0.5_wp * umask(ji,jj) * r1_e12u(ji,jj)                 &
              &                       * (   e12t(ji  ,jj) * ( e3t(ji  ,jj) - e3t_0(ji  ,jj) )    &
              &                           + e12t(ji+1,jj) * ( e3t(ji+1,jj) - e3t_0(ji+1,jj) ) )
           e3v(ji,jj) = e3v_0(ji,jj) + 0.5_wp * vmask(ji,jj) * r1_e12v(ji,jj)                 &
              &                       * (   e12t(ji,jj  ) * ( e3t(ji,jj  ) - e3t_0(ji,jj  ) )    &
              &                           + e12t(ji,jj+1) * ( e3t(ji,jj+1) - e3t_0(ji,jj+1) ) )
        END DO
     END DO

     !-- Load variables --
     IF ( jk == 1 ) THEN
        !- variable -
        un(:,:,jkkm1) = 0._wp
        vn(:,:,jkkm1) = 0._wp
        wn(:,:,jkkm1) = 0._wp
        un(:,:,jkk  ) = getvar(cf_uu, cn_vozocrtx, jk  , jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkk  ) = getvar(cf_vv, cn_vomecrty, jk  , jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkk  ) = getvar(cf_ww, cn_vovecrtz, jk  , jpiglo, jpjglo, ktime=jt )
        un(:,:,jkkp1) = getvar(cf_uu, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkkp1) = getvar(cf_vv, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkkp1) = getvar(cf_ww, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
        !- ensemble mean -
        unm(:,:,jkkm1) = 0._wp
        vnm(:,:,jkkm1) = 0._wp
        wnm(:,:,jkkm1) = 0._wp
        unm(:,:,jkk  ) = getvar(cf_um, cn_vozocrtx, jk  , jpiglo, jpjglo, ktime=jt )
        vnm(:,:,jkk  ) = getvar(cf_vm, cn_vomecrty, jk  , jpiglo, jpjglo, ktime=jt )
        wnm(:,:,jkk  ) = getvar(cf_wm, cn_vovecrtz, jk  , jpiglo, jpjglo, ktime=jt )
        unm(:,:,jkkp1) = getvar(cf_um, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vnm(:,:,jkkp1) = getvar(cf_vm, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wnm(:,:,jkkp1) = getvar(cf_wm, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
     ELSE
        !- variable -
        un(:,:,jkkm1) = un(:,:,jkk  )
        vn(:,:,jkkm1) = vn(:,:,jkk  )
        wn(:,:,jkkm1) = wn(:,:,jkk  )
        un(:,:,jkk  ) = un(:,:,jkkp1)
        vn(:,:,jkk  ) = vn(:,:,jkkp1)
        wn(:,:,jkk  ) = wn(:,:,jkkp1)
        un(:,:,jkkp1) = getvar(cf_uu, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkkp1) = getvar(cf_vv, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkkp1) = getvar(cf_ww, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
        !- ensemble mean -
        unm(:,:,jkkm1) = unm(:,:,jkk  )
        vnm(:,:,jkkm1) = vnm(:,:,jkk  )
        wnm(:,:,jkkm1) = wnm(:,:,jkk  )
        unm(:,:,jkk  ) = unm(:,:,jkkp1)
        vnm(:,:,jkk  ) = vnm(:,:,jkkp1)
        wnm(:,:,jkk  ) = wnm(:,:,jkkp1)
        unm(:,:,jkkp1) = getvar(cf_um, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vnm(:,:,jkkp1) = getvar(cf_vm, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wnm(:,:,jkkp1) = getvar(cf_wm, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
     ENDIF

     tmpu(:,:,:) = 0._wp       ! to remove the dissipation term in the UBS scheme
     tmpv(:,:,:) = 0._wp       ! (care must be taken to not double counting it if considered)

     !-- Advection and KE trends --
     ! full
     !CALL dyn_adv_ubs( jt, jk, tmpu, tmpv, un, vn, wn, un, vn )
     !CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un, vn ) 
     !CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un, vn ) 
     !ierr = putvar(ncout_u , id_varout_u(5) , adv_h_u(:,:) +adv_z_u(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     !ierr = putvar(ncout_v , id_varout_v(5) , adv_h_v(:,:) +adv_z_v(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     !ierr = putvar(ncout_ke, id_varout_ke(5), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     ! ubub 
     CALL dyn_adv_ubs( jt, jk, tmpu, tmpv, unm, vnm, wnm, unm, vnm )
     ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) +adv_z_u(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) +adv_z_v(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )       !ubar*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )       !ubar*
     ierr = putvar(ncout_ke1, id_varout_ke1(1), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )     !uprime*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )     !uprime*
     ierr = putvar(ncout_ke2, id_varout_ke2(1), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     ! ubup
     CALL dyn_adv_ubs( jt, jk, tmpu, tmpv, unm, vnm, wnm, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
     ierr = putvar(ncout_u , id_varout_u(2) , adv_h_u(:,:) +adv_z_u(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     ierr = putvar(ncout_v , id_varout_v(2) , adv_h_v(:,:) +adv_z_v(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )     !uprime* 
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )     !uprime*
     ierr = putvar(ncout_ke1, id_varout_ke1(2), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )      !ubar*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )      !ubar*
     ierr = putvar(ncout_ke2, id_varout_ke2(2), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     ! upub
     CALL dyn_adv_ubs( jt, jk, tmpu, tmpv, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:), wn(:,:,:)-wnm(:,:,:), unm, vnm )
     ierr = putvar(ncout_u , id_varout_u(3) , adv_h_u(:,:) +adv_z_u(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     ierr = putvar(ncout_v , id_varout_v(3) , adv_h_v(:,:) +adv_z_v(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )      !ubar* 
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )      !ubar*
     ierr = putvar(ncout_ke1, id_varout_ke1(3), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )    !uprime*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )    !uprime*
     ierr = putvar(ncout_ke2, id_varout_ke2(3), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     ! upup
     CALL dyn_adv_ubs( jt, jk, tmpu, tmpv, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:), wn(:,:,:)-wnm(:,:,:), &
             & un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
     ierr = putvar(ncout_u , id_varout_u(4) , adv_h_u(:,:) +adv_z_u(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     ierr = putvar(ncout_v , id_varout_v(4) , adv_h_v(:,:) +adv_z_v(:,:) , jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )    !uprime*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )    !uprime*
     ierr = putvar(ncout_ke1, id_varout_ke1(4), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )       !ubar*
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )       !ubar*
     ierr = putvar(ncout_ke2, id_varout_ke2(4), adv_h_ke(:,:)+adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=1 )

   ENDDO       !jt-loop

  ENDDO        !jk-loop
 
  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke1)
  ierr = closeout(ncout_ke2)

CONTAINS

   SUBROUTINE dyn_adv_ubs( kt, jk, u0, v0, u1, v1, w1, u2, v2 )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE dyn_adv_ubs  ***
!!
!! ** Purpose :   Compute the now momentum advection trend in flux form
!!              and the general trend of the momentum equation.
!!
!! ** Method  :   The scheme is the one implemeted in ROMS. It depends
!!      on two parameter gamma1 and gamma2. The former control the
!!      upstream baised part of the scheme and the later the centred
!!      part:     gamma1 = 0    pure centered  (no diffusive part)
!!                       = 1/4  Quick scheme
!!                       = 1/3  3rd order Upstream biased scheme
!!                gamma2 = 0    2nd order finite differencing
!!                       = 1/32 4th order finite differencing
!!      For stability reasons, the first term of the fluxes which cor-
!!      responds to a second order centered scheme is evaluated using
!!      the now velocity (centered in time) while the second term which
!!      is the diffusive part of the scheme, is evaluated using the
!!      before velocity (forward in time).
!!      Default value (hard coded in the begining of the module) are
!!      gamma1=1/3 and gamma2=1/32.
!!
!! ** Action : - (ua,va) updated with the 3D advective momentum trends
!!
!! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling.
!!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step index
      INTEGER, INTENT(in) ::   jk     ! ocean vertical level
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u0, v0   ! U and V ocean velocities used to estimate the 
                                                              ! upstream diffusive fluxes (should be full velocities)
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u1, v1, w1   ! U, V and W ocean advecting velocities
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u2, v2       ! U and V    ocean advected  velocities
!
      INTEGER  ::   ji, jj              ! dummy loop indices
      REAL(wp) ::   zbu, zbv    ! temporary scalars
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu, zfv, zfu2, zfv2
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfw
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu_t, zfv_t, zfu_f, zfv_f
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfu_uw, zfv_vw
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zlu_uu, zlv_vv, zlu_uv, zlv_vu
!!----------------------------------------------------------------------
!
      ALLOCATE( zfu(jpiglo, jpjglo)       , zfv(jpiglo, jpjglo)       )
      ALLOCATE( zfu2(jpiglo, jpjglo)      , zfv2(jpiglo, jpjglo)      )
      ALLOCATE( zfw(jpiglo, jpjglo, jpkk)                             )
      ALLOCATE( zfu_t(jpiglo, jpjglo)     , zfv_t(jpiglo, jpjglo)     )
      ALLOCATE( zfu_f(jpiglo, jpjglo)     , zfv_f(jpiglo, jpjglo)     )
      ALLOCATE( zfu_uw(jpiglo, jpjglo, jpkk)    , zfv_vw(jpiglo, jpjglo, jpkk)    )
      !
      ALLOCATE( zlu_uu(jpiglo, jpjglo, 2) , zlv_vv(jpiglo, jpjglo, 2) )
      ALLOCATE( zlu_uv(jpiglo, jpjglo, 2) , zlv_vu(jpiglo, jpjglo, 2) )

!
      zfu_t(:,:)    = 0._wp
      zfv_t(:,:)    = 0._wp
      zfu_f(:,:)    = 0._wp
      zfv_f(:,:)    = 0._wp
      zfu_uw(:,:,:) = 0._wp
      zfv_vw(:,:,:) = 0._wp
!
      zlu_uu(:,:,:) = 0._wp
      zlv_vv(:,:,:) = 0._wp
      zlu_uv(:,:,:) = 0._wp
      zlv_vu(:,:,:) = 0._wp

!                                      ! =========================== !
                                       !  Laplacian of the velocity  !
!                                      ! =========================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = e2u(:,:) * e3u(:,:) * u1(:,:,jkk)
         zfv(:,:) = e1v(:,:) * e3v(:,:) * v1(:,:,jkk)
!
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = 2, jpim1   ! vector opt.
!
               zlu_uu(ji,jj,1) = ( u0 (ji+1,jj  ,jkk) - 2.*u0 (ji,jj,jkk) + u0 (ji-1,jj  ,jkk) ) * umask(ji,jj)
               zlv_vv(ji,jj,1) = ( v0 (ji  ,jj+1,jkk) - 2.*v0 (ji,jj,jkk) + v0 (ji  ,jj-1,jkk) ) * vmask(ji,jj)
               zlu_uv(ji,jj,1) = ( u0 (ji  ,jj+1,jkk) - u0 (ji  ,jj  ,jkk) ) * fmask(ji  ,jj  )   &
                  &               - ( u0 (ji  ,jj  ,jkk) - u0 (ji  ,jj-1,jkk) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,1) = ( v0 (ji+1,jj  ,jkk) - v0 (ji  ,jj  ,jkk) ) * fmask(ji  ,jj  )   &
                  &               - ( v0 (ji  ,jj  ,jkk) - v0 (ji-1,jj  ,jkk) ) * fmask(ji-1,jj  )
!
               zlu_uu(ji,jj,2) = ( zfu(ji+1,jj  ) - 2.*zfu(ji,jj) + zfu(ji-1,jj  ) ) * umask(ji,jj)
               zlv_vv(ji,jj,2) = ( zfv(ji  ,jj+1) - 2.*zfv(ji,jj) + zfv(ji  ,jj-1) ) * vmask(ji,jj)
               zlu_uv(ji,jj,2) = ( zfu(ji  ,jj+1) - zfu(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfu(ji  ,jj  ) - zfu(ji  ,jj-1) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,2) = ( zfv(ji+1,jj  ) - zfv(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfv(ji  ,jj  ) - zfv(ji-1,jj  ) ) * fmask(ji-1,jj  )
            END DO
         END DO

!                                      ! ====================== !
!                                      !  Horizontal advection  !
                                       ! ====================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = 0.25 * e2u(:,:) * e3u(:,:) * u1(:,:,jkk)
         zfv(:,:) = 0.25 * e1v(:,:) * e3v(:,:) * v1(:,:,jkk)
         zfu2(:,:) = 0.25 * e2u(:,:) * e3u(:,:) * u0(:,:,jkk)
         zfv2(:,:) = 0.25 * e1v(:,:) * e3v(:,:) * v0(:,:,jkk)
!
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, jpim1   ! vector opt.
               zui = ( u0(ji,jj,jkk) + u0(ji+1,jj  ,jkk) )
               zvj = ( v0(ji,jj,jkk) + v0(ji  ,jj+1,jkk) )
               IF (zui > 0) THEN   ;   zl_u = zlu_uu(ji  ,jj,1)
               ELSE                ;   zl_u = zlu_uu(ji+1,jj,1)
               ENDIF
               IF (zvj > 0) THEN   ;   zl_v = zlv_vv(ji,jj  ,1)
               ELSE                ;   zl_v = zlv_vv(ji,jj+1,1)
               ENDIF
               zui = ( u2(ji,jj,jkk) + u2(ji+1,jj  ,jkk) )
               zvj = ( v2(ji,jj,jkk) + v2(ji  ,jj+1,jkk) )
!
               zfu_t(ji+1,jj  ) = ( zfu(ji,jj) + zfu(ji+1,jj  )                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,2) + zlu_uu(ji+1,jj  ,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1) = ( zfv(ji,jj) + zfv(ji  ,jj+1)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,2) + zlv_vv(ji  ,jj+1,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
!
               zfuj = ( zfu2(ji,jj) + zfu2(ji  ,jj+1) )
               zfvi = ( zfv2(ji,jj) + zfv2(ji+1,jj  ) )
               IF (zfuj > 0) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,1)
               ELSE                 ;    zl_v = zlv_vu( ji+1,jj,1)
               ENDIF
               IF (zfvi > 0) THEN   ;    zl_u = zlu_uv( ji,jj  ,1)
               ELSE                 ;    zl_u = zlu_uv( ji,jj+1,1)
               ENDIF
               zfuj = ( zfu(ji,jj) + zfu(ji  ,jj+1) )
               zfvi = ( zfv(ji,jj) + zfv(ji+1,jj  ) )
!
               zfv_f(ji  ,jj  ) = ( zfvi - gamma2 * ( zlv_vu(ji,jj,2) + zlv_vu(ji+1,jj  ,2) )  )   &
                  &                * ( u2(ji,jj,jkk) + u2(ji  ,jj+1,jkk) - gamma1 * zl_u )
               zfu_f(ji  ,jj  ) = ( zfuj - gamma2 * ( zlu_uv(ji,jj,2) + zlu_uv(ji  ,jj+1,2) )  )   &
                  &                * ( v2(ji,jj,jkk) + v2(ji+1,jj  ,jkk) - gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = 2, jpim1   ! vector opt.
               zbu = e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj)
               zbv = e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj)
!
               adv_h_u(ji,jj)  = - (  zfu_t(ji+1,jj  ) - zfu_t(ji  ,jj  )    &
                  &                   + zfv_f(ji  ,jj  ) - zfv_f(ji  ,jj-1)  ) / zbu
               adv_h_v(ji,jj) = - (  zfu_f(ji  ,jj  ) - zfu_f(ji-1,jj  )    &
                  &                  + zfv_t(ji  ,jj+1) - zfv_t(ji  ,jj  )  ) / zbv
!
               adv_h_u(ji,jj) = adv_h_u(ji,jj) * umask(ji,jj)
               adv_h_v(ji,jj) = adv_h_v(ji,jj) * vmask(ji,jj)
            END DO
         END DO

!                                      ! ==================== !
!                                      !  Vertical advection  !
                                       ! ==================== !
!                                         ! Vertical volume fluxes
         zfw(:,:,jkk) = 0.25 * e1t(:,:) * e2t(:,:) * w1(:,:,jkk)
         zfw(:,:,jkkp1) = 0.25 * e1t(:,:) * e2t(:,:) * w1(:,:,jkkp1)
!
   !      IF( jk == 1 ) THEN                        ! surface/bottom advective fluxes
   !       ... moved after interior fluxes ...
   !      ELSE                                      ! interior fluxes
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zfu_uw(ji,jj,jkk) = ( zfw(ji,jj,jkk)+ zfw(ji+1,jj  ,jkk) ) * ( u2(ji,jj,jkk) + u2(ji,jj,jkkm1) )
                  zfv_vw(ji,jj,jkk) = ( zfw(ji,jj,jkk)+ zfw(ji  ,jj+1,jkk) ) * ( v2(ji,jj,jkk) + v2(ji,jj,jkkm1) )
                  zfu_uw(ji,jj,jkkp1) = ( zfw(ji,jj,jkkp1)+ zfw(ji+1,jj  ,jkkp1) ) * ( u2(ji,jj,jkkp1) + u2(ji,jj,jkk) )
                  zfv_vw(ji,jj,jkkp1) = ( zfw(ji,jj,jkkp1)+ zfw(ji  ,jj+1,jkkp1) ) * ( v2(ji,jj,jkkp1) + v2(ji,jj,jkk) )
               END DO
            END DO
    !     ENDIF
         IF( jk == jpkm1 ) THEN
            zfu_uw(:,:,jkkp1) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vw(:,:,jkkp1) = 0.e0
         ENDIF
         IF ( jk == 1 ) THEN                        ! Surface value :
!            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uw(:,:, jkk ) = 0.e0
               zfv_vw(:,:, jkk ) = 0.e0
!            ELSE                                             ! constant volume : advection through the surface
!               ...
!            ENDIF
         ENDIF

         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               adv_z_u(ji,jj) =  - ( zfu_uw(ji,jj,jkk) - zfu_uw(ji,jj,jkkp1) )    &
                  &  / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj) )
               adv_z_v(ji,jj) =  - ( zfv_vw(ji,jj,jkk) - zfv_vw(ji,jj,jkkp1) )    &
                  &  / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj) )
!
               adv_z_u(ji,jj) = adv_z_u(ji,jj) * umask(ji,jj)
               adv_z_v(ji,jj) = adv_z_v(ji,jj) * vmask(ji,jj)
            END DO
         END DO
!
      DEALLOCATE( zfu_t , zfv_t  )
      DEALLOCATE( zfu_f , zfv_f  )
      DEALLOCATE( zfu   , zfv    )
      DEALLOCATE( zfw            )
      DEALLOCATE( zlu_uu, zlv_vv )
      DEALLOCATE( zlu_uv, zlv_vu )

   END SUBROUTINE dyn_adv_ubs


   SUBROUTINE trd_ken( putrd, pvtrd, ktrd, u0, v0 )
!!---------------------------------------------------------------------
!!                  ***  ROUTINE trd_ken  ***
!!
!! ** Purpose :   output 3D Kinetic Energy trends using IOM
!!
!! ** Method  : - apply lbc to the input masked velocity trends
!!              - compute the associated KE trend:
!!          zke = 0.5 * (  mi-1[ un * putrd * bu ] + mj-1[ vn * pvtrd * bv]  ) / bt
!!      where bu, bv, bt are the volume of u-, v- and t-boxes.
!!              - vertical diffusion case (jpdyn_zdf):
!!          diagnose separately the KE trend associated with wind stress
!!              - bottom friction case (jpdyn_bfr):
!!          explicit case (ln_bfrimp=F): bottom trend put in the 1st level
!!                                       of putrd, pvtrd
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   putrd, pvtrd   ! U and V masked trends
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u0, v0       ! U and V 
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ktrd           ! KE trend
!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
!
      REAL(wp)                                :: rau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   bu, bv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   r1_bt    ! inverse of t-box volume
!!----------------------------------------------------------------------
!
!
      ALLOCATE( bu(jpiglo,jpjglo) , bv(jpiglo,jpjglo) , r1_bt(jpiglo,jpjglo) )
!
!      IF ( lk_vvl .AND. kt /= nkstp ) THEN   ! Variable volume: set box volume at the 1st call of kt time step
!         nkstp = kt
      bu   (:,:) =  e1u(:,:) * e2u(:,:) * e3u(:,:)
      bv   (:,:) =  e1v(:,:) * e2v(:,:) * e3v(:,:)
      r1_bt(:,:) = 1._wp / ( e12t(:,:) * e3t(:,:) ) * tmask(:,:)
!
      ktrd(:,:) = 0._wp
      ktrd(:,:) = 0._wp
      DO jj = 2, jpjglo
         DO ji = 2, jpiglo
            ktrd(ji,jj) = 0.5_wp * rau0 *( u0(ji  ,jj,jkk) * putrd(ji  ,jj) * bu(ji  ,jj)  &
               &                            + u0(ji-1,jj,jkk) * putrd(ji-1,jj) * bu(ji-1,jj)  &
               &                            + v0(ji,jj  ,jkk) * pvtrd(ji,jj  ) * bv(ji,jj  )  &
               &                            + v0(ji,jj-1,jkk) * pvtrd(ji,jj-1) * bv(ji,jj-1)  ) * r1_bt(ji,jj)
         END DO
      END DO

   END SUBROUTINE trd_ken



  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s)
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! define new variables for output
    ipk(:)                         = jpk
    stypvar1(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar1(1)%cname              = 'adv_ubub'
    stypvar1(1)%cunits             = 'm/s^2'
    stypvar1(1)%rmissing_value     = 99999.
    stypvar1(1)%valid_min          = -1.
    stypvar1(1)%valid_max          = 1.
    stypvar1(1)%clong_name         = 'Advection of mean flow by the mean flow'
    stypvar1(1)%cshort_name        = 'adv_ubub'
    stypvar1(1)%conline_operation  = 'On u-grid'
    stypvar1(1)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar1(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar1(2)%cname              = 'adv_ubup'
    stypvar1(2)%cunits             = 'm/s^2'
    stypvar1(2)%rmissing_value     = 99999.
    stypvar1(2)%valid_min          = -1.
    stypvar1(2)%valid_max          = 1.
    stypvar1(2)%clong_name         = 'Advection of the perturbations by the mean flow'
    stypvar1(2)%cshort_name        = 'adv_ubup'
    stypvar1(2)%conline_operation  = 'On u-grid'
    stypvar1(2)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar1(3)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar1(3)%cname              = 'adv_upub'
    stypvar1(3)%cunits             = 'm/s^2'
    stypvar1(3)%rmissing_value     = 99999.
    stypvar1(3)%valid_min          = -1.
    stypvar1(3)%valid_max          = 1.
    stypvar1(3)%clong_name         = 'Advection of the the mean flow by the perturbations'
    stypvar1(3)%cshort_name        = 'adv_upub'
    stypvar1(3)%conline_operation  = 'On u-grid'
    stypvar1(3)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar1(4)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar1(4)%cname              = 'adv_upup'
    stypvar1(4)%cunits             = 'm/s^2'
    stypvar1(4)%rmissing_value     = 99999.
    stypvar1(4)%valid_min          = -1.
    stypvar1(4)%valid_max          = 1.
    stypvar1(4)%clong_name         = 'Advection of the perturbations by perturbations'
    stypvar1(4)%cshort_name        = 'adv_upup'
    stypvar1(4)%conline_operation  = 'On u-grid'
    stypvar1(4)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    ipk(:)                         = jpk
    stypvar2(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname              = 'adv_vbvb'
    stypvar2(1)%cunits             = 'm/s^2'
    stypvar2(1)%rmissing_value     = 99999.
    stypvar2(1)%valid_min          = -1.
    stypvar2(1)%valid_max          = 1.
    stypvar2(1)%clong_name         = 'Advection of mean flow by the mean flow'
    stypvar2(1)%cshort_name        = 'adv_ubub'
    stypvar2(1)%conline_operation  = 'On v-grid'
    stypvar2(1)%caxis              = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(2)%cname              = 'adv_vbvp'
    stypvar2(2)%cunits             = 'm/s^2'
    stypvar2(2)%rmissing_value     = 99999.
    stypvar2(2)%valid_min          = -1.
    stypvar2(2)%valid_max          = 1.
    stypvar2(2)%clong_name         = 'Advection of the perturbations by the mean flow'
    stypvar2(2)%cshort_name        = 'adv_vbvp'
    stypvar2(2)%conline_operation  = 'On v-grid'
    stypvar2(2)%caxis              = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(3)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(3)%cname              = 'adv_vpvb'
    stypvar2(3)%cunits             = 'm/s^2'
    stypvar2(3)%rmissing_value     = 99999.
    stypvar2(3)%valid_min          = -1.
    stypvar2(3)%valid_max          = 1.
    stypvar2(3)%clong_name         = 'Advection of the the mean flow by the perturbations'
    stypvar2(3)%cshort_name        = 'adv_vpvb'
    stypvar2(3)%conline_operation  = 'On v-grid'
    stypvar2(3)%caxis              = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(4)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(4)%cname              = 'adv_vpvp'
    stypvar2(4)%cunits             = 'm/s^2'
    stypvar2(4)%rmissing_value     = 99999.
    stypvar2(4)%valid_min          = -1.
    stypvar2(4)%valid_max          = 1.
    stypvar2(4)%clong_name         = 'Advection of the perturbations by perturbations'
    stypvar2(4)%cshort_name        = 'adv_vpvp'
    stypvar2(4)%conline_operation  = 'On v-grid'
    stypvar2(4)%caxis              = 'time depthv nav_lon_v nav_lat_v'

    ipk(:)                         = jpk
    stypvar3(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname              = 'advm_mke'
    stypvar3(1)%cunits             = 'm/s^2'
    stypvar3(1)%rmissing_value     = 99999.
    stypvar3(1)%valid_min          = -1.
    stypvar3(1)%valid_max          = 1.
    stypvar3(1)%clong_name         = 'Advection of MKE by the mean flow'
    stypvar3(1)%cshort_name        = 'advm_mke'
    stypvar3(1)%conline_operation  = 'On t-grid'
    stypvar3(1)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(2)%cname              = 'advm_eke'
    stypvar3(2)%cunits             = 'm/s^2'
    stypvar3(2)%rmissing_value     = 99999.
    stypvar3(2)%valid_min          = -1.
    stypvar3(2)%valid_max          = 1.
    stypvar3(2)%clong_name         = 'Advection of EKE by the mean flow'
    stypvar3(2)%cshort_name        = 'advm_eke'
    stypvar3(2)%conline_operation  = 'On t-grid'
    stypvar3(2)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(3)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(3)%cname              = 'advp_mke'
    stypvar3(3)%cunits             = 'm/s^2'
    stypvar3(3)%rmissing_value     = 99999.
    stypvar3(3)%valid_min          = -1.
    stypvar3(3)%valid_max          = 1.
    stypvar3(3)%clong_name         = 'Advection of MKE by the perturbations'
    stypvar3(3)%cshort_name        = 'advp_mke'
    stypvar3(3)%conline_operation  = 'On t-grid'
    stypvar3(3)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(4)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(4)%cname              = 'advp_eke'
    stypvar3(4)%cunits             = 'm/s^2'
    stypvar3(4)%rmissing_value     = 99999.
    stypvar3(4)%valid_min          = -1.
    stypvar3(4)%valid_max          = 1.
    stypvar3(4)%clong_name         = 'Advection of EKE by perturbations'
    stypvar3(4)%cshort_name        = 'advp_eke'
    stypvar3(4)%conline_operation  = 'On t-grid'
    stypvar3(4)%caxis              = 'time deptht nav_lon_t nav_lat_t'

    ipk(:)                         = jpk
    stypvar4(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(1)%cname              = 'adv_eddy1'
    stypvar4(1)%cunits             = 'm/s^2'
    stypvar4(1)%rmissing_value     = 99999.
    stypvar4(1)%valid_min          = -1.
    stypvar4(1)%valid_max          = 1.
    stypvar4(1)%clong_name         = 'ubar * nabla cdot (\mathbf{ubar} * uprime )'
    stypvar4(1)%cshort_name        = 'advm_eddy1'
    stypvar4(1)%conline_operation  = 'On t-grid'
    stypvar4(1)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(2)%cname              = 'advm_eddy2'
    stypvar4(2)%cunits             = 'm/s^2'
    stypvar4(2)%rmissing_value     = 99999.
    stypvar4(2)%valid_min          = -1.
    stypvar4(2)%valid_max          = 1.
    stypvar4(2)%clong_name         = 'uprime * nabla cdot (\mathbf{ubar} * ubar )'
    stypvar4(2)%cshort_name        = 'advm_eddy2'
    stypvar4(2)%conline_operation  = 'On t-grid'
    stypvar4(2)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(3)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(3)%cname              = 'adv_eddyFlx'
    stypvar4(3)%cunits             = 'm/s^2'
    stypvar4(3)%rmissing_value     = 99999.
    stypvar4(3)%valid_min          = -1.
    stypvar4(3)%valid_max          = 1.
    stypvar4(3)%clong_name         = 'Eddy momentum flux (uprime*uprime*d(ubar)/dx)'
    stypvar4(3)%cshort_name        = 'adv_eddyFlx'
    stypvar4(3)%conline_operation  = 'On t-grid'
    stypvar4(3)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(4)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(4)%cname              = 'adv_mec'
    stypvar4(4)%cunits             = 'm/s^2'
    stypvar4(4)%rmissing_value     = 99999.
    stypvar4(4)%valid_min          = -1.
    stypvar4(4)%valid_max          = 1.
    stypvar4(4)%clong_name         = 'Mean-to-Eddy energy Conversion (ubar*d(upr*upr)/dx)'
    stypvar4(4)%cshort_name        = 'adv_mec'
    stypvar4(4)%conline_operation  = 'On t-grid'
    stypvar4(4)%caxis              = 'time deptht nav_lon_t nav_lat_t'


    ! create output fileset
    ncout_u = create      (cf_out_u, cf_uu ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar1,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_uu ,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u, depthu   )

    ncout_v = create      (cf_out_v, cf_vv ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  pnvarout, ipk , id_varout_v           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vv ,  jpiglo, jpjglo, jpk, nav_lon_v, nav_lat_v, depthv   )

    ncout_ke1= create     (cf_out_ke1, cf_tt ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke1 , stypvar3,  pnvarout, ipk , id_varout_ke1          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke1 , cf_tt ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    ncout_ke2= create     (cf_out_ke2, cf_tt ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke2 , stypvar4,  pnvarout, ipk , id_varout_ke2         , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke2 , cf_tt ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )


    dtim = getvar1d(cf_uu , cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u , dtim,        1, 'T')
    ierr = putvar1d(ncout_v , dtim,        1, 'T')
    ierr = putvar1d(ncout_ke1, dtim,       1, 'T')
    ierr = putvar1d(ncout_ke2, dtim,       1, 'T')


  END SUBROUTINE CreateOutput


END PROGRAM
