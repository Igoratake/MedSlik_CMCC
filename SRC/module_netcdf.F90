!-----------------------------------------------------------------------------------
!  module_netcdf.F90
!-----------------------------------------------------------------------------------
!
!
!   --|-------------------------------------------------------------|--
!     | National Institute of Geopyisics and Volcanology (I.N.G.V)  |
!     | Bologna Section                                             |
!     | Via Aldo Moro 44, Bologna, Italy                            |
!     |                                                             |
!     | Programmer: Diego Bruciaferri                               |
!   --|-------------------------------------------------------------|--
!
!
!---------------------------------------------------------------------------------
!
!  This routine module specifies netcdf variables and subroutines needed to 
!  read/write netcdf files involved in MEDSLIK_II netcdf output creation. 
!  
!---------------------------------------------------------------------------------
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!---------------------------------------------------------------------------------

MODULE MODULE_NETCDF

      USE netcdf
      IMPLICIT NONE

! -------------------------------- SEA/LAND MASK NC VARIABLES ------------------------------            

      CHARACTER(LEN=*),PARAMETER            :: MaskFile='../RUN/data/LandSeaMask_150m.nc'
      INTEGER                               :: MaskFileid, maskid, lon_mask_id, lat_mask_id
      INTEGER                               :: xid, nx, yid, ny, STATUS
      CHARACTER(LEN=*),PARAMETER            :: MaskName = "LSM", lon_mask = "Lon", lat_mask = "Lat"
      REAL,DIMENSION(:),ALLOCATABLE         :: lon_msk, lat_msk
      REAL,DIMENSION(:,:),ALLOCATABLE       :: msk

! ------------------------------ MEDSLIK-II NetCDF OUTPUT VARIABLES -------------------------

! COMMON VARIABLES

      CHARACTER(LEN=150)                                  :: nc_file
      INTEGER                                             :: ncid,rec
      INTEGER,PARAMETER                                   :: ndims = 3
      INTEGER                                             :: ntim, tim_dimid, tim_dims(1)
      INTEGER                                             :: count_tim(1), start_tim(1)
      INTEGER                                             :: tim_varid
      REAL                                                :: time_out
      CHARACTER(LEN=*),PARAMETER                          :: TIM_NAME = "time"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_TIM_NAME = "hours_from_the_start_of_&
                                                                              the_simulation"
      CHARACTER(LEN=*),PARAMETER                          :: TIM_UNITS = "hr"

      CHARACTER(LEN=*),PARAMETER                          :: LONG_NAME = "long_name"
      CHARACTER(LEN=*),PARAMETER                          :: UNITS = "units"
      CHARACTER(LEN=*),PARAMETER                          :: LAT_UNITS = "degrees_north"
      CHARACTER(LEN=*),PARAMETER                          :: LON_UNITS = "degrees_east"
      CHARACTER(LEN=*),PARAMETER                          :: VEL_UNITS = "m/s"
      CHARACTER(LEN=*),PARAMETER                          :: COMM = "comment"
      CHARACTER(LEN=*),PARAMETER                          :: COOR = "coordinate"

! OIL-GRID VARIABLES

      INTEGER                                             :: nx_oil, ny_oil
      INTEGER                                             :: lat_oil_dimid, lon_oil_dimid, pnt_cst_dimid
      INTEGER                                             :: lat_oil_dims(1), lon_oil_dims(1)
      CHARACTER(LEN=*),PARAMETER                          :: X_OIL_NAME = "x_oil"
      CHARACTER(LEN=*),PARAMETER                          :: Y_OIL_NAME = "y_oil"
      CHARACTER(LEN=*),PARAMETER                          :: P_CST_NAME = "coastal_segments"
      CHARACTER(LEN=*),PARAMETER                          :: LAT_OIL_NAME = "lat_oil"
      CHARACTER(LEN=*),PARAMETER                          :: LON_OIL_NAME = "lon_oil"
      INTEGER                                             :: start_oil(ndims), count_oil(ndims)
      INTEGER                                             :: start_cst(2), count_cst(2)
      INTEGER                                             :: lon_oil_varid, lat_oil_varid

      CHARACTER(LEN=*),PARAMETER                          :: SURF_NAME="surface"
      CHARACTER(LEN=*),PARAMETER                          :: COAST_NAME="coast"
      CHARACTER(LEN=*),PARAMETER                          :: DISPL_NAME="displaced"
      CHARACTER(LEN=*),PARAMETER                          :: U_OIL_NAME = "u_oil"
      CHARACTER(LEN=*),PARAMETER                          :: V_OIL_NAME = "v_oil"
      CHARACTER(LEN=*),PARAMETER                          :: UW_OIL_NAME = "uw_oil"
      CHARACTER(LEN=*),PARAMETER                          :: VW_OIL_NAME = "vw_oil"

      INTEGER                                             :: surf_varid, coast_varid, displ_varid
      INTEGER                                             :: u_oil_varid, v_oil_varid
      INTEGER                                             :: uw_oil_varid, vw_oil_varid
      INTEGER                                             :: surf_dims(ndims), coast_dims(2), displ_dims(ndims)
      INTEGER                                             :: u_oil_dims(ndims), v_oil_dims(ndims)
      INTEGER                                             :: uw_oil_dims(ndims), vw_oil_dims(ndims)

      CHARACTER(LEN=*),PARAMETER                          :: LONG_LON_OIL = "Longitude_of_oil_map"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_LAT_OIL = "Latitude_of_oil_map"

      CHARACTER(LEN=*),PARAMETER                          :: LONG_surf = "Concentration_of_oil_on_the_surface"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_coast = "Concentration_of_oil_on_coastal_segments"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_displ = "Concentration_of_oil_displaced_in_&
                                                                           the_water_column"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_U_OIL = "Current_zonal_velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_V_OIL = "Current_meridional_velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_UW_OIL = "Stokes_Drift_zonal_velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_VW_OIL = "Stokes_Drift_meridional_velocity"

      CHARACTER(LEN=*),PARAMETER                          :: OIL_UNITS = "kg/m2"
      CHARACTER(LEN=*),PARAMETER                          :: CST_UNITS = "kg/m2"
    
      CHARACTER(LEN=*),PARAMETER                          :: COMM_LON_OIL = "Longitude_of_the_regular_oil_map_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_LAT_OIL = "Latitude_of_the_regular_oil_map_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_U_OIL = "Current_zonal_velocity_defined_on_the_&
                                                                           oil_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_V_OIL = "Current_meridional_velocity_defined_on_&
                                                                           the_oil_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_UW_OIL = "Stokes_Drift_zonal_velocity_defined_on_&
                                                                            the_oil_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_VW_OIL = "Stokes_Drift_meridional_velocity_defined_&
                                                                            on_the_oil_grid"


      CHARACTER(LEN=*),PARAMETER                          :: COOR_oil = "lat_oil lon_oil"
      CHARACTER(LEN=*),PARAMETER                          :: COOR_cst = "coastal_segments"

      REAL,DIMENSION(:,:),ALLOCATABLE                     :: surf_out, displ_out
      REAL,DIMENSION(:),ALLOCATABLE                       :: coast_out
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: u_oil_out, v_oil_out, uw_oil_out, vw_oil_out 


! CURRENT-GRID VARIABLES

      INTEGER                                             :: nx_cur, ny_cur
      INTEGER                                             :: lat_cur_dimid, lon_cur_dimid
      INTEGER                                             :: lat_cur_dims(1), lon_cur_dims(1)
      CHARACTER(LEN=*),PARAMETER                          :: X_CUR_NAME = "x_cur"
      CHARACTER(LEN=*),PARAMETER                          :: Y_CUR_NAME = "y_cur"
      CHARACTER(LEN=*),PARAMETER                          :: LAT_CUR_NAME = "lat_cur"
      CHARACTER(LEN=*),PARAMETER                          :: LON_CUR_NAME = "lon_cur"
      INTEGER                                             :: start_cur(ndims), count_cur(ndims)
      INTEGER                                             :: lon_cur_varid, lat_cur_varid

      CHARACTER(LEN=*),PARAMETER                          :: U_CUR_NAME = "u_cur"
      CHARACTER(LEN=*),PARAMETER                          :: V_CUR_NAME = "v_cur"

      INTEGER                                             :: u_cur_varid, v_cur_varid
      INTEGER                                             :: u_cur_dims(ndims), v_cur_dims(ndims)

      CHARACTER(LEN=*),PARAMETER                          :: LONG_LON_CUR = "Longitude_of_current_map"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_LAT_CUR = "Latitude_of_current_map"

      CHARACTER(LEN=*),PARAMETER                          :: LONG_U_CUR = "Current_zonal_velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_V_CUR = "Current_meridional_velocity"

      CHARACTER(LEN=*),PARAMETER                          :: COMM_U_CUR = "Current_zonal_velocity_defined_on_&
                                                                           the_current_original_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_V_CUR = "Current_meridional_velocity_defined_&
                                                                           on_the_current_original_grid"

      CHARACTER(LEN=*),PARAMETER                          :: COOR_cur = "lat_cur lon_cur"

      REAL,DIMENSION(:,:),ALLOCATABLE                     :: u_cur_out, v_cur_out

! WAVE-GRID VARIABLES

      INTEGER                                             :: nx_wav, ny_wav
      INTEGER                                             :: lat_wav_dimid, lon_wav_dimid
      INTEGER                                             :: lat_wav_dims(1), lon_wav_dims(1)
      CHARACTER(LEN=*),PARAMETER                          :: X_WAV_NAME = "x_wav"
      CHARACTER(LEN=*),PARAMETER                          :: Y_WAV_NAME = "y_wav"
      CHARACTER(LEN=*),PARAMETER                          :: LAT_WAV_NAME = "lat_wav"
      CHARACTER(LEN=*),PARAMETER                          :: LON_WAV_NAME = "lon_wav"
      INTEGER                                             :: start_wav(ndims), count_wav(ndims)
      INTEGER                                             :: lon_wav_varid, lat_wav_varid

      CHARACTER(LEN=*),PARAMETER                          :: U_WAV_NAME = "u_wav"
      CHARACTER(LEN=*),PARAMETER                          :: V_WAV_NAME = "v_wav"

      INTEGER                                             :: u_wav_varid, v_wav_varid
      INTEGER                                             :: u_wav_dims(ndims), v_wav_dims(ndims)

      CHARACTER(LEN=*),PARAMETER                          :: LONG_LON_WAV = "Longitude_of_wave_map"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_LAT_WAV = "Latitude_of_wave_map"

      CHARACTER(LEN=*),PARAMETER                          :: LONG_U_WAV = "Stokes_Drift_zonal_velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_V_WAV = "Stokes_Drift_meridional_velocity"

      CHARACTER(LEN=*),PARAMETER                          :: COMM_U_WAV = "Stokes_Drift_zonal_velocity_defined_&
                                                                           on_the_current_original_grid"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_V_WAV = "Stokes_Drift_meridional_velocity_defined_&
                                                                           on_the_current_original_grid"

      CHARACTER(LEN=*),PARAMETER                          :: COOR_wav = "lat_wav lon_wav"
      
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: u_wav_out, v_wav_out

! WIND POINTS VARIABLES

      INTEGER                                             :: n_cog
      INTEGER                                             :: cog_dimid
      CHARACTER(LEN=*),PARAMETER                          :: N_COG_NAME = "n_cog"

      INTEGER                                             :: lon_cog_dims(2),lat_cog_dims(2) 
      CHARACTER(LEN=*),PARAMETER                          :: LAT_COG_NAME = "lat_cog"
      CHARACTER(LEN=*),PARAMETER                          :: LON_COG_NAME = "lon_cog"
      INTEGER                                             :: lon_cog_varid, lat_cog_varid
      INTEGER                                             :: start_cog(2), count_cog(2)
      CHARACTER(LEN=*),PARAMETER                          :: LONG_LON_COG ="Longitude_of_slicks_gravity_centre_points"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_LAT_COG ="Latitude_of_slicks_gravity_centre_points"

      CHARACTER(LEN=*),PARAMETER                          :: X_WND_NAME = "x_wnd"
      CHARACTER(LEN=*),PARAMETER                          :: Y_WND_NAME = "y_wnd"
      INTEGER                                             :: x_wnd_varid, y_wnd_varid
      INTEGER                                             :: x_wnd_dims(2), y_wnd_dims(2)
      CHARACTER(LEN=*),PARAMETER                          :: LONG_X_WND = "Wind_zonal velocity"
      CHARACTER(LEN=*),PARAMETER                          :: LONG_Y_WND = "Wind_meridional_velocity"

      CHARACTER(LEN=*),PARAMETER                          :: COMM_X_WND = "Wind_zonal_velocity_defined_&
                                                                           on_slicks_gravity_centre_points"
      CHARACTER(LEN=*),PARAMETER                          :: COMM_Y_WND = "Wind_meridional_velocity_defined_&
                                                                           on_slicks_gravity_centre_points"

      CHARACTER(LEN=*),PARAMETER                          :: COOR_wnd = "lat_cog lon_cog"

      REAL,DIMENSION(:),ALLOCATABLE                       :: lon_cog_out, lat_cog_out, x_wnd_out, y_wnd_out


CONTAINS


!=====================================================================================
!                           SUBROUTINE HANDLE_ERR(STATUS)
!=====================================================================================

     SUBROUTINE HANDLE_ERR(status)

!         USE netcdf
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: status

         IF ( STATUS /= nf90_noerr ) THEN

            WRITE(*,*) nf90_strerror(status)
            STOP 'Stopped'

         ENDIF

     END SUBROUTINE HANDLE_ERR


!===========================================================================================
!                                  SUBROUTINE READ_MASK150
!===========================================================================================
     SUBROUTINE READ_MASK150
 

         STATUS = NF90_OPEN(TRIM(MaskFile), NF90_NOWRITE, MaskFileid)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

!--- Get the ID of the DIMENSIONS of the variables

         STATUS = NF90_INQ_DIMID(MaskFileid, "X", xid)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_INQ_DIMID(MaskFileid, "Y", yid)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

!--- Get the SIZE of the DIMENSIONS of the variables

         STATUS = NF90_INQUIRE_DIMENSION(MaskFileid, xid, len = nx)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_INQUIRE_DIMENSION(MaskFileid, yid, len = ny)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

!---  Get the ID of the variables  --------------

         STATUS = NF90_INQ_VARID(MaskFileid, lon_mask, lon_mask_id)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_INQ_VARID(MaskFileid, lat_mask, lat_mask_id)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_INQ_VARID(MaskFileid, MaskName, maskid)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

!---  Get the  values  -------------------------
         PRINT*, "nx: ", nx, "ny: ", ny
         ALLOCATE( lon_msk(nx) ) 
         ALLOCATE( lat_msk(ny) )
         ALLOCATE( msk(nx,ny) )

         STATUS = NF90_GET_VAR(MaskFileid, lon_mask_id, lon_msk)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_GET_VAR(MaskFileid, lat_mask_id, lat_msk)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

         STATUS = NF90_GET_VAR(MaskFileid, maskid, msk)
         IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)

!---  Close netCDF file  -----------------------

        STATUS = NF90_CLOSE(MaskFileid)
        IF (STATUS /= NF90_NOERR) CALL HANDLE_ERR(STATUS)
        PRINT*,"-----------------------------------------"
        PRINT*,'READ MASK FOR 150mx150m GRID'
        PRINT*,"-----------------------------------------"
      
     END SUBROUTINE READ_MASK150

!===========================================================================================
!   SUBROUTINE WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,ntime,&
!                           map_lon,map_lat,lon_cur,lat_cur,lon_wav,lat_wav,&
!                           surf_map,coast_map,displ_map,uhr_int,vhr_int,uwhr_int,vwhr_int,&
!                           uhr,vhr,uwhr,vwhr)

!===========================================================================================
     SUBROUTINE WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,nslick,ntime,&
                            map_lon,map_lat,lon_cur,lat_cur,cog_lon,cog_lat,time_vec,surf_map,&
                            displ_map,uhr,vhr,wdxhr,wdyhr,lon_wav,lat_wav,uwhr,vwhr,coast_map)

        IMPLICIT NONE
        INTEGER                                          :: cimx,cjmx,map_lon_n,map_lat_n,&
                                                            nseg,np_cst,ntime,i,nslick
        INTEGER                                          :: wimx,wjmx
        REAL,DIMENSION(map_lon_n+1)                      :: map_lon
        REAL,DIMENSION(map_lat_n+1)                      :: map_lat
        REAL,DIMENSION(cimx)                             :: lon_cur 
        REAL,DIMENSION(cjmx)                             :: lat_cur
        REAL,DIMENSION(wimx),OPTIONAL                    :: lon_wav
        REAL,DIMENSION(wjmx),OPTIONAL                    :: lat_wav
        REAL,DIMENSION(ntime,nslick)                     :: cog_lat, cog_lon
        REAL,DIMENSION(:),ALLOCATABLE                    :: pnt_vec
        REAL,DIMENSION(ntime)                            :: time_vec
        REAL,DIMENSION(map_lon_n+1,map_lat_n+1,ntime)    :: surf_map,displ_map
        REAL,DIMENSION(nseg,ntime),OPTIONAL              :: coast_map
        REAL,DIMENSION(cimx,cjmx,ntime)                  :: uhr,vhr
        REAL,DIMENSION(wimx,wjmx,ntime),OPTIONAL         :: uwhr,vwhr
        REAL,DIMENSION(nslick,ntime)                     :: wdxhr, wdyhr   

        print*, coast_map

        nx_oil = map_lon_n+1
        ny_oil = map_lat_n+1
        nx_cur = cimx 
        ny_cur = cjmx
        n_cog  = nslick
        ntim = ntime

        IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
           nx_wav = wimx
           ny_wav = wjmx
        ENDIF 

        IF ( PRESENT(coast_map) ) THEN

           np_cst = nseg
           ALLOCATE( pnt_vec(np_cst) )
           ALLOCATE( coast_out(np_cst) )
           pnt_vec = (/(i,i=1,np_cst)/)
        ENDIF

        ALLOCATE( surf_out(nx_oil,ny_oil) ) 
        ALLOCATE( displ_out(nx_oil,ny_oil) )
        
        ALLOCATE( u_cur_out(nx_cur,ny_cur) )
        ALLOCATE( v_cur_out(nx_cur,ny_cur) )

        IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
           ALLOCATE( u_wav_out(nx_wav,ny_wav) )
           ALLOCATE( v_wav_out(nx_wav,ny_wav) )
        ENDIF

        ALLOCATE( lon_cog_out(n_cog) )
        ALLOCATE( lat_cog_out(n_cog) )
        ALLOCATE( x_wnd_out(n_cog) )
        ALLOCATE( y_wnd_out(n_cog) )

! WRITING THE NetCDF OUTPUT

        print*, nc_file

        status=NF90_CREATE(nc_file,NF90_NOCLOBBER,ncid)
        if (status /= nf90_noerr) call handle_err(status)

! Define the dimensions. The time dimension is defined to have
! unlimited length - it can grow as needed. 

        status=NF90_DEF_DIM(ncid, X_OIL_NAME, nx_oil, lon_oil_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_DEF_DIM(ncid, Y_OIL_NAME, ny_oil, lat_oil_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        IF ( PRESENT(coast_map) ) THEN
           status=NF90_DEF_DIM(ncid, P_CST_NAME, np_cst, pnt_cst_dimid)
           if (status /= nf90_noerr) call handle_err(status)
        ENDIF

        status=NF90_DEF_DIM(ncid, X_CUR_NAME, nx_cur, lon_cur_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_DEF_DIM(ncid, Y_CUR_NAME, ny_cur, lat_cur_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
           status=NF90_DEF_DIM(ncid, X_WAV_NAME, nx_wav, lon_wav_dimid)
           if (status /= nf90_noerr) call handle_err(status)

           status=NF90_DEF_DIM(ncid, Y_WAV_NAME, ny_wav, lat_wav_dimid)
           if (status /= nf90_noerr) call handle_err(status)
        ENDIF

        status=NF90_DEF_DIM(ncid, N_COG_NAME, n_cog, cog_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_DEF_DIM(ncid, TIM_NAME, NF90_UNLIMITED, tim_dimid)
        if (status /= nf90_noerr) call handle_err(status)

        lon_oil_dims = (/ lon_oil_dimid /)
        lat_oil_dims = (/ lat_oil_dimid /)
        lon_cur_dims = (/ lon_cur_dimid /)
        lat_cur_dims = (/ lat_cur_dimid /)
        IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
           lon_wav_dims = (/ lon_wav_dimid /)
           lat_wav_dims = (/ lat_wav_dimid /)
        ENDIF
        tim_dims     = (/ tim_dimid /)

! Define the coordinate variables. 

     status=NF90_DEF_VAR(ncid, LON_OIL_NAME, NF90_REAL, lon_oil_dims, lon_oil_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, LAT_OIL_NAME, NF90_REAL, lat_oil_dims, lat_oil_varid)
     if (status /= nf90_noerr) call handle_err(status)


     status=NF90_DEF_VAR(ncid, LON_CUR_NAME, NF90_REAL, lon_cur_dims, lon_cur_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, LAT_CUR_NAME, NF90_REAL, lat_cur_dims, lat_cur_varid)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        status=NF90_DEF_VAR(ncid, LON_WAV_NAME, NF90_REAL, lon_wav_dims, lon_wav_varid)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_DEF_VAR(ncid, LAT_WAV_NAME, NF90_REAL, lat_wav_dims, lat_wav_varid)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_DEF_VAR(ncid, TIM_NAME, NF90_REAL, tim_dims, tim_varid)
     if (status /= nf90_noerr) call handle_err(status)

! Assign attributes to coordinate variables.

     status=NF90_PUT_ATT(ncid, lon_oil_varid, LONG_NAME, LONG_LON_OIL)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lat_oil_varid, LONG_NAME, LONG_LAT_OIL)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lon_oil_varid, UNITS, LON_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lat_oil_varid, UNITS, LAT_UNITS)
     if (status /= nf90_noerr) call handle_err(status)


     status=NF90_PUT_ATT(ncid, lon_cur_varid, LONG_NAME, LONG_LON_CUR)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lat_cur_varid, LONG_NAME, LONG_LAT_CUR)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lon_cur_varid, UNITS, LON_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lat_cur_varid, UNITS, LAT_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        status=NF90_PUT_ATT(ncid, lon_wav_varid, LONG_NAME, LONG_LON_WAV)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, lat_wav_varid, LONG_NAME, LONG_LAT_WAV)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, lon_wav_varid, UNITS, LON_UNITS)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, lat_wav_varid, UNITS, LAT_UNITS)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_PUT_ATT(ncid, tim_varid, LONG_NAME, LONG_TIM_NAME)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, tim_varid, UNITS, TIM_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.

     IF ( PRESENT(coast_map) ) THEN
        coast_dims  = (/ pnt_cst_dimid, tim_dimid /)
     ENDIF

     surf_dims   = (/ lon_oil_dimid, lat_oil_dimid, tim_dimid /)
     displ_dims  = (/ lon_oil_dimid, lat_oil_dimid, tim_dimid /)
    
     u_cur_dims  = (/ lon_cur_dimid, lat_cur_dimid, tim_dimid /)
     v_cur_dims  = (/ lon_cur_dimid, lat_cur_dimid, tim_dimid /)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        u_wav_dims  = (/ lon_wav_dimid, lat_wav_dimid, tim_dimid /)
        v_wav_dims  = (/ lon_wav_dimid, lat_wav_dimid, tim_dimid /)
     ENDIF

     lon_cog_dims = (/ cog_dimid, tim_dimid /)
     lat_cog_dims = (/ cog_dimid, tim_dimid /)
     x_wnd_dims   = (/ cog_dimid, tim_dimid /)
     y_wnd_dims   = (/ cog_dimid, tim_dimid /)

! Define the netCDF variables.

     IF ( PRESENT(coast_map) ) THEN
        status=NF90_DEF_VAR(ncid, COAST_NAME, NF90_REAL, coast_dims, coast_varid)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_DEF_VAR(ncid, SURF_NAME, NF90_REAL, surf_dims, surf_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, DISPL_NAME, NF90_REAL, displ_dims, displ_varid)
     if (status /= nf90_noerr) call handle_err(status)


     status=NF90_DEF_VAR(ncid, U_CUR_NAME, NF90_REAL, u_cur_dims, u_cur_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, V_CUR_NAME, NF90_REAL, v_cur_dims, v_cur_varid)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        status=NF90_DEF_VAR(ncid, U_WAV_NAME, NF90_REAL, u_wav_dims, u_wav_varid)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_DEF_VAR(ncid, V_WAV_NAME, NF90_REAL, v_wav_dims, v_wav_varid)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_DEF_VAR(ncid, LON_COG_NAME, NF90_REAL, lon_cog_dims, lon_cog_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, LAT_COG_NAME, NF90_REAL, lat_cog_dims, lat_cog_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, X_WND_NAME, NF90_REAL, x_wnd_dims, x_wnd_varid)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_DEF_VAR(ncid, Y_WND_NAME, NF90_REAL, y_wnd_dims, y_wnd_varid)
     if (status /= nf90_noerr) call handle_err(status)
      
! Assign attributes to the netCDF variables.

!------------------ LONG NAME ATTRIBUTE

     IF ( PRESENT(coast_map) ) THEN
        status=NF90_PUT_ATT(ncid, coast_varid, LONG_NAME, LONG_coast)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_PUT_ATT(ncid, surf_varid, LONG_NAME, LONG_surf)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, displ_varid, LONG_NAME, LONG_displ)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lon_cog_varid, LONG_NAME, LONG_LON_COG)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, lat_cog_varid, LONG_NAME, LONG_LAT_COG)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, x_wnd_varid, LONG_NAME, LONG_X_WND)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, y_wnd_varid, LONG_NAME, LONG_y_WND)
     if (status /= nf90_noerr) call handle_err(status)

!------------------ UNITS ATTRIBUTE

     status=NF90_PUT_ATT(ncid, surf_varid, UNITS, OIL_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( PRESENT(coast_map) ) THEN
        status=NF90_PUT_ATT(ncid, coast_varid, UNITS, CST_UNITS)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_PUT_ATT(ncid, displ_varid, UNITS, OIL_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, surf_varid, COOR, COOR_oil)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( PRESENT(coast_map) ) THEN
        status=NF90_PUT_ATT(ncid, coast_varid, COOR, COOR_cst)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

     status=NF90_PUT_ATT(ncid, displ_varid, COOR, COOR_oil)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, u_cur_varid, LONG_NAME, LONG_U_CUR)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, u_cur_varid, LONG_NAME, LONG_V_CUR)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, u_cur_varid, UNITS, VEL_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, v_cur_varid, UNITS, VEL_UNITS)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, u_cur_varid, COOR, COOR_cur)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_ATT(ncid, v_cur_varid, COOR, COOR_cur)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        status=NF90_PUT_ATT(ncid, u_wav_varid, LONG_NAME, LONG_U_WAV)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, u_wav_varid, LONG_NAME, LONG_V_WAV)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, u_wav_varid, UNITS, VEL_UNITS)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, v_wav_varid, UNITS, VEL_UNITS)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, u_wav_varid, COOR, COOR_wav)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_ATT(ncid, v_wav_varid, COOR, COOR_wav)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF
! End define mode.

     status=NF90_ENDDEF(ncid)
     if (status /= nf90_noerr) call handle_err(status)

! Write the coordinate variable data.

     status=NF90_PUT_VAR(ncid, lon_oil_varid, map_lon)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_VAR(ncid, lat_oil_varid, map_lat)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_VAR(ncid, lon_cur_varid, lon_cur)
     if (status /= nf90_noerr) call handle_err(status)

     status=NF90_PUT_VAR(ncid, lat_cur_varid, lat_cur)
     if (status /= nf90_noerr) call handle_err(status)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        status=NF90_PUT_VAR(ncid, lon_wav_varid, lon_wav)
        if (status /= nf90_noerr) call handle_err(status)

        status=NF90_PUT_VAR(ncid, lat_wav_varid, lat_wav)
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF

! These settings tell netcdf to write one timestep of data. (The
! setting of start(3) inside the loop below tells netCDF which
! timestep to write.)

     count_oil = (/ nx_oil, ny_oil, 1 /)
     start_oil = (/ 1, 1, 1 /)
     IF ( PRESENT(coast_map) ) THEN
        count_cst = (/ np_cst, 1/)
        start_cst = (/ 1, 1 /)
     ENDIF

     count_cur = (/ nx_cur, ny_cur, 1 /)
     start_cur = (/ 1, 1, 1 /)

     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        count_wav = (/ nx_wav, ny_wav, 1 /)
        start_wav = (/ 1, 1, 1 /)
     ENDIF

     count_cog = (/ n_cog, 1 /)
     start_cog = (/ 1, 1 /)

     count_tim = (/ 1 /)
     start_tim = (/ 1 /)

     DO rec = 1,ntim

          start_oil(3) = rec
          IF (nseg /= 0) start_cst(2) = rec
          start_cur(3) = rec
          IF ( wimx /= 0 .AND. wjmx /= 0 ) start_wav(3) = rec
          start_cog(2) = rec
          start_tim(1) = rec

          time_out = time_vec(rec)
 
          surf_out   = surf_map(:,:,rec)
          IF ( PRESENT(coast_map) ) coast_out  = coast_map(:,rec)
          displ_out  = displ_map(:,:,rec)
                 
          u_cur_out  = uhr(:,:,rec)
          v_cur_out  = vhr(:,:,rec)

          IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
             u_wav_out  = uwhr(:,:,rec)
             v_wav_out  = vwhr(:,:,rec)
          ENDIF

          lon_cog_out = cog_lon(rec,:)
          lat_cog_out = cog_lat(rec,:)
          x_wnd_out   = wdxhr(:,rec)
          y_wnd_out   = wdyhr(:,rec)

          status=NF90_PUT_VAR(ncid, tim_varid, time_out, start = start_tim)!, &
!                              count = count_tim)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, surf_varid, surf_out, start = start_oil, &
                              count = count_oil)
          if (status /= nf90_noerr) call handle_err(status)

          IF ( PRESENT(coast_map) ) THEN

             status=NF90_PUT_VAR(ncid, coast_varid, coast_out, start = start_cst, &
                                 count = count_cst)
             if (status /= nf90_noerr) call handle_err(status)
          ENDIF

          status=NF90_PUT_VAR(ncid, displ_varid, displ_out, start = start_oil, &
                              count = count_oil)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, u_cur_varid, u_cur_out, start = start_cur, &
                              count = count_cur)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, v_cur_varid, v_cur_out, start = start_cur, &
                              count = count_cur)
          if (status /= nf90_noerr) call handle_err(status)

          IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
             status=NF90_PUT_VAR(ncid, u_wav_varid, u_wav_out, start = start_wav, &
                                 count = count_wav)
             if (status /= nf90_noerr) call handle_err(status)

             status=NF90_PUT_VAR(ncid, v_wav_varid, v_wav_out, start = start_wav, &
                              count = count_wav)
             if (status /= nf90_noerr) call handle_err(status)
          ENDIF

          status=NF90_PUT_VAR(ncid, lon_cog_varid, lon_cog_out, start = start_cog, &
                              count = count_cog)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, lat_cog_varid, lat_cog_out, start = start_cog, &
                              count = count_cog)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, x_wnd_varid, x_wnd_out, start = start_cog, &
                              count = count_cog)
          if (status /= nf90_noerr) call handle_err(status)

          status=NF90_PUT_VAR(ncid, y_wnd_varid, y_wnd_out, start = start_cog, &
                              count = count_cog)
          if (status /= nf90_noerr) call handle_err(status)

     ENDDO
! -------------------------------------------------------------------------------------

     status=NF90_CLOSE(ncid)
     if (status /= nf90_noerr) call handle_err(status)
     print*, '   '
     print *,"***  SUCCESS WRITING ", nc_file, "!!!   ***"

! --------------------------------------------------------------------------------------

     DEALLOCATE( surf_out )
     IF ( PRESENT(coast_map) ) DEALLOCATE( coast_out )
     DEALLOCATE( displ_out )

     DEALLOCATE( u_cur_out )
     DEALLOCATE( v_cur_out )
     IF ( wimx /= 0 .AND. wjmx /= 0 ) THEN
        DEALLOCATE( u_wav_out )
        DEALLOCATE( v_wav_out )
     ENDIF
     END SUBROUTINE WRITE_OUTPUT
! =============================================================================================

END MODULE MODULE_NETCDF
