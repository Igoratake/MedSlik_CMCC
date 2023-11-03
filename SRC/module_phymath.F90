!-----------------------------------------------------------------------------------
!  module_phymath.F90
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
!  This routine module specifies physical variables and mathematical functions and 
!  subroutines needed to create the netcdf MEDSLIK_II output.  
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

MODULE MODULE_PHYMATH

      IMPLICIT NONE

! -------------------------------- GENERAL VARIABLES ------------------------------

      REAL,PARAMETER                                      :: g=9.81
      REAL,PARAMETER                                      :: PI=4.0*ATAN(1.0)      
      INTEGER                                             :: i, j, k, mc, nc, t, nst, err
      INTEGER                                             :: nfile, n
      CHARACTER(LEN=150)                                  :: EMPTY
      REAL                                                :: DUMMY
! ---------------------------------- INPUT VARIABLES ------------------------------      
! UNIT 11: out_txt
! UNIT 12: plgn
! UNIT 13: inp
! UNIT 14: par
! UNIT 22: nc_file

      CHARACTER(LEN=150)                                  :: out_dir,outname,oce_frq
      CHARACTER(LEN=*),PARAMETER                          :: inp = "medslik5.inp"
      CHARACTER(LEN=*),PARAMETER                          :: tmp = "medslik.tmp"
      CHARACTER(LEN=*),PARAMETER                          :: par = "medslik5.par"
      CHARACTER(LEN=150)                                  :: inp_file, par_file, tmp_file
      CHARACTER(LEN=150)                                  :: ocefile,metfile,wavfile
      CHARACTER(LEN=*),PARAMETER                          :: p_coord = "points_coord.txt"
      CHARACTER(LEN=*),PARAMETER                          :: map_cst = "coast.map"
      CHARACTER(LEN=*),PARAMETER                          :: plgs = "medslik_plgs.inp"
      CHARACTER(LEN=*),PARAMETER                          :: pnts = "medslik_pnts.inp"
      CHARACTER(LEN=150)                                  :: plg_file, pnt_file
      CHARACTER(LEN=150)                                  :: cur_dir, wav_dir, wnd_dir
      CHARACTER(LEN=3)                                    :: type_cu,type_wd,type_wv
      CHARACTER(LEN=4)                                    :: a0
      CHARACTER(LEN=150)                                   :: out_file
      CHARACTER(LEN=4)                                    :: a4
      INTEGER                                             :: nrec_cu,nrec_wd,nrec_wv

! ---------------------------------- READED VARIABLES ------------------------------


      INTEGER                                             :: icurrents, iwind, iwave, istoke
      INTEGER                                             :: kount, isat, sim_length
      INTEGER                                             :: numspills, ntime, ini_point 
      INTEGER                                             :: tot_point, sea_point, cst_point
      INTEGER,DIMENSION(3)                                :: npoint
      INTEGER,DIMENSION(:),ALLOCATABLE                    :: status_p

      REAL,DIMENSION(:),ALLOCATABLE                       :: lat_oil, lon_oil, conc_oil
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: lat_cst, lon_cst
      REAL                                                :: max_lat_oil, min_lat_oil
      REAL                                                :: max_lon_oil, min_lon_oil
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: cog_lat, cog_lon
     
! FROM input_file

      INTEGER                                             :: DD, MM, YR, HR ! day, month, year and hour of the spill
      INTEGER                                             :: tcomp
      INTEGER                                             :: cur_num, wnd_num, wav_num
      INTEGER                                             :: delt_out
      CHARACTER(LEN=8),DIMENSION(:),ALLOCATABLE           :: cur_file
      CHARACTER(LEN=8),DIMENSION(:),ALLOCATABLE           :: wnd_file
      CHARACTER(LEN=8),DIMENSION(:),ALLOCATABLE           :: wav_file

! FROM par_file

      CHARACTER(LEN=150),DIMENSION(:),ALLOCATABLE         :: cur_name
      CHARACTER(LEN=150),DIMENSION(:),ALLOCATABLE         :: wav_name
      CHARACTER(LEN=150),DIMENSION(:),ALLOCATABLE         :: wnd_name      
      INTEGER                                             :: iwav
      CHARACTER(LEN=2)                                    :: ihr_str, a2
            
! FROM cur_name

      INTEGER                                             :: cimx, cjmx, cndata
      INTEGER                                             :: wimx, wjmx, wndata
      INTEGER                                             :: wdimx, wdjmx, wdndata
      REAL,DIMENSION(:),ALLOCATABLE                       :: raw_lon, raw_lat
      REAL                                                :: clon1, clon2, clat1, clat2, cdlon, cdlat
      REAL                                                :: wlon1, wlon2, wlat1, wlat2, wdlon, wdlat
      REAL                                                :: wdlon1, wdlon2, wdlat1, wdlat2, wddlon, wddlat
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: us, vs, dus, dvs
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: uw, vw, duw, dvw
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: wndx, wndy, dwndx, dwndy
      REAL,DIMENSION(:,:,:),ALLOCATABLE                   :: uhr, vhr, uwhr, vwhr
      REAL,DIMENSION(:,:,:),ALLOCATABLE                   :: wndxhr, wndyhr

! -------------------------------- CALCULATED VARIABLES -------------------------------

      INTEGER                                             :: idd, imm, iyr, ihr, mxst
      INTEGER                                             :: id, im, iy, thr
      INTEGER                                             :: out_len, nday
      INTEGER                                             :: cfile, wfile, wdfile,wndrec
      REAL,DIMENSION(:),ALLOCATABLE                       :: time_vec
      REAL                                                :: time, sim_start, nstph, timehr
      REAL                                                :: fac_c, fac_w, frac_c, frac_w, fac_wd, frac_wd
      REAL                                                :: sim_time, dt_cur, dt_wav, dt_wnd
      REAL,DIMENSION(:),ALLOCATABLE                       :: lon_cur, lat_cur
      REAL,DIMENSION(:),ALLOCATABLE                       :: lon_wav, lat_wav
      REAL,DIMENSION(:),ALLOCATABLE                       :: lon_wnd, lat_wnd
      REAL,DIMENSION(:),ALLOCATABLE                       :: cur_time, wav_time
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: wnd_time
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: mask_c, mask_w
      REAL,DIMENSION(:,:,:),ALLOCATABLE                   :: usseol, vsseol, uwseol, vwseol
      REAL,DIMENSION(:,:,:),ALLOCATABLE                   :: uhr_int, vhr_int, uwhr_int, vwhr_int
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: wndx_int, wndy_int 
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: map_lon2, map_lat2
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: lon_cur2,lat_cur2,lon_wav2,lat_wav2
      LOGICAL                                             :: readdata      

! REGULAR OUTPUT MAP VARIABLES DECLARATION

      REAL,PARAMETER                                      :: lat_margin=8., lon_margin=10.        ! boundaries margins of the regular
                                                                                                  ! output map in minutes
      REAL,PARAMETER                                      :: lat_cell=3., lon_cell=4.
      REAL,PARAMETER                                      :: pixel=150.                           ! grid resolution of the regular 
                                                                                                  ! output map in meters
      REAL                                                :: d, ref_lat                           ! round factor
      REAL                                                :: max_lat, min_lat, max_lon, min_lon
      REAL                                                :: min_map_lat,min_map_lon
      REAL                                                :: max_map_lat,max_map_lon
      REAL                                                :: map_dy, map_dx                       ! grid resolution of the regular 
                                                                                                  ! output map in minutes
      REAL,PARAMETER                                      :: mile=1852.                           ! nautical mile in meters
      REAL,PARAMETER                                      :: deg2met = 60. * mile                 ! factor to convert distance in 
                                                                                                  ! meter to distance in degree
      REAL,PARAMETER                                      :: deg2rad = PI / 180. 
      INTEGER                                             :: map_lat_n, map_lon_n, nseg
      REAL,DIMENSION(:),ALLOCATABLE                       :: map_lat, map_lon
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: lon_sort,lat_sort
      INTEGER,DIMENSION(:),ALLOCATABLE                    :: oil_lat_indx, oil_lon_indx
      INTEGER,DIMENSION(:,:),ALLOCATABLE                  :: oil_indx
      REAL,DIMENSION(:,:,:),ALLOCATABLE                   :: surf_map, displ_map
      REAL,DIMENSION(:,:),ALLOCATABLE                     :: cst_seg, coast_map
      REAL                                                :: distance
!      REAL,PARAMETER                                      :: pixel_lat = pixel/deg2met
!      REAL                                                :: pixel_lon, eps
!      REAL,PARAMETER                                      :: eps = (sqrt(2*(pixel**2)) / 2.) / deg2met
      REAL,PARAMETER                                      :: eps = pixel / deg2met
      INTEGER                                             :: lon_start, lon_end, lat_start, lat_end
      INTEGER,DIMENSION(:,:),ALLOCATABLE                  :: mask_medslik        
      
CONTAINS

!=================================================================================
!                            FUNCTION  FindMin
!---------------------------------------------------------------------------------
!
!    This function returns the array index of the minimum in the section
!    between Start and End of the array x.
!
!=================================================================================

    FUNCTION  FindMin(x, first, last)

         IMPLICIT  NONE
         REAL, DIMENSION(1:), INTENT(IN)    :: x
         INTEGER, INTENT(IN)                :: first, last
         REAL                               :: Minimum
         INTEGER                            :: Location
         INTEGER                            :: i
         INTEGER                            :: FindMin

         Minimum  = x(first)          ! assume the first is the min
         Location = first             ! record its position

         DO i = first+1,last          ! start with next elements

            IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
               Minimum  = x(i)        !      Yes, a new minimum found
               Location = i           !      record its position
            END IF

         END DO

         FindMin = Location            ! return the index position

    END FUNCTION  FindMin

!=================================================================================
!                             SUBROUTINE  Swap():
!
!    This subroutine swaps the values of its two formal arguments.
!---------------------------------------------------------------------------------

    SUBROUTINE  Swap(a, b)

         IMPLICIT  NONE
         REAL, INTENT(INOUT) :: a, b
         REAL                :: Temp

         Temp = a
         a    = b
         b    = Temp

    END SUBROUTINE  Swap

!=====================================================================================
!                                  FUNCTION NEAR
!-------------------------------------------------------------------------------------
!
!   NEAR finds the index of x element that is closest to the point x0.
!   x is an array of dimension n, x0 is a point.  
!   Distance is the abs(x-x0).
!
!=====================================================================================

    FUNCTION Near(n,x,x0)

         INTEGER,INTENT(IN)              :: n
         REAL,DIMENSION(n),INTENT(IN)    :: x
         REAL,INTENT(IN)                 :: x0
         INTEGER                         :: Near
         REAL,DIMENSION(n)               :: distance

         distance = ABS(x-x0)
         Near = FindMin(distance,1,n)

    END FUNCTION Near

!=====================================================================================
!                                 SUBROUTINE NEAR2D
!-------------------------------------------------------------------------------------
!
!   NEAR2D finds the index of the (k)th element of the matrix xy (n,2) that is 
!   closest to the point (x0,y0).
!   x and y are arrays of dimension n, (x0,y0) is a point.  
!   Distance is the sqrt((x-x0)^2+(y-y0)^2).
!
!=====================================================================================

    SUBROUTINE Near2D(n,x,y,x0,y0,indx,dist_nearest)

         INTEGER,INTENT(IN)              :: n
         REAL,DIMENSION(n),INTENT(IN)    :: x,y
         REAL,INTENT(IN)                 :: x0,y0
         INTEGER                         :: indx
         REAL,DIMENSION(n)               :: dist
         REAL                            :: dist_nearest
         

         dist = SQRT((x-x0)**2 + (y-y0)**2)
         indx = FindMin(dist,1,n)
         dist_nearest = dist(indx)

    END SUBROUTINE Near2D

!=================================================================================
!                           FUNCTION  Sort(x):
!
!    This FUNCTION receives a 1D array x(:) and sorts it into ascending order.
!---------------------------------------------------------------------------------

    FUNCTION  Sort(x)

         IMPLICIT  NONE
         REAL,DIMENSION(1:),INTENT(IN)         :: x
         INTEGER                               :: x_size
         INTEGER                               :: i
         INTEGER                               :: Location
         REAL,DIMENSION(SIZE(x))               :: Sort

         x_size = SIZE(x)
         Sort = x

         DO i = 1, x_size-1                          ! except for the last

            Location = FindMin(Sort,i,x_size)           ! find min from this to last
            CALL  Swap(Sort(i), Sort(Location))      ! swap this and the minimum

         END DO

    END FUNCTION  Sort

!=================================================================================
!                           FUNCTION  m_lldist(lon1,lat1,lon2,lat2):
!
!  This function calculates the distance in meters between point1 (lon1,lat1) 
!  and point2 (lon2,lat2) using the Haversine formula on a spherical earth of 
!  radius 6378.137km.
!---------------------------------------------------------------------------------

    FUNCTION  m_lldist(lon1,lat1,lon2,lat2)

         IMPLICIT NONE
         REAL,PARAMETER                        :: PI180 = PI / 180.
         REAL,PARAMETER                        :: EARTH_RADIUS = 6378.137
         REAL                                  :: lon1,lon2,lat1,lat2
         REAL                                  :: dlon, dlat, a, angle, dist
         REAL                                  :: m_lldist

         lon1 = lon1 * PI180
         lon2 = lon2 * PI180
         lat1 = lat1 * PI180
         lat2 = lat2 * PI180

         dlon = lon2 - lon1 
         dlat = lat2 - lat1 

         a = (SIN(dlat/2.))**2 + COS(lat1) * COS(lat2) * (SIN(dlon/2.))**2

         angle = 2. * ATAN2(SQRT(a),SQRT(1. - a))

         dist = EARTH_RADIUS * angle

         m_lldist = 1000. * dist

    END FUNCTION m_lldist

!===================================================================================
!                                FUNCTION JULDAY 
!                    computes the julian day for a given date
!===================================================================================

    FUNCTION JULDAY(IDD,IMM,IYR)
      
         IMPLICIT NONE
         INTEGER                     :: IDD,IMM,IYR
         INTEGER,DIMENSION(12)       :: JS,JE
         INTEGER                     :: JULDAY,LEAP_YR,k,i

         LEAP_YR = 0
         IF ( ( ((IYR/4)*4.eq.IYR).and.((IYR/100)*100.ne.IYR) ).or.&
            ( (IYR/400)*400.eq.IYR) ) LEAP_YR = 1

         JS(1) = 1
         JE(1) = 31
         JS(2) = 32
         JE(2) = 59 + LEAP_YR

         k=1

         DO i=3,12

            JS(i) = JE(i-1) + 1
            JE(i) = JS(i) + 29 + k
            k = 1 - k
            IF (i.eq.7) k=1

         ENDDO

         JULDAY = JS(IMM) + IDD - 1

         RETURN

    END FUNCTION JULDAY


!====================================================================================
!                                 FUNCTION JDIFF 
!              gives the no of days from id1/im1/iy1 to id2/im2/iy2 
!           the day-difference may be negative and the years differ by 1
!====================================================================================
 
    FUNCTION  JDIFF(ID1,IM1,IY1,ID2,IM2,IY2)

         IMPLICIT NONE
         INTEGER            :: IY1,IM1,ID1,IY2,IM2,ID2
         INTEGER            :: JDIFF

         IF (IY2.eq.IY1)   JDIFF = JULDAY(ID2,IM2,IY2) - JULDAY(ID1,IM1,IY1)

         IF (IY1.eq.IY1+1) JDIFF = JULDAY(ID2,IM2,IY2) + JULDAY(31,12,IY1)&
                                   - JULDAY(ID1,IM1,IY1)
         IF (IY2.eq.IY1-1) JDIFF = JULDAY(ID2,IM2,IY2) - JULDAY(31,12,IY2)&
                                   - JULDAY(ID1,IM1,IY1)

         RETURN

     END FUNCTION JDIFF

!=====================================================================================
!                              FUNCTION READFCST_1HR 
!                         read forecast hourly data fields 
!=====================================================================================

     FUNCTION READFCST_1HR(dir,file_name,lon,lat,field_type,field_comp,ndata,mmax,nmax)

         IMPLICIT NONE

         CHARACTER(LEN=150)                    :: dir
         CHARACTER(LEN=150)                    :: file_name
         CHARACTER(LEN=1)                      :: field_type,field_comp
         REAL,DIMENSION(mmax)                  :: lon
         REAL,DIMENSION(nmax)                  :: lat
         REAL,DIMENSION(:,:),ALLOCATABLE       :: raw
         REAL,DIMENSION(mmax,nmax)             :: field, readfcst_1hr
         INTEGER                               :: mmax, nmax, ndata, i, j, k

         ALLOCATE( raw(ndata,4) )
!         PRINT*, FILE_NAME       
         OPEN(UNIT=71,FILE=TRIM(dir)//TRIM(file_name))
         DO i=1,5
            READ(71,*) EMPTY
         ENDDO
         DO i=1,ndata
            IF ( field_type == "c" ) THEN
               READ(71,*) raw(i,1), raw(i,2), EMPTY, raw(i,3), raw(i,4)
            ENDIF
            IF ( field_type == "w" ) THEN
               READ(71,*) raw(i,1), raw(i,2), raw(i,3), raw(i,4), EMPTY 
            ENDIF 
!            PRINT*, RAW(i,:)
         ENDDO

         field = 9999.
         DO i=1,mmax
            DO j=1,nmax
               DO k=1,ndata        
!                  PRINT*, RAW(k,2), LON(i), RAW(k,1), LAT(j)                 
                  IF ( (lon(i) == raw(k,2)).and.(lat(j) == raw(k,1)) ) THEN

                     IF ( field_comp.eq."u" ) THEN
                        field(i,j) = raw(k,3)
                     ENDIF
                     IF ( field_comp.eq."v" ) THEN
                        field(i,j) = raw(k,4)
                     ENDIF

                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CLOSE(71) 

         readfcst_1hr = field
         DEALLOCATE(RAW) 
         
     END FUNCTION READFCST_1HR    

!=====================================================================================
!                              FUNCTION READFCST_NEW 
!                         read forecast hourly data fields 
!=====================================================================================

     FUNCTION READFCST_NEW(dir,file_name,lon,lat,field_type,field_comp,ndata,mmax,nmax)

         IMPLICIT NONE

         CHARACTER(LEN=150)                    :: dir
         CHARACTER(LEN=150)                    :: file_name
         CHARACTER(LEN=1)                      :: field_type,field_comp
         REAL,DIMENSION(mmax)                  :: lon
         REAL,DIMENSION(nmax)                  :: lat
         REAL,DIMENSION(:,:),ALLOCATABLE       :: raw
         REAL,DIMENSION(mmax,nmax)             :: field, readfcst_new
         INTEGER                               :: mmax, nmax, ndata, i, m, n
         REAL                                  :: dlon, dlat

         ALLOCATE( raw(ndata,4) )
!         PRINT*, FILE_NAME       
         OPEN(UNIT=71,FILE=TRIM(dir)//TRIM(file_name))
         DO i=1,5
            READ(71,*) EMPTY
         ENDDO

         dlon = lon(2) - lon(1) 
         dlat = lat(2) - lat(1)
         field = 9999.

         DO i=1,ndata
            IF ( field_type == "c" ) THEN
               READ(71,*) raw(i,1), raw(i,2), EMPTY, raw(i,3), raw(i,4)
            ENDIF
            IF ( field_type == "w" ) THEN
               READ(71,*) raw(i,1), raw(i,2), raw(i,3), raw(i,4), EMPTY
            ENDIF
            
            m = INT((raw(i,2) - lon(1))/dlon+1)
            n = INT((raw(i,1) - lat(1))/dlat+1)
            
            IF ((m > mmax).or.(n > nmax).or.(m < 1).or.(n < 1)) CYCLE
            IF ( field_comp.eq."u" ) THEN
               field(m,n) = raw(i,3)
            ENDIF
            IF ( field_comp.eq."v" ) THEN
               field(m,n) = raw(i,4)
            ENDIF            
         ENDDO
 
         CLOSE(71)

         readfcst_new = field
         DEALLOCATE(RAW)

     END FUNCTION READFCST_NEW

!=====================================================================================
!                              FUNCTION READMET 
!                         read meteo forecast data fields 
!=====================================================================================

     FUNCTION READMET(dir,file_name,wndrec,lon,lat,field_comp,ndata,mmax,nmax)

         IMPLICIT NONE

         CHARACTER(LEN=150)                    :: dir
         CHARACTER(LEN=14)                     :: file_name
         CHARACTER(LEN=1)                      :: field_comp
         CHARACTER(LEN=2)                      :: nout
         CHARACTER(LEN=11)                     :: formt
         REAL,DIMENSION(mmax)                  :: lon
         REAL,DIMENSION(nmax)                  :: lat
         REAL,DIMENSION(:,:),ALLOCATABLE       :: raw
         REAL,DIMENSION(mmax,nmax)             :: field, readmet
         INTEGER                               :: mmax, nmax, ndata, i, j, k
         INTEGER                               :: wndrec

         ALLOCATE( raw(ndata,2*(nrec_wd+1)) )
         OPEN(UNIT=71,FILE=TRIM(dir)//TRIM(file_name))

         DO i=1,5
            READ(71,*) EMPTY
         ENDDO

         WRITE(nout,'(i2)') 2*nrec_wd+2
         formt = "("//nout//"(F12.8)"//")"
         
         DO i=1,ndata
            READ(71,formt) raw(i,1),raw(i,2),(raw(i,k), k=3,2*nrec_wd+2)
         ENDDO

         field = 9999.

         DO i=1,mmax
            DO j=1,nmax
               DO k=1,ndata
                  IF ( (lon(i) == raw(k,2)).and.(lat(j) == raw(k,1)) ) THEN

                     IF ( field_comp.eq."x" ) THEN
                        field(i,j) = raw(k,2*wndrec+1)
                     ENDIF
                     IF ( field_comp.eq."y" ) THEN
                        field(i,j) = raw(k,2*wndrec+2)
                     ENDIF

                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         CLOSE(71)

         readmet = field

         DEALLOCATE( raw )

     END FUNCTION READMET

!=====================================================================================
!                           SUBROUTINE OIL_GRID(lon_msk,lat_msk)
!=====================================================================================
     SUBROUTINE OIL_GRID(msk,lon_msk,lat_msk,nx,ny)

         IMPLICIT NONE
         INTEGER                    :: nx,ny
         REAL,DIMENSION(nx)         :: lon_msk
         REAL,DIMENSION(ny)         :: lat_msk
         REAL,DIMENSION(nx,ny)      :: msk
         REAL                       :: conc, dens
         REAL,PARAMETER             :: area_grid = pixel**2
         
         
         
         ntime = mxst + 1
         ALLOCATE ( time_vec(ntime) )
         time_vec(1) = 0.

         ALLOCATE ( cog_lon(ntime,numspills) )
         ALLOCATE ( cog_lat(ntime,numspills) )

         OPEN(UNIT=21,FILE=p_coord)

         DO t=1,ntime

            IF ( t == 1 ) THEN
               sea_point = ini_point
               cst_point = 0

               ALLOCATE(conc_oil(sea_point))
               ALLOCATE(status_p(sea_point))

               conc_oil = 1000.
               status_p = 5

               WRITE(21,*) "TIME    ", "NPOINTS     "
               WRITE(21,*) time_vec(t), sea_point
               WRITE(21,*) "lat   ", "lon   ", "status   ", "conc   "
               DO k=1,sea_point
                  WRITE(21,'(2(F10.5),I4,3X,F14.5)') lat_oil(k), lon_oil(k), status_p(k), conc_oil(k)
               END DO

               min_lat = MINVAL(lat_oil)
               max_lat = MAXVAL(lat_oil)
               min_lon = MINVAL(lon_oil)
               max_lon = MAXVAL(lon_oil)

               DO i=1,numspills
                  cog_lat(t,i) = (min_lat + max_lat)/2.
                  cog_lon(t,i) = (min_lon + max_lon)/2.
               ENDDO

               DEALLOCATE(lat_oil) 
               DEALLOCATE(lon_oil)
               DEALLOCATE(conc_oil)
               DEALLOCATE(status_p)
               
            ELSE

               npoint = 0
               DO n=1,3
               
                  IF (n==1) a0 = ".srf"
                  IF (n==2) a0 = ".cst"
                  IF (n==3) a0 = ".dsp"

                  WRITE(a4,'(i4)') t-1
 
                  IF ( t < 11 ) THEN
                     a4 = "000"//a4(4:4)
                  ELSE IF (( t >= 11 ).AND.( t < 101 )) THEN
                     a4 = "00"//a4(3:4)
                  ELSE IF ( t >= 101 ) THEN
                     a4 = "0"//a4(2:4)
                  ENDIF
                  out_file = TRIM(out_dir)//TRIM(outname)//a4//a0
                  nfile=80+n

                  OPEN (nfile,file=TRIM(out_file))
                  IF ( n==1 ) THEN
                     READ(nfile,*) EMPTY
                     READ(nfile,*) time_vec(t) 
                     DO i=1,3
                         READ(nfile,*) EMPTY
                     ENDDO
                  ELSE
                     DO i=1,5
                         READ(nfile,*) EMPTY
                     ENDDO
                  ENDIF
 
                  DO i=1,numspills
                     IF ( n==1 ) THEN     
                        READ(nfile,*) empty, cog_lat(t,i), cog_lon(t,i) ! baricentro line
                     ELSE
                        READ(nfile,*) EMPTY
                     ENDIF
                  ENDDO
                                        
                  READ(nfile,*) EMPTY
                  READ(nfile,*) EMPTY
                  READ(nfile,*) npoint(n)
                  CLOSE(nfile)
               ENDDO

               tot_point = npoint(1)+npoint(2)+npoint(3)
               sea_point = npoint(1)+npoint(3)
               cst_point = npoint(2)

               WRITE(21,*) "TIME    ", "NPOINTS     "
               WRITE(21,*) time_vec(t), sea_point
               WRITE(21,*) "lat   ", "lon   ", "status   ", "conc   "
             
               ALLOCATE(lon_oil(tot_point))
               ALLOCATE(lat_oil(tot_point))
               ALLOCATE(lat_cst(cst_point,2))
               ALLOCATE(lon_cst(cst_point,2))
               ALLOCATE(conc_oil(tot_point))
               ALLOCATE(status_p(tot_point))
               
               kount = 0

               DO n=1,3
                  IF (n==1) a0 = ".srf"
                  IF (n==2) a0 = ".cst"
                  IF (n==3) a0 = ".dsp"

                  IF ( t < 11 ) THEN
                     a4 = "000"//a4(4:4)
                  ELSE IF (( t >= 11 ).AND.( t < 101 )) THEN
                     a4 = "00"//a4(3:4)
                  ELSE IF ( t >= 101 ) THEN
                     a4 = "0"//a4(2:4)
                  ENDIF
                  out_file = TRIM(out_dir)//TRIM(outname)//a4//a0

                  nfile=80+n

                  OPEN(nfile,file=out_file)

                  DO i=1,5
                     READ(nfile,*) EMPTY
                  ENDDO
                  DO i=1,numspills
                     READ(nfile,*) EMPTY
                  ENDDO
                  DO i=1,4
                     READ(nfile,*) EMPTY
                  ENDDO
                  
                  DO i=1,npoint(n)
                     
                     conc = 0.
                     dens = 0.
                     k = kount+i
                     
                     IF (n /= 2) THEN
                        READ(nfile,*) lat_oil(k), lon_oil(k), conc_oil(k)
                        IF (n == 1) status_p(k) = 5
                        IF (n == 3) status_p(k) = 1
                        WRITE(21,'(2(F10.5),I4,3X,F14.5)') lat_oil(k), lon_oil(k),&
                                                        status_p(k), conc_oil(k)
                     ELSE
                        IF (npoint(n) /= 0) THEN
                           READ(nfile,"(5(F11.7))") lat_cst(i,1), lon_cst(i,1), lat_cst(i,2), &
                                         lon_cst(i,2), conc_oil(k)
                           lat_oil(k) = (lat_cst(i,1)+lat_cst(i,2))/2.
                           lon_oil(k) = (lon_cst(i,1)+lon_cst(i,2))/2.
                        ENDIF
                     ENDIF
                  ENDDO
                  CLOSE(nfile)
                  kount = kount + npoint(n)
               ENDDO 

               IF ( MINVAL(lat_oil) < min_lat ) THEN
                  min_lat = MINVAL(lat_oil)
               END IF
               IF ( MAXVAL(lat_oil) > max_lat ) THEN
                  max_lat = MAXVAL(lat_oil)
               END IF
               IF ( MINVAL(lon_oil) < min_lon ) THEN
                  min_lon = MINVAL(lon_oil)
               END IF
               IF ( MAXVAL(lon_oil) > max_lon ) THEN
                  max_lon = MAXVAL(lon_oil)
               END IF
               
               do i=1,tot_point
                  if (lon_oil(i) > max_lon) then
                     max_lon = lon_oil(i)
                  endif
               enddo 

               DEALLOCATE(lat_oil)
               DEALLOCATE(lon_oil)
               DEALLOCATE(conc_oil)
               DEALLOCATE(lat_cst)
               DEALLOCATE(lon_cst)
               DEALLOCATE(status_p)

            END IF
         END DO

         PRINT*, "MIN LAT: ", min_lat, "MAX LAT: ", max_lat
         PRINT*, "MIN LON: ", min_lon, "MAX LON: ", max_lon

         CLOSE(21)

         min_map_lon = lon_msk(Near(nx,lon_msk,min_lon))
         max_map_lon = lon_msk(Near(nx,lon_msk,max_lon))
         min_map_lat = lat_msk(Near(ny,lat_msk,min_lat))
         max_map_lat = lat_msk(Near(ny,lat_msk,max_lat))

         PRINT*, "MIN LAT MAP: ", min_map_lat, "MAX LAT MAP: ", max_map_lat
         PRINT*, "MIN LON MAP: ", min_map_lon, "MAX LON MAP: ", max_map_lon

         DO i=1,nx

            IF (lon_msk(i) == min_map_lon ) THEN
                lon_start = i
                PRINT*, "lon_start: ", lon_start
            ENDIF
            IF (lon_msk(i) == max_map_lon ) THEN
               lon_end = i
               PRINT*, "lon_end: ", lon_end
            ENDIF

         ENDDO
         DO i=1,ny

            IF (lat_msk(i) == min_map_lat ) THEN
               lat_start = i
               PRINT*, "lat_start: ", lat_start
            ENDIF
            IF (lat_msk(i) == max_map_lat ) THEN
               lat_end = i
               PRINT*, "lat_end: ", lat_end
            ENDIF

         ENDDO

         map_lon_n = lon_end-lon_start
         map_lat_n = lat_end-lat_start

         ALLOCATE( mask_medslik(map_lon_n+1,map_lat_n+1) )
         ALLOCATE( map_lat(map_lat_n+1) )
         ALLOCATE( map_lon(map_lon_n+1) )

         mask_medslik = msk(lon_start:lon_end,lat_start:lat_end)
         
         DO i=0,map_lat_n
            map_lat(i+1) = lat_msk(lat_start+i)
         ENDDO
         DO i=0,map_lon_n
            map_lon(i+1) = lon_msk(lon_start+i)
         ENDDO

     END SUBROUTINE OIL_GRID

!=====================================================================================
!                           SUBROUTINE COAST_SEG(mnlon,mxlon,mnlat,mxlat)
!=====================================================================================
     SUBROUTINE COAST_SEG(mnlon,mxlon,mnlat,mxlat)
 
         IMPLICIT NONE
         REAL                               :: mnlon, mxlon
         REAL                               :: mnlat, mxlat
         CHARACTER(LEN=*),PARAMETER         :: cstfile = "data/plot_medf.map"
         INTEGER                            :: var
         INTEGER                            :: totseg, tyseg
         REAL                               :: loncst, latcst, lon1, lat1
         REAL                               :: lon2, lat2
         LOGICAL                            :: presence, wrte

         OPEN(UNIT=57,FILE=cstfile)
         OPEN(UNIT=56,FILE=map_cst)
         var = 0
         nseg = 0
         
         DO WHILE (var == 0 )

            READ(UNIT=57,FMT=*,IOSTAT=var) totseg, tyseg

            presence = .FALSE.
            wrte = .FALSE.           
 
            DO k=1,totseg
               READ(UNIT=57,FMT=*,IOSTAT=var) loncst, latcst
               IF ( loncst>=mnlon.AND.loncst<=mxlon.AND.latcst>=mnlat.AND.latcst<=mxlat) THEN
                  presence = .TRUE.
                  IF (wrte .EQV. .FALSE.) THEN
                      wrte = .TRUE.
                      lon1 = loncst
                      lat1 = latcst
                      CYCLE
                  ENDIF
                  IF (wrte .EQV. .TRUE.) THEN
                      wrte = .FALSE.
                      lon2 = loncst
                      lat2 = latcst 
                      WRITE(56,"(4(F10.5))") lon1, lat1, lon2, lat2
                  ENDIF
                  nseg = nseg + 1
               ENDIF
            ENDDO

            IF (presence .EQV. .TRUE.) WRITE(56,*) "999999"

         ENDDO

         PRINT*, "NUMBER OF COASTAL SEGMENTS INTO THE DOMAIN: ", nseg

         CLOSE(57)
         CLOSE(56)

         IF ( nseg /= 0 ) THEN

            ALLOCATE(cst_seg(nseg,4))

            OPEN(UNIT=56,FILE=map_cst)

            var = 0
            i=1

            DO WHILE (var == 0 .AND. i <= nseg)

               READ(UNIT=56,FMT=*,IOSTAT=var) lon1, lat1, lon2, lat2 

               IF (loncst == 999999) CYCLE
            
               cst_seg(i,1) = lon1     ! lon point 1
               cst_seg(i,2) = lat1     ! lat point 1
               cst_seg(i,3) = lon2     ! lon point 2
               cst_seg(i,4) = lat2     ! lat point 2
 
               i = i + 1     
         
            ENDDO

            CLOSE(56)

         ELSE
 
            ALLOCATE(cst_seg(1,4))
            cst_seg = 0.

         ENDIF

     END SUBROUTINE COAST_SEG


END MODULE MODULE_PHYMATH
