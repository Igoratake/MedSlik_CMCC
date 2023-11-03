!================================================================================
!  MEDSLIK-II_1.02                                                               |
!                                                                                |
!  Oil spill fate and transport numerical model                                  |
!--------------------------------------------------------------------------------|
!  Extract_II.F90                                                                |
!                                                                                |
!  This routine reads winds and currents from                                    |
!  meteo-oceanogrpahic model output (NetCDF files)                               |
!  and write them into MEDSLIK-II ascii formatted input files                    |
!                                                                                |
!--------------------------------------------------------------------------------|
!                                                                                |
!  Copyright (C) <2012>                                                          |
!                                                                                |
!  This program was originally written by Robin Lardner and George Zodiatis.     |
!                                                                                |
!  Subsequent additions and modifications have been made by Michela De Dominicis |
!  and Diego Bruciaferri.                                                        |
!                                                                                |
!--------------------------------------------------------------------------------|
!  The development of the MEDSLIK-II model is supported by a formal agreement    |
!  Memorandum of Agreement for the Operation and Continued Development of        |
!  MEDSLIK-II signed by the following institutions:                              |
!                                                                                |
!  INGV     - Istituto Nazionale di Geofisica e Vulcanologia                     |
!  OC-UCY   - Oceanography Center at the University of Cyprus                    |
!  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per                  |
!             lo Studio dell’Ambiente Marino Costiero                            |
!  CMCC     - Centro Euro-Mediterraneo sui Cambiamenti Climatici                 |
!--------------------------------------------------------------------------------|  
!  This program is free software: you can redistribute it and/or modify          |
!  it under the terms of the GNU General Public License as published by          |
!  the Free Software Foundation, either version 3 of the License, or             |
!  any later version.                                                            |
!                                                                                |
!  This program is distributed in the hope that it will be useful,               |
!  but WITHOUT ANY WARRANTY; without even the implied warranty of                |
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 |
!  GNU General Public License for more details.                                  |
!  You should have received a copy of the GNU General Public License             |
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.         |
!================================================================================

 MODULE env

         IMPLICIT NONE

         CHARACTER        :: regn*4, indate_cu(30)*8, indate_wd(30)*8, indate_wv(30)*8
         CHARACTER        :: fcdir*120, type_cu*2, type_wd*2, type_wv*2
         REAL             :: alon1, alon2, alat1, alat2
         REAL             :: blon1, blon2, blat1, blat2, blon, blat
         REAL             :: clon1, clon2, clat1, clat2
         INTEGER          :: numfiles_cu, numfiles_wd, numfiles_wv, iviod
         INTEGER          :: numf_cu
         INTEGER          :: icurrent, iwind, iwave
         REAL             :: oplon0, oplat0, op_dlon, op_dlat
         INTEGER          :: i_first, i_last, j_first, j_last
         INTEGER          :: imax, imax1, jmax, jmax1
         INTEGER          :: imx, jmx, kmx, ktmx
         INTEGER          :: i1, i2, j1, j2
         INTEGER          :: kount, nore
         CHARACTER        :: prdate*50, outfile*50, infile_WI*120, infile_WV*120
         CHARACTER        :: infile_T*120, infile_U*120, infile_V*120
         CHARACTER        :: heads*150, empty*80, ora*2, ore*2, mm*2, m*2
         LOGICAL          :: ex

 END MODULE env

!-------------------------------------------------------------------------------------

 MODULE utils

         IMPLICIT NONE

         INTEGER          :: i, j, k, n, t, tt
         INTEGER          :: nx, ny, nz, nt, km
         REAL             :: udef, udef2, rhoa, miss
         REAL,PARAMETER   :: g=9.81
         REAL,PARAMETER   :: PI=4.0*ATAN(1.0)
         
 CONTAINS         

!=====================================================================================
!                  Extract medslik files from OCEAN data
!=====================================================================================

    SUBROUTINE ExtractOCE(fc_dir,indx)
          
          USE env
          USE netcdf

          IMPLICIT NONE

          CHARACTER (LEN = *)                 :: fc_dir
          CHARACTER (LEN = 120)               :: xvar, yvar, zvar, lonvar, latvar, timvar
          CHARACTER (LEN = 120)               :: tvar, uvar, vvar, mess
          CHARACTER (LEN = 2)                 :: scfac, adfac
          INTEGER                             :: indx, startfile, numfile
          INTEGER                             :: start(4), count(4)
          INTEGER                             :: T_id, U_id, V_id, idLON, idLAT
          INTEGER                             :: idT, idU, idV, nwp
          INTEGER                             :: Status
          REAL                                :: sst,us,vs,u10,v10,u30,v30,&
                                                 u120,v120
          REAL                                :: scf_u, off_u, scf_v, off_v,&
                                                 scf_t, off_t
          REAL,DIMENSION(:),ALLOCATABLE       :: lon, lat, oplon, oplat
          REAL,DIMENSION(:,:),ALLOCATABLE     :: lon2D, lat2D
          INTEGER,DIMENSION(:,:),ALLOCATABLE  :: msk
          REAL,DIMENSION(:,:,:),ALLOCATABLE   :: u_tmp, v_tmp, t_tmp
          REAL,DIMENSION(:,:,:),ALLOCATABLE   :: u, v, ts
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: uvel, vvel, potemp
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: u_24, v_24, t_24

!-------------------------------------------------------------------

          udef2 = -9999.
          rhoa = 1.19
          
          OPEN(2,file='files-oce')
          READ(2,*) type_cu
          CLOSE(2)

!         MFS-24dm
          if (indx == 14) then
             xvar = "lon"
             yvar = "lat"
             zvar = "depth"
             lonvar = xvar
             latvar = yvar
             timvar = "time"
             tvar = "votemper"
             uvar = "vozocrtx"
             vvar = "vomecrty"
             scfac = "no"
             adfac = "no"
             startfile = 1
             numfile = numfiles_cu
             miss = 1.e+20
          endif
          

!         MFS-01hm
          IF ( indx == 76 ) THEN
             xvar = "x"
             yvar = "y"
             zvar = "deptht"
             lonvar = "nav_lon"
             latvar = "nav_lat"
             timvar = "time_counter"
             tvar = "votemper"
             uvar = "vozocrtx"
             vvar = "vomecrty"
             startfile = 1
             numfile = numfiles_cu
             miss = 1.e+20
          ENDIF
!         MFS-COP01hm
          IF ( indx == 77 ) THEN
             xvar = "lon"
             yvar = "lat"
             zvar = "depth"
             lonvar = "lon"
             latvar = "lat"
             timvar = "time"
             tvar = "votemper"
             uvar = "vozocrtx"
             vvar = "vomecrty"
             startfile = 2
             numfile = numfiles_cu + 1
             miss = 1.e+20
          ENDIF


!---------------- READING OCEAN DIMENSIONS ------------------

          Status = 0
!          PRINT*, fc_dir
          infile_T = fc_dir//'/OCE/'//indate_cu(startfile)(1:6)//'_T.nc'

          Status = nf90_open(infile_T, nf90_nowrite, T_id)
          IF (Status /= nf90_noerr) call handle_err(Status)

          Status = nf90_inq_dimid(T_id, TRIM(xvar), nx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(T_id, nx, len = imx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inq_dimid(T_id, TRIM(yvar), ny)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(T_id, ny, len = jmx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inq_dimid(T_id, TRIM(zvar), nz)
          IF (status /= nf90_noerr) call handle_err(status)
          
          Status = nf90_inquire_dimension(T_id, nz, len = kmx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inq_dimid(T_id, TRIM(timvar), nt)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(T_id, nt, len = ktmx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_close(T_id)
          IF (Status /= nf90_noerr) call handle_err(Status)

          PRINT*, 'OCE-FILES dimensions: ', imx, jmx, kmx, ktmx
          OPEN(2,file='files-oce')
          WRITE(2,*) type_cu
          WRITE(2,'(i2)') ktmx
          CLOSE(2)

          IF ( (indx.eq.76) ) THEN
             ALLOCATE ( lon2D(imx,jmx) )
             ALLOCATE ( lat2D(imx,jmx) )
          ELSE
             ALLOCATE ( lon(imx) )
             ALLOCATE ( lat(jmx) )
          ENDIF

          ALLOCATE ( oplon(imx) )
          ALLOCATE ( oplat(jmx) )
          ALLOCATE ( msk(imx,jmx) )
          ALLOCATE ( u_tmp(imx,jmx,kmx) )
          ALLOCATE ( v_tmp(imx,jmx,kmx) )
          ALLOCATE ( t_tmp(imx,jmx,kmx) )
          ALLOCATE ( u(imx,jmx,kmx) )
          ALLOCATE ( v(imx,jmx,kmx) )
          ALLOCATE ( ts(imx,jmx,kmx) )
          ALLOCATE ( uvel(imx,jmx,kmx,ktmx) )
          ALLOCATE ( vvel(imx,jmx,kmx,ktmx) )
          ALLOCATE ( potemp(imx,jmx,kmx,ktmx) )
          ALLOCATE ( u_24(imx,jmx,kmx,ktmx) )
          ALLOCATE ( v_24(imx,jmx,kmx,ktmx) )
          ALLOCATE ( t_24(imx,jmx,kmx,ktmx) )

!####################################################################
!                     READ OCEAN DATA FILES
!####################################################################

          DO 60 n = startfile,numfile

              IF ((indx == 76).OR.(indx == 14)) THEN
                 PRINT*, "DAY= ", n
              ELSE
                 PRINT*, "DAY= ", n-1
              ENDIF
               
              Status = 0
              infile_T = fc_dir//'/OCE/'//indate_cu(n)(1:6)//'_T.nc'
              infile_U = fc_dir//'/OCE/'//indate_cu(n)(1:6)//'_U.nc'
              infile_V = fc_dir//'/OCE/'//indate_cu(n)(1:6)//'_V.nc'
              PRINT*, infile_T
              PRINT*, infile_U
              PRINT*, infile_V

!-----------------  OCEAN MEDESS horizontal grid -------------------

!    LONGITUDE FIELD
      
              Status = nf90_open(infile_T, nf90_nowrite, T_id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(T_id, TRIM(lonvar), idLON)
              IF (Status /= nf90_noerr) call handle_err(Status)

              IF (indx.eq.76) THEN
                 Status = nf90_get_var(T_id, idLON, lon2D)
                 IF (Status /= nf90_noerr) call handle_err(Status)
              ELSE
                 Status = nf90_get_var(T_id, idLON, lon)
                 IF (Status /= nf90_noerr) call handle_err(Status)
              ENDIF

              Status = nf90_close(T_id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Longitude in degrees east done'

!    LATITUDE FIELD

              Status = nf90_open(infile_T, nf90_nowrite, T_id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(T_id, TRIM(latvar), idLAT)
              IF (Status /= nf90_noerr) call handle_err(Status)

              IF (indx.eq.76) THEN
                 Status = nf90_get_var(T_id, idLAT, lat2D)
                 IF (Status /= nf90_noerr) call handle_err(Status)
              ELSE
                 Status = nf90_get_var(T_id, idLAT, lat)
                 IF (Status /= nf90_noerr) call handle_err(Status)
              ENDIF

              Status = nf90_close(T_id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Latitude in degrees east done'
          

              IF ( indx.eq.76 ) THEN
                 DO i=1,imx
                    oplon(i) = lon2D(i,1)
                 ENDDO
                
                 DO j=1,jmx
                    oplat(j) = lat2D(1,j)
                 ENDDO
              ELSE
                 DO i=1,imx
                    oplon(i) = lon(i)
                 ENDDO

                 DO j=1,jmx
                    oplat(j) = lat(j)
                 ENDDO
              ENDIF

              op_dlon = oplon(2) - oplon(1)
              op_dlat = oplat(2) - oplat(1)
              oplon0 = oplon(1)
              oplat0 = oplat(1)
             
              i_first = int( (alon1 - oplon0) / op_dlon ) + 1
              i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
              j_first = int( (alat1 - oplat0) / op_dlat ) + 1
              j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

              IF (i_first.lt.1) i_first = 1
              IF (i_last.gt.imx) i_last = imx
              IF (j_first.lt.1) j_first = 1
              IF (j_last.gt.jmx) j_last = jmx

              clon1 = oplon0 + (i_first - 1) * op_dlon
              clon2 = oplon0 + (i_last  - 1) * op_dlon
              clat1 = oplat0 + (j_first - 1) * op_dlat
              clat2 = oplat0 + (j_last  - 1) * op_dlat

              imax = i_last - i_first + 1
              jmax = j_last - j_first + 1

              WRITE(99,*) 'i-limits   = ',i_first,i_last,imax
              WRITE(99,*) 'j-limits   = ',j_first,j_last,jmax
              WRITE(99,*) 'lon-limits = ',clon1,clon2,(clon2-clon1)*(1/op_dlon)
              WRITE(99,*) 'lat-limits = ',clat1,clat2,(clat2-clat1)*(1/op_dlat)
       
!-----------  Reading OCEAN CURRENT and TEMPERATURE fields ---------------

!     ZONAL U VELOCITY FIELD

              status = nf90_open(infile_U, nf90_nowrite, U_id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(U_id, TRIM(uvar), idU)
              IF (status /= nf90_noerr) call handle_err(status)

              count = (/ imx, jmx, kmx, 1 /)
              start = (/ 1, 1, 1, 1 /)

              DO t = 1,ktmx

                 start(4)=t
                 status = nf90_get_var(U_id, idU, u_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 uvel(:,:,:,t) = u_tmp(:,:,:)

              ENDDO

              status = nf90_close(U_id)
              IF (status /= nf90_noerr) call handle_err(status)
              PRINT*,'read zonal u current velocity field done'

!     MERIDIONAL V VELOCITY FIELD

              status = nf90_open(infile_V, nf90_nowrite, V_id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(V_id, TRIM(vvar), idV)
              IF (status /= nf90_noerr) call handle_err(status)

              count = (/ imx, jmx, kmx, 1 /)
              start = (/ 1, 1, 1, 1 /)

              DO t = 1,ktmx

                 start(4)=t
                 status = nf90_get_var(V_id, idV, v_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 vvel(:,:,:,t) = v_tmp(:,:,:)

              ENDDO

              status = nf90_close(V_id)
              IF (status /= nf90_noerr) call handle_err(status)
              PRINT*,'read zonal v current field done'

!     TEMPERATURE FIELD

              status = nf90_open(infile_T, nf90_nowrite, T_id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(T_id, TRIM(tvar), idT)
              IF (status /= nf90_noerr) call handle_err(status)

              count = (/ imx, jmx, kmx, 1 /)
              start = (/ 1, 1, 1, 1 /)

              DO t = 1,ktmx

                 start(4)=t
                 status = nf90_get_var(T_id, idT, t_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 potemp(:,:,:,t) = t_tmp(:,:,:)

              ENDDO

              status = nf90_close(T_id)
              IF (status /= nf90_noerr) call handle_err(status)
              PRINT*,'read temperature field DOne'

              DO t=1,ktmx
                 DO k=1,kmx
                    DO j=1,jmx
                       DO i=1,imx

                          IF (uvel(i,j,k,t).ne.miss) THEN
                             u_24(i,j,k,t)=uvel(i,j,k,t)
                          ELSE
                             u_24(i,j,k,t)=udef2
                          ENDIF
 
                          IF (vvel(i,j,k,t).ne.miss) THEN
                             v_24(i,j,k,t)=vvel(i,j,k,t)
                          ELSE
                             v_24(i,j,k,t)=udef2
                          ENDIF
 
                          IF (potemp(i,j,k,t).ne.miss) THEN
                             t_24(i,j,k,t)=potemp(i,j,k,t)
                          ELSE
                             t_24(i,j,k,t)=udef2
                          ENDIF

                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO

!##########################################################################
!                    CHECK IF THE FILE ALREADY EXISTS
!##########################################################################

              DO t = 1,ktmx

                 IF ( indx == 76 ) THEN

                    mess = 'Writing medslik MFS HOURLY MEAN CURRENT file for date '
                    
                    IF(t.le.12) THEN
                       kount=n
                       nore=t+12
                       WRITE(ore,'(i2)') nore
                       ora=ore(1:2)
                    ENDIF
                    IF (t.gt.12.and.t.lt.22) THEN
                       kount=n+1
                       nore=t-12
                       WRITE(ore,'(i2)') nore
                       ora='0'//ore(2:2)
                    ENDIF
                    IF (t.ge.22.and.t.le.24) THEN
                       kount=n+1
                       nore=t-12
                       WRITE(ore,'(i2)') nore
                       ora=ore(1:2)
                    ENDIF

                 ENDIF

                 IF ( indx == 77 ) THEN

                    mess = 'Writing medslik MFS HOURLY MEAN CURRENT file for date '

                    IF(t.le.12) THEN
                       kount=n-1
                       nore=t+12
                       WRITE(ore,'(i2)') nore
                       ora=ore(1:2)
                    ENDIF
                    IF (t.gt.12.and.t.lt.22) THEN
                       kount=n
                       nore=t-12
                       WRITE(ore,'(i2)') nore
                       ora='0'//ore(2:2)
                    ENDIF
                    IF (t.ge.22.and.t.le.24) THEN
                       kount=n
                       nore=t-12
                       WRITE(ore,'(i2)') nore
                       ora=ore(1:2)
                    ENDIF

                 ENDIF

                 IF (indx == 14) THEN

                    mess = 'Writing medslik MFS-MyO CURRENT file for date '
                    kount=n
                    nore=0
                    WRITE(ore,'(i2)') nore
                    ora='0'//ore(2:2)

                 ENDIF

                 prdate = indate_cu(kount)(5:6)//'/'//indate_cu(kount)(3:4)//&
                          '/20'//indate_cu(kount)(1:2)//' '//ora//':00'
                 WRITE(6,*) TRIM(mess)//" "//TRIM(prdate)
                 WRITE(99,*) TRIM(mess)//" "//TRIM(prdate)

                 outfile = 'INP_DATA/OCE/'//'meds'//indate_cu(kount)(1:6)//ora(1:2)//'.med'

                 INQUIRE(file = TRIM(outfile), EXIST = ex)

                 IF (ex) THEN

                    OPEN(20,file = TRIM(outfile))
                    READ(20,*) empty
                    READ(20,*) empty
                    READ(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1

                    IF (blon1.eq.clon1.and.blon2.eq.clon2.and.blat1.eq.&
                        clat1.and.blat2.eq.clat2.and.imax1.eq.imax.and.&
                        jmax1.eq.jmax) THEN

                        WRITE(6,*) TRIM(outfile)//' already exists for this  subregion'
                        go to 60

                    ENDIF

                    CLOSE(20)

                 ENDIF
!--------------------------------------------------------------------
!                         READ OCEAN data files
!--------------------------------------------------------------------

                 u(:,:,:) = u_24(:,:,:,t)
                 v(:,:,:) = v_24(:,:,:,t)
                 ts(:,:,:) = t_24(:,:,:,t)

        
!####################################################################
!                      MASK & NWP CALCULATION
!####################################################################

! -------------- Mask -------------------

                 DO i=1,imx
                    DO j=1,jmx

                       msk(i,j) = 0
                       IF(ts(i,j,1).gt.udef2) msk(i,j) = 1

                    ENDDO
                 ENDDO

!--------------  nwp = Number of Wet Points -------------

                 nwp = 0

                 DO i=i_first,i_last
                    DO j=j_first,j_last

                       IF (msk(i,j).eq.1) THEN

                           nwp = nwp + 1

                       ENDIF

                    ENDDO
                 ENDDO

                 WRITE(99,*) 'nwp = ',nwp
                 WRITE(99,*) 'mask: '

                 DO j=j_last,j_first,-1

                    WRITE(99,'(300i1)') (msk(i,j),i=i_first,i_last)

                 ENDDO


!##############################################################################
!                  SEA OVER LAND APPLICATION ON ORIGINAL FIELDS
!##############################################################################

                 i1 = i_first-2
                 i2 = i_last+2
                 j1 = j_first-2
                 j2 = j_last+2

                 IF(i1.lt.1) i1 = 1
                 IF(i2.gt.imx) i2 = imx
                 IF(j1.lt.1) j1 = 1
                 IF(j2.gt.jmx) j2 = jmx

                 call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
                 call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
                 call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)

                 IF (indx == 76) THEN

!    MFS old format output are written on staggared grid for U and V volocity
!    components.

                    DO i=i_first,i_last
                       DO j=j_first,j_last
                          IF (msk(i,j).eq.1) THEN
                             DO k=1,kmx
                                IF(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
                                IF(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
                             ENDDO
                          ENDIF
                       ENDDO
                    ENDDO

                 ENDIF

!##############################################################################
!                             WRITE MEDSLIK FILES
!##############################################################################
 
                 outfile = 'INP_DATA/OCE/'//'meds'//indate_cu(kount)(1:6)//&
                            ora(1:2)//'.med'
                 OPEN(20,file = TRIM(outfile))

                 WRITE(20,*) 'OCEAN forecast data for '//TRIM(prdate)
                 WRITE(20,*) 'Mediterranean Sea'
                 WRITE(20,'(4f9.5,2i5,''   Geog. limits'')')&
                                      clon1,clon2,clat1,clat2,imax,jmax
                 WRITE(20,'(i6,''   0.0'')') nwp
                 heads ='    lat        lon        SST        '//&
                        'u_srf      v_srf      u_10m      v_10m'//&
                        '       u_30m      v_30m      u_120m     v_120m'

                 WRITE(20,'(a150)') heads

                 DO i=i_first,i_last
                     DO j=j_first,j_last

                        IF ( msk(i,j).eq.1 ) THEN

                            blon = oplon(i)
                            blat = oplat(j)
 
                            sst = ts(i,j,1)
                            us = u(i,j,1)
                            vs = v(i,j,1)

                            km = 1
                            DO k=2,kmx
                               IF(u(i,j,k).gt.udef2.and.v(i,j,k).gt.udef2) km=k
                            ENDDO

                            IF (km.ge.4) THEN
                               u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615
                               v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615
                            ELSE
                               u10 = u(i,j,km)
                               v10 = v(i,j,km)
                            ENDIF

                            IF (km.ge.9) THEN
                               u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196
                               v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196
                            ELSE
                               u30 = u(i,j,km)
                               v30 = v(i,j,km)
                            ENDIF

                            IF (km.ge.20) THEN
                               u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212
                               v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
                            ELSE
                               u120 = u(i,j,km)
                               v120 = v(i,j,km)
                            ENDIF


                            WRITE(20,'(11f12.8,1i2)') blat,blon,sst,us,vs,u10,v10,&
                                                      u30,v30,u120,v120,msk(i,j)
 !                           WRITE (1,*) i_first,i_last,j_first,j_last

                        ENDIF

                     ENDDO
                  ENDDO
                  CLOSE(20)    
              ENDDO
!              CLOSE(20)
   60     continue
          
          IF ( indx.eq.76 ) THEN
             DEALLOCATE ( lon2D )
             DEALLOCATE ( lat2D )
          ELSE
             DEALLOCATE ( lon )
             DEALLOCATE ( lat )
          ENDIF

          DEALLOCATE ( oplon )
          DEALLOCATE ( oplat )
          DEALLOCATE ( msk )
          DEALLOCATE ( u_tmp )
          DEALLOCATE ( v_tmp )
          DEALLOCATE ( t_tmp )
          DEALLOCATE ( u )
          DEALLOCATE ( v )
          DEALLOCATE ( ts )
          DEALLOCATE ( uvel )
          DEALLOCATE ( vvel )
          DEALLOCATE ( potemp )
          DEALLOCATE ( u_24 )
          DEALLOCATE ( v_24 )
          DEALLOCATE ( t_24 )

          return

    END SUBROUTINE ExtractOCE

!=====================================================================================
!                   Extract MEDSLIK-II files from METEO files
!=====================================================================================

    SUBROUTINE ExtractMET(fc_dir,indx)

          USE env
          USE netcdf

          IMPLICIT NONE

          CHARACTER (LEN = *)                 :: fc_dir
          CHARACTER (LEN = 120)               :: xvar, yvar, wxvar, wyvar, mess
          CHARACTER (LEN = 2)                 :: scfac, adfac
          INTEGER                             :: indx, jstart, jend, step
          INTEGER                             :: start(3), count(3)
          INTEGER                             :: id, idU10, idV10, nwp
          INTEGER                             :: Status
          INTEGER                             :: idLON, idLAT
          REAL                                :: scf_x, off_x, scf_y, off_y                        
          CHARACTER                           :: head1*500, fmw*20
          REAL,DIMENSION(:),      ALLOCATABLE :: oplon, oplat, lon, lat
          INTEGER,DIMENSION(:,:), ALLOCATABLE :: msk
          REAL,DIMENSION(:,:),    ALLOCATABLE :: wx_tmp, wy_tmp
          REAL,DIMENSION(:,:,:),  ALLOCATABLE :: x_wind10, y_wind10
          REAL,DIMENSION(:,:,:),  ALLOCATABLE :: wx, wy
          REAL,DIMENSION(:),      ALLOCATABLE :: wxm, wym


!-------------------------------------------------------------------

          udef = -32766.
          udef2 = -9999.
          rhoa = 1.19
           
          OPEN(3,file='files-met')
          READ(3,*) type_wd
          CLOSE(3)

          IF ( (indx == 11).or.(indx == 12) ) THEN
             xvar  = "lon"
             yvar  = "lat"
             wxvar = "U10M"
             wyvar = "V10M"
          ENDIF
         
          DO 60 n = 1,numfiles_wd

!---------------- READING METEO DIMENSIONS ------------------
         
             Status = 0
             infile_WI = fc_dir//'/MET/'//indate_wd(n)(1:8)//'_WIND.nc'
          
             Status = nf90_open(infile_WI, nf90_nowrite, id)
             IF (Status /= nf90_noerr) call handle_err(Status)

             Status = nf90_inq_dimid(id, xvar, nx)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_inquire_dimension(id, nx, len = imx)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_inq_dimid(id, yvar, ny)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_inquire_dimension(id, ny, len = jmx)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_inq_dimid(id, "time", nt)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_inquire_dimension(id, nt, len = ktmx)
             IF (status /= nf90_noerr) call handle_err(status)

             Status = nf90_close(id)
             IF (Status /= nf90_noerr) call handle_err(Status)

             PRINT*, 'MET-FILES dimensions: ', imx, jmx, ktmx
             OPEN(2,file='files-met')
             WRITE(2,*) type_wd
             WRITE(2,'(i2)') ktmx
             CLOSE(2)

             ALLOCATE ( lon(imx) )
             ALLOCATE ( lat(jmx) )
             ALLOCATE ( oplon(imx) )
             ALLOCATE ( oplat(jmx) )
             ALLOCATE ( msk(imx,jmx) )
             ALLOCATE ( wx_tmp(imx,jmx) )
             ALLOCATE ( wy_tmp(imx,jmx) )
             ALLOCATE ( x_wind10(imx,jmx,ktmx) )
             ALLOCATE ( y_wind10(imx,jmx,ktmx) )
             ALLOCATE ( wx(imx,jmx,ktmx) )
             ALLOCATE ( wy(imx,jmx,ktmx) )
             ALLOCATE ( wxm(ktmx) )
             ALLOCATE ( wym(ktmx) )

!####################################################################
!                     READ METEO DATA FILES
!####################################################################

              Status = 0
              infile_WI = fc_dir//'/MET/'//indate_wd(n)(1:8)//'_WIND.nc'
              PRINT*, infile_WI

!    LONGITUDE FIELD

              Status = nf90_open(infile_WI, nf90_nowrite, id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(id, TRIM(xvar), idLON)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_get_var(id, idLON, lon)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_close(id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Longitude in degrees east done'

!    LATITUDE FIELD

              Status = nf90_open(infile_WI, nf90_nowrite, id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(id, TRIM(yvar), idLAT)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_get_var(id, idLAT, lat)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_close(id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Latitude in degrees east done'


              DO i=1,imx
                 oplon(i) = lon(i)
              ENDDO

              DO j=1,jmx
                 oplat(j) = lat(j)
              ENDDO
             
              op_dlon = oplon(2) - oplon(1)
              IF ( (indx == 11).or.(indx == 12) ) THEN
                 op_dlat = oplat(1) - oplat(2)
              ELSE
                 op_dlat = oplat(2) - oplat(1)
              ENDIF

              oplon0 = oplon(1)
              oplat0 = oplat(1)
        
              i_first = int( (alon1 - oplon0) / op_dlon ) + 1
              i_last  = int( (alon2 - oplon0) / op_dlon ) + 2

              IF ( (indx == 11).or.(indx == 12) ) THEN
                 j_first = int( (oplat0 - alat2) / op_dlat ) + 1
                 j_last  = int( (oplat0 - alat1) / op_dlat ) + 2
              ELSE
                 j_first = int( (alat1 - oplat0) / op_dlat ) + 1
                 j_last  = int( (alat2 - oplat0) / op_dlat ) + 2
              ENDIF

              IF (i_first.lt.1) i_first = 1
              IF (i_last.gt.imx) i_last = imx
              IF (j_first.lt.1) j_first = 1
              IF (j_last.gt.jmx) j_last = jmx

              IF ( (indx == 11).or.(indx == 12) ) THEN
                 clon1 = oplon0 + (i_first - 1) * op_dlon
                 clon2 = oplon0 + (i_last  - 1) * op_dlon
                 clat1 = oplat0 - (j_first - 1) * op_dlat
                 clat2 = oplat0 - (j_last  - 1) * op_dlat              
              ELSE
                 clon1 = oplon0 + (i_first - 1) * op_dlon
                 clon2 = oplon0 + (i_last  - 1) * op_dlon
                 clat1 = oplat0 + (j_first - 1) * op_dlat
                 clat2 = oplat0 + (j_last  - 1) * op_dlat
              ENDIF

              imax = i_last - i_first + 1
              jmax = j_last - j_first + 1

              WRITE(99,*) 'i-limits   = ',i_first,i_last,imax
              WRITE(99,*) 'j-limits   = ',j_first,j_last,jmax
              WRITE(99,*) 'lon-limits = ',clon1,clon2,(clon2-clon1)*(1/op_dlon)
              WRITE(99,*) 'lat-limits = ',clat1,clat2,(clat2-clat1)*(1/op_dlat)

!---------------- Reading METEO wind velocity field -------------------

!     ZONAL U10 VELOCITY FIELD

              status = nf90_open(infile_WI, nf90_nowrite, id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(id, TRIM(wxvar), idU10)
              IF (status /= nf90_noerr) call handle_err(status)


              count = (/ imx, jmx, 1 /)
              start = (/ 1, 1, 1 /)

              DO t = 1,ktmx

                 start(3)=t
                 status = nf90_get_var(id, idU10, wx_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 x_wind10(:,:,t) = wx_tmp(:,:)

              ENDDO

              status = nf90_close(id)
              IF (status /= nf90_noerr) call handle_err(status)
              PRINT*,'read zonal U10 field done'

!     MERIDIONAL V10 VELOCITY FIELD

              status = nf90_open(infile_WI, nf90_nowrite, id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(id, TRIM(wyvar), idV10)
              IF (status /= nf90_noerr) call handle_err(status)

              count = (/ imx, jmx, 1 /)
              start = (/ 1, 1, 1 /)

              DO t = 1,ktmx

                 start(3)=t
                 status = nf90_get_var(id, idV10, wy_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 y_wind10(:,:,t) = wy_tmp(:,:)

              ENDDO

              status = nf90_close(id)
              IF (status /= nf90_noerr) call handle_err(status)
              PRINT*,'read zonal V10M field done'

              DO t=1,ktmx
                 DO j=1,jmx
                    DO i=1,imx

                       IF (x_wind10(i,j,t).ne.miss) THEN
                           wx(i,j,t)=x_wind10(i,j,t)
                       ELSE
                           wx(i,j,t)=udef2
                       ENDIF
                       
                       IF (y_wind10(i,j,t).ne.miss) THEN
                           wy(i,j,t)=y_wind10(i,j,t)
                       ELSE
                           wy(i,j,t)=udef2
                       ENDIF
                      
                    ENDDO
                 ENDDO
              ENDDO

!##########################################################################
!                    CHECK IF THE FILE ALREADY EXISTS
!##########################################################################

              prdate = indate_wd(n)(7:8)//'/'//indate_wd(n)(5:6)//&
                          '/20'//indate_wd(n)(3:4)

              IF ( (indx == 11).or.(indx == 12) ) THEN
                 mess = "Writing medslik ECMWF METEO file for date"
              ENDIF

              WRITE(6,*) "------------------------------------------------"
              WRITE(6,*) TRIM(mess)//" "//TRIM(prdate)
              WRITE(6,*) "------------------------------------------------"
              WRITE(99,*) TRIM(mess)//" "//TRIM(prdate)

              outfile = 'INP_DATA/MET/'//'met_'//indate_wd(n)(3:8)//'.met'

              INQUIRE(file = TRIM(outfile), EXIST = ex)
              IF (ex) THEN

                 OPEN(20,file = TRIM(outfile))
                 READ(20,*) empty
                 READ(20,*) empty
                 READ(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1

                 IF (blon1.eq.clon1.and.blon2.eq.clon2.and.blat1.eq.&
                     clat1.and.blat2.eq.clat2.and.imax1.eq.imax.and.&
                     jmax1.eq.jmax) THEN

                     WRITE(6,*) TRIM(outfile)//' already exists for this  subregion'
                     go to 60

                 ENDIF

                 CLOSE(20)

              ENDIF

!####################################################################
!                      MASK & NWP CALCULATION
!####################################################################

! -------------- Mask -------------------
            
              DO i=1,imx
                 DO j=1,jmx

                    msk(i,j) = 0
                    IF(wx(i,j,1).gt.(udef2)) msk(i,j) = 1

                 ENDDO
              ENDDO

!--------------  nwp = Number of Wet Points -------------

              nwp = 0

              DO i=i_first,i_last
                 DO j=j_first,j_last

                    IF (msk(i,j).eq.1) THEN

                        nwp = nwp + 1

                    ENDIF

                 ENDDO
              ENDDO

              WRITE(99,*) 'nwp = ',nwp

              WRITE(99,*) 'mask: '

              DO j=j_last,j_first,-1

                 WRITE(99,'(300i1)') (msk(i,j),i=i_first,i_last)

              ENDDO

!###############################################################################
!                  SEA OVER LAND APPLICATION ON ORIGINAL FIELDS
!##############################################################################

              i1 = i_first-2
              i2 = i_last+2
              j1 = j_first-2
              j2 = j_last+2

              IF(i1.lt.1) i1 = 1
              IF(i2.gt.imx) i2 = imx
              IF(j1.lt.1) j1 = 1
              IF(j2.gt.jmx) j2 = jmx

              DO t = 1,ktmx

                 call extrap2d(wx(:,:,t), i1, i2, j1, j2, imx, jmx)
                 call extrap2d(wy(:,:,t), i1, i2, j1, j2, imx, jmx)

              ENDDO

!##############################################################################
!                             WRITE MEDSLIK FILES
!##############################################################################

              outfile = 'INP_DATA/MET/'//'met_'//indate_wd(n)(3:8)//'.met'
              OPEN(20,file = TRIM(outfile))

              WRITE(20,*) 'METEO forecast data for '//TRIM(prdate)
              WRITE(20,*) 'Subregion of the Mediterranean with limits:'
              WRITE(20,'(4f9.5,2i5,''   Geog. limits'')')&
                        clon1,clon2,clat1,clat2,imax,jmax
              WRITE(20,'(i6,''   0.0'')') nwp

              IF (ktmx.eq.4) THEN
                  head1 = '     lat       lon       wx00       wy00       wx06       wy06       &
                                wx12       wy12       wx18       wy18'

                  fmw = '(10f12.8)'
              ENDIF
               
              WRITE(20,'(a500)') head1

              IF ( (indx >= 10).and.(indx <= 12) ) THEN
                 jstart = j_last
                 jend   = j_first
                 step   = -1
              ELSE
                 jstart = j_first
                 jend   = j_last
                 step   = +1
              ENDIF
              
              DO i=i_first,i_last
                 DO j=jstart,jend,step

                     blon = oplon(i)
                     blat = oplat(j)

                     DO t = 1,ktmx

                        wxm(t) = wx(i,j,t)
                        wym(t) = wy(i,j,t)
 
                     ENDDO
                     
                     WRITE (20,fmw) blat, blon, (wxm(t), wym(t), t = 1, ktmx)   

                 ENDDO
              ENDDO

              CLOSE(20)

              DEALLOCATE ( lon )
              DEALLOCATE ( lat )
              DEALLOCATE ( oplon )
              DEALLOCATE ( oplat )
              DEALLOCATE ( msk )
              DEALLOCATE ( wx_tmp )
              DEALLOCATE ( wy_tmp )
              DEALLOCATE ( x_wind10 )
              DEALLOCATE ( y_wind10 )
              DEALLOCATE ( wx )
              DEALLOCATE ( wy )
              DEALLOCATE ( wxm )
              DEALLOCATE ( wym )

   60     continue

          return

    END SUBROUTINE ExtractMET

!=====================================================================================
!                  Extract medslik files from WAVE data
!=====================================================================================

    SUBROUTINE ExtractWAV(fc_dir,indx)
          
          USE env
          USE netcdf

          IMPLICIT NONE

          CHARACTER (LEN = *)                 :: fc_dir
          CHARACTER (LEN = 120)               :: xvar, yvar, uvar, vvar, mess
          CHARACTER (LEN = 120)               :: latvar, lonvar
          CHARACTER (LEN = 2)                 :: scfac, adfac
          INTEGER                             :: indx,startfile,finishfile
          INTEGER                             :: start(3), count(3)
          INTEGER                             :: WV_id, idLON, idLAT
          INTEGER                             :: idD, idUSK, idVSK, nwp
          INTEGER                             :: Status
          REAL                                :: scf_u, off_u, scf_v, off_v,&
                                                 scf_d, off_d
          REAL                                :: stoku00, stokv00
          REAL,DIMENSION(:),      ALLOCATABLE :: oplon, oplat, lon, lat
          REAL,DIMENSION(:,:),    ALLOCATABLE :: lonswa, latswa
          INTEGER,DIMENSION(:,:), ALLOCATABLE :: msk, fifmsk
          REAL,DIMENSION(:,:),    ALLOCATABLE :: d_tmp, ustk_tmp, vstk_tmp
          REAL,DIMENSION(:,:,:),  ALLOCATABLE :: wdir, ustk, vstk
          REAL,DIMENSION(:,:,:),  ALLOCATABLE :: dir, dir_rad, dir_x, dir_y
          REAL,DIMENSION(:,:,:),  ALLOCATABLE :: ustok, vstok

!-------------------------------------------------------------------

          IF (indx == 000) THEN
             PRINT*, ''
             PRINT*, '------- NO WAVE MODEL OUTPUT CHOSEN: -------'
             PRINT*, 'the simulation will not use stokes drift or'
             PRINT*, 'it will be calculated directly by MEDSLIK-II'
             PRINT*, 'using JONSWAP spectra. See medslik5.par to '
             PRINT*, 'check what it will be done.'
             PRINT*, '--------------------------------------------'
             PRINT*, ''
             return
          ENDIF


          udef = -32766.
          udef2 = -9999.
          rhoa = 1.19
                    
          OPEN(4,file='files-wav')
          READ(4,*) type_wv
          CLOSE(4)

          IF (indx == 101) THEN
             xvar = "lon"
             yvar = "lat"
             uvar = "sozostdx"
             vvar = "somestdy"
             startfile  = 2
             finishfile = numfiles_wv
             scfac = "no"
             adfac = "no"
             miss = 1.e+20
          ENDIF
        
!---------------- READING WAVE DIMENSIONS ------------------

          Status = 0
!          PRINT*, fc_dir
          infile_WV = fc_dir//'/WAV/'//indate_wv(startfile)(1:8)//'_WAVE.nc'

          Status = nf90_open(infile_WV, nf90_nowrite, WV_id)
          IF (Status /= nf90_noerr) call handle_err(Status)

          Status = nf90_inq_dimid(WV_id, TRIM(xvar), nx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(WV_id, nx, len = imx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inq_dimid(WV_id, TRIM(yvar), ny)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(WV_id, ny, len = jmx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inq_dimid(WV_id, "time", nt)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_inquire_dimension(WV_id, nt, len = ktmx)
          IF (status /= nf90_noerr) call handle_err(status)

          Status = nf90_close(WV_id)
          IF (Status /= nf90_noerr) call handle_err(Status)

          PRINT*, 'WAV-FILES dimensions: ', imx, jmx, ktmx
          OPEN(2,file='files-wav')
          WRITE(2,*) type_wv
          WRITE(2,'(i2)') ktmx
          CLOSE(2)

          ALLOCATE ( lon(imx) )
          ALLOCATE ( lat(jmx) )

          ALLOCATE ( oplon(imx) )
          ALLOCATE ( oplat(jmx) )
          ALLOCATE ( msk(imx,jmx) )
          ALLOCATE ( fifmsk(imx,jmx) )
          ALLOCATE ( ustk(imx,jmx,ktmx) )
          ALLOCATE ( vstk(imx,jmx,ktmx) )
          ALLOCATE ( ustk_tmp(imx,jmx) )
          ALLOCATE ( vstk_tmp(imx,jmx) )
          
!####################################################################
!                     READ WAVE DATA FILES
!####################################################################

          DO 60 n = startfile,finishfile

              IF (indx == 101) THEN
                 PRINT*, "DAY= ", n-1
              ELSE
                 PRINT*, "DAY= ", n
              ENDIF

              Status = 0
!              PRINT*, fc_dir
              infile_WV = fc_dir//'/WAV/'//indate_wv(n)(1:8)//'_WAVE.nc'

              PRINT*, infile_WV

!-----------------  WAVE horizontal grid -------------------

!    LONGITUDE FIELD

              Status = nf90_open(infile_WV, nf90_nowrite, WV_id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(WV_id, TRIM(xvar), idLON)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_get_var(WV_id, idLON, lon)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_close(WV_id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Longitude in degrees east done'

!    LATITUDE FIELD

              Status = nf90_open(infile_WV, nf90_nowrite, WV_id)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_inq_varid(WV_id, TRIM(yvar), idLAT)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_get_var(WV_id, idLAT, lat)
              IF (Status /= nf90_noerr) call handle_err(Status)

              Status = nf90_close(WV_id)
              IF (Status /= nf90_noerr) call handle_err(Status)
              PRINT*,'read Latitude in degrees east done'

              DO i=1,imx
                 oplon(i) = lon(i)
              ENDDO
              DO j=1,jmx
                 oplat(j) = lat(j)
              ENDDO

              op_dlon = oplon(2) - oplon(1)
              op_dlat = oplat(2) - oplat(1)
              oplon0 = oplon(1)
              oplat0 = oplat(1)
              PRINT*,'oplon0', oplon0
              PRINT*,'oplat0', oplat0

              i_first = int( (alon1 - oplon0) / op_dlon ) + 1
              i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
              j_first = int( (alat1 - oplat0) / op_dlat ) + 1
              j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

              IF (i_first.lt.1) i_first = 1
              IF (i_last.gt.imx) i_last = imx
              IF (j_first.lt.1) j_first = 1
              IF (j_last.gt.jmx) j_last = jmx

              clon1 = oplon0 + (i_first - 1) * op_dlon
              clon2 = oplon0 + (i_last  - 1) * op_dlon
              clat1 = oplat0 + (j_first - 1) * op_dlat
              clat2 = oplat0 + (j_last  - 1) * op_dlat

              imax = i_last - i_first + 1
              jmax = j_last - j_first + 1

              WRITE(99,*) 'i-limits   = ',i_first,i_last,imax
              WRITE(99,*) 'j-limits   = ',j_first,j_last,jmax
              WRITE(99,*) 'lon-limits = ',clon1,clon2,(clon2-clon1)*(1/op_dlon)
              WRITE(99,*) 'lat-limits = ',clat1,clat2,(clat2-clat1)*(1/op_dlat)


! -- Reading WAVE stokes_drift fields --


!     ZONAL U-STOKES_DRIFT VELOCITY FIELD 

              status = nf90_open(infile_WV, nf90_nowrite, WV_id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(WV_id, uvar, idUSK)
              IF (status /= nf90_noerr) call handle_err(status)

              PRINT *, 'u-stoke scale-factor and offset', scf_u, off_u
              
              count = (/ imx, jmx, 1 /)
              start = (/ 1, 1, 1 /)

              DO t = 1,ktmx

                 start(3)=t
                 status = nf90_get_var(WV_id, idUSK, ustk_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 ustk(:,:,t) = ustk_tmp(:,:)

              ENDDO

              status = nf90_close(WV_id)
              IF (status /= nf90_noerr) call handle_err(status)

              PRINT*,'READ U-STOKE FIELD DONE'

!     ZONAL V-STOKES_DRIFT VELOCITY FIELD

              status = nf90_open(infile_WV, nf90_nowrite, WV_id)
              IF (status /= nf90_noerr) call handle_err(status)

              status = nf90_inq_varid(WV_id, vvar, idVSK)
              IF (status /= nf90_noerr) call handle_err(status)

              PRINT *, 'v-stoke scale-factor and offset', scf_v, off_v
               
              count = (/ imx, jmx, 1 /)
              start = (/ 1, 1, 1 /)

              DO t = 1,ktmx

                 start(3)=t
                 status = nf90_get_var(WV_id, idVSK, vstk_tmp, start = start, &
                                 count = count)
                 call handle_err(Status)
                 vstk(:,:,t) = vstk_tmp(:,:)

              ENDDO

              status = nf90_close(WV_id)
              IF (status /= nf90_noerr) call handle_err(status)

              PRINT*,'READ V-STOKE FIELD DONE'
 
!----------------------------------------------------------------------------------------

              DO t=1,ktmx
                 DO j=1,jmx
                    DO i=1,imx
                       
                       IF (ustk(i,j,t).ne.miss) THEN
                           ustk(i,j,t)=ustk(i,j,t)
                       ELSE
                           ustk(i,j,t)=udef2
                       ENDIF
 
                       IF (vstk(i,j,t).ne.miss) THEN
                           vstk(i,j,t)=vstk(i,j,t)
                       ELSE
                           vstk(i,j,t)=udef2
                       ENDIF

                    ENDDO
                 ENDDO
              ENDDO
            
              fifmsk(:,:) = ustk(:,:,1)
!###########################################################################

             DO t=1,ktmx

                IF (indx == 101) THEN
                   mess = 'Writing medslik WW3-INGV STOKES DRIFT file for date '
                   WRITE(m,'(i2)') 0
                   mm = '0'//m(2:2)

                   IF (t.le.12) THEN
                      kount=n-1
                      nore=t+12
                      WRITE(ore,'(i2)') nore
                      ora=ore(1:2)
                   ENDIF
                   IF (t.gt.12.and.t.lt.22) THEN
                      kount=n
                      nore=t-12
                      WRITE(ore,'(i2)') nore
                      ora='0'//ore(2:2)
                   ENDIF
                   IF (t.ge.22.and.t.le.24) THEN
                      kount=n
                      nore=t-12
                      WRITE(ore,'(i2)') nore
                      ora=ore(1:2)
                   ENDIF
                ENDIF

                prdate = indate_wv(kount)(7:8)//'/'//indate_wv(kount)(5:6)//&
                         '/20'//indate_wv(kount)(3:4)//' '//ora(1:2)//':'//mm

                WRITE(6,*) TRIM(mess)//' '//TRIM(prdate)
                WRITE(99,*) TRIM(mess)//' '//TRIM(prdate)

                outfile = 'INP_DATA/WAV/'//'wav_'//indate_wv(kount)(3:8)//ora(1:2)//mm(1:2)//'.wav'

                 INQUIRE(file = TRIM(outfile), EXIST = ex)

                 IF (ex) THEN

                    OPEN(20,file = TRIM(outfile))
                    READ(20,*) empty
                    READ(20,*) empty
                    READ(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1

                    IF (blon1.eq.clon1.and.blon2.eq.clon2.and.blat1.eq.&
                       clat1.and.blat2.eq.clat2.and.imax1.eq.imax.and.&
                       jmax1.eq.jmax) THEN

                       WRITE(6,*) TRIM(outfile)//' already exists for this  subregion'
                       go to 60

                    ENDIF

                    CLOSE(20)

                 ENDIF
        
!####################################################################
!                      MASK & NWP CALCULATION
!####################################################################

! -------------- Mask -------------------

                 DO i=1,imx
                    DO j=1,jmx

                       msk(i,j) = 0
                       IF(fifmsk(i,j).gt.udef2+1) msk(i,j) = 1

                    ENDDO
                 ENDDO

!--------------  nwp = Number of Wet Points -------------

                 nwp = 0

                 DO i=i_first,i_last
                    DO j=j_first,j_last

                       IF (msk(i,j).eq.1) THEN

                           nwp = nwp + 1

                       ENDIF

                    ENDDO
                 ENDDO

                 WRITE(99,*) 'nwp = ',nwp
                 WRITE(99,*) 'mask: '

                 DO j=j_last,j_first,-1

                    WRITE(99,'(300i1)') (msk(i,j),i=i_first,i_last)

                 ENDDO

!##############################################################################
!                  SEA OVER LAND APPLICATION ON ORIGINAL FIELDS
!##############################################################################

                 i1 = i_first-2
                 i2 = i_last+2
                 j1 = j_first-2
                 j2 = j_last+2

                 IF(i1.lt.1) i1 = 1
                 IF(i2.gt.imx) i2 = imx
                 IF(j1.lt.1) j1 = 1
                 IF(j2.gt.jmx) j2 = jmx

                 ustk_tmp(:,:) = ustk(:,:,t)
                 vstk_tmp(:,:) = vstk(:,:,t)
                 
                 call extrap2d(ustk_tmp, i1, i2, j1, j2, imx, jmx)
                 call extrap2d(vstk_tmp, i1, i2, j1, j2, imx, jmx)
                 
!##############################################################################
!                             WRITE MEDSLIK FILES
!##############################################################################
 
                 outfile = 'INP_DATA/WAV/'//'wav_'//indate_wv(kount)(3:8)//&
                            ora(1:2)//mm(1:2)//'.wav'
                 OPEN(20,file = TRIM(outfile))

                 WRITE(20,*) 'WAVE forecast data for '//TRIM(prdate)
                 WRITE(20,*) 'Mediterranean Sea'
                 WRITE(20,'(4f9.5,2i5,''   Geog. limits'')')&
                                      clon1,clon2,clat1,clat2,imax,jmax
                 WRITE(20,'(i6,''   0.0'')') nwp
                 heads ='        lat        lon        stokx0        stoky0'

                 WRITE(20,'(a150)') heads

                 DO i=i_first,i_last
                     DO j=j_first,j_last

                        IF ( msk(i,j).eq.1 ) THEN

                            blon = oplon(i)
                            blat = oplat(j)

                            stoku00 = ustk_tmp(i,j)
                            stokv00 = vstk_tmp(i,j) 

                            WRITE(20,'(4f12.8,1i2)') blat,blon,stoku00,stokv00,msk(i,j)
                         
                        ENDIF

                     ENDDO
                  ENDDO
                  CLOSE(20)    
              ENDDO

   60     continue
          
          DEALLOCATE ( lon )
          DEALLOCATE ( lat )
        
          DEALLOCATE ( oplon )
          DEALLOCATE ( oplat )
          DEALLOCATE ( msk )
          DEALLOCATE ( fifmsk )
          DEALLOCATE ( ustk )
          DEALLOCATE ( vstk )
          DEALLOCATE ( ustk_tmp )
          DEALLOCATE ( vstk_tmp )

          return

    END SUBROUTINE ExtractWAV

!************************************************************************
!      Extrapolation of 2-D fields over land points
!************************************************************************

    SUBROUTINE extrap2d(field, i_first, i_last, j_first, j_last, imx, jmx)
 
          IMPLICIT NONE
               
          REAL,DIMENSION(imx,jmx)              :: field
          INTEGER,INTENT(IN)                   :: i_first, i_last, j_first,&
                                                  j_last, imx, jmx
          INTEGER                              :: ngridpts, iter, knt
          INTEGER                              :: im1, ip1, jm1, jp1, jcn,&
                                                  ii, jj
          REAL                                 :: dat
!------------------------------------------------------------------------

          udef = -9999.
          ngridpts = (i_last - i_first + 1) * (j_last - j_first + 1)

          DO iter=1,100

             knt=0

             DO j = j_first, j_last

                DO i = i_first, i_last

                   IF(field(i,j).le.udef) THEN

                      knt = knt + 1
                      im1=i-1
                      ip1=i+1
                      jm1=j-1
                      jp1=j+1

                      IF(im1.lt.i_first) im1 = i_first
                      IF(ip1.gt.i_last)  ip1 = i_last
                      IF(jm1.lt.j_first) jm1 = j_first
                      IF(jp1.gt.j_last)  jp1 = j_last

                      dat=0.
                      jcn=0

                      DO jj=jm1,jp1

                         DO ii=im1,ip1

                            IF (field(ii,jj).gt.udef) THEN

                                dat = dat + field(ii,jj)
                                jcn = jcn + 1

                            ENDIF

                         ENDDO 

                      ENDDO 

                      IF (jcn.gt.2) THEN

                          field(i,j) = dat / float(jcn)

                      ENDIF

                   ENDIF

                ENDDO 

             ENDDO

!            WRITE(99,*) iter,knt

             IF (knt.eq.0.or.knt.eq.ngridpts) return 
        
          ENDDO
                    
          return

    END SUBROUTINE extrap2d        

!************************************************************************
!      Extrapolation of 3-D fields over land points, 7 levels only
!------------------------------------------------------------------------

    SUBROUTINE extrap3d (field, i_first, i_last, j_first, j_last, imx, jmx, kmx)

          IMPLICIT NONE

          REAL,DIMENSION(imx,jmx,kmx)            :: field 
          INTEGER,INTENT(IN)                     :: i_first, i_last, j_first,&
                                                    j_last, imx, jmx, kmx

          REAL,DIMENSION(imx,jmx)                :: tmp
          INTEGER                                :: kd(7), k
          REAL                                   :: udef
          DATA udef /-9999./
          DATA kd /1,3,4,8,9,19,20/

          DO k=1,kmx

             DO i=1,imx

                DO j=1,jmx

                   tmp(i,j) = field(i,j,k)

                ENDDO

             ENDDO
        
             call extrap2d(tmp, i_first, i_last, j_first, j_last, imx, jmx)

             DO i=1,imx

                DO j=1,jmx

                   field(i,j,k) = tmp(i,j)

                ENDDO

             ENDDO
        

          ENDDO   !k

          return
           
    END SUBROUTINE extrap3d

!=====================================================================================
!                           SUBROUTINE HANDLE_ERR(STATUS)
!=====================================================================================

    SUBROUTINE handle_err(status)

         USE netcdf
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: status
         CHARACTER (LEN=80) :: nf_90_strerror

         IF (status /= nf90_noerr) THEN
            WRITE (*,*) nf90_strerror(status)
            stop 'Stopped'
         end IF

    end SUBROUTINE handle_err

!=======================================================================================

 END MODULE utils

!=======================================================================================

 PROGRAM extract

      USE netcdf
      USE env
      USE utils

      IMPLICIT NONE
 
      call getarg(1,fcdir)
      PRINT*, trim(fcdir)
      OPEN(1,file='medslik.tmp')

      READ(1,*) regn, icurrent, iwind, iwave
      READ(1,*) alon1,alon2
      READ(1,*) alat1,alat2
      PRINT*, "alon1= ", alon1, "alon2= ", alon2
      PRINT*, "alat1= ", alat1, "alat2= ", alat2

      READ(1,*) numfiles_cu

      IF (icurrent.eq.14) THEN
         numf_cu=numfiles_cu
      ELSE
         numf_cu=numfiles_cu+1
      ENDIF

      DO n=1,numf_cu
         READ(1,'(a8)') indate_cu(n)
      ENDDO

      READ(1,*) numfiles_wd

      DO n=1,numfiles_wd
         READ(1,'(a8)') indate_wd(n)
      ENDDO
      
      IF (iwave.ne.000) THEN

         READ(1,*) numfiles_wv

         IF (iwave.eq.101) THEN
            numfiles_wv=numfiles_wv+1
         ENDIF
  
         DO n=1,numfiles_wv
            READ(1,'(a8)') indate_wv(n)
            PRINT*, indate_wv(n)
         ENDDO

      ENDIF

      READ(1,*) iviod
      CLOSE(1)
      OPEN(99,file='Extract.log')

      call ExtractOCE(trim(fcdir),icurrent)
      PRINT*,''
      PRINT*,''
      call ExtractMET(trim(fcdir),iwind)
      PRINT*,''
      PRINT*,''
      call ExtractWAV(trim(fcdir),iwave)
      PRINT*,''
      PRINT*,''

      stop

 END PROGRAM extract
