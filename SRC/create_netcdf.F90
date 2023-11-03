
!  create_netcdf.F90
!-----------------------------------------------------------------------------------
!
!
!   --|-------------------------------------------------------------|--
!     | National Institute of Geophyisics and Volcanology (I.N.G.V) |
!     | Bologna Section                                             |
!     | Via Franceschini 11, Bologna, Italy                         |
!     |                                                             |
!     | Programmer: Diego Bruciaferri                               |
!   --|-------------------------------------------------------------|--
!
!
!---------------------------------------------------------------------------------
!  This routine reads MEDSLIK-II ascii formatted outputs and creates,  
!  by the bilinear interpolation and linear interpolation methods, netcdf outputs 
!  of the Oil Concentration, the Ocean Currents, the Stokes' Drift and the Winds 
!  velocities fields.
!
!  Wind velocities are defined for the gravity centre points of the slicks.
!
!  Ocean forcing fields are defined on the original grids  while concentration 
!  is defined only on the 150mx150m grid. 
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
!------------------------------------------------------------------------------------

PROGRAM create_netcdf

      USE module_phymath
      USE module_interpolation
      USE module_netcdf
      USE netcdf
      IMPLICIT NONE
!-------------------------------------------------------------------------------------
!                                 READING INPUT FILES

      CALL getarg(1,out_dir)
      
      IF (TRIM(out_dir) == "") THEN
         PRINT*, " ----------- WARNING !!! ---------- "
         PRINT*, " USAGE: ./create_output.exe OUT_DIR "
         PRINT*, "------------------------------------"
         STOP
      ENDIF

      plg_file = TRIM(out_dir)//plgs                  ! UNIT 12
      pnt_file = TRIM(out_dir)//pnts                  ! UNIT 13
      inp_file = TRIM(out_dir)//inp                   ! UNIT 14
      tmp_file = TRIM(out_dir)//tmp                   ! UNIT 14
      par_file = TRIM(out_dir)//par                   ! UNIT 15
      ocefile  = TRIM(out_dir)//"files-oce"           ! UNIT 16
      metfile  = TRIM(out_dir)//"files-met"           ! UNIT 16
      wavfile  = TRIM(out_dir)//"files-wav"           ! UNIT 16
      cur_dir  = TRIM(out_dir)//"OCE/"
      wav_dir  = TRIM(out_dir)//"WAV/"
      wnd_dir  = TRIM(out_dir)//"MET/"


! READING medslik5.par -----------------------------------------

      OPEN(UNIT=15,FILE=par_file)

      DO i=1,4
         READ(15,*) EMPTY
      ENDDO
      READ(15,*) istoke
      PRINT*, "ISTOKE: ", istoke

      CLOSE(15)

! READING medslik5.inp -----------------------------------------

     OPEN(14,file = TRIM(inp_file))

     DO i=1,4
        READ(14,*) empty
     ENDDO
     READ(14,*) numspills
     READ(14,'(I2,1X,I2,1X,I4)') dd,mm,yr
     READ(14,'(I4)') hr
     READ(14,'(A4)') outname
     READ(14,*) sim_length
     READ(14,'(I3)') delt_out
     READ(14,'(I2)') icurrents
     READ(14,'(I2)') iwind
     READ(14,'(I3)') iwave
     DO i=1,3
        READ(14,*) empty
     ENDDO
     READ(14,*) isat
     CLOSE(14)

     nc_file  = TRIM(out_dir)//TRIM(outname)//".nc"

! READING medslik.tmp -----------------------------------------

     OPEN(14,file = TRIM(tmp_file))
     DO i=1,3
        READ(14,*) empty
     ENDDO

     READ(14,*) cur_num

     IF (icurrents == 14) THEN
        cur_num = cur_num
     ELSE
        cur_num = cur_num + 1
     ENDIF 

     ALLOCATE( cur_file(cur_num) )

     DO n=1,cur_num
        READ(14,'(a8)') cur_file(n)
        print*, cur_file(n)
     ENDDO
     
     IF (icurrents /= 14) THEN
        cur_num = cur_num - 1
     ENDIF

     READ(14,*) wnd_num
     ALLOCATE( wnd_file(wnd_num) )

     DO n=1,wnd_num
        READ(14,'(a8)') wnd_file(n)
     ENDDO

     IF (iwave /= 000) THEN
        READ(14,*) wav_num
       
        ALLOCATE( wav_file(wav_num+1) )

        DO n=1,wav_num+1
           READ(14,'(a8)') wav_file(n)
           wav_file(n) = wav_file(n)(3:8)
           print*, wav_file(n)
        ENDDO
     ENDIF

     CLOSE(14)

! READING ocefile, metfile, wavfile

      OPEN(UNIT=16,FILE=ocefile)
      READ(16,'(a3)') type_cu
      READ(16,'(i2)') nrec_cu
      CLOSE(16)

      OPEN(UNIT=16,FILE=metfile)
      READ(16,'(a3)') type_wd
      READ(16,'(i2)') nrec_wd
      CLOSE(16)

      IF (iwave /= 000) THEN
         OPEN(UNIT=16,FILE=wavfile)
         READ(16,'(a3)') type_wv
         READ(16,'(i2)') nrec_wv
      ENDIF

! ------------------------------------------------------------------------

      nstph = 1.d0 / delt_out
      sim_start = DFLOAT(hr / 100) + DFLOAT(hr - 100 * (hr / 100)) / 60.d0
      PRINT*, "SIM_START: ", sim_start
      mxst = INT(sim_length * nstph)
      PRINT*, "MXST: ", mxst
      ntime = mxst + 1
 
!------------------------------------------------------------------------
!                          CURRENT FILES

      ALLOCATE( cur_name(nrec_cu*cur_num) )
      ALLOCATE( cur_time(nrec_cu*cur_num) )
      print*, "nrec_cu: ", nrec_cu
      print*, "cur_num: ", cur_num
     
      DO i = 1,cur_num

         IF (icurrents == 14) THEN
             cur_name(i) = "meds"//cur_file(i)(1:6)//"00.med"
             print*, cur_name(i)
         ENDIF

         DO ihr = 1,nrec_cu

            IF ((icurrents /= 76).AND.(icurrents /= 77)) THEN 
               IF (nrec_cu.eq.24) THEN
! ------ 1hr resolution (forecast) -------
                  IF (type_cu.eq.' hm') THEN
! Hourly mean fields
                     WRITE (ihr_str,'(i2)') ihr
                     IF (ihr.lt.10) THEN
                        ihr_str='0'//ihr_str(2:2)
                     ENDIF
                  ELSE
! Instantaneous fields
                     WRITE(ihr_str,'(i2)') ihr-1                                      
                     IF (ihr.le.10) THEN
                        ihr_str='0'//ihr_str(2:2)
                     ENDIF
                  ENDIF
               ELSE
! ------ Different time resolution (analysis/historical data) -------
                  thr = ((24 / nrec_cu)*(ihr-1))
! Hourly mean fields
                  IF (type_cu.eq.' hm') THEN
                     WRITE(ihr_str,'(i2)') thr+1
                     IF (thr.lt.9) then
                        ihr_str='0'//ihr_str(2:2)
                     ENDIF
! Instantaneous fields
                  ELSE
                     WRITE(ihr_str,'(i2)') thr
                     IF (thr.lt.10) THEN
                        ihr_str='0'//ihr_str(2:2)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

            IF (icurrents >= 76) THEN
               IF (ihr <= 12) THEN
                  kount=i
                  WRITE(ihr_str,'(i2)') ihr+12
               ENDIF
               IF (ihr > 12 .AND.  ihr < 22) THEN
                  kount=i+1
                  WRITE(ihr_str,'(i2)') ihr-12
                  ihr_str = '0'//ihr_str(2:2)
               ENDIF
               IF (ihr >= 22 .AND. ihr <= 24) THEN
                  kount=i+1
                  WRITE(ihr_str,'(i2)') ihr-12
               ENDIF
            ENDIF

            IF (icurrents >= 76) THEN
               cur_name((i-1)*nrec_cu+ihr) = "meds"//cur_file(kount)(1:6)//ihr_str//".med"
            ELSE
               cur_name((i-1)*nrec_cu+ihr) = "meds"//cur_file(i)(1:6)//ihr_str//".med"
            ENDIF
         ENDDO
      ENDDO

! COMPUTATION OF CURRENT FORCING TIME RESPECT TO SIMULATION START

      DO i=1,cur_num*nrec_cu

         a2 = cur_name(i)(5:6)
         READ(a2,'(i2)') iyr
         a2 = cur_name(i)(7:8)
         READ(a2,'(i2)') imm
         a2 = cur_name(i)(9:10)
         READ(a2,'(i2)') idd
         a2 = cur_name(i)(11:12)
         READ(a2,'(i2)') ihr

         nday = JDIFF(dd,mm,yr,idd,imm,iyr+2000)
         IF (type_cu.eq.' hi') THEN
            cur_time(i) = nday*24.d0 + DFLOAT(ihr) !snapshot
         ENDIF
         IF (type_cu.eq.' hm') THEN
            cur_time(i) = nday*24.d0 + (DFLOAT(ihr)-0.5) !hourly mean
         ENDIF
         IF (type_cu.eq.' dm') THEN
            cur_time(i) = nday*24.d0 + (DFLOAT(ihr)) !daily mean
            print*, cur_time(i)
         ENDIF

      ENDDO

! READ COORDINATES OF CURRENTS FIELDS 

      OPEN(UNIT=17,FILE=TRIM(cur_dir)//TRIM(cur_name(1)))
      print*, TRIM(cur_dir)//TRIM(cur_name(1))
      READ(17,*) EMPTY
      READ(17,*) EMPTY
      READ(17,*) EMPTY
      READ(17,*) cndata
      READ(17,*) EMPTY

      ALLOCATE( raw_lat(cndata) )
      ALLOCATE( raw_lon(cndata) )
      ALLOCATE( lat_sort(cndata,2) )
      ALLOCATE( lon_sort(cndata,2) )

      DO i=1,cndata
         READ(17,*) raw_lat(i), raw_lon(i)
      ENDDO
      print*, raw_lat(1), raw_lon(1)

      lat_sort = 0.
      lon_sort = 0.

      lat_sort(:,1) = Sort(raw_lat)
      lon_sort(:,1) = Sort(raw_lon)

      nc = 0
      mc = 0

      DO i=2,cndata
         IF (lat_sort(i,1) == lat_sort(i-1,1)) THEN
            nc = nc + 1
            lat_sort(i,2) = 1.
         ENDIF
         IF (lon_sort(i,1) == lon_sort(i-1,1)) THEN
            mc = mc + 1
            lon_sort(i,2) = 1.
         ENDIF
      ENDDO

      cjmx = cndata - nc
      cimx = cndata - mc
      print*, cimx, cjmx
      ALLOCATE( lon_cur(cimx) )
      ALLOCATE( lat_cur(cjmx) )

      lat_cur = PACK(lat_sort(:,1), lat_sort(:,2) == 0.)
      lon_cur = PACK(lon_sort(:,1), lon_sort(:,2) == 0.)
      print*, lat_cur
      print*, size(lat_cur,1)
      print*, size(lon_cur,1)

      clat1 = lat_cur(1)
      clat2 = lat_cur(cjmx)
      clon1 = lon_cur(1)
      clon2 = lon_cur(cimx)

      cdlon = (clon2 - clon1) / DFLOAT(cimx - 1)
      cdlat = (clat2 - clat1) / DFLOAT(cjmx - 1)

      cdlon = lon_cur(2) - lon_cur(1)
      cdlat = lat_cur(2) - lat_cur(1)

      PRINT*, "clon1, ",clon1,"clon2, ",clon2
      PRINT*, "clat1, ",clat1,"clat2, ",clat2
      PRINT*, "cdlat, ",cdlat,"cdlon, ",cdlon

      DEALLOCATE(raw_lat)
      DEALLOCATE(raw_lon)
      DEALLOCATE(lat_sort)
      DEALLOCATE(lon_sort)

      CLOSE(17)

      ALLOCATE( us(cimx,cjmx) )
      ALLOCATE( vs(cimx,cjmx) )
      ALLOCATE( dus(cimx,cjmx) )
      ALLOCATE( dvs(cimx,cjmx) )
      ALLOCATE( uhr(cimx,cjmx,mxst+1) )
      ALLOCATE( vhr(cimx,cjmx,mxst+1) )
 
! ------------------------------------------------------------------------
!                                METEO FILES

      ALLOCATE( wnd_name(wnd_num) )
      ALLOCATE( wnd_time(wnd_num,nrec_wd) )

      DO i=1,wnd_num
             wnd_name(i) = "met_"//wnd_file(i)(3:8)//".met"
      ENDDO

      dt_wnd = 24.d0 / DFLOAT(nrec_wd)

!----- COMPUTATION OF WIND FORCING TIME RESPECT TO SIMULATION START -----

      DO i = 1,wnd_num

         a2 = wnd_name(i)(5:6)
         READ(a2,'(i2)') iyr
         a2 = wnd_name(i)(7:8)
         READ(a2,'(i2)') imm
         a2 = wnd_name(i)(9:10)
         READ(a2,'(i2)') idd

         nday = JDIFF(dd,mm,yr,idd,imm,iyr+2000)

         DO k = 1,nrec_wd
            wnd_time(i,k) = nday * 24.d0 + (k-1) * dt_wnd
         END DO

      ENDDO

! READ COORDINATES OF WIND FIELDS -------------

      OPEN(UNIT=19,FILE=TRIM(wnd_dir)//TRIM(wnd_name(1)))

      READ(19,*) EMPTY
      READ(19,*) EMPTY
      READ(19,*) EMPTY
      READ(19,*) wdndata
      READ(19,*) EMPTY

      ALLOCATE( raw_lat(wdndata) )
      ALLOCATE( raw_lon(wdndata) )
      ALLOCATE( lat_sort(wdndata,2) )
      ALLOCATE( lon_sort(wdndata,2) )

      DO i=1,wdndata
         READ(19,"(2(F12.8))") raw_lat(i), raw_lon(i)       
      ENDDO

      lat_sort = 0.
      lon_sort = 0.

      lat_sort(:,1) = Sort(raw_lat)
      lon_sort(:,1) = Sort(raw_lon)

      nc = 0
      mc = 0

      DO i=2,wdndata
         IF (lat_sort(i,1) == lat_sort(i-1,1)) THEN
             nc = nc + 1
             lat_sort(i,2) = 1.
         ENDIF
         IF (lon_sort(i,1) == lon_sort(i-1,1)) THEN
            mc = mc + 1
            lon_sort(i,2) = 1.
         ENDIF
      ENDDO

      wdjmx = wdndata - nc
      wdimx = wdndata - mc

      ALLOCATE( lon_wnd(wdimx) )
      ALLOCATE( lat_wnd(wdjmx) )

      lat_wnd = PACK(lat_sort(:,1), lat_sort(:,2) == 0.)
      lon_wnd = PACK(lon_sort(:,1), lon_sort(:,2) == 0.)

      wdlat1 = lat_wnd(1)
      wdlat2 = lat_wnd(wdjmx)
      wdlon1 = lon_wnd(1)
      wdlon2 = lon_wnd(wdimx)

      wddlon = (wdlon2 - wdlon1) / DFLOAT(wdimx - 1)
      wddlat = (wdlat2 - wdlat1) / DFLOAT(wdjmx - 1)
      PRINT*, "wdlon1, ",wdlon1,"wdlon2, ",wdlon2
      PRINT*, "wdlat1, ",wdlat1,"wdlat2, ",wdlat2
      PRINT*, "wddlat, ",wddlat,"wddlon, ",wddlon

      DEALLOCATE(raw_lat)
      DEALLOCATE(raw_lon)
      DEALLOCATE(lat_sort)
      DEALLOCATE(lon_sort)

      CLOSE(19)

      ALLOCATE( wndx(wdimx,wdjmx) )
      ALLOCATE( wndy(wdimx,wdjmx) )
      ALLOCATE( dwndx(wdimx,wdjmx) )
      ALLOCATE( dwndy(wdimx,wdjmx) )
      ALLOCATE( wndxhr(wdimx,wdjmx,mxst+1) )
      ALLOCATE( wndyhr(wdimx,wdjmx,mxst+1) )

! ------------------------------------------------------------------------
!                                 WAVE FILES

      IF ( iwave /= 000 ) THEN

         ALLOCATE( wav_name(nrec_wv*wav_num) )
         ALLOCATE( wav_time(nrec_wv*wav_num) )

         DO i = 1,wav_num

            DO ihr=1,nrec_wv

               IF (iwave /= 101) THEN 
                  IF (nrec_wv.eq.24) THEN
! ------ 1hr resolution (forecast) -------
                     IF (type_wv.eq.' hm') THEN
! Hourly mean fields
                         WRITE(ihr_str,'(i2)') ihr
                         IF (ihr.lt.10) THEN
                            ihr_str='0'//ihr_str(2:2)
                         ENDIF
                     ELSE
! Instantaneous fields
                         WRITE(ihr_str,'(i2)') ihr-1
                         IF (ihr.le.10) THEN
                            ihr_str='0'//ihr_str(2:2)
                         ENDIF
                     ENDIF
                  ELSE
                     thr = ((24 / nrec_wv)*(ihr-1))
! Hourly mean fields
                     IF (type_cu.eq.' hm') THEN
                        WRITE(ihr_str,'(i2)') thr+1
                        IF (thr.lt.9) then
                           ihr_str='0'//ihr_str(2:2)
                        ENDIF
! Instantaneous fields
                     ELSE
                        WRITE(ihr_str,'(i2)') thr
                        IF (thr.lt.10) THEN
                           ihr_str='0'//ihr_str(2:2)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF

               IF (iwave == 101) THEN
                  IF (ihr <= 12) THEN
                     kount=i
                     WRITE(ihr_str,'(i2)') ihr+12
                     ihr_str=ihr_str(1:2)
                  ENDIF
                  IF (ihr > 12 .AND. ihr < 22) THEN
                     kount=i+1
                     WRITE(ihr_str,'(i2)') ihr-12
                     ihr_str='0'//ihr_str(2:2)
                  ENDIF
                  IF (ihr >= 22 .AND. ihr <= 24) THEN
                     kount=i+1
                     WRITE(ihr_str,'(i2)') ihr-12
                     ihr_str=ihr_str(1:2)
                  ENDIF
               ENDIF
               IF (iwave == 101) THEN
                  wav_name((i-1)*nrec_wv+ihr) = "wav_"//wav_file(kount)(1:6)//ihr_str//"00.wav"
               ELSE
                  wav_name((i-1)*nrec_wv+ihr) = "wav_"//wav_file(i)(1:6)//ihr_str//"00.wav"
               ENDIF
     
            ENDDO
         ENDDO
         

!----- COMPUTATION OF WAVE FORCING TIME RESPECT TO SIMULATION START -----

         DO i=1,wav_num*nrec_wv

            a2 = wav_name(i)(5:6)
            READ(a2,'(i2)') iyr
            a2 = wav_name(i)(7:8)
            READ(a2,'(i2)') imm
            a2 = wav_name(i)(9:10)
            READ(a2,'(i2)') idd
            a2 = wav_name(i)(11:12)
            READ(a2,'(i2)') ihr

            nday = JDIFF(dd,mm,yr,idd,imm,iyr+2000)
            IF (type_wv.eq.' hi') THEN
               wav_time(i) = nday*24.d0 + DFLOAT(ihr) !snapshot
            ENDIF
            IF (type_wv.eq.' hm') THEN
               wav_time(i) = nday*24.d0 + (DFLOAT(ihr)-0.5) !hourly mean
            ENDIF

         ENDDO


! READ COORDINATES OF WAVE FIELDS -------------

         OPEN(UNIT=18,FILE=TRIM(wav_dir)//TRIM(wav_name(1)))

         READ(18,*) EMPTY
         READ(18,*) EMPTY
         READ(18,*) EMPTY
         READ(18,*) wndata
         READ(18,*) EMPTY

         ALLOCATE( raw_lat(wndata) )
         ALLOCATE( raw_lon(wndata) )
         ALLOCATE( lat_sort(wndata,2) )
         ALLOCATE( lon_sort(wndata,2) )

         DO i=1,wndata
            READ(18,*) raw_lat(i), raw_lon(i)
         ENDDO

         lat_sort = 0.
         lon_sort = 0.

         lat_sort(:,1) = Sort(raw_lat)
         lon_sort(:,1) = Sort(raw_lon)

         nc = 0
         mc = 0

         DO i=2,wndata
            IF (lat_sort(i,1) == lat_sort(i-1,1)) THEN
               nc = nc + 1
               lat_sort(i,2) = 1.
            ENDIF
            IF (lon_sort(i,1) == lon_sort(i-1,1)) THEN
               mc = mc + 1
               lon_sort(i,2) = 1.
            ENDIF
         ENDDO

         wjmx = wndata - nc
         wimx = wndata - mc

         ALLOCATE( lon_wav(wimx) )
         ALLOCATE( lat_wav(wjmx) )

         lat_wav = PACK(lat_sort(:,1), lat_sort(:,2) == 0.)
         lon_wav = PACK(lon_sort(:,1), lon_sort(:,2) == 0.)

         wlat1 = lat_wav(1)
         wlat2 = lat_wav(wjmx)
         wlon1 = lon_wav(1)
         wlon2 = lon_wav(wimx)
         wdlon = (wlon2 - wlon1) / DFLOAT(wimx - 1)
         wdlat = (wlat2 - wlat1) / DFLOAT(wjmx - 1)

         PRINT*, "wlon1, ",wlon1,"wlon2, ",wlon2
         PRINT*, "wlat1, ",wlat1,"wlat2, ",wlat2
         PRINT*, "wdlat, ",wdlat,"wdlon, ",wdlon

         DEALLOCATE(raw_lat)
         DEALLOCATE(raw_lon)
         DEALLOCATE(lat_sort)
         DEALLOCATE(lon_sort)

         CLOSE(18)

         ALLOCATE( uw(wimx,wjmx) )
         ALLOCATE( vw(wimx,wjmx) )
         ALLOCATE( duw(wimx,wjmx) )
         ALLOCATE( dvw(wimx,wjmx) )
         ALLOCATE( uwhr(wimx,wjmx,mxst+1) )
         ALLOCATE( vwhr(wimx,wjmx,mxst+1) )
     
      ENDIF
!------------------------------------------------------------------------
! READING OIL SLICK/SPILL INITIAL CONDITION

      IF (numspills == 1) THEN
         ! SINGLE POINT SOURCE
         IF (isat == 0) THEN

           ALLOCATE( lat_oil(1) )
           ALLOCATE( lon_oil(1) )

           OPEN(13,file = TRIM(pnt_file))

           DO i=1,5
              READ(13,*) empty
           ENDDO
           READ(13,'(F8.4)') lat_oil(1)
           READ(13,'(F8.4)') lon_oil(1)
           ini_point = 1
           
           CLOSE(13)
         ! SINGLE POLYGON SLICK
         ELSE IF (isat == 1) THEN

           OPEN(12,file = TRIM(plg_file))

           DO i=1,5
              READ(12,*) empty
           ENDDO
           READ(12,*) ini_point
           READ(12,*) empty 

           ALLOCATE( lat_oil(ini_point) )
           ALLOCATE( lon_oil(ini_point) )

           DO i=1,ini_point
              READ(12,'(F8.4)') lat_oil(i), lon_oil(i)
           ENDDO

           CLOSE(12)
         ENDIF

      ELSE IF (numspills > 1) THEN
         ! MULTIPLE POINT SOURCES
         IF (isat == 0) THEN

           ALLOCATE( lat_oil(numspills) )
           ALLOCATE( lon_oil(numspills) )

           OPEN(13,file = TRIM(pnt_file))

!           READ(13,*) empty

           DO i=1,numspills
              DO j=1,5
                 READ(13,*) empty
              ENDDO
              READ(13,'(F8.4)') lat_oil(i)
              READ(13,'(F8.4)') lon_oil(i)
           ENDDO

           ini_point = numspills 

           CLOSE(13)

        ! MULTIPLE POLYGON SLICKS
         ELSE IF (isat == 1) THEN

           OPEN(12,file = TRIM(plg_file))

           DO i=1,4
              READ(12,*) empty
           ENDDO
           DO i=1,numspills
              READ(12,*) empty
           ENDDO
           READ(12,*) ini_point
           READ(12,*) empty

           ALLOCATE( lat_oil(ini_point) )
           ALLOCATE( lon_oil(ini_point) )

           DO i=1,ini_point
              READ(12,'(2F8.4)') lat_oil(i), lon_oil(i)
           ENDDO

           CLOSE(12)

         ENDIF

      ENDIF

! -------------------------------------------------------------------------
! REGULAR OIL GRID DEFINITION

      ref_lat = lat_oil(1)

      CALL READ_MASK150
      CALL OIL_GRID(msk,lon_msk,lat_msk,nx,ny)
      CALL COAST_SEG(min_map_lon,max_map_lon,min_map_lat,max_map_lat)

      PRINT*, "MAP_LON_N: ", map_lon_n + 1
      PRINT*, "MAP_LAT_N: ", map_lat_n + 1

! MATRIXES FOR OIL CONCENTRATION

      ALLOCATE(surf_map(map_lon_n+1,map_lat_n+1,ntime))
      ALLOCATE(displ_map(map_lon_n+1,map_lat_n+1,ntime))

      IF ( nseg /= 0 ) THEN
          ALLOCATE(coast_map(nseg,ntime))
      ELSE
          ALLOCATE(coast_map(1,ntime))
      ENDIF

      surf_map  = 0.
      displ_map = 0.
      coast_map = 0.

      OPEN(UNIT=21,FILE=p_coord)

! TIME LOOP **************************************
      
      DO t=0,mxst

! OIL CONCENTRATION OUTPUT -----------

        PRINT*, 't= ', t
        READ(21,*) EMPTY
        READ(21,*) EMPTY, sea_point
        READ(21,*) EMPTY

        IF ( t == 0 ) THEN
           PRINT*, "--------------------------"
           PRINT*, "OIL CONC INITIAL CONDITION"
        ELSE
           PRINT*, "OIL CONC on REGULAR FIXED GRID INTERPOLATED for TIMESTEP ", t
        ENDIF

        ALLOCATE(lat_oil(sea_point))
        ALLOCATE(lon_oil(sea_point))
        ALLOCATE(conc_oil(sea_point))
        ALLOCATE(status_p(sea_point))

        DO i=1,sea_point
           READ(21,*) lat_oil(i), lon_oil(i), status_p(i), conc_oil(i)
        ENDDO

        IF ( t /= 0 ) THEN

           PRINT*,""
           PRINT*, "--- SEARCHING THE OIL CONC NEAREST POINT TO EACH FINAL GRID POINT ----"
 
           ALLOCATE(oil_indx(map_lon_n+1,map_lat_n+1))
            
           DO i=1,map_lon_n+1
              DO j=1,map_lat_n+1
                 CALL Near2D(sea_point, lon_oil, lat_oil, map_lon(i), &
                            map_lat(j), oil_indx(i,j), distance)
                 k = oil_indx(i,j)
                 IF ( distance > eps ) k = -1
                 IF ( k /= -1 ) THEN
                    IF (status_p(k) == 5) surf_map(i,j,t+1) = conc_oil(k)    ! surface
                    IF (status_p(k) == 1) displ_map(i,j,t+1) = conc_oil(k)  ! dispersed
                 ENDIF
              ENDDO
           ENDDO

! OIL CONCENTRATION ON COASTAL SEGMENTS 

           IF ( nseg /= 0 ) THEN

              a0 = ".cst"
              WRITE(a4,'(i4)') t
              IF ( t < 10 ) THEN
                 a4 = "000"//a4(4:4)
              ELSE IF (( t >= 10 ).AND.(t < 100)) THEN
                 a4 = "00"//a4(3:4)
              ELSE IF (t >= 100) THEN
                 a4 = "0"//a4(2:4)
              ENDIF
              
              out_file = TRIM(out_dir)//TRIM(outname)//a4//a0
              OPEN (83,file=TRIM(out_file))
              DO i=1,5
                 READ(83,*) EMPTY
              ENDDO
              DO i=1,numspills     
                 READ(83,*) EMPTY
              ENDDO
              READ(83,*) EMPTY
              READ(83,*) EMPTY
              READ(83,*) cst_point
              READ(83,*) EMPTY

              ALLOCATE(lat_cst(cst_point,2))
              ALLOCATE(lon_cst(cst_point,2))
              
              DO i=1,cst_point
                 READ(83,*) lat_cst(i,1), lon_cst(i,1), lat_cst(i,2), &
                               lon_cst(i,2), conc_oil(i)
              ENDDO

              DO i=1,nseg
                 DO k=1,cst_point
                    PRINT*, i
                    
                       IF ( (lon_cst(k,1) == cst_seg(i,1)).AND.&
                            (lat_cst(k,1) == cst_seg(i,2)).AND.&
                            (lon_cst(k,2) == cst_seg(i,3)).AND.&
                            (lat_cst(k,2) == cst_seg(i,4)) ) THEN
                          coast_map(i,t+1) = conc_oil(k) 
                       ENDIF
                    
                 ENDDO
              ENDDO
              DEALLOCATE(lat_cst)
              DEALLOCATE(lon_cst)
              CLOSE(83)
           ENDIF

           DEALLOCATE(oil_indx)

           PRINT*, "--- REGRIDDING OF OIL ON THE NEW REGULAR GRID DONE ----"
           PRINT*, ""

        ELSE
! INITIAL CONDITION
           ALLOCATE(oil_lat_indx(ini_point))
           ALLOCATE(oil_lon_indx(ini_point))

           DO k=1,ini_point

              PRINT*, 'lat_oil= ', lat_oil(k), 'lon_oil= ', lon_oil(k)

              oil_lat_indx(k) = Near(map_lat_n,map_lat,lat_oil(k))
              oil_lon_indx(k) = Near(map_lon_n,map_lon,lon_oil(k))

              PRINT*, 'oil_lat_indx= ', oil_lat_indx(k), 'oil_lon_indx= ', oil_lon_indx(k)

           ENDDO

           DO k=1,ini_point
              DO i=1,map_lon_n
                 DO j=1,map_lat_n

                    IF ( (i == oil_lon_indx(k)).AND.(j == oil_lat_indx(k)).AND.(status_p(k) == 5) ) THEN
                       surf_map(i,j,t+1) = conc_oil(k)
                    ENDIF

                 ENDDO
              ENDDO
           ENDDO

           DEALLOCATE(oil_lat_indx)
           DEALLOCATE(oil_lon_indx)

        ENDIF

        DEALLOCATE(lat_oil)
        DEALLOCATE(lon_oil)
        DEALLOCATE(conc_oil)
        DEALLOCATE(status_p)

! ==========================================================================

! METEO-OCEAN FORCING FIELDS LINEAR INTERPOLATION IN TIME

        timehr = t * delt_out
        time   = sim_start + timehr

        PRINT*, ""
        PRINT*, "METEO-OCEAN FORCING FIELDS LINEARLY INTERPOLATED IN TIME"
        PRINT*, ""
        PRINT*, ", time: ", time

        dt_cur = 24.d0 / DFLOAT(nrec_cu)
        dt_wav = 24.d0 / DFLOAT(nrec_wv)
        dt_wnd = 24.d0 / DFLOAT(nrec_wd)

        frac_c  = dt_cur / delt_out
        frac_w  = dt_wav / delt_out
        frac_wd = dt_wnd / delt_out

        IF ( t == 0 ) THEN

           DO i=1,cur_num*nrec_cu
              IF(time.ge.cur_time(i)) cfile=i
           ENDDO

           DO i=1,wav_num*nrec_wv
              IF(time.ge.wav_time(i)) wfile=i
           ENDDO

           readdata=.TRUE.
           wdfile=1
           DO k=1,nrec_wd
              IF (time.ge.wnd_time(wdfile,k)) THEN
                 wndrec=k
                 EXIT
              ENDIF
           ENDDO
        ENDIF

! OCEAN CURRENT TIME INTERPOLATION

        IF (cfile <= cur_num*nrec_cu) THEN

           IF ( time.ge.cur_time(cfile) ) THEN

              fac_c = (time - cur_time(cfile)) / delt_out

              PRINT*, "Reading forecast currents from file ", TRIM(cur_name(cfile))
              us  =READFCST_NEW(cur_dir,cur_name(cfile),lon_cur,lat_cur,"c","u",cndata,cimx,cjmx)
              vs = READFCST_NEW(cur_dir,cur_name(cfile),lon_cur,lat_cur,"c","v",cndata,cimx,cjmx)

              cfile = cfile + 1

              IF (cfile <= cur_num*nrec_cu) THEN
                 PRINT*, "Reading forecast currents from file ", TRIM(cur_name(cfile))
                 dus = READFCST_NEW(cur_dir,cur_name(cfile),lon_cur,lat_cur,"c","u",cndata,cimx,cjmx)
                 dvs = READFCST_NEW(cur_dir,cur_name(cfile),lon_cur,lat_cur,"c","v",cndata,cimx,cjmx)

                 PRINT*, "LINEAR TIME INTERP. OF FILES ", TRIM(cur_name(cfile-1)), " AND ", TRIM(cur_name(cfile))

                 dus(:,:) = ( dus(:,:) - us(:,:) ) / frac_c
                 dvs(:,:) = ( dvs(:,:) - vs(:,:) ) / frac_c

                 us(:,:)  = us(:,:) + dus(:,:) * fac_c
                 vs(:,:)  = vs(:,:) + dvs(:,:) * fac_c
              ELSE
                 us(:,:)  =  us(:,:) + dus(:,:)
                 vs(:,:)  =  vs(:,:) + dvs(:,:)          
              ENDIF

           ELSE

              us(:,:)  =  us(:,:) + dus(:,:)
              vs(:,:)  =  vs(:,:) + dvs(:,:)

           ENDIF

        ELSE

           PRINT*, "DAILY MEAN DATA, PERSISTING WITH FILE ",TRIM(cur_name(cfile-1))
           us(:,:)  = us(:,:)
           vs(:,:)  = vs(:,:)

        ENDIF

! WAVE TIME INTERPOLATION

        IF (iwave /= 000) THEN

        IF ( time.ge.wav_time(wfile) ) THEN

            fac_w = (time - wav_time(wfile)) / delt_out

            PRINT*, "Reading forecast waves from file ", TRIM(wav_name(wfile))
            uw = READFCST_NEW(wav_dir,wav_name(wfile),lon_wav,lat_wav,"w","u",wndata,wimx,wjmx)
            vw = READFCST_NEW(wav_dir,wav_name(wfile),lon_wav,lat_wav,"w","v",wndata,wimx,wjmx)
            wfile = wfile + 1
            PRINT*, "Reading forecast waves from file ", TRIM(wav_name(wfile))
            duw = READFCST_NEW(wav_dir,wav_name(wfile),lon_wav,lat_wav,"w","u",wndata,wimx,wjmx)
            dvw = READFCST_NEW(wav_dir,wav_name(wfile),lon_wav,lat_wav,"w","v",wndata,wimx,wjmx)

            PRINT*, "LINEAR TIME INTERP. OF FILES ", TRIM(wav_name(wfile-1)), " AND ", TRIM(wav_name(wfile))

            duw(:,:) = ( duw(:,:) - uw(:,:) ) / frac_w
            dvw(:,:) = ( dvw(:,:) - vw(:,:) ) / frac_w

            uw(:,:)  = uw(:,:) + duw(:,:) * fac_w
            vw(:,:)  = vw(:,:) + dvw(:,:) * fac_w

        ELSE
            uw(:,:)  =  uw(:,:) + duw(:,:)
            vw(:,:)  =  vw(:,:) + dvw(:,:)

        ENDIF

        ENDIF

! WIND TIME INTERPOLATION

        IF ( time.ge.wnd_time(wdfile,wndrec) ) THEN

           readdata = .TRUE.
           fac_wd = (time - wnd_time(wdfile,wndrec)) / delt_out

        ENDIF
        IF (readdata) THEN
           readdata = .FALSE.
           PRINT*, "Reading forecast wind from file ", wnd_name(wdfile), " record ", wndrec
           wndx = READMET(wnd_dir,wnd_name(wdfile),wndrec,lon_wnd,lat_wnd,"x",wdndata,wdimx,wdjmx)
           wndy = READMET(wnd_dir,wnd_name(wdfile),wndrec,lon_wnd,lat_wnd,"y",wdndata,wdimx,wdjmx)
           IF (wndrec.eq.nrec_wd) THEN
              wndrec = 1
              wdfile = wdfile + 1
           ELSE
              wndrec = wndrec + 1
           ENDIF
           PRINT*, "Reading forecast wind from file ", wnd_name(wdfile), " record ", wndrec
           dwndx = READMET(wnd_dir,wnd_name(wdfile),wndrec,lon_wnd,lat_wnd,"x",wdndata,wdimx,wdjmx)
           dwndy = READMET(wnd_dir,wnd_name(wdfile),wndrec,lon_wnd,lat_wnd,"y",wdndata,wdimx,wdjmx)

           dwndx(:,:) = ( dwndx(:,:) - wndx(:,:) ) / frac_wd
           dwndy(:,:) = ( dwndy(:,:) - wndy(:,:) ) / frac_wd

           wndx(:,:)  = wndx(:,:) + dwndx(:,:) * fac_wd
           wndy(:,:)  = wndy(:,:) + dwndy(:,:) * fac_wd

        ELSE
           wndx(:,:)  =  wndx(:,:) + dwndx(:,:)
           wndy(:,:)  =  wndy(:,:) + dwndy(:,:)

        ENDIF


        PRINT*,"---------------------------------------------------------"

        uhr(:,:,t+1) = us(:,:)
        vhr(:,:,t+1) = vs(:,:)

        IF (iwave /= 000) THEN
           uwhr(:,:,t+1) = uw(:,:)
           vwhr(:,:,t+1) = vw(:,:)
        ENDIF

        wndxhr(:,:,t+1) = wndx(:,:)
        wndyhr(:,:,t+1) = wndy(:,:)

!***************************************************************************************************
! END TIME LOOP
      ENDDO

      CLOSE(21)
 
      DEALLOCATE( cst_seg )
      DEALLOCATE( cur_name )
      DEALLOCATE( cur_time )
      IF (iwave /= 000) THEN
         DEALLOCATE( wav_name )
         DEALLOCATE( wav_time )
         DEALLOCATE( uw )
         DEALLOCATE( vw )
         DEALLOCATE( duw )
         DEALLOCATE( dvw )
      ENDIF
      DEALLOCATE( wnd_time )
      DEALLOCATE( wnd_name )
      DEALLOCATE( us )
      DEALLOCATE( vs )
      DEALLOCATE( dus )
      DEALLOCATE( dvs )
      DEALLOCATE( wndx )
      DEALLOCATE( wndy )
      DEALLOCATE( dwndx )
      DEALLOCATE( dwndy )
! ------------------------------------------------------------------------

      print*, coast_map

! SPATIAL INTERPOLATION of FORCING FIELDS

      ALLOCATE( wndx_int(numspills,ntime) )
      ALLOCATE( wndy_int(numspills,ntime) )

      DO t=1,ntime

         PRINT*, ""
         PRINT*, "-----------------------------------------------"

         wndx_int(:,t) = POINT_INTERP(wdimx,lon_wnd,wdjmx,lat_wnd,wndxhr(:,:,t),numspills,cog_lon(t,:),cog_lat(t,:))
         wndy_int(:,t) = POINT_INTERP(wdimx,lon_wnd,wdjmx,lat_wnd,wndyhr(:,:,t),numspills,cog_lon(t,:),cog_lat(t,:))

         PRINT*,"POINT INTERPOLATION DONE FOR RECORD ", t
      ENDDO

      IF ( nseg /= 0 ) THEN
         IF (iwave /= 000) THEN
            print*, "nseg /=0 and iwave /=0"
            CALL WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,numspills,ntime,&
                              map_lon,map_lat,lon_cur,lat_cur,cog_lon,cog_lat,time_vec,surf_map,&
                              displ_map,uhr,vhr,wndx_int,wndy_int,lon_wav,lat_wav,uwhr,vwhr,&
                              coast_map)
         ELSE
            wimx = 0
            wjmx = 0
            print*, "nseg /=0 and iwave == 0"
            CALL WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,numspills,ntime,&
                              map_lon,map_lat,lon_cur,lat_cur,cog_lon,cog_lat,&
                              time_vec,surf_map,displ_map,uhr,vhr,wndx_int,wndy_int,&
                              coast_map=coast_map)
         ENDIF
      ELSE
         IF (iwave /= 000) THEN
            print*, "nseg ==0 and iwave /=0"
            CALL WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,numspills,ntime,&
                              map_lon,map_lat,lon_cur,lat_cur,cog_lon,cog_lat,time_vec,surf_map,&
                              displ_map,uhr,vhr,wndx_int,wndy_int,lon_wav,lat_wav,uwhr,vwhr)
         ELSE
            wimx = 0
            wjmx = 0
            print*, "nseg =0 and iwave =0"
            CALL WRITE_OUTPUT(cimx,cjmx,wimx,wjmx,map_lon_n,map_lat_n,nseg,numspills,ntime,&
                              map_lon,map_lat,lon_cur,lat_cur,cog_lon,cog_lat,&
                              time_vec,surf_map,displ_map,uhr,vhr,wndx_int,wndy_int)
         ENDIF
      ENDIF


END PROGRAM create_netcdf
