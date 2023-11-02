c-----------------------------------------------------------------------------------
c  MEDSLIK-II_1.01 
c  oil spill fate and transport model 
c-----------------------------------------------------------------------------------
c  Extract_II.for
c  This routine reads winds and currents from 
c  meteo-oceanogrpahic model output (NetCDF files)
c-----------------------------------------------------------------------------------
c  Copyright (C) <2012>
c  This program was originally written
c  by Robin Lardner and George Zodiatis.
c  Subsequent additions and modifications
c  have been made by Michela De Dominicis. 
c----------------------------------------------------------------------------------
c  The development of the MEDSLIK-II model is supported by a formal agreement
c  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
c  signed by the following institutions:
c  INGV - Istituto Nazionale di Geofisica e Vulcanologia
c  OC-UCY - Oceanography Center at the University of Cyprus
c  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per 
c  lo Studio dell’Ambiente Marino Costiero
c  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
c 
c  This program is free software: you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation, either version 3 of the License, or
c  any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c  You should have received a copy of the GNU General Public License
c  along with this program.  If not, see <http://www.gnu.org/licenses/>.
c-----------------------------------------------------------------------------------

      character regn*4, indate(30)*8,indate_wind(30)*8, fc_dir*120
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       numfiles_wind,indate_wind,iviod, icurrents,fc_dir
      integer len_dir
c--------------------------------------------------------------------
c     read subregion limits & file dates. Adjust limits to lie on OPA grid
c--------------------------------------------------------------------
      call getarg(1,fc_dir)
       
      len_dir=120
      do while(fc_dir(len_dir:len_dir).eq.' ')
      len_dir=len_dir-1
      enddo


      open(1,file='medslik.tmp')
	read(1,*) regn, icurrents, iwind
	read(1,*) alon1,alon2
	read(1,*) alat1,alat2
	read(1,*) numfiles
	do n=1,numfiles+1
	  read(1,'(a8)') indate(n)

	enddo
	read(1,*) numfiles_wind
	do n=1,numfiles_wind
	  read(1,'(a8)') indate_wind(n)

	enddo
	read(1,*) iviod
	close(1)
    
      open(99,file='Extract.log')
      if(icurrents.eq.10) call ExtractOPA(fc_dir,len_dir)
      if(icurrents.eq.11) call ExtractADRI24(fc_dir,len_dir)
      if(icurrents.eq.12) call ExtractSICI24(fc_dir,len_dir)
      if(icurrents.eq.13) call ExtractTYRR24(fc_dir,len_dir)
      
      if(icurrents.eq.70) call ExtractOPA_1hr(fc_dir,len_dir)
      if(icurrents.eq.71) call ExtractSICI(fc_dir,len_dir)
      if(icurrents.eq.72) call ExtractADRI(fc_dir,len_dir)
      if(icurrents.eq.73) call ExtractTYRR(fc_dir,len_dir)
      if(icurrents.eq.74) call ExtractRELO(fc_dir,len_dir)
      if(icurrents.eq.75) call ExtractWESTMED(fc_dir,len_dir)

      if(iwind.eq.6)      call ExtractECMWF(fc_dir,len_dir)
      if(iwind.eq.25)     call ExtractECMWF25(fc_dir,len_dir)
      stop
	end
c********************************************************************
c     Extract medslik files from OPA data
c********************************************************************
      subroutine ExtractOPA(fc_dir,len_dir)
      
      parameter(imx=677, jmx=253, kmx=10)
      
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &           ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx)
      character indate(30)*8, prdate*16, outfile*30,
     &          infile*120,heads*150, empty*80, 
     &          regn*4, fc_dir*120
      integer len_dir
      logical ex
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c-------------------------------------------------------------------
c     OPA horizontal grid
c--------------------------------------------------------------------

      oplon0 = -6.0
      oplat0 = 30.25
      op_dlon = 1./16.
      op_dlat = 1./16.

      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
      do 60 n=1,numfiles
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//indate(n)(7:8)//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate

        outfile = 'fcst_data/OPA/'//'medf'//indate(n)//'.opa'
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read OPA data files
c--------------------------------------------------------------------
      ! Open *_U.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/OPA/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idU, start, count, u)
      call handle_err(Status)

      Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/OPA/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idV, start, count, v)
      call handle_err(Status)
 
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/OPA/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idT, start, count, ts)
      call handle_err(Status)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)

c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile = 'fcst_data/OPA/'//'medf'//indate(n)//'.opa'
      open(20,file = outfile)
        write(20,*) 'MFS/OPA forecast data for '//prdate 
        write(20,*) 'Subregion of the Mediterranean with limits:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  '       u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
 
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

	    blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then     
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
                                   
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
          endif
        enddo
        enddo
        
        close(20)
   60 continue
        
      return
      end   
c********************************************************************
c     Extract medslik files from OPA data 1hr
c********************************************************************
      subroutine ExtractOPA_1hr(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=677, jmx=253, kmx=10)
      
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx),v(imx,jmx,kmx),
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
      character indate(30)*8, prdate*16, outfile*40,infile*120,                
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120
      logical ex
      integer t,kount,nore, len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     OPA horizontal grid
c--------------------------------------------------------------------
      oplon0 = -6.0
      oplat0 = 30.25
      op_dlon = 1./16.
      op_dlat = 1./16.
      
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo

       
	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
      print *, 'DAY=', n 
c--------------------------------------------------------------------
c     read OPA data files
c--------------------------------------------------------------------
    
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/O1h/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      

      
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/O1h/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/O1h/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
         
       do t=1,ktmx

          if(t.le.12) then
          kount=n
	  nore=t+12
	  write(ore,'(i2)') nore
	  ora=ore(1:2)
	  endif 
          if(t.gt.12.and.t.lt.22) then
	  kount=n+1
	  nore=t-12
	  write(ore,'(i2)') nore
	  ora='0'//ore(2:2)
	  endif
	  if(t.ge.22.and.t.le.24) then
	  kount=n+1
	  nore=t-12
	  write(ore,'(i2)') nore
	  ora=ore(1:2)
	  endif

	  prdate = indate(kount)(5:6)//'/'//indate(kount)(3:4)//'/20'//
     &                  indate(kount)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/O1h/'//'medf'//indate(kount)(1:6)//ora(1:2)//'.opa'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read OPA data files
c--------------------------------------------------------------------

	u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	
        call extrap3d(ts, i1, i2, j1, j2, imx, jmx,kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   

          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/O1h/'//'medf'//indate(kount)(1:6)//ora(1:2)//'.opa'
      open(20,file = outfile)
        write(20,*) 'MFS/OPA forecast data for '//prdate 
        write(20,*) 'Subregion of the Mediterranean with limits:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)

            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
c********************************************************************
c     Extract medslik files from Sicily data 1hr
c********************************************************************
      subroutine ExtractSICI(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=257, jmx=273, kmx=10)
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx), 
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),             
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
      character indate(30)*8, prdate*16, outfile*40,infile*120,
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120
      logical ex
      integer t, len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     SICI horizontal grid
c--------------------------------------------------------------------
        
      oplon0=8.9844
      oplat0=30.984
      op_dlon = 1./32.
      op_dlat = 1./32.

      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read SICILY data files
c--------------------------------------------------------------------
      
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/S1h/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      
      
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'U', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/S1h/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'V', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/S1h/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'TEM', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
         
      
          do t=1,ktmx
          if(t.lt.10) then
           write(ore,'(i2)') t
	   ora='0'//ore(2:2)
	 else
	 write(ora,'(i2)') t
	 endif 
       
      
	   
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/S1h/'//'sici'//indate(n)(1:6)//ora(1:2)//'.sic'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read SICILY data files
c--------------------------------------------------------------------
      

	u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       


        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/S1h/'//'sici'//indate(n)(1:6)//ora(1:2)//'.sic'
      open(20,file = outfile)
        write(20,*) 'SCRM forecast data for '//prdate 
        write(20,*) 'Subregion of the Sicily Strait:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
      
      
c********************************************************************
c     Extract medslik files from Sicily data 24hr
c********************************************************************
      subroutine ExtractSICI24(fc_dir,len_dir)
      
      parameter(imx=257, jmx=273, kmx=10)

      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx)              
      character indate(30)*8, prdate*16, outfile*40,
     &          infile*120,heads*150, empty*80,
     &          regn*4,fc_dir*120
      logical ex
      integer len_dir
	common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
	data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     Sicily model horizontal grid
c--------------------------------------------------------------------
 
      oplon0=8.9844
      oplat0=30.984

      op_dlon = 1./32.
      op_dlat = 1./32.
      
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read Sicily data files
c--------------------------------------------------------------------
      ! Open *_U.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/S24/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'U', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idU, start, count, u)
      call handle_err(Status)
      Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/S24/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'V', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idV, start, count, v)
      call handle_err(Status)
   
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/S24/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'TEM', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idT, start, count, ts)
      call handle_err(Status)
     
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//indate(n)(7:8)//':00'

        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate

        outfile = 'fcst_data/S24/'//'sici'//indate(n)//'.sic'
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------
        outfile = 'fcst_data/S24/'//'sici'//indate(n)//'.sic'

      open(20,file = outfile)
        write(20,*) 'Sicily Channel Regional Model
     & forecast data for '//prdate 
        write(20,*) 'Subregion of the Sicily Channel Regional Model:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo

        close(20)
   60 continue
        
      return
      end         
      
      
      
c********************************************************************
c     Extract medslik files from Adriatic data 1hr
c********************************************************************
      subroutine ExtractADRI(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=287, jmx=311, kmx=10)

      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx), 
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),             
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
    
    
      character indate(30)*8, prdate*16, outfile*40,infile*120,
     &        heads*150, empty*80, regn*4, ora*2, ore*2,fc_dir*120
      logical ex
      integer t,kount,nore,len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     Adriatic horizontal grid
c--------------------------------------------------------------------
          
      oplon0=12.2
      oplat0=39.0

      op_dlon = 0.03
      op_dlat = 0.022
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read Adriatic data files
c--------------------------------------------------------------------

      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/A1h/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      
      
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/A1h/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/A1h/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
      
	  
          do t=1,ktmx
          if(t.le.12) then
          kount=n
	  nore=t+12
	  write(ore,'(i2)') nore
	  ora=ore(1:2)
	  endif 
          if(t.gt.12.and.t.lt.22) then
	  kount=n+1
	  nore=t-12
	  write(ore,'(i2)') nore
	  ora='0'//ore(2:2)
	  endif
	  if(t.ge.22.and.t.le.24) then
	  kount=n+1
	  nore=t-12
	  write(ore,'(i2)') nore
	  ora=ore(1:2)
	  endif
     
 	   
      prdate = indate(kount)(5:6)//'/'//indate(kount)(3:4)//'/20'//
     &                  indate(kount)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/A1h/'//'adri'//indate(kount)(1:6)//ora(1:2)//'.adr'      
     
   

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read Adriatic data files
c--------------------------------------------------------------------
      

	u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/A1h/'//'adri'//indate(kount)(1:6)//ora(1:2)//'.adr'
      open(20,file = outfile)
        write(20,*) 'Adriatic forecast data for '//prdate 
        write(20,*) 'Subregion of the Adriatic Sea:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
	    us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.6)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
c********************************************************************
c     Extract medslik files from Adriatic data 24hr
c********************************************************************
      subroutine ExtractADRI24(fc_dir,len_dir)
      
      parameter(imx=287, jmx=311, kmx=10)
     
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx)              
    
      character indate(30)*8, prdate*16, outfile*40,infile*120,
     &          heads*150, empty*80, regn*4,fc_dir*120
      logical ex
      integer len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     Adriatic horizontal grid
c--------------------------------------------------------------------
        
      oplon0=12.2
      oplat0=39.0

      op_dlon = 0.03
      op_dlat = 0.022
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
      	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//indate(n)(7:8)//':00'

        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate

        outfile = 'fcst_data/A24/'//'adri'//indate(n)//'.adr'
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif   
c--------------------------------------------------------------------
c     read Adriatic data files
c--------------------------------------------------------------------
    ! Open *_U.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/A24/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idU, start, count, u)
      call handle_err(Status)

      Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/A24/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idV, start, count, v)
      call handle_err(Status)
    
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/A24/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idT, start, count, ts)
      call handle_err(Status)
       

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------
        outfile = 'fcst_data/A24/'//'adri'//indate(n)//'.adr'
      open(20,file = outfile)
        write(20,*) 'Adriatic forecast data for '//prdate 
        write(20,*) 'Subregion of the Adriatic Sea:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
c********************************************************************
c     Extract medslik files from Tyrrhenian data 1hr
c********************************************************************
      subroutine ExtractTYRR(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=360, jmx=376, kmx=10)

      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx), 
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),             
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
    
      character indate(30)*8, prdate*16, outfile*40,infile*120,                
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120
      logical ex
      integer t,len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /999.999/,      rhoa /1.19/

c--------------------------------------------------------------------
c     Tyrrhenian horizontal grid
c--------------------------------------------------------------------
        
      oplon0=8.8125
      oplat0=36.6875
      op_dlon = 1./48.
      op_dlat = 1./48.

      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read Tyrrhenian data files
c--------------------------------------------------------------------


      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/T1h/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
       
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo
       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/T1h/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/T1h/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
         
          do t=1,ktmx
          if(t.lt.10) then
           write(ore,'(i2)') t
	   ora='0'//ore(2:2)
	 else
	 write(ora,'(i2)') t
	 endif 
       
      
	   
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/T1h/'//'tyrr'//indate(n)(1:6)//ora(1:2)//'.tyr'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read Tyrrhenian data files
c--------------------------------------------------------------------
      

	u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/T1h/'//'tyrr'//indate(n)(1:6)//ora(1:2)//'.tyr'
      open(20,file = outfile)
        write(20,*) 'Tyrrhenian forecast data for '//prdate 
        write(20,*) 'Subregion of the Tyrrhenian Sea:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
c********************************************************************
c     Extract medslik files from Tyrrhenian data 24hr
c********************************************************************
 
      subroutine ExtractTYRR24(fc_dir,len_dir)
      
      parameter(imx=360, jmx=376, kmx=10)
    
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx)              
    
	character indate(30)*8, prdate*16, outfile*40,
     &          infile*120, heads*150, empty*80, 
     &          regn*4, fc_dir*120
    
        logical ex
        integer len_dir
	common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
	data udef /9999./,      rhoa /1.19/
c--------------------------------------------------------------------
c     Tyrrhenian horizontal grid
c--------------------------------------------------------------------
        
      oplon0=8.81
      oplat0=36.68
      op_dlon = 1./48.
      op_dlat = 1./48. 
      
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
      	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//indate(n)(7:8)//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile = 'fcst_data/T24/'//'tyrr'//indate(n)//'.tyr'
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            

c--------------------------------------------------------------------
c     read Tyrr data files
c--------------------------------------------------------------------
    ! Open *_U.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/T24/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idU, start, count, u)
      call handle_err(Status)

      Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/T24/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idV, start, count, v)
      call handle_err(Status)

    ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/T24/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1
      Status = nf_get_vara_real (id, idT, start, count, ts)
      call handle_err(Status)


c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------
        outfile = 'fcst_data/T24/'//'tyrr'//indate(n)//'.tyr'

      open(20,file = outfile)
        write(20,*) 'Tyrrhenian Regional Model
     & forecast data for '//prdate 
        write(20,*) 'Subregion of the Tyrrhenian Regional Model:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end         
      
      
      
c********************************************************************
c     Extract medslik files from Relocatable model output (1hr-3km)
c********************************************************************
      subroutine ExtractRELO(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=100, jmx=100, kmx=19)
c     parameter(ktmx=24, imx=130, jmx=88, kmx=19)  
      
      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx),
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),              
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
      character indate(30)*8, prdate*16, outfile*40,infile*120,                
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120,
     &          start_lat*19
      logical ex
      integer t,len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     Relocatable horizontal grid
c--------------------------------------------------------------------
c      open(1,file=fc_dir(1:len_dir)//'/fcst_data/H3k//'
c     &                                //indate(1)(1:6)//'_T.txt')
	  
c	  do i=1,200
c          read(1,'(a11)') start
c	  if(start.eq.'		nav_lat:l') then

c	  read(1,'(a20,1x,f9.5)') empty,oplat0
c	  print *,'oplat0',oplat0
c	  endif
c	  if(start.eq.'		nav_lon:l') then

c	  read(1,'(a20,1x,f9.5)') empty,oplon0
c	  print *,'oplon0',oplon0
c	  endif
c        enddo
c      close(1)
       oplon0=9.195901
       oplat0=40.91457

c      oplon0=6.998218
c      oplat0=42.14644
c      oplon0=3.629855
c      oplat0=36.55457

       op_dlon = 0.036447
c      op_dlon = 0.03708
c      op_dlon = 0.03418
      op_dlat = 0.02698

      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read Relocatable data files
c--------------------------------------------------------------------

    
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo

      
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vozocrtx', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'vomecrty', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'votemper', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
          do t=1,ktmx
          if(t.lt.10) then
           write(ore,'(i2)') t
	   ora='0'//ore(2:2)
	 else
	 write(ora,'(i2)') t
	 endif 
       
      
	   
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/H3k//'//'relo'//indate(n)(1:6)//ora(1:2)//'.rel'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read Relocatable data files
c--------------------------------------------------------------------
        u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/H3k//'//'relo'//indate(n)(1:6)//ora(1:2)//'.rel'
      open(20,file = outfile)
        write(20,*) 'Relocatable model 3 km forecast data for '//prdate 
        write(20,*) 'Subregion of the Mediterranean Sea:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.5) then     
              u10 = (u(i,j,4) * 1.559 + u(i,j,5) * 2.056 ) / 3.615 
              v10 = (v(i,j,4) * 1.559 + v(i,j,5) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.10) then     
              u30 = (u(i,j,9) * 4.164 + u(i,j,10) * 1.032 ) / 5.196 
              v30 = (v(i,j,9) * 4.164 + v(i,j,10) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.21) then   
              u120 = (u(i,j,20) * 3.459 + u(i,j,21) * 7.753 ) / 11.212 
              v120 = (v(i,j,20) * 3.459 + v(i,j,21) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
            write (1,*) i_first,i_last,j_first,j_last
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
                  
c********************************************************************
c     Extract medslik files from West Mediterranean data 1hr
c********************************************************************
      subroutine ExtractWESTMED(fc_dir,len_dir)
      
      parameter(ktmx=24, imx=432, jmx=250, kmx=10)

      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx),
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),              
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
      character indate(30)*8, prdate*16, outfile*40,infile*120,                
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120
      logical ex
      integer t,len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c     WESTMED horizontal grid
c--------------------------------------------------------------------
        
      oplon0=3
      oplat0=36.700

      op_dlon = 1./32.
      op_dlat = 1./32.
      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
        oplat(j) = oplat0 + (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
	alat1 = oplat0 + (j_first - 1) * op_dlat
	alat2 = oplat0 + (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1

	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read West Med data files
c--------------------------------------------------------------------

    
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/WME/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo

      
       ! Open *_U.nc File 
      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'U', idU)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idU, start, count, u_tmp)
       call handle_err(Status)
       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)

      ! Open *_V.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/WME/'
     &                                //indate(n)(1:6)//'_V.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'V', idV)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idV, start, count, v_tmp)
       call handle_err(Status)
       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
   
      
      ! Open *_T.nc File 
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/WME/'
     &                                //indate(n)(1:6)//'_T.nc'
      Status = nf_open(infile,nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'TEM', idT)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     
 
       do t = 1,ktmx
       start(4)=t
       Status = nf_get_vara_real (id, idT, start, count, ts_tmp)
       call handle_err(Status)
       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)
       enddo

       Status = nf_close (id)
          
         
      
          do t=1,ktmx
          if(t.lt.10) then
           write(ore,'(i2)') t
	   ora='0'//ore(2:2)
	 else
	 write(ora,'(i2)') t
	 endif 
       
      
	   
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/WME/'//'wmed'//indate(n)(1:6)//ora(1:2)//'.wme'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read West Med data files
c--------------------------------------------------------------------
      

	u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	
c	call extrap2d(taux, i1, i2, j1, j2, imx, jmx)
c        call extrap2d(tauy, i1, i2, j1, j2, imx, jmx)
        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            do k=1,kmx
              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
            enddo   
          endif
        enddo
        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/WME/'//'wmed'//indate(n)(1:6)//ora(1:2)//'.wme'
      open(20,file = outfile)
        write(20,*) 'West Med forecast data for '//prdate 
        write(20,*) 'Subregion of the West Mediterranean Sea:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)

	      km = 1
            do k=2,kmx
              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
            enddo

            if(km.ge.4) then     
              u10 = (u(i,j,3) * 1.559 + u(i,j,4) * 2.056 ) / 3.615 
              v10 = (v(i,j,3) * 1.559 + v(i,j,4) * 2.056 ) / 3.615 
            else
              u10 = u(i,j,km)
              v10 = v(i,j,km) 
            endif
            if(km.ge.9) then     
              u30 = (u(i,j,8) * 4.164 + u(i,j,9) * 1.032 ) / 5.196 
              v30 = (v(i,j,8) * 4.164 + v(i,j,9) * 1.032 ) / 5.196 
            else
              u30 = u(i,j,km)
              v30 = v(i,j,km) 
            endif
            if(km.ge.20) then   
              u120 = (u(i,j,19) * 3.459 + u(i,j,20) * 7.753 ) / 11.212 
              v120 = (v(i,j,19) * 3.459 + v(i,j,20) * 7.753 ) / 11.212
            else
              u120 = u(i,j,km)
              v120 = v(i,j,km) 
            endif
             
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   


 
                  

 
                  
c********************************************************************
c     Extract medslik files from ECMWF
c********************************************************************
      subroutine ExtractECMWF(fc_dir,len_dir)
      

       parameter(ktmx=4, imx=123, jmx=37)
       real fmis !netcdf
       parameter(fmis=0.) !netcdf
       integer start(3), count(3) !netcdf  
       integer id, idWX, idWY !netcdf
       integer Status !netcdf
       dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          wx(imx,jmx,ktmx),wy(imx,jmx,ktmx),
     &          wx_tmp(imx,jmx),wy_tmp(imx,jmx)

       character indate(30)*8, prdate*16, outfile*40, infile*120,
     &          heads*150, empty*80, regn*4, ora*2, ore*2,fc_dir*120,
     &          indate_wind(30)*8
       logical ex
       integer t,len_dir
       common regn, alon1, alon2, alat1, alat2, numfiles,indate,  
     &       numfiles_wind,indate_wind, iviod, icurrents
       data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c    ECWMF horizontal grid
c--------------------------------------------------------------------
      numfiles=numfiles_wind
      indate=indate_wind
      oplon0 = -19
      oplat0 = 48

      op_dlon = 1./2.
      op_dlat = 1./2.

      do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
c        oplat(j) = oplat0 + (j-1) * op_dlat
         oplat(j) = oplat0 - (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2 ! +2??
c	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
c	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2
        j_first = int( (oplat0 - alat2) / op_dlat ) + 1
	j_last  = int( (oplat0 - alat1) / op_dlat ) + 2 ! +2?

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
c	alat1 = oplat0 + (j_first - 1) * op_dlat
c	alat2 = oplat0 + (j_last  - 1) * op_dlat
        alat2 = oplat0 - (j_first - 1) * op_dlat
	alat1 = oplat0 - (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1
c        print *, 'i-limits   = ',i_first,i_last,imax
c	print *, 'j-limits   = ',j_first,j_last,jmax
c	print *, 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
c	print *, 'lat-limits = ',alat1,alat2,(alat2-alat1)*16
	
	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read ECMWF data files
c--------------------------------------------------------------------

      ! Open *.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/ECM/'
     &                                //indate_wind(n)(1:8)//'.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'U10M', idWX)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
   
     
      count(1) = imx
      count(2) = jmx
      count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWX, start, count, wx_tmp)
       call handle_err(Status)
       wx(1:imx,1:jmx,t) = wx_tmp(1:imx,1:jmx)
       enddo

       Status = nf_close (id)
       
       Status = nf_open(infile(1:len_file),nf_nowrite,id)
       call handle_err(Status)
       Status = nf_inq_varid (id, 'V10M', idWY)
       call handle_err(Status)
       start(1) = 1
       start(2) = 1
       start(3) = 1
        
       count(1) = imx
       count(2) = jmx
       count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWY, start, count, wy_tmp)
       call handle_err(Status)
       wy(1:imx,1:jmx,t) = wy_tmp(1:imx,1:jmx)
       enddo

       Status = nf_close (id)
    
  	   
	  prdate = indate_wind(n)(7:8)//'/'//indate(n)(5:6)//'/20'//
     &                  indate(n)(3:4)
        write(6,*) 'Writing medslik ECMWF file for date '//prdate
        write(99,*) 'Writing medslik ECMWF file for date '//prdate
        outfile =
     & 'fcst_data/ECM/'//'ecm_'//indate(n)(3:8)//'.ecm'
  
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
	
c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(wx(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
     
	
       call extrap2d(wx, i1, i2, j1, j2, imx, jmx)
       call extrap2d(wy, i1, i2, j1, j2, imx, jmx)


c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
  	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/ECM/'//'ecm_'//indate(n)(3:8)//'.ecm'
      open(20,file = outfile)
        write(20,*) 'ECMWF forecast data for '//prdate 
        write(20,*) 'Subregion of the Mediterranean with limits:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
	heads = '    lat        lon        winx0   '//
     &  '    winy0      winx6      winy6      winx12     winy12'//
     &  '      winx18     winy18'
        write(20,'(a150)') heads

	
        do i=i_first,i_last
c        do j=j_first,j_last
         do j=j_last,j_first,-1
c           if(msk(i,j).eq.1) then

	    blon = oplon(i)
            blat = oplat(j)
            
c	    wx0 = wx(i,j,1)
c	    wy0 = wy(i,j,1)
c            endif

	     
              wx00 = wx(i,j,1)
              wy00 = wy(i,j,1)
            
	      
              wx06 = wx(i,j,2)
              wy06 = wy(i,j,2)
            
	     
              wx12 = wx(i,j,3)
              wy12 = wy(i,j,3)
          
	 
              wx18 = wx(i,j,4)
              wy18 = wy(i,j,4)
           
           
             
            write(20,'(13f11.4)') blat,blon,wx00,wy00,wx06,wy06
     &      ,wx12,wy12,wx18,wy18
          
        
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
     


c********************************************************************
c     Extract medslik files from ECMWF
c********************************************************************
      subroutine ExtractECMWF25(fc_dir,len_dir)
      

       parameter(ktmx=4, imx=245, jmx=73)

       real fmis !netcdf
       parameter(fmis=0.) !netcdf
       integer start(3), count(3) !netcdf  
       integer id, idWX, idWY !netcdf
       integer Status !netcdf
       dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          wx(imx,jmx,ktmx),wy(imx,jmx,ktmx),
     &          wx_tmp(imx,jmx),wy_tmp(imx,jmx)

       character indate(30)*8, prdate*16, outfile*40, infile*120,
     &          heads*150, empty*80, regn*4, ora*2, ore*2,fc_dir*120,
     &          indate_wind(30)*8
       logical ex
       integer t,len_dir
       common regn, alon1, alon2, alat1, alat2, numfiles,indate,  
     &       numfiles_wind,indate_wind, iviod, icurrents
       data udef /9999./,      rhoa /1.19/

c--------------------------------------------------------------------
c    ECWMF horizontal grid
c--------------------------------------------------------------------
      numfiles=numfiles_wind
      indate=indate_wind
      oplon0 = -19
c     oplat0 = 30
      oplat0 = 48
      op_dlon = 1./4.
      op_dlat = 1./4.
       do i=1,imx
        oplon(i) = oplon0 + (i-1) * op_dlon
      enddo
      do j=1,jmx
c        oplat(j) = oplat0 + (j-1) * op_dlat
         oplat(j) = oplat0 - (j-1) * op_dlat
      enddo


	i_first = int( (alon1 - oplon0) / op_dlon ) + 1
	i_last  = int( (alon2 - oplon0) / op_dlon ) + 2 ! +2??
c	j_first = int( (alat1 - oplat0) / op_dlat ) + 1
c	j_last  = int( (alat2 - oplat0) / op_dlat ) + 2
        j_first = int( (oplat0 - alat2) / op_dlat ) + 1
	j_last  = int( (oplat0 - alat1) / op_dlat ) + 2 ! +2?

	if(i_first.lt.1) i_first = 1
	if(i_last.gt.imx) i_last = imx
	if(j_first.lt.1) j_first = 1
	if(j_last.gt.jmx) j_last = jmx

	alon1 = oplon0 + (i_first - 1) * op_dlon
	alon2 = oplon0 + (i_last  - 1) * op_dlon
c	alat1 = oplat0 + (j_first - 1) * op_dlat
c	alat2 = oplat0 + (j_last  - 1) * op_dlat
        alat2 = oplat0 - (j_first - 1) * op_dlat
	alat1 = oplat0 - (j_last  - 1) * op_dlat

	imax = i_last - i_first + 1
	jmax = j_last - j_first + 1
c        print *, 'i-limits   = ',i_first,i_last,imax
c	print *, 'j-limits   = ',j_first,j_last,jmax
c	print *, 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
c	print *, 'lat-limits = ',alat1,alat2,(alat2-alat1)*16
	
	write(99,*) 'i-limits   = ',i_first,i_last,imax
	write(99,*) 'j-limits   = ',j_first,j_last,jmax
	write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
	write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16

c--------------------------------------------------------------------
c     Begin main loop. 
c     First check if the files already exist for the current subregion
c--------------------------------------------------------------------
    
      do 60 n=1,numfiles
c--------------------------------------------------------------------
c     read ECMWF data files
c--------------------------------------------------------------------

      ! Open *.nc File 
      Status = 0
      
      infile=fc_dir(1:len_dir)//'/fcst_data/E25/'
     &                                //indate_wind(n)(1:8)//'.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
      call handle_err(Status)
      Status = nf_inq_varid (id, 'U10M', idWX)
      call handle_err(Status)
      start(1) = 1
      start(2) = 1
      start(3) = 1
   
     
      count(1) = imx
      count(2) = jmx
      count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWX, start, count, wx_tmp)
       call handle_err(Status)
       wx(1:imx,1:jmx,t) = wx_tmp(1:imx,1:jmx)
       enddo

       Status = nf_close (id)
       
       Status = nf_open(infile(1:len_file),nf_nowrite,id)
       call handle_err(Status)
       Status = nf_inq_varid (id, 'V10M', idWY)
       call handle_err(Status)
       start(1) = 1
       start(2) = 1
       start(3) = 1
        
       count(1) = imx
       count(2) = jmx
       count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWY, start, count, wy_tmp)
       call handle_err(Status)
       wy(1:imx,1:jmx,t) = wy_tmp(1:imx,1:jmx)
       enddo

       Status = nf_close (id)
    
  	   
	  prdate = indate_wind(n)(7:8)//'/'//indate(n)(5:6)//'/20'//
     &                  indate(n)(3:4)
        write(6,*) 'Writing medslik ECMWF file for date '//prdate
        write(99,*) 'Writing medslik ECMWF file for date '//prdate
        outfile =
     & 'fcst_data/E25/'//'ecm_'//indate(n)(3:8)//'.ecm'
  
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
	
c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(wx(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
     
	
       call extrap2d(wx, i1, i2, j1, j2, imx, jmx)
       call extrap2d(wy, i1, i2, j1, j2, imx, jmx)


c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
  	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/E25/'//'ecm_'//indate(n)(3:8)//'.ecm'
      open(20,file = outfile)
        write(20,*) 'ECMWF forecast data for '//prdate 
        write(20,*) 'Subregion of the Mediterranean with limits:' 
        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
	heads = '    lat        lon        winx0   '//
     &  '    winy0      winx6      winy6      winx12     winy12'//
     &  '      winx18     winy18'
        write(20,'(a150)') heads

	
        do i=i_first,i_last
c        do j=j_first,j_last
         do j=j_last,j_first,-1
c           if(msk(i,j).eq.1) then

	    blon = oplon(i)
            blat = oplat(j)
            
c	    wx0 = wx(i,j,1)
c	    wy0 = wy(i,j,1)
c            endif

	     
              wx00 = wx(i,j,1)
              wy00 = wy(i,j,1)
            
	      
              wx06 = wx(i,j,2)
              wy06 = wy(i,j,2)
            
	     
              wx12 = wx(i,j,3)
              wy12 = wy(i,j,3)
          
	 
              wx18 = wx(i,j,4)
              wy18 = wy(i,j,4)
           
           
             
            write(20,'(13f11.4)') blat,blon,wx00,wy00,wx06,wy06
     &      ,wx12,wy12,wx18,wy18
          
        
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
     




c************************************************************************
c      Extrapolation of 2-D fields over land points
c------------------------------------------------------------------------

      subroutine extrap2d(data, i_first, i_last, j_first, j_last, 
     &                          imx, jmx)

      real data(imx,jmx)
      
      data udef /999.999/

	ngridpts = (i_last - i_first + 1) * (j_last - j_first + 1)

      do iter=1,100
        knt=0
        do j = j_first, j_last
        do i = i_first, i_last

          if(data(i,j). ge. udef) then
             knt = knt + 1
               im1=i-1
               ip1=i+1
               jm1=j-1
               jp1=j+1
               if(im1.lt.i_first) im1 = i_first
               if(ip1.gt.i_last)  ip1 = i_last
               if(jm1.lt.j_first) jm1 = j_first
               if(jp1.gt.j_last)  jp1 = j_last

               datan=0.
               jcn=0
               do jj=jm1,jp1
               do ii=im1,ip1
                    if(data(ii,jj). lt. udef) then
                      datan = datan + data(ii,jj)
                      jcn = jcn + 1
                   end if
               enddo 
               enddo 

               if(jcn. gt. 2) then
                 data(i,j) = datan / float(jcn)
               end if

          end if

        enddo 
        enddo
c        write(99,*) iter,knt
        if(knt.eq.0.or.knt.eq.ngridpts) return 
        
      enddo
                    
      return
      end      

c************************************************************************
c      Extrapolation of 3-D fields over land points, 7 levels only
c------------------------------------------------------------------------

      subroutine extrap3d (data, i_first, i_last, j_first, j_last, 
     &                           imx, jmx, kmx)


      real data(imx,jmx,kmx), tmp(imx,jmx)
      integer kd(7)

      data udef /999.999/
      data kd /1,3,4,8,9,19,20/
c
      do k=1,kmx
c        k = kd(n)

        do i=1,imx
        do j=1,jmx
          tmp(i,j) = data(i,j,k)
        enddo
        enddo
        
        call extrap2d(tmp, i_first, i_last, j_first, j_last, imx, jmx)

        do i=1,imx
        do j=1,jmx
          data(i,j,k) = tmp(i,j)
        enddo
        enddo
        
      enddo   !k

      return
      end      

c****************************************************************************
c****************************************************************************
c****************************************************************************
      SUBROUTINE HANDLE_ERR(Status)
      include 'netcdf.inc'
c      include '/usr/local/include/netcdf.inc'
      INTEGER Status
      IF (Status .NE. NF_NOERR) THEN
        PRINT *, NF_STRERROR(Status)
        STOP 'Stopped'
      ENDIF
      END
