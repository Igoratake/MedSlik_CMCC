      PROGRAM lat_lon      

      IMPLICIT NONE

      CHARACTER                           :: dummy*20, filename*16
      REAL,DIMENSION(:),ALLOCATABLE       :: lat
      REAL,DIMENSION(:),ALLOCATABLE       :: lon
      REAL*8                              :: dist,Glon,Glat,dep,length
      REAL*8                              :: sum_lat, sum_lon
      REAL                                :: alon1,alon2,alat1,alat2
      REAL,PARAMETER                      :: pi = 4.*datan(1.d0)
      REAL,PARAMETER                      :: degrad = 180.d0/pi
      INTEGER                             :: i, j, npoint, nslicks      

      CALL getarg(1,filename)
      
      open(102,file="medslik5.inp",status='old')
      do i=1,8
            read(102,*) dummy
      enddo
      read(102,*) length
      close(102)

      open(100,file=filename,status='old')

      if ( filename == "medslik_plgs.inp" ) then     
         do i=1,3
            read(100,*) dummy
         enddo
         read(100,*) nslicks
         do i=1,nslicks
            read(100,*) dummy
         enddo  
         read(100,*) npoint
         read(100,*) dummy

         allocate( lat(npoint) )
         allocate( lon(npoint) )

         do i=1,npoint      
            read(100,*) lat(i), lon(i)
         enddo
      endif

      if ( filename == "medslik_pnts.inp" ) then
         read(100,*) npoint

         allocate( lat(npoint) )
         allocate( lon(npoint) )

         do i=1,npoint
            do j=1,4
               read(100,*) dummy
            enddo
            read(100,*) lat(i)
            read(100,*) lon(i)
            read(100,*) dummy
         enddo
      endif

      close(100)

      dist = 1.5d0*length
      dist = dist / 60.d0

      sum_lat=0.0
      sum_lon=0.0

      do i=1,npoint
         sum_lat=sum_lat + lat(i)
         sum_lon=sum_lon + lon(i)
      enddo

      Glat=sum_lat/npoint
      Glon=sum_lon/npoint

      dep = dist / dcos(Glat/degrad)

      alon1 = Glon - dep
      alon2 = Glon + dep
      alat1 = Glat - dist
      alat2 = Glat + dist

      open(101,file='medslik1.tmp')

      read(101,*) dummy

      write(101,*) alon1,'   ',alon2,'        ','Min & Max Longitudes'
      write(101,*) alat1,'   ',alat2,'        ','Min & Max Latitudes'

      close(101)
           
      END PROGRAM lat_lon
