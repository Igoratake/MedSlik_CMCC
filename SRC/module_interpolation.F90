MODULE MODULE_INTERPOLATION

CONTAINS

! =======================================================================
     FUNCTION BINARYSEARCH(mode, length, array, value, delta)
!
! Given an array and a value, returns the index of the element that
! is closest to, but less than, the given value.
! Uses a binary search algorithm.
! "delta" is the tolerance used to determine if two values are equal
! if ( abs(x1 - x2) <= delta) then
!     assume x1 = x2
! endif
!=========================================================================

       IMPLICIT NONE
       INTEGER,INTENT(in)                   :: mode, length
       REAL,DIMENSION(length),INTENT(in)    :: array
       REAL,INTENT(in)                      :: value
       REAL,INTENT(in),OPTIONAL             :: delta
       INTEGER                              :: binarysearch
       INTEGER                              :: left, middle, right
       REAL                                 :: d

       IF ( present(delta) .EQV. .true.) THEN
          d = delta
       ELSE
          d = 1e-9
       ENDIF


       left = 1
       right = length

       DO
          IF ( left > right ) THEN
             EXIT
          ENDIF

          middle = nint( (left+right) / 2.0 )

          IF ( abs(array(middle) - value) <= d ) THEN
                  binarySearch = middle
                  RETURN
          ELSE IF ( array(middle) > value ) THEN
                  right = middle - 1
          ELSE
                  left = middle + 1
          ENDIF

       ENDDO
       
       IF (mode == 1) THEN
           binarySearch = -999
       ELSE
           binarySearch = right
       ENDIF
          
     END FUNCTION binarysearch

!=====================================================================================================
     FUNCTION POINT_INTERP(lon_in_len, lon_in_array, lat_in_len, lat_in_array,&
                           input_field, npoint, lon_out_array, lat_out_array)
!
! This function uses bilinear interpolation to estimate the value of a function input_field 
! at grid (lon_out_array,lat_out_array). Input_field is assumed to be sampled on a regular grid, 
! with the grid lon values specified by lon_in_array and the grid lat values specified by lat_in_array
!
!======================================================================================================

       IMPLICIT NONE

       INTEGER,INTENT(in)                                 :: lon_in_len, lat_in_len
       REAL,DIMENSION(lon_in_len),INTENT(in)              :: lon_in_array
       REAL,DIMENSION(lat_in_len),INTENT(in)              :: lat_in_array
       REAL,DIMENSION(lon_in_len,lat_in_len),INTENT(in)   :: input_field

       INTEGER,INTENT(in)                                 :: npoint
       REAL,DIMENSION(npoint),INTENT(in)                  :: lon_out_array
       REAL,DIMENSION(npoint),INTENT(in)                  :: lat_out_array
       REAL,DIMENSION(npoint)                             :: POINT_INTERP

       REAL                                               :: denom, x, y, x1, x2, y1, y2
       REAL                                               :: q1, q2, q3, q4
       INTEGER                                            :: i,j,m,n

       DO i = 1,npoint

          x = lon_out_array(i)

          m = binarysearch(0, lon_in_len, lon_in_array, x)

          y = lat_out_array(i)

          n = binarysearch(0, lat_in_len, lat_in_array, y)

          IF ( ( (m >= 1).AND.(m < lon_in_len) ).AND.( (n >= 1).AND.(n < lat_in_len) ) ) THEN

             x1 = lon_in_array(m)
             x2 = lon_in_array(m+1)

             q1 = input_field(m,n)
             q2 = input_field(m+1,n)

             y1 = lat_in_array(n)
             y2 = lat_in_array(n+1)

             q3 = input_field(m,n+1)
             q4 = input_field(m+1,n+1)

             denom = (x2 - x1)*(y2 - y1)

             POINT_INTERP(i) = ( ( q1 * (x2-x) + q2 * (x-x1) ) * (y2-y) + &
                                 ( q3 * (x2-x) + q4 * (x-x1) ) * (y-y1) )&
                                 / denom

          ELSE

             POINT_INTERP(i) = 0.

          ENDIF

       ENDDO


     END FUNCTION POINT_INTERP

!=======================================================================================================

END MODULE MODULE_INTERPOLATION
