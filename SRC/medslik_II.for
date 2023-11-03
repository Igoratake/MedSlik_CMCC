!================================================================================
!  MEDSLIK-II_1.02                                                               |
!                                                                                |
!  Oil spill fate and transport numerical model                                  |
!--------------------------------------------------------------------------------|
!  medslik_II.F90                                                                |
!                                                                                |
!  This routine reads winds and currents from                                    |
!  meteo-oceanogrpahic model output (NetCDF files)                               |
!  and write them into MEDSLIK-II ascii formatted input files                    |
!                                                                                |
!--------------------------------------------------------------------------------|
!                                                                                |
!  Copyright (C) <2012>                                                          |
!                                                                                |
!  This program was originally written by Robin Lardner and George  Zodiatis.    |
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
!                                                                                |
!  This program is distributed in the hope that it will be useful,               |
!  but WITHOUT ANY WARRANTY; without even the implied warranty of                |
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 |
!  GNU General Public License for more details.                                  |
!  You should have received a copy of the GNU General Public License             |
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.         |
!================================================================================

      PROGRAM medslik_II

      implicit real*8(a-h,o-z)
      character dummy*80


      open(1,file='medslik5.par')
      do k=1,43
         read(1,'(a80)') dummy
      enddo

      read(1,'(a80)') dummy
      if(dummy(1:6).eq.'Comput') then
         read(1,*) nstph
         read(1,*) npl
         delt = 1.d0 / dfloat(nstph) 
      else
         delt = 0.5d0
         npl = 2000
      endif
      close(1)     

      print*, "delta t: ", delt     
      call main(delt,npl)
      
      END PROGRAM medslik_II
    
c======================================================================
c                        BEGIN THE MAIN PROGRAM
c======================================================================
      
      subroutine main(delt, npl)
 

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=200000,msp=1200,
     &          npol=10)

c----------------------------------------------------------------------
c     Dimension declarations for PHYSICAL-ENVIRONMENTAL variables
c----------------------------------------------------------------------

c    Bathymetry variables

      dimension 
     &      h(mm,nm),itype(mm,nm)

c    Current and wind advective velocity variables

      dimension
     &      uadv(mm,nm),vadv(mm,nm),uadv1(mm,nm),vadv1(mm,nm),
     &      uadv2(mm,nm),vadv2(mm,nm),wdrftx(mm,nm),wdrfty(mm,nm),
     &      wdrftx1(mm,nm),wdrfty1(mm,nm),wdrftx2(mm,nm),
     &      wdrfty2(mm,nm)


c    Jonswap spectra variables

      dimension
     &      wvel_vec(720),wvel_mean(720),wdir_vec(720),wdir_mean(720),
     &      freq_sp(700),fang(700),erre(700),sp_exp1(700), spectra(700),
     &      wave_num(700),stoke_sp(700),hwave_d(700),stoke_d(700),
     &      sp_exp2(700)


c    Corrections to oil drift velocity from slick observation variables

      dimension 
     &      crntim(5),crnu(5),crnv(5)

c    Forcing input data related variables

      dimension
     &      ifcstfn(720),fcsttim(720),iwfcstfn(30),wfcsttim(30,24),
     &      iwavfcstfn(720), wavfcsttim(720)
 
c    Booms related variables

      dimension 
     &      bmtim(20),bmlat1(20),bmlon1(20),bmlat2(20),bmlon2(20),
     &      bmx1(20),bmy1(20),bmx2(20),bmy2(20),ibmeff(20)

c----------------------------------------------------------------------
c     Dimension declarations for oil slick properties variables
c----------------------------------------------------------------------

      dimension 
     &      avgdep(npl,npl),sum_avgdep(npl,npl),n_disp(npl,npl),
     &      zdispl(npc),c1p(npc),c2p(npc),poll(npl,npl),
     &      px(npc),py(npc),pz(npc),
     &      is(npc),ib(npc),itmp(npc),bd1(nss,4),ibd(nss),
     &      seg(nss,4),alngt(nss),segx(nss),segy(nss),ns0(nss),
     &      seep(nss),prel(nss),sfrac(nss),vcst(nss),
     &      den(msp),vis(msp),visem(msp),tre(msp),c1ms(msp),c2ms(msp),
     &      atn(msp),atk(msp),ato(msp),ttk(msp),tto(msp),
     &      vtn(msp),vtk(msp),vto(msp),xcl(msp),xcs(msp),
     &      vtne(msp),vtke(msp),vte(msp),vtnd(msp),vtkd(msp),vtd(msp),
     &      ftk(msp),ftn(msp),fw(msp),pcte(msp),pctd(msp)


c----------------------------------------------------------------------
c     Dimension declarations for multiple slick properties variables
c----------------------------------------------------------------------

      dimension
     &      iddv(20),immv(20),iyrv(20),istartv(20),ispillv(20),
     &      ce1_vec(20),tstartv(20),tspillv(20),splatv(20),
     &      splonv(20),splrtev(20),x0v(20),y0v(20),totbblv(20),
     &      nppmsv(20),nstrel(20),vmsplv(20),bblpmsv(20),
     &      rtk0v(20),rtn0v(20),xavg_multi(20),yavg_multi(20),
     &      yavg_lat_multi(20),xavg_lon_multi(20),tkv(20)

c----------------------------------------------------------------------
c     COMMON DECLARATION
c----------------------------------------------------------------------

c     Geographical and geometrical variables

      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk2/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

c     Correction from slick observations variables

      common /encr/ aloncorr,iaencr,ibencr,icencr

c     Input forcing varables

      common /wind/ winx(mm,nm),winy(mm,nm),dwinx(mm,nm),dwiny(mm,nm)
      common /curr/ sst(mm,nm),usrf(mm,nm),vsrf(mm,nm),u10(mm,nm),
     &              v10(mm,nm),u30(mm,nm),v30(mm,nm),u120(mm,nm),
     &              v120(mm,nm),dsst(mm,nm),dusrf(mm,nm),dvsrf(mm,nm),
     &              du10(mm,nm),dv10(mm,nm),du30(mm,nm),dv30(mm,nm),
     &              du120(mm,nm),dv120(mm,nm)
      common /wave/ stoku(mm,nm),stokv(mm,nm),dstoku(mm,nm),
     &              dstokv(mm,nm)

c     Runge-Kutta 4th order scheme variables

      common /rng1/ ucur1(mm,nm,4),vcur1(mm,nm,4),usto1(mm,nm),
     &              vsto1(mm,nm),ducr1(mm,nm,4),dvcr1(mm,nm,4),
     &              dust1(mm,nm),dvst1(mm,nm),wx1(mm,nm),wy1(mm,nm),
     &              dwx1(mm,nm),dwy1(mm,nm),irngc1,irngw1,irngwd1,
     &              isubc1,isubw1,isubwd1
      common /rng2/ ucur2(mm,nm,4),vcur2(mm,nm,4),usto2(mm,nm),
     &              vsto2(mm,nm),ducr2(mm,nm,4),dvcr2(mm,nm,4),
     &              dust2(mm,nm),dvst2(mm,nm),wx2(mm,nm),wy2(mm,nm),
     &              dwx2(mm,nm),dwy2(mm,nm),irngc2,irngw2,irngwd2,
     &              isubc2,isubw2,isubwd2
      common /rng3/ urng1(mm,nm),vrng1(mm,nm),urng2(mm,nm),
     &              vrng2(mm,nm)

c     Oil & spill/slick properties variables

      common /spill/ idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /evap/  ce,vappr,fmaxe
      common /disp/  cd1,cd3,cd4,cd5,vl,vs1,um,stk,stn,fmaxd
      common /emul/  cm1,cm2,cm3,visemx
      common /sprd/  cs1,cs2,cs3
      common /phys/  deno,denk,cdt,viso,visk,cvt,tvk0,denw,tk,tk0

c----------------------------------------------------------------------
c     CHARACTER declarations
c---------------------------------------------------------------------- 

      character empty*80, a0(3)*4, a1*1, a2*2, a3*3, a4*4, 
     &          ay*2, am*2, ad*2, ah*2, d24*3, d06*3, vlunit*4, pref*4,
     &          nore*2, nora*2,ora*2,ore*2,indate(30)*8,
     &          dummy*80, indate_wv(30)*8, indate_wd(30)*8,
     &          min_w*2, min_wav*2, wora*2, wore*2,date*10
 
      character regn1*4, seas1*1, seas2*1, fcstcurdir*13,list*11 
      character fn(3)*12,fcstfn(720)*16, wfcstfn(30)*14,wavfn*16,
     &          wavfcstfn(720)*18     
      character type_cu*3, type_wd*3, type_wv*3
      character outfile_traj*21,ahm*4,text*57

c----------------------------------------------------------------------

      integer   tt,ktmx_cu, ktmx_wd, ktmx_wv
      integer   nt,ix,count,is1

      logical   ex

c----------------------------------------------------------------------
c           FUNCTION: transformations between 
c                     the medslik grid coords and lat/lon
c----------------------------------------------------------------------

      xgrid(alon)=(alon-along1)/dlong+1.d0
      ygrid(alat)=(alat-alatg1)/dlatg+1.d0
      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      data hmin /1.d0/,    hmax /2500.d0/
      data a0 /'.cst','.srf','.dsp'/

      pi=4.*datan(1.d0)
      degrad=180.d0/pi
      open(90,file='medslik.log')
c      open(91,file='smag.log')


c-----------------------------------------------------------------------
c     LIST OF INPUT & OUTPUT files and UNIT numbers
c-----------------------------------------------------------------------
c      open(39,file='medslik5.par',status='old')
c      open(40,file='medslik5.inp',status='old')
c      open(41,file='medslik.tmp',status='old')
c      open(46,file='medslik.crn',status='old')
c      open(48,file='medslik.bms',status='old')
c
c      open(50,file='data/'//regn1//'.bath',status='old')
c      open(51,file='data/'//regn1//'.map',status='old')
c      open(52,file='data/'//regn1//'cst1.d',dispose='delete')
c      open(55,file='data/'//regn1//seas1//'.vel',status='old')
c      open(56,file='data/'//regn1//seas2//'.vel',status='old')

c      open(57,file='files-oce',status='old')
c      open(58,file='files-met',status='old')
c      open(59,file='files-wav',status='old')
c      open(71,file=fcstcurdir//fcstfile)
c      open(72,file=fcstwinddirr//fcstfile)
c
c      open(90,file='medslik.log')
c      open(99,file='output/medslik.fte')
c      open(81,file='output/'//outhhhh.cst)
c      open(82,file='output/'//outhhhh.srf)
c      open(83,file='output/'//outhhhh.dsp)
c
c      open(38,file='medslik_plgs.inp') satellite or user defined 
c                                       slick data file
c-------------------------------------------------------------------------


c============================================================================================
c      READ medslik5.par (model options) 
c============================================================================================
c ----- Physical Parameters -----
c
c       istoke=      Stokes drift: 0=not_computed, 1=computed_using_Jonswap, 
c                                  2=computed_using_wave_model_outputs
c       alpha=       Wind correction (Drift Factor)
c       beta0=       Wind correction (Drift Angle at zero wind speed)
c       ibetared=    0/1, variable Drift Angle: if drift angle reduces 
c                    with wind speed (beta reduces by 50% at wind speed = halfspeed)
c       halfspeed=   Wind speed at which drift angle is halved (used if previous entry is 1)
c                    beta=beta0*(1.-ibetared*wvel/(wvel+halfspeed))
c       iwindred=    0/1 if effective wind speed is reduced 
c                    (fraction of forecast wind is subtracted in drift formula)
c       wredfrac=    Reduction Fraction (used if previous entry is 1), fraction of fcst wind 
c                    to be subtracted when fcst currents are used
c       currincr=    Factor added to increment currents velocity 
c                    (FOR REASEARCH PURPOSES IT HAS TO BE 0.0)
c       ismag=       0/1 if Smagorinski scheme is to be used for horiz diffus 
c                    (only for forecast currents)
c       horizk=      Horizontal Diffusivity
c       vertk1=      Vertical Diffusivity above the thermocline
c       vertk2=      Vertical Diffusivity below the thermocline
c       thermocl=    Depth of the thermocline
c       ntot=        No of parcels used to model diffusion and dispersion
c       fcstdepth1|  
c       fcstdepth2|= Depths (3) at which velocities are printed in the forecast files   
c       fcstdepth3|  
c       idepth=      Depth level of water velocity used for slick convection - surface = 0
c----------------------------------------------------------------------
      open(39,file='medslik5.par',status='old')

      do i=1,4
         read(39,*) empty
      enddo 
      read(39,*) istoke   
      read(39,*) alpha    
      read(39,*) beta0    
      read(39,*) ibetared 
      read(39,*) halfspeed
      read(39,*) iwindred 
      read(39,*) wredfrac 
      read(39,*) currincr    
      read(39,*) ismag    
      read(39,*) horizk   
      read(39,*) vertk1   
      read(39,*) vertk2   
      read(39,*) thermocl 
      read(39,*) ntot     
      read(39,*) fcstdep1,fcstdep2,fcstdep3 
      read(39,*) idepth

      if(ntot.gt.100000) then
         write(6,*) 'Total number of parcels cannot exceed 100,000'
         write(6,*) 'Return to Parameters Form and change this number'
         stop
      endif
c---------------------------------------------------------------------------------
c ----- Evaporation Parameters -----
c Evaporation constants from Mackay et al.
c
c       ce=    12.0   Controls change in vapour pressure with evaporated fraction, 
c                     coeff accounts for drop in vapour pressure with evaporation 
c                     (ce=10-20)
c       akew=  3.3e-5 Controls overall evaporation rate
c       gamma= 0.78   Exponent of wind speed in evaporation rate
c       visk=  4.0    Controls effect of evaporated fraction on oil viscosity 
c
c       Use of akew and gamma:
c
c              ce1=akew*(wvel*3.6)**gamma (evaporative exposure to wind)
c
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) ce     
      read(39,*) akew   
      read(39,*) gamma  
      read(39,*) visk   
c----------------------------------------------------------------------
c ----- Emulsification Parameters -----
c       Emulsion constants from Mackay et al
c
c       cm1=   0.65    controls the effect of water fraction on mousse viscosity
c       cm2=   1.6e-06 controls the rate of water absorption
c       cm3=   1.333   controls maximum water fraction in the mousse 
c                      (decreased for heavy oils), reciprocal of maximum water fraction
c       icm3=  1       1/0 if maximum water fraction does/does not increase with API
c       visemx=        maximum mousse viscosity after which emulsification stops
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) cm1   
      read(39,*) cm2   
      read(39,*) cm3  
      read(39,*) icm3 
      visemx=100000.d0
c----------------------------------------------------------------------
c ----- Dispersion Parameters -----
c       Dispersion constants from Mackay et al
c
c       cd1=   0.001   downward diffusion velocity of small droplets (m/s)
c       cd3=   0.8e-05 controls the rate of dispersion of all droplets by waves
c       cd4=   50.0    controls the fraction of droplets below the critical size
c       cd5=   2000.0  controls the dispersion from the thin slick (sheen)
c       vs1=   0.08    rising velocity of small droplets (m/s)
c       vl=    0.0003  rising velocity of large droplets (m/s)
c       um=    0.5     controls depth of well-mixed surface layer (m)
c       st=    24.0    interfacial tension between oil and water
c       fmaxd= 1.0     max dispersive fraction
c----------------------------------------------------------------------
      read(39,*) empty
      read(39,*) cd1   
      read(39,*) cd3   
      read(39,*) cd4   
      read(39,*) cd5   
      read(39,*) vl    
      read(39,*) vs1   
      read(39,*) um    
      read(39,*) st    
      read(39,*) fmaxd 
      stk=st
      stn=st
c--------------------------------------------------------------------------------
c ----- Spreading Parameters -----
c       Spreading constants from Mackay et al
c
c       cs1= 1.0    Controls rate of spreading of thin slick
c       cs2= 150.0  Controls rate of spreading of thick slick
c       cs3= 0.0015 Controls dependence of spreading of thin slick on thickness
c       sprdmx =    max time for spreading of spill (hours)
c--------------------------------------------------------------------------------

      read(39,*) empty
      read(39,*) cs1    
      read(39,*) cs2    
      read(39,*) cs3    
      sprdmx=24.d0

c-----------------------------------------------------------------------------
c       seepmx=   5000.0 limiting concentration on coast (bbls/km)
c       apicoeff= 0.0    coeff of reduction of coastalretention rate 
c                        for heavy oils
c-----------------------------------------------------------------------------

      data seepmx /5000.d0/

       read(39,*) empty
       read(39,*) seepmx     ! 5000.0
       read(39,*) apicoeff   ! 0.0

       close(39)

c======================================================================
c    TIME STEP RELATED VARIABLES COMPUTATION
c 
c       delt = time increment for spill re-computations (hrs)
c	stph, nstph = no of steps per hour
c======================================================================

       deltsc=delt*3600.d0
       stph=1.d0/delt
       nstph=stph+0.001d0

c============================================================================================
c      READ medslik5.inp (user inputs) 
c============================================================================================
c
c       regn1=      4-character name of the region where the spill is located
c       irestart= 1 if a previous run is to be restarted
c       ihrestart=  number of hours of this previous run
c       numscheme=  0/1 for Eulero Forward/Fourth Order Runge-Kutta
c                   numerical schemes
c       trackmode=  Lagrangian_path_simulated=1, 
c                   Oil-Slick/Spill_fate_and_transport_simulated=0
c       numspills=  number of simulated slicks/spills
c----------------------------------------------------------------------

       open(40,file='medslik5.inp',status='old')

       read(40,'(a4)') regn1
       read(40,*) irestart,ihrestart
       read(40,*) numscheme
       read(40,*) trackmode
       read(40,*) numspills

       if (regn1.eq.'medf') then
              iregn=0
       elseif (regn1.eq.'emed'.or.regn1.eq.'cyba'.or.regn1.eq.'cyse'.or.
     &        regn1.eq.'syri') then 
              iregn=1
       elseif (regn1.eq.'adri'.or.regn1.eq.'adno'.or.regn1.eq.'anza')
     &        then
              iregn=2
       elseif(regn1.eq.'sici'.or.regn1.eq.'mlts'.or.regn1.eq.'mltc') 
     &        then
              iregn=3
       elseif(regn1.eq.'tyrr') then
              iregn=4
       else
              iregn=10
       endif

       if(irestart.eq.0) ihrestart=0

       call setcodes(regn1)

c----------------------------------------------------------------------
c     --- Read Variables ---
c
c       idd/imm/iyr= date of the first spill/slick specified in 
c                    the medslik_inputfile.txt
c       istart=      Start hour and minutes of the first spill/slick 
c                    specified in the medslik_inputfile.txt 
c                    (format = hhmm)
c       pref=        3-letter prefix for labelling output files
c       icomp=       duration of computation from start of spill (hrs)
c       iprs=        interval for output (hrs)
c       icurrents=   index that specifies source of ocean currents data
c       iwind=       index that specifies source of wind data
c       iwave=       index that specifies source of ocean waves data
c       icrn =       1/0 if observation of spill posn is (is not) to be applied
c
c     --- Derived Variables ---
c
c       tstart = nearest hour to start of spill (time in decimal hours)
c       tspill = ispill
c       tcomp = icomp
c       splat = Latitude in decimal grads
c       splon = Longitude in decimal grads
c----------------------------------------------------------------------

      read(40,'(i2,1x,i2,1x,i4)') idd,imm,iyr
      read(40,'(i4)') istart
      read(40,'(a4)') pref
      read(40,'(i4)') icomp
      read(40,'(i3)') iprs
      read(40,'(i2)') icurrents
      read(40,'(i2)') iwind
      read(40,'(i3)') iwave 
      read(40,'(i2)') icrn
      
      tstart=dfloat(istart/100)+dfloat(istart-100*(istart/100))/60.d0
      tspill=dfloat(ispill)
      tcomp=dfloat(icomp)
c----------------------------------------------------------------------
c       --- Read Variables ---
c       vlunit=  unit of volume (tons, cu.m, bbls, gals)
c       iage=    age of the slick (in hours, 0, 24, 48)
c       isat=    0: point source;
c                1: areal source of spill (for manual slick contour or from 
c                                          satellite data)
c       api=     API number of the oil
c       deno=    oil density  (g/cm**3)
c       den2=    density of residual part (g/cm**3)
c       respc=   residual percentage of oil
c       viso=    initial oil viscosity
c       tem0=    Temperature at which Viscosity determined
c       vappr=   oil vapour pressure of oil (bar)
c
c     max water content reduced/increased for heavy/light oils
c
c     --- Derived Variables ---
c
c     tvk0 = reference temperature for viscosity in Kelvin grads
c     po = vappr
c----------------------------------------------------------------------


      read(40,'(a4)') vlunit
      read(40,*) iage       
      read(40,*) isat
      print *, 'isat= ', isat, 'iage', iage
      read(40,'(a80)') empty
      read(40,'(f5.2)') api
      read(40,'(f5.3)') deno
      read(40,'(f5.3)') den2
      read(40,'(f5.2)') respc
      read(40,'(f5.1)') viso
      read(40,'(f5.1)') tem0
      read(40,'(f5.1)') vappr

      tvk0=tem0+273.d0
      if(icm3.eq.1.and.cm3.eq.1.333d0) then
        cm3=(10.d0/9.d0)*(2.d0-1.d0/(1.d0+4.d0**(1.d0-api/17.d0)))
      end if
      
      po=vappr

c----------------------------------------------------------------------
c     --- Read Variables --- 
c
c       isst=   index that specifies source of SST data
c       tc=     SST in deg C - read only if isst = 8, otherwise zero
c       ibooms= 1/0 if booms are (are not) deployed
c       al5=    output pixel size (m) for counting up slick parcels
c
c     --- Derived Variables ---
c
c     gridkm = output pixel size (km)
c     area = area (km^2) of the output horizontal grid cell
c
c----------------------------------------------------------------------

      read(40,'(i2)') isst
      read(40,'(f4.1)') tc
      read(40,'(i2)') ibooms
      read(40,'(f5.1)') al5

      gridkm=al5/1000.d0
      area=gridkm*gridkm

c---------------------------------------------------------------------------------
c
c       fcstfn()=   name of forecast file of currents
c       ifcstfn()=  1/0 if the file is (is not) available
c                   array of n elements where  n = 30(dd) * num_timestep_cur_data
c       wfcstfn()=  name of forecast file of winds
c       iwfcstfn()= 1/0 if the wind file is (is not) available
c       wavfcstfn=  name of forecast file of waves (Added by D. Bruciaferri)
c       iwavfcstfn= 1/0 if the file is (is not) available
c                   array of n elements where  n = 30(dd) * num_timestep_wav_data
c
c---------------------------------------------------------------------------------

c  *** OCEAN DATA FILES ***

      if ( (icurrents.eq.76).or.(icurrents.eq.77)
     &      .or.(icurrents.eq.14) ) then

         iprod=1
         read(40,'(i4)') nfcst

         open(57,file='files-oce')
         read(57,'(a3)') type_cu
         read(57,'(i2)') ktmx_cu
         close(57)
       
         open(41,file='medslik.tmp')
         do n=1,3
            read(41,*) dummy
         enddo
         read(41,*) numfiles

         if (icurrents.eq.14) then
            do n=1,numfiles
               read(41,'(a8)') indate(n)
            enddo
         else
            do n=1,numfiles+1
               read(41,'(a8)') indate(n)
            enddo
         endif

         do i=1,nfcst

c        Daily forecast files
 
            if (icurrents.eq.14) then

              do k=1,3
                    read(40,'(a11,i2)') fn(k),ifcstfn(i)
              enddo
              fcstfn(i) = 'meds'//fn(1)(1:6)//'00.med'

            endif

c        Hourly forecast files
c---------------------------------------------------------------------------
c 1. type_cu = hm -> hourly mean fields (for example from 8:00 to 9:00, 
c                    centered at 8:30) So, we write files referred to 30 min 
c                    after and then medslik-II subtract 30 when reads the  
c                    file in input.
c 2. type_cu = hi -> Instantaneous fields
c----------------------------------------------------------------------------

            if (icurrents.ge.76) then 

               do k=1,3
                  read(40,'(a11,i2)') fn(k),ifcstfn(i)
               enddo

               list=fn(1)

               do nt=1,ktmx_cu

c DATE LIST FOR HOURLY CURRENTS FILES STARTING FROM 00:00

                  if (ktmx_cu.eq.24) then
! ------ 1hr resolution (forecast) -------
                     if (type_cu.eq.' hm') then
! Hourly mean fields
                        if (nt.lt.10) then
                           write(nore,'(i2)') nt
                           nora='0'//nore(2:2)
                        else
                           write(nora,'(i2)') nt
                        endif

                     else
! Instantaneous fields                                         
                        if (nt.le.10) then
                           write(nore,'(i2)') nt-1
                           nora='0'//nore(2:2)
                        else
                           write(nora,'(i2)') nt-1
                        endif

                     endif
                  else
! ------ Different time resolution (analysis/historical data) -------
                     tt = ((24 / ktmx_cu)*(nt-1))
                     if (type_cu.eq.' hm') then
! Hourly mean fields
                        if (tt.lt.9) then
                           write(nore,'(i2)') tt+1
                           nora='0'//nore(2:2)
                        else
                           write(nora,'(i2)') tt+1
                        endif
! Instantaneous fields
                     else
                        if (tt.lt.10) then
                           write(nore,'(i2)') tt
                           nora='0'//nore(2:2)
                        else
                           write(nora,'(i2)') tt
                        endif
                     endif
                  endif

c DATE LIST FOR HOURLY CURRENTS FILES STARTING FROM 12:00 

                  if ((icurrents.eq.76).or.(icurrents.eq.77)) then
                     if (nt.le.12) then
                        count=i
                        write(ore,'(i2)') nt+12
                        ora=ore(1:2)
                     endif
                     if (nt.gt.12.and.nt.lt.22) then
                        count=i+1
                        write(ore,'(i2)') nt-12
                        ora='0'//ore(2:2)
                     endif
                     if (nt.ge.22.and.nt.le.24) then
                        count=i+1
                        write(ore,'(i2)') nt-12
                        ora=ore(1:2)
                     endif
                  endif

c-----------------------------------------------------------------------
c OCEAN DATA FILENAME 

                  if ( (icurrents.eq.76).or.(icurrents.eq.77) ) then
                     fcstfn((i-1)*ktmx_cu+nt)=
     &               'meds'//indate(count)(1:6)//ora(1:2)//'.med'
                  else
                     fcstfn((i-1)*ktmx_cu+nt)=
     &               'meds'//list(1:6)//nora(1:2)//'.med'
                  endif

                  ifcstfn((i-1)*ktmx_cu+nt)=ifcstfn(i)

               enddo
            endif
            iprod=iprod*ifcstfn(i)
         enddo
c************************************************************************

        if (icurrents.ge.1.and.icurrents.le.80) then
           if ( (icurrents.ge.1.and.icurrents.le.80).and.ifcstfn(1)
     &           .eq.0.or.ifcstfn(2).eq.0 ) then

          write(6,*) '*************************************************'
          write(6,*) 'The first two forecast files MUST be present'
          write(6,*) 'The required files are: '
          write(6,*) '               ',fcstfn(1),', ',fcstfn(2)
          write(6,*) '*************************************************'
          stop

           end if
        endif

        if(iprod.eq.0) then

          write(6,*) '*************************************************'
          write(6,*) 'WARNING: The following forecast files are missing'

          do i=1,nfcst
             if (ifcstfn(i).eq.0.and.icurrents.ge.70) then 
                write(6,*) fn(1),', ',fn(2),', ',fn(3)
             end if

             if (ifcstfn(i).eq.0.and.icurrents.ge.70) then
                write(6,*) fcstfn(i)
             end if
          end do

          write(6,*) 'The simulation will extrapolate from the last'
          write(6,*) 'available file and so will lose reliability.'
          write(6,*) '*************************************************'

        end if

      else
      
        if(ismag.eq.1) then
          write(6,*) '*************************************************'
          write(6,*) 'Smagorinsky Diffusivity Model can ONLY be used'
          write(6,*) '    when FORECAST Water Currents are used.'
          write(6,*) 'Return to the Startup Screen and enter a value' 
          write(6,*) '         for Horizontal Diffusivity' 
          write(6,*) '*************************************************'
          stop
        end if   
      
      end if

c *** METEO DATA FILES ***

      if ( (iwind.eq.11).or.(iwind.eq.12) ) then

         open(58,file='files-met')
         read(58,'(a2)') type_wd
         read(58,'(i2)') ktmx_wd
         close(58)
  
         read(41,*) numfiles_wd

         do n=1,numfiles_wd
            read(41,'(a8)') indate_wd(n)
         enddo

         iprod=1
         read(40,'(i3)') nwfcst
         do i=1,nwfcst
            read(40,*) wfcstfn(i),iwfcstfn(i)
            iprod=iprod*iwfcstfn(i)
         end do

         if (iwfcstfn(1).eq.0.or.(nwfcst.gt.1.and.
     &      iwfcstfn(2).eq.0)) then

          write(6,*) '*************************************************'
          write(6,*) 'The first two wind forecast files MUST be present'
          write(6,*) 'The required files are: '
          write(6,*) '               ',wfcstfn(1),', ',wfcstfn(2)
          write(6,*) '*************************************************'
          stop
         end if

         if (iprod.eq.0) then
          write(6,*) '*************************************************'
          write(6,*) 'WARNING: The following forecast files are missing'
          do i=1,nwfcst
             if(iwfcstfn(i).eq.0) write(6,*) wfcstfn(i)
          end do
          write(6,*) 'The simulation will extrapolate from the last'
          write(6,*) 'available file and so will lose reliability.'
          write(6,*) '*************************************************'
         end if
 
      end if 

c *** WAVE DATA FILES ***

      if (iwave.gt.0) then

         open(59,file='files-wav')
         read(59,'(a3)') type_wv
         read(59,'(i2)') ktmx_wv
         close(59)       

         read(41,*) numfiles_wv

         if (iwave.eq.101) then
            numfiles_wv=numfiles_wv+1
         endif

         do n=1,numfiles_wv
            read(41,'(a8)') indate_wv(n)
         enddo

         close(41)


         iprod=1
         read(40,'(i3)') nwavfcst

         do i=1,nwavfcst

c        Hourly forecast files
c---------------------------------------------------------------------------
c 1. type_cu = hm -> average hourly fields (for example from 8:00 to 9:00, 
c                    centered at 8:30) So, we write files referred to 30 min 
c                    after and then medslik-II subtract 30 when reads the  
c                    file in input.
c 2. type_cu = hi -> Instantaneous fields
c----------------------------------------------------------------------------

            read(40,'(a16,i2)') wavfn,iwavfcstfn(i)
            write(6,'(a16,i2)') wavfn,iwavfcstfn(i)

            do nt=1,ktmx_wv

c DATE LIST FOR HOURLY WAVES FILES STARTING FROM 00:00

               if (ktmx_wv.eq.24) then
! Instantaneous fields
                  if (type_wv.eq.' hi') then
                      if (nt.le.10) then
                           write(nore,'(i2)') nt-1
                           nora='0'//nore(2:2)
                       else
                           write(nora,'(i2)') nt-1
                       endif
                  else 
                       if (nt.lt.10) then
                           write(nore,'(i2)') nt
                           nora='0'//nore(2:2)
                       else
                           write(nora,'(i2)') nt
                       endif
                  endif
               else
! ------ Different time  resolutions (forecast/analysis/historical data) -------
                  tt = ((24 / ktmx_wv)*(nt-1))
! Instantaneous fields
                  if (type_wv.eq.' hi') then
                     if (tt.lt.10) then
                         write(nore,'(i2)') tt
                         nora='0'//nore(2:2)
                      else
                         write(nora,'(i2)') tt
                      endif
                  else
! Hourly mean fields
                      if (tt.lt.9) then
                          write(nore,'(i2)') tt+1
                          nora='0'//nore(2:2)
                      else
                           write(nora,'(i2)') tt+1
                      endif
                  endif
               endif

c DATE LIST FOR HOURLY WAVE FILES STARTING FROM 12:00 

               if (nt.le.12) then
                  count=i
                  write(wore,'(i2)') nt+12
                  wora=wore(1:2)
               endif
               if (nt.gt.12.and.nt.lt.22) then
                  count=i+1
                  write(wore,'(i2)') nt-12
                  wora='0'//wore(2:2)
               endif
               if (nt.ge.22.and.nt.le.24) then
                  count=i+1
                  write(wore,'(i2)') nt-12
                  wora=wore(1:2)
               endif

c-----------------------------------------------------------------------
c OCEAN WAVES DATA FILENAME 
 
               if (iwave.eq.101) then
                  wavfcstfn((i-1)*ktmx_wv+nt)='wav_'//
     &            indate_wv(count)(3:8)//wora(1:2)//'00.wav'
               else
                  wavfcstfn((i-1)*ktmx_wv+nt)='wav_'//wavfn(3:8)//
     &            nora(1:2)//'.wav'
               endif
               
               iwavfcstfn((i-1)*ktmx_wv+nt)=iwavfcstfn(i)

            enddo
         enddo

         if (iwavfcstfn(1).eq.0.or.(nwavfcst.gt.1.and.
     &       iwavfcstfn(2).eq.0)) then

             write(6,*) '***************************
     &                   **********************'
             write(6,*) 'The first two wave forecast files MUST 
     &                   be present'
             write(6,*) 'The required files are: '
             write(6,*) '               ',wavfcstfn(1),', ',wavfcstfn(2)
             write(6,*) '****************************
     &                   *********************'
             stop

          end if

          if(iprod.eq.0) then

             write(6,*) '***************************
     &                   **********************'
             write(6,*) 'WARNING: The following forecast 
     &                   files are missing'
             do i=1,nwavfcst
                if(iwavfcstfn(i).eq.0) write(6,*) wavfcstfn(i)
             end do
             write(6,*) 'The simulation will extrapolate from the last'
             write(6,*) 'available file and so will lose reliability.'
             write(6,*) '***************************
     &                   ********************'
          end if
        
        end if      
        close(40)

c=========================================================================
c                      READ OIL SLICKS/SPILLS DATA
c=========================================================================
c     --- Read Variables ---
c
c FOR SINGLE POINT SOURCE
c
c     ispill=       duration of the spill
c     splat splon=  geographical location of the spill (decimal degrees)
c     splrte=       spill rate of the spill (no units per hour) or 
c                   total volume of the spill for ispill(i)=0
c
c FOR MULTPILE POINT SOURCES OR 
c USER SPECIFIED SINGLE/MULTIPLE POLYGONAL SLICKS    
c
c     iddv(i)/immv(i)/iyrv(i)= date of the ith spill
c     istartv(i)=              time of day at start of ith spill 
c                              (hhmm of the ith spill)
c     ispillv(i)=              duration of the ith spill (hours)
c     splatv(i) splonv(i)=     geographical location of ith 
c                              the spill (decimal degree)
c     splrtv(i)=               spill rate of the ith spill 
c                              (no units per hour) or total volume 
c                              of the ith spill for ispillv(i)=0
c
c     --- Derived Variables ---
c
c FOR SINGLE POINT SOURCE
c
c     tspill = ispill
c
c     splatmx = Max Latitude of the spills
c     splatmn = Min Latitude of the spills
c     splonmx = Max Longitude of the spills
c     splonmn = Min Longitude of the spills 
c
c FOR MULTPILE POINT SOURCES OR 
c USER SPECIFIED SINGLE/MULTIPLE POLYGONAL SLICKS
c     
c     tstartv(i)=  nearest hour to start of the ith spill 
c                  (time in decimal hours)
c     tspillv(i)=  ispillv(i)
c-------------------------------------------------------------------------
      
      if (numspills.eq.1) then

         if (isat.eq.0) then
            
            open(41,file='medslik_pnts.inp',status='old')
            do i=1,4
               read(41,*) empty
            enddo
            read(41,*) ispill 
            read(41,*) splat
            read(41,*) splon
            read(41,*) splrte
            tspill=dfloat(ispill)

            splatmx = splat
            splatmn = splat
            splonmx = splon
            splonmn = splon

            print*,'=================================='
            print*, "POINT 1"
            print*,''
            print*, "Date:", idd,imm,iyr
            print*, "Time: ", istart 
            print*, "Spill duration: ", ispill 
            print*, "Spill rate: ", splrte
            print*, "Spill latitude: ", splat
            print*, "Spill longitude: ", splon
     
         endif

         if (isat.eq.1) then

            iddv(1) = idd
            immv(1) = imm
            iyrv(1) = iyr
            istartv(1) = istart
            tstartv(1) = tstart

            open(41,file='medslik_plgs.inp',status='old')

            do i=1,2
               read(41,*) empty
            enddo

            read(41,*) ispillv(1)

            if (ispillv(1).ne.0) then
               print*, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
               print*, '           !!! ERROR !!!'
               print*, ''
               print*, 'For satellite detected or hand builded'
               print*, 'oil slicks the duration of the release'
               print*, '            MUST BE ZERO'
               print*, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
               stop
            endif
            read(41,*) empty
            read(41,*) splrtev(1)
            tspillv(1)=dfloat(ispillv(1)) 
            tstartv(1)=tstart

            print*,'=================================='
            print*, "POLYGON 1"
            print*,''
            print*, "Date:", iddv(1),immv(1),iyrv(1)
            print*, "Time: ", istartv(1)
            print*, "Spill duration: ", ispillv(1)
            print*, "Spill rate: ", splrtev(1)

            read(41,*) empty
            read(41,*) empty
            read(41,*) splat, splon
            splatmx = splat
            splatmn = splat
            splonmx = splon
            splonmn = splon 
                      
         endif

      endif

      if (numspills.gt.1) then

         if (isat.eq.0) then

            open(41,file='medslik_pnts.inp',status='old')
            read(41,*) empty

            do i=1,numspills
               read(41,*) empty
               read(41,'(i2,2x,i2,1x,i4)') iddv(i), immv(i), iyrv(i)
               read(41,*) istartv(i)
               read(41,*) ispillv(i)
               read(41,*) splatv(i)
               read(41,*) splonv(i)
               read(41,*) splrtev(i)

               iyrv(i) = iyrv(i) + 2000
               tspill = dfloat(ispill)
               tstartv(i) = dfloat(istartv(i)/100) +
     &                      dfloat(istartv(i) - 100*
     &                      (istartv(i)/100)) / 60.d0
               if (i.eq.1) then
                  splatmx = splatv(1)
                  splatmn = splatv(1)
                  splonmx = splonv(1)
                  splonmn = splonv(1)
               else
                  if(splatv(i).gt.splatmx) splatmx = splatv(i)
                  if(splatv(i).lt.splatmn) splatmn = splatv(i)
                  if(splonv(i).gt.splonmx) splonmx = splonv(i)
                  if(splonv(i).lt.splonmn) splonmn = splonv(i)
               endif

               print*,'=================================='
               print*, "POINT ", i
               print*,''
               print*, "Date:", iddv(i),immv(i),iyrv(i)
               print*, "Time: ", istartv(i)
               print*, "Spill duration: ", ispillv(i)
               print*, "Spill rate: ", splrtev(i)
               print*, "Spill latitude: ", splatv(i)
               print*, "Spill longitude: ", splonv(i)
            enddo

         endif

         if (isat.eq.1) then

            open(41,file='medslik_plgs.inp',status='old')

            do i=1,2
               read(41,*) empty
            enddo

            read(41,*) ispill

            if (ispill.ne.0) then
               print*, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
               print*, '           !!! ERROR !!!'
               print*, ''
               print*, 'For satellite detected or hand builded'
               print*, 'oil slicks the duration of the release'
               print*, '            MUST BE ZERO'
               print*, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
               stop
            endif

            read(41,*) empty

            do i=1,numspills
               read(41,*) splrtev(i)
               iddv(i) = idd
               immv(i) = imm
               iyrv(i) = iyr
               istartv(i) = istart
               tstartv(i) = tstart
               ispillv(i) = ispill
               tspillv(i) = dfloat(ispillv(i))

               print*,'=================================='
               print*, "POLYGON ", i
               print*,''
               print*, "Date:", iddv(i),immv(i),iyrv(i)
               print*, "Time: ", istartv(i)
               print*, "Spill duration: ", ispillv(i)
               print*, "Spill rate: ", splrtev(i)
            enddo
            read(41,*) empty
            read(41,*) empty
            read(41,*) splat, splon
            splatmx = splat
            splatmn = splat
            splonmx = splon
            splonmn = splon 
         endif
      endif

      close(41)

c====================================================================
c     INITIAL PARTICLE POSITIONS FOR POLYGONAL SLICKS
c====================================================================

      if (isat.eq.1) then

         call seedmedslik(ix)
         call readsat(ix,ntot,numspills,nppmsv,px,py)

         write(6,*) '========================================='
         write(6,*) 'POLYGON SLICK: READING SLICK CONTOUR DATA'

         ntot = 0

         do i=1,numspills
            ntot = ntot + nppmsv(i)
         enddo

         do i=1,ntot

            if(py(i).gt.splatmx) splatmx = py(i)
            if(py(i).lt.splatmn) splatmn = py(i)
            if(px(i).gt.splonmx) splonmx = px(i)
            if(px(i).lt.splonmn) splonmn = px(i)

         enddo

      endif
c-------------------------------------------------------------------------
c START - TRAJECTORY output file - HEADING
c-------------------------------------------------------------------------
      if (trackmode.eq.1) then

          write(a2,'(a2)') '_F'
          if(iyr.ge.2000) write(ay,'(i2.2)') iyr-2000
          if(iyr.lt.2000) write(ay,'(i2.2)') iyr-1900
          write(am,'(i2.2)') imm
          write(ad,'(i2.2)') idd
          write(ahm,'(i4.4)') istart
          write(a3,'(i3.3)') ihrestart
          outfile_traj=pref//ay//am//ad//'_'//ahm//a2//'.trj'
          open(93,file = outfile_traj)
          ndefrec = 0
          write(93,*) "MEDSLIK-II trajectory simulation output:"
          write(93,'(a25)')'Time  GC_lat  GC_lon'
      endif
c-------------------------------------------------------------------------
c END - TRAJECTORY output file - HEADING
c-------------------------------------------------------------------------       
      
c----------------------------------------------------------------------
c     Assign directory for forecast current data
c     first check if the region is an added region
c----------------------------------------------------------------------
      if (icurrents.ge.1.and.icurrents.le.80) then

          inquire(file='data/newregions.txt',EXIST=ex)
          if (ex) then
              open(1,file='data/newregions.txt')
              read(1,*) nnew
              do k=1,nnew
                 read(1,'(a4)') a4
                 read(1,'(a80)') empty
                 read(1,'(a80)') empty
                 if (a4.eq.regn1) then
                    read(1,*) d24,d06
                    iregn = 70 + k
                    go to 5 
                 else
                    read(1,'(a80)') empty
                 endif 
              enddo   
   5          continue
              close(1)
          endif

          a4 = fcstfn(1)(1:4)
          if ( (iregn.eq.1.or.iregn.eq.0).and.a4.eq.'meds' ) 
     &         fcstcurdir='INP_DATA/OCE/'
      endif

c======================================================================
c     CONSTRUCT BASIC GRID PARAMETERS AND BATHYMETRY
C======================================================================
c     delx, dely= bathymetry grid spacing in metres
c                 (set grid spacing delx, at mid-latitude of region)
c                 In case iregn = 0 (whole Mediterranean) restrict 
c                 region: assuming slick travels at less than 1.5 kts 
c                 from spill site.
c     rax, ray=   ratios of grid spacing to output pixel size
c======================================================================
            
      if (icencr.ne.0) then
         open(50,file='data/'//regn1//'_.bath',status='old')
      else
         open(50,file='data/'//regn1//'.bath',status='old')
      endif

      read(50,'(a80)') empty 
      read(50,*) along1,along2,alatg1,alatg2
      read(50,*) mmax,nmax
            
      dlong = (along2 - along1) / dfloat(mmax-1)
      dlatg = (alatg2 - alatg1) / dfloat(nmax-1)
      avlat = (alatg1 + alatg2) / 2.d0

      if (icencr.ne.0.or.iregn.eq.4.or.iregn.eq.2) then
 
         if (iregn.ne.0.and.iregn.ne.4.and.iregn.ne.2) then
          
            do n=nmax,1,-1
               do m=1,mmax 
                   read(50,'(i4)') ih1     
                  if (ih1.ne.9999) ih1 = xor(ih1, icencr)
                  itype(m,n) = ih1
               enddo
            enddo   

         else if (iregn.eq.0.or.iregn.eq.4.or.iregn.eq.2) then

            dist = 1.5d0 * tcomp                ! distance in NM
            dist = dist / 60.d0                 ! convert to deg of latitude

            avlat = (splatmx + splatmn) / 2.d0

            dep = dist / dcos(avlat/degrad)    ! same dist in deg longitude

            dsplon = (splonmx - splonmn)

            print*, '--------------------------------------------------'
            print*, '      GEOMETRIC PROPERTIES OF THE DOMAIN          '
            print*, '                                                  '
            print*, '           Max Lat = ', splatmx
            print*, '           Min Lat = ', splatmn
            print*, '           Max Lon = ', splonmx
            print*, '           Min Lon = ', splonmn
            print*, '           Mean Lat = ', avlat
            print*, '           Delta Lon = ', dsplon
            print*, '--------------------------------------------------'


            if(dep.gt.(dlong * (mm-4) - dsplon) / 2.d0) then

               dep = dlong * (mm-4) / 2.d0 - dsplon
               dist = dep * dcos(avlat/degrad)

            endif

c     ----------------------------------        

            alon1 = splonmn - dep
            alon2 = splonmx + dep

            if (alon1.lt.along1) alon1 = along1
            if (alon2.gt.along2) alon2 = along2

            m1 = int( xgrid(alon1) )
            m2 = int( xgrid(alon2) + 1.d0 )

            mmax1 = m2 - m1 + 1
            alon1 = glon(dfloat(m1))

c     ----------------------------------        

            alat1 = splatmn - dist
            alat2 = splatmx + dist
  
            if(alat1.lt.alatg1) alat1 = alatg1
            if(alat2.gt.alatg2) alat2 = alatg2

            n1 = int( ygrid(alat1) )
            n2 = int( ygrid(alat2) + 1.d0 )

            nmax1 = n2 - n1 + 1
            alat1 = glat(dfloat(n1)) 

c     ----------------------------------

            do n=nmax,1,-1
               do m=1,mmax

                   read(50,'(i4)') ih1     

                   if (n.ge.n1.and.n.le.n2.and.m.ge.
     &                m1.and.m.le.m2) then

                      mp = m - m1 + 1 
                      np = n - n1 + 1

                      if(ih1.ne.9999.and.iregn.ne.4.and.
     &                   iregn.ne.2) then
  
                         ih1 = xor(ih1, icencr)
                         itype(mp,np) = ih1 

                      endif

                      if(ih1.ne.9999.and.iregn.eq.4.and.
     &                   iregn.eq.2) ih1 = ih1

                      itype(mp,np) = ih1 

                   endif  

               enddo
            enddo
 
            mmax = mmax1
            nmax = nmax1
            along1 = alon1      
            along2 = alon1 + (mmax-1) * dlong      
            alatg1 = alat1      
            alatg2 = alat1 + (nmax-1) * dlatg      

        endif
        
      else            

        do n=nmax,1,-1
           read(50,*) (itype(m,n),m=1,mmax)
        enddo

      endif

      close(50)

      avlat = (alatg1 + alatg2) / 2.d0
      dely = dlatg*60.d0*1849.d0
      delx = dlong*60.d0*1849.d0*dcos(avlat/degrad)

      rax = delx / al5
      ray = dely / al5

      hmax = 0.
      do n=nmax,1,-1
        do m=1,mmax
          if(itype(m,n).eq.9999) itype(m,n) = 0
          h(m,n) = float( itype(m,n) )
          if(itype(m,n).gt.0) itype(m,n) = 1
          if(h(m,n).lt.hmin.and.itype(m,n).ne.0) h(m,n) = hmin
          if(h(m,n).gt.hmax.and.itype(m,n).ne.0) hmax = h(m,n)
        enddo
      enddo

c      call hsmoot(h)

c==========================================================================
c     vertd1,2= vertical  diffusion distance during time delt
c     horizd=   horizontal diffusion distance during time delt
c
c     check stability of horizontal diffusion simulation (ismag=0)
c
c     ismag = 1 if Smagorinski scheme used for horizd (forecast data only)
c==========================================================================

      vertd1=dsqrt(6.d0*vertk1*delt*3600.d0)
      vertd2=dsqrt(6.d0*vertk2*delt*3600.d0)

      if (ismag.eq.0) then
         horizd=dsqrt(6.d0*horizk*delt*3600.d0)
         dtmx1=delx**2/(6.d0*horizk*3600.d0)
         dtmx2=dely**2/(6.d0*horizk*3600.d0)
         if (delt.gt.dtmx1.or.delt.gt.dtmx2) then
            write(6,*)'delt must not exceed ',dtmx1,' or ',dtmx2,' hrs'
            stop
         end if
      end if

c======================================================================
c     x0,y0 = location of spill in coords of the bathymetry grid
c
c     Rough check that spill is inside the given water body
c
c     Write headings to screen and log file
c======================================================================
      
      if (isat.eq.0.and.numspills.gt.1) then

         do k=1,numspills

            x0v(k)=xgrid(splonv(k))
            y0v(k)=ygrid(splatv(k))

         enddo

         m0=int(x0v(1))
         n0=int(y0v(1))

      else

         x0=xgrid(splon)
         y0=ygrid(splat)

         m0=int(x0)
         n0=int(y0)

      endif

      isum=itype(m0,n0)+itype(m0+1,n0)+itype(m0,n0+1)
     &     +itype(m0+1,n0+1)

      if(m0.lt.1.or.m0.gt.mmax-1.or.n0.lt.
     &   1.or.n0.gt.nmax-1.or.isum.eq.0) then

        write(6,*) '************************************************'
        write(6,*) ' Spill location is outside the given water body'
        write(6,*) ' You should continue only if using restart file'
        write(6,*) '************************************************'

      else if(isum.le.2) then

        write(6,*) '************************************************'
        write(6,*) 'WARNING:  Spill location is very close to a coast'
        write(6,*) '   and  may be outside the given water body.'
        write(6,*) 'Check its latitude and longitude and if necessary'
        write(6,*) '      move it further into the water body'
        write(6,*) '************************************************'

      endif

      write(6,'(/'' Welcome to the MEDSLIK run module''/)')
      write(6,*) 'The spill simulated has the following parameters:'
      write(6,'('' The spill is located in the region '',a4)') regn1
      write(6,600) idd,imm,iyr,istart
      write(6,*) 'Geog, and grid coordinates of the spill points:'

      if (isat.eq.0.and.numspills.gt.1) then

          do k=1,numspills
             write(6,605) k,splatv(k),splonv(k),x0v(k),y0v(k),
     &                    splrtev(k),vlunit
          enddo

      else 

          write(6,601) splat,splon,x0,y0,tspill,icomp 

      endif

      if (ispill.eq.0) write(6,603) splrte,vlunit
      if (ismag.eq.0) write(6,604) horizd
      
      if ((icurrents.eq.76).or.(icurrents.eq.77).or.(icurrents.eq.14)) 
     &   write(6,645)

      if ( (iwind.eq.11).or.(iwind.eq.12) ) write(6,703)

  600 format(' Date:           ',i2,'/',i2,'/',i4,'  hour: ',i4.4)    
  601 format(' Location:  lat  ',f8.4,',    lon ',f8.4/
     &       ' Grid coords: x  ',f8.4,',      y ',f8.4/
     &       ' Spill duration: ',f7.0,' hours'/
     &       ' Length of run: ',i4,' hours')
  602 format(' Spill rate:    ',f9.2,' ',a4,'/hr')    
  603 format(' Total volume:  ',f9.2,' ',a4)    
  604 format(' Horiz diffusion distance:  ',f12.5)    
  605 format(' Spill no:       ',i3/
     &       ' Location:  lat  ',f8.4,',    lon ',f8.4/
     &       ' Grid coords: x  ',f8.4,',      y ',f8.4/
     &       ' Total volume:   ',f9.2,' ',a4)
  645 format(/' Water currents will be read directly from the'
     &        ' Mediterranean Forecasting System netcdf output files'/)
  703 format(' Winds will be read from ECMWF 6-hourly '
     &       ' netcdf output files'/)

  770 format(/'Spill positions will be corrected from observations')

      write(90,'(/'' Welcome to the MEDSLIK run module''/)')
      write(90,*) 'The spill simulated has the following parameters:'
      write(90,'('' The spill is located in the region '',a4)') regn1
      write(90,600) idd,imm,iyr,istart
      write(90,601) splat,splon,x0,y0,tspill,icomp 
      if(ispill.gt.0) write(90,602) splrte,vlunit
      if(ispill.eq.0) write(90,603) splrte,vlunit
      if(ismag.eq.0) write(90,604) horizd
      
      if (icurrents.eq.76.or.icurrents.eq.77.or.icurrents.eq.14)
     &   write(90,645) 
      if ( (iwind.eq.11).or.(iwind.eq.12) ) write(90,703)

      if(icrn.eq.1) then
         write(6,770)
         write(90,770)
      endif

c--------------------------------------------------------------------
c       Gravity Centre computation for the initial condition
c       in satellite detected oil slicks
c---------------------------------------------------------------------
 
       if(isat.eq.1) then

          ilast = 0
          do k=1,numspills

             sumx=0.d0
             sumy=0.d0
             ifirst = ilast + 1
             ilast = ilast + nppmsv(k)             

             do i=ifirst,ilast
c
                px(i) = xgrid(px(i))
                py(i) = ygrid(py(i))
                sumx=sumx+px(i)
                sumy=sumy+py(i)

             enddo

             xavg_multi(k) = sumx/nppmsv(k)
             yavg_multi(k) = sumy/nppmsv(k)

          enddo

          nst0 = iage / delt   !no of timesteps for aging the oil properties

       endif

c----------------------------------------------------------------------
c     For case when region is whole mediterranean, print grid details 

c     for the selected subregion +/- 'dist' from spill.
c     subroutine setcst selectes coastal segments within the subregion 
c----------------------------------------------------------------------

       if(iregn.eq.0) then

          write(90,*) ' '
          write(90,*) 'Grid details in selected region'
          write(90,750) mmax,nmax
          write(90,751) along1,along2,alatg1,alatg2
          write(90,*) ' '
          write(90,*) 'Land/water mask:'

          do n=nmax,1,-1
             write(90,'(400i1)') (itype(m,n),m=1,mmax)
          enddo

          write(90,*) ' '
          write(90,*) 'Bathymetry in selected region'

          do n=nmax,1,-1
             write(90,'(400i5)') (int(h(m,n)),m=1,mmax)
          enddo

       endif

  750  format(' Dimensions: ',2i5)
  751  format(' Longitude limits:',2f9.5/' Latitude limits: ',2f9.5)

c----------------------------------------------------------------------
c     convert spill volume to barrels
c     rbm3 = cu m per barrel, gpm3 = gals per cu m.
c	ix=random 6-digit seed for random number generator
c----------------------------------------------------------------------

       rbm3=0.158987d0
       gpm3=219.970d0

       if(vlunit.eq.'bbls') volfac=1.d0
       if(vlunit.eq.'cu.m') volfac=1.d0/rbm3
       if(vlunit.eq.'tons') volfac=1.d0/(rbm3*deno)
       if(vlunit.eq.'gals') volfac=1.d0/(rbm3*gpm3)
 
       call seedmedslik(ix)

c----------------------------------------------------------------------
c
c     volfac = factor to convert spill volume to barrels
c     totbbl = total no of barrels spilt
c     totbblv(k) = total no of barrels spilt in the kth spill
c
c     ntot = initial value of total no of lagrangian parcels spilt read 
c            from par-file
c     nppms = no of parcels per mini-spill
c     nppmsv(k) = no of parcels of the kth spill
c     ntot1 = real total no of lagrangian parcels spilt calculated considering
c             nppmsv(k), totbblv(k), ntot and totbbl
c 
c     nmini = number of incremental mini-spills (one mini-spill per delt)
c             In polygon (single or multiple) case nmini=1 -> the code is not 
c             for polygons with a continuos release
c     bpp = bbls per lagrangian parcel
c     vmspl = volume per mini-spill (in cu m)
c     vmsplv(k) = volume of the kth spill (in cu m)
c
c     sscale=scale factor for spill printout (= 1 if totbbl < 500000)
c
c----------------------------------------------------------------------

      if (isat.eq.0) then

         if (numspills.gt.1) then   ! PNTM

             totbbl=0.d0

             do k=1,numspills
                splrtev(k) = splrtev(k) * volfac
                if(ispillv(k).gt.0) totbblv(k) = splrtev(k) * tspillv(k)
                if(ispillv(k).eq.0) totbblv(k) = splrtev(k)
                totbbl = totbbl + totbblv(k)
             enddo
c      
             ntot1 = 0

             do k=1,numspills
                nppmsv(k) = totbblv(k) * ntot / totbbl
                ntot1 = ntot1 + nppmsv(k)
             enddo

             ntot = ntot1

c            In the present implementation Multiple k Point sources
c            can be only instantaneus (tspill(k)=0)

             nmini=numspills
             bpp = totbbl/ntot

             do k=1,numspills
                bblpmsv(k) = nppmsv(k) * bpp
                vmsplv(k) = totbblv(k) * rbm3
             enddo

c
c     step for release of each minispill
c

             write(90,*) ' '
             do k=1,numspills
                print*, iddv(1),immv(1),iyrv(1)
                print*, iddv(k),immv(k),iyrv(k) 
                ddays = jdiff(iddv(1),immv(1),iyrv(1),iddv(k),immv(k),
     &                  iyrv(k))
c For leap years ----------------------------------------------------

                ly = 0
             
                if ( ( ((iyr/4)*4.eq.iyr).and.((iyr/100)*100.ne.iyr))
     &               .or.( (iyr/400)*400.eq.iyr) ) ly=1

                ndays1 = 365 + ly

                if( iyrv(k).gt.iyrv(1) ) ddays = ddays + ndays1
c-------------------------------------------------------------------

                dhrs = ddays * 24.d0 + tstartv(k) - tstartv(1)
                nstrel(k) = int( dhrs / delt + 1.49d0)
                print*, "=============================================="
                print*, "point: ", k
                print*, "ddays: ", ddays
                print*, "tspillv: ", tspillv(k)
                print*, "tstartv: ", tstartv(k)
                print*, "dhrs: ", dhrs
                print*, "nstrel: ", nstrel(k)
                print*," =============================================="
             enddo

         else     ! PNTS

             splrte=splrte*volfac
             if(ispill.gt.0) totbbl=splrte*tspill
             if(ispill.eq.0) totbbl=splrte
c      
             nmini=int(tspill / delt + 0.0001d0) + 1
             if(tspill.gt.(nmini-1)*delt) nmini=nmini+1
             vmspl=totbbl*rbm3/nmini
c
             nppms=int( dfloat(ntot)/dfloat(nmini) + 0.01d0 )
             ntot=nppms*nmini
             bpp=totbbl/ntot
             bblpms=nppms*bpp

         endif

      else if (isat.eq.1) then   ! PGNS / PGNM

         nmini = 1
         totbbl=0.d0
          
         do k=1,numspills

            splrtev(k) = splrtev(k) * volfac

            if(ispillv(k).eq.0) totbblv(k) = splrtev(k)

            totbbl = totbbl + totbblv(k) ! Total bbls (barrels) released

         enddo  

         bpp = totbbl/ntot ! barrels for each particle

         do k=1,numspills

            bblpmsv(k) = nppmsv(k) * bpp
            vmsplv(k) = totbblv(k) * rbm3

         enddo

      endif

      sscale=1.0d0
      i=int(1.d0+dlog10(totbbl/500000.d0))

      if (i.ge.1) sscale=10**i

      if ((isat.eq.0.and.numspills.gt.1).or.(isat.eq.1)) then

          write(90,*) 'Total bbls released = ',totbbl
          write(90,*) 'Total no of parcels = ',ntot
          write(90,*) 'Bbls per parcel     = ',bpp
          write(90,*) 'No of minispills    = ',nmini
          write(90,*) 'Parcels per m_spill = ',(nppmsv(k),k=1,
     &                                            numspills)
          write(90,*) 'Release step        = ',(nstrel(k),k=1,
     &                                            numspills)
          write(90,*) 'Spill duration      = ',(ispillv(k),k=1,
     &                                            numspills)
          write(90,*) 'Spill rate (bbls)   = ',(splrtev(k),k=1,
     &                                            numspills)
          write(90,*) 'bbls per minispill  = ',(bblpmsv(k),k=1,
     &                                            numspills)
            
      else

          write(90,*) 'Total bbls released = ',totbbl
          write(90,*) 'Total no of parcels = ',ntot
          write(90,*) 'Bbls per parcel     = ',bpp
          write(90,*) 'No of minispills    = ',nmini
          write(90,*) 'Parcels per m_spill = ',nppms

      endif

      
c----------------------------------------------------------------------
c     c1i = fraction of evaporative part
c     c2i = fraction of residual part,   c1i + c2i = 1
c     den1 = density of evaporative part (kg/m**3)
c     den2 = density of residual part (kg/m**3)
c     deno = oil density (kg/m**3):    deno = c1i*den1 + c2i*den2
c     tsat = saturation temp for evaporative component
c     fmaxe=max evaporative fraction
c
c     denk = coeff for change of density with evaporation
c     cdt, cvt = coeffs for change of density & viscosity with temperature
c     denw0 = sea water density at temp tk0 (kg/m3)
c     tk0 = reference temperature (deg Kelvin)
c----------------------------------------------------------------------

      deno=deno*1000.d0
      den2=den2*1000.d0
      c2i=respc/100.d0
      c1i=1.d0-c2i
      den1=(deno-c2i*den2)/c1i
      tsat=400.d0-3.d0*api
      fmaxe=c1i
      fmaxd=c2i

      denk=(den2-den1)/deno   !0.18
      cdt=0.008d0
      cvt=5000.d0
      denw0=1026.d0
      denw=1026.d0
      row=(denw0-deno)/deno
      tk0=293.d0

c----------------------------------------------------------------------
c
c     maxst = no of steps of computation of length delt each
c     ned = no of calls to evap/disp subroutine per delt
c     tedsec = time interval in secs between evap/disp calls
c     tiprs = time interval (hours) for printing output
c     nprs = step number at which output printing begins
c     iprs = no of steps between output printing
c
c----------------------------------------------------------------------

      maxst = tcomp*stph+0.001d0
      if (irestart.eq.1) maxst = maxst - ihrestart*stph+0.001d0
c
      ned=30
      tedsec=deltsc/ned
c
      tiprs=dfloat(iprs)
      nprs=tiprs*stph+0.001d0
      iprs=tiprs*stph+0.001d0

c----------------------------------------------------------------------
c
c	initialise wind forecast velocity
c
c----------------------------------------------------------------------

      do m=1,mm
         do n=1,nm
            winx(m,n)=0.d0
            winy(m,n)=0.d0
         end do
      end do               

c----------------------------------------------------------------------
c     if booms deployed, read boom data
c     bmtim(k) = hour of kth boom deployment relative to 0 hrs on nday 0
c     bmx1(k),bmy1(k),bmx2(k),bmy2(k) = 5 km grid coords of boom ends 
c     ibmeff(k) = efficiency (%) of kth boom 
c----------------------------------------------------------------------


c(M. De Dominicis, booms deployment has not been used and tested in the MEDSLIK-II_1.01)
c(D. Bruciaferri, booms deployment has not been used and tested in the MEDSLIK-II_1.02)

      nbooms=0


      if(ibooms.eq.1) then
         open(48,file='medslik.bms',status='old')
         read(48,'(a80)') empty
         read(48,'(i2)') nbooms
         read(48,'(a80)') empty
         do k=1,nbooms
            read(48,'(i2,1x,i2,1x,i4,2x,i2,3x,i5,4(3x,i2,2x,f5.2),
     &           2x,i3)') id,im,iy,ihb,lenbm,latdg1,alatm1,londg1,
     &           alonm1,latdg2,alatm2,londg2,alonm2,ibmeff(k)
          
            nday = jdiff(idd,imm,iyr,id,im,iy)
            bmtim(k)=nday*24.+ihb
            bmlat1(k)=dfloat(latdg1)+alatm1/60.d0
            bmlon1(k)=dfloat(londg1)+alonm1/60.d0
            bmlat2(k)=dfloat(latdg2)+alatm2/60.d0
            bmlon2(k)=dfloat(londg2)+alonm2/60.d0

            bmx1(k)=xgrid(bmlon1(k))
            bmy1(k)=ygrid(bmlat1(k))
            bmx2(k)=xgrid(bmlon2(k))
            bmy2(k)=ygrid(bmlat2(k))
         enddo    
         close(48)
      end if

c----------------------------------------------------------------------
c     If icrn = 1, read spill correction data and construct correction 
c     velocities in m/s to shift centre of slick from the previous
c     computed position (xold,yold) to the observed position (xnew,ynew)
c     For 1st correction, t0avg = mean time of oil release
c----------------------------------------------------------------------

      if(icrn.eq.1) then

         open(46,file='medslik.crn',status='old')
         read(46,'(a80)') empty
         read(46,'(i2)') ncrns

         do k=1,ncrns

            read(46,'(a80)') empty 
            read(46,'(i4)') kkk
            crntim(k)=dfloat(kkk)
            read(46,'(i2)') nvert 

            do j=1,nvert+1
               read(46,'(a80)') empty
            end do
          
            read(46,'(f6.3,3x,f6.3)') alat,alon
            xnew=xgrid(alon)
            ynew=ygrid(alat)
                   
            read(46,'(f6.3,3x,f6.3)') alat,alon
            xold=xgrid(alon)
            yold=ygrid(alat)

            deltax=(xnew-xold)*delx
            deltay=(ynew-yold)*dely

            if(k.eq.1) then

               t0avg = tspill/2.d0

               if(crntim(1).lt.tspill) t0avg = crntim(1)/2.d0
                 deltat = (crntim(1) - t0avg)*3600.d0

            endif

            if(k.ge.2) deltat = (crntim(k) - crntim(k-1))*3600.d0

            if(deltat.ne.0) then
               crnu(k)=deltax/deltat
               crnv(k)=deltay/deltat
               if(k.gt.1) then
                  crnu(k)=crnu(k)+crnu(k-1)
                  crnv(k)=crnv(k)+crnv(k-1)
               end if
            else
               write(6,*) 'Two spill corrections for the same times'
               write(6,*) 'Times are: ',crntim(k)
               stop
            end if
            read(46,'(a80)') empty 
         end do
         close(46)
      end if  

c----------------------------------------------------------------------
c     ttn=thickness of thin slick (m) (= 10 microns)
c     ttki=initial thickness of thick slick (m) (= 2 cm)
c     afac=initial area ratio of thin to thick slick areas (afac=4-8)
c     compute initial areas, volumes & radii of thick and thin slicks
c----------------------------------------------------------------------

      ttn=0.00001d0
      ttki=0.02d0
      afac=4.0d0

c----------------------------------------------------------------------
c     INITIAL ASSIGNMENT of mini-spill PROPERTIES
c----------------------------------------------------------------------
c     tre = time from release
c     ttk/atk/vtk = thickness/area/volume of thick slick
c     ttn/atn/vtn = thickness/area/volume of thick slick
c     tt0/at0/vt0 = average thickness/total area/total volume of 2 slicks
c
c     den = average density of oil in minispill i
c     vis/visem = viscosity of oil/mousse in minispill i
c     vtke/vtne/vte = evaporated volumes in thick/thin/total slick
c     vtkd/vtnd/vtd = dispersed volumes in thick/thin/total slick
c     ftk/ftn = evaporated fractions
c
c     fw = water fraction in mousse
c     xcl/xcs = volumes of large/small oil droplets dispersed below slick
c     pcte/pctd = percentages evaporated/dispersed
c
c----------------------------------------------------------------------

      write(90,*) ' '

      if (isat.eq.0) then

         do i=1,nmini

            if (numspills.gt.1) then   ! PNTM

                atk0=vmsplv(i)/(ttki+afac*ttn)
                atn0=atk0*afac
                vtk0=atk0*ttki
                vtn0=atn0*ttn
                rtk0v(i)=dsqrt(atk0/pi)
                rtn0v(i)=dsqrt((atn0+atk0)/pi)
                write(90,*) 'initial mini-slick radii = ',
     &                    rtk0v(i),rtn0v(i)

            else                       ! PNTS

                atk0=vmspl/(ttki+afac*ttn)
                atn0=atk0*afac
                vtk0=atk0*ttki
                vtn0=atn0*ttn
                rtk0=dsqrt(atk0/pi)
                rtn0=dsqrt((atn0+atk0)/pi)
                write(90,*) 'initial mini-slick radii = ',rtk0,rtn0

            endif

            write(90,*) ' '
            tre(i)=0.d0
            ttk(i)=ttki
            atk(i)=atk0
            atn(i)=atn0
            ato(i)=atk(i)+atn(i)
            vtk(i)=vtk0
            vtn(i)=vtn0
            vto(i)=vtk(i)+vtn(i)
            tto(i)=vto(i)/ato(i)
            den(i)=deno
            vis(i)=viso
            visem(i)=viso

            vtke(i) =0.d0
            vtne(i) =0.d0
            vte(i)  =0.d0
            vtkd(i) =0.d0
            vtnd(i) =0.d0
            vtd(i)  =0.d0
            ftk(i)  =0.d0
            ftn(i)  =0.d0
            fw(i)   =0.d0
            xcl(i)  =0.d0
            xcs(i)  =0.d0
            pcte(i) =0.d0
            pctd(i) =0.d0

            write(90,*) 'Initial minispill properties:'
            write(90,*) tre(i)
            write(90,*) ttk(i),atk(i),vtk(i)
            write(90,*) ttn   ,atn(i),vtn(i)
            write(90,*) tto(i),ato(i),vto(i)
            write(90,*) ftk(i),ftn(i),fw(i)
            write(90,*) ttk(i),atk(i),vtk(i)
            write(90,*) vtke(i),vtne(i),vte(i)
            write(90,*) pcte(i),pctd(i)

         enddo

      else         ! PGNS / PGNM

         do i=1,numspills

            atk0=vmsplv(i)/(ttki+afac*ttn)
            atn0=atk0*afac
            vtk0=atk0*ttki
            vtn0=atn0*ttn
            rtk0v(i)=dsqrt(atk0/pi)
            rtn0v(i)=dsqrt((atn0+atk0)/pi)
            write(90,*) 'Slick: ',i,', Initial radii: ',
     &                    rtk0v(i),rtn0v(i)

            write(90,*) ' '
            tre(i)=0.d0
            ttk(i)=ttki
            atk(i)=atk0
            atn(i)=atn0
            ato(i)=atk(i)+atn(i)
            vtk(i)=vtk0
            vtn(i)=vtn0
            vto(i)=vtk(i)+vtn(i)
            tto(i)=vto(i)/ato(i)
            den(i)=deno
            vis(i)=viso
            visem(i)=viso

            vtke(i) =0.d0
            vtne(i) =0.d0
            vte(i)  =0.d0
            vtkd(i) =0.d0
            vtnd(i) =0.d0
            vtd(i)  =0.d0
            ftk(i)  =0.d0
            ftn(i)  =0.d0
            fw(i)   =0.d0
            xcl(i)  =0.d0
            xcs(i)  =0.d0
            pcte(i) =0.d0
            pctd(i) =0.d0

            write(90,*) 'Initial slick ',i,' properties:'
            write(90,*) tre(i)
            write(90,*) ttk(i),atk(i),vtk(i)
            write(90,*) ttn   ,atn(i),vtn(i)
            write(90,*) tto(i),ato(i),vto(i)
            write(90,*) ftk(i),ftn(i),fw(i)
            write(90,*) ttk(i),atk(i),vtk(i)
            write(90,*) vtke(i),vtne(i),vte(i)
            write(90,*) pcte(i),pctd(i)

         enddo

      endif

c-----------------------------------------------------------------------
c    INITIAL ASSIGNMENT of STATUS PARAMETER "is" for each parcel 
c-----------------------------------------------------------------------
c
c        is=0 parcel not released
c
c        is=1 in the spreading surface slick 
c
c        is=2 on surface but not spreading
c
c        is=3 dispersed into water column
c
c        is=-nsg stuck on shore segment number nsg
c
c        is=9 beyond boundary limit
c
c        ib(i) = k or -k indicates parcel is stuck on kth boom
c
c        c1p(i) = bbls of evaporative component left in ith parcel
c
c        c2p(i) = bbls of non-evaporative component left in ith parcel
c
c        px(i), py(i) = horizontal grid coordinates of ith parcel
c
c        pz(i) = vertical sigma-coordinate of ith parcel (0/1 on bottom/surface)
c
c----------------------------------------------------------------------

      if (isat.eq.1) then

         do i=1,ntot
            is(i)=1
            ib(i)=0
            c1p(i)=c1i*bpp
            c2p(i)=(1.d0-c1i)*bpp
            pz(i)=1.d0
         enddo

      elseif(isat.eq.0) then 

         do i=1,ntot
            is(i)=0 
            ib(i)=0
            c1p(i)=c1i*bpp
            c2p(i)=(1.d0-c1i)*bpp
         enddo

         if (numspills.gt.1) then

            ilast = 0
            do k=1,numspills
         
               ifirst = ilast + 1
               ilast = ilast + nppmsv(k)
         
               do i=ifirst,ilast
                  px(i) = x0v(k)
                  py(i) = y0v(k)
                  pz(i) = 1.d0
               enddo
      
            enddo   
         
         else

            do i=1,ntot

               px(i)=x0
               py(i)=y0
               pz(i)=1.d0

            enddo

         endif

      endif

c----------------------------------------------------------------------
c    open file for fate parameters and write headings
c----------------------------------------------------------------------

      if (irestart.eq.0) then

         open(99,file='medslik.fte')

         if (isat.eq.0.and.numspills.gt.1) then

             write(99,670) (iddv(k),immv(k),iyrv(k),istartv(k),
     &                      k=1,numspills)
             write(99,671) (splatv(k),splonv(k),k=1,numspills)
             write(99,672) (x0v(k),y0v(k),k=1,numspills)
             write(99,673) (splrtev(k),vlunit,k=1,numspills)
             write(99,674) tspill,icomp
             write(99,680)
  670        format(' Dates:hours   ',20('     ',i2,'/',i2,'/',i4,' : ',
     &                i4.4))
  671        format(' Location:',20('    lat  ',f8.4,'  lon ',f8.4))
  672        format(' Grid coords:',20('    x  ',f8.4,'   y ',f8.4))
  673        format(' Total volumes:  ',20(f9.2,' ',a4,'   '))
  674        format(' Spill durations: ',f7.0,' hours'/
     &         ' Length of run: ',i4,' hours')

         else

             write(99,600) idd,imm,iyr,istart
             write(99,601) splat,splon,x0,y0,tspill,icomp 
             if(ispill.gt.0) write(99,602) splrte,vlunit
             if(ispill.eq.0) write(99,603) splrte,vlunit
             write(99,680)
                
         endif

680      format('   time   vol spilt  %evap    %srf   %srftot   ',
     &              '%disp   %cstfxd   %csttot   visem1   visem2  ',
     &              'visoil1  visoil2  denoil1  denoil2   wfrac1  ',
     &              'wfrac2  volratio')
c    &         '   wfrac1  wfrac2  volratio denw')

      elseif (irestart.eq.1) then

         open(99,file='medslik.fte',access='APPEND')

      endif

c----------------------------------------------------------------------
c
c     npcl = number of parcels released so far
c
c     nspill = number of minispills which have occured
c
c     nseg = number of coastal segments in vicinity of spill
c
c     icur = current date number if icurrents = 9
c
c     totucrn,totvcrn = total observational corrections to u,v up to current time 
c
c     utot, vtot, nuv used to compute an average slick velocity for diagnostics
c
c----------------------------------------------------------------------

      npcl=0
      nspill = 0
      nseg=0
      icur=0
      totucrn=0.d0
      totvcrn=0.d0
c
      utot=0.d0
      vtot=0.d0
      nuv=0 

c----------------------------------------------------------------------
c     For hot restart read restart file
c----------------------------------------------------------------------

      if (irestart.eq.1) then

         if(iyr.ge.2000) write(ay,'(i2.2)') iyr-2000
         if(iyr.lt.2000) write(ay,'(i2.2)') iyr-1900
         write(am,'(i2.2)') imm
         write(ad,'(i2.2)') idd
         write(ah,'(i2.2)') istart/100
         write(a3,'(i3.3)') ihrestart
         open(98,file=ay//am//ad//ah//'_'//a3//'.rso',
     &           form='unformatted')

         read(98) jtime,npcl,nspill,pcevp,pcsrf,pcdsp,pccst,deno,viso

         if(jtime.ne.ihrestart) then

            write(6,*) 'Restart file time ',jtime,' does not agree ',
     &                 'with restart time in input file ',ihrestart
            stop
         endif 

         do i=1,npcl

            read(98) is(i),ib(i),c1p(i),c2p(i),alon,alat,pz(i),seep(i) 
            px(i) = xgrid(alon)
            py(i) = ygrid(alat)

         enddo

         do i=1,nspill

            read(98) den(i),vis(i),visem(i),tre(i),c1ms(i),c2ms(i),
     &      atn(i),atk(i),ato(i),ttk(i),tto(i),
     &      vtn(i),vtk(i),vto(i),xcl(i),xcs(i),
     &      vtne(i),vtke(i),vte(i),vtnd(i),vtkd(i),vtd(i),
     &      ftk(i),ftn(i),fw(i),pcte(i),pctd(i) 

         enddo

         nsgr=1
    1    continue
         read(98,end=2) ns0(nsgr),xx1,yy1,vcst(nsgr)
         segx(nsgr) = (xx1 - along1) / dlong + 1.d0
         segy(nsgr) = (yy1 - alatg1) / dlatg + 1.d0 
         nsgr = nsgr + 1
         go to 1 
    2    continue
         close(98)
         nsgr = nsgr - 1 
        
         write(90,*) 'Impacted coastal segments '
         write(90,*) 'Num impacted segments = ',nsgr
         write(90,*) '    No Old-num Numpcls   Seep'

         do ns=1,nsgr

            num=0

            do i=1,npcl
               if(is(i).eq.-ns0(ns)) num=num+1
            enddo

            write(90,'(3i7,f13.6)') ns,ns0(ns),num,vcst(ns)

         enddo     
        
         totseep=0.

         do ns=1,nsgr
            totseep =totseep + vcst(ns)
         enddo

         write(90,*) 'total beached = ',totseep
         write(90,*) ' '

         pccst1=100.*(totseep)/(npcl*bpp)
         write(90,*) 'percentages beached = ',pccst,pccst1
           
      endif

c----------------------------------------------------------------------
c      ********* BEGIN THE MAIN COMPUTATIONAL LOOP *********
c----------------------------------------------------------------------
c
c     nst = time step
c     timehr = time in hours from start of spill
c     tsec = time in seconds
c
c----------------------------------------------------------------------
      print*, '========================================================'

      write(6,161)

      if(irestart.eq.1) write(6,162) ihrestart
       
  161 format(//' Program initialization is now complete. '/
     &       ' Spill simulation will start now.'/)
  162 format(' Run will restart from ',i4,' hours after spill'//)

      print*, '========================================================'
c

c-------------------------------------------------
      do 94 nst=1,maxst
c-------------------------------------------------
       
         timehr=nst*delt
      
         if(irestart.eq.1) timehr = timehr + ihrestart

         tsec=timehr*3600.d0

c----------------------------------------------------------------------
c
c 1.  Select the area that contains all the parcels: (mzb,nzb)-(mzf,nzf) 
c     (used also for selection of coastal segments).
c
c
c       y  |
c          |
c          |
c          - nzf       ##########################
c          |           #                        #
c          |           #                        #
c          |           #                        #
c     ymax -           #       *D               #
c          |           #                        #
c          |           #                    *C  #
c          - ymcntr    #            +           #
c          |           #                        #
c          |           #              *B        #
c     ymin -           #  *A                    #
c          |           #                        #
c          |           #                        #
c          - nzb       ##########################
c          |
c          |
c          |
c          |                           
c          |          mzb                       mzf 
c          |-----------|--|---------|--------|---|-----------
c                        xmin    xmcntr     xmax             x
c
c
c	(xcntr,ycntr) is the centre of the grid (mzb,nzb)-(mzf,nzf).
c
c 2.  Select the area (mzb2,nzb2)-(mzf2,nzf2) for interpolation of forecast
c     data, assuming spill does not travel faster than 1.5 nm/hr on average
c     -- extending distance = 'dist' nm  from spill --
c
c----------------------------------------------------------------------
         if (isat.eq.0) then !a

            if (npcl.eq.0) then !b (START TIME)

c for multiple/single point source of spill at start time

                if (numspills.gt.1) then !c
        
                   xmin=x0v(1)
                   xmax=x0v(1)
                   ymin=y0v(1)
                   ymax=y0v(1)

                   do k=1,numspills
 
                      xavg_multi(k) = x0v(k)
                      yavg_multi(k) = y0v(k)
 
                      if (x0v(k).lt.xmin) xmin = x0v(k)
                      if (x0v(k).gt.xmax) xmax = x0v(k)
                      if (y0v(k).lt.ymin) ymin = y0v(k)
                      if (y0v(k).gt.ymax) ymax = y0v(k)

                   enddo
 
                else    !c
 
                   xavg=x0
                   yavg=y0

                   xmin=x0
                   xmax=x0
                   ymin=y0
                   ymax=y0

                endif !c

            else   !b (THE FOLLOWING TIME STEPS)

c for single/multiple point source of spill at the following time step

                xmin=10000.d0
                xmax=-10000.d0
                ymin=10000.d0
                ymax=-10000.d0

                do 34 i=1,npcl

                   if(px(i).lt.xmin) xmin=px(i)
                   if(px(i).gt.xmax) xmax=px(i)
                   if(py(i).lt.ymin) ymin=py(i)
                   if(py(i).gt.ymax) ymax=py(i)

   34           continue

            endif !b

         else !a

c for single/multiple areal source of spill (added by M. De Dominicis) 

            xmin=10000.d0
            xmax=-10000.d0
            ymin=10000.d0
            ymax=-10000.d0

            do 35 i=1,ntot
               if(px(i).lt.xmin) xmin=px(i)
               if(px(i).gt.xmax) xmax=px(i)
               if(py(i).lt.ymin) ymin=py(i)
               if(py(i).gt.ymax) ymax=py(i)
   35       continue

         endif !a

c 1. Select the area that contains all the parcels: (mzb,nzb)-(mzf,nzf)

         mzb=int(xmin)-2
         mzf=int(xmax+1)+2
         nzb=int(ymin)-2
         nzf=int(ymax+1)+2

         if(mzb.lt.1) mzb=1
         if(nzb.lt.1) nzb=1
         if(mzf.gt.mmax) mzf=mmax
         if(nzf.gt.nmax) nzf=nmax

c    (xcntr,ycntr) is the centre of the grid (mzb,nzb)-(mzf,nzf)

         xcntr=dfloat(mzb+mzf)/2.d0
         ycntr=dfloat(nzb+nzf)/2.d0

c--------------------------------------------------------------------------
c 2. Select the area (mzb2,nzb2)-(mzf2,nzf2) for interpolation of forecast
c    data, assuming spill does not travel faster than 1.5 nm/hr on average
c     -- extending distance = 'dist' nm  from spill --

         dist = 1.5 * (tcomp - timehr)
         mex = int( dist * 1849. / delx )
         nex = int( dist * 1849. / dely )

         mzb2 = int(xmin) - mex
         mzf2 = int(xmax+1) + mex
         nzb2 = int(ymin) - nex
         nzf2 = int(ymax+1) + nex

         if(mzb2.lt.1) mzb2=1
         if(nzb2.lt.1) nzb2=1
         if(mzf2.gt.mmax) mzf2=mmax
         if(nzf2.gt.nmax) nzf2=nmax

c----------------------------------------------------------------------
c     compute water currents
c
c     icurrents=14: currents directly read from MFS Copernicus format 
c                   daily netcdf output files
c     icurrents=76: currents directly read from MFS old format hourly 
c                   netcdf output files
c     icurrents=77: currents directly read from MFS Copernicus format 
c                   hourly netcdf output files
c     if ismag = 1, compute horiz diffusivity from Smagorinsky model
c         (done later for each minispill)
c----------------------------------------------------------------------
                    
        if (icurrents.ge.14) then

          call fcstcur_intrpl(nst,timehr+tstart,delt,nfcst,fcstfn,
     &       ifcstfn,fcsttim,icurrents,iregn,fcstcurdir,ktmx_cu,type_cu)

c............ Runge Kutta 4th order scheme ..............................

          if ( numscheme.eq.1 ) then

c            CURRENTS FIELD INTERPOLATED AT TIME= t + delt/2

             trunge = delt/2.d0
             call curr_runge(nst,timehr+tstart,trunge,delt,nfcst,fcstfn,
     &                       ifcstfn,fcsttim,icurrents,iregn,ktmx_cu,
     &                       type_cu)

c            CURRENTS FIELD INTERPOLATED AT TIME= t + delt

             trunge = delt
             call curr_runge(nst,timehr+tstart,trunge,delt,nfcst,fcstfn,
     &                       ifcstfn,fcsttim,icurrents,iregn,ktmx_cu,
     &                       type_cu)

          endif

c.........................................................................
        end if

c----------------------------------------------------------------------
c     Select water current field to be used for slick advective transport; 
c----------------------------------------------------------------------
        do m=1,mmax
           do n=1,nmax
              if (idepth.eq.0) then
                 uadv(m,n)=usrf(m,n)
                 vadv(m,n)=vsrf(m,n)
              else if(idepth.eq.1) then
                 uadv(m,n)=u10(m,n)
                 vadv(m,n)=v10(m,n)
              else if(idepth.eq.2) then
                 uadv(m,n)=u30(m,n)
                 vadv(m,n)=v30(m,n)
              else if(idepth.eq.3) then
                 uadv(m,n)=u120(m,n)
                 vadv(m,n)=v120(m,n)
              end if   
           end do
        end do

c----------------------------------------------------------------------
c   Select water current field to be used for slick advective transport 
c   for the Runge Kutta 4th order numerical scheme 
c----------------------------------------------------------------------
       
       if (numscheme.eq.1) then
          do m=1,mmax
             do n=1,nmax
                if (idepth.eq.0) then
                   uadv1(m,n)=ucur1(m,n,1)
                   vadv1(m,n)=vcur1(m,n,1)
                   uadv2(m,n)=ucur2(m,n,1)
                   vadv2(m,n)=vcur2(m,n,1)
                else if(idepth.eq.1) then
                   uadv1(m,n)=ucur1(m,n,2)
                   vadv1(m,n)=vcur1(m,n,2)
                   uadv2(m,n)=ucur2(m,n,2)
                   vadv2(m,n)=vcur2(m,n,2)
                else if(idepth.eq.2) then
                   uadv1(m,n)=ucur1(m,n,3)
                   vadv1(m,n)=vcur1(m,n,3)
                   uadv2(m,n)=ucur2(m,n,3)
                   vadv2(m,n)=vcur2(m,n,3)
                else if(idepth.eq.3) then
                   uadv1(m,n)=ucur1(m,n,4)
                   vadv1(m,n)=vcur1(m,n,4)
                   uadv2(m,n)=ucur2(m,n,4)
                   vadv2(m,n)=vcur2(m,n,4)
                end if
             end do
          end do
       endif
        
c----------------------------------------------------------------------
c     Compute wind speed and direction
c	(winx,winy) = wind velocity components (m/s)
c
c     Set evaporation rate parameter that depends on wind speed.
c----------------------------------------------------------------------

        if ( (isat.eq.0.and.numspills.gt.1).or.(isat.eq.1) ) then

           do k=1,numspills

               if ((iwind.eq.11).or.(iwind.eq.12)) then
              
                     call fcstwnd_intrpl(xavg_multi(k),yavg_multi(k),
     &                         nst,timehr+tstart,delt,nwfcst,wfcstfn,
     &                         iwfcstfn,wfcsttim,iwind,ktmx_wd,
     &                         wvel_vec(k),wdir_vec(k))
                                                                        
               end if
               ce1_vec(k)=akew*(wvel_vec(k)*3.6d0)**gamma
           enddo
           
        else

           if ((iwind.eq.11).or.(iwind.eq.12)) then 

                call fcstwnd_intrpl(xavg,yavg,nst,timehr+tstart,delt,
     &                              nwfcst,wfcstfn,iwfcstfn,wfcsttim,
     &                              iwind,ktmx_wd,wvel,wdir)
    
           end if
           ce1=akew*(wvel*3.6d0)**gamma 
        endif

c............ Runge Kutta 4th order scheme ........................

          if ( numscheme.eq.1 ) then

c            WINDS FIELD INTERPOLATED AT TIME= t + delt/2

             trunge = delt/2.d0
             call wind_runge(nst,timehr+tstart,trunge,delt,nwfcst,
     &                       wfcstfn,iwfcstfn,wfcsttim,iwind,ktmx_wd)

c            WINDS FIELD INTERPOLATED AT TIME= t + delt

             trunge = delt
             call wind_runge(nst,timehr+tstart,trunge,delt,nwfcst,
     &                       wfcstfn,iwfcstfn,wfcsttim,iwind,ktmx_wd)

          endif

c----------------------------------------------------------------------
c     Compute the wind-induced drift velocity for the oil parcels.
c     slick speed = alpha * wind speed at angle beta to right of wind.
c
c     If ibetared = 1, beta reduces with wind velocity, else const. 
c
c     Subtract a fraction wredfrac of the wind already incorporated in 
c     the forecast current in case when iwindred = 1.
c----------------------------------------------------------------------

        do m=1,mm
           do n=1,nm

              winu=winx(m,n)
              winv=winy(m,n)
              winmod=dsqrt((winu*winu)+(winv*winv))

              if(winmod.eq.0.d0) then
                beta=beta0
              else
                beta=beta0*(1.d0-ibetared*winmod/(winmod+halfspeed))
              end if 

              csbeta=dcos(beta/degrad)
              snbeta=dsin(beta/degrad)
               
              wdrftx(m,n)=alpha*( winu*csbeta+winv*snbeta)
              wdrfty(m,n)=alpha*(-winu*snbeta+winv*csbeta)
  
c..................... Runge Kutta wdrft ...........................

              if ( numscheme.eq.1 ) then
                 winu1=wx1(m,n)
                 winv1=wy1(m,n)
                 winu2=wx2(m,n)
                 winv2=wy2(m,n)
                 winmod1=dsqrt((winu1*winu1)+(winv1*winv1))
                 winmod2=dsqrt((winu2*winu2)+(winv2*winv2))

                 if (winmod1.eq.0.d0) then
                    betar=beta0
                 else
                    betar=beta0*(1.d0-ibetared*winmod1/
     &                          (winmod1+halfspeed))
                 end if
                 
                 csbetar=dcos(betar/degrad)
                 snbetar=dsin(betar/degrad)

                 wdrftx1(m,n)=alpha*( winu1*csbetar+winv1*snbetar)
                 wdrfty1(m,n)=alpha*(-winu1*snbetar+winv1*csbetar)

                 if (winmod2.eq.0.d0) then
                    betar=beta0
                 else
                    betar=beta0*(1.d0-ibetared*winmod2/
     &                          (winmod2+halfspeed))
                 end if

                 csbetar=dcos(betar/degrad)
                 snbetar=dsin(betar/degrad)

                 wdrftx2(m,n)=alpha*( winu2*csbetar+winv2*snbetar)
                 wdrfty2(m,n)=alpha*(-winu2*snbetar+winv2*csbetar)
              endif
c......................................................................
           end do
        end do
 
c======================================================================
c                    STOKES DRIFT VELOCITY
c======================================================================
c
c    istoke = 0:     No stokes drift calculation
c
c    istoke = 1:     Stokes drift velocity calculation using JONSWAP 
c                    spectrum parameterization and instantaneous wind 
c                    fields
c
c    istoke = 2:     Stokes drift velocity calculated using WAVE model
c                    outputs
c----------------------------------------------------------------------         

         write(6,*) '**********************************************'
         if (istoke.eq.0) then
            
            write(6,*) 'ISTOKE = 0: NO STOKE DRIFT CALCULATION'

         elseif(istoke.eq.2) then

            write(6,*) 'STOKE DRIFT CALCULATION from WAVE MODEL OUTPUT'

         elseif(istoke.eq.1) then

            wdirstoke=wdir
            write(6,*) 'STOKE DRIFT CALCULATION USING JONSWAP SPECTRUM'

         endif
         write(6,*) '**********************************************'
        
         if (istoke.eq.1) then !c

            if ( (isat.eq.0.and.numspills.gt.1).or.(isat.eq.1) )  then !b

               do l=1,numspills !1

                  call calcfetch(xavg_multi(l),yavg_multi(l),
     &                           wdirstoke,fetch)
                  xavg_lon=glon(xavg_multi(l))
                  yavg_lat=glat(yavg_multi(l))

                  grav=9.8
                  pi=4.*datan(1.d0)
         
                  m0=int(xavg_multi(l)+0.5d0)
                  n0=int(yavg_multi(l)+0.5d0)
                  wxx=winx(m0,n0)
                  wyy=winy(m0,n0)
   
                  wvel=dsqrt(wxx*wxx+wyy*wyy)
                   
                  wwdir=0.d0

                  if (wxx.eq.0.) then

                     wwdir=0.d0
                     if(wyy.gt.0) wwdir=180.d0
                     if(wyy.le.0) wwdir=0.d0

                  else

                     wwdir=datan(wyy/wxx) * degrad

                     if (wxx.lt.0.d0) wwdir=wwdir+180.d0
                     wwdir=270.d0-wwdir

                  endif

                  do m=1,mm !2

                     do n=1,nm !3

                        ww=dsqrt((winx(m,n)*winx(m,n))+
     &                           (winy(m,n)*winy(m,n))) !wind velocity module

                        if (ww.eq.0) then !a

                           stoke=0
                           stoku(m,n)=0
                           stokv(m,n)=0

                        else !a

                           gamma=3.3
                           alfa=0.076*((ww**2)/(fetch*grav))**(0.22)
                           fang_peak=22*((grav**2)/(ww*fetch))
     &                               **(0.3333333333)
                           freq_sp(1)=0.001

                           do k=2,700      ! to be changed

                              freq_sp(k)=freq_sp(k-1)+0.001
                              fang(k)=2*pi*freq_sp(k)
                              if (fang(k).ge.fang_peak) sigma=0.09
                              if (fang(k).lt.fang_peak) sigma=0.07
                              erre(k)=exp(-((fang(k)-fang_peak)**2)
     &                                   /(2*(sigma**2)*(fang_peak**2)))
                              sp_exp2(k)=exp(-1.25*(fang_peak/
     &                                       fang(k))**4)
                              sp_exp1(k)=alfa*grav**2*fang(k)**(-5)
                              spectra(k)=sp_exp1(k)*sp_exp2(k)*
     &                                    (gamma**erre(k))
                              wave_num(k)=fang(k)*fang(k)/grav   ! deep water
                              stoke_sp(k)=2*spectra(k)*fang(k)*
     &                                    wave_num(k)

                           enddo

                           hwave_tot=0
                           stoke_tot=0
 
                           do k=1,699

                              hwave_d(k)=0.001*2*pi*
     &                                   (spectra(k)+spectra(k+1))/2
                              stoke_d(k)=0.001*2*pi*(stoke_sp(k)+
     &                                   stoke_sp(k+1))/2
                              stoke_tot=stoke_tot+stoke_d(k)
                              hwave_tot=hwave_tot+hwave_d(k)

                           enddo

                           hwave=4*sqrt(hwave_tot)
                           stoke=stoke_tot

                           wwangle = (270.d0-wwdir) / degrad
                           cswwdir = dcos(wwangle)

                           snwwdir = dsin(wwangle)

                           stoku(m,n)=stoke*cswwdir
                           stokv(m,n)=stoke*snwwdir

                        endif !a

                     enddo !3

                  enddo !2

               enddo!1

            else !b

               call calcfetch(xavg,yavg,wdirstoke,fetch)
               xavg_lon=glon(xavg)
               yavg_lat=glat(yavg)

               grav=9.8
               pi=4.*datan(1.d0)
         
               m0=int(xavg+0.5d0)
               n0=int(yavg+0.5d0)
               wxx=winx(m0,n0)
               wyy=winy(m0,n0)
   
               wvel=dsqrt(wxx*wxx+wyy*wyy)
               wwdir=0.d0
 
               if (wxx.eq.0.) then

                  wwdir=0.d0
                  if (wyy.gt.0) wwdir=180.d0
                  if (wyy.le.0) wwdir=0.d0

               else

                  wwdir=datan(wyy/wxx) * degrad

                  if(wxx.lt.0.d0) wwdir=wwdir+180.d0
                  wwdir=270.d0-wwdir
         
               endif 
 
        
               do m=1,mm

                  do n=1,nm

                     ww=dsqrt((winx(m,n)*winx(m,n))+(winy(m,n)*
     &                         winy(m,n))) !wind velocity module
        
                     if (ww.eq.0) then !a

                         stoke=0
                         stoku(m,n)=0
                         stokv(m,n)=0

                     else !a

                         gamma=3.3
                         alfa=0.076*((ww**2)/(fetch*grav))**(0.22)
                         fang_peak=22*((grav**2)/(ww*fetch))
     &                             **(0.3333333333)
                         freq_sp(1)=0.001

                         do k=2,700      

                            freq_sp(k)=freq_sp(k-1)+0.001
                            fang(k)=2*pi*freq_sp(k)
                            if (fang(k).ge.fang_peak) sigma=0.09
                            if (fang(k).lt.fang_peak) sigma=0.07
                            erre(k)=exp(-((fang(k)-fang_peak)**2)
     &                                  /(2*(sigma**2)*(fang_peak**2)))
                            sp_exp2(k)=exp(-1.25*(fang_peak/fang(k))**4)
                            sp_exp1(k)=alfa*grav**2*fang(k)**(-5)
                            spectra(k)=sp_exp1(k)*sp_exp2(k)*
     &                                 (gamma**erre(k))
                            wave_num(k)=fang(k)*fang(k)/grav  ! deep water
                            stoke_sp(k)=2*spectra(k)*fang(k)*wave_num(k)

                         enddo

                         hwave_tot=0
                         stoke_tot=0

                         do k=1,699

                            hwave_d(k)=0.001*2*pi*
     &                                 (spectra(k)+spectra(k+1))/2
                            stoke_d(k)=0.001*2*pi*(stoke_sp(k)+
     &                                  stoke_sp(k+1))/2
                            stoke_tot=stoke_tot+stoke_d(k)
                            hwave_tot=hwave_tot+hwave_d(k)

                         enddo

                         hwave=4*sqrt(hwave_tot)
                         stoke=stoke_tot

                         wwangle = (270.d0-wwdir) / degrad
                         cswwdir = dcos(wwangle)
                         snwwdir = dsin(wwangle)

                         stoku(m,n)=stoke*cswwdir 
                         stokv(m,n)=stoke*snwwdir
 
                     endif !a

                  enddo

               enddo

            endif !b

         endif !c

         if (istoke.eq.2) then 

            if (iwave.eq.101) then

               call fcstwav_intrpl(nst,timehr+tstart,delt,nwavfcst,
     &                             wavfcstfn,iwavfcstfn,wavfcsttim,
     &                             iwave,ktmx_wv,type_wv)

c........ Runge Kutta 4th order scheme Stokes Drift velocities .......

               if ( numscheme.eq.1 ) then
c              WAVES FIELD INTERPOLATED AT TIME= t + delt/2
                  trunge = delt/2.d0
                  call wave_runge(nst,timehr+tstart,trunge,delt,
     &                 nwavfcst,wavfcstfn,iwavfcstfn,wavfcsttim,
     &                 iwave,iregn,ktmx_wv,type_wv)
c              WAVES FIELD INTERPOLATED AT TIME= t + delt
                  trunge = delt
                  call wave_runge(nst,timehr+tstart,trunge,delt,
     &                 nwavfcst,wavfcstfn,iwavfcstfn,wavfcsttim,
     &                 iwave,iregn,ktmx_wv,type_wv)
               endif
c.........................................................................

            endif

         endif
       
c--------------------------------------------------------------------------
c  Compute advection velocities for the Runge-Kutta 4th order numerical
c  scheme
c--------------------------------------------------------------------------

         if (numscheme.eq.1) then

            do m=1,mmax
               do n=1,nmax
                  if (istoke.eq.0) then
                     urng1(m,n)=uadv1(m,n)+wdrftx1(m,n)
                     vrng1(m,n)=vadv1(m,n)+wdrfty1(m,n)
                     urng2(m,n)=uadv2(m,n)+wdrftx2(m,n)
                     vrng2(m,n)=vadv2(m,n)+wdrfty2(m,n)
                  else if (istoke.ne.0) then
                     urng1(m,n)=uadv1(m,n)+usto1(m,n)+wdrftx1(m,n)
                     vrng1(m,n)=vadv1(m,n)+vsto1(m,n)+wdrfty1(m,n)
                     urng2(m,n)=uadv2(m,n)+usto2(m,n)+wdrftx2(m,n)
                     vrng2(m,n)=vadv2(m,n)+vsto2(m,n)+wdrfty2(m,n)
                  endif
               enddo
            enddo

         endif
c============================================================================
c
c  *) Compute Sea Surface Temperature
c
c     tc, tk = sea surface temp in deg centigrade and kelvin
c     if isst = 0 tc is obtaines from climatology
c     if isst = 1 tc = sst(m,n) is obtained in subroutine fcstcur
c     if isst = 8 tc is already read from file medslik.inp
c
c  *) Set Evap and Disp properties that depend on sst
c
c     denw = density of sea water at temp tk
c     viso = oil viscosity at temp tk; tvk0 = ref temp for initial viscosity
c
c  *) Re-initialize minispill oil properties that depend on temp
c
c  *) Write initial data into fate output file for time 0
c
c==============================================================================

         if ((isat.eq.0).and.(numspills.eq.1)) then

            m1=int(xavg+0.5d0)
            n1=int(yavg+0.5d0)
         
            if (isst.eq.1) then
                tc=sst(m1,n1)
            end if  
    
            tk=tc+273.d0

         else
           
            do k=1,numspills
          
               m1 = int(xavg_multi(k)+0.5d0)
               n1 = int(yavg_multi(k)+0.5d0)

               if (isst.eq.1) then

                  tc = sst(m1,n1)

               end if

               tkv(k)=tc+273.d0
 
            enddo

         endif 

         if (nst.eq.1.and.irestart.eq.0) then

            if (isat.eq.0.and.numspills.eq.1) then      ! PNTS
               fac = dexp(cvt*((1.d0/tk)-(1.d0/tvk0)))
               viso = viso * fac

               do i=1,nmini
                  den(i)=deno
                  vis(i)=viso
                  visem(i)=viso
               enddo

            else                                         ! PNTM / PGNS / PGNM

               do i=1,numspills
                  fac = dexp(cvt*((1.d0/tkv(i))-(1.d0/tvk0)))
                  v_swap = viso
                  v_swap = v_swap * fac
                  den(i) = deno
                  vis(i) = v_swap
                  visem(i) = v_swap
               enddo

            endif 

            timout=0.
            nbblr=0.

            if (ispill.eq.0) nbblr = totbbl

            pcevap   =0.d0
            pcsrf    =100.d0
            pcsrtot  =100.d0
            pcdsp    =0.d0
            pccstfxd =0.d0
            pccsttot =0.d0
            wfr      =0.d0
            volratio =1.d0

            write(99,'(f9.2,f9.0,6f9.4,6f9.2,2f8.3,f10.5)') 

     &      timehr,nbblr/volfac,pcevp,pcsrf,pcsrftot,pcdsp,
     &      pccstfxd,pccsttot,visem(1),viso,viso,viso,deno,
     &      deno,wfr,wfr,volratio

         endif

               
c	  denw=denw0*(1.-cdt*(tk-tk0))
c	  m0=x0+0.5
c	  n0=y0+0.5
c	  write(90,*) tc,wvel,wdir,winx(m0,n0),winy(m0,n0)
c	  if((nst/8)*8.eq.nst) write(90,978) timehr,tk,denw
c  978   format('time = ',f6.1,' hrs    water temp = ',f5.1,' deg K'/
c     &         '           water dens = ',f6.1)


c----------------------------------------------------------------------
c     Compute corrections to oil drift velocity from slick observation
c----------------------------------------------------------------------

         if(icrn.eq.1) then
            if(timehr.le.crntim(1)) k1=1
            do k=2,ncrns
            if(timehr.le.crntim(k).and.timehr.gt.crntim(k-1)) k1=k
            end do
            if(timehr.gt.crntim(ncrns)) k1=ncrns
            totucrn=crnu(k1)
            totvcrn=crnv(k1)
         end if

c----------------------------------------------------------------------
c     Subroutine coast reads coast segments and coastal types and 
c     selects those segments that lie within the spill region
c     
c     alngt(i) = length of coastal segment i in metres
c     prel(i) = prob of oil being released in interval delt
c     sfrac(i) = frac of oil seeping into segment i in delt
c
c     For restart locate beached parcels on their new coastal segments
c----------------------------------------------------------------------

        call coast(delt,seg,sfrac,prel,nseg,alngt,regn1,api,apicoeff)
        
        if(irestart.eq.1.and.nst.eq.1) then
          write(90,*) 'Reassign beached parcels to ',nseg,' new segmnts'
          do ns=1,nss
            seep(ns) = 0.d0
          enddo
          do i=1,npcl
            itmp(i) = 0
          enddo
          
          do n=1,nsgr
            dmin=100000.d0
            do ns=1,nseg
              dx = segx(n) - seg(ns,1)
              dy = segy(n) - seg(ns,2)
              d = dx*dx + dy*dy
              if(d.lt.dmin) then
                dmin = d
                nsmin = ns
              endif
            enddo
            seep(nsmin) = vcst(n)
c	write(90,'(i5,4d18.8)')n,segx(n),segy(n),seg(nsmin,1),seg(nsmin,2)
c 
            num=0
            do i=1,npcl
              if(is(i).eq.-ns0(n)) then
                itmp(i) = -nsmin
                num=num+1
              endif   
            enddo
c	      write(90,*) 'Num parcels: ',num,vcst(n)
          enddo
           
          do i=1,npcl
            if(is(i).lt.0) is(i) = itmp(i)
          enddo

          write(90,*)
          write(90,*) 'Newly impacted coastal segments '
          write(90,*) 'Segment Numpcls Seep'

          totseep = 0.d0
          do ns=1,nss
            if(seep(ns).gt.0.d0) then
              totseep =totseep + seep(ns)
              num=0
              do i=1,npcl
                if(is(i).eq.-ns) num=num+1
              enddo
              write(90,'(2i7,f13.6)') ns,num,seep(ns)
            endif
          enddo     
          write(90,*) 'total beached 1 = ',totseep
          write(90,*) ' '
        endif

c======================================================================
c
c    - release new mini-spill
c
c    - initially locate new parcels at spill site
c
c    - initialize totals
c
c=======================================================================

        npcl0=npcl

        if (nspill .lt. nmini) then

           if (isat.eq.0) then

              if (numspills.gt.1) then

                 do k=nspill+1,nmini

                    if (nst.eq.nstrel(k)) then

                        nspill = nspill + 1
                        npcl2 = npcl + nppmsv(k)

                        do i=npcl + 1,npcl2

                           is(i)=1
                           px(i)=x0v(k)
                           py(i)=y0v(k)

                        enddo 

                        npcl = npcl2

                    endif

                 enddo

              else

                 nspill = nspill + 1
                 npcl = npcl + nppms

                 do 99 i=npcl0+1,npcl

                       is(i) = 1
                       px(i) = x0
                       py(i) = y0

   99            continue

              endif

           elseif (isat.eq.1) then

              nppms = ntot
              nspill = nspill + 1 ! isat=1 case: nspill=nspill+1=1
              npcl = npcl + nppms !              nppms = ntot for polygons

              do 54 i = npcl0+1,npcl

                    is(i) = 1
 
   54         continue

           endif

        endif

        en1ps=0.d0
        en2ps=0.d0
        en3ps=0.d0
        en4ps=0.d0
        en5ps=0.d0
        volsrf=0.d0
        tvolsrf=0.d0

c-----------------------------------------------------------------------------
c     START mini-spill LOOP (for each parcel released with the mini-spill ns)
c
c     mini-spill ns contains parcels l1 through l2
c-----------------------------------------------------------------------------
c
c     tre(ns): time in hours since release of mini spill ns
c
c     nsps: no of parcels still spreading from mini spill ns
c
c     xcm, ycm: centre of spreading part of slick mini spill ns
c
c     xcmall, ycmall: centre of whole mini spill ns
c-----------------------------------------------------------------------------

       if (nst.eq.1.and.irestart.eq.0) then

          write(90,*) 
          write(90,*) 'Release of parcels:'

       endif

       if ( timehr.le.(tspill+delt).or.nst.eq.nstrel(nspill) ) 
     &     write(90,*) 'Time hrs = ',timehr
c
       if ((isat.eq.0.and.numspills.gt.1).or.(isat.eq.1)) then

          fracetot = 0.d0
          l2 = 0

       endif

       if (isat.eq.1) nloop = numspills
       if (isat.eq.0) nloop = nspill

       do 70 ns = 1,nloop

          if ((isat.eq.0.and.numspills.gt.1).or.(isat.eq.1)) then

             l1 = l2 + 1
             l2 = l2 + nppmsv(ns)
             tre(ns) = (nst - nstrel(ns)) * delt

          else

             l1=(ns-1)*nppms+1
             l2=ns*nppms
             tre(ns)=tre(ns)+delt

          endif
                   
          nsps=0
          c1ms(ns)=0.d0
          c2ms(ns)=0.d0
          xcm=0.d0
          ycm=0.d0
          nspsall = 0
          xcmall = 0.d0
          ycmall = 0.d0

          if (timehr.le.(tspill+delt).or.nst.eq.nstrel(nspill)) then

              write(90,*) '    Minispill no. ',ns
              write(90,*) '    Parcels ',l1,' to ',l2

          endif

c----------------------------------------------------------------------
c    - stop minispill spreading after time sprdmx
c
c    - compute centres of spreading part and whole mini-slick
c
c    - compute Smagorinsky diffusivity at centre of each mini-slick
c----------------------------------------------------------------------

          do 58 i=l1,l2

                 if (tre(ns).ge.sprdmx.and.is(i).eq.1) is(i)=2

                 if (is(i).eq.1) then
                     xcm=xcm+px(i)
                     ycm=ycm+py(i)
                     nsps=nsps+1
                 end if

                 xcmall = xcmall + px(i)
                 ycmall = ycmall + py(i)
                 nspsall = nspsall + 1
   58     continue

          if (nsps.ge.1) then
             xcm=xcm/nsps
             ycm=ycm/nsps
          else
             if (isat.eq.0.and.numspills.eq.1) then
                xcm = x0
                ycm = y0
             else
                xcm = x0v(ns)
                ycm = y0v(ns)
             endif
          endif

          if (nspsall.ge.1) then
             xcmall = xcmall / nspsall 
             ycmall = ycmall / nspsall 
          else
             if (isat.eq.0.and.numspills.eq.1) then
                xcmall = xavg
                ycmall = yavg
             else
                xcmall = xavg_multi(ns)
                ycmall = yavg_multi(ns)
             endif
          end if
c
          if (icurrents.ge.1.and.icurrents.le.80.and.ismag.eq.1) then

              call smag(xcmall, ycmall, usrf, vsrf, horizk)
              horizd=sqrt(6.d0*horizk*delt*3600.d0)

c	    if(mod(nst,4).eq.0) 
c     &        write(91,*) 'Smagorinsky Horiz Diffus: ',nst/2,ns,horizk

          end if   

c----------------------------------------------------------------------
c
c     - subroutine ed computes evaporating and dispersing fractions
c       and rate of mousse formation; it is called ned times per step
c
c     - evaporation, dispersion and emulsification by mckay
c
c     - pctd,probd = percent dispersed and probability of disperion
c
c     - pcte,frace = percent and fraction evaporated
c
c     - fw(ns) = fraction of water in minispill ns
c
c     - rratio = ratio of radii of thick slick before and after time step
c
c----------------------------------------------------------------------

          pctd0=pctd(ns)
          atkns0 = atk(ns)
          rtkns0 = dsqrt(atk(ns) / pi) 

          if (trackmode.eq.0) then 
             if ((isat.eq.0.and.numspills.gt.1).or.(isat.eq.1)) then

                do k=1,ned
 
                   call ed(tedsec,vmsplv(ns),den(ns),vis(ns),visem(ns),
     &                     ttk(ns),ttn,tto(ns),atk(ns),atn(ns),ato(ns),
     &                     vtk(ns),vtn(ns),vto(ns),ftk(ns),ftn(ns),ft,
     &                     fw(ns),vtke(ns),vtne(ns),vte(ns),vtkd(ns),
     &                     vtnd(ns),vtd(ns),xcl(ns),xcs(ns),pcte(ns),
     &                     pctd(ns),ce1_vec(ns),wvel_vec(ns))
                enddo
             else

                do k=1,ned

                   call ed(tedsec,vmspl,den(ns),vis(ns),visem(ns),
     &                     ttk(ns),ttn,tto(ns),atk(ns),atn(ns),ato(ns),
     &                     vtk(ns),vtn(ns),vto(ns),ftk(ns),ftn(ns),ft,
     &                     fw(ns),vtke(ns),vtne(ns),vte(ns),vtkd(ns),
     &                     vtnd(ns),vtd(ns),xcl(ns),xcs(ns),pcte(ns),
     &                     pctd(ns),ce1,vwel)

                enddo
             endif          
          endif
         
          if(isat.eq.1.and.nst.eq.nst0) then 

              write(90,*) 'Aged minispill properties:'
              write(90,*) tre(ns)
              write(90,*) ttk(ns),atk(ns),vtk(ns) 
              write(90,*) ttn   ,atn(ns),vtn(ns) 
              write(90,*) tto(ns),ato(ns),vto(ns) 
              write(90,*) ftk(ns),ftn(ns),fw(ns) 
              write(90,*) ttk(ns),atk(ns),vtk(ns) 
              write(90,*) vtke(ns),vtne(ns),vte(ns) 
              write(90,*) pcte(ns),pctd(ns) 
 
          endif

          
          probd=(pctd(ns)-pctd0)/(100.d0-pctd0)
          frace=pcte(ns)/100.d0
          rtkns = dsqrt(atk(ns) / pi) 
          rratio = rtkns / rtkns0

c----------------------------------------------------------------------
c
c    - Displace and transform the lagrangian parcels
c    - First evaporation and dispersion check
c 
c----------------------------------------------------------------------

          do 66 i=l1,l2

              c1p(i)=(c1i-frace)*bpp

              if(isat.eq.1.and.nst.le.nst0) go to 65

              if(is(i).eq.9) go to 66

              rrr=randmedslik(ix)

              if ((is(i).eq.1.or.is(i).eq.2).and.rrr.lt.probd) is(i)=3  

c----------------------------------------------------------------------
c
c  - COMPUTE SPREADING DISPLACEMENT (thick slick contains most of the oil
c    since mechanical spreading is stopped after sprdmx = 24 hrs)
c  - new parcels randomly & uniformly distributed within initial thick
c    slick radius rtk0. (Thin slick volume is initially v small)
c
c----------------------------------------------------------------------

              xds=0.d0
              yds=0.d0

              if (is(i).eq.1) then

                  if (i.le.npcl0) then
                     xds = (rratio - 1.d0) * (px(i) -xcm) 
                     yds = (rratio - 1.d0) * (py(i) -ycm) 

                  else

                     rnd=randmedslik(ix)
                     radd=rtk0*dsqrt(rnd)
                     phi=2.d0*pi*randmedslik(ix)
                     xds=radd*dcos(phi)
                     yds=radd*dsin(phi)

                  end if

              end if    

c----------------------------------------------------------------------
c
c     COMPUTE DIFFUSION displacement: (m,n) = nearest grid pt to parcel
c
c----------------------------------------------------------------------

              m=int(px(i)+0.5d0)
              n=int(py(i)+0.5d0)

              xdd=(2.d0*randmedslik(ix)-1.d0)*horizd
              ydd=(2.d0*randmedslik(ix)-1.d0)*horizd 

              zdd=0.d0

              if (is(i).eq.3) then

                 call intrpl0(px(i),py(i),h,hint)
                 dep=(1.d0-pz(i))*hint
                 zdd=(2.d0*randmedslik(ix)-1.d0)*
     &                vertd(dep,vertd1,vertd2,thermocl)
              end if  

c--------------------------------------------------------------------------
c
c     COMPUTE ADVECTIVE displacements of surface & dispersed parcels
c     (including wind drift & correction term from observation(s) of spill)
c
c     Bilinear interpolation for wind & stoke drift velocity (added by 
c     M. De Dominicis) 
c
c---------------------------------------------------------------------------
        
c Interpolation of uadv,vadv,wdrftx,wdrfty,stoku,stokv at current simulation time
c on each particle position. It is needed by both the numerical schemes.

              if (is(i).ne.3) then

                 call intrpl(px(i),py(i),uadv,itype,ui)
                 call intrpl(px(i),py(i),vadv,itype,vi)
                 call intrpl(px(i),py(i),wdrftx,itype,wui)
                 call intrpl(px(i),py(i),wdrfty,itype,wvi)

                 zdispl(i)=0.d0

                 if (istoke.eq.0) then

                        ui=ui+wui
                        vi=vi+wvi

                 elseif (istoke.ne.0) then

                        call intrpl(px(i),py(i),stoku,itype,sui)
                        call intrpl(px(i),py(i),stokv,itype,svi)

                        ui=ui+sui+wui
                        vi=vi+svi+wvi

                 endif

c******************************************************                

                 if (is(i).eq.1.or.is(i).eq.2) then

                    utot=utot+ui+totucrn
                    vtot=vtot+vi+totvcrn
                    nuv=nuv+1

                 endif

c******************************************************                
              elseif (is(i).eq.3) then

                 dep=(1.d0-pz(i))*hint      
                 zdispl(i)=dep              

                 if (dep.le.fcstdep1) then

                    ui=(usrf(m,n)*(fcstdep1-dep)+u10(m,n)*dep)/10.d0
                    vi=(vsrf(m,n)*(fcstdep1-dep)+v10(m,n)*dep)/10.d0

                 elseif (dep.gt.fcstdep1.and.dep.le.fcstdep2) then

                    if (hint.ge.fcstdep2) then
                       denom=fcstdep2-fcstdep1
                       fac1=(dep-fcstdep1)
                       fac2=(fcstdep2-dep) 
                       ui=(u30(m,n)*fac1+u10(m,n)*fac2)/denom
                       vi=(v30(m,n)*fac1+v10(m,n)*fac2)/denom
                    else
                       denom=hint-fcstdep1
                       fac1=(dep-fcstdep1)
                       fac2=(hint-dep) 
                       ui=(u30(m,n)*fac1+u10(m,n)*fac2)/denom
                       vi=(v30(m,n)*fac1+v10(m,n)*fac2)/denom
                    endif

                 elseif (dep.gt.fcstdep2.and.dep.le.fcstdep3) then

                    if (hint.ge.fcstdep3) then
                       denom=fcstdep3-fcstdep2
                       fac1=(dep-fcstdep2)
                       fac2=(fcstdep3-dep) 
                       ui=(u120(m,n)*fac1+u30(m,n)*fac2)/denom
                       vi=(v120(m,n)*fac1+v30(m,n)*fac2)/denom
                    elseif(hint.lt.fcstdep2) then
                       ui=u30(m,n)
                       vi=v30(m,n)
                    else
                       denom=hint-fcstdep2
                       fac1=(dep-fcstdep2)
                       fac2=(hint-dep) 
                       ui=(u120(m,n)*fac1+u30(m,n)*fac2)/denom
                       vi=(v120(m,n)*fac1+v30(m,n)*fac2)/denom
                    endif

                 elseif (dep.gt.fcstdep3) then

                    ui=u120(m,n)
                    vi=v120(m,n)

                 end if

              end if

              if ( numscheme.eq.0 ) then

c                EULERO FORWARD NUMERICAL SCHEME IS USED

                 xdc=(ui+totucrn)*delt*3600.d0
                 ydc=(vi+totvcrn)*delt*3600.d0

              elseif ( numscheme.eq.1 ) then

c                FOURTH ORDER RUNGE-KUTTA NUMERICAL SCHEME IS USED

                 call runge(is(i),px(i),py(i),pz(i),hint,fcstdep1,
     &                      fcstdep2,fcstdep3,delt,ui,vi,xdc,ydc)

              endif

c----------------------------------------------------------------------
c
c     - displace parcels on surface and dispersed parcels
c
c     - (m21,n21) = old grid coordinates on bathy grid (nearest grid point)
c
c     - (ppx,ppy) = new grid coordinates of parcel i
c
c----------------------------------------------------------------------

              m21=int(px(i)+0.5d0)
              n21=int(py(i)+0.5d0)
             
              if (trackmode.eq.0) then
                 xdispl = (xds+xdc+xdd)
                 ydispl = (yds+ydc+ydd)
              else
                 xdispl = (xdc+xdd)
                 ydispl = (ydc+ydd)
              endif

              ppx = px(i) + xdispl/delx
              ppy = py(i) + ydispl/dely

c----------------------------------------------------------------------
c
c     check if beached parcels are released
c
c     sfrac = fraction of beached oil seeping into sand per step at 'ns'
c
c     prel = probability of beached oil being washed off at 'ns'
c
c     seep(nsg) = volume of oil at the impacted coastal segment 'nsg'
c
c     seepage slows when it approaches carrying capacity seepmx (bbls/km)
c
c     (xdispl,ydispl) must be to left of (dxseg,dyseg)
c
c----------------------------------------------------------------------

              if (is(i).lt.0.and.ib(i).eq.0) then

                 nsg = -is(i)
                 seepge = c2p(i)*sfrac(nsg)

                 fac = seep(nsg) / ( seepmx * alngt(nsg) )
                 if (fac.gt.10.d0) write(90,*) nst,i,fac

                 reducn = dexp( - fac )
                 seepge = seepge * reducn

                 c2p(i)=c2p(i)-seepge
                 seep(nsg)=seep(nsg)+seepge
                 probr=prel(nsg)

      
                 dxseg = seg(nsg,3) - seg(nsg,1)
                 dyseg = seg(nsg,4) - seg(nsg,2)
                 cross = dxseg * (ppy - py(i)) - dyseg * (ppx - px(i))

                 if ((randmedslik(ix).lt.probr).and. cross.gt.0.d0) then

                    pz(i)=1.d0
                    is(i)=2

                 end if

              end if

c----------------------------------------------------------------------
c     check if parcels stuck on booms are released
c----------------------------------------------------------------------
c
c            if((is(i).eq.1.or.is(i).eq.2).and.ib(i).ne.0) then
c              k=iabs(ib(i))
c              dxbm = bmx2(k) - bmx1(k)            
c              dybm = bmy2(k) - bmy1(k)
c              cross = dxbm * (ppy - py(i)) - dybm * (ppx - px(i))            
c              if(cross*ib(i).gt.0) ib(i)=0
c            end if
c----------------------------------------------------------------------
c     check if surface parcel hits one of the booms 
c     ib(i) = +k indicates parcel is on right side of boom k
c     ib(i) = -k indicates parcel is on left side of boom k
c----------------------------------------------------------------------

              if ((is(i).eq.1.or.is(i).eq.2).and.ib(i).eq.0) then

                   xi1 = px(i)            
                   eta1 = py(i)
                   xi2 = ppx
                   eta2 = ppy
                   dmin = 100.d0
                   nbmin=0
                   ibmin=0

                   do 64 k=1,nbooms
 
                      if (timehr+tstart.lt.bmtim(k)) go to 64

                      xx1 = bmx1(k)
                      yy1 = bmy1(k)
                      xx2 = bmx2(k)
                      yy2 = bmy2(k)
                      ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
                      ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
                      ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)

                      if (ddel.eq.0.d0) go to 64
                      if ( 100*randmedslik(ix).gt.ibmeff(k) ) go to 64
                 
                      alam=ddel1/ddel
                      alamp=ddel2/ddel

                      if (alam.ge.0.d0.and.alam.le.1.d0.and.
     &                    alamp.ge.0.d0.and.alamp.le.1.d0) then

                         xx=xx1+alam*(xx2-xx1)
                         yy=yy1+alam*(yy2-yy1)
                         dd1=dabs(xx-xi1)+dabs(yy-eta1)

                         if(dd1.lt.dmin) then

                            dmin=dd1
                            pppx=xx
                            pppy=yy
                            nbmin=k
                            if(ddel.gt.0.d0) ibmin=1
                            if(ddel.lt.0.d0) ibmin=-1

                         end if  

                      end if

   64              continue
c              
                   if (nbmin.ne.0) then

                      px(i)=pppx
                      py(i)=pppy
                      ib(i)=ibmin*nbmin

                   end if

              end if      

c----------------------------------------------------------------------
c     check if parcel hits any coastal segment
c     if parcel 'i' hits coastal segment  'ns', set is(i) = -ns
c----------------------------------------------------------------------
              if (is(i).gt.0.and.ib(i).eq.0) then

                xi1 = px(i)            
                eta1 = py(i)
                xi2 = ppx
                eta2 = ppy
                dmin = 100.d0
                nsmn=0
c            
                do 62 nsg = 1,nseg

                      xx1 = seg(nsg,1)
                      yy1 = seg(nsg,2)
                      xx2 = seg(nsg,3)
                      yy2 = seg(nsg,4)

                      ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
                      ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
                      ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)

                      if (ddel.eq.0.d0) go to 62
                 
                      alam=ddel1/ddel
                      alamp=ddel2/ddel

                      if (alam.ge.0.d0.and.alam.lt.1.d0.and.
     &                       alamp.ge.0.d0.and.alamp.lt.1.d0) then

                          xx=xx1+alam*(xx2-xx1)
                          yy=yy1+alam*(yy2-yy1)
c                          dd1=dabs(xx-xi1)+dabs(yy-eta1)
                          dd1=dsqrt( (xx-xi1)*(xx-xi1) + (yy-eta1)
     &                             *(yy-eta1) )
                          if (dd1.lt.dmin) then

                              dmin=dd1
                              pppx=xx
                              pppy=yy
                              nsmn=nsg

                          end if  

                      end if
   62           continue
              
                if (nsmn.ne.0) then

                   px(i)=pppx
                   py(i)=pppy  
                   is(i)=-nsmn

                else

                   px(i)=ppx
                   py(i)=ppy

                end if

            end if   

c----------------------------------------------------------------------
c
c     VERTICAL DIFFUSION OF DISPERSED PARCELS
c
c     ppz = temporary new sigma coordinate of parcel (1 at surface)
c
c     if ppz > 1 reflect displacement from surface
c
c     if ppz*caph < 0.2 (< 20 cm from bottom) parcel is sedimented
c
c----------------------------------------------------------------------

            if (is(i).eq.3) then

               caph=h(m21,n21)
               ppz=1.0d0
             
               if (caph.gt.0.2d0) ppz=pz(i)+zdd/caph

               if (ppz.gt.1.) ppz=2. - ppz

               if (ppz*caph.lt.0.2) then

                   ppz = 0.
                   is(i) = 4

               endif   

               pz(i)=ppz
               px(i)=ppx
               py(i)=ppy

            end if

c----------------------------------------------------------------------
c
c     count fractions of light and heavy components left in minispill ns
c
c----------------------------------------------------------------------

   65       continue
            c1ms(ns)=c1ms(ns)+c1p(i)
            c2ms(ns)=c2ms(ns)+c2p(i)
            
   66     continue

          if ((isat.eq.0.and.numspills.gt.1).or.(isat.eq.1)) then

             write(90,*) ns,c1ms(ns),bblpmsv(ns),c1ms(ns)/bblpmsv(ns)
             c1ms(ns)=c1ms(ns)/bblpmsv(ns)
             c2ms(ns)=c2ms(ns)/bblpmsv(ns)

          else

             c1ms(ns)=c1ms(ns)/bblpms
             c2ms(ns)=c2ms(ns)/bblpms

          endif

c----------------------------------------------------------------------
c
c  Compute totals for output:
c
c    - en1ps = vol of oil still spreading 
c    - en2ps = vol of oil on surface but no longer spreading 
c    - en3ps = vol of oil dispersed 
c    - en4ps = vol of oil on the coast but not permanently there
c    - volsrf = vol of oil-water mousse = (vol of oil)/(1-fw) 
c    - tvolsrf = total vol of oil-water mousse incl releasable oil on coast
c
c  Compute spill centre and shift to nearby al5-grid point
c----------------------------------------------------------------------

          do i=l1,l2

             if (is(i).eq.1) en1ps=en1ps+c1p(i)+c2p(i)
             if (is(i).eq.2) en2ps=en2ps+c1p(i)+c2p(i)
             if (is(i).eq.3) en3ps=en3ps+c1p(i)+c2p(i)
             if (is(i).lt.0) en4ps=en4ps+c1p(i)+c2p(i)
             if (is(i).eq.1.or.is(i).eq.2) 
     &           volsrf=volsrf+(c1p(i)+c2p(i))/(1.d0-fw(ns))
             if (is(i).ne.0.and.is(i).le.2) 
     &           tvolsrf=tvolsrf+(c1p(i)+c2p(i))/(1.d0-fw(ns))
          enddo

          if (isat.eq.0.and.numspills.gt.1) fracetot = fracetot + 
     &                                      frace * nppmsv(ns)
c		volsrf  =volsrf+(en1ps+en2ps)/(1.d0-fw(ns))
c		tvolsrf =tvolsrf+(en1ps+en2ps+en4ps)/(1.d0-fw(ns))
   70   continue

        if (isat.eq.0.and.numspills.gt.1) fracetot = fracetot / npcl

        volratio=tvolsrf/(npcl*bpp)

c----------------------------------------------------------------------
c                   GRAVITY CENTER COMPUTATION
c----------------------------------------------------------------------    

        if ( (isat.eq.0.and.numspills.gt.1).or.(isat.eq.1) ) then

            ilast = 0
            do k=1,numspills

               ifirst = ilast + 1
               ilast = ilast + nppmsv(k)
               sumx=0
               sumy=0

               do i=ifirst,ilast
                  sumx=sumx+px(i)
                  sumy=sumy+py(i)
               enddo

               xavg_multi(k)=sumx/nppmsv(k)
               yavg_multi(k)=sumy/nppmsv(k)
               xavg_lon_multi(k)=glon(xavg_multi(k))
               yavg_lat_multi(k)=glat(yavg_multi(k))

            enddo
        
        else ! for PNTS 

            sumx=0.d0
            sumy=0.d0

            do i=1,npcl

               sumx=sumx+px(i)
               sumy=sumy+py(i)

            end do

            xavg=sumx/npcl
            yavg=sumy/npcl
            xavg_lon=glon(xavg)
            yavg_lat=glat(yavg)

        endif

c----------------------------------------------------------------------
c
c     Print spill output. First summary statistics of fate percentages:
c
c     - ttseep = vol of oil permanently on the coast -> pccstfxd = %
c
c     - pcsrftot & tvolsrf include oil on coast but not permanently there
c     - pcsrf &  volsrf exclude such oil -> oil only on sea surface away 
c                                           from coasts
c     - pccsttot includes all oil on coast both permanent and temporary
c
c     - pccstfxd incudes only oil permanently attached to coast
c
c----------------------------------------------------------------------

        ttseep=0.d0

        do 72 nsg=1,nseg
           ttseep=ttseep+seep(nsg)
   72   continue
c
        c1tot=0.d0
        c2tot=0.d0

        if (isat.eq.0.and.numspills.gt.1) then

            do ns=1,nspill
           
                  c1tot=c1tot+c1ms(ns) * nppmsv(ns)
                  c2tot=c2tot+c2ms(ns) * nppmsv(ns)

            enddo

            c1tot=c1tot / dfloat(ntot)
            c2tot=c2tot / dfloat(ntot)
c
            pcevp    =100.d0* fracetot   !(c1i-c1tot)

        else         

            do ns=1,nspill

               c1tot=c1tot+c1ms(ns)
               c2tot=c2tot+c2ms(ns)

            enddo
            c1tot=c1tot/dfloat(nspill)
            c2tot=c2tot/dfloat(nspill)
c
            pcevp    =100.d0*(c1i-c1tot)
            
        endif

        pcsrftot =100.d0*(en1ps+en2ps+en4ps)/(npcl*bpp)
        pcsrf    =100.d0*(en1ps+en2ps)/(npcl*bpp)
        pcdsp    =100.d0*en3ps/(npcl*bpp)
        pccstfxd =100.d0*(ttseep)/(npcl*bpp)
        pccsttot =100.d0*(en4ps+ttseep)/(npcl*bpp)
        nbblr=int(npcl*bpp+0.5d0)

c----------------------------------------------------------------------
c
c    write fate parameters
c
c    time  volrel  %evap   %srf  %disp   %cst  visc1   visc2
c
c----------------------------------------------------------------------

        write(99,'(f9.2,f9.0,6f9.4,6f9.2,2f8.3,f10.5)') timehr,
     &     nbblr/volfac,pcevp,pcsrf,pcsrftot,pcdsp,pccstfxd,pccsttot,
     &     visem(1),visem(nmini),vis(1),vis(nmini),
     &     den(1),den(nmini),fw(1),fw(nmini),volratio

c----------------------------------------------------------------------
c
c     print spill distributions on coast, surface and dispersed in water
c
c     first open files and write headings
c
c----------------------------------------------------------------------

        if (nst.lt.nprs) go to 92

          jtime=timehr+0.001d0

          do 75 i=1,3

             write(a4,'(i4.4)') jtime
             fn(i)=pref//a4//a0(i)
             nfile=80+i

             open (nfile,file=fn(i))

             if (i.eq.1) write(nfile,801) vlunit  
             if (i.eq.2) write(nfile,802) vlunit   
             if (i.eq.3) write(nfile,803) vlunit 

             write(nfile,810) jtime,nbblr/volfac,vlunit,al5
             write(nfile,'(a27)')'  Spill_num  GC_lat  GC_lon'

             if ( (isat.eq.0.and.numspills.gt.1).or.(isat.eq.1) ) then
                do k=1,numspills
                  write(nfile,681) k,yavg_lat_multi(k),xavg_lon_multi(k)
                enddo
             else
                 ndefrec = 1
                 write(nfile,681) ndefrec,yavg_lat,xavg_lon
             endif

             write(nfile,811) pcevp,pcsrf,pcdsp,pccsttot,
     &                        volsrf/volfac,tvolsrf/volfac

                   
   75     continue   

  801     format(a4,' of oil on coast')
  802     format(a4,' of oil on surface')
  803     format(a4,' of dispersed oil')
  810     format(i4,'     : Hours after start of spill'/
     &           f8.0,' : ',a4,' of oil released so far'/
     &           f6.1,'   : Pixel size (m) for spill plotting'/)
  811  format('  %evap   %srf   %disp   %cst   srf_vol   Tot srf vol'/ 
     &           4f9.4,2e14.5)

c----------------------------------------------------------------------
c
c     count up and print barrels permanently stuck on coast
c
c     use the coastal segments to count & print parcel densities
c
c     vcst(nsg) = vol of oil per km at coastal segment 'nsg'
c
c----------------------------------------------------------------------

          kount=0        
          sum=0.d0

          do 76 nsg=1,nseg

             vcst(nsg)=0.d0
             if(seep(nsg).eq.0.d0) go to 76
             sum=sum+seep(nsg)
             dist=alngt(nsg)
             vcst(nsg)=seep(nsg)/dist
             kount=kount+1        

   76     continue
c
          write(81,839) kount
          write(81,840) vlunit

          do 79 nsg=1,nseg

              if(vcst(nsg).eq.0.) go to 79

              blon1=glon(seg(nsg,1))
              blat1=glat(seg(nsg,2))
              blon2=glon(seg(nsg,3))
              blat2=glat(seg(nsg,4))
              write(81,842) blat1,blon1,blat2,blon2,vcst(nsg)/volfac
       
   79     continue

          write(81,843) vlunit,sum/volfac
          
          kbooms=0

          do 80 k=1,nbooms

            if(timehr+tstart.lt.bmtim(k)) go to 80
            kbooms=kbooms+1

   80     continue  

          write(81,845) kbooms

          do 81 k=1,nbooms

            if(timehr+tstart.lt.bmtim(k)) go to 81
            write(81,846) bmlat1(k),bmlon1(k),bmlat2(k),bmlon2(k),
     &                    ibmeff(k)
   81     continue  

          close(81)

  839     format(i6,'     : Number of data points')
  840     format('   Lat      Lon      Lat      Lon        ',a4,'/km')

  841     format('    Lat      Lon         ',a4,'/km2')
  852     format('    Lat      Lon             m3/km2')
  842     format(4f11.7,e20.5)
  843     format('Total ',a4,' = ',f12.3)
  844     format(2f9.4,e20.5)
  845     format(i2,'       : Number of booms')
  846     format(4f9.4,i5)
  847     format(2f9.4)
  850     format('   Lat      Lon        u        v   ',
     &            2f9.6,'   dlat,dlon')
  851     format(2f12.6,2f9.4)
c----------------------------------------------------------------------
c     COMPUTATION OF CONCENTRATION
c
c     The grid with pixel size al5 is used to count & print
c
c     j=2: count up parcels on surface 
c          and print the corresponding barrels
c
c     j=3: count up parcels dispersed 
c          and print the corresponding barrels 
c
c----------------------------------------------------------------------

          do 89 j=2,3

             do m=1,npl

                do n=1,npl
 
                   poll(m,n)=0.d0
                   sum_avgdep(m,n)=0.d0
                   n_disp(m,n)=0.d0

                enddo

             enddo
c   
             sum=0.d0

             do 83 i=1,npcl

                if(is(i).eq.0) go to 83
                if(j.eq.2.and.is(i).gt.2) go to 83
                if(j.eq.3.and.is(i).ne.3) go to 83              

                mp=npl/2 + nint( (px(i)-xcntr) * rax )
                np=npl/2 + nint( (py(i)-ycntr) * ray )
 
                if(mp.lt.1) mp=1
                if(np.lt.1) np=1
                if(mp.gt.npl) mp=npl
                if(np.gt.npl) np=npl
                poll(mp,np)=poll(mp,np)+c1p(i)+c2p(i)
                sum_avgdep(mp,np)=sum_avgdep(mp,np)+zdispl(i)
                n_disp(mp,np)=n_disp(mp,np)+1
            
                sum=sum+c1p(i)+c2p(i)

   83        continue
  
             avgdep=sum_avgdep/n_disp

             kount=0

             do 84 m=1,npl
                do 84 n=1,npl

                   if (poll(m,n).eq.0.) go to 84
                   poll(m,n)=poll(m,n)/area
                   kount=kount+1

   84        continue
c
             nfile=80+j   
 
             write(nfile,839) kount
             write(nfile,841) vlunit
             kkount=0

             do 85 mp=1,npl
                do 85 np=1,npl

                   if(poll(mp,np).eq.0.d0) go to 85

                   xx = xcntr + (mp - npl/2) / rax
                   yy = ycntr + (np - npl/2) / ray

                   blon = glon(xx)
                   blat = glat(yy)

                   pden = poll(mp,np)

                   write(nfile,844) blat,blon,pden/volfac

   85           continue
                write(nfile,843) vlunit,sum/volfac

                write(nfile,845) kbooms

                do 87 k=1,nbooms

                      if(timehr+tstart.lt.bmtim(k)) go to 87
                      write(nfile,846) bmlat1(k),bmlon1(k),bmlat2(k),
     &                                 bmlon2(k),ibmeff(k)
   87            continue  
   
       
c----------------------------------------------------------------------
c     write wind vel and water currents in case of surface oil (j=2)
c----------------------------------------------------------------------

                 if (j.eq.2) then
                 
                     write(nfile,850) dlatg,dlong

                     mavg=int(xavg+0.5d0)
                     navg=int(yavg+0.5d0)
                     winx1=winx(mavg,navg)
                     winy1=winy(mavg,navg)
                     write(nfile,851) glat(yavg),glon(xavg),winx1,winy1

                     m1=mzb-5
                     m2=mzf+5
                     n1=nzb-5
                     n2=nzf+5

                     if(mod(m1,2).gt.0) m1 = m1 - mod(m1,2) 
                     if(mod(n1,2).gt.0) n1 = n1 - mod(n1,2) 

                     if(m1.lt.2) m1=2
                     if(n1.lt.2) n1=2
                     if(m2.gt.mmax) m2=mmax
                     if(n2.gt.nmax) n2=nmax

                     nstep=2
 
                     if(dlatg.gt.0.0125d0) nstep=1
 
                     do m=m1,m2,nstep
                        do n=n1,n2,nstep

                         write(nfile,851) glat(dfloat(n)),
     &                                    glon(dfloat(m)),
     &                                    uadv(m,n),vadv(m,n)
                        end do
                     end do

                 end if
    
                 close(nfile)
   89     continue
c----------------------------------------------------------------------
c     Trajectory output file - FINALIZING: write Gravity Centre data 
c----------------------------------------------------------------------

          if ( trackmode.eq.1) then

             write(93,990) nst*delt,yavg_lat,xavg_lon

          endif

 990      format(f5.1,4x,f10.4,2x,f10.4)
 681      format(i5,4x,f8.4,1x,f8.4)

          if(nst.gt.1) nprs=nprs+iprs 
c          endif
c		nprs=nprs+iprs
c
c----------------------------------------------------------------------

   92   continue

        jtime=nst/nstph
        if(jtime*nstph.eq.nst) then
           if(irestart.eq.0) then
              write(6,*) 'Simulation has now completed ',jtime,' hours'
           else 
              write(6,*) 'Simulation has now completed ',jtime,' hours',
     &               ' after restarting'   
           end if 
        end if 
        
c----------------------------------------------------------------------
c     next step
c----------------------------------------------------------------------
   94 continue

c**********************************************************************

      write(6,*) 'Simulation has completed successfully'

c**********************************************************************                

      if(nuv.ne.0) then
        write(6,*) 'Average velocity components of surface slick (m/s):'
        write(6,*) 'East ',utot/nuv,'   North ',vtot/nuv
      end if

c******************************************************                

      close(37)


c----------------------------------------------------------------------
c     Write restart file
c----------------------------------------------------------------------
      if(iyr.ge.2000) write(ay,'(i2.2)') iyr-2000
      if(iyr.lt.2000) write(ay,'(i2.2)') iyr-1900
      write(am,'(i2.2)') imm
      write(ad,'(i2.2)') idd
      write(ah,'(i2.2)') istart/100
      write(a3,'(i3.3)') jtime + ihrestart
      open(98,file=ay//am//ad//ah//'_'//a3//'.rso',form='unformatted')

      write(98) jtime+ihrestart,npcl,nspill,pcevp,pcsrf,pcdsp,pccst,
     &          deno,viso

      do i=1,npcl
          vcst(i)=0.d0
          if(is(i).lt.0) vcst(i) = seep(-is(i))
          alon = glon(px(i))
          alat = glat(py(i))
          write(98) is(i),ib(i),c1p(i),c2p(i),alon,alat,pz(i),vcst(i) 
      enddo

      do i=1,nspill
        write(98) den(i),vis(i),visem(i),tre(i),c1ms(i),c2ms(i),
     &    atn(i),atk(i),ato(i),ttk(i),tto(i),
     &    vtn(i),vtk(i),vto(i),xcl(i),xcs(i),

     &    vtne(i),vtke(i),vte(i),vtnd(i),vtkd(i),vtd(i),
     &    ftk(i),ftn(i),fw(i),pcte(i),pctd(i) 
      enddo

      do ns=1,nseg
         if(seep(ns).gt.0.d0) then
            xx1 = glon(seg(ns,1)) 
            yy1 = glat(seg(ns,2)) 
            write(98) ns,xx1,yy1,seep(ns)
         endif   
      enddo

      close(98)

      end subroutine main

c**********************************************************************
c
c                        LIST OF SUBROUTINES
c
c**********************************************************************
c
c------------------------
c  OIL FATE routines
c------------------------
c      subroutine setcst
c      subroutine coast
c      function hlseep
c      function hlwash
c      subroutine ed
c--------------------------
c  DIFFUSIVITY routines
c--------------------------
c      function vertd
c      subroutine smag
c--------------------------
c  WIND routines
c--------------------------
c      subroutine fcstwnd_intrpl 
c      subroutine readwind
c---------------------------
c  WATER CURRENT routines
c----------------------------
c      subroutine fcstcur_intrpl
c      subroutine readfcst
c----------------------------
c  WAVE routines
c----------------------------
c      subroutine fcstwav_intrpl
c      subroutine readstoke
c----------------------------
c  NUMERICAL SCHEME routines
c----------------------------
c      subroutine runge
c----------------------------
c  UTILITY routines
c----------------------------
c      subroutine interpol
c      subroutine intrpl
c      subroutine hsmoot
c      function julday
c      subroutine  date
c      subroutine merparam
c      subroutine ll2mer
c      subroutine mer2ll
c      function randmedslik
c      subroutine seedmedslik    
c      subroutine readsat  
c      subroutine calcfetch  
c
c**********************************************************************
c     OIL FATE routines
c**********************************************************************

c---------------------------------------------------------------------
      subroutine setcst(regn1)
c     constructs coastal segments from regional map 
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(ndim=1000, nss=200000)
c
      dimension bd1(nss,4),ibd(nss),
     &          cstlat(ndim),cstlon(ndim),icst(ndim)

      character regn1*4, empty*80
      logical ex

      common /blk2/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /encr/ aloncorr,iaencr,ibencr,icencr

      if(regn1.ne.'adri'.and.regn1.ne.'tyrr') then
        if(iaencr.eq.0) then
          open(51,file='data/'//regn1//'.map',status='old')

        else
          open(51,file='data/'//regn1//'_.map',status='old')
        endif 
      endif
      if (regn1.eq.'adri'.or.regn1.eq.'tyrr') then
         open(51,file='medf.map',status='old')
      endif

      open(52,file='data/'//regn1//'cst1.d')

      read(51,*) ncontours

      k=1
      do ni=1,ncontours
         read(51,*) isle
         if(iaencr.eq.0) then
            read(51,*) x1,y1
         else 
            read(51,*) lon
            read(51,*) lat
            x1 = dfloat( xor(lon, iaencr) ) / 1.d6 - aloncorr 
            y1 = dfloat( xor(lat, ibencr) ) / 1.d6
         endif  

         do i=2,isle
            if(iaencr.eq.0) then
               read(51,*) x2,y2
            else 
               read(51,*) lon
               read(51,*) lat
               x2 = dfloat( xor(lon, iaencr) ) / 1.d6 - aloncorr 
               y2 = dfloat( xor(lat, ibencr) ) / 1.d6
            endif  

            if( (x1.ge.along1.and.x1.le.along2.and.
     &          y1.ge.alatg1.and.y1.le.alatg2).or. 
     &         (x2.ge.along1.and.x2.le.along2.and.
     &          y2.ge.alatg1.and.y2.le.alatg2) ) then
                bd1(k,1) = x1 
                bd1(k,2) = y1 
                bd1(k,3) = x2 
                bd1(k,4) = y2

                k = k + 1
            endif           
            x1 = x2
            y1 = y2
         enddo 
      enddo

      nsegt = k-1
      write(90,*) 'No of coastal segments = ',nseg
      close(51)      

      do ns=1,nsegt
         ibd(ns) = 1
      enddo
c
c     read beach types
c
      ncstty = 0
      inquire(file='data/'//regn1//'.cst',EXIST=ex)
      if(ex) then
        write(6,*) 'Setting beach type data'
        write(90,*) 'Setting beach type data'
        open(80,file='data/'//regn1//'.cst')
        do k=1,11
          read(80,'(a80)') empty
        end do
        k = 0    
   95   continue
          read(80,*,end=96) clon,clat,ind
          k = k + 1
          cstlon(k) = clon
          cstlat(k) = clat
          icst(k) = ind
          go to 95
   96   continue 
        close(80)
        ncstty = k
c
c     set beach types for each segment
c
        do ns=1,nsegt
          ibd(ns) = 1
          x2 = (bd1(ns,1) + bd1(ns,3)) / 2.d0
          y2 = (bd1(ns,2) + bd1(ns,4)) / 2.d0
          dmin = 1000.d0
          do k=1,ncstty
            dist = dabs(x2 - cstlon(k)) + dabs(y2 - cstlat(k))
            if(dist.lt.dmin) then
              dmin=dist
              kmin=k
            end if  
          end do 
c          write(90,'(i5,2f8.4,i5,3f8.4,i5)') 
c     &     ns,alon,alat,kmin,dmin,cstlon(kmin),cstlat(kmin),icst(kmin)
          if(dmin.lt.0.01d0) ibd(ns) = icst(kmin)
        end do
      endif
c
c      output results
c
      write(52,'(''no. of segments = '',i6)') nsegt
      do ns=1,nsegt
        write(52,100) (bd1(ns,j),j=1,4),ibd(ns)
      end do
  100 format(4f11.7,i2) 
      rewind(52)
      
      return

      end subroutine setcst

c----------------------------------------------------------------------
      subroutine coast(dt,seg,sfrac,prel,nseg,alngt,regn1,api,apicoeff)
c     reads coastal segment data and selects those coastal segments 
c     that are close enough to spill
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=200000,msp=1200)
c	
      save ind,iused,ibd,bd1,nsegt,dkm
c
      integer*2 iused(nss)
      character regn1*4
      dimension bd1(nss,4),ibd(nss),dkm(nss),itype(mm,nm),
     &          seg(nss,4),prel(nss),sfrac(nss),alngt(nss)

      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk2/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

      data ind /0/

      xgrid(alon)=(alon - along1) / dlong + 1.d0
      ygrid(alat)=(alat - alatg1) / dlatg + 1.d0
      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg
 
c     coeff for reduction of coastal retention rate for heavy oils

      cseep = 1.d0
      if(api.lt.30.) cseep = 1. + (30.-api) * apicoeff
c
c   ....convert the boundary segments to grid coords....
c   ....ibd(n) = shore type of segment n (e.g. sand, rock, mangrove etc.)
c   ....dkm(n) = length of segment n in km.
c
      if(ind.eq.0) then

        call setcst(regn1)

        read(52,'(18x,i6)') nsegt
        do ns=1,nsegt
          read(52,100) (bd1(ns,j),j=1,4),ibd(ns)
          do j=1,3,2
            bd1(ns,j)   = xgrid(bd1(ns,j))
            bd1(ns,j+1) = ygrid(bd1(ns,j+1))
          enddo      

          dx = ( bd1(ns,3) - bd1(ns,1) ) * delx
          dy = ( bd1(ns,4) - bd1(ns,2) ) * dely
          dkm(ns) = dsqrt(dx*dx + dy*dy) / 1000.d0
        enddo
        close(52)
      end if
      ind=1
  100 format(4f11.7,i2) 
c
c   ....xb - xf, yb - yf: expected limits of spill....
c
      xb=mzb
      xf=mzf
      yb=nzb
      yf=nzf
c
c     select boundary segments that are close enough to spill site....
c     prel(j) = prob of release from coastal segment j in time step dt
c     sfrac(j) = fraction absorbed onto coast in time step dt
c     alngt(j) = length of segment j in km
c
      j=nseg
      do 20 ns=1,nsegt
        if(iused(ns).eq.1) go to 20
        x1=bd1(ns,1)
        y1=bd1(ns,2)
        x2=bd1(ns,3)
        y2=bd1(ns,4)
        if((x1.ge.xb.and.x1.le.xf.and.y1.gt.yb.and.y1.le.yf).or.
     &     (x2.ge.xb.and.x2.le.xf.and.y2.gt.yb.and.y2.le.yf)) then
          j=j+1
          if(j.gt.100000) then
            write(6,*) 'Too many coastal segments are impacted:' 
            write(6,*) '           reduce length of computation'
            stop
          end if

          iused(ns)=1  
          seg(j,1)=x1
          seg(j,2)=y1
          seg(j,3)=x2
          seg(j,4)=y2
c      write(90,'(4d25.20)') glon(seg(j,1)),glat(seg(j,2)),
c     &                      glon(seg(j,3)),glat(seg(j,4))
          ib=ibd(ns)
          tseep=hlseep(ib) * cseep
          twash=hlwash(ib)
          if(twash.eq.0.) then
            prel(j)=1.0d0
            sfrac(j)=0.0d0
          else   
            prel(j)  = 1.0d0 - 0.5d0**(dt/twash)
            sfrac(j) = 1.0d0 - 0.5d0**(dt/tseep)
          end if
          alngt(j)=dkm(ns)
        end if
   20 continue
      nseg=j  
  199 format(i3,4f10.3,i6,2f10.5)     
c

      return
      end  subroutine coast           

c----------------------------------------------------------------------
c   ....half-life for seepage into the beach; beach type = ib
c   ....coastal types ib:                   
c     1   sand beach                                                            
c     2   sand and gravel beach                                                 
c     3   cobble beach                                                          
c     4   rocky shore                                                           
c     5   seawall; concrete wharf etc                                                    
c     6   exposed headland                                                      
c     7   sheltered sand or gravel beach                                        
c     8   sheltered rocky shore                                                 
c     9   sheltered marsh or mud flats                                          
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      function hlseep(ib)
c----------------------------------------------------------------------
        implicit real*8(a-h,o-z)

        if(ib.eq.1) hlseep= 24.d0 
        if(ib.eq.2) hlseep= 36.d0
        if(ib.eq.3) hlseep= 48.d0
        if(ib.eq.4) hlseep= 96.d0 
        if(ib.eq.5) hlseep= 96.d0
        if(ib.eq.6) hlseep= 96.d0 
        if(ib.eq.7) hlseep= 24.d0
        if(ib.eq.8) hlseep= 96.d0
        if(ib.eq.9) hlseep= 24.d0

        return

      end function hlseep

c----------------------------------------------------------------------
c   ....half-life for washing off the beach; beach type = ib
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      function hlwash(ib)
c----------------------------------------------------------------------

        implicit real*8(a-h,o-z)

        if(ib.eq.1) hlwash= 24.d0
        if(ib.eq.2) hlwash= 24.d0
        if(ib.eq.3) hlwash= 24.d0
        if(ib.eq.4) hlwash= 18.d0
        if(ib.eq.5) hlwash=  0.d0
        if(ib.eq.6) hlwash=  1.d0 
        if(ib.eq.7) hlwash=120.d0
        if(ib.eq.8) hlwash=120.d0
        if(ib.eq.9) hlwash=120.d0
        return
      end function hlwash      

c-------------------------------------------------------------------------
      subroutine ed(dtim,vmspl,den,vis,visem,ttk,ttn,tto,
     &              atk,atn,ato,vtk,vtn,vto,ftk,ftn,ft,fw,
     &              vtke,vtne,vte,vtkd,vtnd,vtd,
     &              xcl,xcs,pcte,pctd,ce1,wsms)
c     Computes evaporation probability (eprob) and dispersion probability 
c     (dprob) and emulsification rate.
c
c     Program is modified from Mackay, Paterson and Trudel, 'Mathematical 
c     Model of Oil Spill Behaviour', Dec 1980
c-------------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      common /evap/ ce,vappr,fmaxe
      common /disp/ cd1,cd3,cd4,cd5,vl,vs1,um,stk,stn,fmaxd
      common /emul/ cm1,cm2,cm3,visemx
      common /sprd/ cs1,cs2,cs3
      common /phys/ deno,denk,cdt,viso,visk,cvt,tvk0,denw,tk,tk0

c----------------------------------------------------------------------
c     evaporation; ftk/ftn = fraction of oil evaporated in thick/thin slicks
c----------------------------------------------------------------------
      dvtke=0.d0
      dftk=0.d0
c
      if(vtk.gt.0.d0) then
        dextk=ce1*atk*dtim*(1.d0-ftk)/(tk*vtk)
        vpoil=vappr*exp(-ce*ftk)
        dftk=dextk*vpoil
        dftkmax=fmaxe-ftk
        if(dftk.gt.dftkmax/5.d0) dftk = dftkmax / 5.d0
        dvtke=vtk*dftk/(1.d0-ftk)
      end if
      vtke=vtke+dvtke
c
      dvtne=0.d0
      if(ftn.le.fmaxe) dvtne=vtn*(fmaxe-ftn)/(1.d0-ftn)
      ftn=fmaxe
      vtne=vtne+dvtne
      vte=vtke+vtne
      pcte=100.d0*vte/vmspl
      if(pcte.ge.fmaxe*100.d0) pcte=fmaxe*100.d0
c----------------------------------------------------------------------
c     check for dispersion if it exceeds max amount
c	initial/old values
c----------------------------------------------------------------------
   10 continue
      if(pctd.ge.fmaxd*100.d0) go to 30
      dvtkd=0.d0
      dvtnd=0.d0
      xclo=xcl
      xcso=xcs
      if(ttk.le.0.d0) go to 20
c----------------------------------------------------------------------
c     dispersion in thick slick
c     f = volume dispersed per per second per unit vol of slick
c     fb = fraction of small droplets
c     rbl = total vol of large drops dispersed per second
c     rbs = total vol of small drops dispersed per second
c----------------------------------------------------------------------
      f=cd3*(wsms+1.0)**2
      if(ttk.le.0.d0) go to 20
      fb=1.d0/(1.d0+cd4*(visem/10.d0)**(0.5d0)*(ttk/0.001d0)**1.d0*
     &                                           (stk/24.d0))
      rb=f*ttk
      rbt=rb*atk
      rbl=rb*(1.d0-fb)
      rbs=rb*fb

c----------------------------------------------------------------------
c     cl, xcl = concentration and amount of large drops dispersed
c     cs, xcs = concentration and amount of small drops dispersed
c     ct, xct = total concentration and amount dispersed  (ppm)
c----------------------------------------------------------------------

      cl=rbl/vl
      xcl=cl*um*atk
      cs=2.d0*rbs/(vs1+cd1)
      xcs=cs*um*atk
      ct=cl+cs
      xct=xcl+xcs
c----------------------------------------------------------------------
c     dxll = amount lost to lower layer per time step dtim
c     dxcl,dxcs = change in amount of large and small drops during dtim
c----------------------------------------------------------------------
      rd=0.5d0*(cd1-vs1)*cs
      rdt=rd*atk
      dxll=rdt*dtim
      dxcl=xcl-xclo
      dxcs=xcs-xcso
c      write(90,*) dxll,dxcl,dxcs
c----------------------------------------------------------------------
c     volume loss from spill: do not include change in large droplets
c	keep old concentrations
c----------------------------------------------------------------------
      dvtkd=dxll+dxcs !+dxcl
      vtkd=vtkd+dvtkd
      xclo=xcl
      xcso=xcs
c----------------------------------------------------------------------
c     dispersion in thin slick
c----------------------------------------------------------------------
   20 continue
      rtn=f*ttn/(1.d0+cd5*stn/24.d0)
      rtnt=rtn*atn
      dvtnd=rtnt*dtim
      vtnd=vtnd+dvtnd
c
c      ratd=rdt+rtnt
      vtd=vtkd+vtnd
      pctd=vtd*100.d0/vmspl
      dxm=dxll+dvtnd
c      write(90,*) dvtkd,vtkd,dvtnd,vtnd
c      write(90,*) vtd,vmspl,pctd
c      write(90,*) ' '
c----------------------------------------------------------------------
c     calculate mousse formation: fw = water fraction
c----------------------------------------------------------------------
   30 continue
      visr=exp(2.5d0*fw/(1.d0-cm1*fw))
      dfw=cm2*(wsms+1.d0)**2*(1.d0-cm3*fw)*dtim
      fw=fw+dfw
      if(fw.gt.1.d0/cm3) fw=1./cm3
      pcw=100.d0*fw
      visem=vis*visr
      if(visem.gt.visemx) visem=visemx
c----------------------------------------------------------------------
c     spreading: use Fay model for areas of both thick and thin slicks
c     calculate new volumes, areas, properties, etc.
c----------------------------------------------------------------------
   40 continue
      datns=0.
      dvtns=0.
      if(vtk.gt.0.) then
        datos=cs1*(ato**0.333d0)*exp(-cs3/(ttk+0.00001d0))*dtim
        dvtns=ttn*datos
        dvtks=-dvtns
        datks=dvtks/ttk + cs2*atk**0.333d0*ttk**1.333d0*dtim
        vtk=vtk-dvtke-dvtkd-dvtns
        atk=atk+datks
        ttk=vtk/atk
      end if
      vtn=vtn-dvtne-dvtnd+dvtns
c----------------------------------------------------------------------
c     transfer thick slick and droplet clouds to thin slick if ttk <= ttn
c----------------------------------------------------------------------
      if(ttk.le.ttn) then
        vtn=vtn+vtk+xcl+xcs
        vtk=0.d0
        atk=0.d0
        ttk=0.d0
        xcl=0.d0
        xcs=0.d0
      end if
      atn=vtn/ttn
c
      vto=vtn+vtk
      ato=atk+atn
      tto=vto/ato
c----------------------------------------------------------------------
c     calculate new compositions of slicks
c----------------------------------------------------------------------
      ftn=fmaxe
      ftk=ftk+dftk
      if(ftk.gt.fmaxe) ftk=fmaxe
      ft=(ftk*vtk+ftn*vtn)/(vtk+vtn)
      ftn=ftn-dvtns*(ftn-ftk)/vtn
c----------------------------------------------------------------------
c     Calculate new oil parameters. Note, denw & den do not contain 
c	temperature expansion effects so only the ratio den/denw is correct.
c	The ratio only is displayed on the output interface.
c     temperature effect on viscosity is included in main program
c----------------------------------------------------------------------
      if(ttk.gt.0.d0) then

c	  den=denw*fw+deno*(1.d0-fw)*(1.d0-cdt*(tk-tk0))*(1.d0+denk*ftk)
        den=denw*fw+deno*(1.d0-fw)*(1.d0+denk*ftk)
c        vis=viso*exp(visk*ftk+cvt*((1.d0/tk)-(1.d0/tvk0)))
        fac = exp(visk * ftk)
        vis = viso * fac
c        write(90,*) ftk,fac,viso,vis,visr,visem

      end if
c
      return
      end subroutine ed

c**********************************************************************
c     DIFFUSIVITY routines
c**********************************************************************

c----------------------------------------------------------------------
      function vertd(dep,vertd1,vertd2,thermocl)
c     Computes the mean diffusion step at depth 'dep'
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
c
      if(dep.lt.thermocl) then
        vertd=vertd1
      else
        vertd=vertd2
      end if

      return
      end function vertd            

c----------------------------------------------------------------------
      subroutine smag(x0,y0,us,vs,horizk)
c     compute horizontal diffusity from Smagorinsky model
c----------------------------------------------------------------------
      
      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370)
      dimension itype(mm,nm),us(mm,nm),vs(mm,nm)
      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad

      data c /0.1d0/

        m=int(x0+0.5d0)
        n=int(y0+0.5d0)

        if(itype(m,n).eq.0) then
          ux=0.d0
          vx=0.d0
        else if(itype(m+1,n).ne.0.and.itype(m-1,n).ne.0) then
          ux=(us(m+1,n)-us(m-1,n))/(2.d0*delx)
          vx=(vs(m+1,n)-vs(m-1,n))/(2.d0*delx)
        else if(itype(m+1,n).ne.0.and.itype(m-1,n).eq.0) then
          ux=(us(m+1,n)-us(m,n))/(delx)
          vx=(vs(m+1,n)-vs(m,n))/(delx)
        else if(itype(m+1,n).eq.0.and.itype(m-1,n).ne.0) then
          ux=(us(m,n)-us(m-1,n))/(delx)
          vx=(vs(m,n)-vs(m-1,n))/(delx)
        else 
          ux=0.d0 
          vx=0.d0
        end if

        if(itype(m,n).eq.0) then
          uy=0.d0
          uy=0.d0
        else if(itype(m,n+1).ne.0.and.itype(m,n-1).ne.0) then
          uy=(us(m,n+1)-us(m,n-1))/(2.d0*dely)
          vy=(vs(m,n+1)-vs(m,n-1))/(2.d0*dely)
        else if(itype(m,n+1).ne.0.and.itype(m,n-1).eq.0) then
          uy=(us(m,n+1)-us(m,n))/(dely)
          vy=(vs(m,n+1)-vs(m,n))/(dely)
        else if(itype(m,n+1).eq.0.and.itype(m,n-1).ne.0) then
          uy=(us(m,n)-us(m,n-1))/(dely)
          vy=(vs(m,n)-vs(m,n-1))/(dely)
        else 
          uy=0.d0
          vy=0.d0
        end if

        fac=ux*ux+0.5d0*(vx+uy)*(vx+uy)+vy*vy
        horizk=c*delx*dely*dsqrt(fac)

        if(horizk.lt.0.5d0) horizk=0.5d0
        if(horizk.gt.15.d0) horizk=15.d0

        return
      end subroutine smag

c**********************************************************************
c     WIND routines
c**********************************************************************

c---------------------------------------------------------------------
      subroutine fcstwnd_intrpl(xavg,yavg,nst,time,delt,nwfcst,wfcstfn,
     &                      iwfcstfn,wfcsttim,iwind,ktmx_wd,wvel,wdir)
c
c     constructs wind velocity field in case of forecast wind
c     (winx,winy) = forecast wind velocity at grid point (m,n)
c     (dwinx,dwiny) = increment per time step
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=100000,npl=2000,msp=1200)
      save isub,readdata,ifile,iwindrec!,dwinx,dwiny
c
c      dimension winx(mm,nm),winy(mm,nm),dwinx(mm,nm),dwiny(mm,nm)
      dimension iwfcstfn(30),wfcsttim(30,24),temp(24),itype(mm,nm)

      common /wind/ winx(mm,nm),winy(mm,nm),dwinx(mm,nm),dwiny(mm,nm)
      common /spill/ idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /blk1/  mmax,nmax,delx,dely,itype,pi,degrad!,ihcst
      common /blk2/  along1,alatg1,along2,alatg2,dlong,dlatg
      common /temp/  m0,n0

      character dirr*13,wfcstfn(30)*14,fn*14
      logical readdata

      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg

      data isub /1/
      if(isub.eq.0) return         ! end of available data files

      m0=int(xavg+0.5d0)           ! centre of slick
      n0=int(yavg+0.5d0)
      slon=glon(dfloat(m0))
      slat=glat(dfloat(n0))

      if((iwind.eq.11).or.(iwind.eq.12)) nrecs=ktmx_wd ! no of wind records in each daily file

      dtwind = 24.d0 / dfloat(nrecs)
      timew = time            

      frac = dtwind / delt  

      if((iwind.eq.11).or.(iwind.eq.12)) dirr='INP_DATA/MET/'
c
c     On 1st step set times of forecast records: wfcsttim(i,k) = time of 
c     the kth record in ith file measured from 0 hrs on spill date
c     Then set the initial wind forecast file and record
c
      if(nst.eq.1) then

          write(90,*) ' '
          write(90,*) 'Wind forecast files and times of wind records:'

          do i=1,nwfcst
             write(6,*) wfcstfn(i)
             read(wfcstfn(i)(5:6),'(i2)') iy
             read(wfcstfn(i)(7:8),'(i2)') im
             read(wfcstfn(i)(9:10),'(i2)') id

             ndaym1 = jdiff(idd,imm,iyr,id,im,iy+2000)
             do k = 1,nrecs
                wfcsttim(i,k) = ndaym1 * 24.d0 + (k-1) * dtwind
                temp(nrecs + 1 - k) = 24.d0 - wfcsttim(i,k)
             end do
             if(ihcst.eq.1) then
                do k = 1,nrecs
                   wfcsttim(i,k) = temp(k)
                end do
             endif
             write(90,'(a14,2x,24f6.0)')wfcstfn(i),
     &            (wfcsttim(i,k),k=1,nrecs)
          end do

          readdata=.true.
          ifile=1
          do k=1,nrecs
             if(timew.ge.wfcsttim(ifile,k)) iwindrec=k
          end do
      end if
c
c       read a fresh pair of wind records
c
      if(timew.ge.wfcsttim(ifile,iwindrec)) then

         readdata=.true.
         write(90,*) timew,wfcsttim(ifile,iwindrec)
         fac=(timew - wfcsttim(ifile,iwindrec))/delt
      endif

      if(readdata) then

         readdata=.false.
         fn=wfcstfn(ifile)
         iw = iwindrec
         if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
         write(6,*) 'Reading wind from file '//fn//', record ',iw
         write(90,*) 'Reading wind from file '//fn//', record ',iw

         call readwind(dirr,fn,iw,nrecs,winx,winy)

         if (iwindrec.eq.nrecs) then
              iwindrec=1
              ifile=ifile+1
         else
              iwindrec=iwindrec+1
         end if

         if(iwfcstfn(ifile).eq.0) then
              isub=0
              return
         end if

         fn=wfcstfn(ifile)
         iw = iwindrec
         if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
         write(6,*) 'Reading wind from file '//fn//', record ',iw
         write(90,*) 'Reading wind from file '//fn//', record ',iw
         call readwind(dirr,fn,iw,nrecs,dwinx,dwiny)
c
c
         do m=1,mmax
            do n=1,nmax
               dwinx(m,n) = (dwinx(m,n)-winx(m,n)) / frac
               dwiny(m,n) = (dwiny(m,n)-winy(m,n)) / frac
               winx(m,n) = winx(m,n) + dwinx(m,n) * fac
               winy(m,n) = winy(m,n) + dwiny(m,n) * fac
            end do
         end do

      else
c
c       increment the values for each re-computation of the spill
c
         do m=1,mmax
            do n=1,nmax
               winx(m,n)  = winx(m,n)  + dwinx(m,n)
               winy(m,n)  = winy(m,n)  + dwiny(m,n)
            end do
         end do

      end if
c
c       compute wind speed and direction at centre of spill
c

      wxx=winx(m0,n0)
      wyy=winy(m0,n0)
      wvel=dsqrt(wxx*wxx+wyy*wyy)
      wdir=0.d0
      if(wxx.eq.0.d0) then
          wdir=0.d0
          if(wyy.gt.0.d0) wdir=180.d0
          if(wyy.le.0.d0) wdir=0.d0
      else
          wdir=datan(wyy/wxx) * degrad
          if(wxx.lt.0.d0) wdir=wdir+180.d0
          wdir=270.d0-wdir
      end if
c
      return

      END SUBROUTINE fcstwnd_intrpl

c-------------------------------------------------------------------------

      SUBROUTINE wind_runge(nst,time,deltrng,delt,nfile,filename,
     &                             ifilename,fcsttim,ifield,nrecs)

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=100000,npl=2000,msp=1200)
      save readdata1,iwindrec1,readdata2,iwindrec2
      dimension ifilename(30),fcsttim(30,24),temp(24),itype(mm,nm)
      character fcst_dir*14,filename(30)*14, fn*16
      logical readdata,readdata1,readdata2
     
      common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

      common /rng1/ ucur1(mm,nm,4),vcur1(mm,nm,4),usto1(mm,nm),
     &              vsto1(mm,nm),ducr1(mm,nm,4),dvcr1(mm,nm,4),
     &              dust1(mm,nm),dvst1(mm,nm),wx1(mm,nm),wy1(mm,nm),
     &              dwx1(mm,nm),dwy1(mm,nm),irngc1,irngw1,irngwd1,
     &              isubc1,isubw1,isubwd1
      common /rng2/ ucur2(mm,nm,4),vcur2(mm,nm,4),usto2(mm,nm),
     &              vsto2(mm,nm),ducr2(mm,nm,4),dvcr2(mm,nm,4),
     &              dust2(mm,nm),dvst2(mm,nm),wx2(mm,nm),wy2(mm,nm),
     &              dwx2(mm,nm),dwy2(mm,nm),irngc2,irngw2,irngwd2,
     &              isubc2,isubw2,isubwd2

     
      if (nst.eq.1) then
         isub = 1
         readdata=.true.
         ifile=1
      else
         if (deltrng.eq.(delt/2.d0)) then
            isub = isubwd1
            ifile = irngwd1
            readdata = readdata1
            iwindrec = iwindrec1
         else if (deltrng.eq.delt) then
            isub = isubwd2
            ifile = irngwd2
            readdata = readdata2
            iwindrec = iwindrec2
         endif
      endif
          
      if(isub.eq.0) return         ! end of available data files

      dtwind = 24.d0 / dfloat(nrecs)
      timew = time + deltrng 
      frac = dtwind / delt   

      fcst_dir='INP_DATA/MET/'

      write(6,*) '    Linear time interpolation of wind'
      write(6,*) 'velocities for 4th order Runge-Kutta scheme'

c
c     On 1st step set times of forecast records: wfcsttim(i,k) = time of 
c     the kth record in ith file measured from 0 hrs on spill date
c     Then set the initial wind forecast file and record
c
      if(nst.eq.1) then
         do i=1,nfile
            write(6,*) filename(i)
            read(filename(i)(5:6),'(i2)') iy
            read(filename(i)(7:8),'(i2)') im
            read(filename(i)(9:10),'(i2)') id
            ndaym1 = jdiff(idd,imm,iyr,id,im,iy+2000)
            do k = 1,nrecs
               fcsttim(i,k) = ndaym1 * 24.d0 + (k-1) * dtwind
               temp(nrecs + 1 - k) = 24.d0 - fcsttim(i,k)
            end do
          end do
          do k=1,nrecs
             if(timew.ge.fcsttim(ifile,k)) iwindrec=k
          end do
      end if

c
c       read a fresh pair of wind records
c
      if(timew.ge.fcsttim(ifile,iwindrec)) then
         readdata=.true.
         fac=(timew - fcsttim(ifile,iwindrec))/delt
      endif

      if (readdata) then
         readdata=.false.
         fn=filename(ifile)
         iw = iwindrec
!         if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
         write(6,*) 'Reading wind from file '//fn//', record ',iw

         if (deltrng.eq.(delt/2.d0)) then
            call readwind(fcst_dir,fn,iw,nrecs,wx1,wy1)
         else if (deltrng.eq.delt) then
            call readwind(fcst_dir,fn,iw,nrecs,wx2,wy2)
         endif

         if (iwindrec.eq.nrecs) then
              iwindrec=1
              ifile=ifile+1
         else
              iwindrec=iwindrec+1
         end if

         if(ifilename(ifile).eq.0) then
              isub=0
              return
         end if

         fn=filename(ifile)
         iw = iwindrec
c         if(ihcst.eq.1) iw = nrecs + 1 - iwindrec
         write(6,*) 'Reading wind from file '//fn//', record ',iw

         if (deltrng.eq.(delt/2.d0)) then
            call readwind(fcst_dir,fn,iw,nrecs,dwx1,dwy1)
         else if (deltrng.eq.delt) then
            call readwind(fcst_dir,fn,iw,nrecs,dwx2,dwy2)
         endif
c
         do m=1,mmax
            do n=1,nmax

               if (deltrng.eq.(delt/2.d0)) then
                  dwx1(m,n) = (dwx1(m,n)-wx1(m,n)) / frac
                  dwy1(m,n) = (dwy1(m,n)-wy1(m,n)) / frac
                  wx1(m,n) = wx1(m,n) + dwx1(m,n) * fac
                  wy1(m,n) = wy1(m,n) + dwy1(m,n) * fac
               else if (deltrng.eq.delt) then
                  dwx2(m,n) = (dwx2(m,n)-wx2(m,n)) / frac
                  dwy2(m,n) = (dwy2(m,n)-wy2(m,n)) / frac
                  wx2(m,n) = wx2(m,n) + dwx2(m,n) * fac
                  wy2(m,n) = wy2(m,n) + dwy2(m,n) * fac
               endif

            end do
         end do

      else
c
c       increment the values for each re-computation of the spill
c
         do m=1,mmax
            do n=1,nmax
               if (deltrng.eq.(delt/2.d0)) then
                  wx1(m,n)  = wx1(m,n)  + dwx1(m,n)
                  wy1(m,n)  = wy1(m,n)  + dwy1(m,n)
               else if (deltrng.eq.delt) then
                  wx2(m,n)  = wx2(m,n)  + dwx2(m,n)
                  wy2(m,n)  = wy2(m,n)  + dwy2(m,n)
               endif
            end do
         end do

      end if

      if (deltrng.eq.(delt/2.d0)) then
         isubwd1 = isub
         irngwd1 = ifile
         readdata1 = readdata
         iwindrec1 = iwindrec
      else if (deltrng.eq.delt) then
         isubwd2 = isub
         irngwd2 = ifile
         readdata2 = readdata
         iwindrec2 = iwindrec
      endif
 
      return

      END SUBROUTINE wind_runge

c----------------------------------------------------------------------
      subroutine readwind(dirr,fn,iwindrec,nrecs,wfx,wfy)
c     read wind forecast data and interpolate to Medslik grid
c     optional interpolation also over land points
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370)
c
      dimension wfx(mm,nm),wfy(mm,nm),wfx2(mm,nm),wfy2(mm,nm),
     &          itype(mm,nm),wx1(24),wy1(24)
      character dirr*13,fn*14,dummy*80,acall*2

      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk2/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /temp/ m0,n0

      data icall /0/
      icall = icall + 1
      write(acall,'(i2.2)') icall

      open(72,file=dirr//fn)
      read(72,*) dummy
      read(72,*) dummy
      read(72,*) alon1,alon2,alat1,alat2,mmaxs,nmaxs
      read(72,*) ndat
      read(72,*) dummy

      dlon=(alon2-alon1)/dfloat(mmaxs-1)
      dlat=(alat2-alat1)/dfloat(nmaxs-1)

      do i=1,ndat

         read(72,'(10(F12.8))') alat,alon,(wx1(k),wy1(k),k=1,nrecs)

         m=int((alon-alon1)/dlon+1.01d0)
         n=int((alat-alat1)/dlat+1.01d0)
  
         wfx(m,n)=wx1(iwindrec)
         wfy(m,n)=wy1(iwindrec)
  
      end do

      close(72)

c     now interpolate data from wind grid onto medslik grid
c     comment out the line (*) to include interpolation to land points

      nwp = 0
      do 10 m=1,mm
         do 10 n=1,nm
            wfx2(m,n)=0.d0
            wfy2(m,n)=0.d0
c        if(itype(m,n).eq.0) go to 10          ! (*) 
            nwp = nwp + 1

            x = along1 + dfloat(m-1) * dlong
            y = alatg1 + dfloat(n-1) * dlatg
            xdata = (x - alon1) / dlon + 1.d0
            ydata = (y - alat1) / dlat + 1.d0
            mdata=int(xdata)
            ndata=int(ydata)
            if(mdata.lt.1.or.ndata.lt.1.or.mdata.gt.mmaxs.or.
     &         ndata.gt.nmaxs) go to 10

            wfx2(m,n)=(wfx(mdata,ndata)*(dfloat(mdata+1)-xdata)+
     &                 wfx(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &                 *(dfloat(ndata+1)-ydata)+(wfx(mdata,ndata+1)
     &                 *(dfloat(mdata+1)-xdata)+wfx(mdata+1,ndata+1)
     &                 *(xdata-dfloat(mdata)))*(ydata-dfloat(ndata))

            wfy2(m,n)=(wfy(mdata,ndata)*(dfloat(mdata+1)-xdata)+
     &                 wfy(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &                 *(dfloat(ndata+1)-ydata)+(wfy(mdata,ndata+1)
     &                 *(dfloat(mdata+1)-xdata)+wfy(mdata+1,ndata+1)
     &                 *(xdata-dfloat(mdata)))*(ydata-dfloat(ndata))

   10 continue
        
      do m=1,mm
         do n=1,nm
            wfx(m,n)=wfx2(m,n)
            wfy(m,n)=wfy2(m,n)
         end do
      end do

      return

      end subroutine readwind

c**********************************************************************
c     WATER CURRENT routines
c**********************************************************************

c----------------------------------------------------------------------
      subroutine season(jd,seas1,seas2,ss)
c     computes season of the year for a given julian date (<=365)
c----------------------------------------------------------------------

        implicit real*8(a-h,o-z)
        dimension ss(2)
        character seas1*1,seas2*1

        if(jd.ge.46.and.jd.le.136) then
          nseas1=1
          nseas2=2
          ss(1)=dfloat(136-jd)/90.d0
          ss(2)=dfloat(jd-46)/90.d0
        else if(jd.ge.137.and.jd.le.227) then
          nseas1=2
          nseas2=3
          ss(1)=dfloat(227-jd)/90.d0
          ss(2)=dfloat(jd-137)/90.d0        
        else if(jd.ge.228.and.jd.le.318) then
          nseas1=3
          nseas2=4
          ss(1)=dfloat(318-jd)/90.d0
          ss(2)=dfloat(jd-228)/90.d0
        else if(jd.le.45.or.jd.ge.319) then
          nseas1=4
          nseas2=1
          if(jd.ge.319) jd=jd-365
          ss(1)=dfloat(45-jd)/91.d0
          ss(2)=dfloat(jd+46)/91.d0
        end if
      
        write(seas1,'(i1)') nseas1
        write(seas2,'(i1)') nseas2

        return

      end subroutine season

c--------------------------------------------------------------------------------
      SUBROUTINE fcstcur_intrpl(nst,time,delt,nfcst,fcstfn,ifcstfn,
     &                            fcsttim,icurrents,iregn,fcstcurdir,
     &                            ktmx_cu,type_cu)
c
c     This subroutine constructs current field for forecast 
c     data, both instantaneuos (type_cu = hi) than hourly mean (type_cu = hm). 
c     time = current time in hours after 0 hrs on date of spill
c
c---------------------------------------------------------------------------------

        implicit real*8(a-h,o-z)
        parameter(mm=370,nm=370,
     &            ntm=2000,npc=100000,nss=200000,msp=1200)
        save ifile,isub
c
        dimension ifcstfn(720), fcsttim(720), itype(mm,nm)
        character fcstcurdir*13,fcstfn(720)*16, fn*16, a2*2
        character type_cu*3

        common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
        common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
        common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

        common /curr/ sst(mm,nm),us(mm,nm),vs(mm,nm),u10(mm,nm),
     &                v10(mm,nm),u30(mm,nm),v30(mm,nm),u120(mm,nm),
     &                v120(mm,nm),dsst(mm,nm),dus(mm,nm),dvs(mm,nm),
     &                du10(mm,nm),dv10(mm,nm),du30(mm,nm),dv30(mm,nm),
     &                du120(mm,nm),dv120(mm,nm)

        common /temp/ m0,n0

        data ifile /1/,   isub /1/

        if(isub.eq.0) return         ! end of available data files

        m0=x0+0.49
        n0=y0+0.49

        nrecs   = ktmx_cu                    ! no of current files for each day
        dtfcst = 24.d0 / dfloat(nrecs)
        timew  = time
    
        frac = dtfcst / delt       

        if (nst.eq.1) then
c
c     on 1st step set the times of the forecast files (hrs after 0 hrs 
c     on spill date) then open the initial forecast file and read initial data
c
           write(90,*) ' '
           write(90,*) 'Frcst files and file times from 0 hrs 
     &                  on spill day'

          nfcst=nfcst*nrecs

          do i=1,nfcst
             a2=fcstfn(i)(5:6)
             read(a2,'(i2)') iy
             a2=fcstfn(i)(7:8)
             read(a2,'(i2)') im
             a2=fcstfn(i)(9:10)
             read(a2,'(i2)') id
             a2=fcstfn(i)(11:12)
             read(a2,'(i2)') ih
  
             nday = jdiff(idd,imm,iyr,id,im,iy+2000)

             if (type_cu.eq.' hi') then
                fcsttim(i)=nday*24.d0 + dfloat(ih) !snapshot 
             endif
             if (type_cu.eq.' hm') then
                fcsttim(i)=nday*24.d0 + dfloat(ih)-0.5 !average
             endif
             if (type_cu.eq.' dm') then
                fcsttim(i)=nday*24.d0 + dfloat(ih) ! dayly average
             endif
             write(90,*) i,'   ',fcstfn(i),fcsttim(i)

          enddo
        
          do i=1,nfcst
             if(timew.ge.fcsttim(i)) ifile=i ! first current data file which can be used in the simulation
          end do

        end if
c
c     open the next forecast file and read new data (on 1st step file #2)
c
        if (timew.ge.fcsttim(ifile)) then

           fac=(timew - fcsttim(ifile)) / delt
     
           fn=fcstfn(ifile)

           write(6,*) ''
           write(6,*) 'Reading forecast currents from file ',fn

           write(90,*) ' '
           write(90,*) 'Forecast current directory = ',fcstcurdir
           write(90,*) 'Reading forecast currents from file ',fn
           write(90,*) ' '

          call readfcst(fcstcurdir,fn,sst,us,vs,u10,v10,
     &                     u30,v30,u120,v120)

           ifile=ifile+1

           if (ifcstfn(ifile).eq.0) then      ! past the last available file
               isub=0
                write(6,*) 'WARNING: Not enough forecast 
     &                      data is available.'
                write(6,*) 'Forecast data will be kept constant 
     &                      from now on'
                go to 22
           end if

           fn=fcstfn(ifile)

           write(6,*) 'Reading forecast currents from file ',fn

           write(90,*) 'Reading forecast currents from file ',fn
           write(90,*) ' '

           call readfcst(fcstcurdir,fn,dsst,dus,dvs,du10,dv10,
     &                        du30,dv30,du120,dv120)
c
c          calculate the increments for each spill re-computation
c
           do m=1,mmax
              do n=1,nmax

                     dsst(m,n)  = (dsst(m,n)  - sst(m,n))  / frac
                     dus(m,n)   = (dus(m,n)   - us(m,n))   / frac
                     dvs(m,n)   = (dvs(m,n)   - vs(m,n))   / frac
                     du10(m,n)  = (du10(m,n)  - u10(m,n))  / frac
                     dv10(m,n)  = (dv10(m,n)  - v10(m,n))  / frac
                     du30(m,n)  = (du30(m,n)  - u30(m,n))  / frac
                     dv30(m,n)  = (dv30(m,n)  - v30(m,n))  / frac
                     du120(m,n) = (du120(m,n) - u120(m,n)) / frac
                     dv120(m,n) = (dv120(m,n) - v120(m,n)) / frac
 
                     sst(m,n)  = sst(m,n)  + dsst(m,n)  * fac
                     us(m,n)   = us(m,n)   + dus(m,n)   * fac
                     vs(m,n)   = vs(m,n)   + dvs(m,n)   * fac
                     u10(m,n)  = u10(m,n)  + du10(m,n)  * fac
                     v10(m,n)  = v10(m,n)  + dv10(m,n)  * fac
                     u30(m,n)  = u30(m,n)  + du30(m,n)  * fac
                     v30(m,n)  = v30(m,n)  + dv30(m,n)  * fac
                     u120(m,n) = u120(m,n) + du120(m,n) * fac
                     v120(m,n) = v120(m,n) + dv120(m,n) * fac

              end do
           end do

        else
c
c          increment the values for each re-computation of the spill except
c          first

           do m=1,mmax
              do n=1,nmax

                 sst(m,n)  = sst(m,n)  + dsst(m,n)
                 us(m,n)   = us(m,n)   + dus(m,n)
                 vs(m,n)   = vs(m,n)   + dvs(m,n)
                 u10(m,n)  = u10(m,n)  + du10(m,n)
                 v10(m,n)  = v10(m,n)  + dv10(m,n)
                 u30(m,n)  = u30(m,n)  + du30(m,n)
                 v30(m,n)  = v30(m,n)  + dv30(m,n)
                 u120(m,n) = u120(m,n) + du120(m,n)
                 v120(m,n) = v120(m,n) + dv120(m,n)

               end do
            end do

        endif

   22   continue

        return


        END SUBROUTINE fcstcur_intrpl


c----------------------------------------------------------------------
      subroutine readfcst(fcstcurdir,fn,sst,us,vs,u10,v10,
     &                    u30,v30,u120,v120)

c     read forecast current data and interpolate to Medslik grid
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370)
c
      dimension sst(mm,nm),
     &          us(mm,nm),vs(mm,nm),u10(mm,nm),v10(mm,nm),
     &          u30(mm,nm),v30(mm,nm),u120(mm,nm),v120(mm,nm)
      character fcstcurdir*13,fn*16,dummy*80
      common /temp/ m0,n0

      open(71,file=fcstcurdir//fn)

      read(71,*) dummy
      read(71,*) dummy
      read(71,*) alon1,alon2,alat1,alat2,mmax,nmax
      read(71,*) ndata
      read(71,*) dummy

      dlon=(alon2-alon1)/dfloat(mmax-1)
      dlat=(alat2-alat1)/dfloat(nmax-1)

      do m=1,mm
         do n=1,nm
            sst(m,n)=0.d0
            us(m,n)=0.d0
            vs(m,n)=0.d0
            u10(m,n)=0.d0
            v10(m,n)=0.d0
            u30(m,n)=0.d0
            v30(m,n)=0.d0
            u120(m,n)=0.d0
            v120(m,n)=0.d0
         enddo
      enddo

      do i=1,ndata
         read(71,*) alat,alon,st,u,v,u1,v1,u3,v3,u4,v4
         m=int((alon-alon1)/dlon+1.1d0)
         n=int((alat-alat1)/dlat+1.1d0)
         if(m.gt.mm.or.n.gt.nm.or.m.lt.1.or.n.lt.1) go to 1
            sst(m,n)=st
            us(m,n)=u
            vs(m,n)=v
            u10(m,n)=u1
            v10(m,n)=v1
            u30(m,n)=u3
            v30(m,n)=v3
            u120(m,n)=u4
            v120(m,n)=v4
    1       continue
      end do

      close(71)
c
c        now interpolate data onto medslik grid
c       
      call interpol(sst,alon1,alat1,dlon,dlat)
      call interpol(us,alon1,alat1,dlon,dlat)
      call interpol(vs,alon1,alat1,dlon,dlat)
      call interpol(u10,alon1,alat1,dlon,dlat)
      call interpol(v10,alon1,alat1,dlon,dlat)
      call interpol(u30,alon1,alat1,dlon,dlat)
      call interpol(v30,alon1,alat1,dlon,dlat)
      call interpol(u120,alon1,alat1,dlon,dlat)
      call interpol(v120,alon1,alat1,dlon,dlat)
c

      return
      end subroutine readfcst



c--------------------------------------------------------------------------------
      SUBROUTINE curr_runge(nst,time,deltrng,delt,nfile,filename,
     &                      ifilename,fcsttim,ifield,iregn,nrecs,
     &                      timetype)
c
c     This subroutine constructs current field for forecast data for the 
c     Runge scheme, both instantaneuos (type_cu = hi) than hourly mean (type_cu = hm). 
c     time = current time in hours after 0 hrs on date of spill
c---------------------------------------------------------------------------------

        implicit real*8(a-h,o-z)
        parameter(mm=370,nm=370,
     &            ntm=2000,npc=100000,nss=200000,msp=1200)
c
        dimension ifilename(720), fcsttim(720), itype(mm,nm)
        character fcst_dir*13,filename(720)*16, fn*16, a2*2
        character timetype*3
        dimension rsst(mm,nm),rus(mm,nm),rvs(mm,nm),ru10(mm,nm),
     &            rv10(mm,nm),ru30(mm,nm),rv30(mm,nm),ru120(mm,nm),
     &            rv120(mm,nm),rdsst(mm,nm),rdus(mm,nm),rdvs(mm,nm),
     &            rdu10(mm,nm),rdv10(mm,nm),rdu30(mm,nm),rdv30(mm,nm),
     &            rdu120(mm,nm),rdv120(mm,nm)

        common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
        common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
        common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

        common /rng1/ ucur1(mm,nm,4),vcur1(mm,nm,4),usto1(mm,nm),
     &              vsto1(mm,nm),ducr1(mm,nm,4),dvcr1(mm,nm,4),
     &              dust1(mm,nm),dvst1(mm,nm),wx1(mm,nm),wy1(mm,nm),
     &              dwx1(mm,nm),dwy1(mm,nm),irngc1,irngw1,irngwd1,
     &              isubc1,isubw1,isubwd1
        common /rng2/ ucur2(mm,nm,4),vcur2(mm,nm,4),usto2(mm,nm),
     &              vsto2(mm,nm),ducr2(mm,nm,4),dvcr2(mm,nm,4),
     &              dust2(mm,nm),dvst2(mm,nm),wx2(mm,nm),wy2(mm,nm),
     &              dwx2(mm,nm),dwy2(mm,nm),irngc2,irngw2,irngwd2,
     &              isubc2,isubw2,isubwd2

        common /temp/ m0,n0

        if ((ifield.eq.76).or.(ifield.eq.77)) then

           if (nst.ne.1) then
           
              if (deltrng.eq.(delt/2.d0)) then
                 isub = isubc1
                 ifile = irngc1
              endif
              if (deltrng.eq.(delt)) then
                 isub = isubc2
                 ifile = irngc2
              endif
           else

              isub = 1

           endif

           if (isub.eq.0) return       ! end of available data files

           dtfcst = 24.d0 / dfloat(nrecs)
           frac = dtfcst / delt
           timew = time + deltrng 

           write(6,*) '   Linear time interpolation of ocean'
           write(6,*) 'currents for 4th order Runge-Kutta scheme'
           fcst_dir='INP_DATA/OCE/'

           if (nst.eq.1) then
c
c     on 1st step set the times of the forecast files (hrs after 0 hrs 
c     on spill date) then open the initial forecast file and read initial data

              do i=1,nfile
                 a2=filename(i)(5:6)
                 read(a2,'(i2)') iy
                 a2=filename(i)(7:8)
                 read(a2,'(i2)') im
                 a2=filename(i)(9:10)
                 read(a2,'(i2)') id
                 a2=filename(i)(11:12)
                 read(a2,'(i2)') ih

                 nday = jdiff(idd,imm,iyr,id,im,iy+2000)

                 if (timetype.eq.' hi') then
                    fcsttim(i)=nday*24.d0 + dfloat(ih) !snapshot 
                 endif
                 if (timetype.eq.' hm') then
                    fcsttim(i)=nday*24.d0 + dfloat(ih)-0.5 !average
                 endif            
              enddo

              do i=1,nfile
                 if (timew.ge.fcsttim(i)) ifile=i ! first current data file which can be used in the simulation
              end do

           endif
c
c     open the next forecast file and read new data (on 1st step file #2)
c
           if (timew.ge.fcsttim(ifile)) then

              fac=(timew - fcsttim(ifile)) / delt
              fn=filename(ifile)

              write(6,*) ''
              write(6,*) 'Reading forecast field from file ',fn

              call readfcst(fcst_dir,fn,rsst,rus,rvs,ru10,rv10,
     &                          ru30,rv30,ru120,rv120)
              ifile=ifile+1

              if (ifilename(ifile).eq.0) then      ! past the last available file
                 isub=0
                 write(6,*) 'WARNING: Not enough forecast 
     &                     data is available.'
                 write(6,*) 'Forecast data will be kept constant 
     &                     from now on'
                 go to 22
              end if

              fn=filename(ifile)

              write(6,*) ''
              write(6,*) 'Reading forecast field from file ',fn

              call readfcst(fcst_dir,fn,rdsst,rdus,rdvs,rdu10,rdv10,
     &                          rdu30,rdv30,rdu120,rdv120)

              do m=1,mmax
                 do n=1,nmax

                     rdus(m,n)   = (rdus(m,n)   - rus(m,n))   / frac
                     rdvs(m,n)   = (rdvs(m,n)   - rvs(m,n))   / frac
                     rdu10(m,n)  = (rdu10(m,n)  - ru10(m,n))  / frac
                     rdv10(m,n)  = (rdv10(m,n)  - rv10(m,n))  / frac
                     rdu30(m,n)  = (rdu30(m,n)  - ru30(m,n))  / frac
                     rdv30(m,n)  = (rdv30(m,n)  - rv30(m,n))  / frac
                     rdu120(m,n) = (rdu120(m,n) - ru120(m,n)) / frac
                     rdv120(m,n) = (rdv120(m,n) - rv120(m,n)) / frac

                     rus(m,n)   = rus(m,n)   + rdus(m,n)   * fac
                     rvs(m,n)   = rvs(m,n)   + rdvs(m,n)   * fac
                     ru10(m,n)  = ru10(m,n)  + rdu10(m,n)  * fac
                     rv10(m,n)  = rv10(m,n)  + rdv10(m,n)  * fac
                     ru30(m,n)  = ru30(m,n)  + rdu30(m,n)  * fac
                     rv30(m,n)  = rv30(m,n)  + rdv30(m,n)  * fac
                     ru120(m,n) = ru120(m,n) + rdu120(m,n) * fac
                     rv120(m,n) = rv120(m,n) + rdv120(m,n) * fac
                   
                     if (deltrng.eq.(delt/2.d0)) then

                         ucur1(m,n,1) = rus(m,n)
                         vcur1(m,n,1) = rvs(m,n)
                         ucur1(m,n,2) = ru10(m,n)
                         vcur1(m,n,2) = rv10(m,n)
                         ucur1(m,n,3) = ru30(m,n)
                         vcur1(m,n,3) = rv30(m,n)
                         ucur1(m,n,4) = ru120(m,n)
                         vcur1(m,n,4) = rv120(m,n)

                         ducr1(m,n,1) = rdus(m,n)
                         dvcr1(m,n,1) = rdvs(m,n)
                         ducr1(m,n,2) = rdu10(m,n)
                         dvcr1(m,n,2) = rdv10(m,n)
                         ducr1(m,n,3) = rdu30(m,n)
                         dvcr1(m,n,3) = rdv30(m,n)
                         ducr1(m,n,4) = rdu120(m,n)
                         dvcr1(m,n,4) = rdv120(m,n)

                     else if (deltrng.eq.delt) then
                          
                         ucur2(m,n,1) = rus(m,n)
                         vcur2(m,n,1) = rvs(m,n)
                         ucur2(m,n,2) = ru10(m,n)
                         vcur2(m,n,2) = rv10(m,n)
                         ucur2(m,n,3) = ru30(m,n)
                         vcur2(m,n,3) = rv30(m,n)
                         ucur2(m,n,4) = ru120(m,n)
                         vcur2(m,n,4) = rv120(m,n)

                         ducr2(m,n,1) = rdus(m,n)
                         dvcr2(m,n,1) = rdvs(m,n)
                         ducr2(m,n,2) = rdu10(m,n)
                         dvcr2(m,n,2) = rdv10(m,n)
                         ducr2(m,n,3) = rdu30(m,n)
                         dvcr2(m,n,3) = rdv30(m,n)
                         ducr2(m,n,4) = rdu120(m,n)
                         dvcr2(m,n,4) = rdv120(m,n) 

                     endif

                 end do
              end do
           else
              do m=1,mmax
                 do n=1,nmax
                    
                    if (deltrng.eq.(delt/2.d0)) then
                                  
                       ucur1(m,n,1) = ucur1(m,n,1) + ducr1(m,n,1)
                       vcur1(m,n,1) = vcur1(m,n,1) + dvcr1(m,n,1)
                       ucur1(m,n,2) = ucur1(m,n,2) + ducr1(m,n,2)
                       vcur1(m,n,2) = vcur1(m,n,2) + dvcr1(m,n,2)
                       ucur1(m,n,3) = ucur1(m,n,3) + ducr1(m,n,3)
                       vcur1(m,n,3) = vcur1(m,n,3) + dvcr1(m,n,3)
                       ucur1(m,n,4) = ucur1(m,n,4) + ducr1(m,n,4)
                       vcur1(m,n,4) = vcur1(m,n,4) + dvcr1(m,n,4)

                    else if (deltrng.eq.delt) then
                   
                       ucur2(m,n,1) = ucur2(m,n,1) + ducr2(m,n,1)
                       vcur2(m,n,1) = vcur2(m,n,1) + dvcr2(m,n,1)
                       ucur2(m,n,2) = ucur2(m,n,2) + ducr2(m,n,2)
                       vcur2(m,n,2) = vcur2(m,n,2) + dvcr2(m,n,2)
                       ucur2(m,n,3) = ucur2(m,n,3) + ducr2(m,n,3)
                       vcur2(m,n,3) = vcur2(m,n,3) + dvcr2(m,n,3)
                       ucur2(m,n,4) = ucur2(m,n,4) + ducr2(m,n,4)
                       vcur2(m,n,4) = vcur2(m,n,4) + dvcr2(m,n,4)

                    endif

                 end do
              end do

           endif
        endif

        if (deltrng.eq.(delt/2.d0)) then
           isubc1 = isub
           irngc1 = ifile
        endif
        if (deltrng.eq.(delt)) then
           isubc2 = isub
           irngc2 = ifile
        endif
                
   22   continue
        return

      END SUBROUTINE curr_runge

c**********************************************************************
c     WAVE subroutines
c**********************************************************************

c----------------------------------------------------------------------
      SUBROUTINE fcstwav_intrpl(nst,time,delt,nfile,filename,
     &                            ifilename,wavfcsttim,iwave,
     &                            ktmx_wv,type_wv)
c
c      This subroutine constructs wave field for forecast data  
c      both instantaneuos (type_cu = hi) than hourly mean
c      (type_cu = hm). 
c      time = current time in hours after 0 hrs on date of spill
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=200000,msp=1200)
      save ifile,isub

      dimension ifilename(720), wavfcsttim(720), itype(mm,nm)
      character fcstwavdir*13, filename(720)*18, fn*18, a2*2
      character type_wv*3

      common /wave/  stoku(mm,nm),stokv(mm,nm),dstoku(mm,nm),
     &               dstokv(mm,nm)
      common /spill/ idd,imm,iyr,ispill,tstart,tcomp,x0,y0
      common /blk1/  mmax,nmax,delx,dely,itype,pi,degrad
      common /blk3/  mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

      common /temp/ m0,n0

      data ifile /1/,   isub /1/



      if(isub.eq.0) return         ! end of available data files

      m0=x0+0.49
      n0=y0+0.49

      nrecs   = ktmx_wv                    ! no of current files for each day
      dtfcst = 24.d0 / dfloat(nrecs)
      timew  = time 
      frac = dtfcst / delt    

      fcstwavdir='INP_DATA/WAV/'

      if (nst.eq.1) then
c
c     on 1st step set the times of the forecast files (hrs after 0 hrs 
c     on spill date) then open the initial forecast file and read
c     initial data
c
         write(90,*) ' '
         write(90,*) 'Frcst files and file times from 
     &                0 hrs on spill day'

         nfile=nfile*ktmx_wv

         do i=1,nfile

            a2=filename(i)(5:6)
            read(a2,'(i2)') iy
            a2=filename(i)(7:8)
            read(a2,'(i2)') im
            a2=filename(i)(9:10)
            read(a2,'(i2)') id
            a2=filename(i)(11:12)
            read(a2,'(i2)') ih

            nday = jdiff(idd,imm,iyr,id,im,iy+2000)

            if (type_wv.eq.' hi') then
                wavfcsttim(i)=nday*24.d0 + dfloat(ih) !snapshot 
            endif
            if (type_wv.eq.' hm') then
                wavfcsttim(i)=nday*24.d0 + dfloat(ih)-0.5 !average
            endif
            write(90,*) i,'   ',filename(i),wavfcsttim(i)

         enddo

         do i=1,nfile
            if(timew.ge.wavfcsttim(i)) ifile=i ! first current data file which can be used in the simulation
         end do

      endif
c
c     open the next forecast file and read new data (on 1st step file#2)
c     time has already been advanced to the end of the current step
c
      if (timew.ge.wavfcsttim(ifile)) then

         fac=(timew - wavfcsttim(ifile)) / delt

         fn=filename(ifile)

         write(6,*) ''
         write(6,*) 'Reading wave forecast from file ',fn

         write(90,*) ' '
         write(90,*) 'Forecast wave directory = ',fcstwavdir
         write(90,*) 'Reading wave forecast from file ',fn
         write(90,*) ' '

         call readstoke(fcstwavdir,fn,stoku,stokv)

         ifile=ifile+1

         if (ifilename(ifile).eq.0) then      ! past the last available file
             isub=0
             write(6,*) 'WARNING: Not enough forecast data 
     &                   is available.'
             write(6,*) 'Forecast data will be kept constant 
     &                   from now on'
             go to 22
         end if

         fn=filename(ifile)

         write(6,*) 'Reading wave forecast from file ',fn

         write(90,*) 'Reading wave forecast from file ',fn
         write(90,*) ' '

         call readstoke(fcstwavdir,fn,dstoku,dstokv)
c
c        calculate the increments for each spill re-computation
c

         do m=1,mmax
            do n=1,nmax

                  dstoku(m,n) = (dstoku(m,n) - stoku(m,n)) / frac
                  dstokv(m,n) = (dstokv(m,n) - stokv(m,n)) / frac

                  stoku(m,n)  = stoku(m,n) + dstoku(m,n) * fac
                  stokv(m,n)  = stokv(m,n) + dstokv(m,n) * fac

            end do
         end do

      else
c
c        increment the values for each re-computation of the spill except 
c        first
c
         do m=1,mmax
            do n=1,nmax

               stoku(m,n)  = stoku(m,n)  + dstoku(m,n)
               stokv(m,n)  = stokv(m,n)  + dstokv(m,n)

            end do
         end do

      end if

   22 continue

      return


      END SUBROUTINE fcstwav_intrpl

c----------------------------------------------------------------------
      subroutine readstoke(wavdir,fn,ustoke,vstoke)

c     Read forecast stokes drift data and interpolate 
c     to Medslik grid.
c     
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370)
c
      dimension ustoke(mm,nm),vstoke(mm,nm)
      character wavdir*13,fn*18,dummy*80
      common /temp/ m0,n0


      open(71,file=wavdir//fn)

      read(71,*) dummy
      read(71,*) dummy
      read(71,*) alon1,alon2,alat1,alat2,mmax,nmax
      read(71,*) ndata
      read(71,*) dummy


      dlon=(alon2-alon1)/dfloat(mmax-1)
      dlat=(alat2-alat1)/dfloat(nmax-1)

      do m=1,mm

         do n=1,nm

            ustoke(m,n)=0.d0
            vstoke(m,n)=0.d0

         enddo

      enddo

      do i=1,ndata

          read(71,*) alat,alon,u,v

          m=int((alon-alon1)/dlon+1.1d0)
          n=int((alat-alat1)/dlat+1.1d0)

          if(m.gt.mm.or.n.gt.nm.or.m.lt.1.or.n.lt.1) go to 1

          ustoke(m,n) = u
          vstoke(m,n) = v

    1     continue

      end do

      close(71)
c
c   now interpolate data onto medslik grid
c       
      call interpol(ustoke,alon1,alat1,dlon,dlat)
      call interpol(vstoke,alon1,alat1,dlon,dlat)

      return

      end subroutine readstoke
 
c---------------------------------------------------------------------------------
      SUBROUTINE wave_runge(nst,time,deltrng,delt,nfile,filename,
     &                      ifilename,fcsttim,ifield,iregn,nrecs,
     &                      timetype)

c
c     This subroutine constructs wave field for forecast data for the
c     Runge scheme, both instantaneuos (type_cu = hi) than hourly mean 
c     (type_cu = hm). 
c     time = current time in hours after 0 hrs on date of spill
c----------------------------------------------------------------------------------

        implicit real*8(a-h,o-z)
        parameter(mm=370,nm=370,
     &            ntm=2000,npc=100000,nss=200000,msp=1200)
c
        dimension ifilename(720), fcsttim(720), itype(mm,nm)
        character fcst_dir*13,filename(720)*18, fn*18, a2*2
        character timetype*3
        dimension rsst(mm,nm),rus(mm,nm),rvs(mm,nm),ru10(mm,nm),
     &            rv10(mm,nm),ru30(mm,nm),rv30(mm,nm),ru120(mm,nm),
     &            rv120(mm,nm),rdsst(mm,nm),rdus(mm,nm),rdvs(mm,nm),
     &            rdu10(mm,nm),rdv10(mm,nm),rdu30(mm,nm),rdv30(mm,nm),
     &            rdu120(mm,nm),rdv120(mm,nm)

        common /spill/  idd,imm,iyr,ispill,tstart,tcomp,x0,y0
        common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
        common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

        common /rng1/ ucur1(mm,nm,4),vcur1(mm,nm,4),usto1(mm,nm),
     &              vsto1(mm,nm),ducr1(mm,nm,4),dvcr1(mm,nm,4),
     &              dust1(mm,nm),dvst1(mm,nm),wx1(mm,nm),wy1(mm,nm),
     &              dwx1(mm,nm),dwy1(mm,nm),irngc1,irngw1,irngwd1,
     &              isubc1,isubw1,isubwd1

        common /rng2/ ucur2(mm,nm,4),vcur2(mm,nm,4),usto2(mm,nm),
     &              vsto2(mm,nm),ducr2(mm,nm,4),dvcr2(mm,nm,4),
     &              dust2(mm,nm),dvst2(mm,nm),wx2(mm,nm),wy2(mm,nm),
     &              dwx2(mm,nm),dwy2(mm,nm),irngc2,irngw2,irngwd2,
     &              isubc2,isubw2,isubwd2

        common /temp/ m0,n0

        if (ifield.eq.101) then

           if (nst.ne.1) then

              if (deltrng.eq.(delt/2.d0)) then
                 isub = isubw1
                 ifile = irngw1
              endif
              if (deltrng.eq.(delt)) then
                 isub = isubw2
                 ifile = irngw2
              endif
           else

              isub = 1

           endif

           if (isub.eq.0) return       ! end of available data files

           dtfcst = 24.d0 / dfloat(nrecs)
           frac = dtfcst / delt
           timew = time + deltrng 

           write(6,*) 'Linear time interpolation of Stokes Drift'
           write(6,*) 'velocities for 4th order Runge-Kutta scheme'
           fcst_dir='INP_DATA/WAV/'
        

           if (nst.eq.1) then
c
c     on 1st step set the times of the forecast files (hrs after 0 hrs 
c     on spill date) then open the initial forecast file and read initial data
c
              do i=1,nfile
                 a2=filename(i)(5:6)
                 read(a2,'(i2)') iy
                 a2=filename(i)(7:8)
                 read(a2,'(i2)') im
                 a2=filename(i)(9:10)
                 read(a2,'(i2)') id
                 a2=filename(i)(11:12)
                 read(a2,'(i2)') ih

                 nday = jdiff(idd,imm,iyr,id,im,iy+2000)

                 if (timetype.eq.' hi') then
                    fcsttim(i)=nday*24.d0 + dfloat(ih) !snapshot 
                 endif

                 if (timetype.eq.' hm') then
                    fcsttim(i)=nday*24.d0 + dfloat(ih)-0.5 !average
                 endif          
              enddo

              do i=1,nfile
                 if (timew.ge.fcsttim(i)) ifile=i ! first current data file which can be used in the simulation
              end do

           endif

c
c     open the next forecast file and read new data (on 1st step file #2)
c     time has already been advanced to the end of the current step
c
           if (timew.ge.fcsttim(ifile)) then


              fac=(timew - fcsttim(ifile)) / delt
              fn=filename(ifile)

              write(6,*) ''
              write(6,*) 'Reading forecast field from file ',fn
               
              call readstoke(fcst_dir,fn,rus,rvs) 
          
              ifile=ifile+1

              if (ifilename(ifile).eq.0) then      ! past the last available file
                 isub=0
                 write(6,*) 'WARNING: Not enough forecast 
     &                       data is available.'
                 write(6,*) 'Forecast data will be kept constant 
     &                       from now on'
                 go to 22
              end if

              fn=filename(ifile)

              write(6,*) ''
              write(6,*) 'Reading forecast field from file ',fn

              call readstoke(fcst_dir,fn,rdus,rdvs)
               
              do m=1,mmax
                 do n=1,nmax

                    rdus(m,n) = (rdus(m,n) - rus(m,n)) / frac
                    rdvs(m,n) = (rdvs(m,n) - rvs(m,n)) / frac

                    rus(m,n)  = rus(m,n) + rdus(m,n) * fac
                    rvs(m,n)  = rvs(m,n) + rdvs(m,n) * fac

                    if (deltrng.eq.(delt/2.d0)) then
                        usto1(m,n) = rus(m,n)
                        vsto1(m,n) = rvs(m,n)
                        dust1(m,n) = rdus(m,n)
                        dvst1(m,n) = rdvs(m,n)
                    elseif (deltrng.eq.delt) then
                        usto2(m,n) = rus(m,n)
                        vsto2(m,n) = rvs(m,n)
                        dust2(m,n) = rdus(m,n)
                        dvst2(m,n) = rdvs(m,n)
                    endif
                 end do
              end do          
           else
              do m=1,mmax
                 do n=1,nmax

                    if (deltrng.eq.(delt/2.d0)) then
                       usto1(m,n) = usto1(m,n) + dust1(m,n)
                       vsto1(m,n) = vsto1(m,n) + dvst1(m,n)
                    elseif (deltrng.eq.delt) then
                       usto2(m,n) = usto2(m,n) + dust2(m,n)
                       vsto2(m,n) = vsto2(m,n) + dvst2(m,n)
                    endif
                 end do
              end do
           endif
        endif

        if (deltrng.eq.(delt/2.d0)) then
              isubw1 = isub
              irngw1 = ifile
        endif
        if (deltrng.eq.(delt)) then
              isubw2 = isub
              irngw2 = ifile
        endif
     
   22   continue

        return

        END SUBROUTINE wave_runge

c******************************************************************************
c                      NUMERICAL SCHEME SUBROUTINE
c******************************************************************************
c------------------------------------------------------------------------------
        SUBROUTINE runge(is,px0,py0,psz,hint,fcstdep1,fcstdep2,fcstdep3,
     &                   delt,ui,vi,xdc,ydc)

c     This subroutine calculates the solution af the linear advection equation
c     by using a 4th order Rung-Kutta numerical scheme
c------------------------------------------------------------------------------

        implicit real*8(a-h,o-z)
        parameter(mm=370,nm=370,
     &            ntm=2000,npc=100000,nss=200000,msp=1200)

        dimension itype(mm,nm),xk(4),yk(4)

        common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
        common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2

        common /rng1/ ucur1(mm,nm,4),vcur1(mm,nm,4),usto1(mm,nm),
     &                vsto1(mm,nm),ducr1(mm,nm,4),dvcr1(mm,nm,4),
     &                dust1(mm,nm),dvst1(mm,nm),wx1(mm,nm),wy1(mm,nm),
     &                dwx1(mm,nm),dwy1(mm,nm),irngc1,irngw1,irngwd1,
     &                isubc1,isubw1,isubwd1
        common /rng2/ ucur2(mm,nm,4),vcur2(mm,nm,4),usto2(mm,nm),
     &                vsto2(mm,nm),ducr2(mm,nm,4),dvcr2(mm,nm,4),
     &                dust2(mm,nm),dvst2(mm,nm),wx2(mm,nm),wy2(mm,nm),
     &                dwx2(mm,nm),dwy2(mm,nm),irngc2,irngw2,irngwd2,
     &                isubc2,isubw2,isubwd2
        common /rng3/ urng1(mm,nm),vrng1(mm,nm),urng2(mm,nm),
     &                vrng2(mm,nm)
 
        hsec = delt * 3600.d0

        coef1 = 1.d0 / 6.d0
        coef2 = 1.d0 / 3.d0
        coef3 = 1.d0 / 3.d0
        coef4 = 1.d0 / 6.d0

        xk(1) = hsec*ui
        yk(1) = hsec*vi

        do k=2,4

           if (k.lt.4) then
              a =0.5
           else if (k.eq.4) then
              a = 1.d0
           endif

           psx = px0 + (xk(k-1)*a)/delx
           psy = py0 + (yk(k-1)*a)/dely

           if (is.ne.3) then

              if ((k.eq.2).or.(k.eq.3)) then
                 call intrpl(psx,psy,urng1,itype,uir)
                 call intrpl(psx,psy,vrng1,itype,vir)
              else if (k.eq.4) then
                 call intrpl(psx,psy,urng2,itype,uir)
                 call intrpl(psx,psy,vrng2,itype,vir)
              endif
              zdispl=0.d0

           elseif (is.eq.3) then

              m = int(psx + 0.5d0)
              n = int(psy + 0.5d0)

              dep=(1.d0-psz)*hint     
              zdispl=dep              

           if (dep.le.fcstdep1) then
               if ((k.eq.2).or.(k.eq.3)) then
                  uir=(ucur1(m,n,1)*(fcstdep1-dep)+ucur1(m,n,2)*dep)/10.d0
                  vir=(vcur1(m,n,1)*(fcstdep1-dep)+vcur1(m,n,2)*dep)/10.d0
               else if (k.eq.4) then
                  uir=(ucur2(m,n,1)*(fcstdep1-dep)+ucur2(m,n,2)*dep)/10.d0
                  vir=(vcur2(m,n,1)*(fcstdep1-dep)+vcur2(m,n,2)*dep)/10.d0
               endif
           elseif (dep.gt.fcstdep1.and.dep.le.fcstdep2) then
               if (hint.ge.fcstdep2) then
                  denom=fcstdep2-fcstdep1
                  fac1=(dep-fcstdep1)
                  fac2=(fcstdep2-dep)
                  if ((k.eq.2).or.(k.eq.3)) then
                     uir=(ucur1(m,n,3)*fac1+ucur1(m,n,2)*fac2)/denom
                     vir=(vcur1(m,n,3)*fac1+vcur1(m,n,2)*fac2)/denom
                  else if (k.eq.4) then
                     uir=(ucur2(m,n,3)*fac1+ucur2(m,n,2)*fac2)/denom
                     vir=(vcur2(m,n,3)*fac1+vcur2(m,n,2)*fac2)/denom
                  endif
               else
                  denom=hint-fcstdep1
                  fac1=(dep-fcstdep1)
                  fac2=(hint-dep)
                  if ((k.eq.2).or.(k.eq.3)) then
                     uir=(ucur1(m,n,3)*fac1+ucur1(m,n,2)*fac2)/denom
                     vir=(vcur1(m,n,3)*fac1+vcur1(m,n,2)*fac2)/denom
                  else if (k.eq.4) then
                     uir=(ucur2(m,n,3)*fac1+ucur2(m,n,2)*fac2)/denom
                     vir=(vcur2(m,n,3)*fac1+vcur2(m,n,2)*fac2)/denom
                  endif
               endif
           elseif (dep.gt.fcstdep2.and.dep.le.fcstdep3) then
               if (hint.ge.fcstdep3) then
                   denom=fcstdep3-fcstdep2
                   fac1=(dep-fcstdep2)
                   fac2=(fcstdep3-dep)
                   if ((k.eq.2).or.(k.eq.3)) then
                      uir=(ucur1(m,n,4)*fac1+ucur1(m,n,3)*fac2)/denom
                      vir=(vcur1(m,n,4)*fac1+vcur1(m,n,3)*fac2)/denom
                   else if (k.eq.4) then
                      uir=(ucur2(m,n,4)*fac1+ucur2(m,n,3)*fac2)/denom
                      vir=(vcur2(m,n,4)*fac1+vcur2(m,n,3)*fac2)/denom
                   endif
               elseif(hint.lt.fcstdep2) then
                   if ((k.eq.2).or.(k.eq.3)) then
                      uir=ucur1(m,n,3)
                      vir=vcur1(m,n,3)
                   else if (k.eq.4) then
                      uir=ucur2(m,n,3)
                      vir=vcur2(m,n,3)
                   endif
               else
                   denom=hint-fcstdep2
                   fac1=(dep-fcstdep2)
                   fac2=(hint-dep)
                   if ((k.eq.2).or.(k.eq.3)) then
                      uir=(ucur1(m,n,4)*fac1+ucur1(m,n,3)*fac2)/denom
                      vir=(vcur1(m,n,4)*fac1+vcur1(m,n,3)*fac2)/denom
                   else if (k.eq.4) then
                      uir=(ucur2(m,n,4)*fac1+ucur2(m,n,3)*fac2)/denom
                      vir=(vcur2(m,n,4)*fac1+vcur2(m,n,3)*fac2)/denom
                   endif
               endif
           elseif (dep.gt.fcstdep3) then
               if ((k.eq.2).or.(k.eq.3)) then
                  uir=ucur1(m,n,4)
                  vir=vcur1(m,n,4)
               else if (k.eq.4) then
                  uir=ucur2(m,n,4)
                  vir=vcur2(m,n,4)
               endif
           end if

           end if

           xk(k) = hsec * uir
           yk(k) = hsec * vir

        enddo

        xdc = coef1*xk(1) + coef2*xk(2) + coef3*xk(3) + coef4*xk(4)
        ydc = coef1*yk(1) + coef2*yk(2) + coef3*yk(3) + coef4*yk(4) 

        return
          
        END SUBROUTINE runge
 
c**********************************************************************
c     UTILITY routines
c**********************************************************************

c----------------------------------------------------------------------
      subroutine interpol(array,alon1,alat1,dlon,dlat)
c     interpolate an array to the medslik grid - dble precision
c----------------------------------------------------------------------

      
      implicit real*8 (a-h,o-z)
      parameter(mm=370,nm=370)
c	
      dimension array(mm,nm),array2(mm,nm)
      dimension itype(mm,nm)

      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      common /blk2/ along1,alatg1,along2,alatg2,dlong,dlatg
      common /blk3/ mzb,mzf,nzb,nzf,mzb2,mzf2,nzb2,nzf2
      common /temp/ m0,n0
c
c	first extrapolate to the whole grid (i.e. over the land points)

      knt0=0
    1 continue
      knt=0
          do m=mzb2,mzf2
             do n=nzb2,nzf2
                array2(m,n)=array(m,n)
                if(array2(m,n).ne.0.d0) go to 5

                knt=knt+1

                sum=0.d0
                no=0
                do i=m-1,m+1
                do j=n-1,n+1

                   if(array(i,j).ne.0.d0.and.i.ge.1.and.j.ge.1.and.
     &                i.le.mm.and.j.le.nm) then
                       sum=sum+array(i,j)
                       no=no+1
                   end if
                end do
                end do
                if(no.gt.0) array2(m,n)=sum/dfloat(no)

    5           continue
             end do
          end do
        
          do m=mzb2,mzf2
             do n=nzb2,nzf2
                array(m,n)=array2(m,n)
             end do
          end do

          knt0=knt0+1

          if(knt.gt.0.and.knt0.le.50) go to 1
c
c	now interpolate from the data grid to the medslik grid 
c

          do m=mzb2,mzf2
             do n=nzb2,nzf2
                array2(m,n)=0.d0
c        if(itype(m,n).eq.0) go to 10
                x = along1 + dfloat(m-1) * dlong
                y = alatg1 + dfloat(n-1) * dlatg
                xdata = (x - alon1) / dlon + 1.d0
                ydata = (y - alat1) / dlat + 1.d0
                mdata=int(xdata)
                ndata=int(ydata)
                if(mdata.lt.1.or.ndata.lt.1.or.mdata.gt.mm-1.
     &             or.ndata.gt.nm-1) go to 10

                array2(m,n)=(array(mdata,ndata)*(dfloat(mdata+1)-
     &          xdata)+array(mdata+1,ndata)*(xdata-dfloat(mdata)))
     &          * (dfloat(ndata+1)-ydata)+(array(mdata,ndata+1)*
     &          (dfloat(mdata+1)-xdata)+array(mdata+1,ndata+1)*
     &          (xdata-dfloat(mdata)))* (ydata-dfloat(ndata))

   10           continue
             end do
          end do
        
          do m=mzb2,mzf2
             do n=nzb2,nzf2
                array(m,n) = array2(m,n)
             enddo
          enddo

      return
      end subroutine interpol        

c----------------------------------------------------------------------
      subroutine intrpl0(x,y,q,qint)
c     subroutine intrpl0 interpolates values of an array q
c     from grid values to a point (x,y)
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370)
c
      dimension q(mm,nm),itype(mm,nm)
      m=int(x)
      n=int(y)
      q1=q(m,n)
      q2=q(m+1,n)
      q3=q(m+1,n+1)
      q4=q(m,n+1)
      qint=(q1*(m+1-x)+q2*(x-m))*(n+1-y)+
     &      (q4*(m+1-x)+q3*(x-m))*(y-n)
      return
      end subroutine intrpl0

c----------------------------------------------------------------------
      subroutine intrpl(x,y,q,itype,qint)
c     subroutine intrpl interpolates values of an array q from grid 
c     values to a point (x,y) taking account of land/water mask
c----------------------------------------------------------------------


      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=200000,msp=1200)
c     &          ntm=2000,npc=100000,nss=200000,npl=2000,msp=1200)
c
      dimension q(mm,nm),itype(mm,nm)
      m=int(x)
      n=int(y)
      it1=itype(m,n)
      it2=itype(m+1,n)
      it3=itype(m+1,n+1)
      it4=itype(m,n+1)
      q1=q(m,n)
      q2=q(m+1,n)
      q3=q(m+1,n+1)
      q4=q(m,n+1)


      if(it1.eq.0.and.it2.ne.0.and.it3.ne.0.and.it4.ne.0) then
        q1=q2+q4-q3
      else if(it1.ne.0.and.it2.eq.0.and.it3.ne.0.and.it4.ne.0) then
        q2=q1+q3-q4
      else if(it1.ne.0.and.it2.ne.0.and.it3.eq.0.and.it4.ne.0) then
        q3=q2+q4-q1
      else if(it1.ne.0.and.it2.ne.0.and.it3.ne.0.and.it4.eq.0) then
        q4=q1+q3-q2

      else if(it1.eq.0.and.it2.eq.0.and.it3.ne.0.and.it4.ne.0) then
        if(itype(m,n+2).ne.0) then
          q1=2.*q4-q(m,n+2)
        else
          q1=q4
        end if
        if(itype(m+1,n+2).ne.0) then  
          q2=2.*q3-q(m+1,n+2)
        else
          q2=q3
        end if
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=q1
        end if  
        if(itype(m-1,n+1).ne.0) then

          q3=2.*q4-q(m-1,n+1)
        else
          q3=q4
        end if
      
      else if(it1.ne.0.and.it2.ne.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m+1,n-1).ne.0) then
          q3=2.*q2-q(m+1,n-1)
        else
          q3=q2
        end if
        if(itype(m,n-1).ne.0) then
          q4=2.*q1-q(m,n-1)
        else
          q4=q1
        end if
      
      else if(it1.eq.0.and.it2.ne.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=q3
        end if        
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=q2
        end if        
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=(q1+q3)/2.d0
        end if        
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=(q1+q3)/2.d0
        end if
      
      else if(it1.eq.0.and.it2.ne.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=(q2+q4)/2.d0
        end if
        if(itype(m-1,n+1).ne.0) then
          q3=2.*q4-q(m-1,n+1)
        else
          q3=(q2+q4)/2.d0
        end if
      
      else if(it1.ne.0.and.it2.eq.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m-1,n).ne.0) then
          q2=2.*q1-q(m-1,n)
        else
          q2=q1
        end if
        if(itype(m,n-1).ne.0) then
          q4=2.*q1-q(m,n-1)
        else
          q4=q1
        end if
        q3=q2+q4-q1

      else if(it1.eq.0.and.it2.ne.0.and.it3.eq.0.and.it4.eq.0) then
        if(itype(m+2,n).ne.0) then
          q1=2.*q2-q(m+2,n)
        else
          q1=q2
        end if
        if(itype(m+1,n-1).ne.0) then
          q3=2.*q2-q(m+1,n-1)
        else
          q3=q2
        end if
        q4=q1+q3-q2

      else if(it1.eq.0.and.it2.eq.0.and.it3.ne.0.and.it4.eq.0) then
        if(itype(m+1,n+2).ne.0) then
          q2=2.*q3-q(m+1,n+2)
        else
          q2=q3
        end if
        if(itype(m+2,n+1).ne.0) then
          q4=2.*q3-q(m+2,n+1)
        else
          q4=q3
        end if
        q1=q2+q4-q3

      else if(it1.eq.0.and.it2.eq.0.and.it3.eq.0.and.it4.ne.0) then
        if(itype(m,n+2).ne.0) then
          q1=2.*q4-q(m,n+2)
        else
          q1=q4
        end if
        if(itype(m-1,n+1).ne.0) then
          q3=2.*q4-q(m-1,n+1)
        else
          q3=q4
        end if
        q2=q1+q3-q4
      end if
c
      qint=(q1*(m+1-x)+q2*(x-m))*(n+1-y)+
     &     (q4*(m+1-x)+q3*(x-m))*(y-n)
      return

      end subroutine intrpl

c----------------------------------------------------------------------
      subroutine hsmoot(h)
c     smooth the bathymetry - slmin = max bottom slope
c---------------------------------------------------------------------- 

      implicit real*8(a-h,o-z)
      parameter(mm=370,nm=370,
     &          ntm=2000,npc=100000,nss=200000,msp=1200)
c     &          ntm=2000,npc=100000,nss=200000,npl=2000,msp=1200)
      dimension itype(mm,nm),h(mm,nm) 
      common /blk1/ mmax,nmax,delx,dely,itype,pi,degrad
      
      slmin=0.2d0
      
      do 10 k=1,10
c         sweep right
        do 3 n=1,nmax-1
          do 1 m=1,mmax-1
            if(itype(m,n).eq.0.or.itype(m+1,n).eq.0) go to 1
            sl=dabs(h(m+1,n)-h(m,n))/(h(m+1,n)+h(m,n))
            if(sl.lt.slmin) go to 1
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m+1,n))
            sn=-1.d0
            if(h(m+1,n).gt.h(m,n)) sn=1
            h(m+1,n)=h(m+1,n)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    1     continue
C           sweep left
          do 2 m=mmax-1,1,-1
            if(itype(m,n).eq.0.or.itype(m+1,n).eq.0) go to 2
            sl=dabs(h(m+1,n)-h(m,n))/(h(m+1,n)+h(m,n))
            if(sl.lt.slmin) go to 2
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m+1,n))
            sn=-1.d0
            if(h(m+1,n).gt.h(m,n)) sn=1
            h(m+1,n)=h(m+1,n)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    2     continue
    3   continue
   
c         sweep up
        do 6 m=1,mmax-1
          do 4 n=1,nmax-1
            if(itype(m,n).eq.0.or.itype(m,n+1).eq.0) go to 4
            sl=dabs(h(m,n+1)-h(m,n))/(h(m,n+1)+h(m,n))
            if(sl.lt.slmin) go to 4
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m,n+1))
            sn=-1.d0
            if(h(m,n+1).gt.h(m,n)) sn=1
            h(m,n+1)=h(m,n+1)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    4     continue
C           sweep down
          do 5 n=nmax-1,1,-1
            if(itype(m,n).eq.0.or.itype(m,n+1).eq.0) go to 5
            sl=dabs(h(m,n+1)-h(m,n))/(h(m,n+1)+h(m,n))
            if(sl.lt.slmin) go to 5
            dh=0.5d0*(sl-slmin)*(h(m,n)+h(m,n+1))
            sn=-1.d0
            if(h(m,n+1).gt.h(m,n)) sn=1
            h(m,n+1)=h(m,n+1)-sn*dh
            h(m,n)=h(m,n)+sn*dh
    5     continue
    6   continue
   
   10 continue

      return
      end subroutine hsmoot

c----------------------------------------------------------------------
      function julday(idd,imm,iyr)
c     function julday computes the julian day for a given date
c----------------------------------------------------------------------

      dimension js(12),je(12)
c
      ly=0
      if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or.
     &    ( (iyr/400)*400.eq.iyr) ) ly=1
      js(1)=1
      je(1)=31
      js(2)=32
      je(2)=59+ly
      k=1
      do 5 i=3,12
        js(i)=je(i-1)+1
        je(i)=js(i)+29+k
        k=1-k
        if(i.eq.7) k=1
    5 continue
    
      julday=js(imm)+idd-1
    
      return
      end function julday 

c----------------------------------------------------------------------
      subroutine  date(jd,iyr,id1,im1)
c     subroutine date computes the date for a given julian day no and year
c----------------------------------------------------------------------
      dimension js(12),je(12)
c
    1 continue
      ly=0
      if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or. 
     &      ((iyr/400)*400.eq.iyr) ) ly=1
      js(1)=1
      je(1)=31
      js(2)=32
      je(2)=59+ly
      k=1
      do 5 i=3,12
        js(i)=je(i-1)+1
        je(i)=js(i)+29+k
        k=1-k
        if (i.eq.7) k=1
    5 continue
c
      iy1=iyr
      if(jd.gt.je(12)) then
        iyr=iyr+1
        jd=jd-je(12)
        go to 1
      else if(jd.lt.0) then
        iyr=iyr-1
        jd=jd+je(12)
        if( ( ((iyr/4)*4.eq.iyr) .and. ((iyr/100)*100.ne.iyr) ) .or. 
     &      ((iyr/400)*400.eq.iyr) ) jd=jd+1
        go to 1
      end if
c          
      do 10 i=1,12
        if(jd.ge.js(i).and.jd.le.je(i)) im1=i
   10 continue
      id1=jd-js(im1)+1
c
      return
      end subroutine  date           

c----------------------------------------------------------------------
      function  jdiff(id1,im1,iy1,id2,im2,iy2)
c     function jdiff gives the no of days from id1/im1/iy1 to id2/im2/iy2 
c     the day-difference may be negative and the years differ by 1
c----------------------------------------------------------------------

      if(iy2.eq.iy1) jdiff = julday(id2,im2,iy2) - julday(id1,im1,iy1)
      if(iy2.eq.iy1+1) jdiff = julday(id2,im2,iy2) + julday(31,12,iy1)
     &                                  - julday(id1,im1,iy1)
      if(iy2.eq.iy1-1) jdiff = julday(id2,im2,iy2) - julday(31,12,iy2)
     &                                  - julday(id1,im1,iy1)

      return
      end function  jdiff

c----------------------------------------------------------------------
      subroutine merparam(avlat)
c    ....parameters for mercator metric coords....
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      pi=3.1415926536d0
      degrad1=pi/180.d0
      fi=avlat*degrad1
      ec=0.082094438d0
      ecc=ec/2.d0
      e1=ec*dsin(fi)
      demod=dsqrt(1.d0-e1*e1)/dcos(fi)
      a=6378137.d0

      return
      end subroutine merparam

c----------------------------------------------------------------------
      subroutine ll2mer(avlat,alat,alon,x,y,nst,maxst)
c    ....converts latitude, longitude to mer coords....
c----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      save ind      
c
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      data ind /0/
      if(ind.eq.0) call merparam(avlat)
      ind=1

      rlat=alat*degrad1
      rlon=alon*degrad1
      x=a*rlon/demod
      ec1=ec*dsin(rlat)
      vy=( (1.d0-ec1)/(1.d0+ec1) ) ** ecc
      vy=vy*dtan(pi/4.d0+rlat/2.d0)
      y=(a/demod)*dlog(vy)
      
      return
      end subroutine ll2mer

c----------------------------------------------------------------------
      subroutine mer2ll(avlat,x,y,alat,alon)
c    ....converts mercator metric coords to latitude, longitude....
c----------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      save ind
c
      common /data1/ pi,degrad1,fi,ec,ecc,demod,a
      
      data ind /0/
      if(ind.eq.0) call merparam(avlat)
      ind=1

      rlon=x*demod/a
      rlat=0.68d0
    1 continue
        
        rlat1=rlat
        ec1=ec*dsin(rlat1)
        vy=( (1.d0-ec1)/(1.d0+ec1) ) ** ecc
        ve=dexp(y*demod/a) / vy
        rlat=2.d0*datan(ve) - pi/2.d0
        if(dabs((rlat-rlat1)/rlat).gt.0.0000001) go to 1
        
      alat=rlat/degrad1
      alon=rlon/degrad1
      return
      end subroutine mer2ll

c----------------------------------------------------------------------
      function randmedslik(ix)
c     function rand computes random numbers between 0.0 and 1.0
c----------------------------------------------------------------------
 
      implicit real*8(a-h,o-z)
    
      integer a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k

      data a/16807/,b15/32768/,b16/65536/,p/2147483647/
      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if(ix.lt.0) ix=ix+p
      randmedslik=dfloat(ix)*4.656612875d-10
      return
      end function randmedslik

c----------------------------------------------------------------------
      subroutine seedmedslik(ix)
c     calculates an approximately random first seed for rand
c----------------------------------------------------------------------

        real*8 hms
c        character date*10
c        call date_and_time(time=date)
c        read(date,'(f10.3)') hms
        character date*8, time*10
        call date_and_time(date,time)
        call date_and_time(DATE=date)
        call date_and_time(TIME=time)
        read(time,'(f10.3)') hms
        ix=int(hms*100)
        iz=ix-(ix/1000000)*1000000
        iy=1000000
        ix=0
        do k=1,6
          ix0=iz
          iz=iz/10
          ix1=ix0-iz*10
          iy=iy/10
          ix=ix+ix1*iy
        end do
        write(90,*) 'Random seed = ',ix
        return
        end subroutine seedmedslik
      include 'codes.h'

c----------------------------------------------------------------------
      subroutine readsat(ix,ntot,numspills,nppmsv,px,py)
c     READS POLYGONS (isat=1) CONTOUR
c----------------------------------------------------------------------

         implicit real*8                 (a-h,o-z)
         integer,parameter               :: npc=100000
         integer                         :: err, npol
         real,dimension(10000,2)         :: segx, segy
         real*8,dimension(npc)           :: px, py
         integer,dimension(20)           :: nppmsv
         integer,dimension(numspills+1)  :: pol_indx
         character                       :: datestamp*10, aslik*2
         character                       :: empty*80

         data datestamp /'0711020357'/, aslik /'01'/


         open(38,file='medslik_plgs.inp')
         do i=1,3
            read(38,*) empty
         enddo
         read(38,*) numspills
         do i=1,numspills
            read(38,*) empty
         enddo
         
         read(38,*) ndata
         read(38,*) empty
         read(38,*) yini,xini

         npol = 1
         pol_indx(npol)=1
         nsegs = 1
         read(38,*) ystart, xstart
         segx(nsegs,1) = xstart
         segy(nsegs,1) = ystart
         xlast = xstart
         ylast = ystart

! HAND MADE polygon/s (to be changed for multpile polygons with ndata > 30)

         if ( ndata.lt.40 ) then
    
            do k=3,ndata

               read(38,*,IOSTAT=err) y,x

               if (err.ne.0) then
                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart
                  exit
               endif
   
               if ( dabs(x-xlast).lt.1.and.dabs(y-ylast).lt.1
     &              .and.( xpre.ne.xini.or.ypre.ne.yini ) ) then

                  segx(nsegs,2) = x
                  segy(nsegs,2) = y

                  nsegs = nsegs + 1
                  segx(nsegs,1) = x
                  segy(nsegs,1) = y
                  xlast = x
                  ylast = y
                  xpre=x 
                  ypre=y

               else   

                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart


                  npol = npol + 1
                  pol_indx(npol) = int(nsegs)
                 
                  xini=x
                  yini=y

                  nsegs = nsegs + 1 
                  read(38,*,IOSTAT=err) y,x
c                  if (err.eq.0) then
                  xstart = x
                  ystart = y
                  segx(nsegs,1) = xstart
                  segy(nsegs,1) = ystart

                  xlast = xstart
                  ylast = ystart

               endif

               if ( k.eq.ndata ) then

                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart
                  
               endif  

            enddo 

         else
      
            do k=3,ndata

               read(38,*,IOSTAT=err) y,x

               if (err.ne.0) then
                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart
                  exit
               endif
 
               if ( dabs(x-xlast).lt.0.01.and.dabs(y-ylast).lt.0.01
     &              .and.( xpre.ne.xini.or.ypre.ne.yini ) ) then



                                 
                  segx(nsegs,2) = x
                  segy(nsegs,2) = y
                  nsegs = nsegs + 1
                  segx(nsegs,1) = x
                  segy(nsegs,1) = y

                  xpre = xlast 
                  ypre = ylast
                  xlast = x
                  ylast = y
               
               else   
                
                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart


                  npol = npol + 1
                  pol_indx(npol) = int(nsegs)

                  xini=x
                  yini=y

                  nsegs = nsegs + 1 
                  read(38,*,IOSTAT=err) y,x

                  xstart = x
                  ystart = y
                  segx(nsegs,1) = xstart
                  segy(nsegs,1) = ystart
                 
                  xlast = xstart
                  ylast = ystart

               endif

               if ( k.eq.ndata ) then

                  segx(nsegs,2) = xstart
                  segy(nsegs,2) = ystart
                  
               endif  

            enddo 

         endif
                 
         pol_indx(npol+1) = int(nsegs)

         npcl_count = 0

         do k=1,npol

            if ( k.eq.1) then
                 istart = pol_indx(k)
                 ilast = pol_indx(k+1)
            else
                 istart = pol_indx(k) + 1
                 ilast = pol_indx(k+1)
            endif

            box_xmax = segx(istart,1)
            box_ymax = segy(istart,1)
            box_xmin = segx(istart,1)
            box_ymin = segy(istart,1)
 
            do i=istart,ilast
               do j=1,2
                  if ( box_xmax.lt.segx(i,j) ) box_xmax = segx(i,j)
                  if ( box_ymax.lt.segy(i,j) ) box_ymax = segy(i,j)
                  if ( box_xmin.gt.segx(i,j) ) box_xmin = segx(i,j)
                  if ( box_ymin.gt.segy(i,j) ) box_ymin = segy(i,j)
               enddo  
            enddo
 
            nppmsv(k) = 0

  107       continue

              randx = randmedslik(ix) * (box_xmax-box_xmin) + box_xmin
              randy = randmedslik(ix) * (box_ymax-box_ymin) + box_ymin
              ind = 0

              do i=pol_indx(k),pol_indx(k+1)

                 if ( ( randx.ge.segx(i,1).and.randx.lt.segx(i,2) ).or.
     &              ( randx.le.segx(i,1).and.randx.gt.segx(i,2) ) ) then

                    y_num = (randx - segx(i,1)) * (segy(i,2) - 
     &                      segy(i,1))+ segy(i,1) * (segx(i,2) - 
     &                      segx(i,1))
                    y_dem = segx(i,2) - segx(i,1)

                    if ( y_dem.ne.0.d0 ) then

                        y_int = y_num / y_dem
                        if(y_int.gt.randy) ind = ind + 1

                    endif 

                 endif

              enddo
             
              if ( mod(ind,2).eq.1 ) then

                  nppmsv(k) = nppmsv(k)+1
                  npcl_count=npcl_count+1

                  px(npcl_count) = randx
                  py(npcl_count) = randy
                  if (nppmsv(k).ge.(ntot/npol)) go to 108

              endif     

            go to 107    
      
  108    continue     
         enddo

        return
       end

c----------------------------------------------------------------------
      subroutine calcfetch(xavg,yavg,wdirstoke,fetch)
c       fetch calculation
c----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      parameter(npts=400000, imx=700, jmx=300)
      dimension alon(npts), alat(npts), seg(npts,5)

      character regn*4, dummy*80, a3*3
      logical ex
      glon(x)=along1+(x-1.d0)*dlong
      glat(y)=alatg1+(y-1.d0)*dlatg


      data regn /'medf'/  istart /1/
      pi = 4.d0 * datan(1.d0)
      degrad = pi / 180.d0

      m0=int(xavg+0.5d0)           ! centre of slick
      n0=int(yavg+0.5d0)
      xavg_lon=glon(dfloat(m0))
      yavg_lat=glat(dfloat(n0))
c
c     read map points
c
      open(100,file='makefetch.log')

      open(1,file='data/'//regn//'.map')
      read(1,*) ncontours

      k = 1
      nseg = 0
      do ni=1,ncontours
         read(1,*) isle
         do i=1,isle
            read(1,*) alon(i), alat(i) 
         enddo 
         if(isle.lt.50) go to 1

         seg(k,1) = alon(1) 
         seg(k,2) = alat(1) 
         do i=2,isle
            seg(k,3) = alon(i) 
            seg(k,4) = alat(i)

            k = k + 1
            seg(k,1) = alon(i) 
            seg(k,2) = alat(i) 
         enddo 
    1    continue
      enddo
      nseg = k - 1
c
      close(1)

      cs2 = dcos(yavg_lat * degrad) **2
      xi1  = xavg_lon
      eta1 = yavg_lat
             
      angledeg = wdirstoke
      angle = angledeg * degrad
      csangle = dcos(angle)
      snangle = dsin(angle)

      xi2 = xi1 + 50.d0 * snangle
      eta2 = eta1 + 50.d0 * csangle
      dmin = 100.d0
      nsmin = 0

      do 62 ns = 1,nseg
           
            xx1 = seg(ns,1)
            yy1 = seg(ns,2)
            xx2 = seg(ns,3)
            yy2 = seg(ns,4)

            ddel  = (xx2-xx1)*(eta2-eta1)-(yy2-yy1)*(xi2-xi1)
            ddel1 = (eta2-eta1)*(xi1-xx1)-(xi2-xi1)*(eta1-yy1)
            ddel2 = (yy2-yy1)*(xi1-xx1)-(xx2-xx1)*(eta1-yy1)
            if(ddel.eq.0.d0) go to 62
                 
            alam=ddel1/ddel
            alamp=ddel2/ddel
            if(alam.ge.0.d0.and.alam.le.1.d0.and.alamp.gt.0.d0) then
              xx=xx1+alam*(xx2-xx1)
              yy=yy1+alam*(yy2-yy1)
              dd1=dsqrt( cs2*(xx-xi1)*(xx-xi1) + (yy-eta1)*(yy-eta1) )
              if(dd1.lt.dmin) then
                dmin=dd1
                nsmn=ns
              end if  
            end if

   62 enddo    !ns
           
      fetch = dmin
      fetch = fetch*60*1852
      if(fetch.gt.20000) fetch=20000              


      return

      end subroutine calcfetch
