!#include <misc.h>
!#include <params.h>



c--------------------------------------------------------------------
	subroutine sunrisesetxx( xlong, ylat, iyear, imonth, iday,
     +		iflag, trise, tset, solardec )
c
c   provides interface to subr srisesetxx without use of common blocks
c
c   input parameters
c	xlong - longitude in degrees (east longitudes are positive)
c	ylat  - latitude in degrees (north latitudes are positive)
c	iyear - year
c	imonth - month
c	iday - day
c   output parameters
c	iflag - status flag
c	    +1 - OK and there is a sunrise and sunset
c	     0 - OK but no sunrise or sunset
c	    -1 = input parameters (date or position) are bad
c	trise - time of sunrise in UT hours
c	tset  - time of sunset  in UT hours
c	solardec - apparent solar declination in degrees
c
c   written 17-aug-93 by r.c.easter

        use shr_kind_mod, only: r8 => shr_kind_r8
	implicit none

c   arguments
	integer iyear, imonth, iday, iflag
	real(r8) xlong, ylat, trise, tset, solardec
	
	integer jday, iriseset
	real(r8) sunrise, sunset, rloc_timehrs, ap_dec
!	common / sundataxx_cmn / jday, iriseset, 
!     +		sunrise, sunset, rloc_timehrs, ap_dec

c   local
	integer i
	real(r8) xlongb


c   need xlong between -180 and +180
	xlongb = xlong
	do 1000 i = 1, 10
	    if (xlongb .lt. -180._r8) then
		xlongb = xlongb + 360._r8
	    else if (xlongb .gt. 180._r8) then
		xlongb = xlongb - 360._r8
	    else
		goto 1050
	    end if
1000	continue
1050	continue
	
	call srisesetxx( iyear, imonth, iday, ylat, xlongb, 
     &   jday, iriseset,sunrise, sunset, rloc_timehrs, ap_dec)

	iflag = iriseset
	if (iflag .eq. 0) then
	    iflag = 1
	    if (abs(sunrise+100._r8) .le. 0.01_r8) iflag = 0
	end if
	trise = sunrise
	tset  = sunset
	solardec = ap_dec

	return
	end


c--------------------------------------------------------------------
	subroutine solarzangxx( xlong, ylat, iyear, imonth, iday,
     +		thours, zenithang )
c
c   provides interface to subr solarxx without use of common blocks
c
c   input parameters
c	xlong - longitude in degrees (east longitudes are positive)
c	ylat  - latitude in degrees (north latitudes are positive)
c	iyear - year
c	imonth - month
c	iday - day
c	thours - time in hours from UT midnight
c   output parameter
c	zenithang - solar zenith angle in degrees
c
c   written 17-aug-93 by r.c.easter
c
        use shr_kind_mod, only: r8 => shr_kind_r8
	implicit none

c   arguments
	integer iyear, imonth, iday
	real(r8) xlong, ylat, thours, zenithang

	integer iy, im, id, idaytime
	real sla, slo, zcoef, zangle, altm
	common /solarxx_cmn/ sla, slo, iy, im, id, zcoef, zangle,
     +		idaytime, altm

c   local
!	real thours


	slo = -xlong
	sla = ylat
	iy = iyear
	im = imonth
	id = iday
	call solarxx( thours )

	zenithang = zangle

	return
	end




c*********************************************************************
	subroutine solarxx (timehrs)
c based on subr sunae of joe michaelsky
c added statements and subr daynum to produce quantities
c needed for joe's routine.  also added statements to yield solar 
c zenith angle (equal to the complement of the elevation) and cosine
c of this angle. it is assumed that the input timehrs refers to 
c greenwich mean time (gmt).
c eg chapman  feb 92
c---------------------------------------------------------------------
c original commenting of joe michaelsky:
c this subr calculates the local azimuth and elevation of the sun
c at a specified location and time using an approximation to the equations
c used to generate the nautical almanac.  refraction correction added
c so sun position is apparent one.  if sun more than 0.5 deg below horizon
c returns e uncorrected for refraction.
c references:  solar energy 20,393; 22,195; 25,91; 27,67
c
c input parameters
c	year=the year number (eg, 1981)
c	day=the day number (eg, jan 1 = 1)
c	tm=time in sec from ut midnight
c	lat=local latitude in degrees
c	long=local longitude in degrees
c output parameters
c	azimuth=azimuth angle of sun (positive is east from north)
c	e=elevation of sun
c	h=hour angle in deg-negative before meridian (equal to suhrang)
c	decl=dec of sun (equal to sundecln)
c additional outputs added by eg chapman
c	zangle = zenith angle, equal to complement of elevation angle
c	zcoef = cosine of zenith angle
c
c-----------------------------------------------------------------------
	real timehrs

	integer iy, im, id, idaytime
	real sla, slo, zcoef, zangle, altm
	common /solarxx_cmn/ sla, slo, iy, im, id, zcoef, zangle,
     +		idaytime, altm

	real lat, long
	data twopi, rad /6.2831853, 0.017453293/

c compute day number and time in sec from ut midnight from available
c model variables
	tm = timehrs * 3600.
	call daynumxx (iy, im, id, jday)
	day = real(jday)
	year = real (iy)     
	lat = sla
	long = slo
c preserve original day number
	day0 = day
	delyr = year - 1980.
	leap = aint(delyr/4.)
	t =tm/3600.
	time = delyr * 365. + leap + day -1. + t/24.
	if (delyr.eq.leap*4.) time = time -1.
	if (delyr.lt.0. .and. delyr.ne.leap*4.)time=time-1.
	theta = (360.*time /365.25) * rad
	g = -0.031271 -4.53963e-7 * time + theta
	el = 4.900968 + 3.67474e-7 * time + (0.033434 - 2.3e-9*time)
     +   * sin(g) + 0.000349 * sin(2.*g) + theta
	eps = 0.409140 - 6.2149e-9 * time
	sel = sin(el)
	a1 = sel * cos(eps)
	a2 = cos(el)
	rangle = atan2(a1,a2)
	if (rangle.lt.0.) rangle = rangle + twopi
	decl = asin (sel * sin(eps))
	st = 1.759335 + twopi * (time/365.25 -delyr) + 3.694e-7 *time
	if (st.ge.twopi)st=st-twopi
	s = st +(t* 15. -long) * rad
	if (s.ge.twopi)s = s -twopi
	h = rangle - s
	phi = lat * rad
	q = sin(phi) * sin(decl) + cos(phi) * cos(decl) * cos(h)
	qq = asin (q)
	if (qq.lt. -0.00872665) then
		e = qq
	        argue = 0.0
		go to 20
	end if
	rc = (1./(0.955 + (20.267 * q))) - 0.047121
	qrc = q + (0.0083 * rc)
	e = asin (qrc)
cec following 4 statements added 5/28/92 by eg chapman to handle cases
cec near equinoxes, where roundoff leads to fatal math
cec error for trying to compute arc sin of a quantity greater than 1.
	argue = (cos(decl) * sin (h)) / cos(qq)
	if (argue .gt. 1.00) argue = 1.00
	if (argue .lt. -1.00) argue = -1.00
  20	azimuth = asin (argue) / rad
c 20 	azimuth = asin(cos(decl) * sin(h) /cos(qq))/rad
c following statement added 2/24/92 by eg chapman to handle case where
c latitude is 0 degrees
	if (lat .eq. 0.0) phi = 0.000001 * rad
	if (sin(e) .ge. sin(decl)/sin(phi)) go to 10
	if (azimuth .lt. 0.) azimuth = azimuth +360.
	azimuth = 180. - azimuth 
 10	e = e /rad
	azimuth = 180. - azimuth 
	decl = decl / rad
	h = h/rad
	if (h .lt. -180.) h = 360 + h
	h = -h
c return original day number
	day = day0
	sundecln = decl
	sunhrang = h
	zangle = 90. - e
	zcoef = cosd (zangle)
	return
	end


c***************************************************************************
	subroutine srisesetxx(iyear, month, iday, rlat, rlong,
     &   jday, iriseset,sunrise, sunset, rloc_timehrs, ap_dec)

        use shr_kind_mod, only: r8 => shr_kind_r8

	integer iyear, month, iday
	real (r8) rlat, rlong

	integer jday, iriseset
	real (r8) sunrise, sunset, rloc_timehrs,apdec
!	common / sundataxx_cmn / jday, iriseset, 
!     +		sunrise, sunset, rloc_timehrs, ap_dec
c--------------------------------------------------------------------------
c "sunriseset.c" contains the integer function sunriseset() for calculating
c  the rising and setting times of the Sun as seen from a place on Earth on a
c  specific date.
c
c  Version 1.0 - April 6, 1992.
c  (This code was adapted from "solarpos.c" Version 3.1.)
c
c  sunriseset() employs the low precision formulas for the Sun's coordinates
c  given in the "Astronomical Almanac" of 1990 to compute the Sun's apparent
c  right ascension, apparent declination, and Greenwich mean sidereal time at
c  0 hours Universal Time, and then the rising and setting times of the Sun.
c  The "Astronomical Almanac" (A. A.) states a precision of 0.01 degree for the
c  apparent coordinates between the years 1950 and 2050.
c
c  The following assumptions and simplifications are made:
c  -> diurnal parallax is ignored, resulting in 0 to 9 arc seconds error in
c     apparent position.
c  -> diurnal aberration is also ignored, resulting in 0 to 0.02 second error
c     in right ascension and 0 to 0.3 arc second error in declination.
c  -> geodetic site coordinates are used, without correction for polar motion
c     (maximum amplitude of 0.3 arc second) and local gravity anomalies.
c  -> the formulas ignore atmospheric refraction, semi-diameter, and changes
c     in right ascension and declination over the course of a day; the
c     accuracies of sunrise and sunset are about 2 and 7 minutes for latitude
c     and longitude of 0 degrees, but accuracy degrades significantly for high
c     latitudes.
c
c
c  The necessary input parameters are:
c  -> the UT date, specified in one of three ways:
c       1) year, month, day.fraction
c       2) year, daynumber.fraction
c       3) days.fraction elapsed since January 0, 1900.
c  Note: in GChM application, only specification #1 is currently valid
c  -> site geodetic (geographic) latitude and longitude.
c
c  Refer to the function declaration for the parameter type specifications and
c  formats.
c
c  sunriseset() returns -1 if an input parameter is out of bounds, or 0 if
c  values were written to the locations specified by the output parameters.
c  Sunrise and sunset times are in UT hours; if there is no sunrise or sunset
c  the values are -1.0.
c
c  Author: Nels Larson
c          Pacific Northwest Lab.
c          P.O. Box 999
c          Richland, WA 99352
c          U.S.A.
c
c--------------------------------------------------------------------------
c modifications for gchm application by eg chapman
c	1. translated from c language to fortran
c	2. input date must be in year, month, day.fraction format; other input 
c	   code eliminated.              
c	3. added indicator iriseset. when equal to -1, indicates location
c	   or date is out of range. 
c	
c---------------------------------------------------------------------------
c math definitions. 
 	parameter (twopi = 6.2831853071795864_r8)
	parameter (deg_rad = 0.017453292519943295_r8)
	parameter (rad_deg = 57.295779513082323_r8)
	real mean_anomaly, mean_longitude, mean_obliquity
c-------------------------------------------------------------------------
c explanation of terms taken from c code
c int    iyear,         Four digit year (Gregorian calendar).
c                       [1950 through 2049; 0 if using days_1900] 
c        month;        Month number.
c                       [1 through 12; 0 if using daynumber for day] 
c 
c day,           /* Calendar day.fraction, or daynumber.fraction.
c                       *   [If month is NOT 0:
c
c                       *      0 through 32; 31st @ 18:10:00 UT = 31.75694
c                       *    If month IS 0:
c                       *      0 through 367; 366 @ 18:10:00 UT = 366.75694] */
c       days_1900,     /* Days since 1900 January 0 @ 00:00:00 UT.
c                       *   [18262.0 (1950/01/00) through 54788.0 (2049/12/32);
c
c                       *    1990/01/01 @ 18:10:00 UT = 32873.75694;
c                       *    0.0 o.k. if using {year, month, day} or
c                       *    {year, daynumber}] */
c       rlat       Observation site geographic latitude.
c                  [degrees.fraction, North positive] 
c       rlong      Observation site geographic longitude.
c                  [degrees.fraction, East positive] 
c       *ap_ra,        /* Apparent solar right ascension.
c                       *   [hours; 0.0 <= *ap_ra < 24.0] */
c       *ap_dec,       /* Apparent solar declination.
c                       *   [degrees; -90.0 <= *ap_dec <= 90.0] */
c
c       *sunrise,      /* Time of sunrise.
c                           [UT hours.fraction; -1.0 if no sunrise or sunset] */
c       *sunset;       /* Time of sunset.
c                           [UT hours.fraction; -1.0 if no sunset or sunrise] */
c  int    daynum();        /* Computes a sequential daynumber during a year. */
c  int    daynumber,       /* Sequential daynumber during a year. */
c         delta_days,      /* Whole days since 2000 January 0. */
c         delta_years;     /* Whole years since 2000. */
c  double cent_J2000,      /* Julian centuries since epoch J2000.0 at 0h UT. */
c         days_J2000,      /* Days since epoch J2000.0. */
c         ecliptic_long,   /* Solar ecliptic longitude. */
c
c         gmst0h,          /* Greenwich mean sidereal time at 0 hours UT. */
c         integral,        /* Integral portion of double precision number. */
c         mean_anomaly,    /* Earth mean anomaly. */
c         mean_longitude,  /* Solar mean longitude. */
c         mean_obliquity,  /* Mean obliquity of the ecliptic. */
c         tan_dec,         /* Tangent of apparent declination. */
c         tan_lat,         /* Tangent of latitude. */
c
c         tangterm,        /* Tangent term of Sun rise/set equation. */
c         timeterm;        /* Time term of Sun rise/set equation. */
c----------------------------------------------------------------------
	iriseset = 0
c check latitude, longitude, dates for proper range before calculating dates.
	 if (((rlat .lt. -90.0_r8) .or. (rlat .gt. 90.0_r8)) .or. 
     +  ((rlong .lt. -180.0_r8) .or. (rlong .gt. 180.0_r8))) then
		iriseset = -1
		return
	end if
  	if (iyear .lt. 1950 .or. iyear .gt. 2049) then
		iriseset = -1
		return
	end if
        if (((month .lt. 1) .or. (month .gt. 12)) .or. 
     +      ((iday .lt. 0) .or. (iday .gt. 32))) then
     		iriseset = -1    
		return
	end if
c determine julian day number
	call daynumxx(iyear, month, iday, jday)
c construct Julian centuries since J2000 at 0 hours UT of date,
c days.fraction since J2000, and UT hours.
 	   delta_years = iyear - 2000._r8
c delta_days is days from 2000/01/00 (1900's are negative). 
  	 delta_days = delta_years * 365._r8 + delta_years / 4. + jday
    	if (iyear .gt. 2000) delta_days = delta_days + 1._r8
c J2000 is 2000/01/01.5 
  	  days_j2000 = delta_days - 1.5_r8
    	cent_j2000 = days_j2000 / 36525.0_r8
c compute solar position parameters.
c    A. A. 1990, C24.
  	f_mean_anomaly = (357.528_r8 + 0.9856003_r8 * days_j2000)
  	f_mean_longitude = (280.460_r8 + 0.9856474_r8 * days_j2000)
c put mean_anomaly and mean_longitude in the range 0 -> 2 pi. 
	mean_anomaly = (f_mean_anomaly / 360._r8 - int(f_mean_anomaly
     +   /360._r8)) * twopi
	mean_longitude = (f_mean_longitude /360. - int(
     +   f_mean_longitude/360._r8)) * twopi
  	mean_obliquity = (23.439_r8 - 4.0e-7_r8 * days_j2000) * deg_rad 
  	ecliptic_long = ((1.915_r8 * sin(mean_anomaly)) +
     +              (0.020_r8 * sin(2.0_r8 * mean_anomaly))) * deg_rad +
     +             mean_longitude
c tangent of ecliptic_long separated into sine and cosine parts for ap_ra. 
  	f_ap_ra = atan2(cos(mean_obliquity) * sin(ecliptic_long), 
     +            cos(ecliptic_long))
c change range of ap_ra from -pi -> pi to 0 -> 2 pi. 
  	if (f_ap_ra .lt. 0.0_r8) f_ap_ra = f_ap_ra + twopi
c put ap_ra in the range 0 -> 24 hours. 
	ap_ra = (f_ap_ra / twopi - int(f_ap_ra /twopi)) * 24.0_r8
  	ap_dec = asin(sin(mean_obliquity) * sin(ecliptic_long))
c calculate local mean sidereal time.
c A. A. 1990, B6-B7.
c horner's method of polynomial exponent expansion used for gmst0h. 
	f_gmst0h = 24110.54841_r8 + cent_j2000 * (8640184.812866_r8 
     +     +cent_j2000 * (0.093104_r8 - cent_j2000 * 6.2e-6_r8))
c convert gmst0h from seconds to hours and put in the range 0 -> 24. 
c 24 hours = 86400 seconds 
	gmst0h = (f_gmst0h / 86400._r8 - int(f_gmst0h / 86400._r8)) * 24.0_r8
  	if (gmst0h .lt. 0.0_r8) gmst0h = gmst0h + 24.0_r8
c convert latitude to radians. 
  	rlat_r =  rlat * deg_rad 
c avoid tangent overflow at +-90 degrees.
c 1.57079615 radians is equal to 89.99999 degrees.
  	if (abs(rlat_r) .lt. 1.57079615_r8) then
    		tan_lat = tan(rlat_r)
  	else
    		tan_lat = 6.0e6_r8
	end if
  	if (abs(ap_dec) .lt. 1.57079615_r8) then
    		tan_dec = tan(ap_dec)
  	else
    		tan_dec = 6.0e6_r8
	end if
c compute UTs of sunrise and sunset.
c A. A. 1990, A12.
  	tangterm = tan_lat * tan_dec
  	if (abs(tangterm) .gt. 1.0_r8) then
    		sunrise = -100.0_r8
    		sunset = -100.0_r8
  	else
c compute angle of tangterm and convert to hours. 
    		tangterm = acos(-tangterm) / twopi * 24.0_r8
    		timeterm = ap_ra - rlong / 15.0_r8 - gmst0h
    		sunrise = timeterm - tangterm
    		sunset = timeterm + tangterm
c put sunrise and sunset in the range 0 to 24 hours. 
cec inserted following statement since in some latitudes timeterm
cec minus tangterm is less than -25
		if (sunrise .le. -24.0_r8) sunrise = sunrise + 48.0_r8
    		if (sunrise .lt. 0.0_r8) sunrise = sunrise + 24.0_r8
    		if (sunrise .ge. 24.0_r8) sunrise = sunrise - 24.0_r8
    		if (sunset .lt. 0.0_r8) sunset = sunset + 24.0_r8
    		if (sunset .ge. 24.0_r8) sunset = sunset - 24.0_r8
c mean sidereal day is 0.99727 mean solar days. 
    		sunrise = sunrise * 0.99727_r8
		sunset =  sunset * 0.99727_r8
  	end if 
c convert ap_dec to degrees. 
  	ap_dec = ap_dec * rad_deg 
	return
	end


c*************************************************************************
	subroutine daynumxx (iyear, month, iday, jday)
c 'daynum()' returns the sequential daynumber of a calendar date during a
c  Gregorian calendar year (for years 1 onward).
c  The integer arguments are the four-digit year, the month number, and
c  the day of month number.
c  (Jan. 1 = 01/01 = 001; Dec. 31 = 12/31 = 365 or 366.)
c  A value of -1 is returned if the year is out of bounds.
c Author: Nels Larson
c         Pacific Northwest Lab.
c         P.O. Box 999
c         Richland, WA 99352
c         U.S.A.
c converted to fortran from c by eg chapman 5-92
c------------------------------------------------------------------------
	integer iyear, month, iday, jday

	dimension iimonth(12)
	data iimonth /0,31,59,90,120,151,181,212,243,273,304,334/
        leapyr = 0
c there is no year 0 in the Gregorian calendar and the leap year cycle
c changes for earlier years. 
	if (iyear .lt. 1) then
		iriseset = -1
		return
	end if
c leap years are divisible by 4, except for centurial years not divisible
c by 400. 
	year = real (iyear)
	if ((amod(year,4.) .eq. 0.0) .and. (amod(year,100.) .ne. 0.0)) 
     +	  leapyr = 1
	if(amod(year,400.) .eq. 0.0) leapyr = 1 
  	jday = iimonth(month) + iday
	if ((leapyr .eq. 1) .and. (month .gt. 2)) jday = jday + 1
  	return 
	end





