/* dpgrdr.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* $Procedure DPGRDR ( Derivative of planetographic w.r.t. rectangular ) */
/* Subroutine */ int dpgrdr_(char *body, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *re, doublereal *f, doublereal *jacobi, 
	ftnlen body_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rnge(char *, integer, 
	    char *, integer);

    /* Local variables */
    integer i__, n;
    extern /* Subroutine */ int chkin_(char *, ftnlen), ucase_(char *, char *,
	     ftnlen, ftnlen), errch_(char *, char *, ftnlen, ftnlen);
    logical found;
    extern /* Subroutine */ int errdp_(char *, doublereal *, ftnlen);
    integer sense;
    extern /* Subroutine */ int repmi_(char *, char *, integer *, char *, 
	    ftnlen, ftnlen, ftnlen), bods2c_(char *, integer *, logical *, 
	    ftnlen), dgeodr_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    integer bodyid;
    extern /* Subroutine */ int gcpool_(char *, integer *, integer *, integer 
	    *, char *, logical *, ftnlen, ftnlen);
    char kvalue[80];
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    char pmkvar[32], pgrlon[4];
    extern /* Subroutine */ int setmsg_(char *, ftnlen), cmprss_(char *, 
	    integer *, char *, char *, ftnlen, ftnlen, ftnlen);
    extern integer plnsns_(integer *);
    extern logical return_(void);
    char tmpstr[32];

/* $ Abstract */

/*     This routine computes the Jacobian matrix of the transformation */
/*     from rectangular to planetographic coordinates. */

/* $ Disclaimer */

/*     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE */
/*     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. */
/*     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE */
/*     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE */
/*     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" */
/*     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY */
/*     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A */
/*     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC */
/*     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE */
/*     SOFTWARE AND RELATED MATERIALS, HOWEVER USED. */

/*     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA */
/*     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT */
/*     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, */
/*     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, */
/*     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE */
/*     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY. */

/*     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF */
/*     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY */
/*     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE */
/*     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE. */

/* $ Required_Reading */

/*     None. */

/* $ Keywords */

/*     COORDINATES */
/*     DERIVATIVES */
/*     MATRIX */

/* $ Declarations */
/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     BODY       I   Body with which coordinate system is associated. */
/*     X          I   X-coordinate of point. */
/*     Y          I   Y-coordinate of point. */
/*     Z          I   Z-coordinate of point. */
/*     RE         I   Equatorial radius of the reference spheroid. */
/*     F          I   Flattening coefficient. */
/*     JACOBI     O   Matrix of partial derivatives. */

/* $ Detailed_Input */

/*     BODY       Name of the body with which the planetographic */
/*                coordinate system is associated. */

/*                BODY is used by this routine to look up from the */
/*                kernel pool the prime meridian rate coefficient giving */
/*                the body's spin sense.  See the Files and Particulars */
/*                header sections below for details. */

/*     X, */
/*     Y, */
/*     Z          are the rectangular coordinates of the point at */
/*                which the Jacobian of the map from rectangular */
/*                to planetographic coordinates is desired. */

/*     RE         Equatorial radius of the reference spheroid. */

/*     F          Flattening coefficient = (RE-RP) / RE,  where RP is */
/*                the polar radius of the spheroid.  (More importantly */
/*                RP = RE*(1-F).) */

/* $ Detailed_Output */

/*     JACOBI     is the matrix of partial derivatives of the conversion */
/*                from rectangular to planetographic coordinates.  It */
/*                has the form */

/*                    .-                               -. */
/*                    |  DLON/DX    DLON/DY   DLON/DZ   | */
/*                    |  DLAT/DX    DLAT/DY   DLAT/DZ   | */
/*                    |  DALT/DX    DALT/DY   DALT/DZ   | */
/*                    `-                               -' */

/*               evaluated at the input values of X, Y, and Z. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the body name BODY cannot be mapped to a NAIF ID code, */
/*        and if BODY is not a string representation of an integer, */
/*        the error SPICE(IDCODENOTFOUND) will be signaled. */

/*     2) If the kernel variable */

/*           BODY<ID code>_PGR_POSITIVE_LON */

/*        is present in the kernel pool but has a value other */
/*        than one of */

/*            'EAST' */
/*            'WEST' */

/*        the error SPICE(INVALIDOPTION) will be signaled.  Case */
/*        and blanks are ignored when these values are interpreted. */

/*     3) If polynomial coefficients for the prime meridian of BODY */
/*        are not available in the kernel pool, and if the kernel */
/*        variable BODY<ID code>_PGR_POSITIVE_LON is not present in */
/*        the kernel pool, the error SPICE(MISSINGDATA) will be signaled. */

/*     4) If the equatorial radius is non-positive, the error */
/*        SPICE(VALUEOUTOFRANGE) is signaled. */

/*     5) If the flattening coefficient is greater than or equal to one, */
/*        the error SPICE(VALUEOUTOFRANGE) is signaled. */

/*     6) If the input point is on the Z-axis (X = 0 and Y = 0), the */
/*        Jacobian matrix is undefined.  The error will be diagnosed */
/*        by routines in the call tree of this routine. */

/* $ Files */

/*     This routine expects a kernel variable giving BODY's prime */
/*     meridian angle as a function of time to be available in the */
/*     kernel pool.  Normally this item is provided by loading a PCK */
/*     file.  The required kernel variable is named */

/*        BODY<body ID>_PM */

/*     where <body ID> represents a string containing the NAIF integer */
/*     ID code for BODY.  For example, if BODY is 'JUPITER', then */
/*     the name of the kernel variable containing the prime meridian */
/*     angle coefficients is */

/*        BODY599_PM */

/*     See the PCK Required Reading for details concerning the prime */
/*     meridian kernel variable. */

/*     The optional kernel variable */

/*        BODY<body ID>_PGR_POSITIVE_LON */

/*     also is normally defined via loading a text kernel. When this */
/*     variable is present in the kernel pool, the prime meridian */
/*     coefficients for BODY are not required by this routine. See the */
/*     Particulars section below for details. */

/* $ Particulars */

/*     When performing vector calculations with velocities it is usually */
/*     most convenient to work in rectangular coordinates. However, once */
/*     the vector manipulations have been performed, it is often */
/*     desirable to convert the rectangular representations into */
/*     planetographic coordinates to gain insights about phenomena in */
/*     this coordinate frame. */

/*     To transform rectangular velocities to derivatives of coordinates */
/*     in a planetographic system, one uses the Jacobian of the */
/*     transformation between the two systems. */

/*     Given a state in rectangular coordinates */

/*        ( x, y, z, dx, dy, dz ) */

/*     the velocity in planetographic coordinates is given by the matrix */
/*     equation: */
/*                          t          |                     t */
/*        (dlon, dlat, dalt)   = JACOBI|       * (dx, dy, dz) */
/*                                     |(x,y,z) */

/*     This routine computes the matrix */

/*              | */
/*        JACOBI| */
/*              |(x, y, z) */


/*     The planetographic definition of latitude is identical to the */
/*     planetodetic (also called "geodetic" in SPICE documentation) */
/*     definition. In the planetographic coordinate system, latitude is */
/*     defined using a reference spheroid.  The spheroid is */
/*     characterized by an equatorial radius and a polar radius. For a */
/*     point P on the spheroid, latitude is defined as the angle between */
/*     the X-Y plane and the outward surface normal at P.  For a point P */
/*     off the spheroid, latitude is defined as the latitude of the */
/*     nearest point to P on the spheroid.  Note if P is an interior */
/*     point, for example, if P is at the center of the spheroid, there */
/*     may not be a unique nearest point to P. */

/*     In the planetographic coordinate system, longitude is defined */
/*     using the spin sense of the body.  Longitude is positive to the */
/*     west if the spin is prograde and positive to the east if the spin */
/*     is retrograde.  The spin sense is given by the sign of the first */
/*     degree term of the time-dependent polynomial for the body's prime */
/*     meridian Euler angle "W":  the spin is retrograde if this term is */
/*     negative and prograde otherwise.  For the sun, planets, most */
/*     natural satellites, and selected asteroids, the polynomial */
/*     expression for W may be found in a SPICE PCK kernel. */

/*     The earth, moon, and sun are exceptions: planetographic longitude */
/*     is measured positive east for these bodies. */

/*     If you wish to override the default sense of positive longitude */
/*     for a particular body, you can do so by defining the kernel */
/*     variable */

/*        BODY<body ID>_PGR_POSITIVE_LON */

/*     where <body ID> represents the NAIF ID code of the body. This */
/*     variable may be assigned either of the values */

/*        'WEST' */
/*        'EAST' */

/*     For example, you can have this routine treat the longitude */
/*     of the earth as increasing to the west using the kernel */
/*     variable assignment */

/*        BODY399_PGR_POSITIVE_LON = 'WEST' */

/*     Normally such assignments are made by placing them in a text */
/*     kernel and loading that kernel via FURNSH. */

/*     The definition of this kernel variable controls the behavior of */
/*     the SPICELIB planetographic routines */

/*        PGRREC */
/*        RECPGR */
/*        DPGRDR */
/*        DRDPGR */

/*     It does not affect the other SPICELIB coordinate conversion */
/*     routines. */

/* $ Examples */

/*     Numerical results shown for this example may differ between */
/*     platforms as the results depend on the SPICE kernels used as */
/*     input and the machine specific arithmetic implementation. */


/*         Find the planetographic state of the earth as seen from */
/*         Mars in the J2000 reference frame at January 1, 2005 TDB. */
/*         Map this state back to rectangular coordinates as a check. */


/*              PROGRAM EX1 */
/*              IMPLICIT NONE */
/*        C */
/*        C     SPICELIB functions */
/*        C */
/*              DOUBLE PRECISION      RPD */
/*        C */
/*        C     Local variables */
/*        C */
/*              DOUBLE PRECISION      ALT */
/*              DOUBLE PRECISION      DRECTN ( 3 ) */
/*              DOUBLE PRECISION      ET */
/*              DOUBLE PRECISION      F */
/*              DOUBLE PRECISION      JACOBI ( 3, 3 ) */
/*              DOUBLE PRECISION      LAT */
/*              DOUBLE PRECISION      LON */
/*              DOUBLE PRECISION      LT */
/*              DOUBLE PRECISION      PGRVEL ( 3 ) */
/*              DOUBLE PRECISION      RADII  ( 3 ) */
/*              DOUBLE PRECISION      RE */
/*              DOUBLE PRECISION      RECTAN ( 3 ) */
/*              DOUBLE PRECISION      RP */
/*              DOUBLE PRECISION      STATE  ( 6 ) */

/*              INTEGER               N */
/*        C */
/*        C     Load a PCK file containing a triaxial */
/*        C     ellipsoidal shape model and orientation */
/*        C     data for Mars. */
/*        C */
/*              CALL FURNSH ( 'pck00008.tpc' ) */

/*        C */
/*        C     Load an SPK file giving ephemerides of earth and Mars. */
/*        C */
/*              CALL FURNSH ( 'de405.bsp' ) */

/*        C */
/*        C     Load a leapseconds kernel to support time conversion. */
/*        C */
/*              CALL FURNSH ( 'naif0007.tls' ) */

/*        C */
/*        C     Look up the radii for Mars.  Although we */
/*        C     omit it here, we could first call BADKPV */
/*        C     to make sure the variable BODY499_RADII */
/*        C     has three elements and numeric data type. */
/*        C     If the variable is not present in the kernel */
/*        C     pool, BODVRD will signal an error. */
/*        C */
/*              CALL BODVRD ( 'MARS', 'RADII', 3, N, RADII ) */

/*        C */
/*        C     Compute flattening coefficient. */
/*        C */
/*              RE  =  RADII(1) */
/*              RP  =  RADII(3) */
/*              F   =  ( RE - RP ) / RE */

/*        C */
/*        C     Look up the geometric state of earth as seen from Mars at */
/*        C     January 1, 2005 TDB, relative to the J2000 reference */
/*        C     frame. */
/*        C */
/*              CALL STR2ET ( 'January 1, 2005 TDB', ET ) */

/*              CALL SPKEZR ( 'Earth', ET,    'J2000', 'LT+S', */
/*             .              'Mars',  STATE, LT               ) */

/*        C */
/*        C     Convert position to planetographic coordinates. */
/*        C */
/*              CALL RECPGR ( 'MARS', STATE, RE, F, LON, LAT, ALT ) */

/*        C */
/*        C     Convert velocity to planetographic coordinates. */
/*        C */

/*              CALL DPGRDR ( 'MARS', STATE(1), STATE(2), STATE(3), */
/*             .               RE,    F,        JACOBI             ) */

/*              CALL MXV ( JACOBI, STATE(4), PGRVEL ) */

/*        C */
/*        C     As a check, convert the planetographic state back to */
/*        C     rectangular coordinates. */
/*        C */
/*              CALL PGRREC ( 'MARS', LON, LAT, ALT, RE, F, RECTAN ) */

/*              CALL DRDPGR ( 'MARS', LON, LAT, ALT, RE, F, JACOBI ) */

/*              CALL MXV ( JACOBI, PGRVEL, DRECTN ) */


/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Rectangular coordinates:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  X (km)                 = ', STATE(1) */
/*              WRITE(*,*) '  Y (km)                 = ', STATE(2) */
/*              WRITE(*,*) '  Z (km)                 = ', STATE(3) */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Rectangular velocity:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  dX/dt (km/s)           = ', STATE(4) */
/*              WRITE(*,*) '  dY/dt (km/s)           = ', STATE(5) */
/*              WRITE(*,*) '  dZ/dt (km/s)           = ', STATE(6) */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Ellipsoid shape parameters: ' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  Equatorial radius (km) = ', RE */
/*              WRITE(*,*) '  Polar radius      (km) = ', RP */
/*              WRITE(*,*) '  Flattening coefficient = ', F */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Planetographic coordinates:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  Longitude (deg)        = ', LON / RPD() */
/*              WRITE(*,*) '  Latitude  (deg)        = ', LAT / RPD() */
/*              WRITE(*,*) '  Altitude  (km)         = ', ALT */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Planetographic velocity:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  d Longitude/dt (deg/s) = ', PGRVEL(1)/RPD() */
/*              WRITE(*,*) '  d Latitude/dt  (deg/s) = ', PGRVEL(2)/RPD() */
/*              WRITE(*,*) '  d Altitude/dt  (km/s)  = ', PGRVEL(3) */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Rectangular coordinates from inverse ' // */
/*             .           'mapping:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  X (km)                 = ', RECTAN(1) */
/*              WRITE(*,*) '  Y (km)                 = ', RECTAN(2) */
/*              WRITE(*,*) '  Z (km)                 = ', RECTAN(3) */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) 'Rectangular velocity from inverse mapping:' */
/*              WRITE(*,*) ' ' */
/*              WRITE(*,*) '  dX/dt (km/s)           = ', DRECTN(1) */
/*              WRITE(*,*) '  dY/dt (km/s)           = ', DRECTN(2) */
/*              WRITE(*,*) '  dZ/dt (km/s)           = ', DRECTN(3) */
/*              WRITE(*,*) ' ' */
/*              END */


/*        Output from this program should be similar to the following */
/*        (rounding and formatting differ across platforms): */


/*           Rectangular coordinates: */

/*             X (km)                 =   146039732. */
/*             Y (km)                 =   278546607. */
/*             Z (km)                 =   119750315. */

/*           Rectangular velocity: */

/*             dX/dt (km/s)           =  -47.0428824 */
/*             dY/dt (km/s)           =   9.07021778 */
/*             dZ/dt (km/s)           =   4.75656274 */

/*           Ellipsoid shape parameters: */

/*             Equatorial radius (km) =   3396.19 */
/*             Polar radius      (km) =   3376.2 */
/*             Flattening coefficient =   0.00588600756 */

/*           Planetographic coordinates: */

/*             Longitude (deg)        =   297.667659 */
/*             Latitude  (deg)        =   20.844504 */
/*             Altitude  (km)         =   336531825. */

/*           Planetographic velocity: */

/*             d Longitude/dt (deg/s) =  -8.35738632E-06 */
/*             d Latitude/dt  (deg/s) =   1.59349355E-06 */
/*             d Altitude/dt  (km/s)  =  -11.2144327 */

/*           Rectangular coordinates from inverse mapping: */

/*             X (km)                 =   146039732. */
/*             Y (km)                 =   278546607. */
/*             Z (km)                 =   119750315. */

/*           Rectangular velocity from inverse mapping: */

/*             dX/dt (km/s)           =  -47.0428824 */
/*             dY/dt (km/s)           =   9.07021778 */
/*             dZ/dt (km/s)           =   4.75656274 */


/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */
/*     W.L. Taber     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 26-DEC-2004 (NJB) (WLT) */

/* -& */
/* $ Index_Entries */

/*     Jacobian of planetographic  w.r.t. rectangular coordinates */

/* -& */
/* $ Revisions */

/*     None. */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DPGRDR", (ftnlen)6);
    }

/*     Convert the body name to an ID code. */

    bods2c_(body, &bodyid, &found, body_len);
    if (! found) {
	setmsg_("The value of the input argument BODY is #, this is not a re"
		"cognized name of an ephemeris object. The cause of this prob"
		"lem may be that you need an updated version of the SPICE Too"
		"lkit. ", (ftnlen)185);
	errch_("#", body, (ftnlen)1, body_len);
	sigerr_("SPICE(IDCODENOTFOUND)", (ftnlen)21);
	chkout_("DPGRDR", (ftnlen)6);
	return 0;
    }

/*     The equatorial radius must be positive. If not, signal an error */
/*     and check out. */

    if (*re <= 0.) {
	setmsg_("Equatorial radius was #.", (ftnlen)24);
	errdp_("#", re, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("DPGRDR", (ftnlen)6);
	return 0;
    }

/*     If the flattening coefficient is greater than 1, the polar radius */
/*     is negative. If F is equal to 1, the polar radius is zero. Either */
/*     case is a problem, so signal an error and check out. */

    if (*f >= 1.) {
	setmsg_("Flattening coefficient was #.", (ftnlen)29);
	errdp_("#", f, (ftnlen)1);
	sigerr_("SPICE(VALUEOUTOFRANGE)", (ftnlen)22);
	chkout_("DPGRDR", (ftnlen)6);
	return 0;
    }

/*     Look up the longitude sense override variable from the */
/*     kernel pool. */

    repmi_("BODY#_PGR_POSITIVE_LON", "#", &bodyid, pmkvar, (ftnlen)22, (
	    ftnlen)1, (ftnlen)32);
    gcpool_(pmkvar, &c__1, &c__1, &n, kvalue, &found, (ftnlen)32, (ftnlen)80);
    if (found) {

/*        Make sure we recognize the value of PGRLON. */

	cmprss_(" ", &c__0, kvalue, tmpstr, (ftnlen)1, (ftnlen)80, (ftnlen)32)
		;
	ucase_(tmpstr, pgrlon, (ftnlen)32, (ftnlen)4);
	if (s_cmp(pgrlon, "EAST", (ftnlen)4, (ftnlen)4) == 0) {
	    sense = 1;
	} else if (s_cmp(pgrlon, "WEST", (ftnlen)4, (ftnlen)4) == 0) {
	    sense = -1;
	} else {
	    setmsg_("Kernel variable # may have the values EAST or WEST.  Ac"
		    "tual value was #.", (ftnlen)72);
	    errch_("#", pmkvar, (ftnlen)1, (ftnlen)32);
	    errch_("#", kvalue, (ftnlen)1, (ftnlen)80);
	    sigerr_("SPICE(INVALIDOPTION)", (ftnlen)20);
	    chkout_("DPGRDR", (ftnlen)6);
	    return 0;
	}
    } else {

/*        Look up the spin sense of the body's prime meridian. */

	sense = plnsns_(&bodyid);

/*        If the required prime meridian rate was not available, */
/*        PLNSNS returns the code 0.  Here we consider this situation */
/*        to be an error. */

	if (sense == 0) {
	    repmi_("BODY#_PM", "#", &bodyid, pmkvar, (ftnlen)8, (ftnlen)1, (
		    ftnlen)32);
	    setmsg_("Prime meridian rate coefficient defined by kernel varia"
		    "ble # is required but not available for body #. ", (
		    ftnlen)103);
	    errch_("#", pmkvar, (ftnlen)1, (ftnlen)32);
	    errch_("#", body, (ftnlen)1, body_len);
	    sigerr_("SPICE(MISSINGDATA)", (ftnlen)18);
	    chkout_("DPGRDR", (ftnlen)6);
	    return 0;
	}

/*        Handle the special cases:  earth, moon, and sun. */

	if (bodyid == 399 || bodyid == 301 || bodyid == 10) {
	    sense = 1;
	}
    }

/*     At this point, SENSE is set to +/- 1. */

/*     To obtain the Jacobian matrix we want, first find the Jacobian */
/*     matrix of rectangular coordinates with respect to geodetic */
/*     coordinates. */

    dgeodr_(x, y, z__, re, f, jacobi);

/*     Letting GLON represent geodetic longitude, the matrix JACOBI is */

/*        .-                             -. */
/*        |  DGLON/DX  DGLON/DY  DGLON/DZ | */
/*        |  DLAT/DX   DLAT/DY   DLAT/DZ  | */
/*        |  DALT/DX   DALT/DY   DALT/DZ  | */
/*        `-                             -' */

/*     evaluated at the input values of X, Y, and Z. */

/*     Since planetographic longitude LON satisfies */

/*        LON = SENSE * GLON */

/*     applying the chain rule to D(*)/DGLON, the above is equivalent to */

/*        .-                                                         -. */
/*        |  (1/SENSE)*DLON/DX  (1/SENSE)*DLON/DY  (1/SENSE)*DLON/DZ  | */
/*        |            DLAT/DX            DLAT/DY            DLAT/DZ  | */
/*        |            DALT/DX            DALT/DY            DALT/DZ  | */
/*        `-                                                         -' */

/*     So, multiplying the first row of JACOBI by SENSE gives us the */
/*     matrix we actually want to compute:  the Jacobian matrix of */
/*     rectangular coordinates with respect to planetographic */
/*     coordinates. */

    for (i__ = 1; i__ <= 3; ++i__) {
	jacobi[(i__1 = i__ * 3 - 3) < 9 && 0 <= i__1 ? i__1 : s_rnge("jacobi",
		 i__1, "dpgrdr_", (ftnlen)712)] = sense * jacobi[(i__2 = i__ *
		 3 - 3) < 9 && 0 <= i__2 ? i__2 : s_rnge("jacobi", i__2, 
		"dpgrdr_", (ftnlen)712)];
    }
    chkout_("DPGRDR", (ftnlen)6);
    return 0;
} /* dpgrdr_ */

