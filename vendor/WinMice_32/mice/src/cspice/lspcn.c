/* lspcn.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__10 = 10;
static integer c__3 = 3;
static integer c__2 = 2;

/* $Procedure    LSPCN  ( Longitude of the sun, planetocentric ) */
doublereal lspcn_(char *body, doublereal *et, char *abcorr, ftnlen body_len, 
	ftnlen abcorr_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    doublereal tipm[9]	/* was [3][3] */;
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    logical found;
    doublereal uavel[3], npole[3], trans[9]	/* was [3][3] */;
    extern /* Subroutine */ int ucrss_(doublereal *, doublereal *, doublereal 
	    *), bods2c_(char *, integer *, logical *, ftnlen);
    extern logical failed_(void);
    integer idcode;
    doublereal lt;
    extern /* Subroutine */ int recrad_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), tipbod_(char *, integer *, 
	    doublereal *, doublereal *, ftnlen);
    doublereal bstate[6], radius;
    extern /* Subroutine */ int spkgeo_(integer *, doublereal *, char *, 
	    integer *, doublereal *, doublereal *, ftnlen), sigerr_(char *, 
	    ftnlen), chkout_(char *, ftnlen), setmsg_(char *, ftnlen);
    doublereal sstate[6];
    extern /* Subroutine */ int twovec_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *), spkezr_(char *, doublereal *, char *, 
	    char *, char *, doublereal *, doublereal *, ftnlen, ftnlen, 
	    ftnlen, ftnlen);
    extern logical return_(void);
    doublereal lat, pos[3];
    extern /* Subroutine */ int mxv_(doublereal *, doublereal *, doublereal *)
	    ;

/* $ Abstract */

/*     Compute L_s, the planetocentric longitude of the sun, as seen */
/*     from a specified body. */

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

/*     NAIF_IDS */
/*     PCK */
/*     TIME */
/*     SPK */

/* $ Keywords */

/*     GEOMETRY */
/*     TIME */

/* $ Declarations */
/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     BODY       I   Name of central body. */
/*     ET         I   Epoch in seconds past J2000 TDB. */
/*     ABCORR     I   Aberration correction. */

/*     The function returns the value of L_s for the specified body */
/*     at the specified time. */

/* $ Detailed_Input */

/*     BODY        is the name of the central body, typically a planet. */

/*     ET          is the epoch at which the longitude of the sun (L_s) */
/*                 is to be computed. ET is expressed as seconds past */
/*                 J2000 TDB (Barycentric Dynamical Time). */

/*     ABCORR      indicates the aberration corrections to be applied */
/*                 when computing the longitude of the sun.  ABCORR may */
/*                 be any of the following. */

/*                    'NONE'     Apply no correction. */

/*                    'LT'       Correct the position of the sun, */
/*                               relative to the central body, for */
/*                               planetary (light time) aberration. */

/*                    'LT+S'     Correct the position of the sun, */
/*                               relative to the central body, for */
/*                               planetary and stellar aberrations. */

/* $ Detailed_Output */

/*     The function returns the planetocentric longitude of the sun, */
/*     often called "L_s," for the specified body at the specified time. */
/*     This is the longitude of the body-sun vector in a right-handed */
/*     frame whose basis vectors are defined as follows: */

/*        - The positive Z direction is given by the instantaneous */
/*          angular velocity vector of the orbit of the body about */
/*          the sun. */

/*        - The positive X direction is that of the cross product of the */
/*          instantaneous north spin axis of the body with the positive */
/*          Z direction. */

/*        - The positive Y direction is Z x X. */

/*     Units are radians; the range is 0 to 2*pi.  Longitudes are */
/*     positive to the east. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If the input body name cannot be translated to an ID code, */
/*        and if the name is not a string representation of an integer */
/*        (for example, '399'), the error SPICE(NOTRANSLATION) is */
/*        signaled. */

/*     2) If no SPK (ephemeris) file has been loaded prior to calling */
/*        this routine, or if the SPK data has insufficient coverage, an */
/*        error will be diagnosed and signaled by a routine in the call */
/*        tree of this routine. */

/*     3) If a PCK file containing rotational elements for the central */
/*        body has not been loaded prior to calling this routine, an */
/*        error will be diagnosed and signaled by a routine called by a */
/*        routine in the call tree of this routine. */

/*     4) If the instantaneous angular velocity and spin axis of BODY */
/*        are parallel, the error will be diagnosed and signaled by a */
/*        routine in the call tree of this routine. */

/* $ Files */

/*     1) An SPK file (or file) containing ephemeris data sufficient to */
/*        compute the geometric state of the central body relative to */
/*        the sun at ET must be loaded before this routine is called. If */
/*        light time correction is used, data must be available that */
/*        enable computation of the state the sun relative to the solar */
/*        system barycenter at the light-time corrected epoch.  If */
/*        stellar aberration correction is used, data must be available */
/*        that enable computation of the state the central body relative */
/*        to the solar system barycenter at ET. */

/*     2) A PCK file containing rotational elements for the central body */
/*        must be loaded before this routine is called. */

/* $ Particulars */

/*     The direction of the vernal equinox for the central body is */
/*     determined from the instantaneous equatorial and orbital planes */
/*     of the central body.  This equinox definition is specified in */
/*     reference [1].  The "instantaneous orbital plane" is interpreted */
/*     in this routine as the plane normal to the cross product of the */
/*     position and velocity of the central body relative to the sun. */
/*     The geometric state of the central body relative to the sun is */
/*     used for this normal vector computation. The "instantaneous */
/*     equatorial plane" is normal to the central body's north pole */
/*     at the requested epoch.  The pole direction is determined from */
/*     rotational elements loaded via a PCK file. */

/*     The result returned by this routine will depend on the */
/*     ephemeris data and rotational elements used.  The result may */
/*     differ from that given in any particular version of the */
/*     Astronomical Almanac, due to differences in these input data, */
/*     and due to differences in precision of the computations. */

/* $ Examples */

/*     1) A simple program that computes L_s for a body and time */
/*        supplied interactively.  The geometric state of the sun is */
/*        used. */


/*            PROGRAM EX1 */
/*            IMPLICIT NONE */

/*            DOUBLE PRECISION      DPR */
/*            DOUBLE PRECISION      LSPCN */

/*            CHARACTER*(*)         ABCORR */
/*            PARAMETER           ( ABCORR = 'NONE' ) */

/*            INTEGER               FILSIZ */
/*            PARAMETER           ( FILSIZ = 255 ) */

/*            INTEGER               NAMLEN */
/*            PARAMETER           ( NAMLEN = 36 ) */

/*            INTEGER               TIMLEN */
/*            PARAMETER           ( TIMLEN = 40 ) */

/*            CHARACTER*(NAMLEN)    BODY */
/*            CHARACTER*(FILSIZ)    LSK */
/*            CHARACTER*(FILSIZ)    PCK */
/*            CHARACTER*(FILSIZ)    SPK */
/*            CHARACTER*(TIMLEN)    TIMSTR */

/*            DOUBLE PRECISION      ET */
/*            DOUBLE PRECISION      LON */


/*            CALL PROMPT ( 'Enter name of leapseconds kernel > ', LSK ) */
/*            CALL PROMPT ( 'Enter name of PCK file           > ', PCK ) */
/*            CALL PROMPT ( 'Enter name of SPK file           > ', SPK ) */

/*            CALL FURNSH ( LSK ) */
/*            CALL FURNSH ( PCK ) */
/*            CALL FURNSH ( SPK ) */

/*            WRITE (*,*) ' ' */
/*            WRITE (*,*) 'Kernels have been loaded.' */
/*            WRITE (*,*) ' ' */

/*            DO WHILE ( .TRUE. ) */

/*               CALL PROMPT ( 'Enter name of central body       > ', */
/*           .                  BODY                                  ) */
/*               CALL PROMPT ( 'Enter calendar, JD, or DOY time  > ', */
/*           .                  TIMSTR                                ) */

/*               CALL STR2ET ( TIMSTR, ET ) */

/*      C */
/*      C        Convert longitude to degrees. */
/*      C */
/*               LON = DPR() * LSPCN ( BODY, ET, ABCORR ) */

/*               WRITE (*,*) ' ' */
/*               WRITE (*,*) 'Central body              = ',  BODY */
/*               WRITE (*,*) 'Time                      = ',  TIMSTR */
/*               WRITE (*,*) 'Planetocentric L_s (deg.) = ',  LON */
/*               WRITE (*,*) ' ' */

/*            END DO */

/*            END */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     [1] "The Astronomical Almanac for the Year 2005." U.S. Government */
/*         Printing Office, Washington, D.C., 1984, page L9. */

/* $ Author_and_Institution */

/*     N.J. Bachman       (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 07-JAN-2005 (NJB) */

/* -& */
/* $ Index_Entries */

/*     planetocentric longitude of sun */
/*     compute L_s */
/*     compute Ls */
/*     compute L_sub_s */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Give the function an initial value. */

    ret_val = 0.;

/*     Standard SPICE error handling. */

    if (return_()) {
	return ret_val;
    }
    chkin_("LSPCN", (ftnlen)5);

/*     Map the body name to an ID code. */

    bods2c_(body, &idcode, &found, body_len);
    if (! found) {
	setmsg_("The body name # could not be translated to a NAIF ID code. "
		" The cause of this problem may be that you need an updated v"
		"ersion of the SPICE Toolkit.", (ftnlen)147);
	errch_("#", body, (ftnlen)1, body_len);
	sigerr_("SPICE(NOTRANSLATION)", (ftnlen)20);
	chkout_("LSPCN", (ftnlen)5);
	return ret_val;
    }

/*     Look up the direction of the North pole of the central body. */
/*     Note that TIPBOD does make use of binary PCK data if available. */

    tipbod_("J2000", &idcode, et, tipm, (ftnlen)5);
    for (i__ = 1; i__ <= 3; ++i__) {
	npole[(i__1 = i__ - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("npole", i__1,
		 "lspcn_", (ftnlen)339)] = tipm[(i__2 = i__ * 3 - 1) < 9 && 0 
		<= i__2 ? i__2 : s_rnge("tipm", i__2, "lspcn_", (ftnlen)339)];
    }

/*     Get the geometric state of the body relative to the sun. */

    spkgeo_(&idcode, et, "J2000", &c__10, bstate, &lt, (ftnlen)5);

/*     Get the unit direction vector parallel to the angular velocity */
/*     vector of the orbit.  This is just the unitized cross product of */
/*     position and velocity. */

    ucrss_(bstate, &bstate[3], uavel);

/*     We want to create a transformation matrix that maps vectors from */
/*     basis REF to the following frame: */
/*        Z  =  UAVEL */

/*        X  =  NPOLE x UAVEL */

/*        Y  =  Z x X */

/*     This is a "two-vector" frame with the unit orbital */
/*     angular velocity vector UAVEL as the primary vector and the */
/*     spin axis NPOLE as the secondary vector.  The primary */
/*     vector is associated with the +Z axis; the secondary vector */
/*     is associated with the +Y axis. */

    twovec_(uavel, &c__3, npole, &c__2, trans);
    if (failed_()) {
	chkout_("LSPCN", (ftnlen)5);
	return ret_val;
    }

/*     We'll find the position of the Sun relative to this frame. */

/*     Get the state of the sun in frame REF.  Since we may be using */
/*     aberration corrections, this is not necessarily the negative of */
/*     the state we've just found. */

    spkezr_("SUN", et, "J2000", abcorr, body, sstate, &lt, (ftnlen)3, (ftnlen)
	    5, abcorr_len, body_len);

/*     Now transform the position of the Sun into the "orbit plane */
/*     and equinox" frame. */

    mxv_(trans, sstate, pos);

/*     Let RECRAD find the longitude LS for us.  RECRAD performs */
/*     the same coordinate transformation as the more commonly used */
/*     RECLAT, but the range of right ascension is 0:2*pi, which is */
/*     what we want for Ls. */

    recrad_(pos, &radius, &ret_val, &lat);
    chkout_("LSPCN", (ftnlen)5);
    return ret_val;
} /* lspcn_ */

