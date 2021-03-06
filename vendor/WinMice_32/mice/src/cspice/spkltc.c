/* spkltc.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static integer c__6 = 6;
static doublereal c_b25 = -1.;

/* $Procedure SPKLTC ( S/P Kernel, light time corrected state ) */
/* Subroutine */ int spkltc_(integer *targ, doublereal *et, char *ref, char *
	abcorr, doublereal *stobs, doublereal *starg, doublereal *lt, 
	doublereal *dlt, ftnlen ref_len, ftnlen abcorr_len)
{
    /* Initialized data */

    static logical first = TRUE_;
    static char prvcor[5] = "     ";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal dist;
    extern doublereal vdot_(doublereal *, doublereal *);
    static logical xmit;
    doublereal a, b, c__;
    integer i__;
    extern /* Subroutine */ int zzprscor_(char *, logical *, ftnlen);
    integer refid;
    extern /* Subroutine */ int chkin_(char *, ftnlen), errch_(char *, char *,
	     ftnlen, ftnlen);
    static logical usecn;
    extern /* Subroutine */ int vlcom_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), vsubg_(doublereal *, doublereal *,
	     integer *, doublereal *);
    doublereal ssblt;
    static logical uselt;
    extern doublereal vnorm_(doublereal *);
    extern logical failed_(void);
    extern doublereal clight_(void);
    logical attblk[15];
    extern /* Subroutine */ int spkgeo_(integer *, doublereal *, char *, 
	    integer *, doublereal *, doublereal *, ftnlen), sigerr_(char *, 
	    ftnlen), chkout_(char *, ftnlen);
    integer ltsign;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), irfnum_(char *, 
	    integer *, ftnlen);
    doublereal ssbtrg[6];
    integer numitr;
    extern logical return_(void);
    logical usestl;

/* $ Abstract */

/*     Return the state (position and velocity) of a target body */
/*     relative to an observer, optionally corrected for light time, */
/*     expressed relative to an inertial reference frame. */

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

/*     SPK */

/* $ Keywords */

/*     EPHEMERIS */

/* $ Declarations */
/* $ Abstract */

/*     Include file zzabcorr.inc */

/*     SPICE private file intended solely for the support of SPICE */
/*     routines.  Users should not include this file directly due */
/*     to the volatile nature of this file */

/*     The parameters below define the structure of an aberration */
/*     correction attribute block. */

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

/* $ Parameters */

/*     An aberration correction attribute block is an array of logical */
/*     flags indicating the attributes of the aberration correction */
/*     specified by an aberration correction string.  The attributes */
/*     are: */

/*        - Is the correction "geometric"? */

/*        - Is light time correction indicated? */

/*        - Is stellar aberration correction indicated? */

/*        - Is the light time correction of the "converged */
/*          Newtonian" variety? */

/*        - Is the correction for the transmission case? */

/*        - Is the correction relativistic? */

/*    The parameters defining the structure of the block are as */
/*    follows: */

/*       NABCOR    Number of aberration correction choices. */

/*       ABATSZ    Number of elements in the aberration correction */
/*                 block. */

/*       GEOIDX    Index in block of geometric correction flag. */

/*       LTIDX     Index of light time flag. */

/*       STLIDX    Index of stellar aberration flag. */

/*       CNVIDX    Index of converged Newtonian flag. */

/*       XMTIDX    Index of transmission flag. */

/*       RELIDX    Index of relativistic flag. */

/*    The following parameter is not required to define the block */
/*    structure, but it is convenient to include it here: */

/*       CORLEN    The maximum string length required by any aberration */
/*                 correction string */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 18-DEC-2004 (NJB) */

/* -& */
/*     Number of aberration correction choices: */


/*     Aberration correction attribute block size */
/*     (number of aberration correction attributes): */


/*     Indices of attributes within an aberration correction */
/*     attribute block: */


/*     Maximum length of an aberration correction string: */


/*     End of include file zzabcorr.inc */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     TARG       I   Target body. */
/*     ET         I   Observer epoch. */
/*     REF        I   Inertial reference frame of output state. */
/*     ABCORR     I   Aberration correction flag. */
/*     STOBS      I   State of the observer relative to the SSB. */
/*     STARG      O   State of target. */
/*     LT         O   One way light time between observer and target. */
/*     DLT        O   Derivative of light time with respect to time. */

/* $ Detailed_Input */

/*     TARG        is the NAIF ID code for a target body.  The target */
/*                 and observer define a state vector whose position */
/*                 component points from the observer to the target. */

/*     ET          is the ephemeris time, expressed as seconds past */
/*                 J2000 TDB, at which the state of the target body */
/*                 relative to the observer is to be computed. ET */
/*                 refers to time at the observer's location. */

/*     REF         is the inertial reference frame with respect to which */
/*                 the input state STOBS and the output state STARG are */
/*                 expressed. REF must be recognized by the SPICE */
/*                 Toolkit. The acceptable frames are listed in the */
/*                 Frames Required Reading, as well as in the SPICELIB */
/*                 routine CHGIRF. */

/*                 Case and blanks are not significant in the string */
/*                 REF. */


/*     ABCORR      indicates the aberration corrections to be applied to */
/*                 the state of the target body to account for one-way */
/*                 light time. See the discussion in the Particulars */
/*                 section for recommendations on how to choose */
/*                 aberration corrections. */

/*                 If ABCORR includes the stellar aberration correction */
/*                 symbol '+S', this flag is simply ignored. Aside from */
/*                 the possible presence of this symbol, ABCORR may be */
/*                 any of the following: */

/*                    'NONE'     Apply no correction. Return the */
/*                               geometric state of the target body */
/*                               relative to the observer. */

/*                 The following values of ABCORR apply to the */
/*                 "reception" case in which photons depart from the */
/*                 target's location at the light-time corrected epoch */
/*                 ET-LT and *arrive* at the observer's location at ET: */

/*                    'LT'       Correct for one-way light time (also */
/*                               called "planetary aberration") using a */
/*                               Newtonian formulation. This correction */
/*                               yields the state of the target at the */
/*                               moment it emitted photons arriving at */
/*                               the observer at ET. */

/*                               The light time correction involves */
/*                               iterative solution of the light time */
/*                               equation (see Particulars for details). */
/*                               The solution invoked by the 'LT' option */
/*                               uses one iteration. */

/*                    'CN'       Converged Newtonian light time */
/*                               correction.  In solving the light time */
/*                               equation, the 'CN' correction iterates */
/*                               until the solution converges (three */
/*                               iterations on all supported platforms). */

/*                               The 'CN' correction typically does not */
/*                               substantially improve accuracy because */
/*                               the errors made by ignoring */
/*                               relativistic effects may be larger than */
/*                               the improvement afforded by obtaining */
/*                               convergence of the light time solution. */
/*                               The 'CN' correction computation also */
/*                               requires a significantly greater number */
/*                               of CPU cycles than does the */
/*                               one-iteration light time correction. */

/*                 The following values of ABCORR apply to the */
/*                 "transmission" case in which photons *depart* from */
/*                 the observer's location at ET and arrive at the */
/*                 target's location at the light-time corrected epoch */
/*                 ET+LT: */

/*                    'XLT'      "Transmission" case:  correct for */
/*                               one-way light time using a Newtonian */
/*                               formulation. This correction yields the */
/*                               state of the target at the moment it */
/*                               receives photons emitted from the */
/*                               observer's location at ET. */

/*                    'XCN'      "Transmission" case:  converged */
/*                               Newtonian light time correction. */


/*                 Neither special nor general relativistic effects are */
/*                 accounted for in the aberration corrections applied */
/*                 by this routine. */

/*                 Case and blanks are not significant in the string */
/*                 ABCORR. */


/*     STOBS       is the geometric (uncorrected) state of the observer */
/*                 relative to the solar system barycenter at epoch ET. */
/*                 STOBS is a 6-vector: the first three components of */
/*                 STOBS represent a Cartesian position vector; the last */
/*                 three components represent the corresponding velocity */
/*                 vector. STOBS is expressed relative to the inertial */
/*                 reference frame designated by REF. */

/*                 Units are always km and km/sec. */

/* $ Detailed_Output */

/*     STARG       is a Cartesian state vector representing the position */
/*                 and velocity of the target body relative to the */
/*                 specified observer. STARG is corrected for the */
/*                 specified aberration, and is expressed with respect */
/*                 to the specified inertial reference frame.  The first */
/*                 three components of STARG represent the x-, y- and */
/*                 z-components of the target's position; last three */
/*                 components form the corresponding velocity vector. */

/*                 The position component of STARG points from the */
/*                 observer's location at ET to the aberration-corrected */
/*                 location of the target. Note that the sense of the */
/*                 position vector is independent of the direction of */
/*                 radiation travel implied by the aberration */
/*                 correction. */

/*                 Units are always km and km/sec. */

/*     LT          is the one-way light time between the observer and */
/*                 target in seconds.  If the target state is corrected */
/*                 for light time, then LT is the one-way light time */
/*                 between the observer and the light time-corrected */
/*                 target location. */

/*     DLT         is the derivative with respect to barycentric */
/*                 dynamical time of the one way light time between */
/*                 target and observer: */

/*                    DLT = d(LT)/d(ET) */

/*                 DLT can also be described as the rate of change of */
/*                 one way light time. DLT is unitless, since LT and */
/*                 ET both have units of TDB seconds. */

/*                 If the observer and target are at the same position, */
/*                 then DLT is set to zero. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) For the convenience of the caller, the input aberration */
/*        correction flag can call for stellar aberration correction via */
/*        inclusion of the '+S' suffix. This portion of the aberration */
/*        correction flag is ignored if present. */

/*     2) If ABCORR calls for stellar aberration but not light */
/*        time corrections, the error SPICE(NOTSUPPORTED) is */
/*        signaled. */

/*     3) If ABCORR calls for relativistic light time corrections, the */
/*        error SPICE(NOTSUPPORTED) is signaled. */

/*     4) If the value of ABCORR is not recognized, the error */
/*        is diagnosed by a routine in the call tree of this */
/*        routine. */

/*     5) If the reference frame requested is not a recognized */
/*        inertial reference frame, the error SPICE(BADFRAME) */
/*        is signaled. */

/*     6) If the state of the target relative to the solar system */
/*        barycenter cannot be computed, the error will be diagnosed */
/*        by routines in the call tree of this routine. */

/*     7) If the observer and target are at the same position, */
/*        then DLT is set to zero. This situation could arise, */
/*        for example, when the observer is Mars and the target */
/*        is the Mars barycenter. */

/*     8) If a division by zero error would occur in the computation */
/*        of DLT, the error SPICE(DIVIDEBYZERO) is signaled. */

/* $ Files */

/*     This routine computes states using SPK files that have been */
/*     loaded into the SPICE system, normally via the kernel loading */
/*     interface routine FURNSH.  Application programs typically load */
/*     kernels once before this routine is called, for example during */
/*     program initialization; kernels need not be loaded repeatedly. */
/*     See the routine FURNSH and the SPK and KERNEL Required Reading */
/*     for further information on loading (and unloading) kernels. */

/*     If any of the ephemeris data used to compute STARG are expressed */
/*     relative to a non-inertial frame in the SPK files providing those */
/*     data, additional kernels may be needed to enable the reference */
/*     frame transformations required to compute the state. Normally */
/*     these additional kernels are PCK files or frame kernels. Any */
/*     such kernels must already be loaded at the time this routine is */
/*     called. */

/* $ Particulars */

/*     This routine supports higher-level SPK API routines that can */
/*     perform both light time and stellar aberration corrections. */
/*     User applications normally will not need to call this routine */
/*     directly. */

/*     See the header of the routine SPKEZR for a detailed discussion */
/*     of aberration corrections. */

/* $ Examples */


/*    1) Look up a sequence of states of the Moon as seen from the */
/*       Earth. Use light time corrections. Compute the first state for */
/*       the epoch 2000 JAN 1 12:00:00 TDB; compute subsequent states at */
/*       intervals of 1 hour. For each epoch, display the states, the */
/*       one way light time between target and observer, and the rate of */
/*       change of the one way light time. */

/*       Use the following meta-kernel to specify the kernels to */
/*       load: */

/*          KPL/MK */

/*          This meta-kernel is intended to support operation of SPICE */
/*          example programs. The kernels shown here should not be */
/*          assumed to contain adequate or correct versions of data */
/*          required by SPICE-based user applications. */

/*          In order for an application to use this meta-kernel, the */
/*          kernels referenced here must be present in the user's */
/*          current working directory. */


/*          \begindata */

/*             KERNELS_TO_LOAD = ( 'de418.bsp', */
/*                                 'pck00008.tpc', */
/*                                 'naif0008.tls'  ) */

/*          \begintext */


/*       The code example follows: */

/*           PROGRAM EX1 */
/*           IMPLICIT NONE */
/*     C */
/*     C     Local constants */
/*     C */
/*     C     The meta-kernel name shown here refers to a file whose */
/*     C     contents are those shown above. This file and the kernels */
/*     C     it references must exist in your current working directory. */
/*     C */
/*           CHARACTER*(*)         META */
/*           PARAMETER           ( META   = 'example.mk' ) */
/*     C */
/*     C     Use a time step of 1 hour; look up 5 states. */
/*     C */
/*           DOUBLE PRECISION      STEP */
/*           PARAMETER           ( STEP   = 3600.0D0 ) */

/*           INTEGER               MAXITR */
/*           PARAMETER           ( MAXITR = 5 ) */
/*     C */
/*     C     Local variables */
/*     C */
/*           DOUBLE PRECISION      DLT */
/*           DOUBLE PRECISION      ET */
/*           DOUBLE PRECISION      ET0 */
/*           DOUBLE PRECISION      LT */
/*           DOUBLE PRECISION      STATE ( 6 ) */
/*           DOUBLE PRECISION      STOBS ( 6 ) */
/*           INTEGER               I */

/*     C */
/*     C     Load the SPK and LSK kernels via the meta-kernel. */
/*     C */
/*           CALL FURNSH ( META ) */
/*     C */
/*     C     Convert the start time to seconds past J2000 TDB. */
/*     C */
/*           CALL STR2ET ( '2000 JAN 1 12:00:00 TDB', ET0 ) */
/*     C */
/*     C     Step through a series of epochs, looking up a */
/*     C     state vector at each one. */
/*     C */
/*           DO I = 1, MAXITR */

/*              ET = ET0 + (I-1)*STEP */

/*     C */
/*     C        Look up a state vector at epoch ET using the */
/*     C        following inputs: */
/*     C */
/*     C           Target:                 Moon (NAIF ID code 301) */
/*     C           Reference frame:        J2000 */
/*     C           Aberration correction:  Light time ('LT') */
/*     C           Observer:               Earth (NAIF ID code 399) */
/*     C */
/*     C        Before we can execute this computation, we'll need the */
/*     C        geometric state of the observer relative to the solar */
/*     C        system barycenter at ET, expressed relative to the */
/*     C        J2000 reference frame: */
/*     C */
/*              CALL SPKSSB ( 399, ET,    'J2000', STOBS ) */
/*     C */
/*     C        Now compute the desired state vector: */
/*     C */
/*              CALL SPKLTC ( 301,   ET,    'J2000', 'LT', */
/*          .                 STOBS, STATE, LT,      DLT     ) */

/*              WRITE (*,*) 'ET = ', ET */
/*              WRITE (*,*) 'J2000 x-position (km):   ', STATE(1) */
/*              WRITE (*,*) 'J2000 y-position (km):   ', STATE(2) */
/*              WRITE (*,*) 'J2000 z-position (km):   ', STATE(3) */
/*              WRITE (*,*) 'J2000 x-velocity (km/s): ', STATE(4) */
/*              WRITE (*,*) 'J2000 y-velocity (km/s): ', STATE(5) */
/*              WRITE (*,*) 'J2000 z-velocity (km/s): ', STATE(6) */
/*              WRITE (*,*) 'One-way light time (s):  ', LT */
/*              WRITE (*,*) 'Light time rate:         ', DLT */
/*              WRITE (*,*) ' ' */

/*           END DO */

/*           END */


/*     The output produced by this program will vary somewhat as */
/*     a function of the platform on which the program is built and */
/*     executed. On a PC/Linux/g77 platform, the following output */
/*     was produced: */

/*        ET =   0. */
/*        J2000 x-position (km):    -291569.265 */
/*        J2000 y-position (km):    -266709.186 */
/*        J2000 z-position (km):    -76099.1551 */
/*        J2000 x-velocity (km/s):   0.643530613 */
/*        J2000 y-velocity (km/s):  -0.666081817 */
/*        J2000 z-velocity (km/s):  -0.301322832 */
/*        One-way light time (s):    1.34231061 */
/*        Light time rate:           1.07316909E-07 */

/*        ET =   3600. */
/*        J2000 x-position (km):    -289240.781 */
/*        J2000 y-position (km):    -269096.441 */
/*        J2000 z-position (km):    -77180.8997 */
/*        J2000 x-velocity (km/s):   0.650062115 */
/*        J2000 y-velocity (km/s):  -0.660162739 */
/*        J2000 z-velocity (km/s):  -0.299642674 */
/*        One-way light time (s):    1.34269395 */
/*        Light time rate:           1.05652599E-07 */

/*        ET =   7200. */
/*        J2000 x-position (km):    -286888.887 */
/*        J2000 y-position (km):    -271462.302 */
/*        J2000 z-position (km):    -78256.5557 */
/*        J2000 x-velocity (km/s):   0.656535992 */
/*        J2000 y-velocity (km/s):  -0.654196577 */
/*        J2000 z-velocity (km/s):  -0.297940273 */
/*        One-way light time (s):    1.34307131 */
/*        Light time rate:           1.03990457E-07 */

/*        ET =   10800. */
/*        J2000 x-position (km):    -284513.792 */
/*        J2000 y-position (km):    -273806.6 */
/*        J2000 z-position (km):    -79326.0432 */
/*        J2000 x-velocity (km/s):   0.662951901 */
/*        J2000 y-velocity (km/s):  -0.648183807 */
/*        J2000 z-velocity (km/s):  -0.296215779 */
/*        One-way light time (s):    1.34344269 */
/*        Light time rate:           1.02330665E-07 */

/*        ET =   14400. */
/*        J2000 x-position (km):    -282115.704 */
/*        J2000 y-position (km):    -276129.17 */
/*        J2000 z-position (km):    -80389.283 */
/*        J2000 x-velocity (km/s):   0.669309504 */
/*        J2000 y-velocity (km/s):  -0.642124908 */
/*        J2000 z-velocity (km/s):  -0.294469343 */
/*        One-way light time (s):    1.3438081 */
/*        Light time rate:           1.00673404E-07 */


/* $ Restrictions */

/*     1) The routine SPKGEO should be used instead of this routine */
/*        to compute geometric states. SPKGEO introduces less */
/*        round-off error when the observer and target have common */
/*        center that is closer to both objects than is the solar */
/*        system barycenter. */

/*     2) The kernel files to be used by SPKLTC must be loaded */
/*        (normally by the SPICELIB kernel loader FURNSH) before */
/*        this routine is called. */

/*     3) Unlike most other SPK state computation routines, this */
/*        routine requires that the output state be relative to an */
/*        inertial reference frame. */

/* $ Literature_References */

/*     SPK Required Reading. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 11-JAN-2008 (NJB) */

/* -& */
/* $ Index_Entries */

/*     low-level light time correction */
/*     light-time corrected state from spk file */
/*     get light-time corrected state */

/* -& */
/* $ Revisions */

/*     None. */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     TOL is the tolerance used for a division-by-zero test */
/*     performed prior to computation of DLT. */


/*     Local variables */


/*     Saved variables */


/*     Initial values */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("SPKLTC", (ftnlen)6);
    }
    if (first || s_cmp(abcorr, prvcor, abcorr_len, (ftnlen)5) != 0) {

/*        The aberration correction flag differs from the value it */
/*        had on the previous call, if any.  Analyze the new flag. */

	zzprscor_(abcorr, attblk, abcorr_len);
	if (failed_()) {
	    chkout_("SPKLTC", (ftnlen)6);
	    return 0;
	}

/*        The aberration correction flag is recognized; save it. */

	s_copy(prvcor, abcorr, (ftnlen)5, abcorr_len);

/*        Set logical flags indicating the attributes of the requested */
/*        correction: */

/*           XMIT is .TRUE. when the correction is for transmitted */
/*           radiation. */

/*           USELT is .TRUE. when any type of light time correction */
/*           (normal or converged Newtonian) is specified. */

/*           USECN indicates converged Newtonian light time correction. */

/*        The above definitions are consistent with those used by */
/*        ZZPRSCOR. */

	xmit = attblk[4];
	uselt = attblk[1];
	usecn = attblk[3];
	usestl = attblk[2];
	if (usestl && ! uselt) {
	    setmsg_("Aberration correction flag # calls for stellar aberrati"
		    "on but not light time corrections. This combination is n"
		    "ot expected.", (ftnlen)123);
	    errch_("#", abcorr, (ftnlen)1, abcorr_len);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    chkout_("SPKLTC", (ftnlen)6);
	    return 0;
	} else if (attblk[5]) {
	    setmsg_("Aberration correction flag # calls for relativistic lig"
		    "ht time correction.", (ftnlen)74);
	    errch_("#", abcorr, (ftnlen)1, abcorr_len);
	    sigerr_("SPICE(NOTSUPPORTED)", (ftnlen)19);
	    chkout_("SPKLTC", (ftnlen)6);
	    return 0;
	}
	first = FALSE_;
    }

/*     See if the reference frame is a recognized inertial frame. */

    irfnum_(ref, &refid, ref_len);
    if (refid == 0) {
	setmsg_("The requested frame '#' is not a recognized inertial frame. "
		, (ftnlen)60);
	errch_("#", ref, (ftnlen)1, ref_len);
	sigerr_("SPICE(BADFRAME)", (ftnlen)15);
	chkout_("SPKLTC", (ftnlen)6);
	return 0;
    }

/*     Find the geometric state of the target body with respect to */
/*     the solar system barycenter. Subtract the state of the */
/*     observer to get the relative state. Use this to compute the */
/*     one-way light time. */

    spkgeo_(targ, et, ref, &c__0, ssbtrg, &ssblt, ref_len);
    vsubg_(ssbtrg, stobs, &c__6, starg);
    dist = vnorm_(starg);
    *lt = dist / clight_();
    if (*lt == 0.) {

/*        This can happen only if the observer and target are at the */
/*        same position. We don't consider this an error, but we're not */
/*        going to compute the light time derivative. */

	*dlt = 0.;
	chkout_("SPKLTC", (ftnlen)6);
	return 0;
    }
    if (! uselt) {

/*        This is a special case: we're not using light time */
/*        corrections, so the derivative */
/*        of light time is just */

/*           (1/c) * d(VNORM(STARG))/dt */

	*dlt = vdot_(starg, &starg[3]) / (dist * clight_());

/*        LT and DLT are both set, so we can return. */

	chkout_("SPKLTC", (ftnlen)6);
	return 0;
    }

/*     To correct for light time, find the state of the target body */
/*     at the current epoch minus the one-way light time. Note that */
/*     the observer remains where it is. */

/*     Determine the sign of the light time offset. */

    if (xmit) {
	ltsign = 1;
    } else {
	ltsign = -1;
    }

/*     Let NUMITR be the number of iterations we'll perform to */
/*     compute the light time. */

    if (usecn) {
	numitr = 3;
    } else {
	numitr = 1;
    }
    i__1 = numitr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = *et + ltsign * *lt;
	spkgeo_(targ, &d__1, ref, &c__0, ssbtrg, &ssblt, ref_len);
	vsubg_(ssbtrg, stobs, &c__6, starg);
	*lt = vnorm_(starg) / clight_();
    }

/*     At this point, STARG contains the light time corrected */
/*     state of the target relative to the observer. */

/*     Compute the derivative of light time with respect */
/*     to time: dLT/dt.  Below we derive the formula for */
/*     this quantity for the reception case. Let */

/*        POBS be the position of the observer relative to the */
/*        solar system barycenter. */

/*        VOBS be the velocity of the observer relative to the */
/*        solar system barycenter. */

/*        PTARG be the position of the target relative to the */
/*        solar system barycenter. */

/*        VTARG be the velocity of the target relative to the */
/*        solar system barycenter. */

/*        S be the sign of the light time correction. S is */
/*        negative for the reception case. */

/*     The light-time corrected position of the target relative to */
/*     the observer at observation time ET, given the one-way */
/*     light time LT is: */

/*         PTARG(ET+S*LT) - POBS(ET) */

/*     The light-time corrected velocity of the target relative to */
/*     the observer at observation time ET is */

/*         VTARG(ET+S*LT)*( 1 + S*d(LT)/d(ET) ) - VOBS(ET) */

/*     We need to compute dLT/dt. Below, we use the facts that, */
/*     for a time-dependent vector X(t), */

/*          ||X||     = <X,X> ** (1/2) */

/*        d(||X||)/dt = (1/2)<X,X>**(-1/2) * 2 * <X,dX/dt> */

/*                    = <X,X>**(-1/2) *  <X,dX/dt> */

/*                    = <X,dX/dt> / ||X|| */

/*     Newtonian light time equation: */

/*        LT     =   (1/c) * || PTARG(ET+S*LT) - POBS(ET)|| */

/*     Differentiate both sides: */

/*        dLT/dt =   (1/c) * ( 1 / || PTARG(ET+S*LT) - POBS(ET) || ) */

/*                  * < PTARG(ET+S*LT) - POBS(ET), */
/*                      VTARG(ET+S*LT)*(1+S*d(LT)/d(ET)) - VOBS(ET) > */


/*               = (1/c) * ( 1 / || PTARG(ET+S*LT) - POBS(ET) || ) */

/*                 * (  < PTARG(ET+S*LT) - POBS(ET), */
/*                        VTARG(ET+S*LT) - VOBS(ET) > */

/*                   +  < PTARG(ET+S*LT) - POBS(ET), */
/*                        VTARG(ET+S*LT)           > * (S*d(LT)/d(ET))  ) */

/*     Let */

/*        A =   (1/c) * ( 1 / || PTARG(ET+S*LT) - POBS(ET) || ) */

/*        B =   < PTARG(ET+S*LT) - POBS(ET), VTARG(ET+S*LT) - VOBS(ET) > */

/*        C =   < PTARG(ET+S*LT) - POBS(ET), VTARG(ET+S*LT) > */

/*     Then */

/*        d(LT)/d(ET) =  A * ( B  +  C * S*d(LT)/d(ET) ) */

/*     which implies */

/*        d(LT)/d(ET) =  A*B / ( 1 - S*C*A ) */



    a = 1. / (clight_() * vnorm_(starg));
    b = vdot_(starg, &starg[3]);
    c__ = vdot_(starg, &ssbtrg[3]);

/*     For physically realistic target velocities, S*C*A cannot equal 1. */
/*     We'll check for this case anyway. */

    if (ltsign * c__ * a > .99999999989999999) {
	setmsg_("Target range rate magnitude is approximately the speed of l"
		"ight. The light time derivative cannot be computed.", (ftnlen)
		110);
	sigerr_("SPICE(DIVIDEBYZERO)", (ftnlen)19);
	chkout_("SPKLTC", (ftnlen)6);
	return 0;
    }

/*     Compute DLT: the rate of change of light time. */

    *dlt = a * b / (1. - ltsign * c__ * a);

/*     Overwrite the velocity portion of the output state */
/*     with the light-time corrected velocity. */

    d__1 = ltsign * *dlt + 1.;
    vlcom_(&d__1, &ssbtrg[3], &c_b25, &stobs[3], &starg[3]);
    chkout_("SPKLTC", (ftnlen)6);
    return 0;
} /* spkltc_ */

