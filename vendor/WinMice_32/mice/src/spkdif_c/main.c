/* spkdiff.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__15 = 15;
static doublereal c_b42 = 3600.;
static doublereal c_b43 = 60.;
static doublereal c_b44 = 1.;
static integer c__0 = 0;

/* $Procedure   SPKDIFF ( Compare two SPK files. ) */
/* Main program */ MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_dnnt(doublereal *), s_cmp(char *, char *, ftnlen, ftnlen), 
	    s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static char line[1024];
    static integer msec;
    static char time[1024*2];
    static doublereal step;
    static integer nitr, hour;
    extern /* Subroutine */ int zzgeoseg_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *);
    static integer i__, j, n, bodid[2], cenid[2];
    extern /* Subroutine */ int etcal_(doublereal *, char *, ftnlen);
    static char frame[32*2];
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    static doublereal epoch[1000000];
    extern /* Subroutine */ int errch_(char *, char *, ftnlen, ftnlen), 
	    repmc_(char *, char *, char *, char *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), repmf_(char *, char *, doublereal *, integer *, char *, 
	    char *, ftnlen, ftnlen, ftnlen, ftnlen), repmi_(char *, char *, 
	    integer *, char *, ftnlen, ftnlen, ftnlen), errdp_(char *, 
	    doublereal *, ftnlen);
    static char hword[32];
    extern integer rtrim_(char *, ftnlen);
    static doublereal state1[6000000]	/* was [6][1000000] */, state2[
	    6000000]	/* was [6][1000000] */, et[2];
    static integer handle;
    static char bodnam[32*2], cennam[32*2];
    static doublereal lt;
    extern /* Subroutine */ int getcml_(char *, ftnlen), rmaind_(doublereal *,
	     doublereal *, doublereal *, doublereal *), chwcml_(char *, char *
	    , char *, integer *, char *, integer *, char *, char *, 
	    doublereal *, integer *, doublereal *, char *, char *, char *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static char infmsg[1024*21];
    extern /* Subroutine */ int spklef_(char *, integer *, ftnlen), stdiff_(
	    doublereal *, doublereal *, integer *, doublereal *, char *, char 
	    *, ftnlen, ftnlen), spkgeo_(integer *, doublereal *, char *, 
	    integer *, doublereal *, doublereal *, ftnlen);
    static integer hanlst[200];
    extern /* Subroutine */ int sigerr_(char *, ftnlen);
    static doublereal dsclst[1000]	/* was [5][200] */;
    static char kernls[1024], timfmt[1024], diftyp[32];
    static integer spkcnt, minuts;
    extern /* Subroutine */ int errprt_(char *, char *, ftnlen, ftnlen), 
	    intstr_(integer *, char *, ftnlen), prefix_(char *, integer *, 
	    char *, ftnlen, ftnlen), tostdo_(char *, ftnlen), ktotal_(char *, 
	    integer *, ftnlen), setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), spkuef_(integer *), chkout_(char *, ftnlen);
    static integer sec, day;
    extern doublereal spd_(void);
    static char spk[1024*2];
    static doublereal hdp1, hdp2, hdp3, hdp4;

/* $ Abstract */

/*     SPKDIFF is a program that finds differences between geometric */
/*     states computed from two SPK files and either displays these */
/*     differences or shows statistics about them. */

/*     For complete information about the program see SPKDIFF User's */
/*     Guide. */

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

/* $ Author_and_Institution */

/*     B.V. Semenov   (JPL) */

/* $ Version */

/* -    Version 1.0.0, 18-JUL-2006 (BVS). */

/* -& */

/*     Global parameters */


/*     Local variables */

/* $ Abstract */

/*     Include Section:  SPKDIFF Global Parameters */

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

/* $ Author_and_Institution */

/*     B.V. Semenov   (JPL) */

/* $ Version */

/* -    Version 1.0.0, 20-JUL-2006 (BVS). */

/* -& */

/*     Program name and version. */


/*     Command line keys. */


/*     Max and min number states that the program can handle. */


/*     Default number states. */


/*     Line size parameters. */


/*     Version/usage display parameters. */


/*     Maximum segment buffer size for Nat's ZZGEOSEG. */


/*     DAF descriptor size and component counts. */


/*     Cell lower boundary. */


/*     Maximum allowed number of coverage windows. */


/*     Smallest allowed step. */


/*     SPICELIB functions. */


/*     Save everything to prevent potential memory problems in f2c'ed */
/*     version. */


/*     Check in. */

    chkin_("spkdiff", (ftnlen)7);

/*     Reset default error messages. */

    errprt_("SET", "NONE, SHORT, LONG, TRACEBACK", (ftnlen)3, (ftnlen)28);

/*     Get command line and call the "big kahuna" deal-with-command-line */
/*     routine and get back all needed setups. */

    getcml_(line, (ftnlen)1024);
    chwcml_(line, spk, bodnam, bodid, cennam, cenid, frame, time, et, &nitr, &
	    step, diftyp, timfmt, kernls, (ftnlen)1024, (ftnlen)1024, (ftnlen)
	    32, (ftnlen)32, (ftnlen)32, (ftnlen)1024, (ftnlen)32, (ftnlen)
	    1024, (ftnlen)1024);

/*     OK. We are done with inputs and report about what we are going to */
/*     compare. */

    s_copy(infmsg, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 1024, "# Comparison of $ '$'-referenced geometric states",
	     (ftnlen)1024, (ftnlen)49);
    s_copy(infmsg + 2048, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 3072, "#    of '$' ($) relative to '$' ($)", (ftnlen)1024,
	     (ftnlen)35);
    s_copy(infmsg + 4096, "#    from SPK '$'", (ftnlen)1024, (ftnlen)17);
    s_copy(infmsg + 5120, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 6144, "# with $ '$'-referenced geometric states", (ftnlen)
	    1024, (ftnlen)40);
    s_copy(infmsg + 7168, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 8192, "#    of '$' ($) relative to '$' ($)", (ftnlen)1024,
	     (ftnlen)35);
    s_copy(infmsg + 9216, "#    from SPK '$'", (ftnlen)1024, (ftnlen)17);
    s_copy(infmsg + 10240, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 11264, "# evenly-spaced with $ second ($) step size", (
	    ftnlen)1024, (ftnlen)43);
    s_copy(infmsg + 12288, "# within the time interval", (ftnlen)1024, (
	    ftnlen)26);
    s_copy(infmsg + 13312, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 14336, "#    from '$' ($ TDB seconds)", (ftnlen)1024, (
	    ftnlen)29);
    s_copy(infmsg + 15360, "#    to   '$' ($ TDB seconds)", (ftnlen)1024, (
	    ftnlen)29);
    s_copy(infmsg + 16384, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 17408, "# using additional data from these kernels", (
	    ftnlen)1024, (ftnlen)42);
    s_copy(infmsg + 18432, "# ", (ftnlen)1024, (ftnlen)2);
    s_copy(infmsg + 19456, "#    '$'", (ftnlen)1024, (ftnlen)8);
    s_copy(infmsg + 20480, "# ", (ftnlen)1024, (ftnlen)2);
    repmi_(infmsg + 1024, "$", &nitr, infmsg + 1024, (ftnlen)1024, (ftnlen)1, 
	    (ftnlen)1024);
    repmc_(infmsg + 1024, "$", frame, infmsg + 1024, (ftnlen)1024, (ftnlen)1, 
	    (ftnlen)32, (ftnlen)1024);
    repmc_(infmsg + 3072, "$", bodnam, infmsg + 3072, (ftnlen)1024, (ftnlen)1,
	     (ftnlen)32, (ftnlen)1024);
    repmi_(infmsg + 3072, "$", bodid, infmsg + 3072, (ftnlen)1024, (ftnlen)1, 
	    (ftnlen)1024);
    repmc_(infmsg + 3072, "$", cennam, infmsg + 3072, (ftnlen)1024, (ftnlen)1,
	     (ftnlen)32, (ftnlen)1024);
    repmi_(infmsg + 3072, "$", cenid, infmsg + 3072, (ftnlen)1024, (ftnlen)1, 
	    (ftnlen)1024);
    repmc_(infmsg + 4096, "$", spk, infmsg + 4096, (ftnlen)1024, (ftnlen)1, (
	    ftnlen)1024, (ftnlen)1024);
    repmi_(infmsg + 6144, "$", &nitr, infmsg + 6144, (ftnlen)1024, (ftnlen)1, 
	    (ftnlen)1024);
    repmc_(infmsg + 6144, "$", frame + 32, infmsg + 6144, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)32, (ftnlen)1024);
    repmc_(infmsg + 8192, "$", bodnam + 32, infmsg + 8192, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)32, (ftnlen)1024);
    repmi_(infmsg + 8192, "$", &bodid[1], infmsg + 8192, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)1024);
    repmc_(infmsg + 8192, "$", cennam + 32, infmsg + 8192, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)32, (ftnlen)1024);
    repmi_(infmsg + 8192, "$", &cenid[1], infmsg + 8192, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)1024);
    repmc_(infmsg + 9216, "$", spk + 1024, infmsg + 9216, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)1024, (ftnlen)1024);

/*     Step in included into the report as the number of seconds and as */
/*     #d #h #m #.######s (for more clarity.) */

    repmf_(infmsg + 11264, "$", &step, &c__15, "F", infmsg + 11264, (ftnlen)
	    1024, (ftnlen)1, (ftnlen)1, (ftnlen)1024);
    d__1 = spd_();
    rmaind_(&step, &d__1, &hdp1, &hdp2);
    day = (integer) hdp1;
    rmaind_(&hdp2, &c_b42, &hdp1, &hdp3);
    hour = (integer) hdp1;
    rmaind_(&hdp3, &c_b43, &hdp1, &hdp4);
    minuts = (integer) hdp1;
    rmaind_(&hdp4, &c_b44, &hdp1, &hdp2);
    sec = (integer) hdp1;
    d__1 = hdp2 * 1e6;
    msec = i_dnnt(&d__1);
    if (msec == 1000000) {
	msec = 999999;
    }
    intstr_(&msec, hword, (ftnlen)32);
    while(rtrim_(hword, (ftnlen)32) < 6) {
	prefix_("0", &c__0, hword, (ftnlen)1, (ftnlen)32);
    }
    s_copy(line, "#d #h #m #.#s", (ftnlen)1024, (ftnlen)13);
    repmi_(line, "#", &day, line, (ftnlen)1024, (ftnlen)1, (ftnlen)1024);
    repmi_(line, "#", &hour, line, (ftnlen)1024, (ftnlen)1, (ftnlen)1024);
    repmi_(line, "#", &minuts, line, (ftnlen)1024, (ftnlen)1, (ftnlen)1024);
    repmi_(line, "#", &sec, line, (ftnlen)1024, (ftnlen)1, (ftnlen)1024);
    repmc_(line, "#", hword, line, (ftnlen)1024, (ftnlen)1, (ftnlen)32, (
	    ftnlen)1024);
    repmc_(infmsg + 11264, "$", line, infmsg + 11264, (ftnlen)1024, (ftnlen)1,
	     (ftnlen)1024, (ftnlen)1024);
    repmc_(infmsg + 14336, "$", time, infmsg + 14336, (ftnlen)1024, (ftnlen)1,
	     (ftnlen)1024, (ftnlen)1024);
    repmf_(infmsg + 14336, "$", et, &c__15, "F", infmsg + 14336, (ftnlen)1024,
	     (ftnlen)1, (ftnlen)1, (ftnlen)1024);
    repmc_(infmsg + 15360, "$", time + 1024, infmsg + 15360, (ftnlen)1024, (
	    ftnlen)1, (ftnlen)1024, (ftnlen)1024);
    repmf_(infmsg + 15360, "$", &et[1], &c__15, "F", infmsg + 15360, (ftnlen)
	    1024, (ftnlen)1, (ftnlen)1, (ftnlen)1024);

/*     If no additional kernels were provided, we don't report them. */

    if (s_cmp(kernls, " ", (ftnlen)1024, (ftnlen)1) != 0) {
	repmc_(infmsg + 19456, "$", kernls, infmsg + 19456, (ftnlen)1024, (
		ftnlen)1, (ftnlen)1024, (ftnlen)1024);
	n = 21;
    } else {
	n = 17;
    }

/*     Report all details about comparison to the screen. */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tostdo_(infmsg + (((i__2 = i__ - 1) < 21 && 0 <= i__2 ? i__2 : s_rnge(
		"infmsg", i__2, "spkdiff_", (ftnlen)233)) << 10), (ftnlen)
		1024);
    }

/*     We are done with command line inputs. Time to get to real */
/*     business. */


/*     Calculate times for which we will compute states. */

    i__1 = nitr - 2;
    for (j = 0; j <= i__1; ++j) {
	epoch[(i__2 = j) < 1000000 && 0 <= i__2 ? i__2 : s_rnge("epoch", i__2,
		 "spkdiff_", (ftnlen)245)] = et[0] + (doublereal) j * step;
    }
    epoch[(i__1 = nitr - 1) < 1000000 && 0 <= i__1 ? i__1 : s_rnge("epoch", 
	    i__1, "spkdiff_", (ftnlen)248)] = et[1];

/*     Before we will go ahead, we will do a small sanity check to */
/*     verify that other kernels provided on the command line and loaded */
/*     earlier don't contain any SPKs that would allow to compute state */
/*     for the first and/or second body-center pair. */

/*     This check needs to be done only if we had any SPKs among other */
/*     kernels. */

    ktotal_("SPK", &spkcnt, (ftnlen)3);
    if (spkcnt > 0) {
	i__1 = nitr;
	for (j = 1; j <= i__1; ++j) {

/*           Use modified Nat's ZZGEOSEG to determine whether we have */
/*           enough SPK data loaded to compute geometric state of the */
/*           first body-center pair at this time without loading the */
/*           first SPK. */

	    zzgeoseg_(bodid, &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? 
		    i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)271)], 
		    cenid, &n, hanlst, dsclst);
	    if (n != 0) {

/*              Yes, we do. We should signal and error and exit. */

		etcal_(&epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
			s_rnge("epoch", i__2, "spkdiff_", (ftnlen)279)], line,
			 (ftnlen)1024);
		setmsg_("Oops ... It looks like we can compute geometric sta"
			"te of # with respect to # in '#' frame at '# TDB' (E"
			"T #), even without loading the first SPK file '#'. I"
			"t means that supporting kernel(s) '#' already contai"
			"n data for this center/body pair.", (ftnlen)240);
		errint_("#", bodid, (ftnlen)1);
		errint_("#", cenid, (ftnlen)1);
		errch_("#", frame, (ftnlen)1, (ftnlen)32);
		errch_("#", line, (ftnlen)1, (ftnlen)1024);
		errdp_("#", &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? 
			i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)292)]
			, (ftnlen)1);
		errch_("#", spk, (ftnlen)1, (ftnlen)1024);
		errch_("#", kernls, (ftnlen)1, (ftnlen)1024);
		sigerr_("SPICE(JEOPARDIZEDRUN1)", (ftnlen)22);
	    }

/*           Same for the second body/center pair. */

	    zzgeoseg_(&bodid[1], &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 
		    ? i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)302)], 
		    &cenid[1], &n, hanlst, dsclst);
	    if (n != 0) {

/*              Yes, we do. We should signal and error and exit. */

		etcal_(&epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
			s_rnge("epoch", i__2, "spkdiff_", (ftnlen)310)], line,
			 (ftnlen)1024);
		setmsg_("Oops ... It looks like we can compute geometric sta"
			"te of # with respect to # in '#' frame at '# TDB' (E"
			"T #), even without loading the second SPK file '#'. "
			"It means that supporting kernel(s) '#' already conta"
			"in data for this center/body pair.", (ftnlen)241);
		errint_("#", &bodid[1], (ftnlen)1);
		errint_("#", &cenid[1], (ftnlen)1);
		errch_("#", frame + 32, (ftnlen)1, (ftnlen)32);
		errch_("#", line, (ftnlen)1, (ftnlen)1024);
		errdp_("#", &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? 
			i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)323)]
			, (ftnlen)1);
		errch_("#", spk + 1024, (ftnlen)1, (ftnlen)1024);
		errch_("#", kernls, (ftnlen)1, (ftnlen)1024);
		sigerr_("SPICE(JEOPARDIZEDRUN2)", (ftnlen)22);
	    }
	}
    }

/*     Load first SPK file and compute NITR states from it. */

    spklef_(spk, &handle, (ftnlen)1024);
    i__1 = nitr;
    for (j = 1; j <= i__1; ++j) {

/*        Use modified Nat's ZZGEOSEG to determine whether we have */
/*        enough SPK data loaded to compute geometric state at this */
/*        time. */

	zzgeoseg_(bodid, &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 :
		 s_rnge("epoch", i__2, "spkdiff_", (ftnlen)346)], cenid, &n, 
		hanlst, dsclst);
	if (n == 0) {

/*           No, we don't. In this version of the program we signal */
/*           and error and exit. */

	    etcal_(&epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
		    s_rnge("epoch", i__2, "spkdiff_", (ftnlen)355)], line, (
		    ftnlen)1024);
	    setmsg_("Geometric state of # with respect to # in '#' frame at "
		    "'# TDB' (ET #), which is within specified interval of in"
		    "terest, cannot be computed due to insufficient data in t"
		    "he first SPK file '#'.", (ftnlen)189);
	    errint_("#", bodid, (ftnlen)1);
	    errint_("#", cenid, (ftnlen)1);
	    errch_("#", frame, (ftnlen)1, (ftnlen)32);
	    errch_("#", line, (ftnlen)1, (ftnlen)1024);
	    errdp_("#", &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
		    s_rnge("epoch", i__2, "spkdiff_", (ftnlen)366)], (ftnlen)
		    1);
	    errch_("#", spk, (ftnlen)1, (ftnlen)1024);
	    sigerr_("SPICE(UNSUFFSPK1DATA)", (ftnlen)21);
	}

/*        Thought the test above was incomplete, we may a feel a */
/*        bit more comfortable when we call SPKGEO now (there can */
/*        still be cases when position in one of the segments found */
/*        ZZGEOSEG was with respect to a frame orientation for which */
/*        cannot be determined using loaded PCK/CK/FRAMES data, but */
/*        we will leave it up to SPKGEO to determine this.) */

	spkgeo_(bodid, &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
		s_rnge("epoch", i__2, "spkdiff_", (ftnlen)380)], frame, cenid,
		 &state1[(i__3 = j * 6 - 6) < 6000000 && 0 <= i__3 ? i__3 : 
		s_rnge("state1", i__3, "spkdiff_", (ftnlen)380)], &lt, (
		ftnlen)32);
    }

/*     We are done with the first file. Unload it. */

    spkuef_(&handle);

/*     Load second SPK file and compute NITR states from it. */

    spklef_(spk + 1024, &handle, (ftnlen)1024);
    i__1 = nitr;
    for (j = 1; j <= i__1; ++j) {

/*        Use modified Nat's ZZGEOSEG to determine whether we have */
/*        enough SPK data loaded to compute geometric state at this */
/*        time. */

	zzgeoseg_(&bodid[1], &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? 
		i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)402)], &
		cenid[1], &n, hanlst, dsclst);
	if (n == 0) {

/*           No, we don't. In this version of the program we signal */
/*           and error and exit. */

	    etcal_(&epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
		    s_rnge("epoch", i__2, "spkdiff_", (ftnlen)411)], line, (
		    ftnlen)1024);
	    setmsg_("Geometric state of # with respect to # in '#' frame at "
		    "'# TDB' (ET #), which is within specified interval of in"
		    "terest, cannot be computed due to insufficient data in t"
		    "he second SPK file '#'.", (ftnlen)190);
	    errint_("#", &bodid[1], (ftnlen)1);
	    errint_("#", &cenid[1], (ftnlen)1);
	    errch_("#", frame + 32, (ftnlen)1, (ftnlen)32);
	    errch_("#", line, (ftnlen)1, (ftnlen)1024);
	    errdp_("#", &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? i__2 : 
		    s_rnge("epoch", i__2, "spkdiff_", (ftnlen)422)], (ftnlen)
		    1);
	    errch_("#", spk + 1024, (ftnlen)1, (ftnlen)1024);
	    sigerr_("SPICE(UNSUFFSPK2DATA)", (ftnlen)21);
	}

/*        Thought the test above was incomplete, we may a feel a */
/*        bit more comfortable when we call SPKGEO now (there can */
/*        still be cases when position in one of the segments found */
/*        ZZGEOSEG was with respect to a frame orientation for which */
/*        cannot be determined using loaded PCK/CK/FRAMES data, but */
/*        we will leave it up to SPKGEO to determine this.) */

	spkgeo_(&bodid[1], &epoch[(i__2 = j - 1) < 1000000 && 0 <= i__2 ? 
		i__2 : s_rnge("epoch", i__2, "spkdiff_", (ftnlen)436)], frame 
		+ 32, &cenid[1], &state2[(i__3 = j * 6 - 6) < 6000000 && 0 <= 
		i__3 ? i__3 : s_rnge("state2", i__3, "spkdiff_", (ftnlen)436)]
		, &lt, (ftnlen)32);
    }

/*     We are done with the second file. Unload it. */

    spkuef_(&handle);

/*     Pass state tables to the routine that will do analysis of the */
/*     differences and will print them to the screen. */

    stdiff_(state1, state2, &nitr, epoch, diftyp, timfmt, (ftnlen)32, (ftnlen)
	    1024);

/*     Check out. */

    chkout_("spkdiff", (ftnlen)7);

/*     We are done. :-) */

    return 0;
} /* MAIN__ */

/* Main program alias */ int spkdiff_ () { MAIN__ (); return 0; }
