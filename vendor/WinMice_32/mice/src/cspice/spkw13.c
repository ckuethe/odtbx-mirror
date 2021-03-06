/* spkw13.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__15 = 15;
static integer c__13 = 13;
static integer c__1 = 1;

/* $Procedure      SPKW13 ( Write SPK segment, type 13 ) */
/* Subroutine */ int spkw13_(integer *handle, integer *body, integer *center, 
	char *frame, doublereal *first, doublereal *last, char *segid, 
	integer *degree, integer *n, doublereal *states, doublereal *epochs, 
	ftnlen frame_len, ftnlen segid_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    extern logical even_(integer *);
    integer i__;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    doublereal descr[5];
    extern /* Subroutine */ int errch_(char *, char *, ftnlen, ftnlen), 
	    errdp_(char *, doublereal *, ftnlen), dafada_(doublereal *, 
	    integer *), dafbna_(integer *, doublereal *, char *, ftnlen), 
	    dafena_(void);
    extern logical failed_(void);
    integer chrcod, refcod;
    extern /* Subroutine */ int namfrm_(char *, integer *, ftnlen);
    extern integer lastnb_(char *, ftnlen);
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    doublereal maxtim;
    extern /* Subroutine */ int setmsg_(char *, ftnlen), errint_(char *, 
	    integer *, ftnlen), spkpds_(integer *, integer *, char *, integer 
	    *, doublereal *, doublereal *, doublereal *, ftnlen);
    extern logical return_(void);
    integer winsiz;

/* $ Abstract */

/*     Write a type 13 segment to an SPK file. */

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
/*     SPC */
/*     SPK */
/*     TIME */

/* $ Keywords */

/*     EPHEMERIS */
/*     FILES */

/* $ Declarations */
/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   Handle of an SPK file open for writing. */
/*     BODY       I   NAIF code for an ephemeris object. */
/*     CENTER     I   NAIF code for center of motion of BODY. */
/*     FRAME      I   Reference frame name. */
/*     FIRST      I   Start time of interval covered by segment. */
/*     LAST       I   End time of interval covered by segment. */
/*     SEGID      I   Segment identifier. */
/*     DEGREE     I   Degree of interpolating polynomials. */
/*     N          I   Number of states. */
/*     STATES     I   Array of states. */
/*     EPOCHS     I   Array of epochs corresponding to states. */
/*     MAXDEG     P   Maximum allowed degree of interpolating polynomial. */

/* $ Detailed_Input */

/*     HANDLE         is the file handle of an SPK file that has been */
/*                    opened for writing. */

/*     BODY           is the NAIF integer code for an ephemeris object */
/*                    whose state relative to another body is described */
/*                    by the segment to be created. */

/*     CENTER         is the NAIF integer code for the center of motion */
/*                    of the object identified by BODY. */

/*     FRAME          is the NAIF name for a reference frame */
/*                    relative to which the state information for BODY */
/*                    is specified. */

/*     FIRST, */
/*     LAST           are, respectively, the start and stop times of */
/*                    the time interval over which the segment defines */
/*                    the state of BODY. */

/*     SEGID          is the segment identifier.  An SPK segment */
/*                    identifier may contain up to 40 characters. */

/*     DEGREE         is the degree of the Hermite polynomials used to */
/*                    interpolate the states.  All components of the */
/*                    state vectors are interpolated by polynomials of */
/*                    fixed degree. */

/*     N              is the number of states in the input state vector */
/*                    array. */

/*     STATES         contains a time-ordered array of geometric states */
/*                    ( x, y, z, dx/dt, dy/dt, dz/dt, in kilometers and */
/*                    kilometers per second ) of BODY relative to CENTER, */
/*                    specified relative to FRAME. */

/*     EPOCHS         is an array of epochs corresponding to the members */
/*                    of the state array.  The epochs are specified as */
/*                    seconds past J2000, TDB. */

/* $ Detailed_Output */

/*     None.  See $Particulars for a description of the effect of this */
/*     routine. */

/* $ Parameters */

/*     MAXDEG         is the maximum allowed degree of the interpolating */
/*                    polynomial.  If the value of MAXDEG is increased, */
/*                    the SPICELIB routine SPKPV must be changed */
/*                    accordingly.  In particular, the size of the */
/*                    record passed to SPKRnn and SPKEnn must be */
/*                    increased, and comments describing the record size */
/*                    must be changed. */

/* $ Exceptions */

/*     If any of the following exceptions occur, this routine will return */
/*     without creating a new segment. */

/*     1) If FRAME is not a recognized name, the error */
/*        SPICE(INVALIDREFFRAME) is signaled. */

/*     2) If the last non-blank character of SEGID occurs past index 40, */
/*        the error SPICE(SEGIDTOOLONG) is signaled. */

/*     3) If SEGID contains any nonprintable characters, the error */
/*        SPICE(NONPRINTABLECHARS) is signaled. */

/*     4) If DEGREE is not at least 1 or is greater than MAXDEG, the */
/*        error SPICE(INVALIDDEGREE) is signaled. */

/*     5) If DEGREE is not odd, the error SPICE(INVALIDDEGREE) is */
/*        signaled. */

/*     6) If the number of states N is not at least (DEGREE+1)/2, */
/*        the error SPICE(TOOFEWSTATES) will be signaled. */

/*     7) If FIRST is greater than or equal to LAST then the error */
/*        SPICE(BADDESCRTIMES) will be signaled. */

/*     8) If the elements of the array EPOCHS are not in strictly */
/*        increasing order, the error SPICE(TIMESOUTOFORDER) will be */
/*        signaled. */

/*     9) If the first epoch EPOCHS(1) is greater than FIRST, the error */
/*        SPICE(BADDESCRTIMES) will be signaled. */

/*    10) If the last epoch EPOCHS(N) is less than LAST, the error */
/*        SPICE(BADDESCRTIMES) will be signaled. */


/* $ Files */

/*     A new type 13 SPK segment is written to the SPK file attached */
/*     to HANDLE. */

/* $ Particulars */

/*     This routine writes an SPK type 13 data segment to the open SPK */
/*     file according to the format described in the type 13 section of */
/*     the SPK Required Reading. The SPK file must have been opened with */
/*     write access. */

/* $ Examples */

/*     Suppose that you have states and are prepared to produce */
/*     a segment of type 13 in an SPK file. */

/*     The following code fragment could be used to add the new segment */
/*     to a previously opened SPK file attached to HANDLE. The file must */
/*     have been opened with write access. */

/*        C */
/*        C     Create a segment identifier. */
/*        C */
/*                  SEGID = 'MY_SAMPLE_SPK_TYPE_13_SEGMENT' */

/*        C */
/*        C     Write the segment. */
/*        C */
/*              CALL SPKW13 (  HANDLE,  BODY,    CENTER,  FRAME, */
/*             .               FIRST,   LAST,    SEGID,   DEGREE, */
/*             .               N,       STATES,  EPOCHS          ) */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.0, 20-MAR-2000 (NJB) */

/* -& */
/* $ Index_Entries */

/*     write spk type_13 ephemeris data segment */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Local variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("SPKW13", (ftnlen)6);
    }

/*     Set the window size corresponding to the input degree.  This */
/*     size will be used in various places below. */

    winsiz = (*degree + 1) / 2;

/*     Get the NAIF integer code for the reference frame. */

    namfrm_(frame, &refcod, frame_len);
    if (refcod == 0) {
	setmsg_("The reference frame # is not supported.", (ftnlen)39);
	errch_("#", frame, (ftnlen)1, frame_len);
	sigerr_("SPICE(INVALIDREFFRAME)", (ftnlen)22);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     Check to see if the segment identifier is too long. */

    if (lastnb_(segid, segid_len) > 40) {
	setmsg_("Segment identifier contains more than 40 characters.", (
		ftnlen)52);
	sigerr_("SPICE(SEGIDTOOLONG)", (ftnlen)19);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     Now check that all the characters in the segment identifier */
/*     can be printed. */

    i__1 = lastnb_(segid, segid_len);
    for (i__ = 1; i__ <= i__1; ++i__) {
	chrcod = *(unsigned char *)&segid[i__ - 1];
	if (chrcod < 32 || chrcod > 126) {
	    setmsg_("The segment identifier contains nonprintable characters",
		     (ftnlen)55);
	    sigerr_("SPICE(NONPRINTABLECHARS)", (ftnlen)24);
	    chkout_("SPKW13", (ftnlen)6);
	    return 0;
	}
    }

/*     Make sure that the degree of the interpolating polynomials is */
/*     in range. */

    if (*degree < 1 || *degree > 15) {
	setmsg_("The interpolating polynomials have degree #; the valid degr"
		"ee range is [1, #]", (ftnlen)77);
	errint_("#", degree, (ftnlen)1);
	errint_("#", &c__15, (ftnlen)1);
	sigerr_("SPICE(INVALIDDEGREE)", (ftnlen)20);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     Make sure that the degree of the interpolating polynomials is odd. */

    if (even_(degree)) {
	setmsg_("The interpolating polynomials have degree #; for SPK type 1"
		"3, the degree must be odd.", (ftnlen)85);
	errint_("#", degree, (ftnlen)1);
	sigerr_("SPICE(INVALIDDEGREE)", (ftnlen)20);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     Make sure that the number of states is sufficient to define a */
/*     polynomial whose degree is DEGREE. */

    if (*n < winsiz) {
	setmsg_("At least # states are required to define a Hermite polynomi"
		"al of degree #.  Number of states supplied:  #", (ftnlen)105);
	errint_("#", &winsiz, (ftnlen)1);
	errint_("#", degree, (ftnlen)1);
	errint_("#", n, (ftnlen)1);
	sigerr_("SPICE(TOOFEWSTATES)", (ftnlen)19);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     The segment stop time should be greater then the begin time. */

    if (*first >= *last) {
	setmsg_("The segment start time: # is greater then the segment end t"
		"ime: #", (ftnlen)65);
	errdp_("#", first, (ftnlen)1);
	errdp_("#", last, (ftnlen)1);
	sigerr_("SPICE(BADDESCRTIMES)", (ftnlen)20);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     Make sure the epochs form a strictly increasing sequence. */

    maxtim = epochs[0];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (epochs[i__ - 1] <= maxtim) {
	    setmsg_("EPOCH # having index # is not greater than its predeces"
		    "sor #.", (ftnlen)61);
	    errdp_("#", &epochs[i__ - 1], (ftnlen)1);
	    errint_("#", &i__, (ftnlen)1);
	    errdp_("#", &epochs[i__ - 2], (ftnlen)1);
	    sigerr_("SPICE(TIMESOUTOFORDER)", (ftnlen)22);
	    chkout_("SPKW13", (ftnlen)6);
	    return 0;
	} else {
	    maxtim = epochs[i__ - 1];
	}
    }

/*     Make sure that the span of the input epochs includes the interval */
/*     defined by the segment descriptor. */

    if (epochs[0] > *first) {
	setmsg_("Segment start time # precedes first epoch #.", (ftnlen)44);
	errdp_("#", first, (ftnlen)1);
	errdp_("#", epochs, (ftnlen)1);
	sigerr_("SPICE(BADDESCRTIMES)", (ftnlen)20);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    } else if (epochs[*n - 1] < *last) {
	setmsg_("Segment end time # follows last epoch #.", (ftnlen)40);
	errdp_("#", last, (ftnlen)1);
	errdp_("#", &epochs[*n - 1], (ftnlen)1);
	sigerr_("SPICE(BADDESCRTIMES)", (ftnlen)20);
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     If we made it this far, we're ready to start writing the segment. */


/*     Create the segment descriptor. */

    spkpds_(body, center, frame, &c__13, first, last, descr, frame_len);

/*     Begin a new segment. */

    dafbna_(handle, descr, segid, segid_len);
    if (failed_()) {
	chkout_("SPKW13", (ftnlen)6);
	return 0;
    }

/*     The type 13 segment structure is eloquently described by this */
/*     diagram from the SPK Required Reading: */

/*        +-----------------------+ */
/*        | State 1               | */
/*        +-----------------------+ */
/*        | State 2               | */
/*        +-----------------------+ */
/*                    . */
/*                    . */
/*                    . */
/*        +-----------------------+ */
/*        | State N               | */
/*        +-----------------------+ */
/*        | Epoch 1               | */
/*        +-----------------------+ */
/*        | Epoch 2               | */
/*        +-----------------------+ */
/*                    . */
/*                    . */
/*                    . */
/*        +-----------------------+ */
/*        | Epoch N               | */
/*        +-----------------------+ */
/*        | Epoch 100             | (First directory) */
/*        +-----------------------+ */
/*                    . */
/*                    . */
/*                    . */
/*        +-----------------------+ */
/*        | Epoch ((N-1)/100)*100 | (Last directory) */
/*        +-----------------------+ */
/*        | Window size - 1       | */
/*        +-----------------------+ */
/*        | Number of states      | */
/*        +-----------------------+ */


    i__1 = *n * 6;
    dafada_(states, &i__1);
    dafada_(epochs, n);
    i__1 = (*n - 1) / 100;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dafada_(&epochs[i__ * 100 - 1], &c__1);
    }
    d__1 = (doublereal) (winsiz - 1);
    dafada_(&d__1, &c__1);
    d__1 = (doublereal) (*n);
    dafada_(&d__1, &c__1);

/*     As long as nothing went wrong, end the segment. */

    if (! failed_()) {
	dafena_();
    }
    chkout_("SPKW13", (ftnlen)6);
    return 0;
} /* spkw13_ */

