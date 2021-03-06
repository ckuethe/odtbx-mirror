/* rotget.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* $Procedure      ROTGET (Frame get transformation) */
/* Subroutine */ int rotget_(integer *infrm, doublereal *et, doublereal *
	rotate, integer *outfrm, logical *found)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    doublereal tipm[9]	/* was [3][3] */;
    integer type__, i__, j;
    extern /* Subroutine */ int zzdynrot_(integer *, integer *, doublereal *, 
	    doublereal *, integer *), chkin_(char *, ftnlen), errch_(char *, 
	    char *, ftnlen, ftnlen);
    char versn[6];
    extern /* Subroutine */ int xpose_(doublereal *, doublereal *);
    extern logical failed_(void);
    integer center;
    extern /* Subroutine */ int tipbod_(char *, integer *, doublereal *, 
	    doublereal *, ftnlen), namfrm_(char *, integer *, ftnlen), 
	    frinfo_(integer *, integer *, integer *, integer *, logical *), 
	    tkfram_(integer *, doublereal *, integer *, logical *), ckfrot_(
	    integer *, doublereal *, doublereal *, integer *, logical *), 
	    sigerr_(char *, ftnlen);
    integer typeid;
    extern /* Subroutine */ int chkout_(char *, ftnlen), setmsg_(char *, 
	    ftnlen), errint_(char *, integer *, ftnlen), irfrot_(integer *, 
	    integer *, doublereal *);
    extern logical return_(void);

/* $ Abstract */

/*     Find the rotation from a user specified frame to */
/*     another frame at a user specified epoch. */

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

/*     FRAMES */

/* $ Declarations */
/* $ Abstract */

/*     The parameters below form an enumerated list of the recognized */
/*     frame types.  They are: INERTL, PCK, CK, TK, DYN.  The meanings */
/*     are outlined below. */

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

/*     INERTL      an inertial frame that is listed in the routine */
/*                 CHGIRF and that requires no external file to */
/*                 compute the transformation from or to any other */
/*                 inertial frame. */

/*     PCK         is a frame that is specified relative to some */
/*                 INERTL frame and that has an IAU model that */
/*                 may be retrieved from the PCK system via a call */
/*                 to the routine TISBOD. */

/*     CK          is a frame defined by a C-kernel. */

/*     TK          is a "text kernel" frame.  These frames are offset */
/*                 from their associated "relative" frames by a */
/*                 constant rotation. */

/*     DYN         is a "dynamic" frame.  These currently are */
/*                 parameterized, built-in frames where the full frame */
/*                 definition depends on parameters supplied via a */
/*                 frame kernel. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 3.0.0, 28-MAY-2004 (NJB) */

/*       The parameter DYN was added to support the dynamic frame class. */

/* -    SPICELIB Version 2.0.0, 12-DEC-1996 (WLT) */

/*        Various unused frames types were removed and the */
/*        frame time TK was added. */

/* -    SPICELIB Version 1.0.0, 10-DEC-1995 (WLT) */

/* -& */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     INFRM      I   The integer code for a SPICE reference frame. */
/*     ET         I   An epoch in seconds past J2000. */
/*     ROTATE     O   A rotation matrix. */
/*     OUTFRM     O   The frame that ROTATE transforms INFRM to. */
/*     FOUND      O   TRUE if a rotation can be found. */

/* $ Detailed_Input */

/*     INFRM       is the SPICE id-code for some reference frame. */

/*     ET          is an epoch in ephemeris seconds past J2000 at */
/*                 which the user wishes to retrieve a transformation */
/*                 matrix. */

/* $ Detailed_Output */

/*     ROTATE      is a 3x3 matrix that transforms positions relative to */
/*                 INFRM to positions relative to OUTFRM.  (Assuming such */
/*                 a rotation can be found.) */

/*     OUTFRM      is a reference frame.  The 3x3 matrix ROTATE rotates */
/*                 positions relative to INFRM to positions relative */
/*                 to OUTFRM. */
/*                 The positions transformation is achieved by */
/*                 multiplying */
/*                 ROTATE on the right by a position relative to INFRM. */
/*                 This */
/*                 is easily accomplished via the subroutine call */
/*                 shown below. */

/*                    CALL MXV  ( ROTATE, INPOS,  OUTPOS ) */

/*     FOUND       is a logical flag indicating whether or not a */
/*                 rotation matrix could be found from INFRM */
/*                 to some other frame.  If a rotation matrix */
/*                 cannot be found OUTFRM will be set to zero, FOUND */
/*                 will be set to FALSE and ROTATE will be returned */
/*                 as the zero matrix. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1) If a rotation matrix cannot be located, then */
/*        FOUND will be set to FALSE, OUTFRM will be set to zero */
/*        and ROTATE will be set to the zero 3x3 matrix. */

/*     2) If the class of the requested frame is not recognized the */
/*        exception 'SPICE(UNKNOWNFRAMETYPE)' will be signalled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     This is a low level routine used for determining a chain of */
/*     position transformation matrices from one frame to another. */

/* $ Examples */

/*     See FRMCHG. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.1.0, 02-MAR-2010 (NJB) */

/*        Bug fix: frame ID rather than frame class ID */
/*        is now passed to dynamic frame evaluation */
/*        routine ZZDYNROT. Order of header sections was */
/*        corrected. */

/* -    SPICELIB Version 2.0.0, 18-DEC-2004 (NJB) */

/*        Added the new frame type 'DYN' to the list of frame */
/*        types recognized by ROTGET. */

/* -    SPICELIB Version 1.0.0, 03-MAR-1999 (WLT) */

/* -& */
/* $ Index_Entries */

/*     Find a rotation matrix from a specified frame */

/* -& */

/*     Spicelib Functions */


/*     Local Variables */

    s_copy(versn, "2.0.0", (ftnlen)6, (ftnlen)5);
    *found = FALSE_;

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    }
    chkin_("ROTGET", (ftnlen)6);

/*     Get all the needed information about this frame. */

    frinfo_(infrm, &center, &type__, &typeid, found);
    if (! (*found)) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		rotate[(i__1 = i__ + j * 3 - 4) < 9 && 0 <= i__1 ? i__1 : 
			s_rnge("rotate", i__1, "rotget_", (ftnlen)203)] = 0.;
	    }
	}
	chkout_("ROTGET", (ftnlen)6);
	return 0;
    }
    if (type__ == 1) {
	irfrot_(infrm, &c__1, rotate);
	*found = TRUE_;
	*outfrm = 1;
    } else if (type__ == 2) {
	tipbod_("J2000", &typeid, et, tipm, (ftnlen)5);
	xpose_(tipm, rotate);
	namfrm_("J2000", outfrm, (ftnlen)5);
	*found = ! failed_();
    } else if (type__ == 3) {
	ckfrot_(&typeid, et, rotate, outfrm, found);
    } else if (type__ == 4) {
	tkfram_(&typeid, rotate, outfrm, found);
    } else if (type__ == 5) {

/*        Unlike the other frame classes, the dynamic frame evaluation */
/*        routine ZZDYNROT requires the input frame ID rather than the */
/*        dynamic frame class ID. ZZDYNROT also requires the center ID */
/*        we found via the FRINFO call. */
	zzdynrot_(infrm, &center, et, rotate, outfrm);

/*        The FOUND flag was set by FRINFO earlier; we don't touch */
/*        it here.  If ZZDYNROT signaled an error, FOUND will be set */
/*        to .FALSE. at end of this routine. */

    } else {
	setmsg_("The reference frame # has class id-code #. This form of ref"
		"erence frame is not supported in version # of ROTGET. You ne"
		"ed to update your version of SPICELIB to the latest version "
		"in order to support this frame. ", (ftnlen)211);
	errint_("#", infrm, (ftnlen)1);
	errint_("#", &type__, (ftnlen)1);
	errch_("#", versn, (ftnlen)1, (ftnlen)6);
	sigerr_("SPICE(UNKNOWNFRAMETYPE)", (ftnlen)23);
	chkout_("ROTGET", (ftnlen)6);
	return 0;
    }
    if (failed_() || ! (*found)) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		rotate[(i__1 = i__ + j * 3 - 4) < 9 && 0 <= i__1 ? i__1 : 
			s_rnge("rotate", i__1, "rotget_", (ftnlen)272)] = 0.;
	    }
	}
	*found = FALSE_;
    }
    chkout_("ROTGET", (ftnlen)6);
    return 0;
} /* rotget_ */

