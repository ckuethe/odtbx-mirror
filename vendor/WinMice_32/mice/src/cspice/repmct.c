/* repmct.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* $Procedure  REPMCT  ( Replace marker with cardinal text ) */
/* Subroutine */ int repmct_(char *in, char *marker, integer *value, char *
	case__, char *out, ftnlen in_len, ftnlen marker_len, ftnlen case_len, 
	ftnlen out_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_indx(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    char card[145];
    extern /* Subroutine */ int lcase_(char *, char *, ftnlen, ftnlen), 
	    chkin_(char *, ftnlen), ucase_(char *, char *, ftnlen, ftnlen), 
	    errch_(char *, char *, ftnlen, ftnlen), ljust_(char *, char *, 
	    ftnlen, ftnlen);
    extern integer lastnb_(char *, ftnlen);
    char tmpcas[1];
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    extern integer frstnb_(char *, ftnlen);
    extern /* Subroutine */ int repsub_(char *, integer *, integer *, char *, 
	    char *, ftnlen, ftnlen, ftnlen), setmsg_(char *, ftnlen);
    integer mrkpos;
    extern logical return_(void);
    extern /* Subroutine */ int inttxt_(integer *, char *, ftnlen);

/* $ Abstract */

/*     Replace a marker with the text representation of a */
/*     cardinal number. */

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

/*     CHARACTER */
/*     CONVERSION */
/*     STRING */

/* $ Declarations */
/* $ Brief_I/O */

/*     VARIABLE  I/O  DESCRIPTION */
/*     --------  ---  -------------------------------------------------- */
/*     IN         I   Input string. */
/*     MARKER     I   Marker to be replaced. */
/*     VALUE      I   Cardinal value. */
/*     CASE       I   Case of replacement text. */
/*     OUT        O   Output string. */
/*     MAXLCN     P   Maximum length of a cardinal number. */

/* $ Detailed_Input */

/*     IN             is an arbitrary character string. */

/*     MARKER         is an arbitrary character string. The first */
/*                    occurrence of MARKER in the input string is */
/*                    to be replaced by the text representation of */
/*                    the cardinal number VALUE. */

/*                    Leading and trailing blanks in MARKER are NOT */
/*                    significant. In particular, no substitution is */
/*                    performed if MARKER is blank. */

/*     VALUE          is an arbitrary integer. */

/*     CASE           indicates the case of the replacement text. */
/*                    CASE may be any of the following: */

/*                       CASE     Meaning        Example */
/*                       ----     -----------    ----------------------- */
/*                       U, u     Uppercase      ONE HUNDRED FIFTY-THREE */

/*                       L, l     Lowercase      one hundred fifty-three */

/*                       C, c     Capitalized    One hundred fifty-three */

/* $ Detailed_Output */

/*     OUT            is the string obtained by substituting the text */
/*                    representation of the cardinal number VALUE for */
/*                    the first occurrence of MARKER in the input string. */

/*                    OUT and IN must be identical or disjoint. */

/* $ Parameters */

/*     MAXLCN         is the maximum expected length of any cardinal */
/*                    text. 145 characters are sufficient to hold the */
/*                    text representing any value in the range */

/*                      ( -10**12, 10**12 ) */

/*                    An example of a number whose text representation */
/*                    is of maximum length is */

/*                       - 777 777 777 777 */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1) If OUT does not have sufficient length to accommodate the */
/*        result of the substitution, the result will be truncated on */
/*        the right. */

/*     2) If MARKER is blank, or if MARKER is not a substring of IN, */
/*        no substitution is performed. (OUT and IN are identical.) */

/*     3) If the value of CASE is not recognized, the error */
/*        SPICE(INVALIDCASE) is signalled. OUT is not changed. */

/* $ Particulars */

/*     This is one of a family of related routines for inserting values */
/*     into strings. They are typically used to construct messages that */
/*     are partly fixed, and partly determined at run time. For example, */
/*     a message like */

/*        'Fifty-one pictures were found in directory [USER.DATA].' */

/*     might be constructed from the fixed string */

/*        '#1 pictures were found in directory #2.' */

/*     by the calls */

/*        CALL REPMCT ( STRING, '#1', NPICS,  'C', STRING ) */
/*        CALL REPMC  ( STRING, '#2', DIRNAM,      STRING ) */

/*     which substitute the cardinal text 'Fifty-one' and the character */
/*     string '[USER.DATA]' for the markers '#1' and '#2' respectively. */

/*     The complete list of routines is shown below. */

/*        REPMC    ( Replace marker with character string value ) */
/*        REPMD    ( Replace marker with double precision value ) */
/*        REPMF    ( Replace marker with formatted d.p. value ) */
/*        REPMI    ( Replace marker with integer value ) */
/*        REPMCT   ( Replace marker with cardinal text) */
/*        REPMOT   ( Replace marker with ordinal text ) */

/* $ Examples */

/*     The following examples illustrate the use of REPMCT to */
/*     replace a marker within a string with the cardinal text */
/*     corresponding to an integer. */

/*     Uppercase */
/*     --------- */

/*        Let */

/*           MARKER = '#' */
/*           IN     = 'INVALID COMMAND.  WORD # WAS NOT RECOGNIZED.' */

/*        Then following the call, */

/*           CALL REPMCT ( IN, '#', 5, 'U', IN  ) */

/*        IN is */

/*           'INVALID COMMAND.  WORD FIVE WAS NOT RECOGNIZED.' */

/*     Lowercase */
/*     --------- */

/*        Let */

/*           MARKER = ' XX ' */
/*           IN     = 'Word XX of the XX sentence was misspelled.' */

/*        Then following the call, */

/*           CALL REPMCT ( IN, '  XX  ', 5, 'L', OUT ) */

/*        OUT is */

/*           'Word five of the XX sentence was misspelled.' */


/*     Capitalized */
/*     ----------- */

/*        Let */

/*           MARKER = ' XX ' */
/*           IN     = 'Name:  YY.  Rank:  XX.' */

/*        Then following the calls, */

/*           CALL REPMC  ( IN,  'YY', 'Moriarty', OUT ) */
/*           CALL REPMCT ( OUT, 'XX',     1, 'C', OUT ) */

/*        OUT is */

/*           'Name:  Moriarty.  Rank:  One.' */

/* $ Restrictions */

/*     1) VALUE must be in the range accepted by subroutine INTTXT. */
/*        This range is currently */

/*           ( -10**12, 10**12 ) */

/*        Note that the endpoints of the interval are excluded. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */
/*     I.M. Underwood (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 1.0.0, 30-AUG-1990 (NJB) (IMU) */

/* -& */
/* $ Index_Entries */

/*     replace marker with cardinal text */

/* -& */

/*     SPICELIB functions */


/*     Local variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("REPMCT", (ftnlen)6);
    }

/*     Bail out if CASE is not recognized. */

    ljust_(case__, tmpcas, (ftnlen)1, (ftnlen)1);
    ucase_(tmpcas, tmpcas, (ftnlen)1, (ftnlen)1);
    if (*(unsigned char *)tmpcas != 'U' && *(unsigned char *)tmpcas != 'L' && 
	    *(unsigned char *)tmpcas != 'C') {
	setmsg_("Case (#) must be U, L, or C.", (ftnlen)28);
	errch_("#", case__, (ftnlen)1, (ftnlen)1);
	sigerr_("SPICE(INVALIDCASE)", (ftnlen)18);
	chkout_("REPMCT", (ftnlen)6);
	return 0;
    }

/*     If MARKER is blank, no substitution is possible. */

    if (s_cmp(marker, " ", marker_len, (ftnlen)1) == 0) {
	s_copy(out, in, out_len, in_len);
	chkout_("REPMCT", (ftnlen)6);
	return 0;
    }

/*     Locate the leftmost occurrence of MARKER, if there is one */
/*     (ignoring leading and trailing blanks). If MARKER is not */
/*     a substring of IN, no substitution can be performed. */

    i__1 = frstnb_(marker, marker_len) - 1;
    mrkpos = i_indx(in, marker + i__1, in_len, lastnb_(marker, marker_len) - 
	    i__1);
    if (mrkpos == 0) {
	s_copy(out, in, out_len, in_len);
	chkout_("REPMCT", (ftnlen)6);
	return 0;
    }

/*     Okay, CASE is recognized and MARKER has been found. */
/*     Generate the cardinal text corresponding to VALUE. */

    inttxt_(value, card, (ftnlen)145);

/*     CARD is always returned in upper case; change to the specified */
/*     case, if required. */

    if (*(unsigned char *)tmpcas == 'L') {
	lcase_(card, card, (ftnlen)145, (ftnlen)145);
    } else if (*(unsigned char *)tmpcas == 'C') {
	lcase_(card + 1, card + 1, (ftnlen)144, (ftnlen)144);
    }

/*     Replace MARKER with CARD. */

    i__1 = mrkpos + lastnb_(marker, marker_len) - frstnb_(marker, marker_len);
    repsub_(in, &mrkpos, &i__1, card, out, in_len, lastnb_(card, (ftnlen)145),
	     out_len);
    chkout_("REPMCT", (ftnlen)6);
    return 0;
} /* repmct_ */

