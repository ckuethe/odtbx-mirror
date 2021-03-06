/* dasudi.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;

/* $Procedure      DASUDI ( DAS, update data, integer ) */
/* Subroutine */ int dasudi_(integer *handle, integer *first, integer *last, 
	integer *data)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer n;
    extern /* Subroutine */ int chkin_(char *, ftnlen);
    integer lastc, lastd, recno, lasti;
    extern /* Subroutine */ int dasa2l_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern logical failed_(void);
    integer clbase;
    extern /* Subroutine */ int daslla_(integer *, integer *, integer *, 
	    integer *), dasuri_(integer *, integer *, integer *, integer *, 
	    integer *);
    integer clsize;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen), setmsg_(char *, ftnlen), errint_(char *, integer *, 
	    ftnlen);
    integer wordno, numint;
    extern logical return_(void);
    integer nwritn;

/* $ Abstract */

/*     Update data in a specified range of integer addresses in a DAS */
/*     file. */

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

/*     DAS */

/* $ Keywords */

/*     ARRAY */
/*     ASSIGNMENT */
/*     DAS */
/*     FILES */

/* $ Declarations */
/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   DAS file handle. */
/*     FIRST, */
/*     LAST       I   Range of integer addresses to write to. */
/*     DATA       I   An array of integers. */

/* $ Detailed_Input */

/*     HANDLE         is a file handle of a DAS file opened for writing. */

/*     FIRST, */
/*     LAST           are the first and last of a range of DAS logical */
/*                    addresses of integers.  These addresses satisfy the */
/*                    inequality */

/*                       1  <   FIRST   <   LAST   <   LASTI */
/*                          _           -          - */

/*                    where LASTI is the last integer logical address in */
/*                    use in the DAS file designated by HANDLE. */

/*     DATA           is an array of integers.  The array elements */
/*                    DATA(1) through DATA(N) will be written to the */
/*                    indicated DAS file, where N is LAST - FIRST + 1. */

/* $ Detailed_Output */

/*     See $Particulars for a description of the effect of this routine. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the input file handle is invalid, the error will be */
/*         diagnosed by routines called by this routine. */

/*     2)  Only logical addresses that already contain data may be */
/*         updated:  if either FIRST or LAST are outside the range */

/*           [ 1,  LASTI ] */

/*         where LASTI is the last integer logical address that */
/*         currently contains data in the indicated DAS file, the error */
/*         SPICE(INVALIDADDRESS) is signalled.  The DAS file will not be */
/*         modified. */

/*     3)  If FIRST > LAST but both addresses are valid, this routine */
/*         will not modify the indicated DAS file.  No error will be */
/*         signalled. */

/*     4)  If an I/O error occurs during the data update attempted */
/*         by this routine, the error will be diagnosed by routines */
/*         called by this routine.  FIRST and LAST will not be modified. */

/* $ Files */

/*     See the description of the argument HANDLE in $Detailed_Input. */

/* $ Particulars */

/*     This routine replaces the integer data in the specified range of */
/*     logical addresses within a DAS file with the contents of the */
/*     input array DATA. */

/*     The actual physical write operations that update the indicated */
/*     DAS file with the contents of the input array DATA may not take */
/*     place before this routine returns, since the DAS system buffers */
/*     data that is written as well as data that is read.  In any case, */
/*     the data will be flushed to the file at the time the file is */
/*     closed, if not earlier.  A physical write of all buffered */
/*     records can be forced by calling the SPICELIB routine DASWUR */
/*     ( DAS, write updated records ). */

/*     In order to append integer data to a DAS file, filling in a range */
/*     of integer logical addresses that starts immediately after the */
/*     last integer logical address currently in use, the SPICELIB */
/*     routine DASADI ( DAS add data, integer ) should be used. */

/* $ Examples */

/*     1)  Write to addresses 1 through 500 in a DAS file in */
/*         random-access fashion by updating the file.  Recall */
/*         that data must be present in the file before it can */
/*         be updated. */


/*                  PROGRAM UP */

/*                  CHARACTER*(4)         TYPE */

/*                  INTEGER               DATA    ( 500 ) */

/*                  INTEGER               HANDLE */
/*                  INTEGER               I */

/*            C */
/*            C     Open the new DAS file RAND.DAS.  Use the file name */
/*            C     as the internal file name. */
/*            C */
/*                  TYPE = 'TEST' */
/*                  CALL DASONW ( 'TEST.DAS', TYPE, 'TEST.DAS', HANDLE ) */

/*            C */
/*            C     Append 500 integers to the file; after the data is */
/*            C     present, we're free to update it in any order we */
/*            C     please.  (CLEARI zeros out an integer array.) */
/*            C */
/*                  CALL CLEARI (           500,  DATA ) */
/*                  CALL DASADI (  HANDLE,  500,  DATA ) */

/*            C */
/*            C     Now the integer logical addresses 1:500 can be */
/*            C     written to in random-access fashion.  We'll fill them */
/*            C     in in reverse order. */
/*            C */
/*                  DO I = 500, 1, -1 */
/*                     CALL DASUDI ( HANDLE, I, I, I ) */
/*                  END DO */

/*            C */
/*            C     Close the file. */
/*            C */
/*                  CALL DASCLS ( HANDLE ) */

/*            C */
/*            C     Now make sure that we updated the file properly. */
/*            C     Open the file for reading and dump the contents */
/*            C     of the integer logical addresses 1:500. */
/*            C */
/*                  CALL DASOPR ( 'RAND.DAS',  HANDLE      ) */

/*                  CALL CLEARI (              500,  DATA  ) */
/*                  CALL DASRDI (  HANDLE,  1, 500,  DATA  ) */

/*                  WRITE (*,*) 'Contents of RAND.DAS:' */
/*                  WRITE (*,*) ' ' */
/*                  WRITE (*,*) DATA */

/*                  END */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     K.R. Gehringer (JPL) */
/*     N.J. Bachman   (JPL) */
/*     W.L. Taber     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.1.1 19-DEC-1995 (NJB) */

/*        Corrected title of permuted index entry section. */

/* -    SPICELIB Version 1.1.0, 12-MAY-1994 (KRG) (NJB) */

/*        Test of FAILED() added to loop termination conditions. */

/*        Removed references to specific DAS file open routines in the */
/*        $ Detailed_Input section of the header. This was done in order */
/*        to minimize documentation changes if the DAS open routines ever */
/*        change. */

/*        Modified the $ Examples section to demonstrate the new ID word */
/*        format which includes a file type and to include a call to the */
/*        new routine DASONW, open new for write, which makes use of the */
/*        file type. Also,  a variable for the type of the file to be */
/*        created was added. */

/* -    SPICELIB Version 1.0.0, 11-NOV-1992 (NJB) (WLT) */

/* -& */
/* $ Index_Entries */

/*     update integer data in a DAS file */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 1.1.0, 12-MAY-1994 (KRG) (NJB) */

/*        Test of FAILED() added to loop termination condition.  Without */
/*        this test, an infinite loop could result if DASA2L or DASURI */
/*        signaled an error inside the loop. */

/*        Removed references to specific DAS file open routines in the */
/*        $ Detailed_Input section of the header. This was done in order */
/*        to minimize documentation changes if the DAS open routines ever */
/*        change. */

/*        Modified the $ Examples section to demonstrate the new ID word */
/*        format which includes a file type and to include a call to the */
/*        new routine DASONW, open new for write, which makes use of the */
/*        file type. Also,  a variable for the type of the file to be */
/*        created was added. */

/* -    SPICELIB Version 1.0.0, 11-NOV-1992 (NJB) (WLT) */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DASUDI", (ftnlen)6);
    }

/*     Get the last logical addresses in use in this DAS file. */

    daslla_(handle, &lastc, &lastd, &lasti);

/*     Validate the input addresses. */

    if (*first < 1 || *first > lasti || *last < 1 || *last > lasti) {
	setmsg_("FIRST was #. LAST was #. Valid range is [1,#].", (ftnlen)46);
	errint_("#", first, (ftnlen)1);
	errint_("#", last, (ftnlen)1);
	errint_("#", &lasti, (ftnlen)1);
	sigerr_("SPICE(INVALIDADDRESS)", (ftnlen)21);
	chkout_("DASUDI", (ftnlen)6);
	return 0;
    }

/*     Let N be the number of addresses to update. */

    n = *last - *first + 1;

/*     We will use the variables RECNO and OFFSET to determine where to */
/*     write data in the DAS file.  RECNO will be the record containing */
/*     the physical location to write to;  WORDNO will be the word */
/*     location that we will write to next. */

/*     Find the first location to write to.  CLBASE and CLSIZE are the */
/*     base record number and size of the cluster of integer records that */
/*     the address FIRST lies within. */

    dasa2l_(handle, &c__3, first, &clbase, &clsize, &recno, &wordno);

/*     Set the number of integer words already written.  Keep */
/*     writing to the file until this number equals the number of */
/*     elements in DATA. */

/*     Note that if N is non-positive, the loop doesn't get exercised. */


    nwritn = 0;
    while(nwritn < n && ! failed_()) {

/*        Write as much data as we can (or need to) into the current */
/*        record.  We assume that CLBASE, RECNO, WORDNO, and NWRITN have */
/*        been set correctly at this point. */

/*        Find out how many words to write into the current record. */
/*        There may be no space left in the current record. */

/* Computing MIN */
	i__1 = n - nwritn, i__2 = 256 - wordno + 1;
	numint = min(i__1,i__2);
	if (numint > 0) {

/*           Write NUMINT words into the current record. */

	    i__1 = wordno + numint - 1;
	    dasuri_(handle, &recno, &wordno, &i__1, &data[nwritn]);
	    nwritn += numint;
	    wordno += numint;
	} else {

/*           It's time to start on a new record.  If the record we */
/*           just finished writing to (or just attempted writing to, */
/*           if it was full) was not the last of the cluster, the next */
/*           record to write to is the immediate successor of the last */
/*           one.  Otherwise, we'll have to look up the location of the */
/*           next integer logical address. */

	    if (recno < clbase + clsize - 1) {
		++recno;
		wordno = 1;
	    } else {
		i__1 = *first + nwritn;
		dasa2l_(handle, &c__3, &i__1, &clbase, &clsize, &recno, &
			wordno);
	    }
	}
    }
    chkout_("DASUDI", (ftnlen)6);
    return 0;
} /* dasudi_ */

