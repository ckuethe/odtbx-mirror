/* daffa.f -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1000 = 1000;
static integer c__1 = 1;
static integer c__128 = 128;

/* $Procedure DAFFA ( DAF, find array ) */
/* Subroutine */ int daffa_0_(int n__, integer *handle, doublereal *sum, char 
	*name__, logical *found, ftnlen name_len)
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical sthvnr[1000] = { FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,FALSE_,
	    FALSE_,FALSE_,FALSE_,FALSE_ };
    static integer stfptr = -1;
    static integer sthead = -1;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer free;
    static doublereal exdc[124];
    static integer exic[250], stfh[1000], prev;
    static char stnr[1000*1000];
    static doublereal stsr[128000]	/* was [128][1000] */;
    static integer i__, p;
    extern logical elemi_(integer *, integer *);
    extern /* Subroutine */ int chkin_(char *, ftnlen), dafps_(integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer bward;
    static doublereal newdc[124];
    extern /* Subroutine */ int dafus_(doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    static integer fward, newic[250];
    extern /* Subroutine */ int errch_(char *, char *, ftnlen, ftnlen), 
	    moved_(doublereal *, integer *, doublereal *), movei_(integer *, 
	    integer *, integer *);
    static integer nextp;
    static doublereal exsum[124];
    static integer nd;
    extern logical failed_(void);
    static char dafnam[255];
    static integer ni;
    extern /* Subroutine */ int dafhof_(integer *), dafhfn_(integer *, char *,
	     ftnlen), dafhsf_(integer *, integer *, integer *), dafsih_(
	    integer *, char *, ftnlen);
    static char ifname[60];
    extern /* Subroutine */ int dafrcr_(integer *, integer *, char *, ftnlen),
	     dafrfr_(integer *, integer *, integer *, char *, integer *, 
	    integer *, integer *, ftnlen), dafgsr_(integer *, integer *, 
	    integer *, integer *, doublereal *, logical *), dafwdr_(integer *,
	     integer *, doublereal *), dafwcr_(integer *, integer *, char *, 
	    ftnlen);
    static integer offset;
    extern /* Subroutine */ int sigerr_(char *, ftnlen), chkout_(char *, 
	    ftnlen);
    static integer namsiz;
    extern /* Subroutine */ int setmsg_(char *, ftnlen);
    static integer stnseg[1000];
    extern /* Subroutine */ int ssizei_(integer *, integer *);
    static integer opnset[1006];
    extern logical return_(void);
    static integer stthis[1000], stpool[1000], stcurr[1000], stprev[1000], 
	    stnext[1000], sumsiz;
    static logical fnd;

/* $ Abstract */

/*     Find arrays in a DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/* $ Abstract */

/*     Parameter declarations for the DAF/DAS handle manager. */

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

/*     DAF, DAS */

/* $ Keywords */

/*     PRIVATE */

/* $ Particulars */

/*     This include file contains parameters defining limits and */
/*     integer codes that are utilized in the DAF/DAS handle manager */
/*     routines. */

/* $ Restrictions */

/*     None. */

/* $ Author_and_Institution */

/*     F.S. Turner       (JPL) */

/* $ Literature_References */

/*     None. */

/* $ Version */

/* -    SPICELIB Version 1.20.0, 13-MAY-2010 (BVS) */

/*        Updated for SUN-SOLARIS-INTEL. */

/* -    SPICELIB Version 1.19.0, 13-MAY-2010 (BVS) */

/*        Updated for SUN-SOLARIS-INTEL-CC_C. */

/* -    SPICELIB Version 1.18.0, 13-MAY-2010 (BVS) */

/*        Updated for SUN-SOLARIS-INTEL-64BIT-CC_C. */

/* -    SPICELIB Version 1.17.0, 13-MAY-2010 (BVS) */

/*        Updated for SUN-SOLARIS-64BIT-NATIVE_C. */

/* -    SPICELIB Version 1.16.0, 13-MAY-2010 (BVS) */

/*        Updated for PC-WINDOWS-64BIT-IFORT. */

/* -    SPICELIB Version 1.15.0, 13-MAY-2010 (BVS) */

/*        Updated for PC-LINUX-64BIT-GFORTRAN. */

/* -    SPICELIB Version 1.14.0, 13-MAY-2010 (BVS) */

/*        Updated for PC-64BIT-MS_C. */

/* -    SPICELIB Version 1.13.0, 13-MAY-2010 (BVS) */

/*        Updated for MAC-OSX-64BIT-INTEL_C. */

/* -    SPICELIB Version 1.12.0, 13-MAY-2010 (BVS) */

/*        Updated for MAC-OSX-64BIT-IFORT. */

/* -    SPICELIB Version 1.11.0, 13-MAY-2010 (BVS) */

/*        Updated for MAC-OSX-64BIT-GFORTRAN. */

/* -    SPICELIB Version 1.10.0, 18-MAR-2009 (BVS) */

/*        Updated for PC-LINUX-GFORTRAN. */

/* -    SPICELIB Version 1.9.0, 18-MAR-2009 (BVS) */

/*        Updated for MAC-OSX-GFORTRAN. */

/* -    SPICELIB Version 1.8.0, 19-FEB-2008 (BVS) */

/*        Updated for PC-LINUX-IFORT. */

/* -    SPICELIB Version 1.7.0, 14-NOV-2006 (BVS) */

/*        Updated for PC-LINUX-64BIT-GCC_C. */

/* -    SPICELIB Version 1.6.0, 14-NOV-2006 (BVS) */

/*        Updated for MAC-OSX-INTEL_C. */

/* -    SPICELIB Version 1.5.0, 14-NOV-2006 (BVS) */

/*        Updated for MAC-OSX-IFORT. */

/* -    SPICELIB Version 1.4.0, 14-NOV-2006 (BVS) */

/*        Updated for PC-WINDOWS-IFORT. */

/* -    SPICELIB Version 1.3.0, 26-OCT-2005 (BVS) */

/*        Updated for SUN-SOLARIS-64BIT-GCC_C. */

/* -    SPICELIB Version 1.2.0, 03-JAN-2005 (BVS) */

/*        Updated for PC-CYGWIN_C. */

/* -    SPICELIB Version 1.1.0, 03-JAN-2005 (BVS) */

/*        Updated for PC-CYGWIN. */

/* -    SPICELIB Version 1.0.1, 17-JUL-2002 */

/*        Added MAC-OSX environments. */

/* -    SPICELIB Version 1.0.0, 07-NOV-2001 */

/* -& */

/*     Unit and file table size parameters. */

/*     FTSIZE     is the maximum number of files (DAS and DAF) that a */
/*                user may have open simultaneously. */


/*     RSVUNT     is the number of units protected from being locked */
/*                to a particular handle by ZZDDHHLU. */


/*     SCRUNT     is the number of units protected for use by scratch */
/*                files. */


/*     UTSIZE     is the maximum number of logical units this manager */
/*                will utilize at one time. */


/*     Access method enumeration.  These parameters are used to */
/*     identify which access method is associated with a particular */
/*     handle.  They need to be synchronized with the STRAMH array */
/*     defined in ZZDDHGSD in the following fashion: */

/*        STRAMH ( READ   ) = 'READ' */
/*        STRAMH ( WRITE  ) = 'WRITE' */
/*        STRAMH ( SCRTCH ) = 'SCRATCH' */
/*        STRAMH ( NEW    ) = 'NEW' */

/*     These values are used in the file table variable FTAMH. */


/*     Binary file format enumeration.  These parameters are used to */
/*     identify which binary file format is associated with a */
/*     particular handle.  They need to be synchronized with the STRBFF */
/*     array defined in ZZDDHGSD in the following fashion: */

/*        STRBFF ( BIGI3E ) = 'BIG-IEEE' */
/*        STRBFF ( LTLI3E ) = 'LTL-IEEE' */
/*        STRBFF ( VAXGFL ) = 'VAX-GFLT' */
/*        STRBFF ( VAXDFL ) = 'VAX-DFLT' */

/*     These values are used in the file table variable FTBFF. */


/*     Some random string lengths... more documentation required. */
/*     For now this will have to suffice. */


/*     Architecture enumeration.  These parameters are used to identify */
/*     which file architecture is associated with a particular handle. */
/*     They need to be synchronized with the STRARC array defined in */
/*     ZZDDHGSD in the following fashion: */

/*        STRARC ( DAF ) = 'DAF' */
/*        STRARC ( DAS ) = 'DAS' */

/*     These values will be used in the file table variable FTARC. */


/*     For the following environments, record length is measured in */
/*     characters (bytes) with eight characters per double precision */
/*     number. */

/*     Environment: Sun, Sun FORTRAN */
/*     Source:      Sun Fortran Programmer's Guide */

/*     Environment: PC, MS FORTRAN */
/*     Source:      Microsoft Fortran Optimizing Compiler User's Guide */

/*     Environment: Macintosh, Language Systems FORTRAN */
/*     Source:      Language Systems FORTRAN Reference Manual, */
/*                  Version 1.2, page 12-7 */

/*     Environment: PC/Linux, g77 */
/*     Source:      Determined by experiment. */

/*     Environment: PC, Lahey F77 EM/32 Version 4.0 */
/*     Source:      Lahey F77 EM/32 Language Reference Manual, */
/*                  page 144 */

/*     Environment: HP-UX 9000/750, FORTRAN/9000 Series 700 computers */
/*     Source:      FORTRAN/9000 Reference-Series 700 Computers, */
/*                  page 5-110 */

/*     Environment: NeXT Mach OS (Black Hardware), */
/*                  Absoft Fortran Version 3.2 */
/*     Source:      NAIF Program */


/*     The following parameter defines the size of a string used */
/*     to store a filenames on this target platform. */


/*     The following parameter controls the size of the character record */
/*     buffer used to read data from non-native files. */

/* $ Brief_I/O */

/*     Variable  I/O  Entry */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE    I,O  DAFBFS, DAFBBS, DAFGH, DAFCS */
/*     SUM       I,O  DAFGS,  DAFRS,  DAFWS */
/*     NAME      I,O  DAFGN,  DAFRN */
/*     FOUND      O   DAFFNA, DAFFPA */

/* $ Detailed_Input */

/*     HANDLE      on input is the handle of the DAF to be searched. */

/*     SUM         on input is an array summary that replaces the */
/*                 summary of the current array in the DAF currently */
/*                 being searched. */

/*     NAME        on input is an array name that replaces the name */
/*                 of the current array in the DAF currently being */
/*                 searched. */

/* $ Detailed_Output */

/*     HANDLE      on output is the handle of the DAF currently being */
/*                 searched. */

/*     SUM         on output is the summary for the array found most */
/*                 recently. */

/*     NAME        on output is the name for the array found */
/*                 most recently. */

/*     FOUND       is true whenever the search for the next or the */
/*                 previous array is successful, and is false otherwise. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     DAFs read by DAFFA and its entry points are opened */
/*     elsewhere, and referred to only by their handles. */

/* $ Exceptions */

/*     1) If DAFFA is called directly, the error SPICE(BOGUSENTRY) */
/*        is signalled. */

/*     2) See entry points DAFBFS, DAFFNA, DAFBBS, DAFFPA, DAFGS, DAFGN, */
/*        DAFGH, DAFRS, DAFWS, DAFRN, and DAFCS for exceptions specific */
/*        to those entry points. */

/* $ Particulars */

/*     DAFFA serves as an umbrella, allowing data to be shared by its */
/*     entry points: */

/*        DAFBFS         Begin forward search. */
/*        DAFFNA         Find next array. */

/*        DAFBBS         Begin backward search. */
/*        DAFFPA         Find previous array. */

/*        DAFGS          Get summary. */
/*        DAFGN          Get name. */
/*        DAFGH          Get handle. */

/*        DAFRS          Replace summary. */
/*        DAFWS          Write summary. */
/*        DAFRN          Replace name. */

/*        DAFCS          Continue search. */

/*     The main function of these entry points is to allow the */
/*     contents of any DAF to be examined on an array-by-array */
/*     basis. */

/*     Conceptually, the arrays in a DAF form a doubly linked list, */
/*     which can be searched in either of two directions: forward or */
/*     backward. It is possible to search multiple DAFs simultaneously. */

/*     DAFBFS (begin forward search) and DAFFNA are used to search the */
/*     arrays in a DAF in forward order.  In applications that search a */
/*     single DAF at a time, the normal usage is */

/*        CALL DAFBFS ( HANDLE ) */
/*        CALL DAFFNA ( FOUND  ) */

/*        DO WHILE ( FOUND ) */
/*           CALL DAFGS ( SUM  ) */
/*           CALL DAFGN ( NAME ) */
/*            . */
/*            . */

/*           CALL DAFFNA ( FOUND ) */
/*        END DO */



/*     DAFBBS (begin backward search) and DAFFPA are used to search the */
/*     arrays in a DAF in backward order.  In applications that search */
/*     a single DAF at a time, the normal usage is */

/*        CALL DAFBBS ( HANDLE ) */
/*        CALL DAFFPA ( FOUND  ) */

/*        DO WHILE ( FOUND ) */
/*           CALL DAFGS ( SUM  ) */
/*           CALL DAFGN ( NAME ) */
/*            . */
/*            . */

/*           CALL DAFFPA ( FOUND ) */
/*        END DO */


/*     In applications that conduct multiple searches simultaneously, */
/*     the above usage must be modified to specify the handle of the */
/*     file to operate on, in any case where the file may not be the */
/*     last one specified by DAFBFS or DAFBBS.  The routine DAFCS */
/*     (DAF, continue search) is used for this purpose.  Below, we */
/*     give an example of an interleaved search of two files specified */
/*     by the handles HANDL1 and HANDL2.  The directions of searches */
/*     in different DAFs are independent; here we conduct a forward */
/*     search on one file and a backward search on the other. */
/*     Throughout, we use DAFCS to specify which file to operate on, */
/*     before calling DAFFNA, DAFFPA, DAFGS, DAFRS, DAFWS, DAFGN, or */
/*     DAFRN. */


/*        CALL DAFBFS ( HANDL1 ) */
/*        CALL DAFBBS ( HANDL2 ) */

/*        CALL DAFCS  ( HANDL1 ) */
/*        CALL DAFFNA ( FOUND1 ) */

/*        CALL DAFCS  ( HANDL2 ) */
/*        CALL DAFFPA ( FOUND2 ) */

/*        DO WHILE ( FOUND1 .OR. FOUND2 ) */

/*           IF ( FOUND1 ) THEN */

/*              CALL DAFCS ( HANDL1 ) */
/*              CALL DAFGS ( SUM    ) */
/*              CALL DAFGN ( NAME   ) */
/*               . */
/*               . */
/*              CALL DAFCS  ( HANDL1 ) */
/*              CALL DAFFNA ( FOUND1 ) */

/*           END IF */

/*           IF ( FOUND2 ) THEN */

/*              CALL DAFCS ( HANDL2 ) */
/*              CALL DAFGS ( SUM    ) */
/*              CALL DAFGN ( NAME   ) */
/*               . */
/*               . */
/*              CALL DAFCS  ( HANDL2 ) */
/*              CALL DAFFPA ( FOUND2 ) */

/*           END IF */

/*        END DO */


/*     At any time, the latest array found (whether by DAFFNA or DAFFPA) */
/*     is regarded as the `current' array for the file in which the */
/*     array was found.  The last DAF in which a search was started, */
/*     executed, or continued by any of DAFBFS, DAFBBS, DAFFNA, DAFFPA */
/*     or DAFCS is regarded as the `current' DAF.  The summary and name */
/*     for the current array in the current DAF can be returned */
/*     separately, as shown above, by calls to DAFGS (get summary) and */
/*     DAFGN (get name).  The handle of the current DAF can also be */
/*     returned by calling DAFGH (get handle). */

/*     The summary and name of the current array in the current DAF can */
/*     be updated (again, separately) by providing new ones through DAFRS */
/*     (replace summary) and DAFRN (replace name). This feature */
/*     should not be used except to correct errors that occurred during */
/*     the creation of a file.  Note that changes can only be made to */
/*     files opened for write access. Also, the addresses of an array */
/*     cannot be changed using these routines. (Another routine, */
/*     DAFWS, is provided for this purpose, but should be used only */
/*     to reorder the arrays in a file.) */

/*     Once a search has been begun, it may be continued in either */
/*     direction. That is, DAFFPA may be used to back up during a */
/*     forward search, and DAFFNA may be used to advance during a */
/*     backward search. */

/* $ Examples */

/*     1) The following code fragment illustrates the way the entry */
/*        points of DAFFA might be used to edit the summaries and names */
/*        for the arrays contained in a DAF. (All subroutines and */
/*        functions are from SPICELIB.) */

/*        In this example, the user begins by supplying the name of */
/*        the file to be edited, followed by any number of the following */
/*        commands. */

/*           NEXT      finds the next array. */

/*           PREV      finds the previous array. */

/*           EDIT      changes the value of an item in the summary or */
/*                     of the entire name. The keyword EDIT is */
/*                     always followed by the name of the item to be */
/*                     edited, */

/*                        DC n */
/*                        IC n */
/*                        NAME */

/*                     and the value, e.g., */

/*                        EDIT IC 2 315 */
/*                        EDIT NAME NAIF test K2905-1 */

/*        The user may terminate the session at any time by typing END. */
/*        Commands other than those listed above are ignored. */

/*           READ (*,FMT='(A)') FNAME */
/*           CALL DAFOPW ( FNAME, HANDLE ) */
/*           CALL DAFBFS ( HANDLE ) */

/*           READ (*,FMT='(A)') COMMAND */

/*           DO WHILE ( COMMAND .NE. 'END' ) */
/*              CALL NEXTWD ( COMMAND, VERB, COMMAND ) */

/*              IF ( VERB .EQ. 'NEXT' ) THEN */
/*                 CALL DAFFNA ( FOUND ) */
/*                 IF ( .NOT. FOUND ) THEN */
/*                    WRITE (*,*) 'At end of array list.' */
/*                 END IF */

/*              IF ( VERB .EQ. 'PREV' ) THEN */
/*                 CALL DAFFPA ( FOUND ) */
/*                 IF ( .NOT. FOUND ) THEN */
/*                    WRITE (*,*) 'At beginning of array list.' */
/*                 END IF */

/*              IF ( VERB .EQ. 'EDIT' ) THEN */
/*                 CALL DAFGS ( SUM ) */
/*                 CALL DAFGN ( NAME ) */
/*                 CALL DAFUS ( SUM, ND, NI, DC, IC ) */

/*                 CALL NEXTWD ( COMMAND, ITEM, VALUE ) */

/*                 IF ( ITEM .EQ. 'DC' ) THEN */
/*                    CALL NEXTWD ( VALUE, INDEX, VALUE ) */
/*                    CALL NPARSI ( INDEX, LOC,     ERR, PTR ) */
/*                    CALL NPARSD ( VALUE, DC(LOC), ERR, PTR ) */

/*                 ELSE IF ( ITEM .EQ. 'IC' ) THEN */
/*                    CALL NEXTWD ( VALUE, INDEX, VALUE ) */
/*                    CALL NPARSI ( INDEX, LOC,     ERR, PTR ) */
/*                    CALL NPARSI ( VALUE, IC(LOC), ERR, PTR ) */

/*                 ELSE IF ( ITEM .EQ. 'NAME' ) THEN */
/*                    NAME = VALUE */
/*                 END IF */

/*                 CALL DAFPS ( ND, NI, DC, IC, SUM ) */
/*                 CALL DAFRS ( SUM ) */
/*                 CALL DAFRN ( NAME ) */
/*              END IF */

/*              READ (*,FMT='(A)') COMMAND */
/*           END DO */


/*     2)  The following program compares data in two DAFs.  The DAFs are */
/*         expected to have the same number of arrays, the same number */
/*         of elements in each corresponding array, and the same summary */
/*         format. */

/*         Each difference whose magnitude exceeds a specified tolerance */
/*         is flagged.  The difference information is written to a file. */


/*                  PROGRAM CMPDAF */

/*            C */
/*            C     Compare data in two DAFs having identical structures. */
/*            C     No array in either DAF is longer than ARRYSZ d.p. */
/*            C     numbers. */
/*            C */

/*            C */
/*            C     Local parameters */
/*            C */
/*                  INTEGER               ARRYSZ */
/*                  PARAMETER           ( ARRYSZ = 1000 ) */

/*                  INTEGER               ERRLEN */
/*                  PARAMETER           ( ERRLEN =  240 ) */

/*                  INTEGER               FILEN */
/*                  PARAMETER           ( FILEN  =  128 ) */

/*                  INTEGER               LINLEN */
/*                  PARAMETER           ( LINLEN =   80 ) */

/*                  INTEGER               MAXND */
/*                  PARAMETER           ( MAXND  =  125 ) */

/*                  INTEGER               MAXNI */
/*                  PARAMETER           ( MAXNI  =  250 ) */

/*                  INTEGER               MAXSUM */
/*                  PARAMETER           ( MAXSUM =  128 ) */

/*                  INTEGER               RLEN */
/*                  PARAMETER           ( RLEN   = 1000 ) */


/*            C */
/*            C     Local variables */
/*            C */
/*                  CHARACTER*(RLEN)      ANAME1 */
/*                  CHARACTER*(RLEN)      ANAME2 */
/*                  CHARACTER*(FILEN)     DAF1 */
/*                  CHARACTER*(FILEN)     DAF2 */
/*                  CHARACTER*(FILEN)     LOG */
/*                  CHARACTER*(ERRLEN)    PRSERR */
/*                  CHARACTER*(LINLEN)    STR */
/*                  CHARACTER*(LINLEN)    TOLCH */

/*                  DOUBLE PRECISION      ARRAY1 ( ARRYSZ ) */
/*                  DOUBLE PRECISION      ARRAY2 ( ARRYSZ ) */
/*                  DOUBLE PRECISION      DC1    ( MAXND ) */
/*                  DOUBLE PRECISION      DC2    ( MAXND ) */
/*                  DOUBLE PRECISION      TOL */
/*                  DOUBLE PRECISION      DIFF */
/*                  DOUBLE PRECISION      SUM1   ( MAXSUM ) */
/*                  DOUBLE PRECISION      SUM2   ( MAXSUM ) */

/*                  INTEGER               FA1 */
/*                  INTEGER               FA2 */
/*                  INTEGER               I */
/*                  INTEGER               IA1 */
/*                  INTEGER               IA2 */
/*                  INTEGER               IC1    ( MAXNI ) */
/*                  INTEGER               IC2    ( MAXNI ) */
/*                  INTEGER               FA */
/*                  INTEGER               HANDL1 */
/*                  INTEGER               HANDL2 */
/*                  INTEGER               LEN1 */
/*                  INTEGER               LEN2 */
/*                  INTEGER               ND1 */
/*                  INTEGER               ND2 */
/*                  INTEGER               NI1 */
/*                  INTEGER               NI2 */
/*                  INTEGER               PTR */

/*                  LOGICAL               FOUND */

/*            C */
/*            C     Start out by obtaining the names of the DAFs to be */
/*            C     compared. */
/*            C */
/*                  WRITE (*,*) 'Enter name of first DAF.' */
/*                  READ  (*,FMT='(A)') DAF1 */

/*                  WRITE (*,*) 'Enter name of second DAF.' */
/*                  READ  (*,FMT='(A)') DAF2 */

/*                  WRITE (*,*) 'Enter name of log file.' */
/*                  READ  (*,FMT='(A)') LOG */

/*                  WRITE (*,*) 'Enter tolerance for data comparison.' */
/*                  READ  (*,FMT='(A)') TOLCH */

/*                  CALL NPARSD ( TOLCH, TOL, PRSERR, PTR ) */

/*                  DO WHILE ( PRSERR .NE. ' ' ) */

/*                     WRITE (*,*) PRSERR */
/*                     WRITE (*,*) 'Enter tolerance for data comparison.' */
/*                     READ  (*,FMT='(A)') TOLCH */

/*                     CALL NPARSD ( TOLCH, TOL, PRSERR, PTR ) */

/*                  END DO */

/*            C */
/*            C     Open both DAFs for reading. */
/*            C */
/*                  CALL DAFOPR ( DAF1, HANDL1 ) */
/*                  CALL DAFOPR ( DAF2, HANDL2 ) */

/*            C */
/*            C     Start forward searches in both DAFS. */
/*            C */
/*                  CALL DAFBFS ( HANDL1 ) */
/*                  CALL DAFBFS ( HANDL2 ) */

/*            C */
/*            C     Obtain the summary formats for each DAF.  Stop now */
/*            C     if the summary formats don't match. */
/*            C */
/*                  CALL DAFHSF ( HANDL1, ND1, NI1 ) */
/*                  CALL DAFHSF ( HANDL2, ND2, NI2 ) */

/*                  IF (  ( ND1 .NE. ND2 ) .OR. ( NI1 .NE. NI2 )  ) THEN */

/*                     STR = 'Summary formats do not match.  NI1 = #, '// */
/*                 .                      'NI2 = #, ND1 = #, ND2 = #.' */

/*                     CALL REPMI  ( STR, '#', NI1, STR ) */
/*                     CALL REPMI  ( STR, '#', NI2, STR ) */
/*                     CALL REPMI  ( STR, '#', ND1, STR ) */
/*                     CALL REPMI  ( STR, '#', ND2, STR ) */

/*                     CALL WRLINE ( LOG,  STR ) */

/*                     CALL SIGERR ( 'Incompatible DAFs' ) */

/*                  END IF */

/*            C */
/*            C     Find the first array in each DAF.  Use DAFCS */
/*            C     (DAF, continue search) to set the handle of the DAF */
/*            C     to search in before calling DAFFNA. */
/*            C */
/*                  CALL DAFCS  ( HANDL1 ) */
/*                  CALL DAFFNA ( FOUND  ) */

/*                  IF ( FOUND ) THEN */
/*                     CALL DAFCS  ( HANDL2 ) */
/*                     CALL DAFFNA ( FOUND  ) */
/*                  END IF */

/*                  DO WHILE ( FOUND ) */

/*            C */
/*            C        Get the summary and name of each array, using */
/*            C        DAFCS to select the DAF to get the information */
/*            C        from.  Unpack the summaries and find the beginning */
/*            C        and ending addresses of the arrays.  Read the */
/*            C        arrays into the variables ARRAY1 and ARRAY2. */
/*            C */
/*                     CALL DAFCS ( HANDL1 ) */
/*                     CALL DAFGN ( ANAME1 ) */
/*                     CALL DAFGS ( SUM1   ) */
/*                     CALL DAFUS ( SUM1, ND1, NI1, DC1, IC1 ) */

/*                     IA1  = IC1 ( NI1 - 1 ) */
/*                     FA1  = IC1 ( NI1     ) */
/*                     LEN1 = FA1 - IA1  + 1 */

/*                     IF (  LEN1  .GT.  ARRYSZ  ) THEN */
/*                        CALL SETMSG ( 'Buffer too small; need # elts.') */
/*                        CALL ERRINT ( '#', LEN1                       ) */
/*                        CALL SIGERR ( 'ARRAYTOOSMALL'                 ) */
/*                     ELSE */
/*                        CALL DAFRDA ( HANDL1, IA1, FA1, ARRAY1 ) */
/*                     END IF */

/*                     CALL DAFCS ( HANDL2 ) */
/*                     CALL DAFGN ( ANAME2 ) */
/*                     CALL DAFGS ( SUM2   ) */
/*                     CALL DAFUS ( SUM2, ND2, NI2, DC2, IC2 ) */

/*                     IA2 = IC2 ( NI2 - 1 ) */
/*                     FA2 = IC2 ( NI2     ) */

/*                     LEN2 = FA2 - IA2  + 1 */

/*                     IF (  LEN1  .GT.  ARRYSZ  ) THEN */

/*                        CALL SETMSG ( 'Buffer too small; need # elts.') */
/*                        CALL ERRINT ( '#', LEN2                       ) */
/*                        CALL SIGERR ( 'ARRAYTOOSMALL'                 ) */

/*                     ELSE IF ( LEN1 .NE. LEN2 ) THEN */

/*                        CALL SETMSG ( 'DAF structures do not match. '// */
/*                    .                 'LEN1 = #, LEN2 = #. ' ) */
/*                        CALL ERRINT ( '#', LEN1              ) */
/*                        CALL ERRINT ( '#', LEN2              ) */
/*                        CALL SIGERR ( 'Incompatible DAFs' ) */

/*                     ELSE */
/*                        CALL DAFRDA ( HANDL2, IA2, FA2, ARRAY2 ) */
/*                     END IF */
/*            C */
/*            C */
/*            C        Compare the data in the two arrays.  Log a message */
/*            C        for every instance of data that differs by more */
/*            C        than the allowed tolerance.  Use the array names */
/*            C        to label the data sources. */
/*            C */
/*                     DO I = 1, LEN1 */

/*                        DIFF  =  ABS( ARRAY1(I) - ARRAY2(I) ) */

/*                        IF (  DIFF  .GT.  TOL  ) THEN */
/*            C */
/*            C              Get the array names. */
/*            C */
/*                           CALL DAFCS ( HANDL1 ) */
/*                           CALL DAFGN ( ANAME1 ) */
/*                           CALL DAFCS ( HANDL2 ) */
/*                           CALL DAFGN ( ANAME2 ) */

/*            C */
/*            C              Construct the report strings.  The number 14 */
/*            C              below is the number of significant digits to */
/*            C              show in the strings representing d.p. */
/*            C              numbers. */
/*            C */

/*                           CALL WRLINE ( LOG, ' ' ) */
/*                           CALL WRLINE ( LOG, 'Difference of array ' // */
/*                    .                         'elements exceeded '   // */
/*                    .                         'tolerance.'            ) */
/*                           CALL WRLINE ( LOG, 'First array:  '//ANAME1) */
/*                           CALL WRLINE ( LOG, 'Second array: '//ANAME2) */

/*                           STR = 'First value:  #' */
/*                           CALL REPMD  ( STR, '#', ARRAY1(I), 14, STR ) */
/*                           CALL WRLINE ( LOG, STR                     ) */

/*                           STR = 'Second value: #' */
/*                           CALL REPMD  ( STR, '#', ARRAY2(I), 14, STR ) */
/*                           CALL WRLINE ( LOG, STR                     ) */

/*                           STR = 'Difference:   #' */
/*                           CALL REPMD  ( STR, '#', DIFF,      14, STR ) */
/*                           CALL WRLINE ( LOG, STR                     ) */
/*                           CALL WRLINE ( LOG, ' '                     ) */

/*                        END IF */

/*                     END DO */

/*            C */
/*            C        Find the next pair of arrays. */
/*            C */
/*                     CALL DAFCS  ( HANDL1 ) */
/*                     CALL DAFFNA ( FOUND  ) */

/*                     IF ( FOUND ) THEN */
/*                        CALL DAFCS  ( HANDL2 ) */
/*                        CALL DAFFNA ( FOUND  ) */
/*                     END IF */

/*                  END DO */

/*            C */
/*            C     Close the DAFs. */
/*            C */
/*                  CALL DAFCLS ( HANDL1 ) */
/*                  CALL DAFCLS ( HANDL2 ) */

/*                  END */


/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 3.0.0, 16-NOV-2001 (FST) */

/*        Updated the entry points of DAFFA to enable its */
/*        internal state table size, TBSIZE, to be smaller */
/*        than the file table maintained by DAFAH: FTSIZE. */

/*        Calls to DAFRDR were replaced with the translation-aware */
/*        interface DAFGSR for retrieving summary records from */
/*        DAFs. */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     find daf array */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 3.0.0, 16-NOV-2001 (FST) */

/*        This umbrella and its entry points were updated to */
/*        work properly with the changes in the DAF system as */
/*        a result of its utilization of the new handle manager. */

/*        Since DAFAH now tracks FTSIZE files as defined in */
/*        the include file 'zzddhman.inc', it was decided that */
/*        in the interest of releasing the toolkit this module */
/*        would undergo simple changes.  As such most previous */
/*        references to FTSIZE in this umbrella have been replaced */
/*        with TBSIZE where appropriate.  DAFBFS and DAFBBS now signal */
/*        errors if there is not enough room to add a new DAF's */
/*        dossier to the state table.  Also, after attempting to */
/*        clean up all files listed in the state table that are */
/*        not currently open, DAFBFS and DAFBBS attempt to locate */
/*        the first dossier with STADDG set to FALSE.  This is then */
/*        freed to make room for the new DAF.  If DAFBNA fails */
/*        to locate such a dossier in the state table, it */
/*        signals the error SPICE(STFULL). */

/*        The parameter FILEN was removed, as it is defined */
/*        on an environmental basis in the include file */
/*        'zzddhman.inc'. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        In previous versions of DAFFA, only one search could be */
/*        conducted at a time.  Therefore, there was no question about */
/*        which DAF was being operated on by any of the DAFFA entry */
/*        points that don't accept file handles as input arguments. */
/*        In the current version of DAFFA, the entry points that don't */
/*        accept file handles as inputs operate on the `current DAF'. */
/*        The current DAF is the last one in which a search was */
/*        started by DAFBFS or DAFBBS, or continued by the new entry */
/*        point DAFCS.  DAFCS was added to allow users to set the */
/*        current DAF, so that searches of multiple DAFs can be */
/*        interleaved. */

/*        Note that the notion of `current DAF' as discussed here applies */
/*        only to DAFs acted upon by entry points of DAFFA.  In DAFANA, */
/*        there is a DAF that is treated as the `current DAF' for */
/*        adding data; there is no connection between the DAFs regarded */
/*        as current by DAFFA and DAFANA. */

/*        The two principal changes to DAFFA are the addition of the */
/*        new entry point DAFCS, and the addition of a data structure */
/*        called the `state table'.  The state table is a collection of */
/*        parallel arrays that maintain information about the state */
/*        of each search that is currently in progress.  The arrays are */
/*        indexed by a singly linked list pool; this mechanism allows */
/*        addition and deletion of information about searches without */
/*        requiring movement of data already in the state table.  The */
/*        linked list pool contains an `active' list and a `free' list. */
/*        Nodes in the active list are used to index elements of the */
/*        state table where data about searches in progress is stored. */
/*        The head node of the active list is of particular significance: */
/*        the state information pointed to by this node is that of the */
/*        current DAF.  Nodes in the free list index elements of the */
/*        state table that are available for use. */

/*        When a search is started on a DAF that is not already `known' */
/*        to DAFFA, information about the DAF is added to the state */
/*        table.  If there are no free elements in the state table, */
/*        the routine starting the search (DAFBFS or DAFBBS) will */
/*        perform garbage collection:  the routine will test the handles */
/*        of each file about which information in stored in the state */
/*        table to see whether that file is still open.  Nodes containing */
/*        information about DAFs that are no longer open will be moved */
/*        to the free list. */

/*        Whenever a DAF becomes the current DAF, the linked list */
/*        that indexes the state table is adjusted so that the */
/*        information about the current DAF is at the head of the list. */
/*        This way, a slight efficiency is gained when repeated search */
/*        accesses are made to the same DAF, since the linear search */
/*        through the state table for information on that DAF will */
/*        be shortened. */

/*        Since the algorithms for maintenance of linked lists are well */
/*        known, they are not documented here.  However, see the */
/*        internals of the SPICELIB routine SPKBSR for a nice diagram */
/*        describing a similar data structure. */

/*        The state table contains two arrays that are quite large: */
/*        there are buffers that contain the last character record */
/*        and summary record read from each DAF.  A parallel situation */
/*        exists in DAFANA, where the name and array summary for each */
/*        array under construction are buffered.  The total storage */
/*        required for these arrays (in DAFANA and DAFFA together) is */
/*        4000 * TBSIZE bytes.  For this reason, it may be a good idea */
/*        to reduce the value of TBSIZE in SPICELIB versions for */
/*        machines where memory is scarce. */

/* -& */

/*     SPICELIB functions */


/*     Local parameters */


/*     Local variables */


/*     State variables. */

/*     These variables define the state of each DAF to which data */
/*     is currently being added.  For each DAF that we're writing to, we */
/*     maintain a copy of: */

/*        STFH           File handle. */

/*        STPREV         Record number of previous array summary. */

/*        STTHIS         Record number of current array summary. */

/*        STNEXT         Record number of next array summary. */

/*        STNSEG         Number of summaries in current summary record. */

/*        STCURR         Index of current summary within summary record. */

/*        STNR           Last name record read. */

/*        STHVNR         Flag indicating whether name record containing */
/*                       name of current array is buffered. */

/*        STSR           Last summary record read. */

/*     These variables are maintained in a table of parallel arrays; */
/*     the size of the table is TBSIZE. */


/*     The table of state variables is indexed by a singly linked list */
/*     of pointers.  This mechanism avoids the work of moving */
/*     the state variable data about as information about DAFs is */
/*     added to or deleted from the table. */

/*     The structure containing the linked list pointers is called a */
/*     `pool'.  The pool contains a list of `active' nodes and a list */
/*     of free nodes.  The head nodes of the active and free lists are */
/*     maintained as the variables STHEAD (`state table head') and */
/*     STFPTR (`state table free pointer'), respectively.  Every node in */
/*     the pool is on exactly one of these lists. */


/*     The pool starts out with all of the nodes on the free list.  The */
/*     first one of DAFBFS or DAFBBS to be called initializes the pool. */
/*     As new DAFs are searched, DAFBFS and DAFBBS add information about */
/*     them to the state table.  Every time a search is started by DAFBFS */
/*     or DAFBBS, the routine in question `moves' the DAF's state */
/*     information to the head of the active list, if the state */
/*     information is not already there.  This re-organization is */
/*     accomplished by deleting the node for the DAF from its current */
/*     position in the active list and inserting the node at the head of */
/*     the list.  Thus, the change is made merely by setting pointers, */
/*     not by moving chunks of data in the state table. */

/*     It may happen that there is no room left in the state table */
/*     to accommodate information about a new DAF.  In this case, */
/*     garbage collection must be performed:  whichever of DAFBFS or */
/*     DAFBBS needs more room frees all nodes in the table that index */
/*     DAFs that are not currently open. */

/*     Note that the routines DAFGS, DAFGN, DAFRS, DAFRN, and DAFWS do */
/*     not modify the state table; they merely act on the current array */
/*     in the DAF that is at the head of the active list. */


/*     Other local variables */


/*     Save everything between calls */


/*     Initial values */

    /* Parameter adjustments */
    if (sum) {
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_dafbfs;
	case 2: goto L_daffna;
	case 3: goto L_dafbbs;
	case 4: goto L_daffpa;
	case 5: goto L_dafgs;
	case 6: goto L_dafgn;
	case 7: goto L_dafgh;
	case 8: goto L_dafrs;
	case 9: goto L_dafrn;
	case 10: goto L_dafws;
	case 11: goto L_dafcs;
	}


/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFFA", (ftnlen)5);
	sigerr_("SPICE(BOGUSENTRY)", (ftnlen)17);
	chkout_("DAFFA", (ftnlen)5);
    }
    return 0;
/* $Procedure DAFBFS ( DAF, begin forward search ) */

L_dafbfs:
/* $ Abstract */

/*     Begin a forward search for arrays in a DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     INTEGER               HANDLE */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   Handle of file to be searched. */

/* $ Detailed_Input */

/*     HANDLE      is the handle of a DAF on which a forward */
/*                 search is to be conducted. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     See argument HANDLE. */

/* $ Exceptions */

/*     1)  If the input handle is invalid, the error will be diagnosed */
/*         by routines called by this routine. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     begin daf forward search */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        Also, the $Exceptions section was filled out. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFBFS", (ftnlen)6);
    }

/*     Check out the file handle before going any further. */

    dafsih_(handle, "READ", (ftnlen)4);
    if (failed_()) {
	chkout_("DAFBFS", (ftnlen)6);
	return 0;
    }

/*     Initialize the state table pool, if this hasn't been done yet. */
/*     Also initialize the cell used to obtain the set of handles of */
/*     open DAFs. */

    if (first) {
	ssizei_(&c__1000, opnset);
	for (i__ = 1; i__ <= 999; ++i__) {
	    stpool[(i__1 = i__ - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stp"
		    "ool", i__1, "daffa_", (ftnlen)1123)] = i__ + 1;
	}
	stpool[999] = -1;
	stfptr = 1;
	first = FALSE_;
    }

/*     See whether we already have an entry for this DAF in the */
/*     state table.  Find the previous node if possible. */

    p = sthead;
    prev = -1;
    fnd = FALSE_;
    while(p != -1 && ! fnd) {
	if (stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
		i__1, "daffa_", (ftnlen)1142)] == *handle) {
	    fnd = TRUE_;
	} else {
	    prev = p;
	    p = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stpool", i__1, "daffa_", (ftnlen)1146)];
	}
    }

/*     At this point, either FND is false, or P points to a */
/*     state table entry describing the DAF indicated by HANDLE. */
/*     In the latter case, PREV is the predecessor of P. */

    if (fnd) {

/*        We already have a dossier on this DAF.  We already have */
/*        the information on the summary format, but we must re-set */
/*        our summary record pointers and our name record availability */
/*        flag. */

/*        Rather than doing the update here, we do it outside of this */
/*        IF block.  That way, the update gets done in just one place. */
/*        This just makes life easier:  if the collection of state */
/*        variables is changed, there are fewer places to forget to */
/*        make the required code changes. */

/*        Move the node for this DAF to the head of the active list, */
/*        if it is not already there: */

/*           - Make the predecessor of P point to the successor of P. */

/*           - Make P point to the head of the active list. */

/*           - Make P the active list head node. */


	if (p != sthead) {

/*           P is in the active list, but is not at the head.  So, */
/*           the predecessor of P is not NIL. */

	    stpool[(i__1 = prev - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stpool", i__1, "daffa_", (ftnlen)1184)] = stpool[(i__2 = 
		    p - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge("stpool", i__2,
		     "daffa_", (ftnlen)1184)];
	    stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool"
		    , i__1, "daffa_", (ftnlen)1185)] = sthead;
	    sthead = p;
	}
    } else {

/*        We don't yet have any information on this DAF.  Make a new */
/*        state table entry for the DAF.  We may need to make room for */
/*        the new information by freeing space allocated to DAFs that */
/*        are no longer open. */

	if (stfptr == -1) {

/*           Oops, we're out of space.  Time for garbage collection. */
/*           Test each file handle to see whether it designates a DAF */
/*           that is still open.  DAFHOF will tell us which handles */
/*           point to open DAFs. */

	    dafhof_(opnset);
	    p = sthead;
	    prev = -1;

/*           For every DAF file represented in the state table, we'll */
/*           delete the corresponding state information if the DAF is */
/*           now closed.  We traverse the active list, examining each */
/*           file handle as we go. */

	    while(p != -1) {
		if (elemi_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
			s_rnge("stfh", i__1, "daffa_", (ftnlen)1217)], opnset)
			) {

/*                 The file is open. Have a look at the next node. */

		    prev = p;
		    p = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
			    s_rnge("stpool", i__1, "daffa_", (ftnlen)1222)];
		} else {

/*                 This file handle is not on the list, so free the */
/*                 node pointing to the information about the DAF it */
/*                 designated: */

/*                    - Save the successor of P. */

/*                    - Link the predecessor of node P to the successor */
/*                      of P, if the predecessor is not NIL. */

/*                    - If it happens that P is the head node of the */
/*                      active list, set the head equal to the */
/*                      successor of P. */

/*                    - Link P into the free list. */

/*                    - Set P equal to its saved successor. */

/*                    - (PREV remains unchanged.) */


		    nextp = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 :
			     s_rnge("stpool", i__1, "daffa_", (ftnlen)1246)];
		    if (p == sthead) {

/*                    Re-assign STHEAD so that we don't lose the head */
/*                    of the active list.  P has no predecessor in this */
/*                    case, so there's no need to set the forward pointer */
/*                    of node PREV. */

			sthead = nextp;
		    } else {

/*                    Since P is not the head node of the active list, */
/*                    PREV is not NIL, so we'll need to set the forward */
/*                    pointer of node PREV. */

			stpool[(i__1 = prev - 1) < 1000 && 0 <= i__1 ? i__1 : 
				s_rnge("stpool", i__1, "daffa_", (ftnlen)1264)
				] = nextp;
		    }
		    stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
			    "stpool", i__1, "daffa_", (ftnlen)1269)] = stfptr;
		    stfptr = p;
		    p = nextp;
		}
	    }

/*           At this point, we've freed all nodes from the active */
/*           list that were used to index information about DAFs that */
/*           are no longer open.  If there's any more room in the state */
/*           table, we have it now. */

	}

/*        If there still is no room, there is a bug in DAFAH, since DAFAH */
/*        should not allow more than TBSIZE DAFs to be open.  So, we */
/*        assume that we've found some room.  The first free node is */
/*        indicated by STFPTR.  We'll allocate this node and use it to */
/*        index the state information for the new DAF. */

	p = stfptr;

/*        Update the free list pointer, link P to the previous head */
/*        of the active list, and make P the head of the active list. */

	stfptr = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stpool", i__1, "daffa_", (ftnlen)1297)];
	stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool", 
		i__1, "daffa_", (ftnlen)1298)] = sthead;
	sthead = p;
    }

/*     At this point, P is the head node of the active list, and P is */
/*     the index in the state table of the information for the current */
/*     DAF. */


/*     Read the file record and first summary record. Do not read the */
/*     corresponding name record until necessary. In most searches, */
/*     names are of no interest. */

    dafrfr_(handle, &nd, &ni, ifname, &fward, &bward, &free, (ftnlen)60);
    dafgsr_(handle, &fward, &c__1, &c__128, &stsr[(i__1 = (p << 7) - 128) < 
	    128000 && 0 <= i__1 ? i__1 : s_rnge("stsr", i__1, "daffa_", (
	    ftnlen)1316)], &fnd);

/*     Set up the state information for this file.  Note that we */
/*     don't have a name record yet, and we have no current array */
/*     yet. */

    stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", i__1, 
	    "daffa_", (ftnlen)1323)] = *handle;
    stthis[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stthis", i__1, 
	    "daffa_", (ftnlen)1324)] = fward;
    stnext[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnext", i__1, 
	    "daffa_", (ftnlen)1325)] = (integer) stsr[(i__2 = (p << 7) - 128) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1325)];
    stprev[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stprev", i__1, 
	    "daffa_", (ftnlen)1326)] = (integer) stsr[(i__2 = (p << 7) - 127) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1326)];
    stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnseg", i__1, 
	    "daffa_", (ftnlen)1327)] = (integer) stsr[(i__2 = (p << 7) - 126) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1327)];
    sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", i__1, 
	    "daffa_", (ftnlen)1328)] = FALSE_;

/*     The arrays are returned in forward order within each summary */
/*     record. */

    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", i__1, 
	    "daffa_", (ftnlen)1333)] = 0;
    chkout_("DAFBFS", (ftnlen)6);
    return 0;
/* $Procedure DAFFNA ( DAF, find next array ) */

L_daffna:
/* $ Abstract */

/*     Find the next (forward) array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     LOGICAL               FOUND */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     FOUND      O   True if an array was found. */

/* $ Detailed_Input */

/*     None. */

/* $ Detailed_Output */

/*     FOUND       is true if an array was found, and is false if, */
/*                 when this routine is called, the current array is */
/*                 the tail of the array list.  (Recall that the */
/*                 arrays in a DAF may be viewed as a doubly linked */
/*                 list, with the tail being the last array in the file.) */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called before a search is begun, the */
/*         error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF to be searched has actually been closed, the error */
/*         will be diagnosed by routines called by this routine. */

/*     3)  If the end of the array list has already been reached when */
/*         this routine is called, this routine has no effect. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     find next daf array */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFFNA", (ftnlen)6);
    }

/*     FOUND will be false until we make it past the error checks. */

    *found = FALSE_;

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFFNA", (ftnlen)6);
	return 0;

/*     Make sure that the `current' DAF is still open. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)1522)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFFNA", (ftnlen)6);
	    return 0;
	}
    }

/*     Now that we know a search is going on, assume that we will find */
/*     an array until proven otherwise. */

    *found = TRUE_;

/*     Either there are more summaries left in this record, or */
/*     there aren't. If there are, just incrementing the pointer */
/*     is sufficient. If there aren't, we have to find the next */
/*     record and point to the first array there. (If that */
/*     record is empty, or doesn't exist, then there are simply */
/*     no more arrays to be found.) */

    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", i__1, 
	    "daffa_", (ftnlen)1548)] = stcurr[(i__2 = p - 1) < 1000 && 0 <= 
	    i__2 ? i__2 : s_rnge("stcurr", i__2, "daffa_", (ftnlen)1548)] + 1;
    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)1550)] > stnseg[(i__2 = p - 1) < 1000 && 
	    0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (ftnlen)1550)]
	    ) {
	if (stnext[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnext"
		, i__1, "daffa_", (ftnlen)1552)] == 0) {

/*           There are no more arrays in the list. */

	    *found = FALSE_;

/*           Make sure that the array pointer stays pointing to */
/*           the position following the end of the list.  Otherwise, */
/*           a call to DAFFPA might fail to find the last array in */
/*           the list. */

	    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr"
		    , i__1, "daffa_", (ftnlen)1563)] = stnseg[(i__2 = p - 1) <
		     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daff"
		    "a_", (ftnlen)1563)] + 1;

/*           The careful reader may note that we're not updating any */
/*           of the pointers */

/*              STTHIS */
/*              STNEXT */
/*              STPREV */

/*           These will not be accessed if there is no current array. */
/*           If the array pointer is backed up again by a call to */
/*           DAFFPA, the values we have right now will be correct. */

	} else {
	    dafgsr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stfh", i__1, "daffa_", (ftnlen)1578)], &stnext[(i__2 = p 
		    - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge("stnext", i__2, 
		    "daffa_", (ftnlen)1578)], &c__1, &c__128, &stsr[(i__3 = (
		    p << 7) - 128) < 128000 && 0 <= i__3 ? i__3 : s_rnge(
		    "stsr", i__3, "daffa_", (ftnlen)1578)], &fnd);

/*           The name (character) record we've saved no longer applies */
/*           to the current summary record.  However, we've just updated */
/*           the summary record, so the summary record remains valid. */

	    sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr"
		    , i__1, "daffa_", (ftnlen)1584)] = FALSE_;
	    stthis[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stthis"
		    , i__1, "daffa_", (ftnlen)1586)] = stnext[(i__2 = p - 1) <
		     1000 && 0 <= i__2 ? i__2 : s_rnge("stnext", i__2, "daff"
		    "a_", (ftnlen)1586)];
	    stnext[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnext"
		    , i__1, "daffa_", (ftnlen)1587)] = (integer) stsr[(i__2 = 
		    (p << 7) - 128) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)1587)];
	    stprev[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stprev"
		    , i__1, "daffa_", (ftnlen)1588)] = (integer) stsr[(i__2 = 
		    (p << 7) - 127) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)1588)];
	    stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnseg"
		    , i__1, "daffa_", (ftnlen)1589)] = (integer) stsr[(i__2 = 
		    (p << 7) - 126) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)1589)];
	    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr"
		    , i__1, "daffa_", (ftnlen)1590)] = 1;
	    *found = stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
		    s_rnge("stnseg", i__1, "daffa_", (ftnlen)1592)] > 0;
	}
    }
    chkout_("DAFFNA", (ftnlen)6);
    return 0;
/* $Procedure DAFBBS ( DAF, begin backward search ) */

L_dafbbs:
/* $ Abstract */

/*     Begin a backward search for arrays in a DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     INTEGER               HANDLE */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   Handle of DAF to be searched. */

/* $ Detailed_Input */

/*     HANDLE      is the handle of a DAF on which a backward */
/*                 search is to be conducted. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     See argument HANDLE. */

/* $ Exceptions */

/*     1)  If the input handle is invalid, the error will be diagnosed */
/*         by routines called by this routine. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     begin daf backward search */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now makes the DAF designated by HANDLE the */
/*        current DAF---the one at the head of the active list.  All */
/*        saved state variables used by this routine are now part of the */
/*        state table, or its associated set of pointers. */

/*        Also, the $Exceptions section was filled out. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFBBS", (ftnlen)6);
    }

/*     Check out the file handle before going any further. */

    dafsih_(handle, "READ", (ftnlen)4);
    if (failed_()) {
	chkout_("DAFBBS", (ftnlen)6);
	return 0;
    }

/*     Initialize the state table pool, if this hasn't been done yet. */
/*     Also initialize the cell used to obtain the set of handles of */
/*     open DAFs. */

    if (first) {
	ssizei_(&c__1000, opnset);
	for (i__ = 1; i__ <= 999; ++i__) {
	    stpool[(i__1 = i__ - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stp"
		    "ool", i__1, "daffa_", (ftnlen)1774)] = i__ + 1;
	}
	stpool[999] = -1;
	stfptr = 1;
	first = FALSE_;
    }

/*     See whether we already have an entry for this DAF in the */
/*     state table.  Find the previous node if possible. */

    p = sthead;
    prev = -1;
    fnd = FALSE_;
    while(p != -1 && ! fnd) {
	if (stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
		i__1, "daffa_", (ftnlen)1793)] == *handle) {
	    fnd = TRUE_;
	} else {
	    prev = p;
	    p = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stpool", i__1, "daffa_", (ftnlen)1797)];
	}
    }

/*     At this point, either FND is false, or P points to a */
/*     state table entry describing the DAF indicated by HANDLE. */
/*     In the latter case, PREV is the predecessor of P. */

    if (fnd) {

/*        We already have a dossier on this DAF.  We already have */
/*        the information on the summary format, but we must re-set */
/*        our summary record pointers and our name record availability */
/*        flag. */

/*        Rather than doing the update here, we do it outside of this */
/*        IF block.  That way, the update gets done in just one place. */
/*        This just makes life easier:  if the collection of state */
/*        variables is changed, there are fewer places to forget to */
/*        make the required code changes. */

/*        Move the node for this DAF to the head of the active list, */
/*        if it is not already there: */

/*           - Make the predecessor of P point to the successor of P. */

/*           - Make P point to the head of the active list. */

/*           - Make P the active list head node. */


	if (p != sthead) {

/*           P is in the active list, but is not at the head.  So, */
/*           the predecessor of P is not NIL. */

	    stpool[(i__1 = prev - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stpool", i__1, "daffa_", (ftnlen)1835)] = stpool[(i__2 = 
		    p - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge("stpool", i__2,
		     "daffa_", (ftnlen)1835)];
	    stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool"
		    , i__1, "daffa_", (ftnlen)1836)] = sthead;
	    sthead = p;
	}
    } else {

/*        We don't yet have any information on this DAF.  Make a new */
/*        state table entry for the DAF.  We may need to make room for */
/*        the new information by freeing space allocated to DAFs that */
/*        are no longer open. */

	if (stfptr == -1) {

/*           Oops, we're out of space.  Time for garbage collection. */
/*           Test each file handle to see whether it designates a DAF */
/*           that is still open.  DAFHOF will tell us which handles */
/*           point to open DAFs. */

	    dafhof_(opnset);
	    p = sthead;
	    prev = -1;

/*           For every DAF file represented in the state table, we'll */
/*           delete the corresponding state information if the DAF is */
/*           now closed.  We traverse the active list, examining each */
/*           file handle as we go. */

	    while(p != -1) {
		if (elemi_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
			s_rnge("stfh", i__1, "daffa_", (ftnlen)1868)], opnset)
			) {

/*                 The file is open. Have a look at the next node. */

		    prev = p;
		    p = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
			    s_rnge("stpool", i__1, "daffa_", (ftnlen)1873)];
		} else {

/*                 This file handle is not on the list, so free the */
/*                 node pointing to the information about the DAF it */
/*                 designated: */

/*                    - Save the successor of P. */

/*                    - Link the predecessor of node P to the successor */
/*                      of P, if the predecessor is not NIL. */

/*                    - If it happens that P is the head node of the */
/*                      active list, set the head equal to the */
/*                      successor of P. */

/*                    - Link P into the free list. */

/*                    - Set P equal to its saved successor. */

/*                    - (PREV remains unchanged.) */


		    nextp = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 :
			     s_rnge("stpool", i__1, "daffa_", (ftnlen)1897)];
		    if (p == sthead) {

/*                    Re-assign STHEAD so that we don't lose the head */
/*                    of the active list.  P has no predecessor in this */
/*                    case, so there's no need to set the forward pointer */
/*                    of node PREV. */

			sthead = nextp;
		    } else {

/*                    Since P is not the head node of the active list, */
/*                    PREV is not NIL, so we'll need to set the forward */
/*                    pointer of node PREV. */

			stpool[(i__1 = prev - 1) < 1000 && 0 <= i__1 ? i__1 : 
				s_rnge("stpool", i__1, "daffa_", (ftnlen)1915)
				] = nextp;
		    }
		    stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
			    "stpool", i__1, "daffa_", (ftnlen)1920)] = stfptr;
		    stfptr = p;
		    p = nextp;
		}
	    }

/*           At this point, we've freed all nodes from the active */
/*           list that were used to index information about DAFs that */
/*           are no longer open.  If there's any more room in the state */
/*           table, we have it now. */

	}

/*        If there still is no room, there is a bug in DAFAH, since DAFAH */
/*        should not allow more than TBSIZE DAFs to be open.  So, we */
/*        assume that we've found some room.  The first free node is */
/*        indicated by STFPTR.  We'll allocate this node and use it to */
/*        index the state information for the new DAF. */

	p = stfptr;

/*        Update the free list pointer, link P to the previous head */
/*        of the active list, and make P the head of the active list. */

	stfptr = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stpool", i__1, "daffa_", (ftnlen)1947)];
	stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool", 
		i__1, "daffa_", (ftnlen)1948)] = sthead;
	sthead = p;
    }

/*     At this point, P is the head node of the active list, and P is */
/*     the index in the state table of the information for the current */
/*     DAF. */


/*     Read the file record and last summary record. Do not read the */
/*     corresponding name record until necessary. In most searches, */
/*     names are of no interest. */

    dafrfr_(handle, &nd, &ni, ifname, &fward, &bward, &free, (ftnlen)60);
    dafgsr_(handle, &bward, &c__1, &c__128, &stsr[(i__1 = (p << 7) - 128) < 
	    128000 && 0 <= i__1 ? i__1 : s_rnge("stsr", i__1, "daffa_", (
	    ftnlen)1965)], &fnd);
    stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", i__1, 
	    "daffa_", (ftnlen)1967)] = *handle;
    stthis[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stthis", i__1, 
	    "daffa_", (ftnlen)1968)] = bward;
    stnext[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnext", i__1, 
	    "daffa_", (ftnlen)1969)] = (integer) stsr[(i__2 = (p << 7) - 128) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1969)];
    stprev[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stprev", i__1, 
	    "daffa_", (ftnlen)1970)] = (integer) stsr[(i__2 = (p << 7) - 127) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1970)];
    stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnseg", i__1, 
	    "daffa_", (ftnlen)1971)] = (integer) stsr[(i__2 = (p << 7) - 126) 
	    < 128000 && 0 <= i__2 ? i__2 : s_rnge("stsr", i__2, "daffa_", (
	    ftnlen)1971)];
    sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", i__1, 
	    "daffa_", (ftnlen)1972)] = FALSE_;

/*     The arrays are returned in backward order from each summary */
/*     record. */

    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", i__1, 
	    "daffa_", (ftnlen)1978)] = stnseg[(i__2 = p - 1) < 1000 && 0 <= 
	    i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (ftnlen)1978)] + 1;
    chkout_("DAFBBS", (ftnlen)6);
    return 0;
/* $Procedure DAFFPA ( DAF, find previous array ) */

L_daffpa:
/* $ Abstract */

/*     Find the previous (backward) array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     LOGICAL               FOUND */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     FOUND      O   True if an array was found. */

/* $ Detailed_Input */

/*     None. */

/* $ Detailed_Output */

/*     FOUND       is true if an array was found, and is false if, */
/*                 when this routine is called, the current array is */
/*                 the head of the array list.  (Recall that the */
/*                 arrays in a DAF may be viewed as a doubly linked */
/*                 list, with the head being the first array in the */
/*                 file.) */


/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1) If this routine is called before a search is begun, the */
/*        error SPICE(DAFNOSEARCH) is signalled. */

/*     2) If the DAF to be searched has actually been closed, the error */
/*        will be diagnosed by routines called by this routine. */

/*     3) If the beginning of the array list has already been reached */
/*        when this routine is called, this routine will not change the */
/*        current array.  FOUND will be false on output. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */
/*        Also, a bug fix was made to the array pointer adjustment */
/*        algorithm. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     find previous daf array */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        Also, a bug fix was made to the array pointer adjustment */
/*        algorithm:  the pointer is no longer decremented if it */
/*        is already less than 1 and the array summary pointer */
/*        is already pointing to the first array summary.  In */
/*        addition, a test made to detect this condition was fixed: */
/*        the test */

/*           CURR .EQ. 0 */

/*        was replaced by */

/*           STCURR(P) .LE. 0 */

/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFFPA", (ftnlen)6);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     FOUND will be false until we make it past the error checks. */

    *found = FALSE_;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFFPA", (ftnlen)6);
	return 0;

/*     Make sure that the `current' DAF is still open. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2189)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFFPA", (ftnlen)6);
	    return 0;
	}
    }

/*     Now that we know a search is going on, assume that we will find */
/*     an array until proven otherwise. */

    *found = TRUE_;

/*     Either there are more summaries left in this record, or */
/*     there aren't. If there are, just decrementing the pointer */
/*     is sufficient. If there aren't, we have to find the previous */
/*     record and point to the last array there. (If that */
/*     record is empty, or doesn't exist, then there are simply */
/*     no more arrays to be found.) */

    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", i__1, 
	    "daffa_", (ftnlen)2212)] = stcurr[(i__2 = p - 1) < 1000 && 0 <= 
	    i__2 ? i__2 : s_rnge("stcurr", i__2, "daffa_", (ftnlen)2212)] - 1;
    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)2214)] <= 0) {
	if (stprev[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stprev"
		, i__1, "daffa_", (ftnlen)2216)] == 0) {

/*           There is no predecessor of the current array in the list. */

	    *found = FALSE_;

/*           Make sure that the array pointer stays pointing to */
/*           the position preceding the front of the list.  Otherwise, */
/*           a call to DAFFNA might fail to find the first array in */
/*           the list. */

	    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr"
		    , i__1, "daffa_", (ftnlen)2227)] = 0;

/*           The careful reader may note that we're not updating any */
/*           of the pointers */

/*              STTHIS */
/*              STNEXT */
/*              STPREV */

/*           These will not be accessed if there is no current array. */
/*           If the array pointer is moved forward again by a call to */
/*           DAFFNA, the values we have right now will be correct. */

	} else {
	    dafgsr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stfh", i__1, "daffa_", (ftnlen)2242)], &stprev[(i__2 = p 
		    - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge("stprev", i__2, 
		    "daffa_", (ftnlen)2242)], &c__1, &c__128, &stsr[(i__3 = (
		    p << 7) - 128) < 128000 && 0 <= i__3 ? i__3 : s_rnge(
		    "stsr", i__3, "daffa_", (ftnlen)2242)], &fnd);

/*           The name (character) record we've saved no longer applies */
/*           to the current summary record.  However, we've just updated */
/*           the summary record, so the summary record remains valid. */

	    sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr"
		    , i__1, "daffa_", (ftnlen)2248)] = FALSE_;
	    stthis[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stthis"
		    , i__1, "daffa_", (ftnlen)2250)] = stprev[(i__2 = p - 1) <
		     1000 && 0 <= i__2 ? i__2 : s_rnge("stprev", i__2, "daff"
		    "a_", (ftnlen)2250)];
	    stnext[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnext"
		    , i__1, "daffa_", (ftnlen)2251)] = (integer) stsr[(i__2 = 
		    (p << 7) - 128) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)2251)];
	    stprev[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stprev"
		    , i__1, "daffa_", (ftnlen)2252)] = (integer) stsr[(i__2 = 
		    (p << 7) - 127) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)2252)];
	    stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnseg"
		    , i__1, "daffa_", (ftnlen)2253)] = (integer) stsr[(i__2 = 
		    (p << 7) - 126) < 128000 && 0 <= i__2 ? i__2 : s_rnge(
		    "stsr", i__2, "daffa_", (ftnlen)2253)];
	    stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr"
		    , i__1, "daffa_", (ftnlen)2254)] = stnseg[(i__2 = p - 1) <
		     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daff"
		    "a_", (ftnlen)2254)];
	    *found = stnseg[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
		    s_rnge("stnseg", i__1, "daffa_", (ftnlen)2256)] > 0;
	}
    }
    chkout_("DAFFPA", (ftnlen)6);
    return 0;
/* $Procedure DAFGS ( DAF, get summary ) */

L_dafgs:
/* $ Abstract */

/*     Return (get) the summary for the current array in the current */
/*     DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     DOUBLE PRECISION      SUM    ( * ) */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     SUM        O   Summary for current array. */

/* $ Detailed_Input */

/*     None. */

/* $ Detailed_Output */

/*     SUM         is the summary for the current array (the array */
/*                 found by the latest call to DAFFNA or DAFFPA). */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF for which the `current' array's summary is to be */
/*         returned has actually been closed, the error will be diagnosed */
/*         by routines called by this routine. */

/*     3)  If no array is current in the current DAF, the error */
/*         SPICE(NOCURRENTARRAY) is signalled.  There is no current */
/*         array when a search is started by DAFBFS or DAFBBS, but no */
/*         calls to DAFFNA or DAFBNA have been made yet, or whenever */
/*         DAFFNA or DAFFPA return the value .FALSE. in the FOUND */
/*         argument. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */
/*        Bug fix made to handle case of having no current array. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     get daf summary */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        In addition, this routine now checks whether an array */
/*        is current before trying to read its summary.  The routine */
/*        previously crashed under these conditions. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFGS", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFGS", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2454)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFGS", (ftnlen)5);
	    return 0;
	}
    }

/*     Check the current pointer position to make sure that it's in */
/*     bounds.  If there is no current array, then we cannot return */
/*     a summary.  This situation occurs if DAFFNA was called when the */
/*     current array was the last, or if DAFFPA was called when the */
/*     current array was the first. */

    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)2470)] == 0) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2472)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `next' array is the first array of"
		" DAF #", (ftnlen)65);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFGS", (ftnlen)5);
	return 0;
    } else if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
	    "stcurr", i__1, "daffa_", (ftnlen)2480)] > stnseg[(i__2 = p - 1) <
	     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (
	    ftnlen)2480)]) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2482)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `previous' array is the last array"
		" of DAF #", (ftnlen)68);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFGS", (ftnlen)5);
	return 0;
    }

/*     The location of the summary depends on the current pointer */
/*     position. */

    dafhsf_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)2496)], &nd, &ni);
    sumsiz = nd + (ni + 1) / 2;
    offset = (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stc"
	    "urr", i__1, "daffa_", (ftnlen)2500)] - 1) * sumsiz + 3;
    moved_(&stsr[(i__1 = offset + 1 + (p << 7) - 129) < 128000 && 0 <= i__1 ? 
	    i__1 : s_rnge("stsr", i__1, "daffa_", (ftnlen)2502)], &sumsiz, 
	    sum);
    chkout_("DAFGS", (ftnlen)5);
    return 0;
/* $Procedure DAFGN ( DAF, get array name ) */

L_dafgn:
/* $ Abstract */

/*     Return (get) the name for the current array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     CHARACTER*(*)         NAME */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     NAME       O   Name of current array. */

/* $ Detailed_Input */

/*     None. */

/* $ Detailed_Output */

/*     NAME        is the name for the current array (the array */
/*                 found by the latest call to DAFFNA or DAFFPA). */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF for which the `current' array's name is to be */
/*         returned has actually been closed, the error will be diagnosed */
/*         by routines called by this routine. */

/*     3)  If no array is current in the current DAF, the error */
/*         SPICE(NOCURRENTARRAY) is signalled.  There is no current */
/*         array when a search is started by DAFBFS or DAFBBS, but no */
/*         calls to DAFFNA or DAFBNA have been made yet, or whenever */
/*         DAFFNA or DAFFPA return the value .FALSE. in the FOUND */
/*         argument. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */
/*        Bug fix made to handle case of having no current array. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     get daf array name */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        In addition, this routine now checks whether an array */
/*        is current before trying to read its summary.  The routine */
/*        previously crashed under these conditions. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFGN", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFGN", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2692)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFGN", (ftnlen)5);
	    return 0;
	}
    }

/*     Check the current pointer position to make sure that it's in */
/*     bounds.  If there is no current array, then we cannot get the */
/*     array's summary's name.  This situation occurs if DAFFNA was */
/*     called when the current array was the last, or if DAFFPA was */
/*     called when the current array was the first. */

    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)2708)] == 0) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2710)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `next' array is the first array of"
		" DAF #", (ftnlen)65);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFGN", (ftnlen)5);
	return 0;
    } else if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
	    "stcurr", i__1, "daffa_", (ftnlen)2718)] > stnseg[(i__2 = p - 1) <
	     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (
	    ftnlen)2718)]) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2720)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `previous' array is the last array"
		" of DAF #", (ftnlen)68);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFGN", (ftnlen)5);
	return 0;
    }

/*     Read the name record for this summary record, if we don't have it */
/*     already. */

    if (! sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", 
	    i__1, "daffa_", (ftnlen)2735)]) {
	i__4 = stthis[(i__2 = p - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge(
		"stthis", i__2, "daffa_", (ftnlen)2737)] + 1;
	dafrcr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2737)], &i__4, stnr + ((i__3 =
		 p - 1) < 1000 && 0 <= i__3 ? i__3 : s_rnge("stnr", i__3, 
		"daffa_", (ftnlen)2737)) * 1000, (ftnlen)1000);
	sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", 
		i__1, "daffa_", (ftnlen)2739)] = TRUE_;
    }

/*     The location of the name depends on the current pointer */
/*     position. */

    dafhsf_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)2748)], &nd, &ni);
    sumsiz = nd + (ni + 1) / 2;
    namsiz = sumsiz << 3;
    offset = (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stc"
	    "urr", i__1, "daffa_", (ftnlen)2754)] - 1) * namsiz;
    i__2 = offset;
    s_copy(name__, stnr + (((i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : 
	    s_rnge("stnr", i__1, "daffa_", (ftnlen)2756)) * 1000 + i__2), 
	    name_len, offset + namsiz - i__2);
    chkout_("DAFGN", (ftnlen)5);
    return 0;
/* $Procedure DAFGH ( DAF, get handle ) */

L_dafgh:
/* $ Abstract */

/*     Return (get) the handle of the DAF currently being searched. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     INTEGER               HANDLE */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     O   Handle for current DAF. */

/* $ Detailed_Input */

/*     None. */

/* $ Detailed_Output */

/*     HANDLE      is the handle for the current DAF (the handle */
/*                 connected to the DAF that is currently being */
/*                 searched). */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF whose handle is to be returned has actually been */
/*         closed, the error will be diagnosed by routines called by */
/*         this routine. */

/* $ Particulars */

/*     Under rare circumstances, it may be necessary to identify */
/*     the particular DAF that is being searched (such as when */
/*     the search is begun by one module and continued by another). */

/* $ Examples */

/*     Consider a program like the following, which examines the */
/*     individual arrays in a DAF and examines the contents of those */
/*     meeting certain criteria. */

/*        CALL DAFOPW ( FNAME, HANDLE ) */
/*        CALL DAFBFS ( HANDLE ) */
/*        CALL DAFFNA ( FOUND  ) */

/*        DO WHILE ( FOUND ) */
/*           CALL CHECK_DAF ( STATUS ) */

/*           IF ( STATUS .EQ. 'EXAMINE' ) THEN */
/*              CALL EXAMINE_DAF */
/*           END IF */

/*           CALL DAFFNA ( FOUND ) */
/*        END DO */

/*     The subroutine CHECK_DAF, which assumes that a search is in */
/*     progress, gets the summary and name for the current array, and */
/*     uses them to decide whether the data in the array merit further */
/*     consideration. */

/*        SUBROUTINE CHECK_DAF ( STATUS ) */

/*        CALL DAFGS ( SUM ) */
/*        CALL DAFGN ( NAME ) */
/*        CALL DAFUS ( SUM, ND, NI, DC, IC ) */
/*         . */
/*         . */

/*     The subroutine EXAMINE_DAF needs to examine the data in */
/*     the array itself. In order to do do, it needs to have access */
/*     not only to the summary, but to the handle of the file */
/*     containing the array. This is provided by DAFGH. */

/*        SUBROUTINE EXAMINE_DAF */

/*        CALL DAFGS ( SUM  ) */
/*        CALL DAFGH ( HANDLE ) */
/*        CALL DAFUS ( SUM, ND, NI, DC, IC ) */

/*        CALL DAFRDA ( HANDLE, BEGIN, END, DATA ) */
/*         . */
/*         . */


/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     get daf handle */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFGH", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFGH", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)2983)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFGH", (ftnlen)5);
	    return 0;
	}
    }
    *handle = stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)2993)];
    chkout_("DAFGH", (ftnlen)5);
    return 0;
/* $Procedure DAFRS ( DAF, replace summary ) */

L_dafrs:
/* $ Abstract */

/*     Change the summary for the current array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     DOUBLE PRECISION      SUM */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     SUM        I   New summary for current array. */

/* $ Detailed_Input */

/*     SUM         is the new summary for the current array. This */
/*                 replaces the existing summary. However, the addresses */
/*                 (the final two integer components) of the original */
/*                 summary are not changed. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF containing the `current' array has actually been */
/*         closed, the error will be diagnosed by routines called by */
/*         this routine. */

/*     3)  If the DAF containing the `current' array is not open for */
/*         writing, the error will be diagnosed by routines called by */
/*         this routine. */

/*     4)  If no array is current in the current DAF, the error */
/*         SPICE(NOCURRENTARRAY) is signalled.  There is no current */
/*         array when a search is started by DAFBFS or DAFBBS, but no */
/*         calls to DAFFNA or DAFBNA have been made yet, or whenever */
/*         DAFFNA or DAFFPA return the value .FALSE. in the FOUND */
/*         argument. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */
/*        Bug fix made to handle case of having no current array. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     replace daf summary */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        In addition, this routine now checks whether an array */
/*        is current before trying to read its summary.  The routine */
/*        previously crashed under these conditions. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFRS", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFRS", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open, and that it */
/*     is open for writing. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3192)], "WRITE", (ftnlen)5);
	if (failed_()) {
	    chkout_("DAFRS", (ftnlen)5);
	    return 0;
	}
    }

/*     Check the current pointer position to make sure that it's in */
/*     bounds.  If there is no current array, then we cannot replace the */
/*     array's  summary.  This situation occurs if DAFFNA was called */
/*     when the current array was the last, or if DAFFPA was called when */
/*     the current array was the first. */

    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)3208)] == 0) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3210)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `next' array is the first array of"
		" DAF #", (ftnlen)65);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFRS", (ftnlen)5);
	return 0;
    } else if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
	    "stcurr", i__1, "daffa_", (ftnlen)3218)] > stnseg[(i__2 = p - 1) <
	     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (
	    ftnlen)3218)]) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3220)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `previous' array is the last array"
		" of DAF #", (ftnlen)68);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFRS", (ftnlen)5);
	return 0;
    }

/*     The location of the summary depends on the current pointer */
/*     position. */

    dafhsf_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3234)], &nd, &ni);
    sumsiz = nd + (ni + 1) / 2;
    offset = (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stc"
	    "urr", i__1, "daffa_", (ftnlen)3238)] - 1) * sumsiz + 3;

/*     Get the existing summary, and unpack it. Replace everything */
/*     but the addresses (the final two integer components), and */
/*     repack. Then replace the existing summary within the record. */

    moved_(&stsr[(i__1 = offset + 1 + (p << 7) - 129) < 128000 && 0 <= i__1 ? 
	    i__1 : s_rnge("stsr", i__1, "daffa_", (ftnlen)3245)], &sumsiz, 
	    exsum);
    dafus_(exsum, &nd, &ni, exdc, exic);
    dafus_(sum, &nd, &ni, newdc, newic);
    moved_(newdc, &nd, exdc);
    i__1 = ni - 2;
    movei_(newic, &i__1, exic);
    dafps_(&nd, &ni, exdc, exic, exsum);
    moved_(exsum, &sumsiz, &stsr[(i__1 = offset + 1 + (p << 7) - 129) < 
	    128000 && 0 <= i__1 ? i__1 : s_rnge("stsr", i__1, "daffa_", (
	    ftnlen)3254)]);

/*     Rewrite the modified summary record. */

    dafwdr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3259)], &stthis[(i__2 = p - 1) < 1000 && 
	    0 <= i__2 ? i__2 : s_rnge("stthis", i__2, "daffa_", (ftnlen)3259)]
	    , &stsr[(i__3 = (p << 7) - 128) < 128000 && 0 <= i__3 ? i__3 : 
	    s_rnge("stsr", i__3, "daffa_", (ftnlen)3259)]);
    chkout_("DAFRS", (ftnlen)5);
    return 0;
/* $Procedure DAFRN ( DAF, change array name ) */

L_dafrn:
/* $ Abstract */

/*     Replace the name for the current array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     CHARACTER*(*)         NAME */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     NAME       I   New name for current array. */

/* $ Detailed_Input */

/*     NAME        is the new name for the current array. */
/*                 This replaces the existing name. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     None. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF containing the `current' array has actually been */
/*         closed, the error will be diagnosed by routines called by */
/*         this routine. */

/*     3)  If the DAF containing the `current' array is not open for */
/*         writing, the error will be diagnosed by routines called by */
/*         this routine. */

/*     4)  If no array is current in the current DAF, the error */
/*         SPICE(NOCURRENTARRAY) is signalled.  There is no current */
/*         array when a search is started by DAFBFS or DAFBBS, but no */
/*         calls to DAFFNA or DAFBNA have been made yet, or whenever */
/*         DAFFNA or DAFFPA return the value .FALSE. in the FOUND */
/*         argument. */

/* $ Particulars */

/*     See DAFFA. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     NAIF Document 167.0, "Double Precision Array Files (DAF) */
/*     Specification and User's Guide" */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/* -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN) */

/*        Literature references added to the header. */

/* -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) */

/* -& */
/* $ Index_Entries */

/*     change daf array name */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        In addition, this routine now checks whether an array */
/*        is current before trying to read its summary.  The routine */
/*        previously crashed under these conditions. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFRN", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFRN", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open, and that it */
/*     is open for writing. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3453)], "WRITE", (ftnlen)5);
	if (failed_()) {
	    chkout_("DAFRN", (ftnlen)5);
	    return 0;
	}
    }

/*     Check the current pointer position to make sure that it's in */
/*     bounds.  If there is no current array, then we cannot replace */
/*     the array's summary's name.  This situation occurs if DAFFNA was */
/*     called when the current array was the last, or if DAFFPA was */
/*     called when the current array was the first. */

    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)3469)] == 0) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3471)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `next' array is the first array of"
		" DAF #", (ftnlen)65);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFRN", (ftnlen)5);
	return 0;
    } else if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
	    "stcurr", i__1, "daffa_", (ftnlen)3479)] > stnseg[(i__2 = p - 1) <
	     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (
	    ftnlen)3479)]) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3481)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `previous' array is the last array"
		" of DAF #", (ftnlen)68);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFRN", (ftnlen)5);
	return 0;
    }

/*     Read the name record for this summary record, if we don't have it */
/*     already. */

    if (! sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", 
	    i__1, "daffa_", (ftnlen)3497)]) {
	i__4 = stthis[(i__2 = p - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge(
		"stthis", i__2, "daffa_", (ftnlen)3499)] + 1;
	dafrcr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3499)], &i__4, stnr + ((i__3 =
		 p - 1) < 1000 && 0 <= i__3 ? i__3 : s_rnge("stnr", i__3, 
		"daffa_", (ftnlen)3499)) * 1000, (ftnlen)1000);
	sthvnr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sthvnr", 
		i__1, "daffa_", (ftnlen)3501)] = TRUE_;
    }

/*     The location of the name depends on the current pointer */
/*     position. */

    dafhsf_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3510)], &nd, &ni);
    sumsiz = nd + (ni + 1) / 2;
    namsiz = sumsiz << 3;
    offset = (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stc"
	    "urr", i__1, "daffa_", (ftnlen)3516)] - 1) * namsiz;
    i__2 = offset;
    s_copy(stnr + (((i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stnr",
	     i__1, "daffa_", (ftnlen)3518)) * 1000 + i__2), name__, offset + 
	    namsiz - i__2, name_len);

/*     Rewrite the character record. */

    i__4 = stthis[(i__2 = p - 1) < 1000 && 0 <= i__2 ? i__2 : s_rnge("stthis",
	     i__2, "daffa_", (ftnlen)3523)] + 1;
    dafwcr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3523)], &i__4, stnr + ((i__3 = p - 1) < 
	    1000 && 0 <= i__3 ? i__3 : s_rnge("stnr", i__3, "daffa_", (ftnlen)
	    3523)) * 1000, (ftnlen)1000);
    chkout_("DAFRN", (ftnlen)5);
    return 0;
/* $Procedure DAFWS ( DAF, write summary ) */

L_dafws:
/* $ Abstract */

/*     Write a new summary for the current array in the current DAF. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     DOUBLE PRECISION      SUM ( * ) */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     SUM        I   New summary for current array in the current DAF. */

/* $ Detailed_Input */

/*     SUM         is the new summary for the current array. This */
/*                 replaces the existing summary, including the */
/*                 addresses (the final two integer components) of */
/*                 the original summary. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Files */

/*     DAFWS updates the DAF currently being searched.  The handle */
/*     of this DAF can be retrieved using the routine DAFGH. */

/* $ Exceptions */

/*     1)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/*     2)  If the DAF containing the `current' array has actually been */
/*         closed, the error will be diagnosed by routines called by */
/*         this routine. */

/*     3)  If the DAF containing the `current' array is not open for */
/*         writing, the error will be diagnosed by routines called by */
/*         this routine. */

/*     4)  If no array is current in the current DAF, the error */
/*         SPICE(NOCURRENTARRAY) is signalled.  There is no current */
/*         array when a search is started by DAFBFS or DAFBBS, but no */
/*         calls to DAFFNA or DAFBNA have been made yet, or whenever */
/*         DAFFNA or DAFFPA return the value .FALSE. in the FOUND */
/*         argument. */

/* $ Particulars */

/*     Unless you are reordering the arrays in the file being searched, */
/*     you should be using DAFRS instead of this routine. */

/*     See also DAFFA, DAFRS. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman    (JPL) */
/*     W.L. Taber      (JPL) */
/*     I.M. Underwood  (JPL) */

/* $ Version */

/* -    SPICELIB Version 2.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */
/*        Bug fix made to handle case of having no current array. */

/* -    SPICELIB Version 1.0.0, 28-MAR-1991 (IMU) */

/* -& */
/* $ Index_Entries */

/*     write daf summary */

/* -& */
/* $ Revisions */

/* -    SPICELIB Version 2.0.0, 04-SEP-1991 (NJB) (WLT) */

/*        Updated to support simultaneous searches of multiple DAFs. */

/*        This routine now operates on the current DAF---the one at */
/*        the head of the active list.  All saved state variables */
/*        used by this routine are now part of the state table, or */
/*        its associated set of pointers. */

/*        In addition, this routine now checks whether an array */
/*        is current before trying to read its summary.  The routine */
/*        previously crashed under these conditions. */
/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFWS", (ftnlen)5);
    }

/*     Operate on the last DAF in which a search has been started. */

    p = sthead;

/*     Make sure that a search has been started in this DAF. */

    if (p == -1) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFWS", (ftnlen)5);
	return 0;

/*     Make sure that the `current' DAF is still open, and that it is */
/*     open for writing. */

    } else {
	dafsih_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3719)], "READ", (ftnlen)4);
	if (failed_()) {
	    chkout_("DAFWS", (ftnlen)5);
	    return 0;
	}
    }

/*     Check the current pointer position to make sure that it's in */
/*     bounds.  If there is no current array, then we cannot write a */
/*     new array summary. This situation occurs if DAFFNA was called */
/*     when the current array was the last, or if DAFFPA was called */
/*     when the current array was the first. */

    if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stcurr", 
	    i__1, "daffa_", (ftnlen)3735)] == 0) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3737)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `next' array is the first array of"
		" DAF #", (ftnlen)65);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFWS", (ftnlen)5);
	return 0;
    } else if (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
	    "stcurr", i__1, "daffa_", (ftnlen)3745)] > stnseg[(i__2 = p - 1) <
	     1000 && 0 <= i__2 ? i__2 : s_rnge("stnseg", i__2, "daffa_", (
	    ftnlen)3745)]) {
	dafhfn_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		"stfh", i__1, "daffa_", (ftnlen)3747)], dafnam, (ftnlen)255);
	setmsg_("No array is current; the `previous' array is the last array"
		" of DAF #", (ftnlen)68);
	errch_("#", dafnam, (ftnlen)1, (ftnlen)255);
	sigerr_("SPICE(NOCURRENTARRAY)", (ftnlen)21);
	chkout_("DAFWS", (ftnlen)5);
	return 0;
    }

/*     The location of the summary depends on the current pointer */
/*     position. */

    dafhsf_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3763)], &nd, &ni);
    sumsiz = nd + (ni + 1) / 2;
    offset = (stcurr[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stc"
	    "urr", i__1, "daffa_", (ftnlen)3767)] - 1) * sumsiz + 3;
    moved_(sum, &sumsiz, &stsr[(i__1 = offset + 1 + (p << 7) - 129) < 128000 
	    && 0 <= i__1 ? i__1 : s_rnge("stsr", i__1, "daffa_", (ftnlen)3769)
	    ]);

/*     Rewrite the modified summary record. */

    dafwdr_(&stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
	    i__1, "daffa_", (ftnlen)3774)], &stthis[(i__2 = p - 1) < 1000 && 
	    0 <= i__2 ? i__2 : s_rnge("stthis", i__2, "daffa_", (ftnlen)3774)]
	    , &stsr[(i__3 = (p << 7) - 128) < 128000 && 0 <= i__3 ? i__3 : 
	    s_rnge("stsr", i__3, "daffa_", (ftnlen)3774)]);
    chkout_("DAFWS", (ftnlen)5);
    return 0;
/* $Procedure DAFCS ( DAF, continue search ) */

L_dafcs:
/* $ Abstract */

/*     Select a DAF that already has a search in progress as the */
/*     one to continue searching. */

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

/*     DAF */

/* $ Keywords */

/*     FILES */

/* $ Declarations */

/*     HANDLE */

/* $ Brief_I/O */

/*     Variable  I/O  Description */
/*     --------  ---  -------------------------------------------------- */
/*     HANDLE     I   Handle of DAF to continue searching. */

/* $ Detailed_Input */

/*     HANDLE         is the handle of a DAF in which either a forward */
/*                    or backward search has already been started by */
/*                    DAFBFS or DAFBBS.  The DAF may be open for read */
/*                    or write access. */

/* $ Detailed_Output */

/*     None. */

/* $ Parameters */

/*     None. */

/* $ Exceptions */

/*     1)  If the input handle is invalid, the error will be diagnosed */
/*         by routines called by this routine. */

/*     2)  If this routine is called when no search is in progress in the */
/*         the current DAF, the error SPICE(DAFNOSEARCH) is signalled. */

/* $ Files */

/*     None. */

/* $ Particulars */

/*     DAFCS supports simultaneous searching of multiple DAFs.  In */
/*     applications that use this capability, DAFCS should be called */
/*     prior to each call to DAFFNA, DAFFPA, DAFGN, DAFGS, DAFRS, or */
/*     DAFWS, to specify which DAF is to be acted upon. */

/* $ Examples */

/*     See DAFFA. */

/* $ Restrictions */

/*     None. */

/* $ Literature_References */

/*     None. */

/* $ Author_and_Institution */

/*     N.J. Bachman   (JPL) */
/*     W.L. Taber     (JPL) */

/* $ Version */

/* -    SPICELIB Version 1.0.1, 10-MAR-1992 (WLT) */

/*        Comment section for permuted index source lines was added */
/*        following the header. */

/* -    SPICELIB Version 1.0.0, 04-SEP-1991 (NJB) (WLT) */

/* -& */
/* $ Index_Entries */

/*     select a daf to continue searching */

/* -& */

/*     Standard SPICE error handling. */

    if (return_()) {
	return 0;
    } else {
	chkin_("DAFCS", (ftnlen)5);
    }

/*     Validate the DAF's handle before going any further.  DAFSIH will */
/*     signal an error if HANDLE doesn't designate an open DAF. */

    dafsih_(handle, "READ", (ftnlen)4);
    if (failed_()) {
	chkout_("DAFCS", (ftnlen)5);
	return 0;
    }

/*     See whether we already have an entry for this DAF in the */
/*     state table.  Find the previous node if possible. */

    p = sthead;
    prev = -1;
    fnd = FALSE_;
    while(p != -1 && ! fnd) {
	if (stfh[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stfh", 
		i__1, "daffa_", (ftnlen)3938)] == *handle) {
	    fnd = TRUE_;
	} else {
	    prev = p;
	    p = stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge(
		    "stpool", i__1, "daffa_", (ftnlen)3942)];
	}
    }

/*     Either FND is false, or P is the index in the state table of */
/*     the DAF specified by HANDLE, and PREV is the predecessor of P. */


/*     You can't continue searching a DAF that you're not already */
/*     searching. */

    if (! fnd) {
	setmsg_("No DAF is currently being searched.", (ftnlen)35);
	sigerr_("SPICE(DAFNOSEARCH)", (ftnlen)18);
	chkout_("DAFCS", (ftnlen)5);
	return 0;
    }

/*     Move the node for this DAF to the head of the active list, */
/*     if it is not already there: */

/*        - Make the predecessor of P point to the successor of P. */

/*        - Make P point to the head of the active list. */

/*        - Make P the active list head node. */


    if (p != sthead) {

/*        P is in the active list, but is not at the head.  So, */
/*        the predecessor of P is not NIL. */

	stpool[(i__1 = prev - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool",
		 i__1, "daffa_", (ftnlen)3983)] = stpool[(i__2 = p - 1) < 
		1000 && 0 <= i__2 ? i__2 : s_rnge("stpool", i__2, "daffa_", (
		ftnlen)3983)];
	stpool[(i__1 = p - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("stpool", 
		i__1, "daffa_", (ftnlen)3984)] = sthead;
	sthead = p;
    }
    chkout_("DAFCS", (ftnlen)5);
    return 0;
} /* daffa_ */

/* Subroutine */ int daffa_(integer *handle, doublereal *sum, char *name__, 
	logical *found, ftnlen name_len)
{
    return daffa_0_(0, handle, sum, name__, found, name_len);
    }

/* Subroutine */ int dafbfs_(integer *handle)
{
    return daffa_0_(1, handle, (doublereal *)0, (char *)0, (logical *)0, (
	    ftnint)0);
    }

/* Subroutine */ int daffna_(logical *found)
{
    return daffa_0_(2, (integer *)0, (doublereal *)0, (char *)0, found, (
	    ftnint)0);
    }

/* Subroutine */ int dafbbs_(integer *handle)
{
    return daffa_0_(3, handle, (doublereal *)0, (char *)0, (logical *)0, (
	    ftnint)0);
    }

/* Subroutine */ int daffpa_(logical *found)
{
    return daffa_0_(4, (integer *)0, (doublereal *)0, (char *)0, found, (
	    ftnint)0);
    }

/* Subroutine */ int dafgs_(doublereal *sum)
{
    return daffa_0_(5, (integer *)0, sum, (char *)0, (logical *)0, (ftnint)0);
    }

/* Subroutine */ int dafgn_(char *name__, ftnlen name_len)
{
    return daffa_0_(6, (integer *)0, (doublereal *)0, name__, (logical *)0, 
	    name_len);
    }

/* Subroutine */ int dafgh_(integer *handle)
{
    return daffa_0_(7, handle, (doublereal *)0, (char *)0, (logical *)0, (
	    ftnint)0);
    }

/* Subroutine */ int dafrs_(doublereal *sum)
{
    return daffa_0_(8, (integer *)0, sum, (char *)0, (logical *)0, (ftnint)0);
    }

/* Subroutine */ int dafrn_(char *name__, ftnlen name_len)
{
    return daffa_0_(9, (integer *)0, (doublereal *)0, name__, (logical *)0, 
	    name_len);
    }

/* Subroutine */ int dafws_(doublereal *sum)
{
    return daffa_0_(10, (integer *)0, sum, (char *)0, (logical *)0, (ftnint)0)
	    ;
    }

/* Subroutine */ int dafcs_(integer *handle)
{
    return daffa_0_(11, handle, (doublereal *)0, (char *)0, (logical *)0, (
	    ftnint)0);
    }

