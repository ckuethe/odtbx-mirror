/*

-Procedure spkw12_c ( Write SPK segment, type 12 )

-Abstract
 
   Write a type 12 segment to an SPK file. 
 
-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading
 
   NAIF_IDS 
   SPC 
   SPK 
   TIME 
 
-Keywords
 
   EPHEMERIS 
   FILES 
 
*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZst.h"
   #include "SpiceZmc.h"
   #undef    spkw12_c
   

   void spkw12_c ( SpiceInt             handle,
                   SpiceInt             body,
                   SpiceInt             center, 
                   ConstSpiceChar     * frame,
                   SpiceDouble          first,
                   SpiceDouble          last,
                   ConstSpiceChar     * segid,
                   SpiceInt             degree,
                   SpiceInt             n,
                   ConstSpiceDouble     states[][6],
                   SpiceDouble          epoch0,
                   SpiceDouble          step        )

/*

-Brief_I/O
 
   Variable  I/O  Description 
   --------  ---  -------------------------------------------------- 
   handle     I   Handle of an SPK file open for writing. 
   body       I   NAIF code for an ephemeris object. 
   center     I   NAIF code for center of motion of body. 
   frame      I   Reference frame name. 
   first      I   Start time of interval covered by segment. 
   last       I   End time of interval covered by segment. 
   segid      I   Segment identifier. 
   degree     I   Degree of interpolating polynomials. 
   n          I   Number of states. 
   states     I   Array of states. 
   epoch0     I   Epoch of first state in states array.
   step       I   Time step separating epochs of states.
   MAXDEG     P   Maximum allowed degree of interpolating polynomial. 
 
-Detailed_Input
 
   handle         is the file handle of an SPK file that has been 
                  opened for writing. 
 
   body           is the NAIF integer code for an ephemeris object 
                  whose state relative to another body is described 
                  by the segment to be created. 
 
   center         is the NAIF integer code for the center of motion 
                  of the object identified by body. 
 
   frame          is the NAIF name for a reference frame 
                  relative to which the state information for body 
                  is specified. 
 
   first, 
   last           are, respectively, the start and stop times of 
                  the time interval over which the segment defines 
                  the state of body. 
 
   segid          is the segment identifier.  An SPK segment 
                  identifier may contain up to 40 characters. 
 
   degree         is the degree of the Hermite polynomials used to 
                  interpolate the states.  All components of the 
                  state vectors are interpolated by polynomials of 
                  fixed degree. 
 
   n              is the number of states in the input state vector 
                  array. 
 
   states         contains a time-ordered array of geometric states 
                  ( x, y, z, dx/dt, dy/dt, dz/dt, in kilometers and 
                  kilometers per second ) of body relative to center, 
                  specified relative to frame. 
 
   epoch0         is the epoch corresponding to the first state in
                  the state array.  Because extra states are needed
                  at the beginning and end of the segment in order
                  for the interpolation method to work, epoch0 will
                  normally precede first.

   step           is the time step separating the epochs of adjacent
                  states in the input state array.  step is specified
                  in seconds.
 
-Detailed_Output
 
   None.  See $Particulars for a description of the effect of this 
   routine. 
 
-Parameters
 
   MAXDEG         is the maximum allowed degree of the interpolating 
                  polynomial.  If the value of MAXDEG is increased, 
                  the SPICELIB routine SPKPVN must be changed 
                  accordingly.  In particular, the size of the 
                  record passed to SPKRnn and SPKEnn must be 
                  increased, and comments describing the record size 
                  must be changed. 
 
-Exceptions
 
   If any of the following exceptions occur, this routine will return 
   without creating a new segment. 
 
   1)  If frame is not a recognized name, the error 
       SPICE(INVALIDREFFRAME) is signaled. 
 
   2)  If the last non-blank character of segid occurs past index 40, 
       the error SPICE(SEGIDTOOLONG) is signaled. 
 
   3)  If segid contains any nonprintable characters, the error 
       SPICE(NONPRINTABLECHARS) is signaled. 
 
   4)  If degree is not at least 1 or is greater than MAXDEG, the 
       error SPICE(INVALIDDEGREE) is signaled. 
 
   5)  If degree is not odd, the error SPICE(INVALIDDEGREE) is  
       signaled. 
 
   6)  If the number of states n is not at least (degree+1)/2,  
       the error SPICE(TOOFEWSTATES) will be signaled. 
 
   7)  If first is greater than or equal to last then the error 
       SPICE(BADDESCRTIMES) will be signaled. 
 
   8)  If step is non-positive, the error SPICE(INVALIDSTEPSIZE) will
       be signaled.

   9)  If the first epoch epoch0 is greater than first, the error
       SPICE(BADDESCRTIMES) will be signaled.

   10) If the last epoch

         first + (n-1)*step

       is less than last, the error SPICE(BADDESCRTIMES) will be
       signaled.
 
   11) If either the input frame or segment ID string pointer is null,
       the error SPICE(NULLPOINTER) is signaled.
   
   12) If either the input frame or segment ID string is empty,
       the error SPICE(EMPTYSTRING) is signaled.
   
 
-Files
 
   A new type 12 SPK segment is written to the SPK file attached 
   to HANDLE. 
 
-Particulars
 
   This routine writes an SPK type 12 data segment to the open SPK 
   file according to the format described in the type 12 section of 
   the SPK Required Reading. The SPK file must have been opened with 
   write access. 
 
-Examples
 
   Suppose that you have states and are prepared to produce 
   a segment of type 12 in an SPK file. 
 
   The following code fragment could be used to add the new segment 
   to a previously opened SPK file attached to handle. The file must 
   have been opened with write access. 
 
      #include "SpiceUsr.h"
           .
           .
           .
        
      /.
      Create a segment identifier. 
      ./
      #define  SEGID  "MY_SAMPLE_SPK_TYPE_12_SEGMENT" 
 
        
      /.
      Write the segment. 
      ./
        
      spkw12_c ( handle,  body,    center,  frame, 
                 first,   last,    segid,   degree, 
                 n,       states,  epoch0,  step   );
 
-Restrictions
 
   None. 
 
-Literature_References
 
   None. 
 
-Author_and_Institution
 
   N.J. Bachman   (JPL) 
 
-Version
 
   -CSPICE Version 1.0.0, 23-FEB-2000 (NJB)

-Index_Entries
 
   write spk type_12 ephemeris data segment 
 
-&
*/

{ /* Begin spkw12_c */


   /*
   Participate in error tracing.
   */
   chkin_c ( "spkw12_c" );

   /*
   Check the input strings to make sure the pointers
   are non-null and the string lengths are non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "spkw12_c", frame );
   CHKFSTR ( CHK_STANDARD, "spkw12_c", segid );
 

   /*
   Write the segment. 
   */
   spkw12_ ( ( integer    * ) &handle,
             ( integer    * ) &body,
             ( integer    * ) &center,
             ( char       * ) frame,
             ( doublereal * ) &first,
             ( doublereal * ) &last,
             ( char       * ) segid,
             ( integer    * ) &degree,
             ( integer    * ) &n,
             ( doublereal * ) states,
             ( doublereal * ) &epoch0,
             ( doublereal * ) &step,
             ( ftnlen       ) strlen(frame),
             ( ftnlen       ) strlen(segid)  );


   chkout_c ( "spkw12_c" );

} /* End spkw12_c */
