/*

-Procedure gipool_c (Get integers from the kernel pool)

-Abstract

   Return the integer value of a kernel variable from the
   kernel pool.

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

   KERNEL

-Keywords

   CONSTANTS
   FILES

*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"

   void gipool_c ( ConstSpiceChar * name,
                   SpiceInt         start,
                   SpiceInt         room,
                   SpiceInt       * n,
                   SpiceInt       * ivals,
                   SpiceBoolean   * found )

/*

-Brief_I/O

   VARIABLE  I/O  DESCRIPTION
   --------  ---  --------------------------------------------------
   name       I   Name of the variable whose value is to be returned.
   start      I   Which component to start retrieving for name
   room       I   The largest number of values to return.
   n          O   Number of values returned for name.
   ivals      O   Values associated with name.
   found      O   True if variable is in pool.

-Detailed_Input

   name       is the name of the variable whose values are to be
              returned. If the variable is not in the pool with
              numeric type, found will be SPICEFALSE.

   start      is the index of the first component of name to return.
              The index follows the C convention of being 0 based.
              If start is less than 0, it will be treated as 0.  If
              start is greater than the total number of components
              available for name, no values will be returned (n will
              be set to zero).  However, found will still be set to
              SPICETRUE

   room       is the maximum number of components that should be
              returned for this variable.  (Usually it is the amount
              of room available in the array ivals). If room is
              less than 1 the error SPICE(BADARRAYSIZE) will be
              signaled.

-Detailed_Output

   n          is the number of values associated with name that
              are returned.  It will always be less than or equal
              to room.

              If name is not in the pool with numeric type, no value
              is given to n.

   ivals      is the array of values associated with name.
              If name is not in the pool with numeric type, no
              values are given to the elements of ivals.

   found      is SPICETRUE if the variable is in the pool and has
              numeric type, SPICEFALSE if it is not.

-Parameters

   None.

-Exceptions

   1) If the value of room is less than one the error
      SPICE(BADARRAYSIZE) is signaled.

   2) If a value requested is outside the valid range
      of integers, the error SPICE(INTOUTOFRANGE) is signaled.

   3) If the input string pointer is null, the error SPICE(NULLPOINTER)
      will be signaled.

   4) If the input string has length zero, the error SPICE(EMPTYSTRING)
      will be signaled.

-Files

   None.

-Particulars

   This routine provides the user interface for retrieving
   integer data stored in the kernel pool.  This interface
   allows you to retrieve the data associated with a variable
   in multiple accesses.  Under some circumstances this alleviates
   the problem of having to know in advance the maximum amount
   of space needed to accommodate all kernel variables.

   However, this method of access does come with a price. It is
   always more efficient to retrieve all of the data associated
   with a kernel pool data in one call than it is to retrieve
   it in sections.

   See also the entry points gdpool_c and gcpool_c.

-Examples


   The following code fragment demonstrates how the data stored
   in a kernel pool variable can be retrieved in pieces.  Using the
   kernel "test.ker" which contains

   \begindata

   CTEST_VAL = ('LARRY', 'MOE', 'CURLY' )

   ITEST_VAL = ( 3141, 186, 282 )

   DTEST_VAL = ( 3.1415, 186.282, .0175 )


   The program...

   #include <stdio.h>
   #include <string.h>

   #include "SpiceUsr.h"
   #include "SpiceZmc.h"

   #define NUMVALS 2


   void main()
      {

      SpiceInt          n;
      SpiceInt          i;

      SpiceBoolean      found;

      SpiceInt          ivals[NUMVALS];


      ldpool_c ( "test.ker" );


      /. Is data available by that name. ./

      gipool_c ( "ITEST_VAL", 0, NUMVALS, &n, ivals, &found );


      /. If so, show me the values. ./

      if ( !found )
         {
         printf ( "No int data available for ITEST_VAL.\n" );
         }

      else
         {

         for ( i=0; i < NUMVALS; i++ )
            {
            gipool_c ( "ITEST_VAL", 1, NUMVALS, &n, ivals, &found );

            printf ( "%d \n", ivals[i] );
            }

         }

      exit(0);
      }


   Output should be

     186
     282


-Restrictions

   None.

-Literature_References

   None.

-Author_and_Institution

   W.L. Taber  (JPL)
   E.D. Wright (JPL)

-Version

   -CSPICE Version 2.1.0 22-JUN-1999   (EDW)

      Re-implemented routine without dynamically allocated, temporary
      strings.

      Added local variable to return boolean/logical values.  This
      fix allows the routine to function if int and long are different
      sizes.

   -CSPICE Version 2.0.1 08-FEB-1998   (EDW)

       The start parameter is now zero based as per C convention.

   -CSPICE Version 1.0.0, 6-JAN-1998    (EDW)

       Replaced example routine.  Included the data for a test kernel.

   -CSPICE Version 1.0.0, 25-OCT-1997   (EDW)

-Index_Entries

   RETURN the integer value of a pooled kernel variable

-&
*/

{ /* Begin gipool_c */

   /*
   Local variables
   */
   logical         yes;


   /* The index is zero based here but not in gipool_. */
   start = start + 1;


   /*
   Participate in error handling
   */

   chkin_c ( "gipool_c");


   /*
   Check the input string name to make sure the pointer is non-null
   and the string length is non-zero.
   */
   CHKFSTR ( CHK_STANDARD, "gipool_c", name );


   /*
   Call the f2c'd routine
   */

   gipool_( ( char    * ) name,
            ( integer * ) &start,
            ( integer * ) &room,
            ( integer * ) n,
            ( integer * ) ivals,
            ( logical * ) &yes,
            ( ftnlen    ) strlen(name) );


   /* Cast back to a SpiceBoolean. */
   *found = yes;


   chkout_c ( "gipool_c");


} /* End gipool_c */
