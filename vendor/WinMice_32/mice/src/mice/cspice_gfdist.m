%-Abstract
%
%   CSPICE_GFDIST determines the time intervals over which a specified 
%   constraint on observer-target distance is met.
%   
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE 
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED 
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE 
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, 
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, 
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF 
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY 
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR 
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL 
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE 
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO 
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING 
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      Parameters-
%
%      All parameters described here are declared in the header file
%      SpiceGF.h. See that file for parameter values.
%
%      SPICE_GF_CNVTOL     
%
%               is the convergence tolerance used for finding endpoints of
%               the intervals comprising the result window.
%               SPICE_GF_CNVTOL is used to determine when binary searches
%               for roots should terminate: when a root is bracketed
%               within an interval of length SPICE_GF_CNVTOL, the root is
%               considered to have been found.
% 
%               The accuracy, as opposed to precision, of roots found 
%               by this routine depends on the accuracy of the input 
%               data. In most cases, the accuracy of solutions will be 
%               inferior to their precision.
%
%
%      Arguments-
%
%      target   the string scalar naming the target body.  Optionally, 
%               you may supply the integer ID code for the object as an
%               integer string.  For example both 'MOON' and '301'
%               are legitimate strings that indicate the moon is the
%               target body.
%
%      abcorr   the string scalar indicating the aberration corrections to apply
%               to the state evaluations to account for one-way light time and 
%               stellar aberration.
%   
%               This routine accepts the same aberration corrections as does 
%               the routine spkezr_c. See the header of spkezr_c for a
%               detailed description of the aberration correction options.
%               For convenience, the options are listed below:
%   
%                  'NONE'     Apply no correction.   
%   
%                  'LT'       "Reception" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%   
%                  'LT+S'     "Reception" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%   
%                  'CN'       "Reception" case:  converged
%                             Newtonian light time correction.
%   
%                  'CN+S'     "Reception" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%   
%                  'XLT'      "Transmission" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%   
%                  'XLT+S'    "Transmission" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%   
%                  'XCN'      "Transmission" case:  converged
%                             Newtonian light time correction.
%   
%                  'XCN+S'    "Transmission" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%   
%               The 'abcorr' string lacks sensitivity to case, and to embedded, 
%               leading and trailing blanks.
%
%      obsrvr   the string scalar naming the observing body. Optionally, you
%               may supply the ID code of the object as an integer
%               string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the 
%               observer is earth.
%
%      relate   the string or character scalar describing the constraint  
%               relational operator on observer-target distance. The result
%               window found  by this routine indicates the time intervals 
%               where the constraint is satisfied. Supported values of 
%               'relate' and corresponding meanings are shown below: 
% 
%                  '>'      Distance is greater than the reference 
%                           value 'refval'. 
% 
%                  '='      Distance is equal to the reference 
%                           value 'refval'. 
% 
%                  '<'      Distance is less than the reference 
%                           value 'refval'. 
% 
% 
%                  'ABSMAX'  Distance is at an absolute maximum. 
% 
%                  'ABSMIN'  Distance is at an absolute  minimum. 
%  
%                  'LOCMAX'  Distance is at a local maximum. 
%  
%                  'LOCMIN'  Distance is at a local minimum. 
% 
%               The caller may indicate that the region of interest 
%               is the set of time intervals where the quantity is 
%               within a specified distance of an absolute extremum. 
%               The argument 'adjust' (described below) is used to 
%               specify this distance. 
% 
%               Local extrema are considered to exist only in the 
%               interiors of the intervals comprising the confinement 
%               window:  a local extremum cannot exist at a boundary 
%               point of the confinement window. 
%
%               The 'relate' string lacks sensitivity to case, leading 
%               and trailing blanks.
%
%      refval   the double precision scalar reference value used together
%               with relate argument to define an equality or inequality to
%               satisfy by the observer-target distance. See the discussion 
%               of relate above for further information.
%
%               The units of 'refval' are km. 
%
%      adjust   a double precision scalar value used to modify searches for
%               absolute extrema: when relate is set to ABSMAX or ABSMIN and
%               adjust is set to a positive value, cspice_gfdist finds times 
%               when the observer-target vector coordinate is within 'adjust'
%               radians/kilometers of the specified extreme value.
%
%               For relate set to ABSMAX, the result window contains
%               time intervals when the observer-target vector coordinate has
%               values between ABSMAX - 'adjust' and ABSMAX.
%
%               For relate set to ABSMIN, the result window contains
%               time intervals when the observer-target distance has
%               values between ABSMIN and ABSMIN + 'adjust'.
%               
%               'adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     the double precision scalar time step size to use in the search.
%               'step' must be short enough for a search using this step
%               size to locate the time intervals where coordinate function
%               of the observer-target vector is monotone increasing or
%               decreasing. However, step must not be *too* short, or the search
%               will take an unreasonable amount of time.
%
%               The choice of 'step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%
%               'step' has units of seconds. 
%
%      nintvls  an integer scalar value specifying the number of intervals in  
%               the internal workspace array used by this routine. 'nintvls'
%               should be at least as large as the number of intervals
%               within the search region on which the specified observer-target
%               vector coordinate function is monotone increasing or decreasing.
%               It does no harm to pick a value of 'nintvls' larger than the
%               minimum required to execute the specified search, but if chosen 
%               too small, the search will fail.
% 
%      cnfine   a double precision SPICE window that confines the time
%               period over which the specified search is conducted.
%               cnfine may consist of a single interval or a collection
%               of intervals. 
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the Particulars section
%               below for further discussion.
%              
%               See the Examples section below for a code example 
%               that shows how to create a confinement window.
%
%   the call:
%
%      result = cspice_gfdist( target, abcorr, obsrvr, relate, refval, ... 
%                               adjust, step, nintvls, cnfine)
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window 'cnfine', on which the specified
%               constraint is satisfied.
%               
%               If the search is for local extrema, or for absolute
%               extrema with adjust set to zero, then normally each
%               interval of result will be a singleton: the left and
%               right endpoints of each interval will be identical.
% 
%               If no times within the confinement window satisfy the
%               constraint, 'result' will return with cardinality zero.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Find times during the first three months of the year 2007  
%      when the Earth-Moon distance is greater than 400000 km. 
%      Display the start and stop times of the time intervals 
%      over which this constraint is met, along with the Earth-Moon 
%      distance at each interval endpoint. 
% 
%      We expect the Earth-Moon distance to be an oscillatory function
%      with extrema roughly two weeks apart. Using a step size of one
%      day will guarantee that the GF system will find all distance
%      extrema. (Recall that a search for distance extrema is an
%      intermediate step in the GF search process.)
%
%      MAXWIN  =  1000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
% 
%      %
%      % Load kernels.
%      % 
%      cspice_furnsh( 'standard.tm' );
%   
%      %
%      % Store the time bounds of our search interval in
%      % the cnfine confinement window.
%      %
%      et = cspice_str2et( { '2007 JAN 01', '2007 APR 01'} );
%   
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      % 
%      % Search using a step size of 1 day (in units of
%      % seconds).  The reference value is 400000 km.
%      % We're not using the adjustment feature, so
%      % we set 'adjust' to zero.
%      %
%      target  = 'MOON';
%      abcorr  = 'NONE';
%      obsrvr  = 'EARTH';
%      relate  = '>';
%      refval  = 4.e5;
%      adjust  = 0.;
%      step    = 1.*cspice_spd;
%      nintvls = MAXWIN;
% 
%      result = cspice_gfdist( target, abcorr, obsrvr, ...
%                              relate, refval, adjust, ...
%                              step,   nintvls,        ...
%                              cnfine);
%
%      %
%      % List the beginning and ending times in each interval
%      % if result contains data.
%      %
%      for i=1:numel(result)/2
%
%         [left, right] = cspice_wnfetd( result, i );
%  
%         output = cspice_timout( [left,right], TIMFMT );
%
%         if( isequal( left, right) )
%
%            disp( ['Event time: ' output(1,:)] )
%
%         else
%
%            disp( ['From : ' output(1,:)] )
%            disp( ['To   : ' output(2,:)] )
%            disp( ' ')
%
%         end
%
%      end
%     
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      From : 2007-JAN-08 00:11:07.676219 (TDB)
%      To   : 2007-JAN-13 06:37:47.949222 (TDB)
% 
%      From : 2007-FEB-04 07:02:35.342186 (TDB)
%      To   : 2007-FEB-10 09:31:01.837531 (TDB)
% 
%      From : 2007-MAR-03 00:20:25.253726 (TDB)
%      To   : 2007-MAR-10 14:04:38.491933 (TDB)
% 
%      From : 2007-MAR-29 22:53:58.203278 (TDB)
%      To   : 2007-APR-01 00:01:05.185655 (TDB)
%
%-Particulars
%
%   This routine provides a simple interface for conducting searches
%   for observer-target distance events. 
%
%   This routine determines a set of one or more time intervals
%   within the confinement window for which the observer-target distance
%   between the two bodies satisfies some defined relationship.
%   The resulting set of intervals is returned as a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%   
%   The Search Process
%   ==================
%
%   Regardless of the type of constraint selected by the caller, this
%   routine starts the search for solutions by determining the time
%   periods, within the confinement window, over which the specified
%   distance function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE window.
%   Having found these windows, all of the distance function's 
%   local extrema within the confinement window are known. Absolute extrema
%   then can be found very easily.
%
%   Within any interval of these "monotone" windows, there will be at
%   most one solution of any equality constraint. Since the boundary
%   of the solution set for any inequality constraint is the set 
%   of points where an equality constraint is met, the solutions of
%   both equality and inequality constraints can be found easily
%   once the monotone windows have been found.
%
%   Step Size
%   =========
%
%   The monotone windows (described above) are found using a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of distance (range rate) will be 
%   sampled. Starting at the left endpoint of an interval, samples 
%   will be taken at each step. If a change of sign is found, a 
%   root has been bracketed; at that point, the time at which the 
%   range rate is zero can be found by a refinement  
%   process, for example, using a binary search. 
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the distance function is monotone:
%   the step size should be shorter than the shortest of these
%   intervals (within the confinement window).
%
%   The optimal step size is *not* necessarily related to the lengths
%   of the intervals comprising the result window. For example, if
%   the shortest monotone interval has length 10 days, and if the
%   shortest result window interval has length 5 minutes, a step size
%   of 9.9 days is still adequate to find all of the intervals in the
%   result window. In situations like this, the technique of using
%   monotone windows yields a dramatic efficiency improvement over a
%   state-based search that simply tests at each step whether the
%   specified constraint is satisfied. The latter type of search can
%   miss solution intervals if the step size is shorter than the
%   shortest solution interval.
%
%   Having some knowledge of the relative geometry of the target and
%   observer can be a valuable aid in picking a reasonable step size.
%   In general, the user can compensate for lack of such knowledge by
%   picking a very short step size; the cost is increased computation
%   time.
%
%   Note that the step size is not related to the precision with which
%   the endpoints of the intervals of the result window are computed.
%   That precision level is controlled by the convergence tolerance.
%
%   Convergence Tolerance
%   =====================
%
%   As described above, the root-finding process used by this routine
%   involves first bracketing roots and then using a search process
%   to locate them. "Roots" are both times when local extrema are
%   attained and times when the distance function is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The convergence tolerance used by this
%   routine is set by the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The Confinement Window
%   ======================
%
%   The simplest use of the confinement window is to specify a time
%   interval within which a solution is sought. However, the
%   confinement window can, in some cases, be used to make searches
%   more efficient. Sometimes it's possible to do an efficient search
%   to reduce the size of the time period over which a relatively
%   slow search of interest must be performed.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gfdist_c.
%
%   MICE.REQ
%   GF.REQ
%   SPK.REQ
%   CK.REQ
%   TIME.REQ
%   WINDOWS.REQ
%   
%-Version
%
%   -Mice Version 1.0.0, 15-APR-2009, EDW (JPL)
%
%-Index_Entries
%
%   GF distance search
%
%-&

function [result] = cspice_gfdist( target, abcorr, obsrvr, relate, refval, ... 
                                   adjust, step, nintvls, cnfine )

   switch nargin

      case 9

         target  = zzmice_str(target);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
 
      otherwise

         error ( [ 'Usage: [result] = cspice_gfdist( `target`, `abcorr`, '  ... 
                                     '`obsrvr`, `relate`, refval, adjust, ' ... 
                                     'step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfdist_c', target, abcorr, obsrvr, relate, ... 
                                  refval, adjust, step, nintvls,  ...
                                  [zeros(6,1); cnfine] );
   catch
      rethrow(lasterror)
   end




