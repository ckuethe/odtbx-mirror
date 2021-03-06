%-Abstract
%
%   CSPICE_UNITIM returns the double precision value of an input epoch converted
%   from one uniform time scale to another.
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
%      epoch    {SCALAR or 1xN ARRAY, DOUBLE PRECISION} epoch(s) relative to
%               the 'insys' time scale.
%
%      insys    {SCALAR, STRING} naming the uniform time scale of 'epoch'.
%                Acceptable values:
%
%                  'TAI'     International Atomic Time.
%
%                  'TDB'     Barycentric Dynamical Time.
%
%                  'TDT'     Terrestrial Dynamical Time.
%
%                  'ET'      Ephemeris time (in the SPICE system, this is
%                            equivalent to TDB).
%
%                  'JDTDB'   Julian Date relative to TDB.
%
%                  'JDTDT'   Julian Date relative to TDT.
%
%                  'JED'     Julian Ephemeris Date (in the SPICE system
%                            this is equivalent to JDTDB).
%
%               The routine is not sensitive to the case of insys;  
%               'tai' 'Tai' and 'TAI' are all equivalent from the point of
%               view of this routine.
%
%      outsys   {SCALAR, STRING} naming the uniform time scale to which 'epoch'
%               should be converted. Acceptable values are the same as for
%               'insys'.
%
%               The routine is not sensitive to the case of 'outsys'.
%
%   the call:
%
%      [unitim] = cspice_unitim( epoch, insys, outsys )
%
%   returns:
%
%      unitim   {SCALAR or 1xN ARRAY, DOUBLE PRECISION} time(s) in the system 
%               specified by 'outsys' equivalent to the 'epoch' in the 'insys'
%               time scale.
%
%               'unitimt' returns with the same vectorization measure (N) 
%               as 'epoch'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapseconds kernel, use a meta kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      et = cspice_str2et( 'Dec 19 2003' );
%
%      converted_et = cspice_unitim(et, 'ET','JED')
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      converted_et =
%
%           2.452992500742865e+06
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine unitim_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-APR-2010, EDW (JPL)
%
%-Index_Entries
%
%   Transform between two uniform numeric time systems
%   Transform between two additive numeric time systems
%   Convert one uniform numeric time system to another
%   Convert one additive numeric time system to another
% 
%-&

function [output] = cspice_unitim( epoch, insys, outsys )

   switch nargin
      case 3

         epoch  = zzmice_dp(epoch);
         insys  = zzmice_str(insys);
         outsys = zzmice_str(outsys);

      otherwise

         error( [ 'Usage: [_output_] = ' ...
                          'cspice_unitim( _epoch_, `insys`, `outsys` )'] )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from 
   % this script.
   %
   try
      [output] = mice( 'unitim_c', epoch, insys, outsys);
   catch
      rethrow(lasterror)
   end
