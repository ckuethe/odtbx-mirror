%-Abstract
%
%   CSPICE_PI returns the double precision value of the constant pi.
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
%      No input required.
%   
%   the call:
%   
%      double = cspice_pi
%   
%   returns:
%   
%      A double precision scalar, the value of pi to machine precision.
%   
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      >> pi_double = cspice_pi
%      
%      pi_double =
%      
%          3.1416
%      
%      >> sprintf( 'PI epoch: %2.11f', cspice_pi )      
%      
%      ans =
%      
%      PI: 3.14159265359
%
%   The MATLAB system system variable "pi" returns a double precision value
%   for PI that equates the value returned by cspice_pi, to machine
%   roundoff.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pi_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   value of pi
% 
%-&

function [return_val] = cspice_pi

   switch nargin
      case 0
         ;
      otherwise
      
         error ( 'Usage: double = cspice_pi' )

   end

   %
   % Call the MEX library.
   %
   try
      [return_val] =  mice('pi_c');
   catch
      rethrow(lasterror)
   end


