%-Abstract
%
%   CSPICE_WNFILD fills small gaps between adjacent intervals of
%   a double precision window. 
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
%      sml      the double precision scalar limiting measure of the small
%               gaps to fill. Adjacent intervals separated by gaps of 
%               measure less than or equal to 'sml' are merged. 
%
%      window   a double precision Mx1 array serving as the Mice version
%               of a SPICE window, containing zero or more intervals
%
%   the call:
%
%      window_f = cspice_wnfild(  sml, window )
%
%   returns:
%
%      window_f   a double precision Nx1 array serving as the Mice version
%                 of a SPICE window, containing zero or more intervals,
%                 representing the original 'window', after adjacent 
%                 intervals separated by 'sml' gaps have been merged.
%                 'window_f' may overwrite 'window'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Let 'window' expand the intervals 
%      %
%      window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%      %   
%      % Apply the following series of calls 
%      %
%      window = cspice_wnfild(  1, window )
%      window = cspice_wnfild(  2, window )
%      window = cspice_wnfild(  3, window )
%      window = cspice_wnfild( 10, window )
%
%   MATLAB outputs:
%
%      window =
%      
%           1
%           3
%           7
%          11
%          23
%          27
%          29
%          29
%
%      Representing the intervals:
%
%         [ 1,  3 ]  [ 7, 11 ]  [ 23, 27 ]  [ 29, 29 ]
%
%      window =
%      
%           1
%           3
%           7
%          11
%          23
%          29
%
%      Representing the intervals:
%
%         [ 1,  3 ]  [ 7, 11 ]  [ 23, 29 ]
%      
%      window =
%      
%           1
%           3
%           7
%          11
%          23
%          29
%
%      Representing the intervals:
%
%         [ 1,  3 ]  [ 7, 11 ]  [ 23, 29 ]
%      
%      window =
%      
%           1
%          11
%          23
%          29
%
%      Representing the intervals:
%      
%         [ 1, 11 ]  [ 23, 29 ]
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnfild_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-SEP-2008, EDW (JPL)
%
%       Edit to I/O return value description.
%
%   -Mice Version 1.0.0, 24-JUL-2007, EDW (JPL)
%
%-Index_Entries
% 
%   fill small gaps in a d.p. window 
% 
%-&

function [window_f] = cspice_wnfild( sml, window)

   switch nargin

      case 2

         sml    = zzmice_dp(sml);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [window_f] = cspice_wnfild( sml, window)' )

   end

%
% Please note, this call does not require addition of space for the window
% control segment as needed by other windows interfaces. The interface
% copies the data in 'window' to a work variable rather than directly 
% pass 'window' into a CSPICE call.
%
   try
      [window_f] = mice('wnfild_c', sml, window );
   catch
      rethrow(lasterror)
   end



