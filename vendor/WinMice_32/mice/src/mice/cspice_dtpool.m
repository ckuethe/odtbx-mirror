%-Abstract
%
%   CSPICE_DTPOOL returns descriptive data about a kernel pool variable
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
%     name   the scalar string or NxM character array of names of variables
%            whose values are to be returned
%
%   the call:
%
%      [found, n, type] = cspice_dtpool(name)
%
%   returns:
%
%      found   a boolean scalar or boolean 1XN array that returns as true if 
%              the variable 'name' exists in the pool; false if not
% 
%      n       a double precision scalar or double precision 1XN array
%              describing the number of values associated with 'name'.
%              If 'name' does not  exist in the pool, 'n' returns
%              with the value 0. 
% 
%      type    a scalar character or 1XN array of characters indicating the 
%              type of the variable associated with 'name'
% 
%                  'C' if the data is character data 
%                  'N' if the data is numeric
%                  'X' if there is no variable name in the pool
%
%              'found', 'n', and 'type' return with the same vectorization
%               measure (N) as 'name'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapsecond kernel.
%      %
%      cspice_furnsh('standard.tm' )
%
%      %
%      % Check for the variables defined in the leapseconds kernel
%      % and a name probably (hopefully) not in the kernel pool.
%      %
%      lmpoolNames  = strvcat(              ...
%                    'DELTET/DELTA_T_A',    ...
%                    'DELTET/K',            ...
%                    'DELTET/EB',           ...
%                    'DELTET/M',            ...
%                    'ECHO419',             ...
%                    'DELTET/DELTA_AT',     ...
%                    'EVERLASTING_GOBSTOPPER' );
%
%      [found, n, dtype] = cspice_dtpool( lmpoolNames );
%
%      for i = 1:size(lmpoolNames,1)
%      
%         name = lmpoolNames(i,:);
%            
%         if (found(i))  
%            fprintf( 'Variable name : %s\n', name       )
%            fprintf( 'Variable size : %d\n', n(i)       )
%            fprintf( 'Variable type : %s\n\n', dtype(i) )
%         else
%            fprintf( 'Unable to find variable name : %s\n\n', name )
%         end
%   
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Variable name : DELTET/DELTA_T_A      
%      Variable size : 1
%      Variable type : N
%      
%      Variable name : DELTET/K              
%      Variable size : 1
%      variable type : N
%      
%      Variable name : DELTET/EB             
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/M              
%      Variable size : 2
%      Variable type : N
%
%      Unable to find variable name : ECHO419               
%      
%      Variable name : DELTET/DELTA_AT       
%      Variable size : 48
%      Variable type : N
%      
%      Unable to find variable name : EVERLASTING_GOBSTOPPER
%
%-Particulars
%
%   A sister version of this routine exists named mice_dtpool that returns 
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dtpool_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%      
%   -Mice Version 1.0.0, 07-MAR-2007, EDW (JPL)
%
%-Index_Entries
% 
%   return summary information about a kernel pool variable
% 
%-&

function [found, n, type] = cspice_dtpool(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise
     
         error ( 'Usage: [_found_, _n_, _`type`_] = cspice_dtpool(_`name`_)' )

   end

%
% Call the MEX library. The "_s" suffix indicates a structure type 
% return argument.
%
   %
   % Call the MEX library.
   %
   try
      [dtpool] = mice('dtpool_s',name);
      found    = reshape( [dtpool.found], 1, [] );
      n        = reshape( [dtpool.n], 1, [] );
      type     = char( dtpool.type );
   catch
      rethrow(lasterror)
   end



