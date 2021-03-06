function options = setOdtbxOptions(varargin)

% SETODTBXOPTIONS  Create/alter OPTIONS structure.
%
%   OPTIONS = setOdtbxOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates
%   an options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names.
%
%   OPTIONS = setOdtbxOptions(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS that was created with odtbxOptions.
%
%   OPTIONS = setOdtbxOptions(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties.  If the two options structures
%   are different types, e.g. 'measurement' and 'estimator', then the 
%   resulting structure will be promoted to include all fields.  
%   Unspecified properties of NEWOPTS are left unchanged in OLDOPTS.
%
%   setOdtbxOptions with no input arguments displays all property names and their
%   possible values.
%
%   See ODTBXOPTIONS for the possible options to set.
%
%    keyword: options
%    See also GETODTBXOPTIONS, ODTBXOPTIONS, VALIDATEODTBXOPTIONS, 
%    CREATEJATWORLD, JATFORCES, ODESET
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2011 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

%  REVISION HISTORY
%   Author      						Date         	Comment
%               						(MM/DD/YYYY)
%   Mark W. Reichelt and Lawrence F. Shampine 	5/6/94		Copyright 1984-2003 The MathWorks, Inc.
%   Date: 2004/04/16 22:05:31
%   Kathryn Bradley    	 				04/17/2006   	Ammended original to work for adaptor options
%   Derek Surka                         07/20/2007      Generalized MathWorks & Kate's work
%   Derek Surka                         09/07/2007      Renamed and revised
% (for additional changes, see the svn repository)
%   Allen Brown                         02/26/2009      Updated documentation.
%   Ravi Mathur                         06/26/2014      Fixes input options structs that don't contain all appropriate fields
%                                                       Correctly interprets inputs

if( nargin==0 )
    odtbxOptions(); % call odtbxOptions with no output args to display help
    return;
elseif( nargin == 1 )
    error('ODTBX:setOdtbxOptions:InvalidArgs', 'Invalid number of arguments.');
end

% check the first argument type
if isstruct(varargin{1})
    options = validateOdtbxOptions(varargin{1});
    i = 2;
elseif ischar(varargin{1})
    % 1st arg is not a struct, but a char, assume a name-value pairing
    % start with the default odtbxOptions struct
    options = odtbxOptions;
    i = 1;
else 
    error('ODTBX:setOdtbxOptions:InvalidArgs', 'The first argument is not a valid ODTBX options structure or the name part of a name-value pair.');    
end
arg = varargin{i};

% Assignments between structs
if(isstruct(arg))
    % Make sure second input is a valid options struct
    arg = validateOdtbxOptions(arg);
    
    % Can't have more than one struct from which to get values
    if(nargin > 2)
        error('ODTBX:setOdtbxOptions:InvalidArgs', 'Cannot assign values from more than one ODTBX options structure.');
    end
    
    % Struct content "promotion":
    % If the options.optionsType is 'all', the promotion is safe since 
    % all the required fields exist.
    % If the options.optionsType is anything else, then we'll have to 
    % promote the options struct to an all-inclusive struct via a recursive
    % call.
    if( ~strcmpi(options.optionsType,arg.optionsType) )
        if( ~strcmpi(options.optionsType,'all') )
            % Promotion is required:
            alloptions = odtbxOptions(); % create an 'all' struct
            % copy the old options data into this new 'all' struct
            newoptions = setOdtbxOptions(alloptions,options);
            % rename and keep going:
            options = newoptions; clear alloptions newoptions;
        end
    end
    
    % Copy fields from new options struct
    names = fieldnames(arg);
    for j = 1:length(names)
        val = arg.(names{j});
        if ~isempty(val)
            options.(names{j}) = val;
        end
    end

    return;
end

% Assigning name-value pairs:

% A finite state machine to parse name-value pairs.
names = fieldnames(options);

if rem(nargin-i+1,2) ~= 0
    error('ODTBX:setOdtbxOptions:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end

expectval = false;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i}; 

    if ~expectval
        if ~ischar(arg)
            error('ODTBX:setOdtbxOptions:NoPropName',...
                'Expected argument %d to be a string property name.', i);
        end

        [j,fullname] = getIndex(arg,names);
        expectval = true;                      % we expect a value next

    else
        options.(fullname) = arg;
        expectval = false;

    end
    i = i + 1;
end
