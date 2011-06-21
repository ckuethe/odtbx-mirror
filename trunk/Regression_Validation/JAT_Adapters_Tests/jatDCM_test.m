function failed = jatDCM_test()
% Regression Test Case
% Function(s) jatDCM
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

failed  = 0;
tol     = 1e-10;
epoch   = datenum('June 1, 2007');

load jatDCM_testdata

validNames = {
    'eci2ecef'
    'ecef2eci'
    'tod'
    'precession'
    'nutation'
    'gha'
    'polar'
    'ecliptic'
    };

n = length(validNames); 

for i=1:n
    d = jatDCM(validNames{i}, epoch);
    if( any( any( abs(d-namevals{i,2}) ) ) > tol )
        failed = 1;
    end
end
