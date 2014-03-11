function [xdot A Q] = pancake_dyn(~,x,wP)
% Get dynamics & estimation values for the pancake tutorial.

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

% Get state information
R = x(1:2,1);
V = x(3:4,1);
mu = x(5,1);
Rs = x(6:7,1);

r3 = norm(R)^3;

wPx = [0 -wP; wP 0];

% Calculate state derivatives
xdot = [V; 
        -mu/norm(R)^3*R; 
        0;
        wPx*Rs];

I = eye(2,2);
    
% Calculate partials
A = zeros(7,7);
A(6:7,6:7) = wPx; 
A(3:4,1:2) = -mu/r3*I + 3*mu*R*R'/r3/norm(R)^2; 
A(3:4,5) = -R/r3; 
A(1:2,3:4) = I;

% Set process noise spectral density
Q = zeros(7,7);
Q(5,5) = 1e-4; % add some noise on mu

end