function [x,rtQd] = clkprop2(tspan,xo,opt)
% A 2-state clock error propagator.

% [x,rtQd] = clkprop2(tspan,xo,opt)
% Units are km, and km/s; set c = 1 for s, and s/s.
%
% See also: gpspseudomeas, clkdyn2
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

% REVISION HISTORY
%   Author      		Date         	Comment
%   Russell Carpenter   08/25/2011      Original

c = JATConstant('c')/1e3;
dt = diff(tspan);
if length(dt) > 1
    error('ODTBX:CLKPROP2:multStep',...
        'CLKPROP2 propagates one step at a time; length(tspan) should be 2')
end
if nargin == 2 || isempty(opt)
    h_0 = 2.4e-22;
    hm2 = 8.0e-28;
    q1 = c^2*h_0/2;
    q2 = c^2*2*pi^2*hm2;
else
    q1 = c^2*opt.q1;
    q2 = c^2*opt.q2;
end
Phi = eye(2) + dt*diag(1,1);
x = Phi*xo;
rtQd(2,2) = sqrt(q2*dt);
rtQd(1,2) = rtQd(2,2)*dt/2;
rtQd(1,1) = (sqrt(q1*dt + q2*dt^3/12));