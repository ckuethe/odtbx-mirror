function [y,H,R,AntLB,dtsv,AntVis] = gpspseudomeas(t,x,opts,qatt,prn)
% Creates GPS pseudo-measurements (one-way measurements) biased by onboard clock errors.

% [y,H,R] = GPSPSEUDOMEAS(t,x,options) creates GPS pseudo-measurements,
% i.e. one-way measurements biased by onboard clock errors.  The clock
% errors must be available as part of the state vector, x.  The function
% requires a field 'clockStateIndices' to be present in the options
% structure that indicates the location in the state vector for the clock
% bias and clock drift states.  For example, if the state vector has clock
% bias and drift in elements 7 and 8, respectively, then
%   <opts>.clockStateIndices = [7;8];
% where <opts> is the name of the options structure.  Note that ODTBX
% estimators require a dynamics or propagation model for the clock error
% states, such as CLKDYN2 or CLKDYN3.  Note also that this function assumes
% that the clock model produces errors with units of km and km/s.
%    Refer to the help for GPSMEAS for information on the various options
% and capabilities for modeling GPS measurements.  Note that not all
% options available in GPSMEAS may be accessible via ODBTXOPTIONS; for such
% options, the user must manually add the desired fields and their contents
% to the options structure.  This is also the case for the field this
% function requires, clockStateIndices.
%
%   keyword: measurement
%   See also: gpsmeas, clkdyn2, clkdyn3
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
%   Author      		Date         	Comment
%   Russell Carpenter   08/25/2011      Original

if nargin < 4
    qatt = [];
end

if nargin < 5
    prn = [];
end

% [y,H,R,AntLB] = gpsmeas(t,x,opts,qatt);
[y,H,R,AntLB,dtsv,AntVis] = gpsmeas(t, x, opts,qatt,prn);
clkind = opts.clockStateIndices; % where in the state the clock states are
if opts.useDoppler
    if isfield(opts,'GPSBand'),
        band = opts.GPSBand;
        switch band
            case 'L1'
                f = 1575.42e6;    % Hz
            case 'L2'
                f = 1227.6e6;     % Hz
            case 'L5'
                f = 1176.45e6;    % Hz
        end
    else
        f = 1575.42e6;
    end
    c = JATConstant('c')/1e3;
    h = f/c;
    x(clkind(2),:) = h*x(clkind(2),:);
else
    h = 1;
end
if opts.useRange && (opts.useRangeRate || opts.useDoppler)
    rngbias = repmat(x(clkind(1),:),[size(y,1)/2 1]);
    rdtbias = repmat(x(clkind(2),:),[size(y,1)/2 1]);
    y(1:2:end,:) = y(1:2:end,:) + rngbias;
    y(2:2:end,:) = y(2:2:end,:) + rdtbias;
    H(1:2:end,clkind(1),:) = ones(size(rngbias));
    H(2:2:end,clkind(2),:) = h*ones(size(rdtbias));
elseif opts.useRange
    bias = repmat(x(clkind(1),:),[size(y,1) 1]);
    y = y + bias;
    H(:,clkind(2),:) = zeros(size(y));
    H(:,clkind(1),:) = ones(size(y));
elseif (opts.useRangeRate || opts.useDoppler)
    bias = repmat(x(clkind(2),:),[size(y,1) 1]);
    y = y + bias;
    H(:,clkind(2),:) = h*ones(size(y));
end
    