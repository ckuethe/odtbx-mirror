function [y,Py] = sigpt(x,Px,f,time,rsite,hbar)
% SIGPT Sigma Point Transformation
%
%    [y,Py] = sigpt(x, Px, f, time, rsite, hbar) Performs a sigma-point
%    transformation through a nonlinear transformation function f. hbar is
%    an optional input with the default value being sqrt(3) for a Gaussian
%    distribution
%
% INPUTS:
% x - Domain mean (6x1 vector of RA and Dec values for 3-observation IOD)
%     (degrees)
% Px - Domain covariance matrix (6x6 with diagonals equal to variance)
% f - Nonlinear transformation function handle
% time - nx6 matrix of n datevectors for n observation times (datevector)
% rsite - 3xn matrix of site position vectors in the intertial frame (km)
% hbar - Divided difference scale size (sqrt(kurtosis)) (Optional)
% 
% OUTPUTS:
% y - Mean nonlinear transformation output formed from sigma points
% Py - Covariance matrix of mean
%
% CALLS:
% None

% REVISION HISTORY:
%   Author             Date (MM/DD/YYYY)    Comment
%   Russell Carpenter  06/23/2015           Original
%   Ryan Willmot       06/24/2015           Modified to work with IOD App.
%   Russell Carpenter  06/29/2015           Fixed mean output, mergesigpts

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
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

if nargin == 5
    hbar = sqrt(3); % kurtosis
end
n = size(Px,1) ;    
sx = spawnsigpts(x,chol(Px)',hbar,n);
sy = [];
for i = 1:2*n+1
    tempdata = [sx(1:2,i)'; sx(3:4,i)'; sx(5:6,i)'];
    [r2, v2] = feval(f,tempdata,time,rsite);
    next = [r2;v2];
    sy = [sy next];
end
y = mergesigpts(sy, hbar, n);
[D1,D2] = diffsigpts(sy,hbar,n);
Py = [D1 D2]*[D1 D2]';


% Subfunctions extracted from ODTBX estspf.m  
% diffsigpts, spawnsigpts, mergesigpts

function [D1,D2] = diffsigpts(X,hbar,n)
% Compute the divided difference matrices containing information about the
% first and second derivatives of a set of sigma points, S1 and S2,
% respectively.
for i = n:-1:1,
    D1(:,i) = 1/2/hbar*(X(:,i+1) - X(:,i+n+1));
    D2(:,i) = sqrt(hbar^2-1)/2/hbar^2*(X(:,i+1) + X(:,i+n+1) - 2*X(:,1));
end

% Checks for NaN elements and replaces them with a 0
n1 = isnan(D1);
n2 = isnan(D2);
if any(n1),
    D1(n1) = 0;
end
if any(n2),
    D2(n2) = 0;
end

function X = spawnsigpts(x,N,hbar,n)
% Creates the set of sigma points from the mean, x, the Cholesky factor
% of the covariance matrix, N, and the divided difference scale factor, hbar.
% The resulting set is returned as a matrix, X, whose first column is x,
% whose 2nd:nth columns are x+hN, and whose n+st1:2nth columns are x-hN.
X = repmat(x,1,2*n+1) + [zeros(size(x)), hbar*N, -hbar*N];

function x = mergesigpts(X,hbar,n)
% Merge a set of sigma points back into a single vector.
x = (hbar^2-n)/hbar^2*X(:,1) + sum(1/2/hbar^2*X(:,2:end),2);

