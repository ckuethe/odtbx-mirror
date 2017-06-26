function [L_Vectors] = GeneralLagrange(time, vectors)
% GENERALLAGRANGE Generic Lagrange Interpolation
% 
%    [L_Vectors] = GeneralLagrange(time, vectors) Uses the Lagrange 
%    interpolation formula to derive an approximate expression for a 
%    particular position vector, velocity vector, and acceleration vector 
%    at the second given time in the time input. 
%
% This function will return the estimation of the vector position,
% velocity, and acceleration vectors at the second observation time
% (time(2, :)) provided in the input. This function is different from
% SimpleLagrange in that it can handle any number of input vectors instead
% of only three.
%
% INPUTS:
% time - nx6 matrix of n datevectors for n observation times (datevector)
% vectors - 3xn matrix containing n vectors to be used in Lagrange
%           Interpolation
%
% OUTPUTS:
% L_Vectors - 3x3 matrix with the first column being the zeroth derivative,
%             second column being the first derivative and third colum 
%             being the second derivative of the second observation
%             (Units: 1,s^-1, s^-2)
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/04/2015                Original

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

L_temp = [];
syms t;

% Initializes symbolic variable array that will contain symbolic variables
% representing the different observation times
tsym = sym(zeros(1,length(vectors)));
for i = 1:length(vectors)
    tsym(i) = sym(sprintf('t%d', i));
end

for i = 1:length(vectors)
    subvector = tsym([1:i-1 i+1:end]);
    subvector2 = time([1:i-1 i+1:end], [1:6]);
    X = 1;
    for j = 1:length(subvector)
        X = X*(t-subvector(j))/etime(time(i,:),subvector2(j,:));
    end
    L_temp = [L_temp X*vectors(:,i)];
end

% Approximate expressions for a vector and its first and second derivatives 
% at any time
L = sum(L_temp, 2);
L_dot = diff(L, t);
L_2dot = diff(L_dot, t);

% Plugging in t2 for the time, t
L = subs(L, t, tsym(2));
L_dot = subs(L_dot, t, tsym(2));
L_2dot = subs(L_2dot, t, tsym(2));

% Plugging in 0 for the time of interest, t2
L = subs(L, tsym(2), 0);
L_dot = subs(L_dot, tsym(2), 0);
L_2dot = subs(L_2dot, tsym(2), 0);

% Skipping second element because a value of 0 was already substituted in
subvector = tsym([1:1 3:end]);

% Substituing in the difference between the observation times
for i = 1:length(subvector)
    if i == 1
        L = subs(L, subvector(i), etime(time(i,:),time(2,:)));
        L_dot = subs(L_dot, subvector(i), etime(time(i,:),time(2,:)));
        L_2dot = subs(L_2dot, subvector(i), etime(time(i,:),time(2,:)));
    else
        L = subs(L, subvector(i), etime(time(i+1,:),time(2,:)));
        L_dot = subs(L_dot, subvector(i), etime(time(i+1,:),time(2,:)));
        L_2dot = subs(L_2dot, subvector(i), etime(time(i+1,:),time(2,:)));
    end   
end

% Output Matrix
L_Vectors = [L L_dot L_2dot];

% Conversion from symbolic back to double
L_Vectors = double(L_Vectors);

end