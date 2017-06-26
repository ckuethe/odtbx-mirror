function [L_Vectors] = SimpleLagrange(time, vectors)
% SIMPLELAGRANGE: Simple Lagrange Interpolation for only three vectors
%
%    [L_Vectors] = SimpleLagrange(time, vectors) Uses the Lagrange 
%    interpolation formula to derive an approximate expression for a 
%    particular position vector, velocity vector, and accelration vector 
%    at the second given time in the time input.
%
% WARNING: This function only handles three vectors. Inputting more
% than three will result in an index error. Use GeneralLagrange for more
% than three input vectors.
%
% INPUTS:
% time - 3x6 matrix of n datevectors for 3 observation times (datevector)
% vectors - 3x3 matrix containing n vectors to be used in Lagrange
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

% Finding the change in observation times (seconds)
t_12 = etime(time(1,:),time(2,:));
t_23 = etime(time(2,:), time(3,:));
t_13 = etime(time(1,:), time(3,:)); 
t_21 = -t_12;
t_32 = -t_23;
t_31 = -t_13;

% Approximate expressions for a vector and its first and second derivatives 
% at the middle observation time

L_gen_dot = t_23*vectors(:,1)/(t_12*t_13) + (t_21 + t_23)*vectors(:,2)/(t_21*t_23) +...
    t_21*vectors(:,3)/(t_31*t_32);    % (s^-1)

L_gen_2dot = 2*vectors(:,1)/(t_12*t_13) + 2*vectors(:,2)/(t_21*t_23) +...
    2*vectors(:,3)/(t_31*t_32);    % (s^-2)

% Output matrix 
L_Vectors = [vectors(:,2) L_gen_dot L_gen_2dot];

end