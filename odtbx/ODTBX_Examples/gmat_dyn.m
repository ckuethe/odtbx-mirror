% GMAT_DYN dynfun for the GMAT API examples

function [xDot, A, Q] = gmat_dyn(t, x, fm)

nt = numel(t);

xDot = zeros(6, nt);
A = zeros(6, 6, nt);

for tIndex=1:nt
    % Get the state into a format GMAT can work with
    state = gmat.gmat.convertJavaDoubleArray(x(:,tIndex));

    % Get derivatives
    fm.GetDerivatives(state, t(tIndex), 1); % Calculate derivatives
    deriv = fm.GetDerivativeArray(); % Get calculated derivatives
    derivArray = gmat.gmat.convertDoubleArray(deriv, 42); % Convert array from GMAT to MATLAB format

    xDot(:,tIndex) = derivArray(1:6);

    % Reshape the vector of state derivatives into the A matrix
    % Need to skip the first 6 elements as they do not contain Jacobian data
    A(:,:,tIndex) = reshape(derivArray(7:42), 6, 6)';
end

% Populate process noise
Q = repmat(diag([0 0 0 1e-9 1e-9 1e-9].^2), 1, 1, nt);
