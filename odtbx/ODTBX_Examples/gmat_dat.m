% GMAT_DAT datfun for the GMAT API examples

function [y, H, R] = gmat_dat(t,x,dat)

lent = numel(t);

% Preallocate
y = zeros(2, lent);
H = zeros(2, 6, lent);
R = repmat(diag([0.05, 5e-5]), 1, 1, lent);

% Get the TrackingFileSet and Spacecraft
tfs = dat{1};
sat = dat{2};

tdas = tfs.GetAdapters();

for tt=1:lent
    % Set the spacecraft state
    sat.SetRealParameter('X', x(1));
    sat.SetRealParameter('Y', x(2));
    sat.SetRealParameter('Z', x(3));
    sat.SetRealParameter('VX', x(4));
    sat.SetRealParameter('VY', x(5));
    sat.SetRealParameter('VZ', x(6));
    
    % Set the spacecraft epoch
    sat.SetEpoch(21545.0 + t(tt)/(3600*24));
    
    for tdaIndex = 1:tdas.size()
        % Get adapter for each measurement
        tda = tdas.get(tdaIndex-1);

        md = tda.CalculateMeasurement(); % Calculate measurement
        y(tdaIndex, tt) = md.getValue().get(0); % Get the measurement

        % Calculate measurement derivatives
        id = sat.GetParameterID('CartesianX');
        tda.CalculateMeasurementDerivatives(sat,id);

        columncount = tda.ApiGetDerivativeValue(0,-1); % Get number of columns in partials

        % Retrieve the measurement's derivatives
        for jj = 1:columncount
            H(tdaIndex, jj, tt) = tda.ApiGetDerivativeValue(0,jj-1);
        end
    end
end
