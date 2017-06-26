function [AntLB, HVIS, AntVis] = calc_linkbudgets(out, options, RX_link, TX_link)

% CALC_LINKBUDGETS  Calculates link budget and visibility between satellite targets
%
% [AntLB, HVIS] = getlinkbudget(out, options, RX_link, TX_link) calculates
% link budget based on information in OPTIONS given information about
% transmitter and receivers and their locations.
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings.  The options parameters
% that are valid for this function are:
%


%   INPUTS
%   VARIABLE        SIZE    DESCRIPTION (Optional/Default)
%      out           (1x1)  Output structure containing the following
%                           fields:
%                              out.epoch      (1 x 1) 
%                                  epoch in Matlab datenum
%                              out.TX_az      (N x GPS_SIZE) 
%                                  transmitter azimuth [deg]
%                              out.TX_el      (N x GPS_SIZE) 
%                                  transmitter elevation [deg]
%                              out.RX_az      (N x GPS_SIZExnum_ant) 
%                                  receiver azimuth [deg]
%                              out.RX_el      (N x GPS_SIZExnum_ant) 
%                                  receiver elevation [deg]
%                              out.range      (GPS_SIZE x N) 
%                                  range [km]
%                              out.rrate      (GPS_SIZE x N) 
%                                  range rate [km/sec]
%                              out.GPS_yaw    (N x GPS_SIZE) 
%                                  GPS SV yaw angle
%                              out.rgps_mag   (N x GPS_SIZE) 
%                                  magnitude of the GPS SV position [km]
%                              out.health     (N x GPS_SIZE) 
%                                  health flag of each GPS SV
%                              out.prn        (1 x H) 
%                                  PRN ID of each healthy GPS SV in this out struct
%
%                              If the params.doH == 1, then the following 
%                              fields are also set:
%                                  out.eciRotation   (9 x 9 x N) 
%                                     rotation matrix from ECI to ECEF
%                                  out.los_3d        (3 x N x GPS_SIZE) 
%                                     LOS 3D vector
%                                  out.gps_vel_tot   (3 x N x GPS_SIZE) 
%                                     total velocity of GPS SVs in ECEF
%                                     frame [km/s]
%                                  out.sat_vel_tot   (3 x N) 
%                                     total velocity of the satellite in
%                                     ECEF frame [km/s]
%                                  out.Rotation2ECI  (1 x 1) 
%                                     name of function that computes 
%                                     rotation from input frame to ECI
%      options      (1x1)   data structure (see above description)
%      RX_link      (1x1)
%      TX_link      (1x1)
%
%   OUTPUTS
%      AntLB        {num_antx1}  Cell array containing link budget for each
%                           antenna.  Note that "masked" refers to an applied
%                           bias due to visibility or blockage.  The fields
%                           are:
%           Halpha_r    (MxN)   The receiver elevation angle (rad)
%           Halpha_t    (MxN)   The transmitter elevation angle (rad)
%           Hvis_beta   (MxN)   Logical array where both transmit and
%                               receive antenna elevation angles are
%                               within antenna mask angle limits
%           Hvis_CN0    (MxN)   Logical array where CN0 is above the
%                               acquisition/tracking threshold
%           HCN0        (MxN)   Masked Signal carrier to noise ratio
%           HAd         (MxN)   Masked attenuation from R^2 losses (dB)
%           HAr         (MxN)   Receive antenna gain (dB)
%           HAP         (MxN)   Masked budget gain before receiver antenna (dB)
%           HRP         (MxN)   Masked budget gain before receiver amplifiers
%                               and conversion (dB)
%           HAt         (MxN)   Masked transmit antenna gain (dB)
%       HVIS
%
%   keyword: measurement
%   See also gpsmeas, odtbxOptions, getgpsmeas, gpslinkbudget
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
%   Author      		Date         	Comment
%   Phillip Anderson    10/16/2013      Abstracted out of gpsmeas.m

%% Process incoming data

% Generic constant
d2r             = pi/180;

% Get Some Constants from JAT
EARTH_RADIUS = JATConstant('rEarth','WGS84') / 1000;  % km Equatorial radius of Earth

% Get link budget information
link_budget = getOdtbxOptions(options, 'linkbudget', []);

%% Error check incoming data
% Generate an error if there isn't link budget information
if isempty(link_budget)
    warning('linkbudget variable was undefined in odtbxOptions structure')
end

% Check that we have the necessary information for a link budget
% calculation. If not, assume a value and warn user.
options = linkbudget_default(options, 'frequencyTransmit', JATConstant('L1Frequency') );  % Hz

link_budget = linkbudget_default(link_budget, 'ReceiverNoise', -3 );  % dB, Noise figure of receiver/LNA
link_budget = linkbudget_default(link_budget, 'RecConversionLoss', -1.5 );  % dB
link_budget = linkbudget_default(link_budget, 'NoiseTemp', 300); % K
link_budget = linkbudget_default(link_budget, 'SystemLoss', 0 ); % dB, System losses, in front of LNA
link_budget = linkbudget_default(link_budget, 'AtmAttenuation', 0.0); % dB
link_budget = linkbudget_default(link_budget, 'TXAntennaMask', pi); % rad
link_budget = linkbudget_default(link_budget, 'RXAntennaMask', pi); % rad
link_budget = linkbudget_default(link_budget, 'TX_AntennaPointing', -1); % nadir-pointing by default
link_budget = linkbudget_default(link_budget, 'AntennaPattern', {'sensysmeas_ant.txt','sensysmeas_ant.txt'});
    %  Specify antenna pattern for each antenna, existing antennas are:
    %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
    %     omni.txt                  - zero dB gain,  180 degree half beamwidth
    %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
    %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
    %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
    
% Reassign the options structure with any changed/default link budget values
options = setOdtbxOptions(options, 'linkbudget', link_budget);

num_ant = length(link_budget.AntennaPattern); %hasn't been tested for >4 antennas

%% Assign variables to link structures
% Set receiver and transmitter structure data from link budget information
RX_link.Nf = link_budget.ReceiverNoise;
RX_link.L = link_budget.RecConversionLoss;
RX_link.freq = getOdtbxOptions(options, 'frequencyTransmit', JATConstant('L1Frequency') );  % Hz;
RX_link.Ts = link_budget.NoiseTemp;
RX_link.As = link_budget.SystemLoss;
RX_link.Ae = link_budget.AtmAttenuation;

% Transmitter and receiver antenna masks
TX_link.beta = link_budget.TXAntennaMask;  % User input from options
RX_link.beta = link_budget.RXAntennaMask;   % User input from options

% Measurement physical parameter results:
TX_link.alpha = out.TX_el*d2r; % The transmitter elevation angle (rad) [nn x GPS_SIZE]
RX_link.alpha = out.RX_el*d2r; % The receiver elevation angle (rad)

% Note, the receiver angles from "out" are below, in the ANT loop
Hrange = out.range;      % [GPS_SIZE x nn]
% Hrrate = out.rrate;      % [GPS_SIZE x nn]
rgps_mag = out.rgps_mag; % [nn x GPS_SIZE]
health = out.health;     % the health indicator, [nn x GPS_SIZE]

%% Calculate link constraint and mask information
% Set alpha_t for non-existent/unhealthy SVs to pi rad
TX_link.alpha(~health) = pi;
% Set alpha_r for non-existent/unhealthy SVs to 180 deg
RX_link.alpha(~health) = pi;

% Compute angle subtended by Earth and Earth mask angles for each SV
r_mask          = EARTH_RADIUS + link_budget.AtmosphereMask;	% Atmosphere mask radius (km)
denom=rgps_mag;
denom(denom==0)=NaN;
gamma = asin(EARTH_RADIUS./denom);   % Angle subtended by Earth at SV (nn,GPS_SIZE)
gamma_mask = asin(r_mask./denom);    % Angle subtended by Earth plus altitude mask (nn,GPS_SIZE)
%  Set prns visible if not blocked by Earth of horizon
if link_budget.TX_AntennaPointing == -1 % nadir pointing
    vis_earth = (TX_link.alpha > gamma) | (Hrange' <= rgps_mag.*cos(gamma));   % (nn,GPS_SIZE)
    %  Set prns visible if not subject to atmosphere mask
    vis_atm = (TX_link.alpha > gamma_mask) | (Hrange' <= rgps_mag.*cos(gamma_mask)); % (nn,GPS_SIZE)
elseif link_budget.TX_AntennaPointing == 1 % zenith pointing
    vis_earth = ((-out.TX_el+90) > getOdtbxOptions(options, 'gsElevationConstraint', 10));
    vis_atm = vis_earth;
end

Hvis_earth = vis_earth';
Hvis_atm = vis_atm';

%% Compute antenna specific data

% Set receiver antenna loop number
loop = max([1,num_ant]);
TARGET_SIZE = size(out.range,1);

% ----------------------------------------
%  Antenna calculation loop
% ----------------------------------------
AntLB = cell(loop,1); % cell array of structs to hold link budget data
                      % for each antenna

for ANT=1:loop
    
    AntLB{ANT} = struct('Halpha_r',[],'Halpha_t',[],'Hvis_beta',[],'Hvis_CN0',[],...
        'HCN0',[],'HAd',[],'HAr',[],'HAP',[],'HRP',[],'HAt',[]);
    
    AntLB_raw = struct('CN0',[],'Ad',[],'Ar',[],'AP',[],'RP',[],'At',[]);
   
    % Determine if the pattern is elevation only (1-D) or azimuth and
    % elevation (2-D) and compute the receiver gain
    for j = 1:TARGET_SIZE
        % Encapsulate RX and TX data
        % Originally set to be 1D receive, 1D transmit patterns
        RX_antenna = struct('pattern', RX_link.pattern{ANT}, ...
            'el', RX_link.alpha(:,j,ANT));
        TX_antenna = struct('pattern', TX_link.pattern, ...
            'el', TX_link.alpha(:,j));
        
        % Change dimensions on transmit patterns from 1D to 2D, if required
        if size(RX_link.pattern{ANT},2) > 2
            % 2D receive antenna
            RX_antenna.az = out.RX_az(:,j)*d2r;
        end
        if size(TX_link.pattern,2) > 2
            % 2D transmit antenna
            TX_antenna.az = out.TX_az(:,j)*d2r;
        end
         
         % Compute gain/attenuation of receiving antenna pattern
        [AntLB_raw.CN0(:,j), AntLB_raw.Ar(:,j), AntLB_raw.At(:,j), ...
            AntLB_raw.Ad(:,j), AntLB_raw.AP(:,j), AntLB_raw.RP(:,j)] = ...
            linkbudget(Hrange(j,:)', RX_link, TX_link, RX_antenna, TX_antenna);
    end

    % Apply the receiver gain penalty for the user-defined mask angle
    AntLB_raw = gainpenalty_mask(RX_antenna, AntLB_raw, RX_link.alpha(:,:,ANT), RX_link.beta);
    
    % Apply the transmit gain penalty for the user-defined mask angle
    AntLB_raw = gainpenalty_mask(TX_antenna, AntLB_raw, TX_link.alpha, TX_link.beta);
    
    %------------------------------------------------------------------------------
    % EVALUATION OF GEOMETRIC CONSTRAINTS

    %  Set prns visible if los within antenna mask angles
    %  So far, alpha_t was computed assuming GPS antenna is nadir pointing
    vis_beta_t = (TX_link.alpha <= TX_link.beta);
    vis_beta = vis_beta_t & (RX_link.alpha(:,:,ANT) <= RX_link.beta);    % (nn,GPS_SIZE)

    %  Set prns visible if CN0 is above acquisition/tracking threshold
    vis_CN0 = AntLB_raw.CN0 >= link_budget.RecTrackThresh;                               % [nn,GPS_SIZE]

    %  OUTPUT PARAMETERS
    AntLB{ANT}.Halpha_r = RX_link.alpha(:,:,ANT)';             % [GPS_SIZE,nn]
    AntLB{ANT}.Halpha_t = TX_link.alpha(:,:)'; % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_beta = vis_beta';             % [GPS_SIZE,nn]
    AntLB{ANT}.Hvis_CN0 = vis_CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HCN0 = AntLB_raw.CN0';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAd = AntLB_raw.Ad';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAr = AntLB_raw.Ar';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAP = AntLB_raw.AP';             % [GPS_SIZE,nn]
    AntLB{ANT}.HRP = AntLB_raw.RP';             % [GPS_SIZE,nn]
    % Note: this variable can be used to pass out other values as well
    AntLB{ANT}.HAt = AntLB_raw.At';             % [GPS_SIZE,nn]

    % Mask undefined values for SV dependent parameters using (Health,
    %  Earth blockage, and xmit antennna masks)
    VIS_sv = vis_beta_t' & Hvis_earth & health';             % [GPS_SIZE,nn]
    AntLB{ANT}.HAd(~VIS_sv) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HRP(~VIS_sv) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HAP(~VIS_sv) = -300;             % [GPS_SIZE,nn]

    % Mask undefined values for Antenna dependent parameters using
    % (Health, Earth blockage, and both antennna masks)
    VIS_ant = vis_beta' & Hvis_earth & health';             % [GPS_SIZE,nn]
    AntLB{ANT}.HCN0(~VIS_ant) = -300;             % [GPS_SIZE,nn]
    AntLB{ANT}.HAt(~VIS_ant) = -300;             % [GPS_SIZE,nn]
    
end % ANT loop

%% Combine results across multiple antennas
[HVIS, AntVis] = visibility_constraints(AntLB, options, health, Hvis_earth, Hvis_atm);
end

