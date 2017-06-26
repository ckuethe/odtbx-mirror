function TX_ID = get_gps_block(options,PRN)
% get_gps_block
% This function finds the correct GPS contellation configuration 
% (PRN, SV, and block type mappings) for the simulation epoch. 
%
% Data Source:    GPS_PRN_Usage_Data.txt
% Data Origin:    Mark E. Strub, The Aerospace Corporation, 
%                 GPS Operations, Schriever AFB CO
% Date Created:   08/21/14
%
% OPTIONS is an OD Toolbox Measurement Options data structure. See
% ODTBXOPTIONS for all available options settings. The options parameters
% that are valid for this function are:
%
%   PARAMETER           VALID VALUES           NOTES
%   epoch               datenum                Time associated with start of simulation
%
%   INPUTS
%   VARIABLE        type    SIZE    DESCRIPTION (Optional/Default)
%   options         struct  1       data structure (see above description)
%   PRN             double  1       (optional) if second argument is used, 
%                                   TX_ID structure for PRN is returned, else, 
%                                   TX_ID{32} is returned indexed by PRN      
%
%
%   OUTPUTS
%   VARIABLE        TYPE    SIZE    DESCRIPTION
%   TX_ID           struct  1       A struct to hold raw GPS measurements
%       its fields:
%   .GPS_ID         double  1       A user-supplied GPS satellite vehicle 
%                                   identifier (note, this is not the GPS PRN)
%                                   (defaulted to -1)
%   .GPS_PRN        double  1       The assigned GPS PRN number for this
%                                   transmitter (defaulted to -1)
%   .block_type     double  1       The GPS Satellite Block type:
%                                               1=II/IIA (default)
%                                               2=IIR
%                                               3=IIR-M
%                                               4=IIF
%                                               5=III
%                                   (This numerical value should be
%                                   consistent with the gpsmeas GPSBlock 
%                                   option.)
%
%   See also odtbxOptions, gpsmeas, makeGpsTXID
%

%  REVISION HISTORY
%   Author      		Date         	Comment
%   Jennifer Valdez     08/21/2014      Original
%   Jennifer Valdez     09/22/2014      Code optimizations


if nargin == 1
    % get all 32 prns
    GPS_SIZE = 32;
    prns = 1:32;
elseif nargin == 2
    % do one prn at a time
    GPS_SIZE = 1;
    prns = PRN;
else
	% error
end

TX_ID = cell(GPS_SIZE,1);

fid = fopen('GPS_PRN_Usage_Data.txt','r');
for ii = 1:8
    % throw out the headers
    hdr = fgetl(fid);
    if ii < 4
        % Print the info?
%         fprintf('%s\n',hdr);
    end
end

% For now we discard the clock information
% SVN #/T PRN Start_Time - Stop_Time   => Clk_Turn-on
[A, COUNT] = fscanf(fid,'%d %*s %d %d %*s %d %*s %d');
Data = reshape(A,5,[])';

% % remove this code later (used for testing only)
% epoch = datenum('Jan 1 2006');
% options = odtbxOptions('measurement');
% options = setOdtbxOptions(options,'epoch',epoch);

for prn = prns
    % Find data by prn
    prn_rows = find(Data(:,2) == prn);
    unixTimeStart = Data(prn_rows,3);
    unixTimeEnd = Data(prn_rows,4);
    
    %% Time of applicability
    
    % WARNING: this code does not handle leap seconds properly
    % times are in seconds since Unix Epoch (1/1/1970 00:00:00)
    % Get the serial date number representing the UNIX epoch
    unixepoch = datenum(1970,1,1,0,0,0);
    
    % Divide the timestamp you want to convert by the number of seconds in a day
    unixDaysStart = unixTimeStart/(60*60*24);
    unixDaysEnd = unixTimeEnd/(60*60*24);
    
    % Add this value to the serial date number representing the epoch
    tstart = unixepoch + unixDaysStart;
    tend = unixepoch + unixDaysEnd;
    
    % Use datestr to format this value as a stringified date to confirm
%     tsStringStart = datestr(tstart);
%     tsStringEnd = datestr(tend);
    
    doa_row = prn_rows(tstart < options.epoch & options.epoch < tend);
   
    if isempty(doa_row)
        % this prn was not transmitting at this epoch, set svn to 0
        svn = 0;
    else
        svn = Data(doa_row,1);
    end
    
%     toaStart = datestr(tstart(doa));
%     toaEnd = datestr(tend(doa));
    
    %% Now which Block is it?
    
    if svn > 0 && svn < 12
        % block I
        block = 6;
    elseif svn > 12 && svn < 41
        % block II/IIA
        block = 1;
    elseif svn > 40 && svn < 47 || svn == 51 || svn == 54 || svn == 56
        % block IIR
        block = 2;
    elseif svn == 47 || svn > 58 && svn < 62
        % block IIR with modernized antenna (behaves like IIR-M)
        block = 3;
    elseif svn > 47 && svn < 51 || svn == 52 || svn == 53 || svn == 55 || svn == 57 || svn == 58
        % block IIR-M
        block = 4;
    elseif svn > 61 && svn < 69
        % block IIF
        block = 5;
    else
        % dummy block
        block = 6;
    end
    
    %% Use TX_ID structure
    % GPS PRN assignments can change but it is the unique vehicle and
    % transmitter that we want to analyze.  We have to keep track of the
    % current PRN against a particular vehicle.  Let's say we're only
    % interested in eight of the GPS satellites we've observed.  Here is the
    % identifying data for our eight satellites.  The ID is our own unique tag
    % for managing the data.
    
    % Allow for overwritting the block definetion if user so desires
    if isfield(options.linkbudget,'GPSBlock')
        block = options.linkbudget.GPSBlock;
    end
    if GPS_SIZE == 32
        TX_ID{prn} = makeGpsTXID(svn, prn, block); % ID/SVN, PRN, block
    elseif GPS_SIZE == 1
        TX_ID = makeGpsTXID(svn, prn, block); % ID/SVN, PRN, block
    end
        
end % for prn
fclose(fid);
end % function get_gps_block
