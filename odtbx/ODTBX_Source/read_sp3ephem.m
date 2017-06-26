function [sp3] = read_sp3ephem(filename)

% Really rough first crack at this
% Jennifer Valdez 9/23/2014

%% Open the file
fid = fopen(filename);

if fid == -1
    error('read_sp3ephem: file does not exist');
end

%% Parse the Header
% There are 22 lines of header
% We want to know the number of PRNs on line 3
for line = 1:22
    current_line = fgetl(fid);
    if line == 1
        %TODO: position or velocity flag?
    elseif line == 3 % Store the number of satellites in the SP3 file
        PRN_SIZE = sscanf(current_line(2:end),'%d');
    elseif line == 13 % Which time system are we using?
        TIME_SYS = current_line(10:12);
    else
        % throw away the other lines for now
        % print to screen?
    end
end

%% Time correction factor
if strcmp(TIME_SYS,'GPS')
    datenum_correction = 723186; % for ODTBX
else
    datenum_correction = 0; 
end

%% Initialize sp3 matrix
sp3 = [];
data = NaN(PRN_SIZE,5);

%% Read in all the data
while ~feof(fid)
    % place holder for epoch marker
    current_epoch = fgetl(fid);
    if strcmp(current_epoch(1:3),'EOF')
        break
    else
        current_epoch = sscanf(current_epoch(2:end),'%f',6);
    time = datenum(current_epoch') - datenum_correction;
    % PRN X Y Z (km)
    for ii = 1:PRN_SIZE
        data(ii,2:end) = fscanf(fid,'%*c%*c%d %f %f %f %*f %*d %*d %*d %*d');
        data(ii,1) = time;
    end
    
    sp3 = [sp3;data];
    end
end