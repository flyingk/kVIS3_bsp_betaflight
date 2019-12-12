% kVIS3 Data Visualisation
%
% Copyright (C) 2012 - present  Kai Lehmkuehler, Matt Anderson and
% contributors
%
% Contact: kvis3@uav-flightresearch.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function fdses = import_betaflight(file)

if ~nargin
    file = 'C:\Users\matt\OneDrive\Aircraft\Rabbit\logs\2019-05-11_Initial_Testing\btfl_001.bbl';
end

%% Check to make sure blackbox-tools is installed
% otherwise it can be downloaded from https://github.com/cleanflight/blackbox-tools/releases
BSP_Path = getpref('kVIS_prefs','bspDir');
blackbox_decode_directory = [BSP_Path,'\ImportFilter\blackbox_decode\'];

if ~exist([blackbox_decode_directory,'\blackbox_decode.exe'],'file')
    fprintf('Can''t find %s',[blackbox_decode_directory,'\blackbox_decode.exe']);
    % This gives an error but really should be a pop-up
    error('Please download the blackbox_decode folder from https://github.com/cleanflight/blackbox-tools/releases and put it in \ImportFilter\blackbox_decode\');
    %     return;
end


if file==0
    warning('Error loading file.')
    return
end

tic

fprintf('Importing BetaFlight .bbl file\n');

%% Read bbl file
% Get file name parts
[file_path,file_name,file_ext] = fileparts(file);
fprintf('\tCoverting log %s to csv...',file);
tic;

% Convert File
if (ispc)
    evalc(['!"',blackbox_decode_directory,'blackbox_decode.exe" --unit-rotation deg/s --unit-acceleration m/s2 --simulate-imu --unit-height m "',file,'"']);
    
else
    fprintf('Sorry, your platform isn''t support yet... file an issue :)\n');
    return;
end

fprintf('done (%.1f s)\n',toc);

%% Import the Data
% Each file is for a power cycle for which we may have several flights...
% We want to store each of these flights seperately

fdses = {};

for jj = 1:100
    csv_file = [file_path,'/',file_name,'.',num2str(jj,'%02d'),'.csv'];
    
    if (exist(csv_file,'file'))
        fprintf('\t\tImporting Flight %02d ...',jj);
        tic
        
        % Create New fds structure
        fds = kVIS_fdsInitNew();
        fds.BoardSupportPackage = 'BetaFlight';
        [fds, parentNode] = kVIS_fdsAddTreeBranch(fds, 0, 'BetaFlight');
        
        % Import the csv file.  Ignore warnings
        w = warning('off','all');
        data = readtable(csv_file);
        warning(w);
        
        
        % Clean up variable names
        % This should work on data.Properties.VariableDescriptions'
        % and should remove anything after and including the [
        % Should replace [ with _
        
        % for kk = 1:size(yy{1},1)
        %     field_names{kk} = yy{1}{kk};
        %     loc = find(field_names{kk} == '('); % Remove anything in brackets
        %     if (~isempty(loc))
        %         field_names{kk} = field_names{kk}(1:loc-1);
        %     end
        % end
        
        % Print Diagnostics
        fprintf('%5.1f s for readtable()\n',toc);
        
        
        %         data.Properties.VariableNames'
        %         data.Properties.VariableDescriptions'
        
        
        % Add leaves
        groupName = 'data';
        varNames = {};
        varUnits  = {};
        DAT = [];
        
        %% Extract what I want
        t = table2array(data(:,2))/1e6; ...
            t = t - t(1); DAT(:,end+1) = t'; varNames(end+1) = {'Time'}; varUnits(end+1) = {'s'};
        
        dA = 1500 - table2array(data(:,14)); ...
            DAT(:,end+1) = dA'; varNames(end+1) = {'Roll_In'}; varUnits(end+1) = {'us'};
        dE = 1500 - table2array(data(:,15)); ...
            DAT(:,end+1) = dE'; varNames(end+1) = {'Pitch_In'}; varUnits(end+1) = {'us'};
        dR = 1500 - table2array(data(:,16)); ...
            DAT(:,end+1) = dR'; varNames(end+1) = {'Yaw_In'}; varUnits(end+1) = {'us'};
        dT = table2array(data(:,17)); ...
            DAT(:,end+1) = dT'; varNames(end+1) = {'Throttle_In'}; varUnits(end+1) = {'us'};
        
        p_c   = table2array(data(:,18)); ...
            DAT(:,end+1) = p_c';   varNames(end+1) = {'p_c'}; varUnits(end+1) = {'deg/s'};
        Gx = table2array(data(:,25)); ...
            DAT(:,end+1) = Gx'; varNames(end+1) = {'p'}; varUnits(end+1) = {'deg/s'};
        q_c   = table2array(data(:,19)); ...
            DAT(:,end+1) = q_c';   varNames(end+1) = {'q_c'}; varUnits(end+1) = {'deg/s'};
        Gy = table2array(data(:,26)); ...
            DAT(:,end+1) = Gy'; varNames(end+1) = {'q'}; varUnits(end+1) = {'deg/s'};
        r_c   = table2array(data(:,20)); ...
            DAT(:,end+1) = r_c';   varNames(end+1) = {'r_c'}; varUnits(end+1) = {'deg/s'};
        Gz = table2array(data(:,27)); ...
            DAT(:,end+1) = Gz'; varNames(end+1) = {'r'}; varUnits(end+1) = {'deg/s'};
        
        roll  = table2array(data(:,35)); ...
            DAT(:,end+1) = roll'; varNames(end+1) = {'Roll'}; varUnits(end+1) = {'deg'};
        pitch = table2array(data(:,36)); ...
            DAT(:,end+1) = pitch'; varNames(end+1) = {'Pitch'}; varUnits(end+1) = {'deg'};
        yaw   = table2array(data(:,37)); ...
            DAT(:,end+1) = yaw'; varNames(end+1) = {'Yaw'}; varUnits(end+1) = {'deg'};
        
       thr_c = table2array(data(:,21))/10; ...
            DAT(:,end+1) = thr_c'; varNames(end+1) = {'Throttle_Out'}; varUnits(end+1) = {'\%'};
        M1 = table2array(data(:,31))/2500*100; ...
            DAT(:,end+1) = M1'; varNames(end+1) = {'M1'}; varUnits(end+1) = {'\%'};
        M2 = table2array(data(:,32))/2500*100; ...
            DAT(:,end+1) = M2'; varNames(end+1) = {'M2'}; varUnits(end+1) = {'\%'};
        M3 = table2array(data(:,33))/2500*100; ...
            DAT(:,end+1) = M3'; varNames(end+1) = {'M3'}; varUnits(end+1) = {'\%'};
        M4 = table2array(data(:,34))/2500*100; ...
            DAT(:,end+1) = M4'; varNames(end+1) = {'M4'}; varUnits(end+1) = {'\%'};
        
        Voltage = table2array(data(:,22)); ...
            DAT(:,end+1) = Voltage'; varNames(end+1) = {'Volts'}; varUnits(end+1) = {'V'};
        Current = table2array(data(:,23)); ...
            DAT(:,end+1) = Current'; varNames(end+1) = {'Curr'}; varUnits(end+1) = {'A'};
        Power = Voltage.*Current; ...
            DAT(:,end+1) = Power'; varNames(end+1) = {'Power'}; varUnits(end+1) = {'W'};
        CurrTotal = table2array(data(:,38)); CurrTotal(1) = CurrTotal(2); CurrTotal = CurrTotal - CurrTotal(1); ...
            DAT(:,end+1) = CurrTotal'; varNames(end+1) = {'CurrTotal'}; varUnits(end+1) = {'mAh'};
        
        Ax = table2array(data(:,28)); ...
            DAT(:,end+1) = Ax'; varNames(end+1) = {'Ax'}; varUnits(end+1) = {'m/s/s'};
        Ay = table2array(data(:,29)); ...
            DAT(:,end+1) = Ay'; varNames(end+1) = {'Ay'}; varUnits(end+1) = {'m/s/s'};
        Az = table2array(data(:,30)); ...
            DAT(:,end+1) = Az'; varNames(end+1) = {'Az'}; varUnits(end+1) = {'m/s/s'};
        
        rssi = table2array(data(:,24)); ...
            DAT(:,end+1) = rssi'; varNames(end+1) = {'RSSI'}; varUnits(end+1) = {''};
        
        n_channels = numel(varNames);
        varFrames = repmat({'N/A'},n_channels,1);
        
        fds = kVIS_fdsAddTreeLeaf(fds, groupName, varNames, varNames, varUnits, varFrames, DAT, parentNode, false);
        
        % Add new file
        fdses{end+1} = fds;
    end
end

% Remove the csv file
delete([file_path,'\',file_name,'*.csv']);
delete([file_path,'\',file_name,'*.event']);

%% Return the fds struct
return



