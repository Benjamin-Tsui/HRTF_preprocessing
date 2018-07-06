function [ bad_folder, bad_sofa, checked_sofa ] = check_hrtf_folder( folder, angle_range, dist_range, plot )
%CHECK_HRTF_FOLDER Summary of this function goes here
% 
% First check if there is abnormal sofa file in input folder or 
% Then check if there is abnormal hrtf measurement in sofa file 
% (abnormal sofa file will be removed when checking for abnomral hrtf measurement)

% INPUT: 
% e.g. check_hrtf_folder('ITA_HRTF_Database/SOFA', 0.01, 0.1, 0 );
%
% 1. folder = 'ITA_HRTF_Database/SOFA' % folder contains sofa file
%    - example:
%      'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', 'CIPIC_hrtf_database 2/sofa' ...
%      'IRCAM_Listen_hrtf_database/hrtf', 'RIEC_HRTF_Database/hrtf', 'SADIE_HRTF_Database/hrtf'
%
% 2. angle_range = 0.01 % error tolerance in measurement angle (in degrees)
%    - default(if empty) = 0
%
% 3. dist_range = 0.1 % error tolerance in measurement distance (in meters)
%    - default(if empty) = 0
% 
% 4. plot = 0 % plot trigger (0 = no plot, else = plot)
%    - plot abnormal result
%    - default(if empty) = 0 (no plot)
%  
%
% OUTPUT:
% 1. bad_folder % table about bad sofa file in the folder
%    - column 1: Folder (input folder)
%    - column 2: Measurement (inconsistance number of measurement, probably
%                duplicate measurement angles, or missing angle in some sofa files)
%    - column 3: Channels (inconsistance channel, probably missing channel in some files) 
%    - column 4: HRIR_Length (inconsistance hrir length in some sofa files)
%    - column 5: Sampling_Rate (inconsistance sampling rate in some sofa files)
%
% 2. bad_sofa % table about bad measurement in sofa file (if there is any)
%    - column 1: File_Name (SOFA file that contains problematic measurements)
%    - column 2: Repeat_Angles (repeated measurement angles)
%    - column 3: Inconsist_Dist (inconsistance measurement distance) 
%    - column 4: Asymm_Angles (Asymmetrical measurement angles on left and right)
%    - column 5: Asymm_Medium (Asymmetricalmeasurement angles on medium plane,
%                compare front and back, removed top (90) and bottom (-90))   
%
% 3. checked_sofa % table about all checked sofa file
%    - column (1-5) configuration is same as output 2.bad_sofa
%

%% intialise input and output data

if nargin == 1
	angle_range = 0;
    dist_range = 0;
	plot = 0;
elseif nargin == 2
    dist_range = 0;
	plot = 0;
elseif nargin == 3
    plot = 0;
end
% catch empty input


file_names = dir (folder);
sofa_file = {file_names([file_names(:).isdir]==0).name}';
% get all file names inside the folder

bad_folder.folder = folder;
bad_folder.measurement = [];
bad_folder.channels = [];
bad_folder.hrirLength = [];
bad_folder.samplingRate = [];
% initialise output error log

sofa_file_new = sofa_file;
% initialise output sofa_file

%% extract variables

total_measurements = zeros(length(sofa_file), 1);
channels = zeros(length(sofa_file), 1);
hrir_length = zeros(length(sofa_file), 1);
Fs = zeros(length(sofa_file), 1);
for n = 1:length(sofa_file)
    hrtf = SOFAload(sofa_file{n});
    total_measurements(n) = size(hrtf.Data.IR, 1);
    channels(n) = size(hrtf.Data.IR, 2);
    hrir_length(n) = size(hrtf.Data.IR, 3);
    Fs(n) = hrtf.Data.SamplingRate;
end

%% check hrtf measurement count

total_measurements_unique = unique(total_measurements);
% unique hrtf measurement count

if length(total_measurements_unique) > 1
    % if there is more than 1 total hrtf measurement count
    
    if plot ~= 0
        figure
        bar(total_measurements)
        title ('hrtf measurements number in each file')
        % plot data 
    end
    
    total_measurements_hist = zeros(size(total_measurements_unique));
    for n = 1:length(total_measurements_unique)
        total_measurements_hist(n) = sum(total_measurements(:) == ...
            total_measurements_unique(n));
    end
    [~, idx] = sort(total_measurements_hist);
    total_measurements_unique = total_measurements_unique(idx);
    total_measurements_hist = total_measurements_hist(idx);
    % measurements histogram (sort by the histogram, so outliner goes first)
    
    for n = 1:length(total_measurements_hist)-1
        loc = find(total_measurements(:) == total_measurements_unique(n));
    end
    % find outliner location
    
    message = sprintf('non-uniform number of hrtf measurement, file name: \n');
    loc_message = sprintf([char(sofa_file(loc)) '\n']);
    sldiagviewer.reportWarning([message loc_message]);
    % print warning
    
    bad_folder.measurement = sofa_file(loc);
    % save outliner file name to log
    
    sofa_file_new(loc) = [];
    % remove outliners
end

%% check channels

channels_unique = unique(channels);
% find unique channel numbers

if length(channels_unique) > 1 
    % if there is more than 1 channel configuration

    if plot ~= 0
        figure
        bar(channels)
        title ('channel number in each file')
        % plot data 
    end
    
    channels_hist = zeros(size(channels_unique));
    for n = 1:length(channels_unique)
        channels_hist(n) = sum(channels(:) == channels_unique(n));
    end
    [~, idx] = sort(channels_hist);
    channels_unique = channels_unique(idx);
    channels_hist = channels_hist(idx);
    % channels histogram (sort by the histogram, so outliner goes first)
    
    for n = 1:length(channels_hist)-1
        loc = find(channels(:) == channels_unique(n));
    end
    % find outliner location
    
    message = sprintf('non-uniform channel number, file name: \n');
    loc_message = sprintf([char(sofa_file(loc)) '\n']);
    sldiagviewer.reportWarning([message loc_message]);
    % print warning
    
    bad_folder.channels = sofa_file(loc);
    % save outliner file name to log
    
    sofa_file_new(loc) = [];
    % remove outliners
end

%% check hrir length (filter size)

hrir_length_unique = unique(hrir_length);
% find unique hrir length

if length(hrir_length_unique) > 1 
    % if there is more than 1 hrir length
    
    if plot ~= 0
        figure
        bar(hrir_length)
        title ('hrir length in each file')
        % plot data
    end
    
    hrir_length_hist = zeros(size(hrir_length_unique));
    for n = 1:length(hrir_length_unique)
        hrir_length_hist(n) = sum(hrir_length(:) == hrir_length_unique(n));
    end
    [~, idx] = sort(hrir_length_hist);
    hrir_length_unique = hrir_length_unique(idx);
    hrir_length_hist = hrir_length_hist(idx);
    % channels histogram (sort by the histogram, so outliner goes first)
    
    for n = 1:length(hrir_length_hist)-1
        loc = find(hrir_length(:) == hrir_length_unique(n));
    end
    % find outliner location
    
    message = sprintf('non-uniform hrir length, file name: \n');
    loc_message = sprintf([char(sofa_file(loc)) '\n']);
    sldiagviewer.reportWarning([message loc_message]);
    % print warning
    
    bad_folder.hrirLength = sofa_file(loc);
    % save outliner file name to log
    
    sofa_file_new(loc) = [];
    % remove outliners
end

%% check sampling rate

Fs_unique = unique(Fs);
% find unique sampling rate

if length(Fs_unique) > 1
    % if there is more than one sampling rate in folder
    
    warning('non-uniform sampling rate')
    
    if plot ~= 0
        figure
        bar(Fs)
        title ('sampling rate of each file')
        % plot data
    end
    
    Fs_hist = zeros(size(Fs_unique));
    for n = 1:length(Fs_unique)
        Fs_hist(n) = sum(Fs(:) == Fs_unique(n));
    end
    [~, idx] = sort(Fs_hist);
    Fs_unique = Fs_unique(idx);
    Fs_hist = Fs_hist(idx);
    % channels histogram (sort by the histogram, so outliner goes first)
    
    for n = 1:length(Fs_hist)-1
        loc = find(Fs(:) == Fs_unique(n));
    end
    % find outliner location
    
    message = sprintf('non-uniform sampling rate, file name: \n');
    loc_message = sprintf([char(sofa_file(loc)) '\n']);
    sldiagviewer.reportWarning([message loc_message]);
    % print warning
    
    bad_folder.samplingRate = sofa_file(loc);
    % save outliner file name to log
    
    sofa_file_new(loc) = [];
    % remove outliners
end

%% check sofa file inside the folder (without outliners) and orgainise output

checked_sofa = cell2table(cell(length(sofa_file_new),5));
checked_sofa.Properties.VariableNames = {'File_Name','Repeat_Angles',...
    'Inconsist_Dist','Asymm_Angles', 'Asymm_Median'};
% initialise output
reverseStr = '';
for n = 1: length(sofa_file_new)
    msg = sprintf(['\nchecking: ' sofa_file_new{n}  ' (%d/%d files)\n'], ...
        n, length(sofa_file_new));
	fprintf([reverseStr, msg]);
	reverseStr = repmat(sprintf('\b'), 1, length(msg));
    % print progress
    
    [ bad_data_temp ] = check_SOFA(sofa_file_new{n}, angle_range, dist_range, plot);
    % check data in the reamining sofa files
    
    checked_sofa(n,:)= struct2table(bad_data_temp,'AsArray',true);  
    % put result in table
end
% check reamining sofa files and put it into a table

bad_sofa = checked_sofa;
idx=all(cellfun(@isempty,bad_sofa{:,2:end}),2);
bad_sofa(idx,:)=[];
% remove empty row (those have nothing wrong) 


bad_folder = struct2table(bad_folder,'AsArray',true); 
% change bad_sofa output to table

end

