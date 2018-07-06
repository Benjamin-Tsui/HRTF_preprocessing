function [ bad_folder, bad_sofa, checked_sofa ] = check_hrtf( input, angle_range, dist_range, plot )
%CHECK_HRTF Summary of this function goes here
% 
% First check if there is abnormal sofa file in input folder or 
% Then check if there is abnormal hrtf measurement in sofa file 
% (abnormal sofa file will be removed when checking for abnomral hrtf measurement)

% INPUT: 
% e.g. [ bad_folder, bad_sofa, checked_sofa ] = check_hrtf({'ITA_HRTF_Database/SOFA', ...
%       'ARI_hrtf_database/hrtf', 'MRT02.sofa', 'RIEC_hrir_subject_008.sofa', ...
%        'hrtf b_nh15.sofa', 'MRT04.sofa'}, 0.01, 0.1, 0 );
%
% 1. input = {'ITA_HRTF_Database/SOFA'} % folder contains sofa file
%    or
%    input = {'hrtf b_nh15.sofa'} % sofa file
%    - cell array with folder directory or sofa file (could mix)
%    - input should be cell array (although the function could convert
%      character into to cell array) 
%    - example:
%      'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', 'CIPIC_hrtf_database 2/sofa' ...
%      'IRCAM_Listen_hrtf_database/hrtf', 'RIEC_HRTF_Database/hrtf', 'SADIE_HRTF_Database/hrtf' ...
%      'MRT02.sofa', 'RIEC_hrir_subject_008.sofa', 'hrtf b_nh15.sofa
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
%    - !! suggest only use on SINGLE sofa file or folder 
%      (otherwise lots of windows may pop-out!!!)
%  
%
% OUTPUT:
% 1. bad_folder % table about bad sofa file in folder (if there is any)
%    - column 1: Folder (folder that contains problematic sofa files)
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
%

%% organise input 

if ischar(input)
    input = {input};
end
% catch if input is char instead of cell array


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
% catch missing inputs

input = reshape(input,[], 1);
% reshape input into vertical array

idx = zeros(1, length(input));
for n = 1 : length(input)
    idx(n) = isdir(input{n});
end
folder = input(idx ==1);
% group folder from input
rest = input(idx ==0);
% non-folder in input

sofa_file = cell(size(rest));
for n = 1 : length(rest)
    [~, ~, ext] = fileparts(rest{n});
    if contains(ext, '.sofa')
        sofa_file{n} = rest{n};
    else
        error(['input ' rest{n} ' is not a folder nor sofa file'])
        % kill if there is invalid input (nither folder nor sofa file)
    end
end
sofa_file = sofa_file(~cellfun('isempty',sofa_file));
% group sofa files in input (throw warning if input in not a folder not sofa file)

warning off % temperaly switch off warning
%% check folders

if ~isempty(folder)
% if there is folder from input
    for n = 1 : length(folder)
        if any(size(dir([folder{n} '/*.sofa' ]),1)) == 0
            error(['no .sofa file in folder ' folder{n}])
        end
    end

    bad_folder = cell2table(cell(length(folder),5));
    bad_folder.Properties.VariableNames = {'Folder','Measurement',...
        'Channels','HRIR_Length', 'Sampling_Rate'};

    fprintf(['\nchecking folder: ' folder{1}  ' (%d/%d folders)'], 1, length(folder));
    % print progress
    [bad_folder(1, :), bad_sofa, checked_sofa] = check_hrtf_folder(folder{1},...
        angle_range, dist_range, plot);
    % check 1st hrtf folders for outliners
    if length(folder) > 1 % if there is more than 1 input folder
        for n = 2 : length(folder)
            fprintf(['\nchecking folder: ' folder{n}  ' (%d/%d folders)'],n,length(folder));
            % print progress
            
            [bad_sofa_temp, bad_data_temp, checked_data_temp] = check_hrtf_folder(folder{n}, ...
                angle_range, dist_range, plot);
            % check hrtf folders for outliners
            
            bad_folder(n, :) = bad_sofa_temp;
            bad_sofa = [bad_sofa ; bad_data_temp];
            checked_sofa = [checked_sofa; checked_data_temp];
            % group output
        end
    end

idx=all(cellfun(@isempty,bad_folder{:,2:end}),2);
bad_folder(idx,:)=[];
% remove empty row (those have nothing wrong)

end

%% check sofa files

if ~isempty(sofa_file)
% if there is sofa file from input     
    checked_data_sofa = cell2table(cell(length(sofa_file),5));
    checked_data_sofa.Properties.VariableNames = {'File_Name','Repeat_Angles',...
        'Inconsist_Dist','Asymm_Angles', 'Asymm_Median'};
    % initialise output
    reverseStr = '';
    for n = 1: length(sofa_file)
        msg = sprintf(['\nchecking: ' sofa_file{n}  ' (%d/%d files)\n'], ...
            n, length(sofa_file));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        % print progress
    
        [ bad_data_temp ] = check_SOFA(sofa_file{n}, angle_range, dist_range, plot);
        % check data in the reamining sofa files
    
        checked_data_sofa(n,:)= struct2table(bad_data_temp,'AsArray',true);  
        % put result in table
    end
    % check reamining sofa files and put it into a table

    bad_data_sofa = checked_data_sofa;
    idx=all(cellfun(@isempty,bad_data_sofa{:,2:end}),2);
    bad_data_sofa(idx,:)=[];
    % remove empty row (those have nothing wrong) 
end

%% orgainise output

if exist('bad_sofa') % if ouput variable already exist
    if exist('bad_data_sofa') % if there is bad data in sofa file from input exist
        bad_sofa = [bad_sofa; bad_data_sofa];
        checked_sofa = [checked_sofa; checked_data_sofa];
    end
else % bad data in sofa file from input become output
    if exist('bad_data_sofa')
        bad_sofa = bad_data_sofa;
        checked_sofa = checked_data_sofa;
    end
end
% group ouput data

%% print warning 

warning on
warning off backtrace
if exist('bad_folder')
    warning([num2str(size(bad_folder,1)) ' problematic folder was found'])
else 
    bad_folder = [];
end
if exist('bad_sofa')
    warning([num2str(size(bad_sofa,1)) ' problematic sofa file was found'])
else
    bad_sofa = [];
end
warning on backtrace
% print warning



end

