function [ Obj, SOFAfilename ] = modify_sofa_measurement(sofa_name, measurement_idx, ...
    keep_or_remove, SOFAdir, SOFAfilename)
%REMOVE_SOFA_MEASUREMENT Summary of this function goes here
%   Modify SOFA file by removing / keep part of the measurements in the SOFA file 

% INPUT: 
% e.g. modify_sofa_measurement('irc_1004.sofa', [1 3 4 5 18], ...
%           'remove', 'good_SOFA/', [])
%
% 1. sofa_name = 'irc_1004.sofa'; % original SOFA file
%
% 2. measurement_idx = [1 3 4 5 18]; % target meausrements row number
%    - if empty or []: keeps all measurements, so it will works like a
%       rename and dulplicate SOFA file function.
%
% 3. keep_or_remove = 'remove'; % decide whether should keep or remove the
%       meansurement_idx
%    - 'keep' (or anything start with 'k' or 'K'or 0 or empty): 
%       keep only meausrements in 'measurement_idx'.
%       (Default choice)
%    - 'remove' (or anything start with 'r' or 'R'or 'n' or 'N' or 1): 
%       remove the measurements in 'measurement_idx' .
% 4. SOFAdir = 'normalised_SOFA/' % new file save location (OPTIONAL)
%    - in char
%    - make sure the folder exist
%    - add / in the end
%    - save in current folder if empty or leave blank
%    
% 5. SOFAfilename = []; or desire file name e.g. 'new_file.sofa'; (OPTIONAL)
%    - in char
%    - add .sofa in the end
%    - save in original input file name + '_normalised.sofa' if empty or leave blank
%    - 0 : no saving SOFA file
%
%

% OUTPUT:
% Obj: SOFA struct (same as SOFAload('new_SOFA_file.sofa')
%
%


%% pre-press inputs and catch missing inputs

if nargin == 1
    measurement_idx = [];
    keep_or_remove = 'remove';
    SOFAdir = [];
    SOFAfilename = [sofa_name(1:end-5) '_modified.sofa'];
elseif nargin == 2
    keep_or_remove = 'keep';
    SOFAdir = [];
    SOFAfilename = [sofa_name(1:end-5) '_modified.sofa'];
elseif nargin == 3
    SOFAdir = [];
    SOFAfilename = [sofa_name(1:end-5) '_modified.sofa'];
elseif nargin == 4
    SOFAfilename = [sofa_name(1:end-5) '_modified.sofa'];
end
% catch missing input

if isempty(keep_or_remove)
    keep_or_remove = 'keep';
end
% catch empty inputs (keep_or_remove)

if isempty(measurement_idx)
    keep_or_remove = 'remove';
    warning('No measurements was removed inside the SOFA file')
end
% catch empty inputs (measurement_idx) 

if isempty(SOFAfilename)
    SOFAfilename = [sofa_name(1:end-5) '_modified.sofa'];
end
% catch empty inputs (SOFAfilename)

%% actual function start here 

if ~isempty(keep_or_remove) % if keep_or_remove variable is not empty
    if keep_or_remove(1) == 'r' || keep_or_remove(1) == 'R' || keep_or_remove(1) == 1 || ...
            keep_or_remove(1) == 'n' || keep_or_remove(1) == 'N'
    % if user wants to remove data (input starts with 'r' or 'R' or 1 or 'n' or 'N'
        sofa_hrtf = SOFAload(sofa_name);
        % load SOFA file
        SOFA_measuremnts = size(sofa_hrtf.Data.IR, 1);
        % find out the number of measurement inside the file
        selected_data = linspace(1, SOFA_measuremnts, SOFA_measuremnts)';
        % create an array that equal the original input measurements
        selected_data(measurement_idx) = [];
        % remove the measurements user wants to remove
        size(selected_data)
    else 
    % if using wants to keep selected data (default)
        sofa_hrtf = SOFAload(sofa_name);
        % find out the number of measurement inside the file
        selected_data = measurement_idx;
        % selected the data user what to keep
    end
else
    % default to keep selected data
    sofa_hrtf = SOFAload(sofa_name);
    % find out the number of measurement inside the file
    selected_data = measurement_idx;
    % selected the data user what to keep  
end

[ Obj ] = save_SOFA(sofa_name, squeeze(sofa_hrtf.Data.IR(selected_data, 1, :)), ...
	squeeze(sofa_hrtf.Data.IR(selected_data, 2, :)), size(sofa_hrtf.Data.IR, 3), ...
	sofa_hrtf.Data.SamplingRate, selected_data, ...
	SOFAdir, SOFAfilename);
% save modified SOFA file in desire location with desire name


end



