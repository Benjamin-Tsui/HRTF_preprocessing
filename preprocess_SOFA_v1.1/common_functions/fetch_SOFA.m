function [ SOFA_cell_array ] = fetch_SOFA(input)
%LIST_SOFA Summary of this function goes here
%   input any folder directory or sofa file in a cell array and it will
%   output 'SOFA_cell_array' which is a cell array that includes all 
%   SOFA file name(s) found in the input(s).

% Example:
%   [ SOFA_cell_array ] = fetch_SOFA({'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', ...
%       'MRT02.sofa', 'RIEC_hrir_subject_008.sofa', 'hrtf b_nh15.sofa', 'MRT04.sofa'})
%

% INPUT: 
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

% OUTPUT:
% SOFA_cell_array: single column cell array that contains all SOFA file name(s) found in the input(s).
%


if ~iscell(input)
    if ischar(input)
        input = {input};
    else
        error('input has to be character or cell array');    
    end
end
% catch if input is not cell array

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

if ischar(folder)
    folder = {folder};
end
sofa_file_in_folder = cell(1,1);
for n = 1:length(folder)
    file_names = dir (folder{n});
    sofa_file_in_folder = [sofa_file_in_folder; {file_names([file_names(:).isdir]==0).name}'];
end
sofa_file_in_folder = sofa_file_in_folder(2:end);
% get all file names inside the folder


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
% group sofa files in input (through warning if input in not a folder nor sofa file


SOFA_cell_array = [sofa_file_in_folder; sofa_file];
% group output

idx = strfind(SOFA_cell_array, '.sofa');
idx = not(cellfun('isempty', idx));
SOFA_cell_array = SOFA_cell_array(idx);
% remove non-sofa files

end

