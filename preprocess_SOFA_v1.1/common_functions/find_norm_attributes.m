function [ max_length, min_fs, min_length, max_fs, input_SOFA ] = find_norm_attributes( input )
%FIND_NORM_ATTRIBUTES Summary of this function goes here
%   Find the normalisation attributes in multiple SOFA file 
%   Main attritubes: 1. maximum hrir length, 2. minimum sampling rate


% INPUT: 
% e.g. find_norm_attributes({'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', ...
%       'MRT02.sofa', 'RIEC_hrir_subject_008.sofa', 'hrtf b_nh15.sofa', 'MRT04.sofa'});
%
% 1. input = {'ITA_HRTF_Database/SOFA'} % folder contains sofa file
%    or
%    input = {'hrtf b_nh15.sofa'} % sofa file
%    - cell array with folder directory or sofa file (could mix)
%    - input should be cell array (although the function could convert
%      character into to cell array) 
%    - example:
%      'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', 'CIPIC_hrtf_database 2/sofa' ...
%      'IRCAM_Listen_hrtf_database/hrtf', 'RIEC_HRTF_Database/hrtf', 'MRT02.sofa', ...
%      'SADIE_HRTF_Database/hrtf', 'RIEC_hrir_subject_008.sofa', 'hrtf b_nh15.sofa
%
%

% OUTPUT:
% 1. max_length: maximum hrir length from input SOFA files
%    - Recommend: using the longer hrir length won't remove any hrir 
%      details, and the shorter hrir will be zero padded when normalise.    
% 
% 2. min_fs: minimum sampling rate from input SOFA files
%    - Recommend: using the lower sampling rate will not introduce any 
%      artifical samples when normalise, although it will remove some data
%      in the higher sample rate file
%
% 3. min_length: minimum hrir length from input SOFA files
%    - just an extra option, not recommend to use
%
% 4. max_fs: maximum sampling rate from input SOFA files 
%    - just an extra option, not recommend to use
%
% 5. input_SOFA: list of input SOFA file (just for checking)
%
%

%% pre-process input

if ischar(input)
    input = {input};
end
% catch if input is char instead of cell array

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


input_SOFA = [sofa_file_in_folder; sofa_file];


%% find min max

warning off % temperaly switch off warning

length_array = zeros(length(input_SOFA), 1);
fs_array = zeros(length(input_SOFA), 1);
for n = 1: length(input_SOFA)
    SOFA_hrtf = SOFAload(input_SOFA{n});
    length_array(n) = size(SOFA_hrtf.Data.IR, 3);
    fs_array(n) = SOFA_hrtf.Data.SamplingRate;
end

[max_length, length_idx] = max(length_array);
[min_fs, fs_idx] = min(fs_array);

disp(['the maximum hrir length is ' num2str(max_length) ' in ' input_SOFA{length_idx}])
disp(['the minimum sampling rate is ' num2str(min_fs) ' in ' input_SOFA{fs_idx}])

min_length = min(length_array);
max_fs = max(fs_array);


warning on

end

