function [feature_array, scaled_feature_array, label_file, label_array, label_vector] = ...
    preprocess_data(folder, desire_feature, sample_method, input_1, input_2)
% PREPROCESS_DATA Summary of this function goes here

% Pre-process sofa files in side a folder by extracting features and labels

% INPUT
%
% 1. folder = 'RIEC_HRTF_Database/hrtf';
%    or 
%    {'RIEC_HRTF_Database/hrtf', 'ITA_HRTF_Database/SOFA'}
%    - e.g. 'ITA_HRTF_Database/SOFA', 'IRCAM_Listen_hrtf_database/hrtf'
%
% 2. desire_feature = {'azi' 'ele' 'hrir' 'octave_mean_dB' 'hrtf' 'P1_freq'};
%    - blank = {'azi' 'ele' 'distance' 'P1_freq' 'N1_freq' 'N2_freq' ...
%      'P1_N1_amp_diff' 'P1_N2_amp_diff' 'hrir' 'hrtf' 'log_hrtf' 'ILD' 'ILD_dB'...
%      'third_octave_mean' 'third_octave_mean_dB' 'octave_mean' 'octave_mean_dB'};
%      % all features
%
% 3. sample_method = 'random' or 'range' or 'all' % set sample method
%
%    3a.if sample_method == 'all' or blank 
%        % process all samples
%    
%    3b.if sample_method == 'random'
%        4. input_1 = 3; % number of sample to generate (blank for suffle all sample)
%        5. input_2 = 6; % random seed (blank = 'suffle' mode)
%
%    3c.if sample_method == 'range'
%        4. input_1 = 3; % start_sample (blank for all samples) 
%        5. input_2 = 6; % end_sample (blank for all samples from start sample) 
%

% OUTPUT
%
% feature_array % feature matrix,
% row: samples (sofa file * hrtf measurement)
% column: feature vector (e.g. [azi, ele, hrir(1,1), hrir(2,1)...])
% note that for 2*n matrix (two channels, left & right) row 1 (left) will
% go first
%
% scaled_feature_array % scaled feature matrix 
% each feature scaled by (x - mean) / standard deviation of that feature
%
% label_file % sofa file name 
% row 1: target sample 1, row 2: target sample 2, etc.
% column 1: file name, column 2: number from the folder
%
% label_target % logical matrix of each sample to identify the file name
%
% label_array % label_target in numeric vector
%

tic
if ischar(folder)
    folder = {folder};
end
sofa_file = cell(1,1);
for n = 1:length(folder)
    file_names = dir (folder{n});
    sofa_file = [sofa_file; {file_names([file_names(:).isdir]==0).name}'];
end
sofa_file = sofa_file(2:end);
% get all file names inside the folder

if nargin == 1
	desire_feature = {'azi' 'ele' 'distance' 'P1_freq' 'N1_freq' 'N2_freq' ...
        'P1_N1_amp_diff' 'P1_N2_amp_diff' 'hrir' 'hrtf' 'log_hrtf' 'ILD' 'ILD_dB'...
        'third_octave_mean' 'third_octave_mean_dB' 'octave_mean' 'octave_mean_dB'};
elseif nargin == 2
    sample_method = 'all';
	input_1 = 1;
	input_2 = length(sofa_file);
end


if strcmpi(sample_method,'random') == 1
    if nargin == 3
        input_1 = length(sofa_file);
    elseif nargin == 4
        input_2 = 'shuffle';
    end
    number_of_sample = input_1;
    random_seed = input_2;
    rng(random_seed);
    sample = randperm(length(sofa_file), number_of_sample);
    % randomly pick sample numbers
    sofa_file = sofa_file(sample);
    start_sample = 1;
    end_sample = length(sofa_file);
    
elseif strcmpi(sample_method,'range') == 1
    if nargin == 3
        input_1 = 1;
    elseif nargin == 4
        input_2 = length(sofa_file);
    end
    
    start_sample = input_1;
    end_sample = input_2;
    
    if start_sample > end_sample 
    error('start_sample must be larger than end_sample')
    % catch error if the start sample is smaller than end sample
    end
    
else
    start_sample = 1;
    end_sample = length(sofa_file);
    
end


%% initialise output

hrtf = SOFAload(sofa_file{1});
feature = feature_extraction(hrtf, 1, 0, desire_feature);
% load and first sofa file and find the features

no_of_file = length(sofa_file);
hrtf_samples = size(hrtf.Data.IR, 1);
% get the number of total files and measurement 

feature_out = cell(length(desire_feature), 1);
for n = 1 : length(desire_feature)
    feature_out{n} = extractfield(feature, desire_feature{n}); 
end
% extract desire features into a cell array

feature_array_test = feature_out{1};
if length(feature_out) > 1
	for j = 2:length(feature_out)
        feature_array_test = [feature_array_test feature_out{j}];
	end
end
% find the length of all features in a single vector

feature_array = zeros(no_of_file * (end_sample - start_sample +1), length(feature_array_test));
label_vector = zeros(no_of_file * (end_sample - start_sample +1), 1);
% initialise output array


%% actual feature extractions

warning('off','all')
reverseStr = '';
j = 0;

k = 1;
for n = start_sample : end_sample 
    hrtf = SOFAload(sofa_file{n}); % load SOFA file
    for m = 1:hrtf_samples
        feature = feature_extraction(hrtf, m, 0, desire_feature);
        % extract feature
        
        for i = 1:length(desire_feature)
            feature_out{i} = extractfield(feature,desire_feature{i}); 
        end
        % extract desire features into cell array
        

        feature_vector = feature_out{1};
        if length(feature_out) > 1
            for i = 2:length(feature_out)
                feature_vector = [feature_vector feature_out{i}];
            end
        end
        % put all features in a single vector
        
        feature_array(k, :) = feature_vector;
        % put the feature vector of each measurement into a matrix
        % rows = samples (measuremt * no of sofa file
        % column = feature details
        
        %feature_array(k, :) = [feature.azi, feature.ele, ...
            %feature.hrir(:,1)', feature.hrir(:,2)', ...
            %feature.octave_mean_dB(:,1)', feature.octave_mean_dB(:,2)'];
        
        label_vector(k) = n;
        % put the number of sofa file in array as a training target
        
        i = k;
        percent = 100 * i / (no_of_file * hrtf_samples); 
        if floor(percent) == j
            %fprintf('Finding match: %d/10 is done \n', j)
            %j = j + 1;
            msg = sprintf('Extracting features: %d/100 is done (in %d x %d = %d samples)', ...
                j, no_of_file, hrtf_samples, ...
                no_of_file * hrtf_samples);
            fprintf([reverseStr, msg, '\n']);
            reverseStr = repmat(sprintf('\b'), 1, length(msg) + 1);
            % calculate and print progress in percentage
            j = j + 1;
        end
        
        k = k + 1;
    end
end

warning('on','all')

%% fearure scaling
 
scaled_feature_array = feature_array;
for n = 1:size(feature_array,2)  
    feature_mean = mean(feature_array(:, n));
    feature_std = std(feature_array(:, n));
    scaled_feature = (feature_array(:, n) - feature_mean) ./ feature_std;
    scaled_feature_array(:, n) = scaled_feature;
    % scale feature by (x - mean) / standard deviation
end

%% label data
if strcmpi(sample_method,'random') == 1
    label_file = [sofa_file(start_sample : end_sample) num2cell(sample')];
else
    label_file = [sofa_file(start_sample : end_sample), ...
        num2cell(linspace(start_sample, end_sample , end_sample - start_sample + 1)')];
end
label_vector = label_vector - start_sample + 1;
label_array = full(ind2vec(label_vector')');
% convert the label to logical matrix
% 1 0 0 = 1
% 0 1 0 = 2
% 0 0 1 = 3

toc
end

