function [ summary ] = preprocess_SOFA( input_files, angle_range, dist_range, output_dir, ...
    plot_trig, target_length, target_fs, target_amplitude)
%PREPROCESS_SOFA Summary of this function goes here
%
%  MAIN SCRIPT:
%   Processing SOFA file(s) for research or machine learning purposes.
%
%  Mainly do 3 things:
%   1. find and fix problematic SOFA files.
%   2. find common angles between all of them.
%   3. normalise the SOFA file to the same length and sampling frequency
%      then save as a new file.
%
%  Pipeline:
%   input -> check SOFA file -> orginaise and fix SOFA files
%    -> find SOFA file that represent a unique angle distribution (speed reason)
%    -> find angle match in those 'unique' SOFA files -> find default normalise attributes
%    -> extract measurments, normalise and save as new SOFA files -> output summary

%  Example: (Declare input variables -> use the function)
%   input_files = {'SADIE_HRTF_Database/hrtf', 'ITA_HRTF_Database/SOFA', 'SADIE_HRTF_Database/hrtf', 'MRT04.sofa', ...
%       'hrtf b_nh15.sofa', 'RIEC_hrir_subject_008.sofa', 'IRCAM_Listen_hrtf_database/hrtf'};
%   angle_range = 3;
%   dist_range = 0.2;
%   output_dir = 'normalised_SOFA/';
%   % Declare input variables
%
%   [ summary ] = preprocess_SOFA( input_files, angle_range, dist_range, output_dir);
%

%  INPUT: 
% 1. input = {'ITA_HRTF_Database/SOFA'};  % folder contains sofa file
%    or
%    input = {'hrtf b_nh15.sofa'};  % sofa file
%    - cell array with folder directory or sofa file (could mix)
%    - input should be cell array (although the function could convert
%      character into to cell array) 
%    - example:
%      'ITA_HRTF_Database/SOFA', 'ARI_hrtf_database/hrtf', 'CIPIC_hrtf_database 2/sofa' ...
%      'IRCAM_Listen_hrtf_database/hrtf', 'RIEC_HRTF_Database/hrtf', 'SADIE_HRTF_Database/hrtf' ...
%      'MRT02.sofa', 'RIEC_hrir_subject_008.sofa', 'hrtf b_nh15.sofa
%
% 2. angle_range = 0.01;  % set angle matching tolerance (in degree)
%    (default angle_range = 0)
%    - this value will be used for:
%       1. check SOFA - make sure all the measurement angles in the SOFA files
%           from the same folder are within the range 
%       2. find represent SOFA - find the SOFA files that share the same
%           measurement distribution with in the range 
%       3. matching angles -  find common angles between different SOFA
%           files within the range
%    - Note: this value could be updated in matching angles process when
%    there is not matching angles could be found between 2 or more SOFA file(s)
%
% 3. dist_range = 0.1;  % set distance tolerance (in meters)  
%    (default dist_range = 0)
%    - is will only use for check SOFA, make sure all meansurements where
%       made within the range of distance.
%
% 4. output_dir = 'normalised_SOFA/';  % set output folder location
%     (Optional but highly recommanded)
%    - input the folder where to save the new SOFA file
%    - if the folder does not exsist, it will create the folder 
%    - if this input is empty, it will save at the current folder
%
%  Optional input:
% 5. plot_trig = 1;  % plot matched angles trigger (1 or 0, default = 1 (on))
%    - (default) when plot_trig = 1 it will plot the result from the matching angle
%      to show the matched angles' distribution 
%    - user can input plot_trig = 0 to supress the plot
% 
% 6. target_length = 512;  % normalise hrir target length
%    - (default) will find the maximum hrir length in all inputs SOFA that  
%      could be found in the 'find default normalise attributes' process
%      (to avoid removing hrir tails)
%
% 7. target_fs = 44100;  % normalise hrir target length
%    - (default) will find the lowest hrir sampling frequency in all inputs SOFA   
%      that could be found in the 'find default normalise attributes' process
%      (to avoid up sampling that may produce artificial data)
%
% 8. target_fs = 0.99;  % normalise hrir maximum amplitude
%    - (default) target_fs = 0.99
%      (to avoid up clipping)
%

%  OUTPUT:
% 1. normalised SOFA files with matched angles in the designated folder
%    with a output_file_log list out the angles and modification of each
%    SOFA file
% 
% 2. summary  % summerise key variables in the process 
%    (a copy will save in the folder in .mat file)
% 


%% catch missing input(s);

if nargin == 1
	angle_range = 0;
    dist_range = 0;
    output_dir = [];
    plot_trig = 1;
    target_length = [];
    target_fs = [];
    target_amplitude = [];
elseif nargin == 2
    dist_range = 0;
    output_dir = [];
    plot_trig = 1;
    target_length = [];
    target_fs = [];
    target_amplitude = [];
elseif nargin == 3
    output_dir = [];
    plot_trig = 1;
    target_length = [];
    target_fs = [];
    target_amplitude = [];
elseif nargin == 4 
    plot_trig = 1;
    target_length = [];
    target_fs = [];
    target_amplitude = [];
elseif nargin == 5 
    target_length = [];
    target_fs = [];
    target_amplitude = [];
elseif nargin == 6 
    target_fs = [];
    target_amplitude = [];
elseif nargin == 7 
    target_amplitude = [];
end
% catch missing inputs


if isempty(angle_range)
    angle_range = 0;
elseif ~isnumeric(angle_range)
    error('input "angle_range" have to be numeric.')
end
if isempty(dist_range)
    dist_range = 0;
elseif ~isnumeric(dist_range)
    error('input "dist_range" have to be numeric.')
end
if isempty(plot_trig)
    plot_trig = 1;
elseif sum(plot_trig == 'Y') > 0 || sum(plot_trig == 'y') > 0
    plot_trig = 1;
elseif sum(plot_trig == 'N') > 0 || sum(plot_trig == 'n') > 0
    plot_trig = 0;
elseif ~isnumeric(plot_trig)
    error('input "plot_trig" have to be numeric.')
end
if isempty(target_length)
    target_length = [];
elseif ~isnumeric(target_length)
    error('input "target_length" have to be numeric.')
end
if isempty(target_fs)
    target_fs = [];
elseif ~isnumeric(target_fs)
    error('input "target_fs" have to be numeric.')
end
if isempty(target_amplitude)
    target_amplitude = [];
elseif ~isnumeric(target_amplitude)
    error('input "target_amplitude" have to be numeric.')
end
% catch empty or odd inputs


if ischar(output_dir)
    if ~strcmpi(output_dir(end), '/')
        output_dir = [output_dir '/'];
    end
end
% catch if / is missing in the end of SOFAdir

if ~isempty(output_dir)
    if ~exist(output_dir, 'dir') % check if output folder exist
        mkdir(output_dir) 
        % if output folder did not exist, create output folder 
        warning([' created output folder ' output_dir(1:end-1) '.'])
    else
        reenter = [];
        while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
        sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
        isempty(reenter)
            reenter = input(['\n ' output_dir(1:end-1) ' folder already exist, '...
                'overright folder? (y/n) \n '], 's');
            if ~isempty(reenter)
                if prod(reenter ~= 'Y') ~= 0 && prod(reenter ~= 'y') ~= 0
                    error([output_dir(1:end-1) ' folder already exist, '...
                        'please remove folder'])
                    % if user did not choose yes (y or Y or anything include 'y')
                    % this fuction will abort and show an error
                else
                    rmdir(output_dir(1:end-1), 's')
                end   
            end
        end
    end
end
% check if the output folder avalible 

%% check SOFA files

warning off % temperaly switch off warning

[ bad_folder, bad_sofa, checked_sofa ] = check_hrtf(input_files, angle_range, dist_range, 0);
 
 
 %% organise SOFA files (fix problematic SOFA file and put all usable SOFA into a new folder)

[ good_hrtf_log ] = group_and_fix_SOFA( bad_folder, bad_sofa, checked_sofa );
 
disp(' Press a key to continue')      % the lower left corner of MATLAB window.
pause;
 
 %% find SOFA file with unique angles (so no need to go through all file in the next step)
 
[ unique_sofa_table ] = find_represent_SOFA ( 'good_SOFA', angle_range, 1 );
 
%% match angle

warning on all
warning on backtrace
% switch on warning

matched_angles = [];
% reset output

sofa_match_cell = table2cell(unique_sofa_table);
matching_angle_range = angle_range;
try
    warning off
    [matched_angles] = common_angle(sofa_match_cell(:,1), matching_angle_range, plot_trig);
    % try to find matched angles was found between previous preprocessed SOFA file table
    warning on
catch ME
    % catch if there is error occured
    if (strcmp(ME.identifier,'common_angle:no2Match'))
        % if there are no matches between 2 spacific SOFA files

        warning(ME.message)
        % print warning message
        
        error_message = ME.message;
        % save message to extract the related file names later
        
        no_match_SOFA = cell(2,1);
        % initialise error file name array
        
        name_start_idx = strfind(error_message,'between') + 8;
        name_end_idx = strfind(error_message,'&') - 2;
        % find the loction of the file name by finding spacific word or symbol 
        no_match_SOFA(1) = {error_message(name_start_idx:name_end_idx)};
        % save the file name in the error file name array

        name_start_idx = strfind(error_message,'&') + 2;
        name_end_idx = strfind(error_message,',') - 1;
        % find the loction of the file name by finding spacific word or symbol 
        no_match_SOFA(2) = {error_message(name_start_idx:name_end_idx)}; 
        % save the file name in the error file name array

        
        fprintf(' Problematic SOFA file (or representative file): ')
        fprintf('\n - %s ',no_match_SOFA{:})
        fprintf('\n')
        % print out the problematic file names

    elseif (strcmp(ME.identifier,'common_angle:noMatch'))
        % if there are no matches was found between different SOFA files
        
        warning(ME.message)
        % print warning message
        
        error_message = ME.message;
        % save message to extract the related file names later
        
        comma_count = strfind(error_message,',');
        % count the comma to find of the number of error files
        % (because in this case 1 comma = 1 error file)
        no_match_SOFA = cell(length(comma_count),1);
        % initialise error file name array by the number just found
        
        name_start_idx = strfind(error_message,'in');
        name_start_idx = name_start_idx(2) + 3;
        name_end_idx = comma_count(1) - 1;
        % find the loction of the first file name by finding spacific word or symbol 
        no_match_SOFA(1) = {error_message(name_start_idx:name_end_idx)};
        % save the first file name in the error file name array
        
        for n = 2:length(comma_count)
            
            name_start_idx = comma_count(n-1) + 2;
            name_end_idx = comma_count(n) - 1;
            % find the loction of the file name by finding spacific word or symbol 
            no_match_SOFA(n) = {error_message(name_start_idx:name_end_idx)};
            % save the first file name in the error file name array

        end
         
        fprintf(' Problematic SOFA file (or representative file): ')
        fprintf('\n - %s ',no_match_SOFA{:})
        fprintf('\n')
        % print out the problematic file names
    else
        error(ME.message)
        % print error message
    end
end

while ~exist('matched_angles', 'var') || isempty(matched_angles)
    % is there is no final result (something is wrong) 
    
    reenter = []; % initialise user input variable
    while sum(reenter == 'R') < 1 && sum(reenter == 'r') < 1 && ...
        isnan(str2double(reenter)) || isempty(reenter)
    
%         sum(reenter == 'R') < 1 && sum(reenter == 'r') < 1 && ...
%         sum(reenter == 'S') < 1 && sum(reenter == 's') < 1 && ...
%         isnan(str2double(reenter)) || isempty(reenter)
    
        % 3 ways to fix the SOFA file with no match:
        % - 1. increase range (by input a new number(in degree))
        % - 2. remove one or mulitple file(s) (by input 'r') 
        %      - the file(s) will be chosen later
%        % - 3. split (split to find a match with each file
%        %      - takes lots of time and the result may not be accurate 
        
        reenter = input(['\n try with a wider range or remove one file ? \n'...
            ' (new range (in number)/ r(remove, will select the file later)) : \n '], 's');
%          reenter = input(['\n try with a wider range or remove one file or split into '...
%             num2str(size(no_match_SOFA, 1)) ' ? \n'...
%             ' (new range (in number)/r(remove, will select the file later)/ s(split)) : \n '],...
%             's');
        % user option to solve the no match problem   
    end
    
	if ~isnan(str2double(reenter))
        % if input is a number, that means this is am update in angle range
        
        matching_angle_range = str2double(reenter);
        % set new angle range
        warning(['changed the angle range to ' reenter ])

	elseif ~isempty(reenter) 
        if sum(reenter == 'R') > 0 || sum(reenter == 'r') > 0
            % if input = r (or include r), that's mean remove file
            
            reenter = []; % reset reenter
            while isempty(reenter) % if reenter is empty (no input)
%                 isnan(str2double(reenter)) || isempty(reenter)...
%                     || str2double(reenter) > size(no_match_SOFA, 1)
%                 
                                 
                fprintf(' Problematic SOFA file (or representative file): ')
                for n = 1:length(no_match_SOFA)
                    fprintf(1, '\n %d - %s ', n, ...
                    no_match_SOFA{n})
                end

                reenter = input(['\n please input index in the first column to select' ...
                    'file you want to remove. \n'],'s');
                % print out the no matched SOFA file names again
                % but this time with a index in the first column
                % user need to input the index to choose which file to delete
                
                remove_idx = str2double(regexp(reenter,'\d*','match')');
                % convert string to double, for more than one numbers.
                
                if max(remove_idx) > size(no_match_SOFA, 1) || min(remove_idx) < 1
                    % make sure user input is in the right range
                    warning('Input Number exceed the maximum number in the First Column.')
                    % print warning if user's input excceded the proper range
                    reenter = []; % reset reenter
                elseif isempty(remove_idx)
                    warning('Input must be number(s). ')
                    % print warning if no number was found in the input
                    reenter = []; % reset reenter
                else
                 
                    sofa_match_cell_old = sofa_match_cell;
                    % duplicate input table for angle match
                    for n = 1: length(remove_idx)
                        matched_idx = strfind(sofa_match_cell(:,1), ...
                            no_match_SOFA(remove_idx(n)));
                        idx = not(cellfun('isempty', matched_idx));
                        sofa_match_cell(idx,:) = [];
                    end
                    % from user's input, find the name of the file to remove
                    % from new_unique_sofa_table.
                    % (dulplicate from unique_sofa_table)
            
                    warning('removed selected files')
                    reenter = 'r';
                end
            end
            
            
%         elseif sum(reenter == 'S') > 0 || sum(reenter == 's') > 0
%             
%             temp_sofa_cell = sofa_match_cell;
%             for n = 1 :length(no_match_SOFA)
%                 matched_idx = strfind(temp_sofa_cell(:,1), ...
%                             no_match_SOFA(n));
%                 idx = find(not(cellfun('isempty', matched_idx)));
%                 temp_sofa_cell(idx,:) = [];
%             end
%             
%             splitted_sofa_match = cell(size(temp_sofa_cell,1)+1, ...
%                 size(sofa_match_cell,2),length(no_match_SOFA));
%             
%             for n = 1 :length(no_match_SOFA)
%                 matched_idx = strfind(sofa_match_cell(:,1), ...
%                             no_match_SOFA(n));
%                 idx = find(not(cellfun('isempty', matched_idx)));
%                 splitted_sofa_match(1,:,n) = sofa_match_cell(idx,:);
%                 splitted_sofa_match(2:end,:,n) = temp_sofa_cell;
%             end
%             
%             for n = 1 : size(splitted_sofa_match, 3)
%                 common_angle(splitted_sofa_match(:,1,n), angle_range, 1);
%             end
%             
%             warning('splitted the no matched SOFA file')
%             reenter = []; 
        else
            warning('please make sure the input is a number / r (remove)')
%             warning(['please make sure the input is a number / r (remove)' ...
%                 ' / s (split)'])
        end
	else
	error('please make sure the input in not empty')
	% if there is no input this fuction will abort and show an error
	end
    
    try
        warning off
        [matched_angles] = common_angle(sofa_match_cell(:,1), matching_angle_range, plot_trig);
        % try to find matched angles was found between previous preprocessed SOFA file table
        warning on
    catch ME
        % catch if there is error occured
        if (strcmp(ME.identifier,'common_angle:no2Match'))
            % if there are no matches between 2 spacific SOFA files

            warning(ME.message)
            % print warning message

            error_message = ME.message;
            % save message to extract the related file names later

            no_match_SOFA = cell(2,1);
            % initialise error file name array

            name_start_idx = strfind(error_message,'between') + 8;
            name_end_idx = strfind(error_message,'&') - 2;
            % find the loction of the file name by finding spacific word or symbol 
            no_match_SOFA(1) = {error_message(name_start_idx:name_end_idx)};
            % save the file name in the error file name array

            name_start_idx = strfind(error_message,'&') + 2;
            name_end_idx = strfind(error_message,',') - 1;
            % find the loction of the file name by finding spacific word or symbol 
            no_match_SOFA(2) = {error_message(name_start_idx:name_end_idx)}; 
            % save the file name in the error file name array


            fprintf(' Problematic SOFA file (or representative file): ')
            fprintf('\n - %s ',no_match_SOFA{:})
            fprintf('\n')
            % print out the problematic file names

        elseif (strcmp(ME.identifier,'common_angle:noMatch'))
            % if there are no matches was found between different SOFA files

            warning(ME.message)
            % print warning message

            error_message = ME.message;
            % save message to extract the related file names later

            comma_count = strfind(error_message,',');
            % count the comma to find of the number of error files
            % (because in this case 1 comma = 1 error file)
            no_match_SOFA = cell(length(comma_count),1);
            % initialise error file name array by the number just found

            name_start_idx = strfind(error_message,'in');
            name_start_idx = name_start_idx(2) + 3;
            name_end_idx = comma_count(1) - 1;
            % find the loction of the first file name by finding spacific word or symbol 
            no_match_SOFA(1) = {error_message(name_start_idx:name_end_idx)};
            % save the first file name in the error file name array

            for n = 2:length(comma_count)

                name_start_idx = comma_count(n-1) + 2;
                name_end_idx = comma_count(n) - 1;
                % find the loction of the file name by finding spacific word or symbol 
                no_match_SOFA(n) = {error_message(name_start_idx:name_end_idx)};
                % save the first file name in the error file name array

            end

            fprintf(' Problematic SOFA file (or representative file): ')
            fprintf('\n - %s ',no_match_SOFA{:})
            fprintf('\n')
            % print out the problematic file names
        else
            error(ME.message)
            % print error message
        end
    end
    
end

matched_angles = reshape(matched_angles,length(matched_angles),[]); 
matching_result = [matched_angles sofa_match_cell];

disp(' Press a key to continue') 
pause;

%% pick azi angels to keep (optional) - future development



%% find normalise attributes

matching_result_reshape = reshape(matching_result(:,2:end),[], 1);
matching_result_reshape = matching_result_reshape(~cellfun('isempty',matching_result_reshape));  

input_SOFA = matching_result(:, 2:end);
input_SOFA =input_SOFA(~cellfun('isempty',matching_result(:, 2:end)));
% extract all SOFA files 

[ max_length, min_fs, min_length, max_fs] = find_norm_attributes(input_SOFA);

%% normalise & save as SOFA

warning off
% temporally swith off warning

if isempty(target_length)
    target_length = max_length;
end
if isempty(target_fs)
    target_fs = min_fs;
end
% set default target length and target sampling frequency if empty 

total_files = sum(sum(~cellfun('isempty',matching_result(:, 2:end))), 2);
% find total number of files to save

output_file_log = cell(total_files, 3);
% initialise output log

[~, id] = unique(good_hrtf_log(:, 1));
good_hrtf_log_unique  = good_hrtf_log(id, :);
reverseStr = '';
k = 1;
% initiallise print variables
for n = 1: size(matching_result, 1)
    % for the number of 'unique' files (2nd column in the matched result)
    temp_cell = matching_result(~cellfun('isempty',matching_result(n,:)));  
    % pick out the file(s) in each row in matching_result (remove empty cell)
    for m = 2: length(temp_cell)
        % for the files in each row in matching_result (saved in 'temp_cell')
        
        msg = sprintf(['\n normalise and saving: ' matching_result{n, m} ' (%d/%d files)\n'], ...
            k, total_files);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        % print progress
        
        normalise_hrir(matching_result{n, m}, matching_result{n, 1}(:,4), ...
            target_length, target_fs, target_amplitude, output_dir, [], 0, 0); 
        % normalise and save file(s)
        
        output_file_log(k,1) = matching_result(n, m);
        output_file_log(k,2) = matching_result(n, 1);        
        idx = ismember(good_hrtf_log_unique, matching_result(n, m));
        output_file_log(k,3) = good_hrtf_log_unique(idx, 2);
        % fill in output log
        
        k = k+1;
        % print count
    end
end

addpath(output_dir)
% add folder to search path
save([output_dir 'output_file_log.mat'],'output_file_log')
% save log

warning on
% turn warning back on

%% export a summary of key variables

summary.input.input_files = input_files;
summary.input.angle_range = angle_range;
summary.input.dist_range = dist_range;
summary.input.output_dir = output_dir;

summary.input_SOFA = input_SOFA;

summary.target_length = target_length;
summary.target_fs = target_fs;

summary.check_hrtf.bad_folder = bad_folder;
summary.check_hrtf.bad_sofa = bad_sofa;
summary.check_hrtf.checked_sofa = checked_sofa;

summary.good_hrtf_log = good_hrtf_log;

summary.unique_sofa_table = unique_sofa_table;

summary.match_angle.sofa_match_cell = sofa_match_cell;
summary.match_angle.matching_angle_range = matching_angle_range;
summary.match_angle.matched_angles = matched_angles;
summary.match_angle.matching_result = matching_result;

summary.find_norm_attributes.matching_result_reshape = matching_result_reshape;
summary.find_norm_attributes.max_length = max_length;
summary.find_norm_attributes.min_fs = min_fs;
summary.find_norm_attributes.min_length = min_length;
summary.find_norm_attributes.max_fs = max_fs;

summary.output_file_log = output_file_log;
% save output summary

save([output_dir 'summary.mat'],'summary')
% save summary

end
