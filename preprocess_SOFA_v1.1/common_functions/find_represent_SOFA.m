function [ unique_sofa_table ] = find_represent_SOFA ( folder, range, print )
%FIND_REPRESENT_SOFA Summary of this function goes here
% 
% Find the SOFA files with the same angle distribution 
% to save computation time when finding common angles in the next stage

% INPUT: 
% e.g. [ unique_sofa_table ] = find_represent_SOFA( 'good_SOFA', 0.5, 1 )

% 1. folder = 'good_SOFA';
%    folder name that holds all the SOFA files
%    (input in character array)
%
% 2. range = 0.5;
%    error tolerance when matching angel
%    (incase the measurement angel were recoded high precision, e.g. rounded to 0.001)
%    - default(if empty or missing) = 0
%
% 3. print = 1; 
%    print out the progress 
%    - default(if empty or missing) = 1
%

% OUTPUT:
% 1. unique_sofa_table 
%    - column 1: the name of the 'unique' SOFA file or the one that could represent 
%                the others SOFA file that has the same angle distribution 
%    - column 2: name of the 'matched' SOFA file that could be represented by the 
%                'unique' SOFA file at the first column
% 
 


if nargin == 1
	range = 0;
    print = 1;
elseif nargin == 2
    print = 1;
end
% catch missing inputs

if isempty(range)
    range = 0;
end
if isempty(print)
    print = 1;
end
% catch empty inputs

if ~ischar(folder)
    error(' the input "folder" should be a character array')
end

warning off

[input_SOFA] = fetch_SOFA(folder);
% find all SOFA file names from input

input_SOFA = unique(input_SOFA);
% incase there is duplicate file

idx = strfind(input_SOFA, '.sofa');
idx = not(cellfun('isempty', idx));
input_SOFA = input_SOFA(idx);
% remove non-sofa files

unique_sofa_log = cell(length(input_SOFA), length(input_SOFA));
% initialise output log

unique_sofa_log(1,1) = input_SOFA(1);
% save the first one as the first 'unique' SOFA file

reverseStr = '';  % for reverse print
j = 1; % intialise for the print loop later
for n = 2: length(input_SOFA)
    unique_length = length(find(~cellfun('isempty',unique_sofa_log(:,1))));
    % length of SOFA files with unqiue meausrments distribution 
    unique_sofa_log(unique_length + 1, 1) = input_SOFA(n);
    % default to assume the up coming SOFA file is unique
    
    for m = 1: unique_length    
        unique_sofa = SOFAload([folder '/' unique_sofa_log{m,1}]);
        % load the 'unqiue' SOFA file that were found previously, one by one
        unique_sofa_angles = SOFAcalculateAPV(unique_sofa);
        % load the angles from the 'unique' SOFA file
        
        new_sofa = SOFAload([folder '/' input_SOFA{n}]);
        % load the next SOFA file to compare with the 'unique' SOFA file later
        sofa_angles = SOFAcalculateAPV(new_sofa);
        % load the angles from the 'new' SOFA file
        
        if size(unique_sofa_angles, 1) == size(sofa_angles, 1) 
            % if the total measured angle btween the 'new' SOFA file and
            % the 'unique' SOFA fil matches
            if max(max(unique_sofa_angles(:, 1:2) - sofa_angles(:, 1:2))) <= range
                % compare the 'new' SOFA file and the 'unique' SOFA fil matches
                % make sure the difference between them is within range
                % and the sequences are the same
                
                matched_idx = length(find(~cellfun('isempty',unique_sofa_log(m,:))));
                % length of SOFA files that matched with the unique file
                % to locate the cell to save the 'new' SOFA file name
                
                unique_sofa_log(m, matched_idx + 1) = input_SOFA(n);
                % add the file name to the furthest column in that 'unique'
                % meausrement row
                unique_sofa_log{unique_length + 1, 1} = [];     
                % remove the file form the unique comlumn
            end
        end   
    end
        
	percent = 100 * ( n / length(input_SOFA)); 
    % percentage of the compared SOFA file
	if floor(percent) >= j
        % make sure it wont run in every time (if there are more than 100)
        % for saving computational power 
        if print ~= 0
            msg = sprintf([' Finding SOFA file with unique angle distribution: '...
                '%d/100 done (in %d SOFA file)'], floor(percent), length(input_SOFA));
            fprintf([reverseStr, msg, '\n']);
            reverseStr = repmat(sprintf('\b'), 1, length(msg) + 1);
            % print progress in percentage
        end
            j = floor(percent);
	end
end


unique_length = length(find(~cellfun('isempty',unique_sofa_log(:,1))));
% find the total number of 'unique' measurements distribution

max_matched_length = 0;
for k = 1:unique_length
    matched_idx = length(find(~cellfun('isempty',unique_sofa_log(k,:))));
    if matched_idx > max_matched_length
        max_matched_length = matched_idx;
    end
end
% find the most number of 'matched' files

unique_sofa_log = unique_sofa_log(1: unique_length, 1: max_matched_length);
% remove empty cell

unique_sofa_table = cell2table(unique_sofa_log);
% convert 'unique_sofa_log' to table
table_name = cell(1,max_matched_length);
table_name(1) = {'Unique_file'};
for n = 2:max_matched_length
    table_name(n) = {['Matched_file_' num2str(n-1)]};
end
% define table column names
unique_sofa_table.Properties.VariableNames = table_name;
% update table column names

if print ~= 0
    fprintf('\n SOFA file(s) that represent a kind of angle distribution: \n')
    disp(unique_sofa_table.Unique_file);
end

warning on

end

