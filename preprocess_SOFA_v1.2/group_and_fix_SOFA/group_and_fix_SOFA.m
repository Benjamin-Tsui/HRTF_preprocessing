function [ good_hrtf_log ] = group_and_fix_SOFA( bad_folder, bad_sofa, checked_sofa )
%GROUP_AND_FIX_SOFA Summary of this function goes here

% Follow after 'check_hrtf' function,
% The function helps user to organise the SOFA files into a new folder.
% Besides just simplily put all perfectly good SOFA files into the folder,
% The function can potentially fixs some abnormal SOFA files and folders 
% that were found by the 'check_hrtf' function.
% 
% By inputing the results from the 'check_hrtf' function,
% this function will: 
% 1. Create a new foler call 'good_SOFA' that inculdes all the good and
%    fixed SOFA files.
% 2. Besides the SOFA files, it will include a 'good_hrtf_log.mat' file,
%    which listed all the files inside the folder and notes about whether
%    the file was originally good, or how it got modified by this function.
%    (This is also the function's output)


% INPUT: 
% e.g. input_SOFA = {'MRT04.sofa', 'hrtf b_nh15.sofa', 'RIEC_hrir_subject_008.sofa'};
%      [bad_folder, bad_sofa, checked_sofa] = check_hrtf(input_SOFA, 0.01, 0.1, 0);
%      [ good_hrtf_log ] = group_and_fix_SOFA( bad_folder, bad_sofa, checked_sofa );

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

% OUTPUT:
% 1. good_hrtf_log % a log about file(s) saved in the 'good_SOFA' folder
%    - column 1: SOFA file name (name of the SOFA files in the folder)
%    - column 2: Note about the SOFA file 
%                (Was it orginally good, or fixed by this function somehow)
% 
%

%% 0.0 Create a new folder call 'good_SOFA' (or check whether it exists)

if exist('good_SOFA', 'dir') % check if 'good_SOFA' folder already exist
    
    reenter = [];
    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
        sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
        isempty(reenter)
    reenter = input('\n good_SOFA folder already exist, overright folder? (y/n) : \n ', 's');
    % user input to decide whether or not to overright existing 'good_SOFA' folder
        if ~isempty(reenter)
            if prod(reenter ~= 'Y') ~= 0 && prod(reenter ~= 'y') ~= 0
                error('folder good_SOFA already exist, please remove folder')
                % if user did not choose yes (y or Y or anything include 'y')
                % this fuction will abort and show an error
            else
                if exist('good_SOFA_old', 'dir')
                    rmdir good_SOFA_old s
                end
                movefile('good_SOFA', 'good_SOFA_old')
                mkdir good_SOFA  
                % if user choose yes (y or Y or anything include 'y')
                % this fuction will 'overwrite' the folder 
                % (read the NOTE below for more details)

                % NOTE: 
                % first time will change the current 'good_SOFA' folder to 'good_SOFA_old'
                % then re-create a new 'good_SOFA' folder
                % * however, when excute the second time 'good_SOFA_old' will be replaced
                %   and files inside will be LOST FOREVER.

            end
        end
    end
else
    mkdir good_SOFA 
    % if 'good_SOFA' folder did not exist, create 'good_SOFA' folder 
end
rmpath('good_SOFA')
% remove 'good_SOFA' folder from search path 
% (avoid dulipate file in following process)

warning off


%% 1.0 Copy good SOFA file into the new 'good_SOFA' folder

good_hrtf_log = {'init', 'init'};
% initialise good_hrtf_log
for n = 1: size(checked_sofa,1)
    if sum(cellfun('isempty',checked_sofa{n, 2:5})) == 4
        copyfile(which(char(checked_sofa{n,1})), 'good_SOFA');
        % copy perfectly good SOFA file into 'good_SOFA' folder
        good_hrtf_log = [good_hrtf_log ; checked_sofa{n,1}, 'good'];
        % update 'good_hrtf_log' with all the good SOFA file names
        % add noted as 'good' in second column
    end
end

good_hrtf_log = good_hrtf_log(2:end,:);
[~,unique_idx]=unique(good_hrtf_log(:,1),'rows');
good_hrtf_log =  good_hrtf_log(unique_idx,:);
% remove initialise and duplicate cell

reenter = []; % initialise user input variable
if ~isempty(bad_folder) || ~isempty(bad_sofa)
    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
        sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || isempty(reenter)
        reenter = input([' Abnormal SOFA file(s) was found, '...
            'want to fix abnormal SOFA file(s)? (y/n) \n '], 's');
    end


    if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
            sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 
            % double check input
        if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0


    %% 2.0 FIX BAD FOLDER ERROR (missing measurement, abnormal HRIR length, abnormal smaple rate)

            if ~isempty(bad_folder) % if there is bad folder
                for n = 1: size(bad_folder,1) % for the folders contain bad measurements
                    if sum(~cellfun('isempty',bad_folder{n, 2:5})) > 0 
                    % double check if all of them have problem

                        bad_folder_struct  = table2struct(bad_folder, 'ToScalar', true);
                        % convert bad_folder from table to struct


        %% 2.1 for Missing Measurement
                        if ~cellfun(@isempty, bad_folder_struct.Measurement) 
                            % FIRSTLY, deal with the ones has missing measurements

                           reenter = []; % initialise user input variable
                           while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                                sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
                                isempty(reenter)
                                reenter = input(['\n There is/ are Missing Measurement(s) ' ...
                                    'in one/some SOFA file(s). \n '...
                                    '(compare with other SOFA files inside the folder ' ...
                                    bad_folder_struct.Folder{n} ') \n '...
                                    'Keep SOFA file with Missing Measurement? (y/n): \n '], 's');
                                    % request user input whether they want to keep data with 
                                    % missing measurements
                            end

                            if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                                    sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 
                                    % double check input
                                if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                                    % if user choose yes (y or Y or anything include 'y')
                                    % this fuction will add the problematic files to 
                                    % 'good_SOFA' folder

                                    warning on
                                    warning('kept data with missing measurement')
                                    warning off
                                    % print warning

                                    problem_Measurement = bad_folder_struct.Measurement(~cellfun(...
                                        'isempty',bad_folder_struct.Measurement));
                                    problem_Measurement = vertcat(problem_Measurement{:});
                                    % create new problematic file only cell array 
                                    for m = 1: length(problem_Measurement)
                                        copyfile(which(char(problem_Measurement{m})), 'good_SOFA');
                                        % copy the problematic SOFA file into 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            problem_Measurement(m), 'missing_measuremnt'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA
                                        % file names, add 'missing_measuremnt' in second column
                                    end     
                                end
                            else
                                error('input must be y or n.')
                            end       
                        end

        %% 2.2 for Abnormal HRIR Length

                        if ~cellfun(@isempty, bad_folder_struct.HRIR_Length)
                            % SECONDLY, deal with the ones with different HRIR_Length

                             reenter = []; % initialise user input variable
                             while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                                 sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
                                 isempty(reenter)
                                reenter = input(['\n There is/are SOFA file(s) with ' ... 
                                    'abnormal HRIR length. \n '...
                                    '(compare with other SOFA files inside the same folder) \n '...
                                    'Keep SOFA file with abnormal HRIR length? (y/n): \n '], 's');
                                    % request user input whether they want to keep data with
                                    % abnormal HRIR length
                             end

                             if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                                    sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 
                                    % double check input
                                if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                                    % if user choose yes (y or Y or anything include 'y')
                                    % this fuction will add the problematic files to 
                                    % 'good_SOFA' folder

                                    warning on
                                    warning('kept data with missing measurement')
                                    warning off
                                    % print warning

                                    problem_HRIR_Length = bad_folder_struct.HRIR_Length(~cellfun(...
                                        'isempty', bad_folder_struct.HRIR_Length));
                                    % create new problematic file only cell array 
                                    for m = 1: length(problem_HRIR_Length)
                                        copyfile(which(char(problem_HRIR_Length{m})), 'good_SOFA');
                                        % copy the problematic SOFA file into 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            problem_HRIR_Length(m), 'abnormal_HRIR_length'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA 
                                        % file names add 'abnormal_HRIR_length' in second column
                                    end     
                                end
                            else
                                error('input must be y or n.')
                            end          
                        end


        %% 2.3 for Abnormal Smaple Rate

                        if ~cellfun(@isempty, bad_folder_struct.Sampling_Rate)
                            % THIRDLY, deal with the ones with different sampling rate       

                            reenter = []; % initialise user input variable
                            while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                                 sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
                                 isempty(reenter)
                                reenter = input(['\n There is/are SOFA file(s) with ' ...
                                    'abnormal Sampling Rate. \n '...
                                    '(compare with other SOFA files inside the same folder) \n '...
                                    'Keep SOFA file with abnormal Sampling Rate? (y/n): \n '], 's');
                                    % request user input whether they want to keep data with
                                    % abnormal Sampling Rate
                            end

                            if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                                sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 
                                % double check input
                                if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                                    % if user choose yes (y or Y or anything include 'y')
                                    % this fuction will add the problematic files to 
                                    % 'good_SOFA' folder

                                    warning on
                                    warning('kept data with missing measurement')
                                    warning off
                                    % print warning

                                    problem_sampling_rate = bad_folder_struct.Sampling_Rate( ...
                                        ~cellfun('isempty', bad_folder_struct.Sampling_Rate));
                                    % create new problematic file only cell array 
                                    for m = 1: length(problem_sampling_rate)
                                        copyfile(which(char(problem_sampling_rate{m})), 'good_SOFA');
                                        % copy the problematic SOFA file into 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            problem_sampling_rate(m), 'abnormal_Sampling_Rate'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA 
                                        % file names, add 'abnormal_Sampling_Rate' in second column
                                    end     
                                end
                            else
                                error('input must be y or n.')
                            end          
                        end       
                    end
                end
            end



    %%  3.0 FIX SOFA FILE ERROR (Duplicate angle, Inconsistent distance, 
     %  asymmetry angle (left vs right), asymmetry median (front (0) vs back (180))

            if ~isempty(bad_sofa) % if there is bad sofa file

                bad_sofa_struct = table2struct(bad_sofa, 'ToScalar', true);
                % convert bad_folder from table to struct

        %% 3.1 for Duplicate Angle(s)

              %%% FIRSTLY, deal with the SOFA files has duplicate measurements
                Repeat_Angles_idx = find(~cellfun(@isempty, bad_sofa_struct.Repeat_Angles));
                % find the files that has duplicated angles
                sofa_Repeat_Angles = bad_sofa_struct.File_Name(Repeat_Angles_idx);
                % get the file names of those files into an array
                [sofa_Repeat_Angles, unique_idx] = unique(sofa_Repeat_Angles);
                % find the unique file names and location inside the file name array
                % just created
                Repeat_Angles_idx = Repeat_Angles_idx(unique_idx);
                % remove duplicated file location in the file name array

                if ~isempty(Repeat_Angles_idx)

                    reenter = []; % initialise user input variable
                    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                         sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                         sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 || isempty(reenter)
                         reenter = input(['\n There is/are SOFA file(s) with '...
                             'Duplicate Angles. \n '...
                            'Keep SOFA file(s) by removing Repeated Angles? (y/n/list): \n '], 's');
                        % request user input whether they want to keep data by removing 
                        % repeated angles
                    end

                    if sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                       % if user types list (l or L or anything include 'l')
                       % this fuction will list out the name of the problematic SOFA file(s)

                        reenter = []; % initialise user input variable
                        while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                              sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || isempty(reenter)

                          fprintf(1, '\n - %s ', bad_sofa_struct.File_Name{Repeat_Angles_idx})
                          % print out problematic file(s)

                          reenter = input(['\n \n Above is/are the SOFA file(s) with '... 
                            'Duplicate Angles. \n '...
                            'Keep these SOFA file(s) by removing Repeated Angles? (y/n): \n '],'s');
                           % request user input whether they want to keep data by removing 
                           % repeated angles
                           % (after they have read the list of problematic SOFA files)
                        end
                    end       

                    if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                            sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 % double check input
                        if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                        % if user choose yes (y or Y or anything include 'y')
                        % this fuction will add the problematic files to 'good_SOFA' folder

                            warning on
                            warning('kept SOFA file by removing Repeated Angles')
                            warning off
                            % print warning

                            for n = 1:length(Repeat_Angles_idx)
                                [~, new_name] = modify_sofa_measurement(...
                                    bad_sofa_struct.File_Name{Repeat_Angles_idx(n)}, ...
                                    bad_sofa_struct.Repeat_Angles{Repeat_Angles_idx(n)}, ...
                                    'remove', 'good_SOFA/', []);
                                % modify, rename and copy the problematic SOFA file into 
                                % 'good_SOFA' folder
                                good_hrtf_log = [good_hrtf_log ; ...
                                    {new_name}, 'removed Repeated Angles'];
                                % update 'good_hrtf_log' by adding the problematic SOFA file names
                                % add noted as 'removed Repeated Angles' in second column
                            end
                        end
                    else
                        error('input must be y or n.')
                    end
                end

        %% 3.2 for Inconsistent Distance

              %%% SECONDLY, deal with the SOFA files has Inconsistent Distance
                Inconsist_Dist_idx = find(~cellfun(@isempty, bad_sofa_struct.Inconsist_Dist));
                % find the files that has Inconsistent Distance
                sofa_Inconsist_Dist = bad_sofa_struct.File_Name(Inconsist_Dist_idx);
                % get the file names of those files into an array
                [sofa_Inconsist_Dist, unique_idx] = unique(sofa_Inconsist_Dist);
                % find the unique file names and location inside the file name array
                % just create
                Inconsist_Dist_idx = Inconsist_Dist_idx(unique_idx);
                % remove duplicated file location in the file name array

                if ~isempty(Inconsist_Dist_idx)

                    reenter = []; % initialise user input variable
                    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                         sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                         sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 || isempty(reenter)
                        reenter = input(['\n There is/are SOFA file(s) with Inconsist '...
                            'Distance(s). \n '...
                            'Keep SOFA file(s) by removing measurements with Inconsist '... 
                            'Distance(s)? (y/n/list): \n '], 's');
                    end
                    % request user input whether they want to keep data by removing repeated angles

                    if sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                       % if user types list (l or L or anything include 'l')
                       % this fuction will list out the name of the problematic SOFA file(s)

                        reenter = []; % initialise user input variable
                        while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                              sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
                              isempty(reenter)


                            fprintf(1, '\n - %s ', bad_sofa_struct.File_Name{Inconsist_Dist_idx})
                            % print out problematic file(s)
                            reenter= input(['\n \n Above is/are the SOFA file(s) with ' ...
                                 'Inconsist Distance(s). \n '...
                                'Keep these SOFA file(s) by removing measurements with '...
                                'Inconsist Distance(s)? (y/n): \n '], 's');
                            % request user input whether they want to keep data by removing 
                            % repeated angles
                            % (after they have read the list of problematic SOFA files)
                        end
                    end       

                     if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                            sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 % double check input
                        if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                        % if user choose yes (y or Y or anything include 'y')
                        % this fuction will add the problematic files to 'good_SOFA' folder

                            warning on
                            warning(['kept SOFA file by removing measurements with '...
                                'Inconsist Distances'])
                            warning off
                            % print warning

                            for n = 1:length(Inconsist_Dist_idx)
                                [~, new_name] = modify_sofa_measurement(...
                                    bad_sofa_struct.File_Name{Inconsist_Dist_idx(n)}, ...
                                    bad_sofa_struct.Inconsist_Dist{n}, 'remove', 'good_SOFA/', []);
                                % modify, rename and copy the problematic SOFA file 
                                % into 'good_SOFA' folder
                                good_hrtf_log = [good_hrtf_log ; ...
                                    {new_name}, 'removed Inconsist Distances'];
                                % update 'good_hrtf_log' by adding the problematic SOFA file names
                                % add noted as 'removed Inconsist Distances' in second column
                            end
                        end
                    else
                        error('input must be y or n.')
                    end
                end

        %% 3.3 for Asymmetry Angle (left vs right without median)

              %%% Thirdly, deal with the SOFA files has Asymmetry Angle (left vs right)
                Asymm_Angles_idx = find(~cellfun(@isempty, bad_sofa_struct.Asymm_Angles));
                % find the files that has Asymmetry Angle
                sofa_Asymm_Angles = bad_sofa_struct.File_Name(Asymm_Angles_idx);
                % get the file names of those files into an array
                [sofa_Asymm_Angles, unique_idx] = unique(sofa_Asymm_Angles);
                % find the unique file names and location inside the file name array
                % just create
                Asymm_Angles_idx = Asymm_Angles_idx(unique_idx);
                % remove duplicated file location in the file name array

                if ~isempty(Asymm_Angles_idx)

                    reenter = []; % initialise user input variable
                    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                         sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                         sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 || isempty(reenter)
                        reenter = input(['\n There is/are SOFA file(s) with Asymmetry '...
                                'Angle(s). \n Keep SOFA file(s) with Asymmetry Angle(s)?'... 
                                '(y/n/list): \n ' '(You can decide whether you want to keep ' ...
                                'the Asymmetry Angle meausrement(s) later on.) \n '], 's');
                    end
                    % request user input whether they want to keep data with 
                    % Asymmetry Angle meausrement(s) or list out all the files that
                    % inculdes Asymmetry Angle meausrement(s).
                    % Use will choose HOW to keep the file later

                    if sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                       % if user types list (l or L or anything include 'l')
                       % this fuction will list out the name of the problematic SOFA file(s)

                        reenter = []; % initialise user input variable
                        while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                              sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || isempty(reenter)

                            for n = 1:length(Asymm_Angles_idx)
                                fprintf(1, '\n %d - %s ', n, ...
                                    bad_sofa_struct.File_Name{Asymm_Angles_idx(n)})
                            end
                            % print out list of problematic file(s)

                            reenter= input(['\n \n Above is/are the SOFA file(s) with ' ...
                                 'Asymmetry Angle meausrement(s). \n '...
                                 'Keep these SOFA file(s) with Asymmetry Angle meausrement(s)? '...
                                 '(y/n): \n ' ...
                                 '(You can decide whether you want to keep ' ...
                                 'the Asymmetry Angle meausrement(s) later on.) \n '], 's');
                            % request user input whether they want to keep data by removing 
                            % Asymmetry Angle meausrement(s
                            % (after they have read the list of problematic SOFA files)
                        end
                    end  

                    if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                            sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 % double check input
                        if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                        % if user choose yes (y or Y or anything include 'y')
                        % which means user wants to keep the data
                        % now user need to decide whether they want to keep the
                        % Asymmetry meausrement(s) or not
                            reenter = []; % initialise user input variable
                            while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                                sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                                sum(reenter == 'D') < 1 && sum(reenter == 'd') < 1 && ...
                                sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 && ...
                                ~isnan(str2double(reenter)) || isempty(reenter)
                                reenter= input(['\n For the SOFA file(s) with ' ...
                                    'Asymmetry Angle meausrement(s). \n '...
                                    'Do you wish to keep the Asymmetry Angle meausrement(s)? \n '...
                                    '(y/n/d(for details)/number(plot asymmetry measurements)'...
                                    '/list): \n \n '...
                                    'NOTE 1: input DETAILS will TRY to show you a reference '...
                                    'of the \n   Minimum, Maximum and Average difference ' ...
                                    'between Asymmetry angles. \n '...
                                    '  This may NOT be accurate depends on the data ' ...
                                    'abnormality. \n '...
                                    'NOTE 2: input ANY NUMBER will plot the asymmetry ' ...
                                    'measurements' ...
                                    'from the list of \n   SOFA file that with asymmetry '...
                                    'measurements. \n '...
                                    '  File number could be check by input LIST. \n '], 's');
                                % Request user input whether they want to keep data by removing 
                                % removing Asymmetry Angle meausrement(s) or keep them.
                                % or
                                % User can also input d or DETAILS to show more details on
                                % the Asymmetry Angle(s), including Minimum, Maximum and
                                % Average in Azimith and Elevation anlge difference. 
                                % or
                                % User can input ANY NUMBER to plot out the Asymmetry Angle
                                % of one of the problematic SOFA file (one at a time)
                                % or
                                % The number of can be problematic SOFA file refercnced 
                                % by input l or LIST, the first column of the list is
                                % the corresponding number.

                                if sum(reenter == 'd') > 0 || sum(reenter == 'D') > 0
                                    % for more details Minimum, Maximum and Average in 
                                    % Azimith and Elevation anlge difference. 

                                    azi_diff_min_max_mean = zeros(length(Asymm_Angles_idx), 3);
                                    ele_diff_min_max_mean = zeros(length(Asymm_Angles_idx), 3);
                                    % inititalise output

                                    % to find out the azimuth and elevation angle
                                    for n = 1 : length(Asymm_Angles_idx) 
                                        % for the number of a abnormal SOFA files
                                        SOFA_angle = SOFAcalculateAPV(SOFAload(...
                                            bad_sofa_struct.File_Name ...
                                            {Asymm_Angles_idx(n)}));
                                        % read all azi and ele angels from the abnormal SOFA file
                                        SOFA_angle_abnormal = [SOFA_angle(...
                                            bad_sofa_struct.Asymm_Angles ...
                                            {Asymm_Angles_idx(n)}, :) ...
                                            bad_sofa_struct.Asymm_Angles{Asymm_Angles_idx(n)}];
                                        % extract only abnormal angles 
                                        % (with the row number added to the last column)
                                        SOFA_angle_abnormal_l = SOFA_angle_abnormal(...
                                            SOFA_angle_abnormal(:,1) >= 0, :);
                                        % extract abnormal angle(s) from the left hand side 
                                        SOFA_angle_abnormal_r = SOFA_angle_abnormal( ...
                                            SOFA_angle_abnormal(:,1) <= 0, :);
                                        % extract abnormal angle(s) from the right hand side
                                        if length(SOFA_angle_abnormal_l) == ...
                                                length(SOFA_angle_abnormal_r)
                                            % check to make sure the number of 
                                            % left abnormal angles is equal to the number of 
                                            % right abnormal angles incase there is outliners.
                                            % if this is not true, the following code won't work, 
                                            % the 'details' of that kind of angel distribution 
                                            % requires huge computational power,
                                            % which should be done by manually 

                                            SOFA_angle_abnormal_l =sortrows(SOFA_angle_abnormal_l);
                                            % sort the left abnormal angle(s) in accending order 
                                            SOFA_angle_abnormal_r =sortrows([abs(...
                                                SOFA_angle_abnormal_r(:,1)) ...
                                                SOFA_angle_abnormal_r(:,2:end)]);
                                            % sort the right abnormal angle(s) in accending order 
                                            % but the take the absolute value of the
                                            % aziuth angle (they should be negative)
                                            SOFA_angle_abnormal_diff = ...
                                                [SOFA_angle_abnormal_l(:,1:3) - ...
                                                SOFA_angle_abnormal_r(:,1:3) ...
                                                SOFA_angle_abnormal_l(:,4) ...
                                                SOFA_angle_abnormal_r(:,4)];
                                            % as both left and right angles got sorted
                                            % out, they should match with each other
                                            % which be the left and right should be the 
                                            % closest angle accroding to row numbers.
                                            % the row number of the left and right angle in the 
                                            % origial SOFA file are added to the end
                                            % just to incase we want to track the location later on
                                            azi_diff_min_max_mean(n,:) = [min(...
                                                SOFA_angle_abnormal_diff(:,1)) ...
                                                max(SOFA_angle_abnormal_diff(:,1)) ...
                                                mean(SOFA_angle_abnormal_diff(:,1))];
                                            ele_diff_min_max_mean(n,:) = [min(...
                                                SOFA_angle_abnormal_diff(:,2)) ...
                                                max(SOFA_angle_abnormal_diff(:,2)) ...
                                                mean(SOFA_angle_abnormal_diff(:,2))];
                                            % calculate the minimum, maximum and average difference
                                            % of the azimuth and elecation angle saperately 
                                        else

                                            warning on
                                            warning(['Unable to calculate details. '... 
                                                ' Asymmetry outliners was found, '...
                                                'number of asymmetry measurements on the left '...
                                                'does not equal to the number of '...
                                                'asymmetry measurements on the right. '])
                                            warning(['Suggest to REMOVE ALL asymmetry angle '... 
                                                'measurements.'])   
                                            warning off
                                            % if left and right does not match, which
                                            % mean there may be an outliner,
                                            % throw a warning
                                        end
                                    end
                                    overall_azi_min_max_mean = [min(azi_diff_min_max_mean(:,1)) ...
                                        max(azi_diff_min_max_mean(:,2)) ...
                                        mean(azi_diff_min_max_mean(:,3))];
                                    overall_ele_min_max_mean = [min(ele_diff_min_max_mean(:,1)) ...
                                        max(ele_diff_min_max_mean(:,2)) ...
                                        mean(ele_diff_min_max_mean(:,3))];
                                    % calculate the overall minimum, maximum and average difference
                                    % of the azimuth and elecation angle in all abnormal SOFA files

                                    fprintf (['\n azimuth difference: \n' ...
                                        ' min: ' num2str(overall_azi_min_max_mean(1)) ...
                                        ' | max: ' num2str(overall_azi_min_max_mean(2)) ...
                                        ' | mean: ' num2str(overall_azi_min_max_mean(3)) ...
                                        '\n elevation difference: \n' ...
                                        ' min: ' num2str(overall_ele_min_max_mean(1)) ...
                                        ' | max: ' num2str(overall_ele_min_max_mean(2)) ...
                                        ' | mean: ' num2str(overall_ele_min_max_mean(3)) '\n'])
                                    % print result

                                    reenter = [];
                                    % reset 'reenter'

                                elseif ~isempty(str2double(regexp(reenter,'\d*','match')'))
                                % if the input is a number then plot angles
                                    plot_idx = str2double(regexp(reenter,'\d*','match')');
                                    % convert the input from string to double

                                    for n = 1: length(plot_idx)

                                        if plot_idx(n) > length(Asymm_Angles_idx) || ...
                                            plot_idx(n) < 1

                                            warning on
                                            warning(['input number ' plot_idx(n) ...
                                                ' excessed number of file limit, '...
                                                'please check the first column of '...
                                                'the list by input "l". '])
                                            warning off
                                        else
                                            figure
                                            plot_3d_angles(bad_sofa_struct.File_Name{...
                                                Asymm_Angles_idx(plot_idx(n))}, ...
                                                bad_sofa_struct.Asymm_Angles{Asymm_Angles_idx(plot_idx(n))})
                                            title(['Asymmetry Angle meausrement(s) in ' ...
                                                bad_sofa_struct.File_Name{Asymm_Angles_idx(plot_idx(n))} ...
                                                '(left and right)'], 'Interpreter', 'none')
                                            % plot the Asymmetry Angle(s) on a new figure
                                        end
                                    end

                                    reenter = [];
                                    % reset 'reenter'
                                end


                                if sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                                % if the input includes l or L, than LIST all abnoramal SOFA files 
                                % with a index number on the first column for PLOT

                                    for n = 1:length(Asymm_Angles_idx)
                                        fprintf(1, '\n %d - %s ', n, ...
                                            bad_sofa_struct.File_Name{Asymm_Angles_idx(n)})
                                    end
                                    % print all abnoramal SOFA files with a index number
                                    % on the first column
                                    fprintf(['\n above is/are the list of SOFA file(s) with '...
                                        'Asymmetry Angle meausrement(s), \n '...
                                        'Input the number in the First Column to plot the '...
                                        'Asymmetry Angle(s) \n '])
                                    % print an extra statment to help user understand

                                    reenter = [];
                                    % reset 'reenter'

                                elseif sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                                % if user choose yes (y or Y or anything include 'y')
                                % to keep all meausrements.
                                % this part will add the problematic files to 
                                % 'good_SOFA' folder and the log
                                    warning on
                                    warning(['kept SOFA file by keeping all measurements even '...
                                        'with Asymmetry Angle(s)'])
                                    warning off
                                    % print warning

                                    for n = 1:length(Asymm_Angles_idx)
                                        copyfile(which(char(...
                                            bad_sofa_struct.File_Name{Asymm_Angles_idx(n)})), ...
                                            'good_SOFA');
                                        % copy the problematic SOFA file into 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            {bad_sofa_struct.File_Name{Asymm_Angles_idx(n)}}, ...
                                            'kept Asymmetry Angles'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA file
                                        %  names, add note 'kept Asymmetry Angles' in second column
                                    end 

                                elseif sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0
                                % if user choose yes (n or N or anything include 'n')
                                % to remove measurements with Asymmetry Angles.
                                % then add to 'good_SOFA' folder and the log.
                                    warning on
                                    warning(['kept SOFA file by removing all measurements '...
                                        'with Asymmetry Angle(s)'])
                                    warning off
                                    % print warning

                                    for n = 1:length(Asymm_Angles_idx)
                                        [~, new_name] = modify_sofa_measurement(...
                                        bad_sofa_struct.File_Name{Asymm_Angles_idx(n)}, ...
                                        bad_sofa_struct.Asymm_Angles{n}, 'remove', 'good_SOFA/',[]);
                                        % modify, rename and copy the problematic SOFA file into 
                                        % 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            {new_name}, 'removed Asymmetry Angles'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA file
                                        % names, add note 'removed Asymmetry Angles' in second column
                                    end 
                                end   
                            end
                        end
                    else
                        error(['input must be Y(es) or N(o) or D(etails) or L(ist) or ' ...
                            'a single number.'])
                    end
                end


        %% 3.4 for Asymmetry Median ( front (0 degree) vs back (180 degree))

              %%% Fourthly, deal with the SOFA files has Asymmetry Median Angle 
                %(front(0) vs back(180))
                Asymm_Median_idx = find(~cellfun(@isempty, bad_sofa_struct.Asymm_Median));
                % find the files that has Asymmetry Median
                sofa_Asymm_Median = bad_sofa_struct.File_Name(Asymm_Median_idx);
                % get the file names of those files into an array
                [sofa_Asymm_Median, unique_idx] = unique(sofa_Asymm_Median);
                % find the unique file names and location inside the file name array
                % just create
                Asymm_Median_idx = Asymm_Median_idx(unique_idx);
                % remove duplicated file location in the file name array

                if ~isempty(Asymm_Median_idx)

                    reenter = []; % initialise user input variable
                    while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                         sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                         sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 || isempty(reenter)
                        reenter = input(['\n There is/are SOFA file(s) with '...
                            'Asymmetry Median(s). (front(0) vs back(180)) \n '...
                            'Keep SOFA file(s) with Asymmetry Angle(s)?'... 
                            '(y/n/list): \n ' '(You can decide whether you want to keep ' ...
                            'the Asymmetry Median meausrement(s) later on.) \n '], 's');
                    end
                    % request user input whether they want to keep data with 
                    % Asymmetry Median meausrement(s) or list out all the files that
                    % inculdes Asymmetry Median meausrement(s).
                    % Use could choose HOW to keep the file later

                    if sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                       % if user types list (l or L or anything include 'l')
                       % this fuction will list out the name of the problematic SOFA file(s)

                        reenter = []; % initialise user input variable
                        while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                              sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 || ...
                              isempty(reenter)

                            for n = 1:length(Asymm_Median_idx)
                                fprintf(1, '\n %d - %s ', n, ...
                                    bad_sofa_struct.File_Name{Asymm_Median_idx(n)})
                            end
                            % print out list of problematic file(s)

                            reenter= input(['\n \n Above is/are the SOFA file(s) with ' ...
                                 'Asymmetry Median meausrement(s). \n '...
                                 'Keep these SOFA file(s) with Asymmetry Medain '...
                                 'meausrement(s)? (y/n): \n ' ...
                                 '(You can decide whether you want to keep ' ...
                                 'the Asymmetry Median meausrement(s) later on.) \n '], 's');
                            % request user input whether they want to keep data by removing 
                            % Asymmetry Median meausrement(s)
                            % (after they have read the list of problematic SOFA files)
                        end
                    end  

                    if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0 || ...
                            sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0 
                        % double check input
                        if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                        % if user choose yes (y or Y or anything include 'y')
                        % which means user wants to keep the data
                        % now user need to decide whether they want to keep the
                        % Asymmetry Median meausrement(s) or not
                            reenter = []; % initialise user input variable
                            while sum(reenter == 'Y') < 1 && sum(reenter == 'y') < 1 && ...
                                sum(reenter == 'N') < 1 && sum(reenter == 'n') < 1 && ...
                                sum(reenter == 'D') < 1 && sum(reenter == 'd') < 1 && ...
                                sum(reenter == 'L') < 1 && sum(reenter == 'l') < 1 && ...
                                ~isnan(str2double(reenter)) || isempty(reenter)
                                reenter= input(['\n For the SOFA file(s) with ' ...
                                    'Asymmetry Median meausrement(s). \n '...
                                    'Do you wish to keep the Asymmetry Median meausrement(s)? \n '...
                                    '(y/n/d(for details)/number(plot asymmetry measurements)'...
                                    '/list): \n \n '...
                                    'NOTE 1: input DETAILS will TRY to show you a reference '...
                                    'of the \n Minimum, Maximum and Average difference' ...
                                    ' between Asymmetry angles. \n '...
                                    'This may NOT be accurate depends on the data '...
                                    'abnormality. \n '...
                                    'NOTE 2: input ANY NUMBER will plot the asymmetry ' ...
                                    'measurements from the list of \n SOFA file that '...
                                    'with asymmetry measurements. \n '...
                                    'File number could be check by input LIST. \n '], 's');
                                % Request user input whether they want to keep data by removing 
                                % removing Asymmetry Median meausrement(s) or keep them.
                                % or
                                % User can also input d or DETAILS to show more details on
                                % the Asymmetry Median, including Minimum, Maximum and
                                % Average in Azimith and Elevation anlge difference. 
                                % or
                                % User can input ANY NUMBER to plot out the Asymmetry Angle
                                % of one of the problematic SOFA file (one at a time)
                                % or
                                % The number of can be problematic SOFA file refercnced 
                                % by input l or LIST, the first column of the list is
                                % the corresponding number.


                                if sum(reenter == 'd') > 0 || sum(reenter == 'D') > 0
                                    % for more details Minimum, Maximum and Average in 
                                    % Azimith and Elevation anlge difference. 

                                    azi_med_diff_min_max_mean = zeros(length(Asymm_Median_idx), 3);
                                    ele_med_diff_min_max_mean = zeros(length(Asymm_Median_idx), 3);
                                    % inititalise output

                                    % to find out the azimuth and elevation angle
                                    for n = 1 : length(Asymm_Median_idx) 
                                        % for the number of a abnormal SOFA files
                                        SOFA_med_angle = SOFAcalculateAPV(SOFAload( ...
                                            bad_sofa_struct.File_Name {Asymm_Median_idx(n)}));
                                        % read all azi and ele angels from the abnormal SOFA file
                                        SOFA_med_angle_abnormal = [SOFA_med_angle( ...
                                            bad_sofa_struct.Asymm_Median {Asymm_Median_idx(n)}, :)...
                                            bad_sofa_struct.Asymm_Median{Asymm_Median_idx(n)}];
                                        % extract only abnormal angles 
                                        % (with the row number added to the last column)
                                        SOFA_med_angle_abnormal_f = SOFA_med_angle_abnormal(...
                                            abs(SOFA_med_angle_abnormal(:,1)) <= 90, :);
                                        % extract abnormal angle(s) at the front
                                        SOFA_med_angle_abnormal_b = SOFA_med_angle_abnormal( ...
                                            abs(SOFA_med_angle_abnormal(:,1)) >= 90, :);
                                        % extract abnormal angle(s) at the back
                                        if length(SOFA_med_angle_abnormal_f) == ...
                                                length(SOFA_med_angle_abnormal_b)
                                            % check to make sure the number of left abnormal angles 
                                            % is equal to the number of right abnormal angles incase  
                                            % there is outliners.
                                            % if this is not true, the following code won't work, 
                                            % the 'details' of that kind of angel distribution  
                                            % requires huge computational power, 
                                            % which should be done by manually 

                                            SOFA_med_angle_abnormal_f = ...
                                                sortrows(SOFA_med_angle_abnormal_f);
                                            % sort the front abnormal angle(s) in accending order 
                                            SOFA_med_angle_abnormal_b = ...
                                                sortrows(SOFA_med_angle_abnormal_b);
                                            % sort the back abnormal angle(s) in accending order 
                                            SOFA_med_angle_abnormal_diff = ...
                                                [SOFA_med_angle_abnormal_f(:,1:3) ...
                                                - SOFA_med_angle_abnormal_b(:,1:3) ...
                                                SOFA_med_angle_abnormal_f(:,4) ...
                                                SOFA_med_angle_abnormal_b(:,4)];
                                            % as both front and back angles got sorted
                                            % out, they should match with each other
                                            % which meand the front and back should be the 
                                            % closest angle accroding to row numbers.
                                            % the row number of the front and back angle in the 
                                            % origial SOFA file are added to the end
                                            % just to incase we want to track the location later on
                                            azi_med_diff_min_max_mean(n,:) = [min(...
                                                SOFA_med_angle_abnormal_diff(:,1)) ...
                                                max(SOFA_med_angle_abnormal_diff(:,1)) ...
                                                mean(SOFA_med_angle_abnormal_diff(:,1))];
                                            ele_med_diff_min_max_mean(n,:) = [min(...
                                                SOFA_med_angle_abnormal_diff(:,2)) ...
                                                max(SOFA_med_angle_abnormal_diff(:,2)) ...
                                                mean(SOFA_med_angle_abnormal_diff(:,2))];
                                            % calculate the minimum, maximum and average difference
                                            % of the azimuth and elecation angle saperately 
                                        else

                                            warning on
                                            warning(['Unable to calculate details. '... 
                                                'Asymmetry outliners was found, '...
                                                'number of asymmetry measurements at the front '...
                                                'does not equal to the number of asymmetry '...
                                                'measurements at the right. '])
                                            warning('Suggest to REMOVE ALL asymmetry measurements.')   
                                            warning off
                                            % if number of front and back median meaurements 
                                            % does not match, which mean there may be an outliner,
                                            % throw a warning
                                        end
                                    end
                                    overall_med_azi_min_max_mean = ...
                                        [min(azi_med_diff_min_max_mean(:,1)) ...
                                        max(azi_med_diff_min_max_mean(:,2)) ...
                                        mean(azi_med_diff_min_max_mean(:,3))];
                                    overall_med_med_ele_min_max_mean = ...
                                        [min(ele_med_diff_min_max_mean(:,1)) ...
                                        max(ele_med_diff_min_max_mean(:,2)) ...
                                        mean(ele_med_diff_min_max_mean(:,3))];
                                    % calculate the overall minimum, maximum and average difference
                                    % of the azimuth and elecation angle in all abnormal SOFA files

                                    fprintf (['\n azimuth difference in Median: \n' ...
                                        ' min: ' num2str(overall_med_azi_min_max_mean(1)) ...
                                        ' | max: ' num2str(overall_med_azi_min_max_mean(2)) ...
                                        ' | mean: ' num2str(overall_med_azi_min_max_mean(3)) ...
                                        '\n elevation difference in Median: \n' ...
                                        ' min: ' num2str(overall_med_ele_min_max_mean(1)) ...
                                        ' | max: ' num2str(overall_med_ele_min_max_mean(2)) ...
                                        ' | mean: ' num2str(overall_med_ele_min_max_mean(3)) '\n'])
                                    % print result

                                    reenter = [];
                                    % reset 'reenter'

                                elseif ~isempty(str2double(regexp(reenter,'\d*','match')'))
                                % if the input is a number then plot angles
                                    plot_idx = str2double(regexp(reenter,'\d*','match')');
                                    % convert the input from string to double

                                    for n = 1: length(plot_idx)

                                        if plot_idx(n) > length(Asymm_Median_idx) || ...
                                            plot_idx(n) < 1

                                            warning on
                                            warning(['input number ' plot_idx(n) ...
                                                ' excessed number of file limit, '...
                                                'please check the first column of '...
                                                'the list by input "l". '])
                                            warning off
                                        else
                                            figure
                                            plot_3d_angles(bad_sofa_struct.File_Name{...
                                                Asymm_Median_idx(plot_idx(n))}, ...
                                                bad_sofa_struct.Asymm_Median{Asymm_Median_idx(plot_idx(n))})
                                            title(['Asymmetry Median meausrement(s) in ' ...
                                                bad_sofa_struct.File_Name{Asymm_Median_idx(plot_idx(n))} ...
                                                '(front and back)'], 'Interpreter', 'none')
                                            % plot the Asymmetry Angle(s) on a new figure
                                        end
                                    end

                                    reenter = [];
                                    % reset 'reenter'

                                elseif sum(reenter == 'L') > 0 || sum(reenter == 'l') > 0
                                % if the input includes l or L, than LIST all abnoramal SOFA files 
                                % with a index number on the first column for PLOT

                                    for n = 1:length(Asymm_Median_idx)
                                        fprintf(1, '\n %d - %s ', n, ...
                                            bad_sofa_struct.File_Name{Asymm_Median_idx(n)})
                                    end
                                    % print all abnoramal SOFA files with a index number
                                    % on the first column
                                    fprintf(['\n above is/are the list of SOFA file(s) with '...
                                        'Asymmetry Median meausrement(s), \n '...
                                        'Input the number in the First Column to plot the '...
                                        'Asymmetry Angle(s) \n '])
                                    % print an extra statment to help user understand

                                    reenter = [];
                                    % reset 'reenter'
                                end

                                if sum(reenter == 'Y') > 0 || sum(reenter == 'y') > 0
                                % if user choose yes (y or Y or anything include 'y')
                                % to keep all meausrements.
                                % this part will add the problematic files to 
                                % 'good_SOFA' folder and the log
                                    warning on
                                    warning(['kept SOFA file by keeping all measurement(s) even '...
                                        'with Asymmetry Medina angle'])
                                    warning off
                                    % print warning

                                    for n = 1:length(Asymm_Median_idx)
                                        copyfile(which(char(...
                                            bad_sofa_struct.File_Name{Asymm_Median_idx(n)})), ...
                                            'good_SOFA');
                                        % copy the problematic SOFA file into 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            {bad_sofa_struct.File_Name{Asymm_Median_idx(n)}}, ...
                                            'kept Asymmetry Median'];
                                        % update 'good_hrtf_log' by adding the problematic
                                        % SOFA file names, with a noted 
                                        % 'kept Asymmetry Median' in second column
                                    end 

                                elseif sum(reenter == 'N') > 0 || sum(reenter == 'n') > 0
                                % if user choose yes (n or N or anything include 'n')
                                % to remove measurements with Asymmetry Angles.
                                % then add to 'good_SOFA' folder and the log.
                                    warning on
                                    warning(['kept SOFA file by removing all measurements '...
                                        'with Asymmetry Median(s)'])
                                    warning off
                                    % print warning

                                    for n = 1:length(Asymm_Median_idx)
                                        [~, new_name] = modify_sofa_measurement(...
                                        bad_sofa_struct.File_Name{Asymm_Median_idx(n)}, ...
                                        bad_sofa_struct.Asymm_Median{n}, 'remove', 'good_SOFA/',[]);
                                        % modify, rename and copy the problematic SOFA file into 
                                        % 'good_SOFA' folder
                                        good_hrtf_log = [good_hrtf_log ; ...
                                            {new_name}, 'removed Asymmetry Median'];
                                        % update 'good_hrtf_log' by adding the problematic SOFA file 
                                        % names add noted 'removed Asymmetry Median' in second column
                                    end 
                                end   
                            end
                        end
                    else
                        error('input must be Y(es) or N(o) or D(etails) or L(ist) or single number.')
                    end
                end
            end

    %%  4.0 Wrap up the function by adding 'good_SOFA' folder to path
     %      and remind user to rename the folder to avoid over-write

            addpath('good_SOFA')
            % add folder to search path
            save('good_SOFA/good_hrtf_log.mat','good_hrtf_log')
            % save log
            message = sprintf(['Kept both good and maybe some abnormal SOFA file(s). \n' ...
                '         Check the good_sofa_log for further infomation. \n         ' ...
                '(column 1 : file name, column 2: note (good or details of the problem)) \n'...
                '         Please remember to RENAME "good_SOFA" folder or it will be ' ...
                'over-written in the next run.' ]);
            warning on
            warning off backtrace
            sldiagviewer.reportWarning(message);
            warning on backtrace
            % print final reminder

        else

            % if user decided to keep all the file
            addpath('good_SOFA')
            % add folder to search path
            save('good_SOFA/good_hrtf_log.mat','good_hrtf_log')
            % save log
            message = sprintf(['Kept only the perfectly fine SOFA file(s) \n' ...
                '         Please remember to RENAME "good_SOFA" folder or it will be ' ...
                'over-written in the next run.' ]);
            warning on
            warning off backtrace
            sldiagviewer.reportWarning(message);
            warning on backtrace
            % print final reminder

        end
    end
end
addpath('good_SOFA')
warning on

end

