function [ final_result ,varargout ] = common_angle( hrtf_array, range, plot, method )
%COMMON_ANGLE Summary of this function goes here
%
% Find common HRTF measurement angles between two HRTF dataset
% Note that the output HRTF angle is FROM A LISTENER POINT OF VIEW
%
% e.g [final_result, HRTF_1, HRTF_2, HRTF_3] = common_angle( {'subject_015.sofa';...
%       'MRT07.sofa';'irc_1007.sofa'}, 2, 0.001);
%

% INPUT:
% e.g. common_angle({'subject_015.sofa'; 'MRT07.sofa';'irc_1007.sofa'}, 2, 0.1);
%
% 1. hrtf_array = {'subject_015.sofa'; 'MRT07.sofa';'irc_1007.sofa'}
%    - cell array of different hrtf file name 
%
% 2. range = 2 % matching range (in degree)
%
% 3. plot = 1 % plot trigger
%    - 1 || 'true' || 'yes' = plot, 
%    - 0.0001 - 0.9999 = plot with an offset between each hrtfs matching result (in meter)
%      (suguest between 0.01 - 0.05)
%    - 'ture distance' || 'real' || 'distance' = real distance
%    - 0 = no plot % or else
%
%    - 'ture' or 'true distance', plot the real distance of the hrtf measurements
%    - if 'plot' is in between 0 - 1, plot with the offset between each hrtfs matching 
%      result which euqal to the 'plot' input number
%    - otherwise the plot will show both hrtf measurement in 1.2 meters for better comparison 
%
% 4. method = 'max' % matching method, 'max', 'min','distance max', 'distance min'  or 'fast'
%    - (max: as many as pssoble, min: as close as possible)
%      the option only avalible between 2 HRTF datasets
%      because using 'max' with HRTF datasets more than 2 will end up with
%      exponential big results
%    - defult method is 'min'
%
%    - max: all matched data possible (will duplicate the data if neccessary)
%    - min: only keep the closest match unless there are more than one in same distance
%    - distance max: similar to 'max' but use great circle distance (more accurate sometimes)
%    - distance min: similar to 'min' but use great circle distance (more
%      accurate sometimes)
%    - fast: fast estimate matching result similar to 'max'
%      (works better with small range, e.g. < 5, and there is no distance details)

%    - the normal mehtod should work the same as 'distance' method but faster. 
%      it is a hybrid method combining fast estimate and great circle distance
%

% OUTPUT:
% 1. final_result; % matched angles in cell arrays
%    - final_result{1} represent hrtf_array{1}... final_result{n} represent hrtf_array{n}
%      cloumn 1: azimuth angle 
%      cloumn 2: elevation angle
%      cloumn 3: distance
%      cloumn 4: data row number in sofa file)
%
% 2. varargout; % variable size ouput
%    - set different output variable e.g. [HRTF_1, HRTF_2, HRTF_3] they will be
%    - equivalent to final_result{1}, final_result{2} and final_result{3} respectively 
%    - make sure the varargout is less than the length of hrtf_array
%
%

%% initialise input and output

if nargin == 1
     range = 0;
     plot = 0;
     method = 'min';
elseif nargin == 2
     plot = 0;
     method = 'min';
elseif nargin == 3
     method = 'min';
end
% catch empty inputs

if isstring(hrtf_array)
    hrtf_array = cellstr(hrtf_array)';
end
% catch if hrtf_array is string array instead of cell array

hrtf_array = reshape(hrtf_array,[],1);
% reshape hrtf_array into single column array

if strncmpi(plot,'true dist', 6) || strncmpi(plot,'real', 4) || strncmpi(plot, 'dist', 4) 
    plot_method = 'true distance';
    plot = 1;
elseif strncmpi(plot, 'plot', 4) || strncmpi(plot,'true', 4) || strncmpi(plot,'yes', 1)
    plot_method = 'assigned distance';
    plot_offset = 0;
    plot = 1;
elseif plot >= 1
    plot_method = 'assigned distance';
    plot_offset = 0;
    plot = 1;
elseif plot < 1 && plot >0
    plot_method = 'assigned distance';
    plot_offset = plot;
    plot = 1;  
else
    plot = 0;
end
% decide plot method

nout = max(nargout,1) - 1;
if nout > length(hrtf_array)
    error('Too many output variables (should be equal or less than inputs)')
end
% catch error when number of output variables are more than number of input

%% Compare hrtf code start here

if length(hrtf_array) < 2
    error('not enough hrtf for comparison')
elseif length(hrtf_array) == 2
    
%% compare between 2 hrtfs

    [hrtf_angle_1_in_2, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(hrtf_array{1}, ...
        hrtf_array{2}, range, 0, method);
    final_result{1} = hrtf_angle_1_in_2;
	final_result{2} = hrtf_angle_2_in_1;  

else
%% compare between more than 2 hrtfs

method = 'min';

data_length = length(hrtf_array);
% compare_count = data_length * (data_length - 1) /2;
results_between_2 = cell(data_length, data_length - 1);
remain = 2;
%figure
for n = 1: data_length - 1
    for m = remain : data_length
        [hrtf_angle_1_in_2, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(hrtf_array{n}, ...
            hrtf_array{m}, range, 0, method);
        %hold on
        results_between_2{n,m-1} = hrtf_angle_1_in_2;
        results_between_2{m,n} = hrtf_angle_2_in_1; 
        
        if isempty(results_between_2{n,m-1})
            error('common_angle:no2Match',['could not find any match between %s & %s, '...
                'please consider remove one of them or set a wider range.'], ...
                hrtf_array{n}, hrtf_array{m})
        end
            
    end
    remain = remain + 1;
end
% find match in all combination between 2 hrtf

common_angles1 = cell(data_length, 1);
for n = 1 : size(results_between_2,1)
    common_angles1{n} = intersect(results_between_2{n, 1}(:,1:3), ...
        results_between_2{n, 2}, 'rows');
    common_angles1{n} = results_between_2{n, 1}(ismember(results_between_2{n, 1}, ...
        common_angles1{n}(:,1:3),'rows'),(1:3));
    for m = 3:size(results_between_2,2)
        common_angles1{n} = intersect(common_angles1{n}, results_between_2{n, m}, 'rows');
        common_angles1{n} = results_between_2{n, 1}(ismember(results_between_2{n, 1}, ...
            common_angles1{n},'rows'),(1:3));
    end
end
% keep the common ones in hrtf dataset bases


data_length2 = length(common_angles1);
%compare_count2 = data_length2 * (data_length2 - 1) /2;
compare_result = cell(data_length2, data_length2 - 1);
remain = 2;
%figure
for n = 1: data_length2 - 1
    for m = remain : data_length2
        [hrtf_angle_1_in_2, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(common_angles1{n}, ...
           common_angles1{m}, range, 0, method);
        %hold on
        compare_result{n,m-1} = hrtf_angle_1_in_2;
        compare_result{m,n} = hrtf_angle_2_in_1;       
    end
    remain = remain + 1;
end

common_angles2 = cell(data_length, 1);
for n = 1 : size(compare_result,1)
    common_angles2{n} =  intersect(compare_result{n, 1}, compare_result{n, 2}, 'rows');
    for m = 3:size(compare_result,2)
        common_angles2{n} =  intersect(common_angles2{n}, compare_result{n, m}, 'rows');
    end
end
% remove left out data by comparing the remain HRTFs again


data_length3 = length(common_angles2);
%compare_count3 = data_length3 * (data_length3 - 1) /2;
compare_result3 = cell(data_length3, data_length3 - 1);
remain = 2;
%figure
for n = 1: data_length3 - 1
    for m = remain : data_length3
        [hrtf_angle_1_in_2, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(common_angles2{n}, ...
            common_angles2{m}, range, 0, method);
        %hold on
        compare_result3{n,m-1} = hrtf_angle_1_in_2;
        compare_result3{m,n} = hrtf_angle_2_in_1;       
    end
    remain = remain + 1;
end

common_angles3 = cell(data_length3, 1);
for n = 1 : size(compare_result3,1)
    common_angles3{n} =  common_angles_2_HRTFs(compare_result3{n, 1}, ...
        compare_result3{n, 2}, range, 0, 'max');
    
    for m = 3:size(compare_result3,2)
        common_angles3{n} =  common_angles_2_HRTFs(common_angles3{n}, ...
            compare_result3{n, m}, range, 0, 'max');
    end
end
% remove left out data by comparing the remain HRTFs again


data_length4 = length(common_angles3);
%compare_count4 = data_length4 * (data_length4 - 1) /2;
final_result = cell(data_length4, 1);
%figure
for n = 1: data_length4 - 1
    [hrtf_angle_1_in_2, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(common_angles3{1}, ...
        common_angles3{2}, range, 0, method);
    final_result{1} = hrtf_angle_1_in_2;
    final_result{2} = hrtf_angle_2_in_1;
    for m = 3:length(final_result)
        [~, hrtf_angle_2_in_1]  = common_angles_2_HRTFs(common_angles3{1}, ...
            common_angles3{m}, range, 0, method);
        final_result{m} = hrtf_angle_2_in_1;
    end
end
% re-organise data to make sure the angles matches between datasets

end

%% if there is no match

empty_result = zeros(length(final_result), 1);
for n = 1:length(final_result)
	empty_result(n) = isempty(final_result{n});
end
% find empty hrtf in final_result

if sum(empty_result) > 0
    empty_result = logical(empty_result);
    empty_hrtf_array = strjoin({hrtf_array{empty_result}}, ', '); 
    % empty hrtf string array 

    error('common_angle:noMatch', ['could not find any match in ' empty_hrtf_array ...
        ', please consider matching different datasets or set a wider range'])
    % print and kill function
end

%% find data numbers

for n = 1:length(hrtf_array)
    [~, ~, row] = fetchHRTFbyAngle(hrtf_array{n}, final_result{n}(:,1), ...
        final_result{n}(:,2), 0, 0);
    final_result{n} = [final_result{n} row];
end
% find the matching data number

%% create variables

nout = max(nargout,1) - 1;
for k = 1:nout
	varargout{k} = final_result{k};
end
% create output variables

%% plot result
if plot ~= 0 
    
    rng(16);
	markers = 'd+o*sxv^p<>h';
	z = 1.2;
    
	figure;
    if strcmp('true distance', plot_method) 
        for  n = 1:length(final_result)
            plot_3d_angles(final_result{n}(:,1), final_result{n}(:,2), final_result{n}(:,3), ...
            markers(mod(n, length(markers))), 'MarkerEdgeColor', rand(1,3), ...
            'DisplayName', [char(hrtf_array{n}) ', d = ' num2str(final_result{n}(1,3))])
            hold on
        end
    else
        for n = 1:length(final_result)
            plot_3d_angles(final_result{n}(:,1), final_result{n}(:,2), z, ...
            markers(mod(n, length(markers))), 'MarkerEdgeColor', rand(1,3), ...
            'DisplayName', [char(hrtf_array{n}) ', d = ' num2str(z)])
            hold on
            z = z + plot_offset;
        end
    end
    % plot common angle on one graph in real distance or assigned distance
    
    if strcmp('fast', method) 
        title(['fast estimate common HRTF angles - ' plot_method])
    elseif strcmp('distance max', method) 
        title(['common HRTF angles by great circle distance (max) - ' plot_method])
    elseif strcmp('distance min', method) || strcmp('distance', method)
        title(['common HRTF angles by great circle distance (min) - ' plot_method])
    elseif strcmp('max', method) 
        title(['common HRTF angles (max) - ' plot_method])
    else 
        title(['common HRTF angles - ' plot_method])
    end
    % set different title by method
    
    legend('Location','southoutside')
    legend('boxoff')
    lgd = legend;
    lgd.FontSize = 12;
    lgd.Interpreter = 'none';
    set(gcf,'Position', [300, 130, 600, 680])
    % configure legend and plot window size
    
end

end

