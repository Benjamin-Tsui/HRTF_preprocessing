function [HRTF_1_common_angles,HRTF_2_common_angles, HRTF_data_desire, HRTF_data_opposite] = ...
    common_angles_2_HRTFs(HRTF_1_or_angle, HRTF_2_or_angle, range, plot, method, print)
%COMMON_ANGLES_IN_2_HRTFS Summary of this function goes here
%
% Find common HRTF measurement angles between two HRTF dataset
% Note that the output HRTF angle is FROM A LISTENER POINT OF VIEW
%
% 
% INPUT:
% e.g. common_angles_2_HRTFs('MRT08.sofa', 'irc_1030.sofa', 0.5, 1, 'min');
%
% HRTF_1 = 'hrtf b_nh2.sofa';  % ARI database
% HRTF_2_or_angle = 'subject_009.sofa';  % CIPIC database
% or
% HRTF_2_or_angle = [azi ele]; % angle matrix
%
% range = 2 % matching range
% plot = 1 % plot trigger, 1 = plot, 0 = no plot, 'ture' = real distance
% 'ture' or 'true distance' plot the real distance of the hrtf measurements
% otherwise the plot will show both hrtf measurement in 1.5 meters
% for better comparison 
%
% method = 'max' % matching method, 'max', 'min','distance max', 'distance min'  or 'fast'
% (max:as many as pssoble, min: as close as possible)
% max: all matched data possible (will duplicate the data if neccessary)
% min: only keep the closest match unless there are more than one in same distance
% distance max: similar to 'max' but use great circle distance (more
% accurate sometimes)
% distance min: similar to 'min' but use great circle distance (more
% accurate sometimes)
% fast: fast estimate matching result similar to 'max'
% (works better with small range, e.g. < 5, and there is no distance details) 
%
% the normal mehtod should work the same as 'distance' method but faster. 
% it is a hybrid method combining fast estimate and great circle distance
%
% defult method is 'min'
%
% print = 0 % if 0: no progress will be shown in command window (advance purpose)
%
%
% OUTPUT:
% HRTF_1_common_angles  % matched angles in hrtf 1 
% HRTF_2_common_angles  % matched angles in hrtf 2 
% (cloumn 1: azimuth angle, 2: elevation angle, 3: distance (not useful))
%
% HRTF_data_desire   % matches between hrtf 1 and hrtf 2 in chosen method(e.g. min)
% HRTF_data_opposite   % matches between hrtf 1 and hrtf 2 in opposite method(e.g. max)
% (cloumn 1: hrtf 1 data number, 2: hrtf 2 data number, 3: distance between two point)
% 
% 

if nargin == 2
     range = 0;
     plot = 0;
     method = 'min';
     print = 1;
elseif nargin == 3
     plot = 0;
     method = 'min';
     print = 1;
elseif nargin == 4
     method = 'min';
     print = 1;
elseif nargin == 5
     print = 1;
end
% Catch empty inputs


if strcmp('true distance', plot) || strcmp('true', plot) || strcmp('true dist', plot) 
    plot_method = 'true distance';
    plot = 1;
elseif strcmp('plot', plot)
    plot_method = 'fixed distance';
    plot = 1;
else
    plot_method = 'fixed distance';
end


% determine HRTF_1_or_angle is sofa file or angles
if ischar(HRTF_1_or_angle)
    HRTF_1 = HRTF_1_or_angle;
    hrtf_1 = SOFAload(HRTF_1); %load hrtf 2 dataset
    hrtf_1_raw_angle = SOFAcalculateAPV(hrtf_1); % load hrtf 2 source position
else
    hrtf_1_raw_angle = HRTF_1_or_angle;
end


% determine HRTF_2_or_angle is sofa file or angles
if ischar(HRTF_2_or_angle)
    HRTF_2 = HRTF_2_or_angle;
    hrtf_2 = SOFAload(HRTF_2); %load hrtf 2 dataset
    hrtf_2_raw_angle = SOFAcalculateAPV(hrtf_2); % load hrtf 2 source position
else
    hrtf_2_raw_angle = HRTF_2_or_angle;
end

[~,idx_1] = sortrows(round(hrtf_1_raw_angle));
hrtf_1_angle = hrtf_1_raw_angle(idx_1,:);    
[~,idx_2] = sortrows(round(hrtf_2_raw_angle));
hrtf_2_angle = hrtf_2_raw_angle(idx_2,:); 
% sort row by column 1

%% find common angles row number between hrtf1 and hrtf2 within a certain range

if print ~= 0
    tic
end

hrtf1_com_2 = cell(size(hrtf_1_angle,1), 1);
hrtf2_com_1 = cell(size(hrtf_2_angle,1), 1);
matched = [0 0 0];
reverseStr = '';
i = 0;
j = 0;
if range < 180
    est_range = atand((tand(range/2)*90)/(2*range)) *2 + range;
    est_range_remain = 90-(2*range);
    % set range to estimate as the great circle distance is not equivalent 
    % to azimuth angle (made it a bit larger to be safe)
    % est_range_remain covers the remain data on the top and bottom
else
    est_range = 360;
    est_range_remain = 0;
    warning('range larger than 180 degrees')
end
for n = 1: size(hrtf_1_angle, 1)
    for m = 1: size(hrtf_2_angle, 1)
        
        if strcmp('fast',method) 
            if (abs(mod((mod(hrtf_1_angle(n,1),360) - mod(hrtf_2_angle(m,1), 360)), 360)) <= range || ...
                    abs(mod((mod(hrtf_2_angle(m,1),360) - mod(hrtf_1_angle(n,1), 360)), 360)) <= range) ...
                    && (abs(hrtf_1_angle(n,2) - hrtf_2_angle(m,2)) <= range || ...
                    abs(hrtf_2_angle(m,2) - hrtf_1_angle(n,2)) <= range) %%%%
                hrtf1_com_2{n} = [hrtf1_com_2{n} m];
                hrtf2_com_1{m} = [hrtf2_com_1{m} n]; 
                matched = [matched; n m 0];
            end 
            % 'fast' method for fast estimate result, by calculting the
            % angle difference in azithum and elevation
            % (only avalive in 'max', 'max' and 'min' are the same result)  
            
        elseif strcmp('distance min',method) || strcmp('distance max',method) || strcmp('distance',method)
            arclen = distance('gc', mod(hrtf_1_angle(n,2),360), hrtf_1_angle(n,1),...
            mod(hrtf_2_angle(m,2), 360), hrtf_2_angle(m,1), 'degrees');
            % calculate the distance between two measurement point
            if  arclen<= range 
                hrtf1_com_2{n} = [hrtf1_com_2{n} m]; % matched hrtfs in cell array 
                % (row number = data number)
                hrtf2_com_1{m} = [hrtf2_com_1{m} n]; % matched hrtfs in cell array 
                % (row number = data number)
                matched = [matched; n m arclen]; % all combination and distance between them
            end
            % 'distance' method for the most accuate result by using great 
            % circle distance
    
        else 
            
            if (abs(mod((mod(hrtf_1_angle(n,1),360) - mod(hrtf_2_angle(m,1), 360)), 360)) <= est_range || ...
                    abs(mod((mod(hrtf_2_angle(m,1),360) - mod(hrtf_1_angle(n,1), 360)), 360)) <= est_range)...
                    && (abs(hrtf_1_angle(n,2) - hrtf_2_angle(m,2)) <= range || ...
                    abs(hrtf_2_angle(m,2) - hrtf_1_angle(n,2)) <= range) || ...
                    abs(hrtf_2_angle(m,2)) >= est_range_remain || abs(hrtf_1_angle(n,2)) >= est_range_remain
                    %abs(hrtf_2_angle(m,2)) >= 80 || abs(hrtf_1_angle(n,2)) >= 80
                
                
                arclen = distance('gc', mod(hrtf_1_angle(n,2),360), hrtf_1_angle(n,1),...
                    mod(hrtf_2_angle(m,2), 360), hrtf_2_angle(m,1), 'degrees');

                if  arclen<= range 
                    hrtf1_com_2{n} = [hrtf1_com_2{n} m]; % matched hrtfs in cell array 
                    % (row number = data number)
                    hrtf2_com_1{m} = [hrtf2_com_1{m} n]; % matched hrtfs in cell array 
                    % (row number = data number)
                    matched = [matched; n m arclen]; % all combination and distance between them
                end
            end
            % combine the 'fast' and 'distance' method for fast and
            % accurate result
            % use 'fast' method to pick out possible match then calculte
            % their distance to confirm result
            %
            % note: the reason use 'range^2 + range' is because on a sphere
            % the azithum angle will become narrower when the elevation
            % increase. e.g. 45-55, 0 have similar distnace to 20-80, 70.
            % And all angels above ele 80 or below ele -80 will use the
            % 'distance' method
            
        end

        
        i = i + 1;
        percent = 100 * i / (size(hrtf_1_angle, 1) * size(hrtf_2_angle, 1)); 
        if floor(percent) == j
            %fprintf('Finding match: %d/10 is done \n', j)
            %j = j + 1;
            if print ~= 0
                msg = sprintf('Finding match: %d/100 is done (in %d x %d = %d data)', ...
                    j, size(hrtf_1_angle, 1), size(hrtf_2_angle, 1), ...
                    (size(hrtf_1_angle, 1) * size(hrtf_2_angle, 1)));
                fprintf([reverseStr, msg, '\n']);
                reverseStr = repmat(sprintf('\b'), 1, length(msg) + 1);
                % calculate and print progress in percentage
            end
            j = j + 1;
        end
        
        %{
        msg = sprintf('Processed percentage: %.2f', percent);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        % calculate and print progress in percentage
        %}
    end
end

if print ~= 0
    toc
end

matched = matched(2: end,:);

%% find all matched data possible (max)
 
hrtf_max_match = matched;
% maximum number matched data (will duplicate all angels within range)


%% find minimum matched data (min)
% (only keep the closest match unless there are more than one in same distance)

hrtf_1_min_val = accumarray(matched(:,1),matched(:,3), [], @(x) {x});  
% group the distances by hrtf 1 into cell array
hrtf_1_min_idx = accumarray(matched(:,1),matched(:,3), [], @(x) {find(x == min(x))});
% find which the location of the minimum distance
hrtf_1_min_loc = accumarray(matched(:,1),matched(:,2), [], @(x) {x});
% group hrtf 2 by hrtf 1 into cell array
% (the order equivalent to hrtf_1_min_val)

hrtf_1_min_detail = [hrtf_1_min_val hrtf_1_min_idx hrtf_1_min_loc];
% group the distance cell array, minimum distance location and hrtf 2 cell
% array together

row_no = linspace(1, size(hrtf_1_min_detail,1), size(hrtf_1_min_detail,1));
hrtf_1_min_detail = [num2cell(row_no') hrtf_1_min_detail];
% add the row number (hrtf 1 data number) number into cell array

hrtf_1_min_no_empty = hrtf_1_min_detail(~cellfun(@isempty,hrtf_1_min_detail(:,2)),:);
% remove empty cell

hrtf_1_new = hrtf_1_min_no_empty(:,1);
hrtf_1_dist = hrtf_1_min_no_empty(:,1);
% initialise array
for n = 1:size(hrtf_1_min_no_empty,1)   
    hrtf_1_new{n} = hrtf_1_min_no_empty{n,4}(hrtf_1_min_no_empty{n,3});
	% find the hrtf 2 that has minimum distance 
	hrtf_1_dist{n} = hrtf_1_min_no_empty{n,2}(hrtf_1_min_no_empty{n,3});
	% record minimum distance
end

hrtf_1_min_no_empty = [hrtf_1_min_no_empty hrtf_1_new hrtf_1_dist];
% add minimum distance hrtf 2 and minium distance to matrix



hrtf_2_min_val = accumarray(matched(:,2),matched(:,3), [], @(x) {x});
% group the distances by hrtf 2 into cell array
hrtf_2_min_idx = accumarray(matched(:,2),matched(:,3), [], @(x) {find(x == min(x))});
% find which the location of the minimum distance
hrtf_2_min_loc = accumarray(matched(:,2),matched(:,1), [], @(x) {x});
% group hrtf 1 by hrtf 2 into cell array
% (the order equivalent to hrtf_1_min_val)

hrtf_2_min_detail = [hrtf_2_min_val hrtf_2_min_idx hrtf_2_min_loc];
% group the distance cell array, minimum distance location and hrtf 2 cell
% array together

row_no = linspace(1, size(hrtf_2_min_detail,1), size(hrtf_2_min_detail,1));
hrtf_2_min_detail = [num2cell(row_no') hrtf_2_min_detail];
% add the row number (hrtf 2 data number) number into cell array

hrtf_2_min_no_empty = hrtf_2_min_detail(~cellfun(@isempty,hrtf_2_min_detail(:,2)),:);
% remove empty cell

hrtf_2_new = hrtf_2_min_no_empty(:,1);
hrtf_2_dist = hrtf_2_min_no_empty(:,1);
% initialise array
for n = 1:size(hrtf_2_min_no_empty,1)   
	hrtf_2_new{n} = hrtf_2_min_no_empty{n,4}(hrtf_2_min_no_empty{n,3});
	% find the hrtf 2 that has minimum distance 
	hrtf_2_dist{n} = hrtf_2_min_no_empty{n,2}(hrtf_2_min_no_empty{n,3});
	% record minimum distance 
end

hrtf_2_min_no_empty = [hrtf_2_min_no_empty hrtf_2_new hrtf_2_dist];
% add minimum distance hrtf 2 and minium distance to matrix


hrtf_1_min = {0 0 0};
% initialise parameter (will remove later)
for n = 1: size(hrtf_1_min_no_empty, 1)
	if length(hrtf_1_min_no_empty{n,5}) > 1
        for m = 1: length(hrtf_1_min_no_empty{n,5})
            hrtf_1_min = [hrtf_1_min; hrtf_1_min_no_empty(n,1) hrtf_1_min_no_empty{n,5}(m) hrtf_1_min_no_empty{n,6}(m)];    
        end
    else
        hrtf_1_min = [hrtf_1_min; hrtf_1_min_no_empty(n,1) hrtf_1_min_no_empty(n,5) hrtf_1_min_no_empty(n,6)];
	end
end
hrtf_1_min = hrtf_1_min(2:size(hrtf_1_min, 1), :);
% remove unwanted detail and duplicate the ones with some distances

hrtf_2_min = {0 0 0};
% initialise parameter (will remove later)
for n = 1: size(hrtf_2_min_no_empty, 1)
	if length(hrtf_2_min_no_empty{n,5}) > 1
        for m = 1: length(hrtf_2_min_no_empty{n,5})
            hrtf_2_min = [hrtf_2_min; hrtf_2_min_no_empty(n,1) hrtf_2_min_no_empty{n,5}(m) hrtf_2_min_no_empty{n,6}(m)];    
        end
    else
        hrtf_2_min = [hrtf_2_min; hrtf_2_min_no_empty(n,1) hrtf_2_min_no_empty(n,5) hrtf_2_min_no_empty(n,6)];
    end
end
hrtf_2_min = hrtf_2_min(2:size(hrtf_2_min, 1), :);
% remove unwanted detail and duplicate the ones with some distances


hrtf_1_min = cell2mat(hrtf_1_min);
hrtf_2_min = cell2mat(hrtf_2_min);
% change cell array to normal numeric array

hrtf_1vs2_min = [[hrtf_1_min;zeros(size(hrtf_2_min,1)-size(hrtf_1_min,1),3)] ...
    [hrtf_2_min;zeros(size(hrtf_1_min,1)-size(hrtf_2_min,1),3)]];
% compare hrtf_1_min and hrtf_2_min by putting them into a singel matrix
% note that the column 1, 2 and 4, 5 are inverted 
% which is column 1, 5 represent hrtf 1, column 2, 4 represent hrtf 2

hrtf_min_match = intersect([hrtf_1vs2_min(:,1) hrtf_1vs2_min(:,2) hrtf_1vs2_min(:,3)], ...
	[hrtf_1vs2_min(:,5) hrtf_1vs2_min(:,4) hrtf_1vs2_min(:,6)], 'rows');
% keep the common result between hrtf_1_min and hrtf_2_min
% this will be the output of minimum matched datas
% column 1 repesnet hrtf 1, column 2 repesnet hrtf 2 and column 3 is the
% distance between the two data locations

if ~strcmp('fast', method)
    if length(hrtf_min_match(:,1)) ~= length(unique(hrtf_min_match(:,1)))
        warning('Duplicate in hrtf 1 (min): probably because there are more than one match with same distance')
    end
    if length(hrtf_min_match(:,2)) ~= length(unique(hrtf_min_match(:,2)))
        warning('Duplicate in hrtf 2 (min): probably because there are more than one match with same distance')
    end
    % disply warning if there is more than one match in output
end

%%  choose method (max or min) and plot result 

if strcmp('max', method) || strcmp('distance max',method)
    hrtf1_data_no = hrtf_max_match(:,1);
    hrtf2_data_no = hrtf_max_match(:,2);  
    
    HRTF_data_desire = hrtf_max_match;
    HRTF_data_opposite = hrtf_min_match;
else   
    hrtf1_data_no = hrtf_min_match(:,1);
    hrtf2_data_no = hrtf_min_match(:,2);
    
    HRTF_data_desire = hrtf_min_match;
    HRTF_data_opposite = hrtf_max_match;
end
% choose method

HRTF_1_common_angles = hrtf_1_angle(hrtf1_data_no, :);
HRTF_2_common_angles = hrtf_2_angle(hrtf2_data_no, :);
% keep common angles from orginal angles matrix

%{
HRTF_1_common_angles = [HRTF_1_common_angles HRTF_data_desire(:, 1)];
HRTF_2_common_angles = [HRTF_2_common_angles HRTF_data_desire(:, 2)];
%}

%% plot data

if plot ~= 0 && ~isempty(HRTF_data_desire)
%    SOFAplotGeometry(hrtf_1, idx_1(hrtf1_data_no)');
%    SOFAplotGeometry(hrtf_2, idx_2(hrtf2_data_no)');
    % plot common angles in HRTF 1 and HRTF 2 (revese back to non-sorted for plot only)

    figure
    if strcmp('true distance', plot_method) 
        plot_3d_angles( HRTF_1_common_angles(:,1), HRTF_1_common_angles(:,2), HRTF_1_common_angles(:,3),'+', 'MarkerEdgeColor', 'b')
        hold on
        plot_3d_angles( HRTF_2_common_angles(:,1), HRTF_2_common_angles(:,2),  HRTF_2_common_angles(:,3), 'MarkerEdgeColor', 'r')
    else
        plot_3d_angles( HRTF_1_common_angles(:,1), HRTF_1_common_angles(:,2), 1.5,'+', 'MarkerEdgeColor', 'b')
        hold on
        plot_3d_angles( HRTF_2_common_angles(:,1), HRTF_2_common_angles(:,2),  1.5, 'MarkerEdgeColor', 'r')
    end
    
    if strcmp('fast', method) 
        title('fast estimate common HRTF angles')
    elseif strcmp('distance max', method) 
        title('common HRTF angles by great circle distance (max)')
    elseif strcmp('distance min', method) || strcmp('distance', method)
        title('common HRTF angles by great circle distance (min)')
    elseif strcmp('max', method) 
        title('common HRTF angles (max)')
    else 
        title('common HRTF angles (min)')
    end
    % plot common angle on one graph in same distance (z=1.5m)
end

if isempty(HRTF_data_desire)
    warning('No matching data')  
end
    
end



