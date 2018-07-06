function [ bad_data ] = check_SOFA( file, angle_range, dist_range, plot)
%CHECK_SOFA Summary of this function goes here
% 
% Check for abnormal hrtf measurement in the input sofa file 

% INPUT: 
% e.g. check_hrtf('MRT02.sofa', 0.01, 0.1, 0 );
%
% 1. file = 'MRT02.sofa' % sofa file
%    - example:
%      'MRT02.sofa' 'RIEC_hrir_subject_008.sofa' 'hrtf b_nh15.sofa' 'MRT04.sofa'
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
%  
%
% OUTPUT:
% 1. bad_data % a struct contains bad measurement from the input sofa file
%    bad_data.file: File_Name (SOFA file that contains problematic measurements)
%    bad_data.repeatAngles: Repeat_Angles (repeated measurement angles)
%    bad_data.distance: Inconsist_Dist (inconsistance measurement distance)
%    bad_data.angles: Asymm_Angles (Asymmetrical measurement angles on left and right)
%    bad_data.median: Asymm_Medium (Asymmetricalmeasurement angles on medium plane,
%                     compare front and back, removed top (90) and bottom (-90)) 


%% initialise input and output

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

hrtf = SOFAload(file);
% load sofa file

apparentSourceVector = SOFAcalculateAPV(hrtf);
% find angles

bad_data.file = file;
bad_data.repeatAngles = [];
bad_data.distance = [];
bad_data.angles = [];
bad_data.median = [];
% initialise error log 'bad_data' (struct)

%% check is there repeated measurement angles

[total_angles, I] = unique(apparentSourceVector,'rows');
% find unique angles (in rows)

if length(apparentSourceVector) ~= length(total_angles)
    warning(['repeated angles in file ' file])
    % print warning
    
    bad_data.repeatAngles = {setdiff(1:size(apparentSourceVector,1), I)'};
    % save repeat angles row number in 'bad_data' log
    
    return
end
% show warning and save the hrtf row number in 'bad_data' if there are repeated angles


%% check distance consistency

total_distance = uniquetol(round(apparentSourceVector(:,3), 4), dist_range);
% find out how many different distances in the sofa file within the
% range ('dist_range')

if length(total_distance) > 1 % if there is more than 1 distance
    distance_hist = zeros(size(total_distance));
    for n = 1:length(total_distance)
        distance_hist(n) = sum(apparentSourceVector(:,3) == total_distance(n));
    end
    [~, idx] = sort(distance_hist);
    total_distance = total_distance(idx);
    distance_hist = distance_hist(idx);
    % distance histogram (sort by the histogram, so outliner goes first)
    
    for n = 1:length(distance_hist)-1
        loc = find(apparentSourceVector(:,3)== total_distance(n));
    end
    % find outliner location
    
    message = sprintf(['inconsistent meausrement distance in file ' file ...
            ' \n' 'Outliners row number: \n']);
    loc_message = sprintf('%2d\n', loc);
    sldiagviewer.reportWarning([message loc_message]);
    % compose and print warning message
    
    bad_data.distance = {loc};
    % save location error row number in 'bad_data' log
    
    if plot ~= 0 
        bar(distance_hist, 0.5)
        plot_dist = cellstr(num2str(total_distance));
        set(gca,'xticklabel',plot_dist)
        text(1:length(distance_hist),distance_hist,num2str(distance_hist), ...
            'vert','bottom','horiz','center','Color','red','FontSize',14); 
        box off
        title('measurement distance histogram')
        xlabel('distance (m)')
        % plot histogram and set all the labels
    end

end

%% check measurement distribution

pos_APV = apparentSourceVector((apparentSourceVector(:,1)>0 + angle_range), :);
pos_APV((pos_APV(:, 1)== 0 & pos_APV(:, 2)== 90), :) = [];
neg_APV = apparentSourceVector((apparentSourceVector(:,1)< 0 - angle_range & ...
    apparentSourceVector(:,1)>-180 + angle_range), :);
neg_APV(:,1) = abs(neg_APV(:,1));
% separate the positive azi angle (left) and negative angle (right)
% and remove the angles on median plane (will check them in the next section)
% convert the negative angle to positive for better comparison 

[common_APV_1, common_APV_2] = common_angles_2_HRTFs(pos_APV(:,1:2), neg_APV(:,1:2),...
    angle_range, 0, 'fast', 0);
% find common angles between positive azi angle (left) and negative angle (right)
% with in a range

if length(common_APV_1) ~= length(pos_APV) || length(common_APV_1) ~= length(neg_APV)||...
        length(common_APV_2) ~= length(pos_APV) || length(common_APV_2) ~= length(neg_APV)
    % if there is outliners (which is not in the common angles
    
    warning('non-symmetric measurement distribution')
    
    pos_diff = setdiff(pos_APV(:,1:2),common_APV_1(:,1:2), 'rows');
    neg_diff = setdiff(neg_APV(:,1:2),common_APV_2(:,1:2), 'rows');
    % find out the outliners in positive azi angle (left) and negative angle (right)
    
    neg_diff_orig = neg_diff;
    neg_diff_orig(:,1) = neg_diff_orig(:,1) .* -1;
    % change the negative angle outliner angles back to negative
    % (was changed to positive value for easy comparison)
    
    [~,bad_data.angles] = intersect(apparentSourceVector(:,1:2), pos_diff(:,1:2),'rows');
    [~,temp_idx] = intersect(apparentSourceVector(:,1:2), neg_diff_orig(:,1:2),'rows') ;
    bad_data.angles = {[bad_data.angles; temp_idx]};
    % save row number in 'bad_data' log positive angle (left) then negative angle (right)
    
    if plot ~= 0 
        figure
        plot_3d_angles(pos_APV(:,1), pos_APV(:,2), pos_APV(:,3), 'o', 'b')
        hold on
        plot_3d_angles(neg_APV(:,1), neg_APV(:,2), neg_APV(:,3), 'x', 'r')
        % plot all angles (flipped the right to the left for easy comparison)
        hold on
        plot_3d_angles(pos_diff(:,1), pos_diff(:,2), pos_APV(1,3), 'o', 'b', ...
            'filled', 'MarkerFaceAlpha',.35)
        hold on
        plot_3d_angles(neg_diff(:,1), neg_diff(:,2), neg_APV(1,3), 's', 'r', ...
            'filled','MarkerFaceAlpha',.35)
        % plot outliners
        title('all measurement angles (flipped the right to the left for comparison)')
        legend({'left', 'right', 'left outliner','right outliner'}, ...
            'Location', 'eastoutside','FontSize',12 )
        set(gcf, 'Position', [300, 500, 880, 700])
        % add texts
    end
end

%% check median distribution

front_APV = apparentSourceVector((apparentSourceVector(:,1) < angle_range & ...
    apparentSourceVector(:,1) > -angle_range), :);
back_APV = apparentSourceVector((apparentSourceVector(:,1) < -180 + angle_range | ...
    apparentSourceVector(:,1) > 180 - angle_range), :);
back_APV(:, 1) = abs(back_APV(:, 1)) - 180;
% find angels on median plan and saparate them in to front and back
% convert the back angles to front for better comparison

median_ele_diff_APV = front_APV(~ismembertol(front_APV(:,2), back_APV(:,2), angle_range), 2);
median_ele_diff_APV = [median_ele_diff_APV back_APV(~ismembertol(back_APV(:,2), ...
    front_APV(:,2), angle_range), 2)];
% find outliners

median_ele_diff_APV((abs(median_ele_diff_APV) > 90 - angle_range) , :) = [];

if ~isempty(median_ele_diff_APV) % if there is outliners
    
    warning('non-symmetric elevation measurement distribution on median plane')
    
    [~, median_ele_diff_front_idx] = intersect(front_APV(:,2), median_ele_diff_APV);
    [~, median_ele_diff_back_idx] = intersect(back_APV(:,2), median_ele_diff_APV);
    % check whether the outline is from the front ot from the back
    
    [~, bad_data.median] = intersect(apparentSourceVector(:,2), median_ele_diff_APV);
    bad_data.median = {bad_data.median};
    % save row number in 'bad_data' log
    
    if plot ~= 0
        figure
        plot_3d_angles(front_APV(:,1), front_APV(:,2), front_APV(:,3),'o','b','LineWidth', 1)
        hold on
        plot_3d_angles(back_APV(:,1), back_APV(:,2), back_APV(:,3),'x','r','LineWidth', 1.25)
        hold on
        % plot all median angles (flipped the back to the front for easy comparison)
        if ~isempty(median_ele_diff_front_idx)
            plot_3d_angles(front_APV(1,1), front_APV(median_ele_diff_front_idx, 2), ...
                front_APV(median_ele_diff_front_idx, 3), 'o', 'b', ...
                'filled', 'MarkerFaceAlpha',.35,'LineWidth', 1)
        else 
            plot_3d_angles(back_APV(1,1), back_APV(median_ele_diff_front_idx, 2), ...
                back_APV(median_ele_diff_front_idx, 3), 's', 'r', ...
                'filled','MarkerFaceAlpha',.35,'LineWidth', 1)
        end
        % plot outliner depends on it's location (from front or back)
        title('all median angles (flipped the back to front for comparison)')
        legend({'left', 'right', 'outliner'}, ...
            'Location', 'eastoutside','FontSize',12 )
        set(gcf, 'Position', [300, 500, 880, 700])
        view([0 0]) 
        % add text and adjust plot window
    end
    
end

end

