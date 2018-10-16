function [HRTF1_selected_angles, HRTF2_selected_angles] = select_azi_angles...
    (HRTF_angle_ref_1, HRTF_angle_ref_2, range, plot_compare_result)
%SELECT_AZI_ANGLES Summary of this function goes here
%
% 1. Find common angle between two database by provided example
% 2. Find out the number of elevation samples in each azimuth angle
%    (for machine learning purpose, we want to pick the azimuth angles with
%    more elevation samples)
% 2a. A histogram will show the spread and count on each azimuth angles
% 2b. Users can choose to print the result in command window
%     (sound from left = positive (+), sound from right = negative (-))
% 3. Users then need to choose the azimuth angles they want to use by
%    inputing the desire angles in vector, e.g. [35 90 50 -175 0 55]
%    for angle 35, 90, 50, -175, 0 and 55, order does not matter
% 4. Based on the selected azimuth angle, the function then will choose the
%    common elevation angles between those azimuth angle
% 5. The function will output a list of azimuth and elevtion angle pairs 
%    for each database (column 1 is azimuth and column 2 is elevtion angle)


% HRTF_angle_ref_1 = 'hrtf b_nh2.sofa';  % a HRTF sample from ARI database
% HRTF_angle_ref_2 = ''MRT01.sofa'';  % a HRTF sample from ITA database 
% range = 2.35  % range for matching common angles between datasets
% plot_compare_result = 0  % no plot for common angles


% output list of azimuth and elevtion angle pairs for each database:
% [HRTF1_selected_angles, HRTF2_selected_angles]
% column 1 is azimuth angle and column 2 is elevtion angle


[HRTF1_angles, HRTF2_angles] = common_angles_in_HRTFs(HRTF_angle_ref_1, HRTF_angle_ref_2, ...
    range, plot_compare_result);
% find common angle between 2 hrtf datasets by feeding in one example
% from each database

[~,idx,~] = unique(round(HRTF1_angles(:, 1)));
HRTF1_azi_angles = HRTF1_angles(idx, 1);    
[~,idx, ~] = unique(round(HRTF2_angles(:, 1)));
HRTF2_azi_angles = HRTF2_angles(idx, 1);     
% find the common azimuth angles
% note that the rsult is rouned off to the most common angle
% when fetching HRTFs, range should set between 0.5 - 1 

figure
subplot(2,1,1);
HRTF1_angles_h = histogram(round(HRTF1_angles(:, 1)), ...
    length(unique(round(HRTF1_angles(:, 1)))));
HRTF1_h_data = HRTF1_angles_h.Values;
title('HRTF 1')
xlabel('azimuth angle') 
ylabel('no. of common samples') 
grid on
grid minor
subplot(2,1,2);
HRTF2_angles_h = histogram(round(HRTF2_angles(:, 1)), ...
    length(unique(round(HRTF2_angles(:, 1)))));
HRTF2_h_data = HRTF2_angles_h.Values;
title('HRTF 2')
xlabel('azimuth angle') 
ylabel('no. of common samples') 
grid on
grid minor
% histogram shows the sample number of each azimuth angle
% for user to reference which angle has enough number of samples for
% further use

print_common_ang = sprintf('Do you want to check the common angles? (Y/n): ');
reenter = input(print_common_ang, 's');
if ~isempty(reenter)
    if reenter == 'Y' || reenter == 'y'
        fprintf('Common angles between the two input datasets as below: \n')
        fprintf(' HRTF1:  HRTF2:  no. of samples:\n')
        fprintf('%7.1f %7.1f %10d\n', [HRTF1_azi_angles, HRTF2_azi_angles, HRTF1_h_data'].')
    end
end
% User choose whether they want to see the common azimuth angles in the
% command window.


azi_angle_choice = input('Choose azimuth angle (vector): ');
% User select the azimuth angle thet want to extract in vector
% e.g. [35 90 50 -175 0]

picked_HRTF1_sample = nan(max(HRTF1_h_data), length(azi_angle_choice));
picked_HRTF2_sample = nan(max(HRTF2_h_data), length(azi_angle_choice));
for n = 1: length(azi_angle_choice)
    temp_HRTF1_angles = HRTF1_angles((round(HRTF1_angles(:, 1)) == azi_angle_choice(n)), 2)';
    picked_HRTF1_sample(1:length(temp_HRTF1_angles), n) = temp_HRTF1_angles;
    
    temp_HRTF2_angles = HRTF2_angles((round(HRTF2_angles(:, 1)) == azi_angle_choice(n)), 2)';
    picked_HRTF2_sample(1:length(temp_HRTF2_angles), n) = temp_HRTF2_angles;
end
% find out all the elevation angles for each azimuth angle

HRTF1_common_ele = picked_HRTF1_sample;
HRTF2_common_ele = picked_HRTF2_sample;
for n = 1: length(azi_angle_choice)-1
    [~, idx, ~] = intersect(round(HRTF1_common_ele(:, 1)), round(HRTF1_common_ele(:, n+1)));
    temp_HRTF1_common_ele = HRTF1_common_ele(idx, 1); 
    HRTF1_common_ele(:,1) = [temp_HRTF1_common_ele; ...
        nan(size(HRTF1_common_ele,1)- length(temp_HRTF1_common_ele),1)];
    
    [~, idx, ~] = intersect(round(HRTF2_common_ele(:, 1)), round(HRTF2_common_ele(:, n+1)));
    temp_HRTF2_common_ele = HRTF2_common_ele(idx, 1); 
    HRTF2_common_ele(:,1) = [temp_HRTF2_common_ele; ...
        nan(size(HRTF2_common_ele,1)- length(temp_HRTF2_common_ele),1)];
end
HRTF1_common_ele = unique(HRTF1_common_ele(:,1));
HRTF2_common_ele = unique(HRTF2_common_ele(:,1));
% find common elevation angles in all azimuth angles

ind = ~isnan(HRTF1_common_ele); 
HRTF1_common_ele = HRTF1_common_ele(ind);
ind = ~isnan(HRTF2_common_ele); 
HRTF2_common_ele = HRTF2_common_ele(ind);
% remove Nan elements in vector

HRTF1_selected_angles = zeros(length(azi_angle_choice) * length(HRTF1_common_ele), 2);
HRTF2_selected_angles = zeros(length(azi_angle_choice) * length(HRTF2_common_ele), 2);
i = 1;
j = 1;
for n = 1:length(azi_angle_choice)
    for m = 1: length(HRTF1_common_ele)
        HRTF1_selected_angles(i, :) = [azi_angle_choice(n), HRTF1_common_ele(m)];
        i = i+1;
    end
    for m = 1: length(HRTF2_common_ele)
        HRTF2_selected_angles(j, :) = [azi_angle_choice(n), HRTF2_common_ele(m)];
        j = j+1;
    end
end
% list out the azimuth and elevation angle of the seleted samples
% column 1 is azimuth and column 2 is elevation

end

