function [hrtf, Fs, row_numbers] = fetchHRTFbyAngle(HRTF_sofa, azi_angle, ele_angle, range, plot)
%FETCHHRTFBYANGLE Summary of this function goes here
% Get single HRTF from a SOFA file by inputing the azimuth and elevation
% angle, user can also set an range

%   Detailed explanation goes here
% e.g. [hrtf, hrtf_left, hrtf_right] = fetchHRTFbyAngle('hrtf b_nh2.sofa', 91, 44, 1, 1);
% HRTF_sofa = 'hrtf b_nh2.sofa';  % ARI database
% azi_angle = 91;  % from a listener point of view
% ele_angle = 44;  % from a listener point of view
% range = 1;  % set angle range on each side 
% (if range = 1, the actual range will be 2)

if ischar(HRTF_sofa)
    hrtf_sofa = SOFAload(HRTF_sofa);  % SOFA file name
else
    hrtf_sofa = HRTF_sofa;  % loaded SOFA variable (for fetching multiple angles)
end
% Load SOFA file or use loaded SOFA file

if nargin == 3
     range = 0;
     plot = 0;
end
if nargin == 4
     plot = 0;
end
% Catch empty inputs

Fs = hrtf_sofa.Data.SamplingRate;

sourceVector = SOFAcalculateAPV(hrtf_sofa);
% Calculate the source position from a listener point of view

if length(azi_angle) ~= length(ele_angle)
    error('azi_angle dimension and ele_angle dimension must agree')
end

%range = ones(length(azi_angle),1) .* range;
%row_numbers = zeros(length(azi_angle), 1);
k = 1;
for n = 1:length(azi_angle)
    [row, ~] = find((sourceVector(:,1) >= (azi_angle(n) - range) & ...
        sourceVector(:,1) <= (azi_angle(n) + range)  & ...
        sourceVector(:,2) >= (ele_angle(n) - range) & ...
        sourceVector(:,2) <= (ele_angle(n) + range)));
    % Find the row number represent the angle
    for m = 1:length(row)
        row_numbers(k) = row(m);
        k = k + 1;
    end
end
row_numbers = row_numbers';

if ~isempty(row_numbers)
    if plot == 1
        SOFAplotGeometry(hrtf_sofa, row_numbers);
    end
    
else
    disp('There is no HRTF was found, please try with another angle or a wider range');
    hrtf = 0;
end

if length(row_numbers) > length(azi_angle)
    disp('There are more than one HRTF was found, please try narrow down the range');
end 
% Get HRTF and plot angle(optional)

hrtf_left = squeeze(hrtf_sofa.Data.IR(row_numbers, 1, :));
hrtf_right = squeeze(hrtf_sofa.Data.IR(row_numbers, 2, :));
if length(row_numbers) > 1
    hrtf(:,1,:) = hrtf_left';
    hrtf(:,2,:) = hrtf_right';
else
    hrtf = [hrtf_left hrtf_right];
end

end

