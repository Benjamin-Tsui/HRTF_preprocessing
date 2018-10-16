function [ feature ] = feature_extraction( hrtf_or_hrir, hrtf_id_or_Fs, plot_fft, special_feature)
%FEATURE_EXTRACTION Summary of this function goes here
%   Detailed explanation goes here
% hrtf_or_hrir = 'hrtf b_nh2.sofa';
% or loaded sofa 
% hrtf_or_hrir = hrtf;
% hrtf_id = 71;
% plot_fft = 1;
% plot fft trigger
% special_feature = {'peaks and notches', 'octave_mean'}
% {'P1_freq','N1_freq','N2_freq','N1_freq', ...
%   'P1_N1_amp_diff','P1_N2_amp_diff','peaks and notches', ...
%   'third_octave_mean', 'third_octave_mean_dB',...
%   'octave_mean', 'octave_mean_dB'}

if nargin == 2
	plot_fft = 0;
    % plot fft trigger
    special_feature = 0;
end
% Catch empty inputs

if nargin == 3
    special_feature = 0;
end
% Catch empty inputs

if ischar(hrtf_or_hrir) % if input is a sofa file
    hrtf = SOFAload(hrtf_or_hrir);
    hrtf_id = hrtf_id_or_Fs;
    hrir = [squeeze(hrtf.Data.IR(hrtf_id, 1, :)), squeeze(hrtf.Data.IR(hrtf_id, 2, :))];
    Fs = hrtf.Data.SamplingRate;
    % load hrtf dataset
    
elseif isstruct(hrtf_or_hrir) % if input is a loaded sofa struct file
    hrtf = hrtf_or_hrir;
    hrtf_id = hrtf_id_or_Fs;
    hrir = [squeeze(hrtf.Data.IR(hrtf_id, 1, :)), squeeze(hrtf.Data.IR(hrtf_id, 2, :))];
    Fs = hrtf.Data.SamplingRate;
    % load hrtf dataset
    
else % if input is a hrir
    hrir = hrtf_or_hrir;
    Fs = hrtf_id_or_Fs;
end
% Input switch (sofa file or hrir)

hrtf_angles = SOFAcalculateAPV(hrtf);
azi = hrtf_angles(hrtf_id, 1);
ele = hrtf_angles(hrtf_id, 2);
distance = hrtf_angles(hrtf_id, 3);

peak_range_lower = 3000;
peak_range_upper = 9000;
notches_range_upper = 17500;
% set range for peak and notches
% notches frequency should be above P1 
% notches_range_upper will be setted after P1 was found


fft_hrtf_left= abs(fft(hrir(:,1))/length(hrir(:,1)));
fft_hrtf_left = fft_hrtf_left(1:length(hrir(:,1))/2+1);
fft_hrtf_left_dB = 20 * log10(fft_hrtf_left);
% measurement frequency response (left)
fft_hrtf_right= abs(fft(hrir(:,2))/length(hrir(:,2)));
fft_hrtf_right = fft_hrtf_right(1:length(hrir(:,2))/2+1);
fft_hrtf_right_dB = 20 * log10(fft_hrtf_right);
% measurement frequency response (right)

if plot_fft == 1
    figure 
    
    fft_hrtf_f = Fs * (0:(length(hrir(:,1))/2))/length(hrir(:,1)); % set plot frequency
    semilogx(fft_hrtf_f, fft_hrtf_left_dB, 'Color', 'b')
    title('Measured HRTF');
    legend('measurement', 'Location','southwest');
    xlabel('frequency (Hz)');
    ylabel('magnitude')
    grid on
end
% plot fft

%% Find peak and notches

if sum(strcmpi(special_feature,'P1_freq')) == 1 || sum(strcmpi(special_feature,'N1_freq')) == 1 || ...
        sum(strcmpi(special_feature,'N2_freq')) == 1 || sum(strcmpi(special_feature,'N1_freq')) == 1 || ...
        sum(strcmpi(special_feature,'P1_N1_amp_diff')) == 1 || ...
        sum(strcmpi(special_feature,'P1_N2_amp_diff')) == 1 || ...
        sum(strcmpi(special_feature,'peaks and notches')) == 1 

	peak_range_lower = round(length(hrir(:,1)) * peak_range_lower / Fs) + 1;
    peak_range_upper = round(length(hrir(:,1)) * peak_range_upper / Fs) + 1;
    % convert to match fft length

    [peak_amp_left, peak_loc_left] = findpeaks(fft_hrtf_left_dB...
        (peak_range_lower :peak_range_upper), 'SortStr', 'descend');
    [peak_amp_right, peak_loc_right] = findpeaks(fft_hrtf_right_dB...
        (peak_range_lower :peak_range_upper), 'SortStr','descend');
    % sort fft desceding in given range (max to min)

    if isempty(peak_amp_left)
        peak_amp_left = 0;
        peak_loc_left = 0;
    end

    if isempty(peak_amp_right)
        peak_amp_right = 0;
        peak_loc_right = 0;
    end

    peak_freq_left = Fs .* (peak_loc_left + peak_range_lower - 2) ./ length(hrir(:,1)) ;
    peak_freq_right = Fs .* (peak_loc_right + peak_range_lower - 2) ./ length(hrir(:,2));
    % convert peak location to frequency

    P1_amp_left = peak_amp_left(1);
    P1_freq_left = peak_freq_left(1);
    P1_amp_right = peak_amp_right(1);
    P1_freq_right = peak_freq_right(1);
    % find the highest peak (P1)

    notches_range_lower_left = P1_freq_left;
    notches_range_lower_right = P1_freq_right;
    % make sure the notches is above P1

    notches_range_lower_left = round(length(hrir(:,1)) * notches_range_lower_left / Fs) + 1;
    notches_range_lower_right = round(length(hrir(:,1)) * notches_range_lower_right / Fs) + 1;
    notches_range_upper = round(length(hrir(:,1)) * notches_range_upper / Fs) + 1;
    % convert to match fft length

    peak_prominence = 30;
    notch_no_left = 0 ;
    while notch_no_left < 2 
        if peak_prominence < 9 && notch_no_left ==1
            break
        end
    [notches_amp_left, notches_loc_left] = findpeaks(-fft_hrtf_left_dB...
        (notches_range_lower_left: notches_range_upper), 'SortStr','descend', ...
        'MinPeakProminence', peak_prominence, 'MinPeakDistance', 10);
    notches_amp_left = -notches_amp_left;

    notch_no_left = length(notches_loc_left);
    peak_prominence = peak_prominence - 1; 
    end

    peak_prominence = 30;
    notch_no_right = 0 ;
    while notch_no_right < 2
        if peak_prominence < 9 && notch_no_right ==1
            break
        end
    [notches_amp_right, notches_loc_right] = findpeaks(-fft_hrtf_right_dB...
        (notches_range_lower_right: notches_range_upper), 'SortStr','descend', ...
        'MinPeakProminence', peak_prominence, 'MinPeakDistance', 10);

    notches_amp_right = -notches_amp_right;

    notch_no_right = length(notches_loc_right);
    peak_prominence = peak_prominence - 1; 
    end
    % sort fft  in given range (min to max)



    notches_freq_left = Fs .* (notches_loc_left + notches_range_lower_left - 2) ./ length(hrir(:,1)) ;
    notches_freq_right = Fs .* (notches_loc_right + notches_range_lower_right - 2) ./ length(hrir(:,2));
    % convert peak location to frequency

    if length(notches_amp_left) > 1
        if notches_freq_left(1) <= notches_freq_left(2)
            N1_amp_left = notches_amp_left(1);
            N1_freq_left = notches_freq_left(1);
            % find the first notch (N1)
            N2_amp_left = notches_amp_left(2);
            N2_freq_left = notches_freq_left(2);
            % find the second notch (N2)
        else
            N1_amp_left = notches_amp_left(2);
            N1_freq_left = notches_freq_left(2);
            % find the first notch (N1)
            N2_amp_left = notches_amp_left(1);
            N2_freq_left = notches_freq_left(1);
            % find the second notch (N2)
        end
        % make sure N1 is always the one closer to P1
    else
        N1_amp_left = notches_amp_left(1);
        N1_freq_left = notches_freq_left(1);
        % find the first notch (N1)
        N2_amp_left = notches_amp_left(1);
        N2_freq_left = notches_freq_left(1);
        % find the first notch (N1)
    end  
    % catch error (when there is only one main notch)

    if length(notches_amp_right) > 1
        if notches_freq_right(1) <= notches_freq_right(2)
            N1_amp_right = notches_amp_right(1);
            N1_freq_right = notches_freq_right(1);
            % find the first notch (N1)
            N2_amp_right = notches_amp_right(2);
            N2_freq_right = notches_freq_right(2);
            % find the second notch (N2)
        else
            N1_amp_right = notches_amp_right(2);
            N1_freq_right = notches_freq_right(2);
            % find the first notch (N1)
            N2_amp_right = notches_amp_right(1);
            N2_freq_right = notches_freq_right(1);
            % find the second notch (N2)
        end
        % make sure N1 is always the one closer to P1
    else
        N1_amp_right = notches_amp_right(1);
        N1_freq_right = notches_freq_right(1);
        % find the first notch (N1)
        N2_amp_right = notches_amp_right(1);
        N2_freq_right = notches_freq_right(1);
        % find the first notch (N1)
    end  
    % catch error (when there is only one main notch)

    if plot_fft == 1
        hold on
        semilogx(fft_hrtf_f(notches_range_lower_left: notches_range_upper),...
            fft_hrtf_left_dB(notches_range_lower_left: notches_range_upper) -1, ...
            'LineWidth', 3 , 'Color', 'g')
        semilogx(fft_hrtf_f(peak_range_lower: peak_range_upper),...
            fft_hrtf_left_dB(peak_range_lower: peak_range_upper) +1, ...
            'LineWidth', 2 , 'Color', 'r')
    
        semilogx(N1_freq_left , N1_amp_left - 5, 'LineWidth', 1.5, ...
            'Marker', '^', 'MarkerSize', 5, 'Color', 'k', 'MarkerFaceColor', 'k');
        text(N1_freq_left, N1_amp_left - 13, ' N1', 'HorizontalAlignment', 'center')
        semilogx(N2_freq_left , N2_amp_left - 5, 'LineWidth', 1.5, ...
            'Marker', '^', 'MarkerSize', 5, 'Color', 'k', 'MarkerFaceColor', 'k');
        text(N2_freq_left , N2_amp_left - 13, ' N2', 'HorizontalAlignment', 'center')
        semilogx(P1_freq_left , P1_amp_left + 5, 'LineWidth', 1.5, ...
            'Marker', 'v', 'MarkerSize', 5, 'Color', 'k', 'MarkerFaceColor', 'k');
        text(P1_freq_left, P1_amp_left + 13, ' P1', 'HorizontalAlignment', 'center')
    
        hold off
    end
    % plot fft, peak search range, notch search range, P1, N1 and N2

    feature.P1_freq = [P1_freq_left P1_freq_right];
    feature.N1_freq = [N1_freq_left N1_freq_right];
    feature.N2_freq = [N2_freq_left N2_freq_right];
    % P1, N1 and N2 frequency output

    P1_N1_amp_diff_left = P1_amp_left - N1_amp_left;
    P1_N1_amp_diff_right = P1_amp_right - N1_amp_right;
    feature.P1_N1_amp_diff = [P1_N1_amp_diff_left P1_N1_amp_diff_right];

    P1_N2_amp_diff_left = P1_amp_left - N2_amp_left;
    P1_N2_amp_diff_right = P1_amp_right - N2_amp_right;
    feature.P1_N2_amp_diff = [P1_N2_amp_diff_left P1_N2_amp_diff_right];
    % magnitude difference between P1 and N1, P1 and N2
end

%% Mean in frequency bands

if sum(strcmpi(special_feature,'third_octave_mean')) == 1 || ...
        sum(strcmpi(special_feature,'third_octave_mean_dB')) == 1 || ...
        sum(strcmpi(special_feature,'octave_mean')) == 1 || ...
        sum(strcmpi(special_feature,'octave_mean_dB')) == 1

    third_octave_value = frequency_bnad_calculation(20, 20000, 'third-octave');
    octave_value = frequency_bnad_calculation(20, 20000, 'octave');
    % find octave bands frequency value (row1 = start, row2 = mid, row3 = end)

    mean_third_octave_left = zeros(size(third_octave_value, 1), 1);
    mean_third_octave_right = zeros(size(third_octave_value, 1), 1);
    mean_third_octave_dB_left = zeros(size(third_octave_value, 1), 1);
    mean_third_octave_dB_right = zeros(size(third_octave_value, 1), 1);

    mean_octave_left = zeros(size(octave_value, 1), 1);
    mean_octave_right = zeros(size(octave_value, 1), 1);
    mean_octave_dB_left = zeros(size(octave_value, 1), 1);
    mean_octave_dB_right = zeros(size(octave_value, 1), 1);
    % initialise output matrix

    for n = 1 :size(third_octave_value, 1)
        band_pass_lower = round(length(hrir(:,1)) * third_octave_value(n,1) / Fs) + 1;
        band_pass_upper = round(length(hrir(:,1)) * third_octave_value(n,3) / Fs) + 1;
        
        if band_pass_upper > length(fft_hrtf_left)
            mean_third_octave_left(n) = mean(fft_hrtf_left(band_pass_lower: end));
            mean_third_octave_right(n) = mean(fft_hrtf_right(band_pass_lower: end));
    
            mean_third_octave_dB_left(n) = mean(fft_hrtf_left_dB(band_pass_lower: end));
            mean_third_octave_dB_right(n) = mean(fft_hrtf_right_dB(band_pass_lower: end)); 
        else
            mean_third_octave_left(n) = mean(fft_hrtf_left(band_pass_lower: band_pass_upper));
            mean_third_octave_right(n) = mean(fft_hrtf_right(band_pass_lower: band_pass_upper));
    
            mean_third_octave_dB_left(n) = mean(fft_hrtf_left_dB(band_pass_lower: band_pass_upper));
            mean_third_octave_dB_right(n) = mean(fft_hrtf_right_dB(band_pass_lower: band_pass_upper));
        end
    end

    for n = 1 :size(octave_value, 1)
        band_pass_lower = round(length(hrir(:,1)) * octave_value(n,1) / Fs) + 1;
        band_pass_upper = round(length(hrir(:,1)) * octave_value(n,3) / Fs) + 1;
    
        if band_pass_upper > length(fft_hrtf_left)
            mean_octave_left(n) = mean(fft_hrtf_left(band_pass_lower: end));
            mean_octave_right(n) = mean(fft_hrtf_right(band_pass_lower: end));
    
            mean_octave_dB_left(n) = mean(fft_hrtf_left_dB(band_pass_lower: end));
            mean_octave_dB_right(n) = mean(fft_hrtf_right_dB(band_pass_lower: end));
        else
            mean_octave_left(n) = mean(fft_hrtf_left(band_pass_lower: band_pass_upper));
            mean_octave_right(n) = mean(fft_hrtf_right(band_pass_lower: band_pass_upper));
    
            mean_octave_dB_left(n) = mean(fft_hrtf_left_dB(band_pass_lower: band_pass_upper));
            mean_octave_dB_right(n) = mean(fft_hrtf_right_dB(band_pass_lower: band_pass_upper));
        end
    end


    feature.third_octave_mean = [mean_third_octave_left mean_third_octave_right];
    % mean megnitude in one-third ovtave bands

    feature.third_octave_mean_dB = [mean_third_octave_dB_left mean_third_octave_dB_right];
    % mean dB in one-third ovtave bands

    feature.octave_mean = [mean_octave_left mean_octave_right];
    % mean megnitude in ovtave bands

    feature.octave_mean_dB = [mean_octave_dB_left mean_octave_dB_right];
    % mean dB in ovtave bands

end

%% itd
[~, location_L] = max(hrir(:,1));
[~, location_R] = max(hrir(:,2));
ITD_sample = location_L - location_R;    % negative value if sound from the left (delay on right)
itd = ITD_sample/Fs;

%% Output features

feature.azi = azi;
feature.ele = ele;
feature.distance = distance;

feature.itd = itd;
% ITD

feature.ild = fft_hrtf_left(2:end) - fft_hrtf_right(2:end);
% ILD

feature.ild_db = fft_hrtf_left_dB(2:end) - fft_hrtf_right_dB(2:end);
% ILD in log scale (dB)

feature.hrtf = [fft_hrtf_left(2:end) fft_hrtf_right(2:end)];
% hrtf

feature.log_hrtf = [fft_hrtf_left_dB(2:end) fft_hrtf_right_dB(2:end)];
% hrtf in log scale (dB)

feature.hrir = hrir;
% hrir


end

