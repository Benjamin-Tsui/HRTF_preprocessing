function [P1_freq, N1_freq, N2_freq, P1_N1_amp_diff, P1_N2_amp_diff] = find_peak_and_notches(hrtf, Fs, plot_fft)
%FIND_PEAK_AND_NOTCHES Summary of this function goes here
% Extract hrtf peak and notches data 

% Input:
% hrtf = hrir left and right in column (n:2)
% Fs = hrir sample frequency
% plot_fft = plot trigger, 1 = plot, 0 = no plot
% (plow will show fft of the left channel)
% (will show P1, N1 and N2 search range with their locations)

% Output:
% P1_freq = peak between 3000-9000Hz
% N1_freq = first notch above peak below 17500Hz
% N2_freq = second notch above peak below 17500Hz
% (N1 is the one its frequency is closer to the peak)
% P1_N1_amp_diff = magnitude difference between P1 and N1
% P1_N2_amp_diff = magnitude difference between P1 and N2

if nargin == 2
     plot_fft = 0;
end
% Catch empty inputs

peak_range_lower = 3000;
peak_range_upper = 9000;
notches_range_upper = 17500;
% set range for peak and notches
% notches frequency should be above P1 
% notches_range_upper will be setted after P1 was found


hrtf = (hrtf ./ max(max(abs(hrtf)))) .* 0.99;
% normalise HRIR

fft_hrtf_left= abs(fft(hrtf(:,1))/length(hrtf(:,1)));
fft_hrtf_left = fft_hrtf_left(1:length(hrtf(:,1))/2+1);
fft_hrtf_left = 20 * log10(fft_hrtf_left);
% measurement frequency response (left)
fft_hrtf_right= abs(fft(hrtf(:,2))/length(hrtf(:,2)));
fft_hrtf_right = fft_hrtf_right(1:length(hrtf(:,2))/2+1);
fft_hrtf_right = 20 * log10(fft_hrtf_right);
% measurement frequency response (right)

if plot_fft == 1
    figure
    fft_hrtf_f = Fs * (0:(length(hrtf(:,1))/2))/length(hrtf(:,1)); % set plot frequency
    semilogx(fft_hrtf_f, fft_hrtf_left, 'Color', 'b')
    title('Measured HRTF');
    legend('measurement', 'Location','southwest');
    xlabel('frequency (Hz)');
    ylabel('magnitude')
    grid on
end
% plot fft

peak_range_lower = round(length(hrtf(:,1)) * peak_range_lower / Fs);
peak_range_upper = round(length(hrtf(:,1)) * peak_range_upper / Fs);
% convert to match fft length

[peak_amp_left, peak_loc_left] = findpeaks(fft_hrtf_left...
    (peak_range_lower :peak_range_upper), 'SortStr', 'descend');
[peak_amp_right, peak_loc_right] = findpeaks(fft_hrtf_right...
    (peak_range_lower :peak_range_upper), 'SortStr','descend');
% sort fft desceding in given range (max to min)

peak_freq_left = Fs .* (peak_loc_left + peak_range_lower - 2) ./ length(hrtf(:,1)) ;
peak_freq_right = Fs .* (peak_loc_right + peak_range_lower - 2) ./ length(hrtf(:,2));
% convert peak location to frequency

P1_amp_left = peak_amp_left(1);
P1_freq_left = peak_freq_left(1);
P1_amp_right = peak_amp_right(1);
P1_freq_right = peak_freq_right(1);
% find the highest peak (P1)



notches_range_lower_left = P1_freq_left;
notches_range_lower_right = P1_freq_right;
% make sure the notches is above P1

notches_range_lower_left = round(length(hrtf(:,1)) * notches_range_lower_left / Fs);
notches_range_lower_right = round(length(hrtf(:,1)) * notches_range_lower_right / Fs);
notches_range_upper = round(length(hrtf(:,1)) * notches_range_upper / Fs);
% convert to match fft length

peak_prominence = 30;
notch_no_left = 0 ;
while notch_no_left < 2 
    if peak_prominence < 9 && notch_no_left ==1
        break
    end
[notches_amp_left, notches_loc_left] = findpeaks(-fft_hrtf_left...
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
[notches_amp_right, notches_loc_right] = findpeaks(-fft_hrtf_right...
    (notches_range_lower_right: notches_range_upper), 'SortStr','descend', ...
    'MinPeakProminence', peak_prominence, 'MinPeakDistance', 10);

notches_amp_right = -notches_amp_right;

notch_no_right = length(notches_loc_right);
peak_prominence = peak_prominence - 1; 
end
% sort fft  in given range (min to max)



notches_freq_left = Fs .* (notches_loc_left + notches_range_lower_left - 2) ./ length(hrtf(:,1)) ;
notches_freq_right = Fs .* (notches_loc_right + notches_range_lower_right - 2) ./ length(hrtf(:,2));
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
        fft_hrtf_left(notches_range_lower_left: notches_range_upper) -1, ...
        'LineWidth', 3 , 'Color', 'g')
    semilogx(fft_hrtf_f(peak_range_lower: peak_range_upper),...
        fft_hrtf_left(peak_range_lower: peak_range_upper) +1, ...
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
end
% plot fft, peak search range, notch search range, P1, N1 and N2


P1_freq = [P1_freq_left P1_freq_right];
N1_freq = [N1_freq_left N1_freq_right];
N2_freq = [N2_freq_left N2_freq_right];
% P1, N1 and N2 frequency output

P1_N1_amp_diff_left = P1_amp_left - N1_amp_left;
P1_N1_amp_diff_right = P1_amp_right - N1_amp_right;
P1_N1_amp_diff = [P1_N1_amp_diff_left P1_N1_amp_diff_right];

P1_N2_amp_diff_left = P1_amp_left - N2_amp_left;
P1_N2_amp_diff_right = P1_amp_right - N2_amp_right;
P1_N2_amp_diff = [P1_N2_amp_diff_left P1_N2_amp_diff_right];
% magnitude difference between P1 and N1, P1 and N2



end

