function [ Obj ] = normalise_hrir(input, hrtf_id, target_length, target_fs, ...
    target_amplitude, SOFAdir, SOFAfilename, plot_trig, print)
%NORMALISE Summary of this function goes here
%   Normalise hrtf and save in SOFA file format, user can normalise the input
%   SOFA file's hrtf to desire length, sampling frequency and maximum amplitude.
%   - user can also decide the file name and save location if they want to
%   - or just output the struct but not saving in sofa file
%   - to check the result, different graph can be plotted, include:
%     'hrtf_waterfall', 'hrtf_heatmap', 'hrir_waterfall', 'hrir_heatmap'
%     or compare the original and normalised of individual hrir 
%

% e.g. normalise_hrir( 'irc_1004.sofa', [1 2 3 4 5] , 1024, 44100, ...
%      [], 'normalised_SOFA/', [], [1 56 78 18 109]);

% INPUT: 
% 1. input = 'irc_1004.sofa' % original SOFA file
%
% 2. hrtf_id = [1 2 3 4 5 6 7 8 9 10]   % sample number vector
%    - hrtf measurement number (row number) from the original SOFA file
%    - if empty (e.g. []), equal to all measurements
%
% 3. target_length = 1024;   % new hrir length
%    - origainl length if empty or missing (target_length = [])
%
% 4. target_fs = 44100;   % new sampling frequency
%    - origainl fs if empty  or missing (target_fs = [])
%
% 5. target_amplitude = 0.99;   % new sampling frequency (OPTIONAL)
%     - default 0.99 if empty or missing
%
% 6. SOFAdir = 'normalised_SOFA/'   % new file save location (OPTIONAL)
%    - in char
%    - make sure the folder exist
%    - add / in the end
%    - save in current folder if empty or leave blank
%
% 7. SOFAfilename = []; or desire file name e.g. 'new_file.sofa'; (OPTIONAL)
%    - in char
%    - add .sofa in the end
%    - save in original input file name + '_normalised.sofa' if empty or leave blank
%    - 0 : no saving SOFA file
%    
% 8. plot_trig = 0   % trigger plot(OPTIONAL)
%
%    - 0 or missing: no plot
%    OR 
%    - 'all': plot everything 'hrtf_waterfall', 'hrtf_heatmap', 'hrir_waterfall', 
%       'hrir_heatmap' plus the 1st hrir comparison (may time a while to process)
%    OR 
%    - {'hrtf_waterfall', 'hrtf_heatmap'}
%      - optional plot in cell array (can be single, e.g. {'hrtf_waterfall'} 
%        or multiple e.g. {'hrtf_waterfall', 'hrtf_heatmap'})
%    OR
%    - [1 3 4 5] : any number or numeric array
%
%    - Plot options:
%      'hrtf_waterfall', 'hrtf_heatmap', 'hrir_waterfall', 'hrir_heatmap', [1 3 4 5]
% 
%      - 'hrtf_waterfall': paired original and normalised hrtf 3D waterfall bar plot
%      - 'hrtf_heatmap': paired original and normalised hrtf heatmap
%        (similar to the top view of the hrtf_waterfall plot)
%      - 'hrir_waterfall': paired original and normalised hrir 3D waterfall bar plot
%      - 'hrir_heatmap': paired original and normalised hrir heatmap
%        (similar to the top view of the hrtf_waterfall plot)
%      - [1 3 4 5] : any number or numeric array 
%        - plot the selected original and normalised hrir on the same figure 
%          to compare the result
%        - NOTE THAT the numbe is base on the original input NOT the new result
%          - if hrtf_id = [4 6 10 41 56 129], input [4 6] will NOT plot 41 
%            and 129, it will plot 4 and 6 instead
%          - if hrtf_id = [] (all measurements), input [4 6] will pot 4 and 6
%            from the data
%      - default 0 (no plot)
%
% 9. print = 1;  % print trigger
%    print out the progress 
%    - default(if empty or missing) = 1
%    

% OUTPUT:
% Obj: new SOFA struct (same as SOFAload('new_SOFA_file.sofa')
%
%

%% pre-processing

SOFA_hrtf = SOFAload(input);
% read sofa file


if nargin == 1
     hrtf_id = [];
     target_length = size(SOFA_hrtf.Data.IR, 3);
     target_fs = SOFA_hrtf.Data.SamplingRate;
     target_amplitude = 0.99;
     SOFAdir = [];
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 2
     target_length = size(SOFA_hrtf.Data.IR, 3);
     target_fs = SOFA_hrtf.Data.SamplingRate;
     target_amplitude = 0.99;
     SOFAdir = [];
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 3
     target_fs = SOFA_hrtf.Data.SamplingRate;
     target_amplitude = 0.99;
     SOFAdir = [];
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 4 
     target_amplitude = 0.99;
     SOFAdir = [];
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 5 
     SOFAdir = [];
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 6 
     SOFAfilename = [];
     plot_trig = 0;
     print = 1;
elseif nargin == 7 
     plot_trig = 0;  
     print = 1;
elseif nargin == 8 
     print = 1;
end
% catch missing inputs


if isempty(hrtf_id)
    hrtf_id = linspace(1,size(SOFA_hrtf.Data.IR, 1),size(SOFA_hrtf.Data.IR, 1));
    % take all hrtf (normalise all hrtf)
end

if isempty(target_length)
    target_length = size(SOFA_hrtf.Data.IR, 3);
    % default target length = same as original input
end

if isempty(target_fs)
    target_fs = SOFA_hrtf.Data.SamplingRate;
    % default target sampling rate = same as original input
end

if isempty(target_amplitude)
    target_amplitude = 0.99;
    % default target amplitude
end

if isempty(print)
    print = 1;
    % default print
end
% catch empty inputs

reshape(hrtf_id,1,[]);
% reshape hrtf_id to horizontal vector if it is a vertical vector


hrtf_angles = SOFAcalculateAPV(SOFA_hrtf);
% get angles

hrir_left = squeeze(SOFA_hrtf.Data.IR(:, 1, :));
hrir_right = squeeze(SOFA_hrtf.Data.IR(:, 2, :));
% extract hrir from sofa file

Fs = SOFA_hrtf.Data.SamplingRate;
% sampling rate

hrir = [hrir_left; hrir_right];
% put left and right hrir into a single matrix


%% resample input data (if neccessery)

if Fs ~= target_fs
    hrir = resample(hrir',target_fs,Fs)';
end
% resample the data to desire Fs, if neccessery


%% fft original and desire hrtf (find desire dB in hrtf)

org_amp_max = max(max(hrir));
% find original maximum amplitude

desire_hrir = hrir / org_amp_max * target_amplitude;
% amplify the orginal hrir to the target amplitude 
% e.g if target amplitude = 0.99, result will be -0.99 to 0.99
fft_desire_hrtf= fft(desire_hrir,target_length,2)/length(desire_hrir(1,:));
% use fft to change it in frequency domain (hrtf)
fft_desire_hrtf_nyquist = (fft_desire_hrtf(:, 1:target_length/2));
% remove mirror part of the fft result
fft_desire_hrtf_nyquist_dB = 20 .* log10(fft_desire_hrtf_nyquist);
% convert it to dB
desire_hrtf_max = max(max(fft_desire_hrtf_nyquist_dB));
% find maxium dB value

hrtf = fft(hrir,[],2);
% find frequency response (hrtf)of the original signal using fft

fft_hrtf = fft(hrir,target_length,2)/length(hrir(1,:));
% find frequency response (hrtf)of the original signal using fft (with
% desire output length)
%fft_hrtf = abs(fft_hrtf);
fft_hrtf_nyquist = fft_hrtf(:, 1:target_length/2);
% remove mirror part of the fft result
fft_hrtf_nyquist_dB = 20 .* log10(fft_hrtf_nyquist);
% convert the nyquist result to dB 
fft_hrtf_dB =  20 .* log10(fft_hrtf);
% convert the original result to dB


%% normalise

amp_max = max(max(fft_hrtf_nyquist_dB));
% find maximum dB value in the nyquist hrtf
amp_diff = desire_hrtf_max - amp_max;
% find the difference between the maximum of the desire hrtf and nyquist hrtf
norm_hrtf_dB = fft_hrtf_dB + amp_diff ;
% add the dufferent to the original hrtf to normaise result
norm_hrtf = (10.^(norm_hrtf_dB ./ 20));
% change it back to magnitude 
norm_hrir = ifft(norm_hrtf,target_length,2)* length(desire_hrir(1,:));
% inverse fft the normalised result
norm_hrir = real(norm_hrir);
% take only real value of the result

org_amp_max = max(max(hrir));
norm_amp_max = max(max(norm_hrir));

if print ~= 0
    disp(['original hrir maximum amplitude: ' num2str(org_amp_max)])
    disp(['normalised hrir maximum amplitude: ' num2str(norm_amp_max)])
    % disp maximum amplitude of the original hrir and the normalised hrir
end


%% saperate result to left and right

norm_hrir_left = norm_hrir(1: size(norm_hrir, 1) / 2, :);
norm_hrir_right = norm_hrir(size(norm_hrir, 1) /2 + 1 : end, :);
% saperate the result back to left and right individually 

norm_hrir_left = norm_hrir_left(hrtf_id, :);
norm_hrir_right = norm_hrir_right(hrtf_id, :);


%% plot normalised result

if ischar(plot_trig)
    plot_trig = {plot_trig};
elseif isnumeric(plot_trig) && plot_trig ~= 0
    plot_trig = {plot_trig};
end
% catch char

if ~iscell(plot_trig)
    if plot_trig == 0
        plot_trig = {};
    else
        error('please input cell array or 0 (no plot) or 1 (plot all)')
    end
end
% catch non cell array

if ~isempty(plot_trig) || iscell(plot_trig)
    
    if strcmpi(plot_trig, 'all') 
        plot_trig = {'hrtf_waterfall', 'hrtf_heatmap', 'hrir_waterfall', ... 
            'hrir_heatmap' 1};
    end
    
    
	if sum(strcmpi(plot_trig, 'hrtf_waterfall')) > 0        
        figure
        b= bar3(real(fft_hrtf_nyquist_dB) - min(min(real(fft_hrtf_nyquist_dB))), 1);
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
            b(k).EdgeAlpha = 0.35;
        end
        view(0, 15)
        set(gca,'xscale','log')
        axis fill
        title('original hrtf')
        set(gcf, 'Position', get(0, 'Screensize'));
        % plot original hrtf waterfall 3d bar chart (added 300 to make all
        % values positive)
    
        figure
        b= bar3(real(norm_hrtf_dB(:, 1:target_length/2)) - min(min(real(norm_hrtf_dB))), 1);
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
            b(k).EdgeAlpha = 0.35;
        end
        view(0, 15)
        set(gca,'xscale','log')
        axis fill
        title('normalised hrtf')
        set(gcf, 'Position', get(0, 'Screensize'));
        % plot normalised hrtf waterfall 3d bar chart (added 300 to make all
        % values positive)       
	end
    
    if sum(strcmpi(plot_trig, 'hrtf_heatmap')) > 0   
        figure
        h = pcolor(real(norm_hrtf_dB(:, 1:target_length/2)));
        set(gcf, 'Position', get(0, 'Screensize'));
        set(h, 'EdgeColor', 'none');
        title('normailsed hrir heat map')
        colorbar
        % plot normalised hrtf heatmap
    
        figure
        h = pcolor(real(fft_hrtf_nyquist_dB));
        set(gcf, 'Position', get(0, 'Screensize'));
        set(h, 'EdgeColor', 'none');
        title('original hrir heat map')
        colorbar
        % plot normalised hrtf heatmap
    end
    
    if sum(strcmpi(plot_trig, 'hrir_waterfall')) > 0
        figure
        b= bar3(hrir, 1);
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
            b(k).EdgeAlpha = 0.35;
        end
        view(0, 15)
        axis fill
        title('original hrir')
        set(gcf, 'Position', get(0, 'Screensize'));
        % plot original hrir waterfall 3d bar chart (added 300 to make all
        % values positive)
    
        figure
        b= bar3(norm_hrir, 1);
        for k = 1:length(b)
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
            b(k).EdgeAlpha = 0.35;
        end
        view(0, 15)
        axis fill
        title ('normalised hrir')
        set(gcf, 'Position', get(0, 'Screensize'));
        % plot normalised hrir waterfall 3d bar chart (added 300 to make all
        % values positive)
    end
 
    if sum(strcmpi(plot_trig, 'hrir_heatmap')) > 0
        figure
        h = pcolor(hrir);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(h, 'EdgeColor', 'none');
        title('original hrir heat map')
        colorbar
        % plot original hrir heatmap

        figure
        h = pcolor(norm_hrir);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(h, 'EdgeColor', 'none');
        title('normailsed hrir heat map')
        colorbar
        % plot normalised hrir heatmap 
    end
    
    if sum(cellfun(@(x) isnumeric(x), plot_trig)) > 0
        plot_id = plot_trig{cellfun(@(x) isnumeric(x), plot_trig)};
        %if max(plot_id) > norm_hrir
        
        for n = 1:length(plot_id)
            figure
            ax1 = subplot(2,1,1);
            plot_L1 = plot(hrir(n,:), 'LineWidth', 0.85);
            hold on
            plot_L2 = plot(ax1, norm_hrir(n,:));
            title(ax1, ['compare hrir number: ' num2str(plot_id(n)) ' (left)'])
            legend(ax1, 'original', 'normalised')
            plot_L1.Color(4) = 0.5;
            plot_L2.Color(4) = 0.9;
            ylim([-target_amplitude target_amplitude])
            % compare the selected hrir (left)
            
            ax2 = subplot(2,1,2);
            plot_R1 = plot(hrir(size(hrir, 1) / 2 + n,:), 'LineWidth', 0.85);
            hold on
            plot_R2 = plot(ax2, norm_hrir(size(norm_hrir, 1) / 2 + n,:));
            title(ax2, ['compare hrir number: ' num2str(plot_id(n)) ' (right)'])
            legend(ax2, 'original', 'normalised')
            plot_R1.Color(4) = 0.5;
            plot_R2.Color(4) = 0.9;
            ylim([-target_amplitude target_amplitude])
            % compare the selected hrir (right)
            
            set(gcf, 'Position', [150, 500 , 1200, 750]);
        end
    end

end


%% save as SOFA file

Obj = save_SOFA(input, norm_hrir_left, norm_hrir_right, target_length, ...
    target_fs, hrtf_id, SOFAdir, SOFAfilename, print);


end

