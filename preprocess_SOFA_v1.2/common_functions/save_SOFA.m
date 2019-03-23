function [ Obj ] = save_SOFA( input_SOFA, hrir_left, hrir_right, target_length, ...
    target_fs, hrtf_id, SOFAdir, SOFAfilename, print )
%SAVE_SOFA Summary of this function goes here
%   Save the updated / normalised hrir in SOFA file format with updated
%   attributes.

% INPUT: 
% 1st e.g. save_SOFA( 'irc_1004.sofa', norm_hrir_left, norm_hrir_right, 1024, 44100, ...
%          [1 2 3], 'normalised_SOFA/', []);
%
% 2nd e.g. sofa_name = 'irc_1004.sofa';
%          selected_data = [1 3 4 5 18];
%
%          sofa_hrtf = SOFAload(sofa_name);
%          save_SOFA(sofa_name, squeeze(sofa_hrtf.Data.IR(selected_data, 1, :)), ...
%               squeeze(sofa_hrtf.Data.IR(selected_data, 1, :)),size(sofa_hrtf.Data.IR, 3),...
%               sofa_hrtf.Data.SamplingRate, selected_data, ...
%               [], ['good_SOFA/' sofa_name(1:end-5) '_edited.sofa']);

%
% 1. input = 'irc_1004.sofa' % original SOFA file
%
% 2. hrir_left = norm_hrir_left % left hrtf matrix
%    - number of measurement x HRIR length matrix, e.g. 3 x 1024 from the example 
%    - left only
%
% 3. hrir_right = norm_hrir_right % left hrtf matrix
%    - number of measurement x HRIR length matrix, e.g. 3 x 1024 from the example 
%    - right only
%
% 4. target_length = 1024; % new hrir length
%
% 5. target_fs = 44100; % new sampling frequency
%
% 6. hrtf_id = [1 2 3 4 5 6 7 8 9 10] % sample number vector
%    - hrtf measurement number (row number) from the original SOFA file
%    - if empty (e.g. []), equal to all measurements
% 
% 7. SOFAdir = 'normalised_SOFA/' % new file save location (OPTIONAL)
%    - in char
%    - make sure the folder exist
%    - add / in the end
%    - save in current folder if empty or leave blank
%    
% 8. SOFAfilename = []; or desire file name e.g. 'new_file.sofa'; (OPTIONAL)
%    - in char
%    - add .sofa in the end
%    - save in original input file name + '_normalised.sofa' if empty or leave blank
%    - 0 : no saving SOFA file
%
% 9. print = 1;  % print trigger
%    print out the progress 
%    - default(if empty or missing) = 1
%    
%

% OUTPUT:
% Obj: SOFA struct (same as SOFAload('new_SOFA_file.sofa')
%
%


%% catch missing input

if nargin == 5
    hrtf_id = [];
    SOFAdir = [];
    SOFAfilename = [];
    print = 1;
elseif nargin == 6
    SOFAdir = [];
    SOFAfilename = [];
    print = 1;
elseif nargin == 7
    SOFAfilename = [];
    print = 1;
elseif nargin == 8
    print = 1;
end     
% catch missing input
    
if isempty(hrtf_id)
    hrtf_id = linspace(1,size(SOFA_hrtf.Data.IR, 1),size(SOFA_hrtf.Data.IR, 1));
    % take all hrtf 
end

if ischar(SOFAdir)
    if ~strcmpi(SOFAdir(end), '/')
        SOFAdir = [SOFAdir '/'];
    end
end
% catch if / is missing in the end of SOFAdir
    
if ischar(SOFAfilename)
    if length(SOFAfilename) > 6
        if ~strcmpi(SOFAfilename(end - 4 : end), '.sofa')
            SOFAfilename = [SOFAfilename '.sofa'];
        end
    else
        SOFAfilename = [SOFAfilename '.sofa'];
    end
end
% catch if .sofa is missing in the SOFAfilename

if ~isempty(SOFAdir)
    if ~exist(SOFAdir, 'dir') % check if output folder exist
        mkdir(SOFAdir) 
        % if output folder did not exist, create output folder 
        warning(['created output folder ' SOFAdir(1:end-1) '.'])
    end
end

if isempty(print)
    print = 1;
    % defult print
end

%% Get an empty conventions structure
Obj = SOFAload(input_SOFA);

%% Fill data with data

Obj.Data.IR = zeros(size(hrir_left,1),2,size(hrir_left,2));
% data alocation
Obj.Data.IR(:,1,:) = hrir_left;
Obj.Data.IR(:,2,:) = hrir_right;
% fill in hrir left and right
Obj.Data.SamplingRate = target_fs; 
% update new sample rate

%% Update attributes

if isfield(Obj,'GLOBAL_ListenerShortName')
    Obj.GLOBAL_ListenerShortName = [Obj.GLOBAL_ListenerShortName '_normalised'];
end
if isfield(Obj,'GLOBAL_History')
    Obj.GLOBAL_History = ['normalised ' input_SOFA ' with Fs: ' num2str(target_fs) ...
    ' and length: ' num2str(target_length) '. ' Obj.GLOBAL_History];
end
if isfield(Obj,'GLOBAL_Comment')
    Obj.GLOBAL_Comment = [Obj.GLOBAL_Comment ...
        ' normalise by Benjamin Tsui ...(bt712@york.ac.uk), UoY Audio Lab'];
end
if isfield(Obj,'GLOBAL_Title')
    Obj.GLOBAL_Title = [Obj.GLOBAL_Title ' normalised with Fs: ' num2str(target_fs) ...
        ' and length: ' num2str(target_length)];
end

t.Format = 'yyyy-MM-dd';
Obj.GLOBAL_DateModified = datestr(datetime('now'));

%% update the mandatory variables

if isfield(Obj,'SourcePosition')
    Obj.SourcePosition = Obj.SourcePosition(hrtf_id, :);
end
if isfield(Obj,'MeasurementSourceAudioChannel')
    Obj.MeasurementSourceAudioChannel = Obj.MeasurementSourceAudioChannel(hrtf_id, :);
end
if isfield(Obj,'MeasurementAudioLatency')
    Obj.MeasurementAudioLatency = Obj.MeasurementAudioLatency(hrtf_id, :);
end
Obj.Data.SamplingRate = target_fs;
Obj.Data.SamplingRate_Units = 'hertz';
Obj.API.N = target_length;

%% Update dimensions
Obj = SOFAupdateDimensions(Obj);

%% save SOFA file
if ischar(SOFAfilename) || isempty(SOFAfilename)
    compression = 1;
    if isempty(SOFAfilename)
        SOFAfilename = [SOFAdir input_SOFA(1 : end-5) '_normalised.sofa' ];
    else
        SOFAfilename = [SOFAdir SOFAfilename];
    end
    
    SOFAsave(SOFAfilename, Obj, compression);
    
    if print ~= 0
        disp([' - Saved: ' SOFAfilename]);
    end
    
elseif SOFAfilename ~= 0
    message = sprintf(['SOFA file unsaved: make sure filename is char or empty []' ...
        '(output will be \n input name + _normalised) or 0 if no file need to be saved']);
    sldiagviewer.reportWarning(message);
end


end

