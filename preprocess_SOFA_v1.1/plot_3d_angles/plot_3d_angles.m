function plot_3d_angles(varargin)
%PLOT_3D_ANGLES Summary of this function goes here
%   Detailed explanation goes here

% extender version of "scatter3" function for hrtf measurement angels
% changed input to degree
% added head model
% fix axis

% DEMO:
% plot_3d_angles('irc_1007.sofa', 'Marker', 'x', 'MarkerEdgeColor', 'blue' )
% hold on
% plot_3d_angles('subject_012.sofa', 'Marker', 'o', 'MarkerEdgeColor', 'red' )
% hold on
% plot_3d_angles('SADIE_004_DFC_256_order_fir_48000.sofa', 'Marker', '+', 'MarkerEdgeColor', 'green' )


% INPUT (3 option):
%
% Option 1:
%   azimuth angle (or array), elevation angle (pr array), distance (or array)
%   e.g. plot_3d_angles(-45, 30, 1.5)
%        or
%        plot_3d_angles([-45 20], [30 15], [1.5 1.2])
%        or
%        plot_3d_angles([-45 20], [30 15], 1.5)
%        or
%        plot_3d_angles(-45, [30 15], [1.5 0.2])
%        or
%        plot_3d_angles([-45 20], 30 , [1.5 0.2])
%
%
% Option 2:
%   loaded sofa file in struct, row number (optional)
%   if 2nd input is empty, it will plot all angles
%   e.g. hrtf = SOFAload('irc_1007.sofa'); 
%        then
%        plot_3d_angles(hrtf) % plot all angles
%        or
%        plot_3d_angles(hrtf, 1) % plot row 1 in the sofa file
%        or
%        plot_3d_angles(hrtf, [1 35 60]) % plot row 1, 35 and 60 in the sofa file
%  
% Option 3 (similar to option 2):
%   name of the sofa file in struct, row number (optional)
%   if 2nd input is empty, it will plot all angles
%   e.g. plot_3d_angles('irc_1007.sofa') % plot all angles
%        or
%        plot_3d_angles('irc_1007.sofa', 1) % plot row 1 in the sofa file
%        or
%        plot_3d_angles('irc_1007.sofa', [1 35 60]) % plot row 1, 35 and 60 in the sofa file
%

% Set Properties:
%   different scatter properties can be set after the input data
%   tested with 'Marker', 'MarkerEdgeColor' and 'Legend'
% 

%% initialise input

if isstruct(varargin{1})  % if input 1 is struct (option 2)
    hrtf_angle = SOFAcalculateAPV(varargin{1});
    
    index = find(cellfun('isclass', varargin(2:end), 'char'));
    if nargin > 1 && ~isempty(index)
        index = index(1) + 1;
        k = length(varargin) - index;
        varargin(4:4+k) = varargin(index:end);
    end
    % move the scatter properties inputs to the end 
    % after azithum angle, elevation angle and distnace (ref. option 1)
        
    if nargin > 1 && isa(varargin{2},'double')
        data_no = varargin{2};
        varargin{1} = hrtf_angle(data_no, 1);
        varargin{2} = hrtf_angle(data_no, 2);
        varargin{3} = hrtf_angle(data_no, 3);
    else
        varargin{1} = hrtf_angle(:, 1);
        varargin{2} = hrtf_angle(:, 2);
        varargin{3} = hrtf_angle(:, 3);
    end
    % find desire ploting angles or all angles if input 2 is empty
    % basically convert option 2  to option 1
        
elseif ischar(varargin{1}) % if input 1 is string (option 3)
    hrtf = SOFAload(varargin{1});
    hrtf_angle = SOFAcalculateAPV(hrtf);
    
    index = find(cellfun('isclass', varargin(2:end), 'char'));
    if nargin > 1 && ~isempty(index)
        index = index(1) + 1;
        k = length(varargin) - index;
        varargin(4:4+k) = varargin(index:end);
    end
    % move the scatter properties inputs to the end 
    % after azithum angle, elevation angle and distnace (ref. option 1)
    
    if nargin > 1 && isa(varargin{2},'double')
        data_no = varargin{2};
        varargin{1} = hrtf_angle(data_no, 1);
        varargin{2} = hrtf_angle(data_no, 2);
        varargin{3} = hrtf_angle(data_no, 3);
    else
        varargin{1} = hrtf_angle(:, 1);
        varargin{2} = hrtf_angle(:, 2);
        varargin{3} = hrtf_angle(:, 3);
    end  
    % find desire ploting angles or all angles if input 3 is empty
    % basically convery option 3  to option 1
end

%% Plot result

z = varargin{3};

[varargin{1}, varargin{2}, varargin{3}] = sph2cart(deg2rad(varargin{1}), deg2rad(varargin{2}), varargin{3});
% convert degree (azimuth, elevation, distance) to cartesian (x, y, z) coordinates
scatter3(varargin{:});
% plot result

ax = gca;
if ax.XLim(2) < max(z)*1.2 || strcmp(ax.XLimMode, 'auto')
    window_size = [-max(z)*1.2 max(z)*1.2];
    set(ax,'XLim',window_size,'YLim',window_size,'ZLim',window_size)
end
% set plot size

hold on
[X,Y,Z] = sphere(64);
x = 0.1*X(:);
y = 0.1*Y(:);
z = 0.1*Z(:);
head = scatter3(x,y,z,'MarkerEdgeColor', [0.9 0.6 0.5], 'MarkerFaceColor' ,[0.9 0.6 0.5]);
% plot head
set(get(get(head,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
    
[x, y, z] = sph2cart(deg2rad(-25), deg2rad(10), 0.115);
left_eye = scatter3(x,y,z,'MarkerEdgeColor', 'k', 'MarkerFaceColor' ,'k');
set(get(get(left_eye,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
[x, y, z] = sph2cart(deg2rad(25), deg2rad(10), 0.115);
right_eye = scatter3(x,y,z,'MarkerEdgeColor', 'k', 'MarkerFaceColor' ,'k');
% plot eyes
set(get(get(right_eye,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend

[x, y, z] = sph2cart(0 , deg2rad(-10), linspace(0.09, 0.13, 8));
nose = scatter3(x,y,z,'MarkerEdgeColor', 'k', 'MarkerFaceColor' ,'k');
set(get(get(nose,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
% plot nose

[x, y, z] = sph2cart(deg2rad(linspace(95,100, 50)), deg2rad(linspace(-15,5, 50)), 0.12);
left_ear = scatter3(x,y,z,'MarkerEdgeColor', 'k', 'MarkerFaceColor' ,'k');
set(get(get(left_ear,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
[x, y, z] = sph2cart(deg2rad(linspace(-95,-100, 50)), deg2rad(linspace(-15,5, 50)), 0.12);
right_ear = scatter3(x,y,z,'MarkerEdgeColor', 'k', 'MarkerFaceColor' ,'k');
set(get(get(right_ear,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
% plot ears

axis square

view([0 90]) 
% change inital angle
  
hold off

dcm_obj = datacursormode;
set(dcm_obj,'UpdateFcn',@disp_ploar)
% display spherical coordinates (in degrees) instead of cartesian coordinates

xlabel('X / metre');
ylabel('Y / metre');
zlabel('Z / metre');
% add labels

%x0=300;
%y0=130;
%width=510;
%height=555;
%set(gcf,'units','points','OuterPosition',[x0,y0,width,height])
% set plot window position and size

end

