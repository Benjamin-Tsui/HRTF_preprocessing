function output_txt = disp_ploar(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
[azimuth,elevation,r] = cart2sph(pos(1), pos(2), pos(3));
output_txt = {['azi: ',num2str(rad2deg(azimuth),4)],...
    ['ele: ',num2str(rad2deg(elevation),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['r: ',num2str(r,4)];
end

end