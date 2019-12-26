function rgb = fgGetRgb(fg,colorBy,colorRange,colormap)
% colorBy: 'T1_std', or 'T1_median'
% colorRange: 1X2 vector

if isempty(colorBy)
%     rgb = [0 0 1];
    rgb = lines(4);
    rgb = rgb(4,:);
    return
end

if isempty(colormap)
    colormap = 'jet';
end

if strcmpi(colorBy,'orientationRgb')
    rgb = fgGetOrientationRgb(fg);
    
    % Add a brightness jitter, to allow for better 3D visualization
    for ii = 1:size(rgb,1)
        fact = 1;%;0.9+0.1*abs(rand);
        rgb{ii} = rgb{ii}.*fact;
    end
else
    
    if isempty(colorRange)
        warning('colorRange empty, using 2nd and 98th percentile values of input data'); 
        vals = fgGetParams(fg,colorBy);
        colorRange = prctile(vals, [1 98]);
    end
    vals = fgGetParams(fg,colorBy);
    vals(vals<colorRange(1)) = colorRange(1);
    vals(vals>colorRange(2)) = colorRange(2);
    vals(end+1:end+2) = colorRange;
    rgb = vals2colormap(vals,colormap);
    rgb(end-1:end,:) = [];
end