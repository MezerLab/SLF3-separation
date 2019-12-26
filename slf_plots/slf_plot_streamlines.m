function slf_plot_streamlines(resultsFile,fg,paramImg,outFileBase,hemi,t1wFile,plotTractSeg,cmapStr)
slf_cmap % Load the colormaps
load(resultsFile)
matplotlib_colormaps; % Load the color maps
figsz = get(0,'ScreenSize');
set(0, 'DefaultFigureVisible', 'on'); % If figures are not shown, they are not exported well

uniColor = 1; % For the graphical abstract
if uniColor==1
   rgb = lines(2);
   cmapTs = [rgb(2,:);rgb(2,:);rgb(2,:)];
end

tubesFlag = 1;
resStr = '-r100'; % Set resolution

if strcmp(hemi,'left')
    a = -71;%-64;%-90;
    b = 31;%23;%40;
    lightPos = [-25 0 10];
    coronalSlice = [-5 0 0];
elseif strcmp(hemi,'right')
    a = 71;%64;%90;
    b = 31;%23;%40;
    lightPos = [25 0 10];
    coronalSlice = [5 0 0];
end

t1w = readFileNifti(t1wFile);
[t1w.data, ~] = mrAnatHistogramClip(t1w.data);
if t1w.qto_xyz(1,1)<0
    % Since the xform here has -1 in the first row, it will cause error downstream in AFQ_AddImageTo3dPlot for plotting an axial slice
    t1w.qto_xyz(1,1) = -t1w.qto_xyz(1,1); 
    t1w.qto_ijk(1,1) = -t1w.qto_ijk(1,1);
end
%% Get xlim, ylim, zlim
if length(fg)==3
    fgTmp = fgMerge(fg(1),fg(2));
    fgTmp = fgMerge(fgTmp,fg(3));
else
    fgTmp = fg;
end
AFQ_RenderFibers(fgTmp,'color',cmapTs(3,:), 'camera', [a,b], 'tubes', 0, 'jittershading',0.2,'radius',[0.1,1], 'newfig',1);
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
zlim = get(gca,'zlim');
close(gcf)

%% TractSeg SLF3 and Entire TractSeg
if plotTractSeg
    classification = [ones(1,2000),2.*ones(1,2000),3.*ones(1,2000)]; % Assume that each SLF sub-bundle is 2000 streamlines
    AFQ_RenderFibers(fg(3),'color',cmapTs(3,:), 'camera', [a,b], 'tubes', tubesFlag, 'jittershading',0.2,'radius',[0.1,1], 'newfig',1);
    ax.Visible = 'off';
    grid off, axis off
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);    
    AFQ_AddImageTo3dPlot(t1w, coronalSlice);
    AFQ_AddImageTo3dPlot(t1w, [0,0,-10]);
    figsz = get(0,'ScreenSize');
    set(gcf,'position',figsz); delete(findall(gcf,'Type','light'));
    set(gcf,'color','k')
    h=light;
    h.Position = lightPos;
    grid off; axis off
    set(gcf,'color','k')
% % % export_fig(gcf,[outFileBase '_tractSeg_SLF3.png'],'-dpng',resStr);

    
    AFQ_RenderFibers(fg(1),'color',cmapTs(1,:), 'camera', [a,b], 'tubes', tubesFlag, 'jittershading',0.2,'radius',[0.1,1], 'newfig',0);
    AFQ_RenderFibers(fg(2),'color',cmapTs(2,:), 'camera', [a,b], 'tubes', tubesFlag, 'jittershading',0.2,'radius',[0.1,1], 'newfig',0);
    AFQ_AddImageTo3dPlot(t1w, coronalSlice);
    AFQ_AddImageTo3dPlot(t1w, [0,0,-10]);
    export_fig(gcf,[outFileBase '_tractSeg_Set.png'],'-dpng',resStr);

    if strcmp(hemi,'left')
        xCoords = min(xlim):1;
    elseif strcmp(hemi,'right')
        xCoords = 1:max(xlim);
    end
    yMost=nanmean(coords(2,:));
    cmapTmp = lines(3);
    zCoords = -10:max(zlim);
    [X,Y,Z] = meshgrid(xCoords,yMost,zCoords);
    scatter3(X(:),yMost*ones(size(Y(:))),Z(:),36,cmapTmp(3,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    export_fig(gcf,[outFileBase '_tractSeg_yellowPlane.png'],'-dpng',resStr);
    close all
end


%%
if length(fg)==3
    fgTmp = fgMerge(fg(1),fg(2));
    fg = fgMerge(fgTmp,fg(3));
end

%% Trim streamlines around the coronal plane
numNodes = 20;
fgTrim = fg;
for fI = 1:length(fgTrim.fibers)
    x = fgTrim.fibers{fI};
    yDist = abs(x(2,:)-yMost);
    idx = find(yDist == min(yDist));
    
    n1 = max(1,idx-numNodes);
    n2 = min(size(fgTrim.fibers{fI},2),idx+numNodes);
    fgTrim.fibers{fI} = fgTrim.fibers{fI}(:,n1:n2);
end

%% Get the paramMedian values
fg = fgRemoveParams(fg,'paramMedian_fgTrim');
perPointFlag = 0;
fgTrim = dtiCreateQuenchStats(fgTrim, 'paramMedian_fgTrim', 'param', perPointFlag, paramImg, 'nanmedian', 1); % this is my version of dtiCreateQuenchStats, where in line 214 a "std" options exists.
fg.params{end+1} = fgTrim.params{end};

%% Color streamlines by parameter
figure('Color','k');

eval(['colorMap = ' cmapStr ';']); % plasamadata (or viridisdata for FA)


vals = fgGetParams(fg,'paramMedian_fgTrim');

colorRange = prctile(vals, [1 98]);
rgb = fgGetRgb(fg,'paramMedian_fgTrim',colorRange, colorMap);

AFQ_RenderFibers(fg,'color',rgb, 'camera', [a,b], 'numfibers',length(fg.fibers), 'tubes', tubesFlag, 'jittershading',0,'radius',[0.1,1], 'newfig',1);
set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
AFQ_AddImageTo3dPlot(t1w, coronalSlice);
AFQ_AddImageTo3dPlot(t1w, [0,0,-10]);
ax.Visible = 'off';
grid off, axis off
figsz = get(0,'ScreenSize');
set(gcf,'position',figsz);
delete(findall(gcf,'Type','light'));
h=light;
h.Position = lightPos;
grid off
axis off
set(gcf,'color','k')
export_fig(gcf,[outFileBase '_paramMedian_range_' num2str(colorRange(1)) '_' num2str(colorRange(2)) '_Set.png'],'-dpng',resStr);

%% Final SLF3
rgb = [0.8 0.2 0.4];
fgSlf3 = fgRetainIndices(fg,slf3IndicesParam);
AFQ_RenderFibers(fgSlf3,'color',rgb, 'camera', [a,b], 'tubes', tubesFlag, 'jittershading',0.2,'radius',[0.1,1], 'newfig',1);
set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
AFQ_AddImageTo3dPlot(t1w, coronalSlice);
AFQ_AddImageTo3dPlot(t1w, [0,0,-10]);
ax.Visible = 'off';
grid off, axis off
figsz = get(0,'ScreenSize');
set(gcf,'position',figsz); delete(findall(gcf,'Type','light'));
set(gcf,'color','k')
h=light;
h.Position = lightPos;
grid off; axis off
set(gcf,'color','k')
export_fig(gcf,[outFileBase '_Final_SLF3.png'],'-dpng',resStr);
close all
end
