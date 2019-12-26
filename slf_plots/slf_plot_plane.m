function slf_plot_plane(resultsFile,fg,paramImg,outFileBase,plotTractSeg,cmapStr)
slf_cmap % Load the colormaps
load(resultsFile)
matplotlib_colormaps; % Load the color maps
figsz = get(0,'ScreenSize');
set(0, 'DefaultFigureVisible', 'on'); % If figures are not shown, they are not exported well
resStr = '-r100';

%% If this is a TractSeg set, merge it
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

%%
if choosingStr=='z(x) is better'
    fitLine = fitLinex;
    %cleanidx = [SLF3x_afterCleaning;SLFnotx_beforeCleaning'];
    division_by_param = ones(length(fg.fibers),1); % SLF12
    division_by_param(slf3IndicesParam) = 2; % SLF3
    if exist('borderx','var')
        border = borderx;
    end
else
    fitLine = fitLinez;
    %cleanidx = [SLF3z_afterCleaning;SLFnotz_beforeCleaning'];
    division_by_param = ones(length(fg.fibers),1);
    division_by_param(slf3IndicesParam) = 2;
    if exist('borderz','var')
        border = borderz;
    end
end

xcoords = coords(1,:);
zcoords = coords(3,:);

%%
figure('Color','k');
eval(['colorMap = ' cmapStr ';']); % plasamadata (or viridisdata for FA)
vals = fgGetParams(fg,'paramMedian_fgTrim');

colorRange = prctile(vals, [1 98]);
rgb = fgGetRgb(fg,'paramMedian_fgTrim',colorRange, colorMap);
h = scatter(xcoords,zcoords,4,rgb,'filled');
hold on
h.SizeData = 49;
h.MarkerEdgeColor = 'none';
h.MarkerFaceAlpha = 0.6;
hLine = plot(fitLine(1,:),fitLine(2,:),'k','LineWidth',6);
hLine.Color = 0.9.*[1 1 1];
if exist('border','var')
    h1 = scatter(fitLine(1,border),fitLine(2,border),20,'w','filled');
    h1.MarkerFaceColor = hLine.Color;
    h1.MarkerEdgeColor = 0.8*[1 1 1];
end
% figsz = get(0,'ScreenSize');
% figsz(4) = figsz(3)/2;
set(gca,'tickdir','out','visible','off','color','k')
axis square
set(gcf,'position',figsz);
export_fig(gcf,[outFileBase, '_paramMedianTrim_ColorRange_' num2str(colorRange(1)) '_' num2str(colorRange(2)) '.png'],'-dpng',resStr);
% export_fig(gcf,[outFileBase, '_paramMedianTrim_ColorRange_' num2str(colorRange(1)) '_' num2str(colorRange(2)) '.png'],'-nocrop','-dpng','-r300');
% NOTE: You have to get the xlim and ylim here, after exporting, otherwise
% it won't work well
x1 = get(gca,'xlim');
y1 = get(gca,'ylim');
close(gcf)


%%
% figure('Color','k');
% colorMap = [0.8 0.2 0.4;0.2 0.2 0.6];
% colorMap(2,:) = colorMap(2,:)*1.2;
% rgb = vals2colormap(division_by_param,colorMap);
% h = scatter(xcoords,zcoords,4,rgb,'filled');
% hold on
% h.SizeData = 49;
% h.MarkerEdgeColor = 'none';
% h.MarkerFaceAlpha = 0.6;
% hLine = plot(fitLine(1,:),fitLine(2,:),'k','LineWidth',6);
% hLine.Color = 0.9.*[1 1 1];
% if exist('border','var')
%     h1 = scatter(fitLine(1,border),fitLine(2,border),20,'w','filled');
%     h1.MarkerFaceColor = hLine.Color;
%     h1.MarkerEdgeColor = 0.8*[1 1 1];
% end
% figsz = get(0,'ScreenSize');
% figsz(4) = figsz(3)/2;
% set(gcf,'position',figsz); %delete(findall(gcf,'Type','light'));
% set(gca,'tickdir','out','visible','off')
% axis square
% set(gca,'xlim',x1, 'ylim' ,y1,'color','k')
% export_fig(gcf,[outFileBase, '_Final_notCleaned.png'],'-dpng','-r300');
%

%%
colorMap = [0.8 0.2 0.4;0.2 0.2 0.6];
colorMap(2,:) = colorMap(2,:)*1.2;
colorMap = [[1 1 1]*0.5;
            0.8 0.2 0.4];
colorMap = [0.2 0.4 0.6;
            0.8 0.2 0.4];

       
rgb = vals2colormap(division_by_param,colorMap);

figure('Color','k');
h = scatter(coords(1,:),coords(3,:),4,rgb,'filled');%
hold on
h.SizeData = 49;
h.MarkerEdgeColor = 'none';
h.MarkerFaceAlpha = 0.6;
hLine = plot(fitLine(1,:),fitLine(2,:),'k','LineWidth',6);
hLine.Color = 0.9.*[1 1 1];
if exist('border','var')
    h1 = scatter(fitLine(1,border),fitLine(2,border),20,'w','filled');
    h1.MarkerFaceColor = hLine.Color;
    h1.MarkerEdgeColor = 0.8*[1 1 1];
end
% figsz = get(0,'ScreenSize');
% figsz(4) = figsz(3)/2;
set(gca,'tickdir','out','visible','off','color','k')
set(gca,'xlim',x1)
set(gca,'ylim',y1)
set(gcf,'position',figsz);
axis square
get(gca,'xlim')
%%
export_fig(gcf,[outFileBase, '_Final.png'],'-dpng',resStr);
% export_fig(gcf,[outFileBase, '_Final.png'],'-nocrop','-dpng','-r300');
close(gcf)

%%
if plotTractSeg
    figure('Color','k');
    rgb = vals2colormap([ones(1,2000),2.*ones(1,2000),3.*ones(1,2000)],cmapTs); % Assume that each SLF sub-bundle is 2000 streamlines
    intvec=1:length(xcoords);
    mixvec=intvec(randperm(length(intvec)));
    h = scatter(xcoords(mixvec),zcoords(mixvec),4,rgb(mixvec,:),'filled');
    hold on
    h.SizeData = 49;
    h.MarkerEdgeColor = 'none';
    h.MarkerFaceAlpha = 0.6;
    hLine = plot(fitLine(1,:),fitLine(2,:),'k','LineWidth',6);
    hLine.Color = 0.9.*[1 1 1];
    if exist('border','var')
        h1 = scatter(fitLine(1,border),fitLine(2,border),20,'w','filled');
        h1.MarkerFaceColor = hLine.Color;
        h1.MarkerEdgeColor = 0.8*[1 1 1];
    end
    % figsz = get(0,'ScreenSize');
    % figsz(4) = figsz(3)/2;
    set(gca,'tickdir','out','visible','off')
    axis square
    set(gca,'xlim',x1, 'ylim' ,y1,'color','k')
    set(gcf,'position',figsz);
    export_fig(gcf,[outFileBase, '_TractSeg_Plane.png'],'-dpng',resStr);
    close(gcf)
end

end
