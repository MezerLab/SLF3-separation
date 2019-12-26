function [border, binEdges] = slfFindBorder(fg, minfit, maxfit, keepIdx,increseDecreaseFlag)
n = 20; % Divie the fitting line to n bins
numOfDots = ceil((maxfit-minfit)/n); % The number of fit-line points in each bin

paramMdnBins = nan(1,n);

fitloc = fgGetParams(fg,'FitLocations');
indicesTmp = 1:length(fitloc);
tooFar = indicesTmp(~ismember(indicesTmp,keepIdx));
fitloc(tooFar) = nan;

binEdges = round(linspace(minfit,maxfit-numOfDots,20));
binEdges = binEdges';
for ii = 1:length(binEdges)-1
    binEdges(ii,2) = binEdges(ii+1,1)-1;
end
binEdges(n,2) = maxfit;

for binI = 1:n
    keep=[];
    for fitI = binEdges(binI,1):binEdges(binI,2)
        keep = [keep;find(fitloc==fitI)]; % Get all streamlines in this bin
    end
    if isempty(keep) && binI~=1
        paramMdnBins(binI) = paramMdnBins(binI-1);
        continue
    end
    keep = unique(keep);
    fgTmp = fgRetainIndices(fg,keep);
    paramMedian = nanmedian(fgGetParams(fgTmp,'paramMedian_fgTrim'));
    paramMdnBins(binI) = paramMedian;
end

% Find maximal change along the fitting line
tmp = diff(paramMdnBins);
tmp(1:ceil(length(tmp)/2)) = nan; % We only want the 2nd half of the fitting line
tmp(end) = nan;

if strcmp(increseDecreaseFlag,'decrease')
    [~,idx] = min(tmp); % Greatest decrease along the axis
    border = binEdges(idx,2);
    slfIIIcandidates = find(fgGetParams(fg,'FitLocations')>=border);
    if length(slfIIIcandidates)<length(fg.fibers)*0.1
        tmp(idx) = nan;
        [~,idx] = min(tmp);% second Greatest decrease
        border = binEdges(idx,2);
    end
elseif strcmp(increseDecreaseFlag,'increase')
    [~,idx] = max(tmp); % Greatest increase along the axis for FA
    border = binEdges(idx,2);
    slfIIIcandidates = find(fgGetParams(fg,'FitLocations')>=border);
    if length(slfIIIcandidates)<length(fg.fibers)*0.1
        tmp(idx) = nan;
        [~,idx] = max(tmp);% second Greatest increase
        border = binEdges(idx,2);
    end
end

%% For sanity checks
% figure('color','w');
% scatter(fgGetParams(fg,'FitLocations'),fgGetParams(fg,'paramMedian_fgTrim'))
% hold all
% plot(border*[1 1], get(gca,'ylim'));
end

