function [fg,minfit,maxfit,keepIdx] = slfReturnFitParams(fg,coords,fitLine,radius)
% slfFitParams finds for ecery streamline its closest point on the fitting
% line. It also returns minfit and maxfit, which mark the first and last
% fit coordinates that have streamlines assigned to them.

Xvec = coords(1,:);
Zvec = coords(3,:);
keep = ~isnan(Zvec) & ~isnan(Xvec);

minDistance = zeros(length(fg.fibers),1);
fitLocation = zeros(length(fg.fibers),1);
fitLocationClose = nan(length(fg.fibers),1); % This will save the locations along the fit of all streamlines closer than "radius" to the fitting line
for fI=1:length(fg.fibers)
    d=[];
    if ismember(fI,find(~keep))
        minDistance(fI)=nan;
        fitLocation(fI)=nan;
        continue
    end
    for ll=1:length(fitLine)
        d(ll) = sqrt((coords(1,fI)-fitLine(1,ll))^2+(coords(3,fI)-fitLine(2,ll))^2); % Get distance between every fiber (its representative coordinate) and the fitting line. Similar to "sum((coords([1,3],fI)-fitLine).^2);"
    end
    
    if min(d)<=radius
        fitLocationClose(fI)=find(d==min(d),1);
    end
    
    minDistance(fI)=min(d);
    fitLocation(fI)=find(d==min(d),1);
end

keep=~isnan(fitLocationClose);
keepIdx = 1:length(Xvec);
keepIdx = keepIdx(keep);

fitLocationClose=sort(unique(fitLocationClose));
fitLocationClose(isnan(fitLocationClose))=[];
minfit = min(fitLocationClose); % The minimal location along the fitting line with close coordinates
maxfit = max(fitLocationClose); % The maximal location along the fitting line with close coordinates

% if fitFlag==2
%     fitLocation=max(fitLocation)-fitLocation;%flip to match the slfSeparateByT1_OVERLAPPINGBINS function
% end

% Add a params fields
fg.params{end+1}.name = 'FitLocations';
fg.params{end}.uid = 0;
fg.params{end}.ile = 1;
fg.params{end}.icpp = 0;
fg.params{end}.ivs = 1;
fg.params{end}.agg = 'FitLocations';
fg.params{end}.lname = 'FitLocations';
fg.params{end}.stat = fitLocation;
fg.params{end+1}.name = 'FitDist';
fg.params{end}.uid = 0;
fg.params{end}.ile = 1;
fg.params{end}.icpp = 0;
fg.params{end}.ivs = 1;
fg.params{end}.agg = 'FitDist';
fg.params{end}.lname = 'FitDist';
fg.params{end}.stat = minDistance;
end