function [SLF3,SLF12] = slfCleanByOutlierParam(fg,border,~,paramBundle,paramBundleStd,slf3candidates,slf12candidates,increseDecreaseFlag)

tmpMedianParam=fgGetParams(fg,'paramMedian_fgTrim');
fitLoc=fgGetParams(fg,'FitLocations');

% addedCandidates=find((fitLoc<=border) & (fitLoc > (border-binEdges(1,2)+binEdges(1,1))) & (tmpMedianT1<=(T1Bundle+T1bundleStd))');
% SLF3Final=[slfIIIcandidates', addedCandidates'];
% notSLF3(ismember(notSLF3,addedCandidates))=[];

if strcmp(increseDecreaseFlag,'decrease')
    removeCandidates=find((fitLoc>=border) & (tmpMedianParam > paramBundle+2*paramBundleStd)');
elseif strcmp(increseDecreaseFlag,'increase')
    removeCandidates=find((fitLoc>=border) & (tmpMedianParam < paramBundle-2*paramBundleStd)');
end
slf3candidates(ismember(slf3candidates,removeCandidates))=[];
SLF3 = slf3candidates; % Return the indices of SLF sub-bundles in term of the original candidate set
SLF12 = [slf12candidates' ,removeCandidates'];


end
