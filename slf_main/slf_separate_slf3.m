function slf_separate_slf3(fg,param_img,hemi,dirFlag,outFile,preAfqClean)
%%The main function for extracting slf3 from the rest of the SLF set.
% fg - the candidate set to separate. Can be an array of 3 fg (e.g., SLF-I, SLF-II & SLF-III from TractSeg)
% param_img - struct of the parameter's Nifti file
% hemi - 'left' or 'right'
% dirFlag - 'increase' or 'decrease' string, depending if we are looking for the biggest increase or decrease in the parameter's
% outFile - The output file for saving results
% preAfqClean - if true, the candidate set will be cleaned with AFQ standard cleaning (useful for probabilistic tractography candidates)

%% (1) Merge candidate streamlines
% If 3 SLF sub-bundles are given as inoput, merge them and keep track of who is who
if length(fg)==3
    origClass=[ones(length(fg(1).fibers),1);2*ones(length(fg(2).fibers),1);3*ones(length(fg(3).fibers),1)];
    fgtmp = fgMerge(fg(1),fg(2),'fgtmp');
    fg = fgMerge(fgtmp,fg(3),'fg');
end

%% (2) Standard AFQ cleaning
if ~exist('preAfqClean','var') || isempty(preAfqClean)
    preAfqClean = false;
end
    
if preAfqClean
    [~,keepPreClean] = AFQ_removeFiberOutliers(fg,4,5,25); %takes out outliers
    warning(['Removing ' num2str((length(fg.fibers)-length(find(keepPreClean)))) '/' num2str(length(fg.fibers)) ' streamlines from candidate set.']);
    fg = fgRetainIndices(fg,keepPreClean);
end

%% (3) Represent the candidates in a coronal plane
[yMost,coords] = slfFindPlane(fg);

%% (4) Add param-median from the trimmed version of the fg
% Trim streamlines around the coronal plane
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

fg = fgRemoveParams(fg,'paramMedian_fgTrim');
perPointFlag = 0;
fgTrim = dtiCreateQuenchStats(fgTrim, 'paramMedian_fgTrim', 'param', perPointFlag, param_img, 'nanmedian', 1); % this is my version of dtiCreateQuenchStats, where in line 214 a "std" options exists.
fg.params{end+1} = fgTrim.params{end};

%% (5) Represent the candidates as a curve in the coronal plane
% fit in middle plane for both fits
fitType='poly5';
[fitByX,fitByZ] = slfFindBothFits(coords,fitType);

%for z(x)
fitLinex = slfCreateFitline(fitByX,1,coords,hemi);

%for x(z)
fitLinez = slfCreateFitline(fitByZ,2,coords,hemi);

%% (6) Algorithm
radius = 3; % In mm. Only streamlines whose coordinate is close to the fitted curve (<radius) will be used for finding the separating border

%for z(x)
[fgx,minfitx,maxfitx,keepIdxx] = slfReturnFitParams(fg,coords,fitLinex,radius);
[borderx, binEdgesx] = slfFindBorder(fgx,minfitx,maxfitx,keepIdxx,dirFlag);

%for x(z)
[fgz,minfitz,maxfitz,keepIdxz] = slfReturnFitParams(fg,coords,fitLinez,radius);
[borderz, binEdgesz] = slfFindBorder(fgz,minfitz,maxfitz,keepIdxz,dirFlag);
disp("Found borders for z(x) and for x(z)")

%% (7) Some calculations for the z(x) fit
%% Clean outlier streamlines (in terms of their parameter, e.g., T1)
slf12candidatesx = find(fgGetParams(fgx,'FitLocations')<borderx);
slf3candidatesx = find(fgGetParams(fgx,'FitLocations')>=borderx);

% Add to the SLF12 group also streamlines whose FitLocations are set to NaN in slfReturnFitParams.m, due to their coords being NaN above
slf12candidatesx = unique([slf12candidatesx; find(isnan(fgGetParams(fgx,'FitLocations')))]);

slf3xTmp = fgRetainIndices(fg,slf3candidatesx);
paramBundle = median(fgGetParams(slf3xTmp,'paramMedian_fgTrim'));
stdBundle = std(fgGetParams(slf3xTmp,'paramMedian_fgTrim'));

% Remove outlier streamlines in terms of their T1 (or other parameter)
[slf3x_preAfqClean,slf12x_preAfqClean] = slfCleanByOutlierParam(fgx,borderx,binEdgesx,paramBundle,stdBundle,slf3candidatesx,slf12candidatesx,dirFlag);

%% Standard AFQ cleaning for each bundle
fgSlf3x = fgRetainIndices(fg,slf3x_preAfqClean);
[~,keep] = AFQ_removeFiberOutliers(fgSlf3x,4,5,25);
slf3x_postAfqClean = slf3x_preAfqClean(keep); % The indices (in terms of the entire SLF system) for calculating Dice later on
fgSlf3x = fgRetainIndices(fgSlf3x,keep);

slf12x = fgRetainIndices(fg,slf12x_preAfqClean);
[~,keep] = AFQ_removeFiberOutliers(slf12x,4,5,25);
slf12x_postAfqClean = slf12x_preAfqClean(keep);
slf12x = fgRetainIndices(slf12x,keep);

%% Calculate variance of param for each sub-bundle
% We will later select the fitting line that gave lower variances (assuming
% that true sub-bundles should have low-variance signatures)
paramMdn_x_slf3 = fgGetParams(fgSlf3x,'paramMedian_fgTrim');
paramMdn_x_slf12 = fgGetParams(slf12x,'paramMedian_fgTrim');
paramVar_x_slf3 = var(paramMdn_x_slf3);
paramVar_x_slf12 = var(paramMdn_x_slf12);
mean_var_x = (paramVar_x_slf3+paramVar_x_slf12)/2;

%% (8) Some calculations for the x(z) fit
%% Clean outlier streamlines (in terms of their parameter, e.g., T1)
slf12candidatesz = find(fgGetParams(fgz,'FitLocations')<borderz);
slf3candidatesz = find(fgGetParams(fgz,'FitLocations')>=borderz);

% Add to the SLF12 group also streamlines whose FitLocations are set to NaN in slfReturnFitParams.m, due to their coords being NaN above
slf12candidatesz = unique([slf12candidatesz; find(isnan(fgGetParams(fgz,'FitLocations')))]);

slf3zTmp = fgRetainIndices(fg,slf3candidatesz);
paramBundle=median(fgGetParams(slf3zTmp,'paramMedian_fgTrim'));
stdBundle=std(fgGetParams(slf3zTmp,'paramMedian_fgTrim'));

% Remove outlier streamlines in terms of their T1 (or other parameter)
[slf3z_preAfqClean,slf12z_preAfqClean] = slfCleanByOutlierParam(fgz,borderz,binEdgesz,paramBundle,stdBundle,slf3candidatesz,slf12candidatesz,dirFlag);

%% Standard AFQ cleaning, each bundle
slf3z = fgRetainIndices(fg,slf3z_preAfqClean);
[~,keep] = AFQ_removeFiberOutliers(slf3z,4,5,25);
slf3z_postAfqClean = slf3z_preAfqClean(keep); % The indices (in terms of the entire SLF system) for calculating Dice later on
slf3z = fgRetainIndices(slf3z,keep);

slf12z = fgRetainIndices(fg,slf12z_preAfqClean);
[~,keep] = AFQ_removeFiberOutliers(slf12z,4,5,25);
slf12z_postAfqClean = slf12z_preAfqClean(keep); % The indices (in terms of the entire SLF system) for calculating Dice later on
slf12z = fgRetainIndices(slf12z,keep);

%% Calculate variance of param for each sub-bundle
paramMdn_z_slf3 = fgGetParams(slf3z,'paramMedian_fgTrim');
paramMdn_z_slf12 = fgGetParams(slf12z,'paramMedian_fgTrim');
paramVar_z_slf3 = var(paramMdn_z_slf3);
paramVar_z_slf12 = var(paramMdn_z_slf12);
mean_var_z = (paramVar_z_slf3+paramVar_z_slf12)/2;
disp("final bundles for x(z) fit")

if mean_var_x<mean_var_z
    choosingStr = 'z(x) is better';
    slf3fg = fgSlf3x;
    slf3IndicesParam = slf3x_postAfqClean; % Indices of SLF3 (identified by the param) in the original candidate set for Dice calculation
    slf3IndicesParam_preAfqClean = slf3x_preAfqClean; % Indices of SLF3 (identified by the param, but before AFQ cleaning) in the original candidate set
    slf12fg = slf12x;
else
    choosingStr = 'x(z) is better';
    slf3fg = slf3z;
    slf3IndicesParam = slf3z_postAfqClean;
    slf3IndicesParam_preAfqClean = slf3z_preAfqClean;
    slf12fg = slf12z;
end

%% Save a results file
% Notice that here we also save the "keepPreClean" variable, so we can
% later identify the relevant streamlines from the original fgMori
if exist('keepPreClean','var')
    save(outFile,'yMost','coords','slf3IndicesParam','slf3IndicesParam_preAfqClean','choosingStr','fitLinex','fitLinez','minfitx','maxfitx','minfitz','maxfitz','borderx','borderz','keepPreClean'); % These variables are used for calculating Dice and for various plots
else
    save(outFile,'yMost','coords','slf3IndicesParam','slf3IndicesParam_preAfqClean','choosingStr','fitLinex','fitLinez','minfitx','maxfitx','minfitz','maxfitz','borderx','borderz'); % These variables are used for calculating Dice and for various plots
end

%% Save SLF3 as a fiber group file
[pathstr, name] = fileparts(outFile);

fileName = ['results_',paramName,'_',hemi,'.mat'];

slf3File = fullfile(pathstr, ['slf3_' name]);
fgWrite(slf3fg,slf3File);

end