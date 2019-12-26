%% Separate SLF3 from TractSeg's SLF complex
exampleDataDir = '/ems/elsc-labs/mezer-a/code/roey/slf_github/ExampleData/WH008';
fgFile = fullfile(exampleDataDir,'fg_SLF_complex_left.mat');
paramName = 'FA';
paramFile = fullfile(exampleDataDir,'FA.nii.gz');
hemi = 'left';
preAfqClean = false;
resultsDir = fullfile(exampleDataDir,'results');

SLF_3_separation_main(fgFile,paramName,paramFile,hemi,preAfqClean,resultsDir)

%% Plot the SLF complex and the resulting SLF3
fgComplex = fgRead(fgFile);
fgSlf3 = fullfile(resultsDir,['slf3_results_',paramName,'_',hemi,'.mat']);

fileName = ['results_',paramName,'_',hemi,'.mat'];
resultsFile = fullfile(resultsDir,fileName);

fg = fgRead(fgFile);
paramImg = readFileNifti(paramFile);
outFileBase = fullfile(resultsDir,'plots',[paramName,'_',hemi]);
if ~exist(fullfile(resultsDir,'plots'),'dir')
    mkdir(fullfile(resultsDir,'plots'))
end
t1wFile = fullfile(exampleDataDir,'T1w_2DTI_brainMasked.nii.gz');
plotTractSeg = true;
cmapStr = 'viridisdata';
slf_plot_streamlines(resultsFile,fg,paramImg,outFileBase,hemi,t1wFile,plotTractSeg,cmapStr)

%% Plot the coronal separation plane
slf_plot_plane(resultsFile,fg,paramImg,outFileBase,plotTractSeg,cmapStr)