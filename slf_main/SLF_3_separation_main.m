function SLF_3_separation_main(fgFile,paramName,paramFile,hemi,preAfqClean,resultsDir)
% This function is used as a wrapper for slf_separate_slf3.
%
% INPUT
% fgFile: a fiber group to be read by fgRead (could be Vistasoft format \ tck
%         format \ more). If created using TractSeg, consider resampling to 1 mm
%         using fg = dtiResampleFiberGroup(fg,1,'L').
%         If fg is an array of 3 fiber groups, it is considered as SLF1,2,3.
% paramName: The parameter to be used for separating SLF3 (FA, T1, MTVm T2m
%            T2star, T2wT1w)
% paramFile: The full path of the nifti file to use (FA file, T1 file etc)
% hemi: 'left' or 'right' hemisphere
% preAfqClean: A logical flag for running standard AFQ cleaning before
%              separating SLF3 (useful for probabilistic tractography)
% resultsDir: The output directory where the results file is saved.
%
% OUTPUT
% Two files are written as output:
% * The results file, which holds many of the important algorithm outputs.
% * The SLF3 fiber group file. You may change it to other formats (like
%   tck), within slf_separate_slf3.

%% Set the sign of change in parameter profile (increase/decrease towrds SLF3)
paramImg = readFileNifti(paramFile);

switch paramName
    case 'FA'
        dirFlag = 'increase';
    case 'T1'
        dirFlag = 'decrease';
    case 'MTV'
        dirFlag = 'decrease';
    case 'T2'
        dirFlag = 'decrease';
    case 'T2star'
        dirFlag = 'decrease';
    case 'T2wT1w'
        dirFlag = 'decrease';
        paramImg.data = 1./paramImg.data;
        paramImg.data(paramImg.data == inf) = 0;
    case 'MD'
        dirFlag = 'increase';% In some datasets, it seems like 'increase' would work better. Nevertheless, using MD to separate SLF-III is not recommended (See paper)
    case 'b0'
        dirFlag = 'decrease';
    otherwise
        error('Unrecognized option for param_imgNum');
end

%% Set ouput file
fileName = ['results_',paramName,'_',hemi,'.mat'];
outFile = fullfile(resultsDir,fileName);
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

%% Run
fg = fgRead(fgFile);
slf_separate_slf3(fg,paramImg,hemi,dirFlag,outFile,preAfqClean);
