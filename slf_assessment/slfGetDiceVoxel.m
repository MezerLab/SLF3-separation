function dice = slfGetDiceVoxel(fgSlf3Param,fgSlf3TractSeg,nii,writeFdImg)
% slfGetDiceStreamline reutrns the Dice coefficient calculated based on
% voxel-wise agreement between SLF3 as identified by TractSeg (fg(3))
% and SLF3 as identified by some local parameter (e.g., FA or T1)
% nii is used as a reference for creating the fiber density image (we use a
% high resolution T1 map)
%
% Note, if using a probabilistic candidate set, consider clipping the fg's
% between the waypoint planes before giving them as input.

if ~exist('writeFdImg','var') || isempty(writeFdImg)
    writeFdImg = false;
end

normalizeMap = false;

% If each fg is actually two fg's, then interpret the 2nd one as the entire
% set, and use it to constrain the comparison only to shared voxels across
% the two initial sets (good for comparing probabilistic tractography with
% TractSeg)
if length(fgSlf3Param)==2 && length(fgSlf3TractSeg)==2
    restrictToShared = 1;
    fgSlf3ParamInitial = fgSlf3Param(2);
    fgSlf3Param = fgSlf3Param(1);
    fgSlf3TractSegInitial = fgSlf3TractSeg(2);
    fgSlf3TractSeg = fgSlf3TractSeg(1);

    fdImgTractSegInitial = dtiComputeFiberDensityNoGUI(fgSlf3TractSegInitial, nii.qto_xyz, size(nii.data), normalizeMap);
    fdImgTractSegInitial = double(logical(fdImgTractSegInitial));
    fdImgParamInitial = dtiComputeFiberDensityNoGUI(fgSlf3ParamInitial, nii.qto_xyz, size(nii.data), normalizeMap);
    fdImgParamInitial = double(logical(fdImgParamInitial));

    fdImgMask = zeros(size(fdImgParamInitial));
    fdImgMask(fdImgTractSegInitial & fdImgParamInitial) = 1;
    fdImgMask = logical(fdImgMask);
else
    restrictToShared = 0;
end

% Calculate the fiber density images
% TractSeg (this is fg without AFQ cleaning)
fdImgTractSeg = dtiComputeFiberDensityNoGUI(fgSlf3TractSeg, nii.qto_xyz, size(nii.data), normalizeMap);
if writeFdImg
    dtiWriteNiftiWrapper(fdImg,nii.qto_xyz,fdImgTractSegFile);
end
fdImgTractSeg = double(logical(fdImgTractSeg));

% T1-based SLF-III (this is the final fg, after AFQ cleaning)
fdImgParam = dtiComputeFiberDensityNoGUI(fgSlf3Param, nii.qto_xyz, size(nii.data), normalizeMap);
if writeFdImg
    dtiWriteNiftiWrapper(fdImg,nii.qto_xyz,fdImgParamFile);
end
fdImgParam = double(logical(fdImgParam));

% Restrict to shared voxels in the initial sets
if restrictToShared
    fdImgParam(~fdImgMask) = 0;
    fdImgTractSeg(~fdImgMask) = 0;
end

% Caculate Dice
overlap = length(find(fdImgTractSeg==fdImgParam & fdImgParam==1));
unite = length(find(fdImgTractSeg==1)) + length(find(fdImgParam==1));
dice = 2*overlap/unite; %not clipped!!!
