%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program name: HCM fusion
% Version 1:
% Signal Processing, Inc
% August 1, 2018
%
% Reference: Please cite the following two papers if one uses this
% package. Also, see those two paper for details
% 1. C. Kwan, B. Budavari, F. Gao, and X. Zhu, “A Hybrid Color
% Mapping Approach to Fusing MODIS and Landsat Images for Forward
% Prediction,” Remote Sensing, vol. 10, no. 4, March 26, 2018.
% DOI: 10.3390/rs10040520.
% 2. C. Kwan, X. Zhu, F. Gao, B. Chou, D. Perez, J. Li, Y. Shen,
% K. Koperski, and G. Marchisio, “Assessment of Spatiotemporal
% Fusion Algorithms for Worldview and Planet Images,” Sensors,
% 18(4), 1051, March 31, 2018. DOI:10.3390/s18041051.
%
% Purpose: The purpose of this package is to fuse two low
% resolution images collected at t1 and t2, and one high
% resolution image collected at t1 and generate a high resolution
% image at t2.
%
% Usage: There are several parameters for performance tuning:
% patchsize: We divide the image into patches (subimages) and then
% apply HCM to each patch.
% useOverlap: 0 means no overlap between the patches; 1 means
% overlap between the patches
% shift: This is used when useOverlap is 1. shift can be any
% integers smaller than the patchsize
% reg_param: It is a regularization parameter to avoid numerical
% problems.
%
% Inputs: Two LR images at two different times; one HR image at an
% earlier time. The images should be in double precision.
%
% Outputs: A HR image at the prediction time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nX = HCM_fusion(HR_image, LR_predictionDate, LR_earlierDate, patchSize, useOverlap, shift, reg_param, addWhiteBand)

if addWhiteBand
    HR_image(:,:,size(HR_image,3)+1) = ones(size(HR_image,1), size(HR_image,2));
    LR_earlierDate(:,:,size(LR_earlierDate,3)+1) = ones(size(LR_earlierDate,1), size(LR_earlierDate,2));
    LR_predictionDate(:,:,size(LR_predictionDate,3)+1) = ones(size(LR_predictionDate,1), size(LR_predictionDate,2));
end
nX = cml_fusion(LR_predictionDate, LR_earlierDate, HR_image, patchSize, useOverlap, shift, reg_param);

end