% This is a program to test the HCM fusion algorithm
clear all; close all; clc;

cbase = 'D:\simImg\C_Time1_30m';
cpre  = 'D:\simImg\C_Time2_30m';
fbase = 'D:\simImg\F_Time1';
[samples, lines, bands, dataType, interleave] = read_envi_header(strcat(fbase,'.hdr'));  

l1 = read_envi_image(fbase);
%l2 = read_envi_image();
m1 = read_envi_image(cbase);
m2 = read_envi_image(cpre);

% for i = 1:3
%     l1(:,:,i) = l1(:,:,1);
%     l2(:,:,i) = l2(:,:,1);
%     m1(:,:,i) = m1(:,:,1);
%     m2(:,:,i) = m2(:,:,1);
% end

patchSize = 7;
useOverlap = 1;
reg_param = 1/1000;
shift = 3;
addWhiteBand = 0;

for i = 1:bands
    output_image(:,:,i) = HCM_fusion(l1(:,:,i), m2(:,:,i), m1(:,:,i), patchSize, useOverlap, shift, reg_param, addWhiteBand);
end
fpre ='D:\simImg\F_Time2_HCM';
multibandwrite(output_image,fpre,'BSQ');
% nX = HCM_fusion(l1, m2, m1, patchSize, useOverlap, shift, reg_param, addWhiteBand);
% [l2_nX, ~] = RMSE(l2,nX);
