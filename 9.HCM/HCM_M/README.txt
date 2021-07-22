Program name: 
	Hybrid Color Mapping (HCM) fusion
	Version 1:
	Signal Processing, Inc
	August 1, 2018
	
Reference: 
	Please cite the following two papers if one uses this 
	package. Also, see those two paper for details
	1. C. Kwan, B. Budavari, F. Gao, and X. Zhu, “A Hybrid Color 
	Mapping Approach to Fusing MODIS and Landsat Images for Forward 
	Prediction,” Remote Sensing, vol. 10, no. 4, March 26, 2018. 
	DOI: 10.3390/rs10040520.
	2. C. Kwan, X. Zhu, F. Gao, B. Chou, D. Perez, J. Li, Y. Shen, 
	K. Koperski, and G. Marchisio, “Assessment of Spatiotemporal 
	Fusion Algorithms for Worldview and Planet Images,” Sensors, 
	18(4), 1051, March 31, 2018. DOI:10.3390/s18041051.
	
Purpose: 
	The purpose of this package is to fuse two low 
	resolution images collected at t1 and t2, and one high 
	resolution image collected at t1 and generate a high resolution 
	image at t2.
	
Usage: 
	There are several parameters for performance tuning:
		patchsize: We divide the image into patches (subimages) and then apply HCM to each patch.
		useOverlap: 0 means no overlap between the patches; 1 means overlap between the patches
		shift: This is used when useOverlap is 1. shift can be any integers smaller than the patchsize
		reg_param: It is a regularization parameter to avoid numerical problems.
		addWhiteBand: input 0 to exclude adding a white band; input 1 to add a white band

Inputs: 
		lr1: LR image at earlier date
		lr2: LR image at later date (prediction date)
		hr1: HR image at earlier date
		patchsize: We divide the image into patches (subimages) and then apply HCM to each patch.
		useOverlap: 0 means no overlap between the patches; 1 means overlap between the patches
		shift: This is used when useOverlap is 1. shift can be any integers smaller than the patchsize
		reg_param: It is a regularization parameter to avoid numerical problems.
		addWhiteBand: input 0 to exclude adding a white band; input 1 to add a white band

Outputs: 
	A HR image at the prediction time.
	
Example:
	-See test_HCM_fusion.m for example usage
	
Notes:
	-To perform image prediction band by band, use a for loop to call HCM_fusion function for each band
		-example:
			for i = 1:num_of_bands
				output_image(:,:,i) = HCM_fusion(l1(:,:,i), m2(:,:,i), m1(:,:,i), patchSize, useOverlap, shift, reg_param, addWhiteBand);
			end
	
	-To perform image prediction on a full color image, only one function call is necessary
		-example
			output_image_full_color = HCM_fusion(l1, m2, m1, patchSize, useOverlap, shift, reg_param, addWhiteBand);
			
	-The test_HCM_fusion.m script can be run as is to generate a predicted image using the example images in the directory.