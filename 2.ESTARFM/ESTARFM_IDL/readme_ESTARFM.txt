;---------------------------------------------------------------------------------
;                            ESTARFM PROGRAM
;               Using two pairs of fine and coarse images
;         the program can be used for whole TM scene and VI index product
;        Developed by (1) Zhu Xiaolin,email: zhuxiaolin55@gmail.com
;             Department of Land Surveying and Geo-Informatics
;             The Hong Kong Polytechnic University
;               
;Please cite the reference: Xiaolin Zhu, Jin Chen, Feng Gao, & Jeffrey G Masek.
;An enhanced spatial and temporal adaptive reflectance fusion model for complex
;heterogeneous regions. Remote Sensing of Environment,2010,114,2610-2623
;
;                     Copyright belongs to Xiaolin Zhu
;---------------------------------------------------------------------------------

1. The program is written in IDL;

2. Before running the program, classic ENVI should be opened first because the program uses some functions of ENVI

3. Parameters for ESTARFM.pro
 w=25.0                 ;set the haif window size, if 25, the window size is 25*2+1=51 fine pixels
 num_class=4.0          ;set the estimated number of classes, please set a larger value if blending images with very few bands
 DN_min=0               ;set the range of DN value of the image,If byte, 0 and 255
 DN_max=10000.0
 background=-9999       ;the value of background and missng pixels in both MODIS and Landsat images
 patch_long=500         ;set the size of block,if process whole ETM scene, set 500
 temp_file='D:\temp'    ;set the temporary file location

(1)w=25         
Set the window size (half),  if 25, the window size is 25*2+1=51;
(2) num_class=4             
Set the estimated number of classes according to the scene; 
(3) DN_min=0  
    DN_max=255             
Set the range of DN value of the image. If it is byte, the range is from 0 to 255.
(4) background=-9999           
the value of background and missng pixels in both MODIS and Landsat images
(5)patch_long=500                 
set the size of block,which is determined by the size of image. If process the whole TM scene, 500 is recommended.(block  is to solve the problem of computer memory limit).
(6)temp_file='D:\temp'            
Set the temporary file address. Please build a folder named "temp" before run the program, and write the address of this folder in program, such as 'D:\temp'.
  
4. Input the images accoring to the name of the pop-up window

5. Output the blended result

The blended result will be saved in the folder of the original images automatically. The names are the original prediction Modis image name followed by '_ESTARFM'


6. All the temporary data in the temporary file will be cleared automatically when the process is finished.


7. The test data:
   1).all the images are preprocessed. They are ready for ESTARFM.

   2).400*400 size, DN value range is 0-10000, in 2001.

   3).only have green (band1), red(band2), and NIR bands(band3),but the code can process image with more bands