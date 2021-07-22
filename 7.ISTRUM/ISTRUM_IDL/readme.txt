;---------------------------------------------------------------------------------
;                            ISTRUM PROGRAM
;        Developed by Jianhang Ma, email: majh@radi.ac.cn
;             Key Laboratory of Digital Earth Science, 
;	   Institute of Remote Sensing and Digital Earth, 
;	            Chinese Academy of Sciences
;            
;Please cite the reference:
;Jianhang Ma, Wenjuan Zhang, Andrea Marinoni, Lianru Gao and Bing Zhang.
;An Improved Spatial and Temporal Reflectance Unmixing Model to synthesize time series Landsat-like images.
;
;                     Copyright belongs to Jianhang Ma
;---------------------------------------------------------------------------------
Files:
1. apply_istrum.pro
   Do ISTRUM with one pair of Landsat and Coarse-resolution (L-C) images.

2. apply_istrum_multi_lc.pro
   Do ISTRUM with multiple L-C pairs.

3. istrum.pro
   Main code of ISTRUM.

4. istrum_multi_lc.pro
   Main code of ISTRUM with multiple L-C pairs.

5. cost_fun.pro
   Cost function for spatial-unmixing.

6. AbundanceCaculateModule.exe , gdal111.dll
   Files for Fully Constrained Least Square (FCLS) method.

7. SVD_Endmembers.csv , SVDEndmember.txt
   Global SVD endmembers shared by Christopher Small.

Usage:
1. The program is written in IDL; IDL and ENVI are required, IDL 7.1.2 is supported;

2. Before running the program, new an IDL project with name of ISTRUM, copy all the files to the project folder;

3. Change file paths in apply_istrum.pro or apply_istrum_multi_lc.pro, like this:

  ;coarse image on predication date
  cpre_fp  = 'F:\test_data\MODIS\MOD09GA_A2004123_500m'
  ;coarse image on base date
  cbase_fp = 'F:\test_data\MODIS\MOD09GA_A2004107_500m'
  ;fine image on base date
  fbase_fp    = 'F:\test_data\Landsat\L7_2004_107_Apr16'
  ;output image
  ds_fp       = 'F:\test_data\LGC_DS\test1\MOD_2004123_E3W'
  ;ratio of spatial-resolution,  coarse/fine
  rs_ratio    = 20
  ;half of sliding-window size
  half_win_size=1
4. Compile and execute