pro apply_istrum_multi_lc
  compile_opt idl2
  ;start ENVI
  envi,/RESTORE_BASE_SAVE_FILES
  envi_batch_init
  
  ;2004187-by-2004107-2004235
  ;coarse image on predication date
  cpre_fp   = 'F:\test_data\MODIS\MOD09GA_A2004187_500m'
  
  ;coarse image on base date 1
  cbase_fp1 = 'F:\test_data\MODIS\MOD09GA_A2004107_500m'
  ;fine image on base date 1
  fbase_fp1 = 'F:\test_data\Landsat\L7_2004_107_Apr16'
  ;coarse image on base date 2
  cbase_fp2 = 'F:\test_data\MODIS\MOD09GA_A2004235_500m'
  ;fine image on base date 2
  fbase_fp2 = 'F:\test_data\Landsat\L7_2004_235_Aug22'
  
  ds_fp     = 'F:\test_data\LGC_DS\test2\MOD_2004187_by_107_235_W3'
  
  cbase_fp     =[cbase_fp1,cbase_fp2]
  fbase_fp     =[fbase_fp1,fbase_fp2]
  rs_ratio     =20
  half_win_size=1
  
  istrum_multi_LC,cbase_fp,cpre_fp,fbase_fp,ds_fp,rs_ratio,half_win_size
  envi_batch_exit
end