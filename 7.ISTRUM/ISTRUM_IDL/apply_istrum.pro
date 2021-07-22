pro apply_istrum
  compile_opt idl2
  ;start ENVI
  envi,/RESTORE_BASE_SAVE_FILES
  envi_batch_init

  ;2004123-by-2004107
  ;coarse image on predication date
  cpre_fp  = 'F:\test_data\MODIS\MOD09GA_A2004123_500m'
  ;coarse image on base date
  cbase_fp = 'F:\test_data\MODIS\MOD09GA_A2004107_500m'
  ;fine image on base date
  fbase_fp    = 'F:\test_data\Landsat\L7_2004_107_Apr16'
  ds_fp       = 'F:\test_data\LGC_DS\test1\MOD_2004123_E3W'
  rs_ratio    = 20

  half_win_size=1;[1,2,3,4,5,6,7,8,9,10]

  istrum,cbase_fp,cpre_fp,fbase_fp,ds_fp+strtrim(string(half_win_size*2+1),2),$
    rs_ratio,half_win_size
  envi_batch_exit
end