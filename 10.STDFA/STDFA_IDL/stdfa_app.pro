pro stdfa_app
  compile_opt idl2
  ;¿ªÆôENVIÊ±¶Î
  envi,/RESTORE_BASE_SAVE_FILES
  envi_batch_init
  rs_ratio    = 20
  ab_threshold= 0.05
  half_win_size=1;[1,2,3,4,5,6,7,8,9,10]
  ;2004187-by-2004123
  ;coarse image on predication date
  cpre_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004187_500m'
  ;coarse image on base date
  cbase_fp = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004123_500m'
  ;fine image on base date
  fbase_fp    = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_123_May02'
  cluster_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_123_May02_isodata_C3'
  
  ds_fp       = 'D:\IrinaEmelyanovaUnmix\LGC\stdfa\MOD_2004187_stdfa_E3W'
  for nw=0,N_ELEMENTS(half_win_size)-1 do begin
    stdfa,cbase_fp,cpre_fp,fbase_fp,cluster_fp,rs_ratio,ab_threshold,half_win_size[nw],$
      ds_fp+strtrim(string(half_win_size[nw]*2+1),2)    
  endfor
  
  ;2004235-by-2004187
  ;coarse image on predication date
  cpre_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004235_500m'
  ;coarse image on base date
  cbase_fp = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004187_500m'
  ;fine image on base date
  fbase_fp    = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_187_Jul05'
  cluster_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_187_Jul05_isodata_C3'

  ds_fp       = 'D:\IrinaEmelyanovaUnmix\LGC\stdfa\MOD_2004235_stdfa_E3W'
  stdfa,cbase_fp,cpre_fp,fbase_fp,cluster_fp,rs_ratio,ab_threshold,half_win_size,$
      ds_fp+strtrim(string(half_win_size*2+1),2)

  ;2004123-by-2004107
  ;coarse image on predication date
  cpre_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004123_500m'
  ;coarse image on base date
  cbase_fp = 'D:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004107_500m'
  ;fine image on base date
  fbase_fp    = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_107_Apr16'
  cluster_fp  = 'D:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_107_Apr16_isodata_C3'

  ds_fp       = 'D:\IrinaEmelyanovaUnmix\LGC\stdfa\MOD_2004123_stdfa_E3W'
  stdfa,cbase_fp,cpre_fp,fbase_fp,cluster_fp,rs_ratio,ab_threshold,half_win_size,$
      ds_fp+strtrim(string(half_win_size*2+1),2)
  
end