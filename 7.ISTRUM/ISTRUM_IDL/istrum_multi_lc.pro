;Improved Spatial and Temporal Reflectance Unmixing Model (ISTRUM)
;               with multiple L-C pairs
;---------------Must input parameters------------------
;cbase_fp     : string array, input file paths for coarse-resolution image on base dates
;cpre_fp      : input file path for coarse-resolution image on prediction date
;fbase_fp     : string array, input file paths for fine-resolution image on base date
;fpre_fp      : ouput file path for synthetic fine-resolution image
;rs_ratio     : ratio of spatial-resolution,  coarse/fine
;half_win_size: half of sliding-window size
;-------------NOTE-------------
;Input images must be reflectance image. The range of DN values is 0<DN<1. Data type is float
;----------------------------------------
;           developed by Ma Jianhang, email: majh@radi.ac.cn
;           Institute of Remote Sensing and Digital Earth, Chinese Academy of Sciences
;Please cite the reference:
pro istrum_multi_LC,cbase_fp,cpre_fp,fbase_fp,fpre_fp,rs_ratio,half_win_size
  compile_opt idl2
  timeBegin=SYSTIME(/SECONDS)
  
  fcnts=N_ELEMENTS(cbase_fp)
  for ii=0,fcnts-1 do begin
    fpre_temp=file_dirname(fpre_fp)+'\temp'+strtrim(string(ii),2)
    istrum,cbase_fp[ii],cpre_fp,fbase_fp[ii],fpre_temp,rs_ratio,half_win_size
  endfor
  ;============================================================
  ;read coarse image on prediction date 
  envi_open_file,cpre_fp,R_FID=cp_fid
  if cp_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_file_query,cp_fid,dims=cp_dims,NL=cp_lines,NS=cp_samples,NB=cp_bands,$
    BNAMES=cp_bnames,DATA_TYPE=cp_dt,WL=cp_WL
  cp_Img=MAKE_ARRAY(cp_samples,cp_lines,cp_bands,TYPE=cp_dt)
  for nb=0,cp_bands-1 do cp_Img[*,*,nb]=envi_get_data(FID=cp_fid,dims=cp_dims,pos=nb)
  ENVI_FILE_MNG,ID=cp_fid,/REMOVE
  ;============================================================
  ;read coarse images on base date 
  cb_Imgs=MAKE_ARRAY(cp_samples,cp_lines,fcnts,cp_bands,TYPE=cp_dt)
  for ii=0,fcnts-1 do begin
    envi_open_file,cbase_fp[ii],R_FID=cb_fid
    if cb_fid eq -1 then begin
      envi_batch_exit
      return
    endif   
    for nb=0,cp_bands-1 do cb_Imgs[*,*,ii,nb]=envi_get_data(FID=cb_fid,dims=cp_dims,pos=nb)
    ENVI_FILE_MNG,ID=cb_fid,/REMOVE
  endfor
  ;============================================================
  ;calculate weight
  cweights=fltarr(cp_samples,cp_lines,fcnts,cp_bands)
  Dijk=fltarr(fcnts)
  for i=0,cp_lines-1 do begin
    ;Top Line
    IsTL = i - half_win_size
    ;Buttom Line
    IsBL = i + half_win_size
    ai   = max([0,IsTL])          ;line begin
    bi   = min([cp_lines-1,IsBL]) ;line end
    for j=0,cp_samples-1 do begin
      ;Left Sample
      IsLS = j - half_win_size
      ;Right Sample
      IsRS = j + half_win_size
      aj   = max([0,IsLS])             ;sample begin
      bj   = min([cp_samples-1,IsRS])  ;sample end 
      for ib=0,cp_bands-1 do begin
        for ii=0,fcnts-1 do $
          Dijk[ii]=1.0/total(abs(reform(cb_Imgs[aj:bj,ai:bi,ii,ib])-reform(cp_Img[aj:bj,ai:bi,ib])))
        
        cweights[j,i,*,ib]=Dijk/total(Dijk)
      endfor
      
    endfor
  endfor
  cb_Imgs=0b
  cp_Img=0b
  ;============================================================
  ;combine each prediction
  envi_open_file,fbase_fp[0],R_FID=fb_fid
  if fb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_file_query,fb_fid,dims=fb_dims,NL=fb_lines,NS=fb_samples,NB=fb_bands,$
    BNAMES=fb_bnames,DATA_TYPE=fb_dt,WL=fb_WL,XSTART=fb_xstart,YSTART=fb_ystart
  fb_mapinfo =ENVI_GET_MAP_INFO(FID=fb_fid)
  ENVI_FILE_MNG,ID=fb_fid,/REMOVE
  
  fp_Img=MAKE_ARRAY(fb_samples,fb_lines,fb_bands,TYPE=fb_dt)
  fpb_Img=MAKE_ARRAY(fb_samples,fb_lines,fb_bands,TYPE=fb_dt)
  for ii=0,fcnts-1 do begin
    fpb_fp=file_dirname(fpre_fp)+'\temp'+strtrim(string(ii),2)
    envi_open_file,fpb_fp,R_FID=fpb_fid
    if fpb_fid eq -1 then begin
      envi_batch_exit
      return
    endif
    for ibr=0,cp_bands-1 do fpb_Img[*,*,ibr]=envi_get_data(FID=fpb_fid,dims=fb_dims,pos=ibr)
    
    for i=0,cp_lines-1 do begin
      fai=i*rs_ratio
      fbi=(i+1)*rs_ratio-1
      for j=0,cp_samples-1 do begin
        faj=j*rs_ratio
        fbj=(j+1)*rs_ratio-1
        for ib=0,cp_bands-1 do $
          fp_Img[faj:fbj,fai:fbi,ib]=fp_Img[faj:fbj,fai:fbi,ib]+(reform(cweights[j,i,ii,ib]))[0]*fpb_Img[faj:fbj,fai:fbi,ib]
        
      endfor 
    endfor
    ENVI_FILE_MNG,ID=fpb_fid,/REMOVE
  endfor
  cweights=0b
  
  openw,lun,fpre_fp,/GET_LUN
  writeu,lun,fp_Img
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=fpre_fp,NB=fb_bands,XSTART=fb_xstart,YSTART=fb_ystart,$
    NS=fb_samples,NL=fb_lines,BNAMES=cb_bnames,WL=fb_WL,$
    INTERLEAVE=0,DATA_TYPE=fb_dt,MAP_INFO=fb_mapinfo,/WRITE
  print,'ISTRUM with multiple L-C pairs finished, time consumed'+$
    string(SYSTIME(/SECONDS)-timeBegin,FORMAT='(f18.6)')+' seconds'  
end