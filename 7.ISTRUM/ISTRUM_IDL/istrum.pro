;Improved Spatial and Temporal Reflectance Unmixing Model (ISTRUM)
;---------------Must input parameters------------------
;cbase_fp     : input file path for coarse-resolution image on base date
;cpre_fp      : input file path for coarse-resolution image on prediction date
;fbase_fp     : input file path for fine-resolution image on base date
;fpre_fp      : ouput file path for synthetic fine-resolution image
;rs_ratio     : ratio of spatial-resolution,  coarse/fine
;half_win_size: half of sliding-window size

;-------------NOTE-------------
;Input images must be reflectance image. The range of DN values is 0<DN<1. Data type is float
;coarse-resolution image in its original resolution, e.g., resolution of MODIS image should be 500m
;---------------Optional Keywords------------------
;EMSpec_fp    : Spectra of endmembers applied in fine resolution image sepctral unmixing
;               If not provide, "SVD_Endmembers.csv" will be used
;               if provide, format should be same with "SVD_Endmembers.csv"
;A_threshold  : For a coarse pixel, if abundance of an endmember is little than the threshold,
;               it will be merged to spectrally similar endmembers. Default value is 0.05
;----------------------------------------
;           developed by Ma Jianhang, email: majh@radi.ac.cn
;           Institute of Remote Sensing and Digital Earth, Chinese Academy of Sciences
;Please cite the reference:

pro istrum,cbase_fp,cpre_fp,fbase_fp,fpre_fp,rs_ratio,half_win_size,EMSpec_fp=EMSpec_fp,A_threshold=A_threshold
  compile_opt idl2
  timeBegin=SYSTIME(/SECONDS)
  ;------------------ open images -------------------
  envi_open_file,cbase_fp,R_FID=cb_fid
  if cb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_open_file,cpre_fp,R_FID=cp_fid
  if cp_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  envi_open_file,fbase_fp,R_FID=fb_fid
  if fb_fid eq -1 then begin
    envi_batch_exit
    return
  endif

  ;query the envi files
  envi_file_query,cb_fid,dims=cb_dims,NL=cb_lines,NS=cb_samples,NB=cb_bands,$
    BNAMES=cb_bnames,DATA_TYPE=cb_dt,WL=cb_WL,XSTART=cb_xstart,YSTART=cb_ystart

  envi_file_query,cp_fid,dims=cp_dims,NL=cp_lines,NS=cp_samples,NB=cp_bands,$
    BNAMES=cp_bnames,DATA_TYPE=cp_dt,WL=cp_WL,XSTART=cp_xstart,YSTART=cp_ystart

  envi_file_query,fb_fid,dims=fb_dims,NL=fb_lines,NS=fb_samples,NB=fb_bands,$
    BNAMES=fb_bnames,DATA_TYPE=fb_dt,WL=fb_WL,XSTART=fb_xstart,YSTART=fb_ystart

  cb_mapinfo =ENVI_GET_MAP_INFO(FID=cb_fid)
  fb_mapinfo =ENVI_GET_MAP_INFO(FID=fb_fid)
  
  ;------------- open endmember spectra file --------------
  routine_dir=file_dirname(ROUTINE_FILEPATH('istrum'))+'\'
  ;open endmember spectra file
  if KEYWORD_SET(EMSpec_fp) then begin
    ;the user provides endmember
    EM_fp=EMSpec_fp
  endif else begin
    ;default SVD endmember
    EM_fp=routine_dir+'SVD_Endmembers.csv'
  endelse
  EM_Spec=read_csv(EM_fp,header=EM_Name)
  em_samples=n_elements(EM_Name)
  em_lines  =n_elements(EM_Spec.(0))
  temp=fltarr(em_samples,em_lines)
  for i=0,em_samples-1 do temp[i,*]=float(EM_Spec.(i))
  EM_Spec=temporary(temp)
  
  ;---------------match the dimensions of the images---------------------
  if (cb_lines ne cp_lines) or (cb_samples ne cp_samples) or (cb_bands ne cp_bands) then begin
    print,'coarse image does not match each other'
    envi_batch_exit
    return
  endif
  if (cb_lines*rs_ratio ne fb_lines) or (cb_samples*rs_ratio ne fb_samples) or (cb_bands ne fb_bands) then begin
    print,'coarse image does not match fine image'
    envi_batch_exit
    return
  endif
  if (em_lines ne fb_bands) then begin
    print,'bands of endmembers do not match bands of fine image'
    envi_batch_exit
    return
  endif
  if (em_samples ge fb_bands) then begin
    print,'Too much endmembers, number of endmembers shold le nb-1'
    envi_batch_exit
    return
  endif
  ;number of pixels in silding-window
  win_size =2*half_win_size+1
  pixelCnts=win_size*win_size
  if pixelCnts le (em_samples) then begin
    print,'moving window is too small'
    envi_batch_exit
    return
  endif
  ;----------------- FCLS spectral-unmixing -----------------
  cd,routine_dir
  fbaseabd_fp=file_dirname(fpre_fp)+'\'+file_basename(fbase_fp)+'_abundance.tif'
  ;print,fbaseabd_fp
  cmdStr='AbundanceCaculateModule.exe '+fbase_fp+' '+EM_fp+' '+fbaseabd_fp
  SPAWN,cmdStr,/hide
  ;open the abundance image
  envi_open_file,fbaseabd_fp,R_FID=fabd_fid
  if fabd_fid eq -1 then begin
    envi_batch_exit
    print,'spectral unmixing error!'
    return
  endif
  print,'Spectral unmixing finished, time from start: '+$
    string(SYSTIME(/SECONDS)-timeBegin,FORMAT='(f13.6)')+' seconds'
  ;----------------- Abundance aggregation ------------------
  ;***** spectral angle (SA) between Endmembers *****
  lut=uintarr(em_samples,em_samples)
  for i=0,em_samples-1 do begin
    cosSAM=fltarr(em_samples)
    for j=0,em_samples-1 do begin
      if j eq i then continue
      cosSAM[j]=total(EM_Spec[i,*]*EM_Spec[j,*])/(SQRT(TOTAL(EM_Spec[i,*]^2))*SQRT(TOTAL(EM_Spec[j,*]^2)))
    endfor
    lut[*,i]=reverse(sort(cosSAM))
  endfor
  ;****************  coarse abundance     *************
  if KEYWORD_SET(A_threshold) then begin
    ab_threshold=A_threshold
  endif else begin
    ab_threshold=0.05
  endelse
  fabd_img=fltarr(fb_samples,fb_lines,em_samples)
  for ib=0,em_samples-1 do fabd_img[*,*,ib]=envi_get_data(FID=fabd_fid,pos=ib,DIMS=fb_dims)
  ENVI_FILE_MNG,ID=fabd_fid,/REMOVE,/DELETE
  cabd_img=fltarr(cb_samples,cb_lines,em_samples)   ;coarse abundance image
  nPixels=float(rs_ratio)*rs_ratio
  for cabd_i=0,cb_lines-1 do begin
    for cabd_j=0,cb_samples-1 do begin
      win_data=fabd_img[(cabd_j*rs_ratio):((cabd_j+1)*rs_ratio-1),(cabd_i*rs_ratio):((cabd_i+1)*rs_ratio-1),*]
      cabd=fltarr(em_samples)
      for ib=0,em_samples-1 do cabd[ib]=total(win_data[*,*,ib])/nPixels
      ;merge endmembers to their similar types if the abundance is lt ab_threshold
      index=where((cabd gt 0),cnts)           ;endmembers with abundance gt 0
      abd_min=min(cabd[index],min_ind)        
      while abd_min lt ab_threshold do begin
        min_ind=index[min_ind]
        abdtemp=cabd[min_ind]
        cabd[min_ind]=0.0
        for sc=0,em_samples-1 do begin
          if cabd[lut[sc,min_ind]] ne 0 then begin
            cabd[lut[sc,min_ind]]+=abdtemp
            break
          endif
        endfor
        index=where((cabd gt 0),cnts)
        abd_min=min(cabd[index],min_ind)
      endwhile
      cabd_img[cabd_j,cabd_i,*]=cabd
    endfor
  endfor
  print,'Abundance aggregation finished, time from start: '+$
    string(SYSTIME(/SECONDS)-timeBegin,FORMAT='(f13.6)')+' seconds'
  
  ;------------------ spatial-unmixing --------------------
  Change_F =MAKE_ARRAY(fb_samples,fb_lines,cb_bands,TYPE=cb_dt)  ;downscaled temporal change image
  cb_Img   =MAKE_ARRAY(cb_samples,cb_lines,cb_bands,TYPE=cb_dt)
  cp_Img   =MAKE_ARRAY(cp_samples,cp_lines,cp_bands,TYPE=cp_dt)
  fb_Img   =MAKE_ARRAY(fb_samples,fb_lines,fb_bands,TYPE=fb_dt)  
  for nb=0,cb_bands-1 do cb_Img[*,*,nb]=envi_get_data(FID=cb_fid,dims=cb_dims,pos=nb)
  for nb=0,cp_bands-1 do cp_Img[*,*,nb]=envi_get_data(FID=cp_fid,dims=cp_dims,pos=nb)
  for nb=0,fb_bands-1 do fb_Img[*,*,nb]=envi_get_data(FID=fb_fid,dims=fb_dims,pos=nb)
  ENVI_FILE_MNG,ID=cb_fid,/REMOVE
  ENVI_FILE_MNG,ID=cp_fid,/REMOVE
  ENVI_FILE_MNG,ID=fb_fid,/REMOVE
  change_C=cp_Img-cb_Img
  cp_Img=0b
  
  unmixing_State=intarr(cb_samples,cb_lines,cb_bands)
  RB=fltarr(cp_samples,cp_lines,cp_bands)
  cf_Img=MAKE_ARRAY(cb_samples,cb_lines,cb_bands,TYPE=cb_dt)
  
  ;allowed change value for each band
  min_allow=fltarr(em_samples,cb_bands)
  max_allow=fltarr(em_samples,cb_bands)
  for ib=0,cb_bands-1,1 do begin
    min_allow[*,ib]=min(change_C[*,*,ib])-stddev(change_C[*,*,ib])
    max_allow[*,ib]=max(change_C[*,*,ib])+stddev(change_C[*,*,ib])
  endfor
  ;settings for CONSTRAINED_MIN
  common V_PUB1, ind_v
  common V_PUB2, dep_v
  gbnd    =[0,100]      ;
  nobj    = 0           ;
  Lcomp   = 'Cost_Fun'  ;name of cost function
  for i=0,cb_lines-1 do begin
    ;line begin of sliding-window  Top Line
    IsTL = i - half_win_size
    ;line end   of sliding-window  Buttom Line
    IsBL = i + half_win_size
    ai   = max([0,IsTL])          ;line begin
    bi   = min([cb_lines-1,IsBL]) ;line end
    fai=i*rs_ratio
    fbi=(i+1)*rs_ratio-1
    for j=0,cb_samples-1 do begin
      ;Sample begin of sliding-window  Left Sample
      IsLS = j - half_win_size
      ;Sample end of sliding-window    Right Sample
      IsRS = j + half_win_size
      aj   = max([0,IsLS])             ;sample begin
      bj   = min([cb_samples-1,IsRS])  ;sample end
      faj=j*rs_ratio
      fbj=(j+1)*rs_ratio-1
      ;counts of coarse pixels in the sliding-window
      ;in case for pixels at the edge of the image
      c_win_pixels=(bi-ai+1)*(bj-aj+1)
      while (c_win_pixels le em_samples) do begin
        if IsTL le 0 then bi++ else ai--
        if IsLS le 0 then bj++ else aj--
        c_win_pixels=(bi-ai+1)*(bj-aj+1)
      endwhile
      ;for fine pixels within the central coarse pixel
      ;identify the endmembers with abundance greater than 0
      fabd_temp = fabd_img[faj:fbj,fai:fbi,*]
      cur_fabd  = total(reform(fabd_temp,nPixels,em_samples),1)
      cur_end_f = where(cur_fabd ne 0)
      ;abundance of coarse pixels within the current sliding-window
      cur_cabd  = total(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,em_samples),1)
      ;endmembers exist on fine image but its coarse image abundance is 0
      zero_Index=where(cur_cabd[cur_end_f] eq 0.0,zero_cnts)
      if zero_cnts gt 0 then begin
        ;merge the endmember exist on fine image but with 0 coarse abundance
        Zero_ID=cur_end_f[zero_Index]
        for zz=0,zero_cnts-1 do begin
          for sc=0,em_samples-1 do begin
            if cur_cabd[lut[sc,Zero_ID[zz]]] ne 0 then begin
              fabd_temp[*,*,lut[sc,Zero_ID[zz]]]+=fabd_temp[*,*,Zero_ID[zz]]
              fabd_temp[*,*,Zero_ID[zz]]=0.0
              break
            endif
          endfor
        endfor
      endif
      ind_v   = transpose(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,em_samples))    ;abundance by number of pixels
      for nb=0,cb_bands-1 do begin
        dep_v = double(reform(change_C[aj:bj,ai:bi,nb],c_win_pixels))   ;reflectance change
        x     = fltarr(1,em_samples)
        xbnd  = [[min_allow[*,nb]], [max_allow[*,nb]]]
        CONSTRAINED_MIN, x, xbnd, gbnd, nobj, Lcomp, inform, NSTOP = 5  ;unmixing

        ;the execution state of CONSTRAINED_MIN 
        unmixing_State[j,i,nb]=inform
        ;linear mixture
        ds_change=fabd_temp*rebin(reform(x,1,1,em_samples),rs_ratio,rs_ratio,em_samples,/SAMPLE)
        Change_F[faj:fbj,fai:fbi,nb]=total(ds_change,3)
        ;residual of spatial-unmixing
        RB[j,i,nb]=change_C[j,i,nb]-total(cabd_img[j,i,*]*x)
        
        cf_Img[j,i,nb]=mean(fb_Img[faj:fbj,fai:fbi,nb])
      endfor
    endfor
  endfor
  print,'spatial unmixing finished, time from start: '+$
    string(SYSTIME(/seconds)-timeBegin,FORMAT='(f13.6)')+' seconds'
  ;-------- Sensor difference adjustment coefficient --------
  nPixels=cb_lines*cb_samples
  fp_Img=MAKE_ARRAY(fb_samples,fb_lines,fb_bands,TYPE=fb_dt)
  print,'** Sensor difference adjustment coefficient **'
  for ib=0,cb_bands-1 do begin
    X =reform(cb_Img[*,*,ib],nPixels)
    Y =reform(cf_Img[*,*,ib],nPixels)
    coefs=LINFIT(X,Y,YFIT=Y_estimate)

    YtrueMean=mean(Y)
    YtrueVariance=TOTAL((Y-YtrueMean)*(Y-YtrueMean))
    YestimateVariance=TOTAL((Y_estimate-Y)*(Y_estimate-Y))
    r2=1-YestimateVariance/YtrueVariance
    RMSE=sqrt(mean((Y-Y_estimate)*(Y-Y_estimate)))
    print,'-->Band'+strtrim(string(ib+1),2)+string(coefs[1],format='(" a: ",f10.6)')+$
      string(r2,format='(" r2: ",f10.6)')+string(RMSE,format='(" rmse: ",f10.6)')  
    fp_Img[*,*,ib]=fb_Img[*,*,ib]+(coefs[1]*Change_F[*,*,ib])
  endfor
  print,'Sensor difference adjustment finished, time from start: '+$
    string(SYSTIME(/SECONDS)-timeBegin,FORMAT='(f13.6)')+' seconds'

  openw,lun,fpre_fp+'_unmixingState',/GET_LUN
  writeu,lun,unmixing_State
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=fpre_fp+'_unmixingState',NB=cb_bands,$
    NS=cb_samples,NL=cb_lines,BNAMES=cb_bnames,XSTART=cb_xstart,YSTART=cb_ystart,$
    INTERLEAVE=0,DATA_TYPE=2,WL=cb_WL,MAP_INFO=cb_mapinfo,/WRITE

  openw,lun,fpre_fp+'_coarseRB',/GET_LUN
  writeu,lun,RB
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=fpre_fp+'_coarseRB',NB=cb_bands,$
    NS=cb_samples,NL=cb_lines,BNAMES=cb_bnames,XSTART=cb_xstart,YSTART=cb_ystart,$
    INTERLEAVE=0,DATA_TYPE=4,WL=cb_WL,MAP_INFO=cb_mapinfo,/WRITE
  
  openw,lun,fpre_fp,/GET_LUN
  writeu,lun,fp_Img
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=fpre_fp,NB=fb_bands,XSTART=fb_xstart,YSTART=fb_ystart,$
    NS=fb_samples,NL=fb_lines,BNAMES=cb_bnames,WL=fb_WL,$
    INTERLEAVE=0,DATA_TYPE=fb_dt,MAP_INFO=fb_mapinfo,/WRITE
  
  print,'ISTRUM finished, time consumed'+string(SYSTIME(/SECONDS)-timeBegin,FORMAT='(f18.6)')+' seconds'
end
