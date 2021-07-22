;基于解混的模拟 STRUM
;cbase_fp     : coarse image文件路径   on base time
;cpre_fp      : coarse image文件路径   on prediction time 2
;fbase_fp     : fine image 路径   on base time
;图像数据必须为反射率数据(0~1)，无比例缩放因子
;fbaseabd_fp  : fine abundance image on base time,不包含残差波段
;endSpec      : Spectra of endmembers applied in fine resolution image sepctral unmixing
;rs_ratio     : resolution 之间的比值  coarse/fine
;half_win_size: 用于解混的窗口大小，与构建方程组的像素个数有关
;ab_threshold : 类别丰度小于阈值时，合并到光谱最相似的类别
;ds_fp:  模拟结果fine image保存路径
;----------------------------------------
pro stdfa,cbase_fp,cpre_fp,fbase_fp,cluster_fp,rs_ratio,ab_threshold,half_win_size,ds_fp
  compile_opt idl2
  time1=SYSTIME(/SECONDS)
  ;----------------------------------------------------------
  ;******************     打     开     数     据     ******************
  ;打开coarse resolution
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
  ;打开fine image
  envi_open_file,fbase_fp,R_FID=fb_fid
  if fb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
  ;打开聚类图像
  envi_open_file,cluster_fp,R_FID=clufid
  if clufid eq -1 then begin
    envi_batch_exit
    return
  endif
  ;获取行列信息
  envi_file_query,cb_fid,dims=cb_dims,NL=cb_lines,NS=cb_samples,NB=cb_bands,BNAMES=cb_bnames,$
    DATA_TYPE=cb_dt,WL=cb_WL,XSTART=cb_xstart,YSTART=cb_ystart
  envi_file_query,cp_fid,dims=cp_dims,NL=cp_lines,NS=cp_samples,NB=cp_bands,BNAMES=cp_bnames,$
    DATA_TYPE=cp_dt,WL=cp_WL,XSTART=cp_xstart,YSTART=cp_ystart
  envi_file_query,fb_fid,dims=fb_dims,NL=fb_lines,NS=fb_samples,NB=fb_bands,BNAMES=fb_bnames,$
    DATA_TYPE=fb_dt,WL=fb_WL,XSTART=fb_xstart,YSTART=fb_ystart
  envi_file_query,clufid,dims=clu_dims,NL=clu_lines,NS=clu_samples,NB=clu_Bands,$
    XSTART=clu_xstart,YSTART=clu_ystart,NUM_CLASSES=num_clu,CLASS_NAMES=clu_names

  cb_mapinfo=ENVI_GET_MAP_INFO(FID=cb_fid)
  fb_mapinfo=ENVI_GET_MAP_INFO(FID=fb_fid)
  num_clu--   ;0为 Unclassified 不做考虑
  clu_ID=indgen(num_clu)+1
  clu_names=clu_names[1:*]
  ;----------------------------------------------------------
  ;******************     数     据     匹     配     ******************
  ;图像行列数不等
  if (cb_lines ne cp_lines) or (cb_samples ne cp_samples) or (cb_bands ne cp_bands) then begin
    print,'coarse image does not match each other'
    envi_batch_exit
    return
  endif
  if (fb_lines ne clu_lines) or (fb_samples ne clu_samples) then begin
    print,'fine image does not match cluster image'
    envi_batch_exit
    return
  endif
  if (cb_lines*rs_ratio ne fb_lines) or (cb_samples*rs_ratio ne fb_samples) then begin
    print,'coarse image does not match cluster image'
    envi_batch_exit
    return
  endif
  ;窗口内的像元数小于端元数,方程无解
  win_size =2*half_win_size+1
  pixelCnts=win_size*win_size
  if pixelCnts le (num_clu) then begin
    print,'moving window to unmixing is too small'
    envi_batch_exit
    return
  endif
  ;----------------------------------------------------------
  ;****************  计算各类别间光谱相似程度 SAM    **************
  clu_img=envi_get_data(FID=clufid,pos=0,DIMS=clu_dims)
  ENVI_FILE_MNG,ID=clufid,/REMOVE
  fb_Img=MAKE_ARRAY(fb_samples,fb_lines,fb_bands,TYPE=fb_dt)
  for i=0,fb_bands-1 do fb_Img[*,*,i]=envi_get_data(FID=fb_fid,pos=i,DIMS=fb_dims)
  ENVI_FILE_MNG,ID=fb_fid,/REMOVE

  avg_allow=MAKE_ARRAY(num_clu,fb_bands,TYPE=fb_dt)
  for i=1,num_clu do begin
    index=where(clu_img eq i)
    for j=0,fb_bands-1 do avg_allow[i-1,j]=mean((fb_Img[*,*,j])[index])
  endfor
  ;计算类别光谱的相似性
  lut=uintarr(num_clu,num_clu)
  for i=0,num_clu-1 do begin
    cosSAM=fltarr(num_clu)
    for j=0,num_clu-1 do begin
      if j eq i then continue
      cosSAM[j]=total(avg_allow[i,*]*avg_allow[j,*])/(SQRT(TOTAL(avg_allow[i,*]^2))*SQRT(TOTAL(avg_allow[j,*]^2)))
    endfor
    lut[*,i]=reverse(sort(cosSAM))
  endfor
  ;----------------------------------------------------------
  ;****************  丰度计算     **************
  cabd_img=fltarr(cb_samples,cb_lines,num_clu)   ;丰度图像
  for clu_i=0,cb_lines-1 do begin
    for clu_j=0,cb_samples-1 do begin
      win_data=clu_img[(clu_j*rs_ratio):((clu_j+1)*rs_ratio-1),(clu_i*rs_ratio):((clu_i+1)*rs_ratio-1)]
      pdf=HISTOGRAM(win_data,BINSIZE=1,LOCATIONS=xlon,min=1,max=num_clu)
      abd=pdf/total(pdf)
      index=where((abd gt 0),cnts)
      abd_min=min(abd[index],min_ind)
      while abd_min lt ab_threshold do begin
        min_ind=index[min_ind]
        abdtemp=abd[min_ind]
        abd[min_ind]=0.0
        for sc=0,num_clu -1 do begin
          if abd[lut[sc,min_ind]] ne 0 then begin
            abd[lut[sc,min_ind]]+=abdtemp
            break
          endif
        endfor
        index=where((abd gt 0),cnts)
        abd_min=min(abd[index],min_ind)
      endwhile
      cabd_img[clu_j,clu_i,*]=abd
    endfor
  endfor
  fd_fp=file_dirname(ds_fp)+'\'+file_basename(cluster_fp)+'_cfd'
  openw,lun,fd_fp,/GET_LUN
  writeu,lun,cabd_img
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=fd_fp,NB=num_clu,NS=cb_samples,NL=cb_lines,$
    XSTART=cb_xstart,YSTART=cb_ystart,BNAMES=clu_names[1:*],$
    INTERLEAVE=0,DATA_TYPE=4,MAP_INFO=cb_mapinfo,/WRITE
  ;----------------------------------------------------------
  ;****************  解混计算     **************  
  cb_Img =MAKE_ARRAY(cb_samples,cb_lines,cb_bands,TYPE=cb_dt)
  cp_Img =MAKE_ARRAY(cp_samples,cp_lines,cp_bands,TYPE=cp_dt)
  for nb=0,cb_bands-1 do cb_Img[*,*,nb]=envi_get_data(FID=cb_fid,dims=cb_dims,pos=nb)
  for nb=0,cp_bands-1 do cp_Img[*,*,nb]=envi_get_data(FID=cp_fid,dims=cp_dims,pos=nb)
  ENVI_FILE_MNG,ID=cb_fid,/REMOVE
  ENVI_FILE_MNG,ID=cp_fid,/REMOVE
  ;downscale的图像
  ds_cb=MAKE_ARRAY(clu_samples,clu_lines,cb_bands,TYPE=cb_dt)  
  ds_cp=MAKE_ARRAY(clu_samples,clu_lines,cb_bands,TYPE=cb_dt)
  ;allowed range value for each band
  min_allow=fltarr(num_clu,nb)
  max_allow=fltarr(num_clu,nb)+1.0
  ;解混参数设置
  common V_PUB1, ind_v
  common V_PUB2, dep_v
  gbnd    =[0,100]
  nobj    = 0
  Lcomp   = 'Cost_Fun'
  nPixels = uint(rs_ratio)*rs_ratio
  for i=0,cb_lines-1 do begin
    ;判断是否超出起始行Top Line
    IsTL = i - half_win_size
    ;判断是否超出最后一行Buttom Line
    IsBL = i + half_win_size
    ai   = max([0,IsTL])          ;line begin
    bi   = min([cb_lines-1,IsBL]) ;line end
    fai=i*rs_ratio
    fbi=(i+1)*rs_ratio-1
    for j=0,cb_samples-1 do begin
      ;判断是否超出起始列 Left Sample
      IsLS = j - half_win_size
      ;判断是否超出最后一列 Right Sample
      IsRS = j + half_win_size
      aj   = max([0,IsLS])             ;sample begin
      bj   = min([cb_samples-1,IsRS])  ;sample end
      faj=j*rs_ratio
      fbj=(j+1)*rs_ratio-1

      ;窗口内coarse pixel的个数
      c_win_pixels=(bi-ai+1)*(bj-aj+1)
      while (c_win_pixels le num_clu) do begin
        if IsTL le 0 then bi++ else ai--
        if IsLS le 0 then bj++ else aj--
        c_win_pixels=(bi-ai+1)*(bj-aj+1)
      endwhile

      ;当前 coarse pixel 下对应分类图像
      clu_temp=clu_Img[faj:fbj,fai:fbi]
      cur_cluster=reform(clu_temp,nPixels)
      ;当前coarse pixel 下的类别
      cur_cluster=cur_cluster[UNIQ(cur_cluster, SORT(cur_cluster))]

      ;当前 WIN_SIZE*WIN_SIZE 窗口内的丰度
      abdTemp=total(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,num_clu),1)
      ;检查是否当前m个端元类别中，是否每个类别的丰度都不为零
      zero_Index=where(abdTemp[cur_cluster-1] eq 0.0,zero_cnts)

      if zero_cnts gt 0 then begin
        ;丰度为0的类别编号，1开始编号
        Zero_ID=cur_cluster[zero_Index]
        for zz=0,zero_cnts-1 do begin
          for sc=0,num_clu-1 do begin
            if abdTemp[lut[sc,Zero_ID[zz]-1]] ne 0 then begin
              clu_temp[where(clu_temp eq Zero_ID[zz])]=clu_ID[lut[sc,Zero_ID[zz]-1]]
              break
            endif
          endfor
        endfor
      endif

      ind_v=transpose(reform(cabd_img[aj:bj,ai:bi,*],c_win_pixels,num_clu))
      for nb=0,cb_bands-1 do begin
        xbnd  = [[min_allow[*,nb]], [max_allow[*,nb]]]
        
        dep_v = double(reform(cb_Img[aj:bj,ai:bi,nb],c_win_pixels))
        x     = fltarr(1,num_clu)        
        CONSTRAINED_MIN, x, xbnd, gbnd, nobj, Lcomp, inform, NSTOP = 5
        ds_cb[faj:fbj,fai:fbi,nb]=x[clu_temp-1]
        
        dep_v = double(reform(cp_Img[aj:bj,ai:bi,nb],c_win_pixels))
        x     = fltarr(1,num_clu)
        CONSTRAINED_MIN, x, xbnd, gbnd, nobj, Lcomp, inform, NSTOP = 5
        ds_cp[faj:fbj,fai:fbi,nb]=x[clu_temp-1]
      endfor
    endfor
  endfor 
  Change_F =ds_cp-ds_cb  ;downscale的图像 
  openw,lun,ds_fp+'_01FineChange',/GET_LUN
  writeu,lun,Change_F
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=ds_fp+'_01FineChange',NB=cb_bands,$
    NS=fb_samples,NL=fb_lines,BNAMES=cb_bnames,XSTART=fb_xstart,YSTART=fb_ystart,$
    INTERLEAVE=0,DATA_TYPE=cb_dt,WL=cb_WL,MAP_INFO=fb_mapinfo,/WRITE

  fp_Img_FC=fb_Img+Change_F
  openw,lun,ds_fp,/GET_LUN
  writeu,lun,fp_Img_FC
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=ds_fp,NB=fb_bands,XSTART=fb_xstart,YSTART=fb_ystart,$
    NS=fb_samples,NL=fb_lines,BNAMES=fb_bnames,WL=fb_WL,$
    INTERLEAVE=0,DATA_TYPE=fb_dt,MAP_INFO=fb_mapinfo,/WRITE
  print,file_basename(ds_fp)+'计算结束，耗时：'+string(SYSTIME(/SECONDS)-time1,FORMAT='(f18.6)')+' seconds'
end

