FUNCTION Cost_Fun, g
  compile_opt idl2
  common V_PUB1
  common V_PUB2
  L=total((ind_v ## g-dep_v)^2)
  RETURN, L
END
;基于解混的模拟  STRUM
;cbase_fp     : coarse image文件路径   on base time
;cpre_fp      : coarse image文件路径   on prediction time 2
;fbase_fp     : fine image 路径   on base time
;图像数据必须为反射率数据(0~1)，无比例缩放因子
;cluster_fp   : 聚类图像路径
;rs_ratio     : resolution 之间的比值  coarse/fine
;half_win_size: 用于解混的窗口大小，与构建方程组的像素个数有关
;ab_threshold : 类别丰度小于阈值时，合并到光谱最相似的类别
;ds_fp:  fine image保存路径
;----------------------------------------
;1--> 计算不同类别在fine image 上的平均光谱，然后计算不同类别之间的之间的光谱相似性
;     在计算丰度图像时，如果某一类的丰度小于阈值，则合并到光谱最相似的类别
;2--> 计算丰度图像，同时注意对丰度小的类别进行合并
;3--> 计算两个时间的coarse image的差值图像，选择窗口像素进行解混
;     将解混后的光谱为按类别赋给原始分类图像，得到各像素的变化
;     ****************************
;     假设coarse pixel与fine pixel的比例系数为s
;     空间解混窗口（w*w）的中心coarse pixel下，对应了s*s个fine pixel
;     假设这s*s个像素中包含m个类别
;     计算coarse pixel的丰度时，如果某一类别丰度小于ab_threshold，就会被合并到相似类别中
;     假设在m个类别中，有一类j,它的丰度小于ab_threshold，对于中心coarse pixel来说，类别j被合并了
;     同时，如果类别j在整个w*w窗口内的其他coarse pixel中丰度也为零
;     构建的解混方程组中，将不包含类别j的信息
;     往分类图像赋值时，将产生错误
;     因此，采用如下解决办法：
;     将中心coarse pixel对应的cluster image 中的类别j也合并到相似类别中
;     ****************************
;4--> 将downscale的变化图像与fine image 相加
pro STRUM_by_Clu,cbase_fp,cpre_fp,fbase_fp,cluster_fp,rs_ratio,ab_threshold,half_win_size,ds_fp
  compile_opt idl2
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
  time1=SYSTIME(/SECONDS)
  
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
  print,'低分辨率丰度图像计算结束，耗时：'+string(SYSTIME(/SECONDS)-time1,FORMAT='(f18.6)')+' seconds'
  ;----------------------------------------------------------
  ;****************  解混计算     **************
  time1=SYSTIME(/SECONDS)
  Change_F =MAKE_ARRAY(clu_samples,clu_lines,cb_bands,TYPE=cb_dt)  ;downscale的图像
  cb_Img =MAKE_ARRAY(cb_samples,cb_lines,cb_bands,TYPE=cb_dt)
  cp_Img =MAKE_ARRAY(cp_samples,cp_lines,cp_bands,TYPE=cp_dt)
  for nb=0,cb_bands-1 do cb_Img[*,*,nb]=envi_get_data(FID=cb_fid,dims=cb_dims,pos=nb)
  for nb=0,cp_bands-1 do cp_Img[*,*,nb]=envi_get_data(FID=cp_fid,dims=cp_dims,pos=nb)
  ENVI_FILE_MNG,ID=cb_fid,/REMOVE
  ENVI_FILE_MNG,ID=cp_fid,/REMOVE
  
  change_C=cp_Img-cb_Img
  ;allowed change value for each band
  min_allow=fltarr(num_clu,nb)
  max_allow=fltarr(num_clu,nb)
  for ib=0,nb-1,1 do begin
    min_allow[*,ib]=min(change_C[*,*,ib])-stddev(change_C[*,*,ib])
    max_allow[*,ib]=max(change_C[*,*,ib])+stddev(change_C[*,*,ib])
  endfor
  
  unmixing_State=intarr(cb_samples,cb_lines,cb_bands)
  RB=fltarr(cp_samples,cp_lines,cp_bands)
  
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
        dep_v = double(reform(change_C[aj:bj,ai:bi,nb],c_win_pixels))
        x     = fltarr(1,num_clu)
        xbnd  = [[min_allow[*,nb]], [max_allow[*,nb]]]
        CONSTRAINED_MIN, x, xbnd, gbnd, nobj, Lcomp, inform, NSTOP = 5
                
        ;解混计算过程的精度
        unmixing_State[j,i,nb]=inform
        ;高分辨率的反射率变化
        Change_F[faj:fbj,fai:fbi,nb]=x[clu_temp-1]
        ;解混结果的残差
        RB[j,i,nb]=change_C[j,i,nb]-total(cabd_img[j,i,*]*x)
      endfor
    endfor
  endfor
  openw,lun,ds_fp+'_unmixingState',/GET_LUN
  writeu,lun,temporary(unmixing_State)
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=ds_fp+'_unmixingState',NB=cb_bands,$
    NS=cb_samples,NL=cb_lines,BNAMES=cb_bnames,XSTART=cb_xstart,YSTART=cb_ystart,$
    INTERLEAVE=0,DATA_TYPE=2,WL=cb_WL,MAP_INFO=cb_mapinfo,/WRITE
  
  openw,lun,ds_fp+'_coarseRB',/GET_LUN
  writeu,lun,RB
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=ds_fp+'_coarseRB',NB=cb_bands,$
    NS=cb_samples,NL=cb_lines,BNAMES=cb_bnames,XSTART=cb_xstart,YSTART=cb_ystart,$
    INTERLEAVE=0,DATA_TYPE=4,WL=cb_WL,MAP_INFO=cb_mapinfo,/WRITE
  
;  openw,lun,ds_fp+'_FineChange',/GET_LUN
;  writeu,lun,Change_F
;  free_lun,lun
;  ENVI_SETUP_HEAD,FNAME=ds_fp+'_FineChange',NB=cb_bands,$
;    NS=fb_samples,NL=fb_lines,BNAMES=cb_bnames,XSTART=fb_xstart,YSTART=fb_ystart,$
;    INTERLEAVE=0,DATA_TYPE=cb_dt,WL=cb_WL,MAP_INFO=fb_mapinfo,/WRITE
  
  print,'空间解混计算结束，耗时：'+string(SYSTIME(/SECONDS)-time1,FORMAT='(f18.6)')+' seconds'
  
  fp_Img_FC=fb_Img+Change_F  
  openw,lun,ds_fp,/GET_LUN
  writeu,lun,fp_Img_FC
  free_lun,lun
  ENVI_SETUP_HEAD,FNAME=ds_fp,NB=fb_bands,XSTART=fb_xstart,YSTART=fb_ystart,$
    NS=fb_samples,NL=fb_lines,BNAMES=fb_bnames,WL=fb_WL,$
    INTERLEAVE=0,DATA_TYPE=fb_dt,MAP_INFO=fb_mapinfo,/WRITE   
end

pro test_STRUM_by_Clu
  compile_opt idl2
  ;开启ENVI时段
  envi,/RESTORE_BASE_SAVE_FILES
  envi_batch_init

  
  cpre_fp  = 'F:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004187_500m'
  ;coarse image on base date
  cbase_fp = 'F:\IrinaEmelyanovaUnmix\LGC\MODIS\MOD09GA_A2004123_500m'
  ;fine image on base date
  fbase_fp   = 'F:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_123_May02'
  cluster_fp = 'F:\IrinaEmelyanovaUnmix\LGC\Landsat\L7_2004_123_May02_ISODATA_C5'
  ds_fp      = 'F:\IrinaEmelyanovaUnmix\LGC\LGC_DS\STRUM_by_Clu\2004187-by-2004123\DiffClu\MOD_2004187_Clu_C'
  
  rs_ratio     = 20
  ab_threshold = 0.05
  half_win_size=2
  number_clu=5
  STRUM_by_Clu,cbase_fp,cpre_fp,fbase_fp,cluster_fp,$
               rs_ratio,ab_threshold,half_win_size,$
               ds_fp+strtrim(string(number_clu),2)+'W'+strtrim(string(half_win_size),2)  
  ENVI_BATCH_EXIT
end