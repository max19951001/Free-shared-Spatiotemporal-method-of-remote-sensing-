
;;; function for open the file
; -------------------------------------------------------------------------------------------------------------------
Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
  FileName = FileName,Map_info = map_Info, Fid = Fid, dims = dims
  Envi_Open_File,FileName,R_Fid = Fid
  Envi_File_Query,Fid,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type
  map_info = envi_get_map_info(fid=Fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case Data_Type Of
    1:ImgData = BytArr(ns,nl,nb)    ;  BYTE  Byte
    2:ImgData = IntArr(ns,nl,nb)    ;  INT  Integer
    3:ImgData = LonArr(ns,nl,nb)    ;  LONG  Longword integer
    4:ImgData = FltArr(ns,nl,nb)    ;  FLOAT  Floating point
    5:ImgData = DblArr(ns,nl,nb)    ;  DOUBLE  Double-precision floating
    6:ImgData = COMPLEXARR(ns,nl,nb); complex, single-precision, floating-point
    9:ImgData = DCOMPLEXARR(ns,nl,nb);complex, double-precision, floating-point
    12:ImgData = UINTARR(ns,nl,nb)   ; unsigned integer vector or array
    13:ImgData = ULONARR(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:ImgData = LON64ARR(ns,nl,nb)   ;a 64-bit integer vector or array
    15:ImgData = ULON64ARR(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  EndCase
  For i = 0,nb-1 Do Begin
    Dt = Envi_Get_Data(Fid = Fid,dims = dims,pos=i)
    ImgData[*,*,i] = Dt[*,*]
  EndFor
End
;-------------------------------------------------------------------------------------------------------


;;; read date and time from Landsat data
; LC8126038201287LGN00_ndvi
;----------------------------------------------------------------
function GetDay, filename, t_deta
  name = file_basename(filename)
  location1 = strpos(name, '_', /REVERSE_SEARCH)        ; return the last location of the sign _
  DOY_tm = strmid(name, location1-8, 3)                 ; get the DOY
  DOY_tm=fix(DOY_tm)
  band_tm=ceil(1.0*DOY_tm/t_deta)                       ; get the sequence

  return, [band_tm, DOY_tm]
end


;;; read date and time from Sentinel data
; L2A_T48QUD_20170305T032611_ndvi
; 要求：要求输入的Landsat数据名称是标准的Landsat数据名，以NDVI或者fmask结尾
;----------------------------------------------------------------------------------
function GetDay2, filename, t_deta, year, month, day
  name = file_basename(filename)
  location1 = strpos(name, '_', /REVERSE_SEARCH)        ; return the first location of underscore
  year = strmid(name, location1-15, 4)
  yea = fix(year)
  month = strmid(name, location1-11, 2)
  mon = fix(month)
  day = strmid(name, location1-9, 2)
  da = fix(day)
  DOY_tm = julday(mon, da, yea) - julday(1, 1, yea) + 1
  ;DOY_tm = strmid(name, location1-8, 3)                 ; get the DOY
  DOY_tm=fix(DOY_tm)
  band_tm=ceil(1.0*DOY_tm/t_deta)                           ; get the sequence
  
  return, [band_tm, DOY_tm]
end
;------------------------------------------------------------------------------------


;;; calculate pixel purity
; ------------------------------------------------------------------------------------------------------
pro GetHindex, Class=Class, Ratio=ratio, Hindex=het_index

  Csize=size(Class)
  ns=Csize(1)
  nl=Csize(2)
  het_index=fltarr(ns,nl)
  scale_d=ratio/2.0               ; 当前像元的纯度是在以当前像元为中心的一个Modis像元大小的窗口内计算的
  for i=0,nl-1, 1 do begin        ; the window location
    for j=0,ns-1, 1 do begin
      ai=max([0,i-scale_d])       ; 窗口的左边界
      bi=min([nl-1,i+scale_d-1])    ; 窗口的右边界
      aj=max([0,j-scale_d])       ; 窗口的上边界
      bj=min([ns-1,j+scale_d-1])    ; 窗口的下边界
      class_t=class[j, i]

      ;select same-class pixels
      ind_same_class=where(class[aj:bj, ai:bi] eq class_t, num_sameclass)
      het_index[j, i]=float(num_sameclass)/((bi-ai+1.0)*(bj-aj+1.0))
    endfor
  endfor

end
;-------------------------------------------------------------------------------------------------------



;;; Constrained Least Square(CLS) combination
;--------------------------------------------------------------------------------------------------------
function CLS_main, dlt_coarse, dlt_c, dlt_t, dlt_loc, ratio, mask_shrink_c
  ; dlt_coarse 表示MODIS像元在两期的差值
  ; dlt_c 表示空间增量
  ; dlt_t 表示时间增量
  
  ;;; common variables for weight conbination
  common V_PUB3, Model
  common V_PUB4, Obs
  
  dlt_com = dlt_t                                              ; give an initial value for the combined increment
  f_size = size(dlt_c)
  ns = f_size[1]
  nl = f_size[2]
  dltc_size = size(dlt_coarse)
  ns_coarse = dltc_size(1)
  nl_coarse = dltc_size(2)
  dlt_c_Coarse = fltarr(ns_coarse, nl_coarse)                  ; 粗像元尺度的空间增量
  dlt_t_Coarse = fltarr(ns_coarse, nl_coarse)                  ; 粗像元尺度的时间增量

  ;;; upscale the fine increment to coarse pixel
  for imsk=0, nl_coarse-1 do begin
    for jmsk=0, ns_coarse-1 do begin
      if mask_shrink_c[jmsk, imsk] eq 0 then begin             ; this coarse pixel is a complete pixel
        dlt_c_coarse[jmsk, imsk]=mean(dlt_c[ratio*jmsk-dlt_loc[0]:ratio*(jmsk+1)-1-dlt_loc[0], ratio*imsk-dlt_loc[1]:ratio*(imsk+1)-1-dlt_loc[1]])
        dlt_t_coarse[jmsk, imsk]=mean(dlt_t[ratio*jmsk-dlt_loc[0]:ratio*(jmsk+1)-1-dlt_loc[0], ratio*imsk-dlt_loc[1]:ratio*(imsk+1)-1-dlt_loc[1]])
      endif
    endfor
  endfor

  ;;; calculate the weights in moving window
  patch_size_c = 7                                           ; patch size of coarse pixels
  n_ns_coarse = ceil(float(ns_coarse)/patch_size_c)
  n_nl_coarse = ceil(float(nl_coarse)/patch_size_c)
  ;;;w_patch = fltarr(2, n_nl_coarse*n_ns_coarse)
  for inl=0,nl_coarse-1 do begin
    for ins=0,ns_coarse-1 do begin

      p_dlt_coarse = dlt_coarse[max([0, ins-patch_size_c/2]):min([ins+patch_size_c/2,ns_coarse-1]), max([0,inl-patch_size_c/2]):min([inl+patch_size_c/2,nl_coarse-1])]
      p_mask_shrink = mask_shrink_c[max([0, ins-patch_size_c/2]):min([ins+patch_size_c/2,ns_coarse-1]), max([0,inl-patch_size_c/2]):min([inl+patch_size_c/2,nl_coarse-1])]
      idx_coarse = where(p_mask_shrink eq 0, num_coarse, /null)
      
      if num_coarse gt 0 then begin                           ; there are coarse pixels in this patch
        if num_coarse gt 20 then begin                       ; greater than 100 coarse pixels, then use BMA
          p_dlt_c_coarse = dlt_c_coarse[max([0, ins-patch_size_c/2]):min([ins+patch_size_c/2,ns_coarse-1]), max([0,inl-patch_size_c/2]):min([inl+patch_size_c/2,nl_coarse-1])]
          p_dlt_t_coarse = dlt_t_coarse[max([0, ins-patch_size_c/2]):min([ins+patch_size_c/2,ns_coarse-1]), max([0,inl-patch_size_c/2]):min([inl+patch_size_c/2,nl_coarse-1])]
          Obs = p_dlt_coarse[idx_coarse]
          Model = fltarr(n_elements(Obs), 2)
          Model[*, 0] = p_dlt_c_coarse[idx_coarse]
          Model[*, 1] = p_dlt_t_coarse[idx_coarse]
          
          ;;; use constrained least square method
          xlb=fltarr(2, 1) 
          xub=fltarr(2, 1) + 1                    
          xbnd = [[xlb], [xub]]
          x = [0.5, 0.5]
          nobj = 1
          gcomp = 'SQUARE_WEIGHT'
          gbnd = [[1,0], [1,0]]
;          gbnd = [[0], [0]]
          CONSTRAINED_MIN, x, xbnd, gbnd, nobj, gcomp, inform, EPSTOP = 1.0e-8, NSTOP = 10 
          result = x
;          result=IMSL_LINLSQ(Obs, Model, c, bc, bc, contype, Xlb = xlb, Xub = xub, /double)
          
          ;;; 二次计算权重
          Model[*, 0] = Model[*, 0] + 0.15*result[0]*stddev(Model[*, 0])*randomn(seed, n_elements(Obs), 1)
          Model[*, 1] = Model[*, 1] + 0.15*result[1]*stddev(Model[*, 0])*randomn(seed, n_elements(Obs), 1)
          CONSTRAINED_MIN, x, xbnd, gbnd, nobj, gcomp, inform, EPSTOP = 1.0e-8, NSTOP = 10
          result = x
;          result=IMSL_LINLSQ(Obs, Model, c, bc, bc, contype, Xlb = xlb, Xub = xub, /double)                           
          w = result  
        endif else begin
          w = [0.5, 0.5]                                      ; use default value
        endelse
      endif else begin                                        ; no coarse pixels in this patch
        continue                                              ; skip the following codes in this loop
      endelse

      ;;; combine the spatial increment and the temporal increment
      p_dlt_c = $
        dlt_c[max([ins*ratio-dlt_loc[0], 0]):min([(ins+1)*ratio-1-dlt_loc[0], ns-1]),$
        max([inl*ratio-dlt_loc[1], 0]):min([(inl+1)*ratio-1-dlt_loc[1], nl-1])]
      p_dlt_t = $
        dlt_t[max([ins*ratio-dlt_loc[0], 0]):min([(ins+1)*ratio-1-dlt_loc[0], ns-1]),$
        max([inl*ratio-dlt_loc[1], 0]):min([(inl+1)*ratio-1-dlt_loc[1], nl-1])]
      dlt_com[max([ins*ratio-dlt_loc[0], 0]):min([(ins+1)*ratio-1-dlt_loc[0], ns-1]),$
        max([inl*ratio-dlt_loc[1], 0]):min([(inl+1)*ratio-1-dlt_loc[1], nl-1])] = w[0]*p_dlt_c+w[1]*p_dlt_t

    endfor
  endfor

  return, dlt_com                                              ; output the combined increment
end
;--------------------------------------------------------------------------------------------------------


;;; add a sub-function for constrained_min in weight combination
; 2019.01.01
;---------------------------------------------------------------
function SQUARE_WEIGHT, x
  common V_PUB3
  common V_PUB4
  
  g = dblarr(2)
  g[0] = x[0] + x[1]
  g[1] = total((Model#x - transpose(Obs))^2)
  return, g
end


;;; Thin Spline interpolation(TPS) 使用划分子块儿的方式来对数据进行TPS插值:保留中间的3*3子块儿
; -------------------------------------------------------------------------------------------------------------------------
function moveWinsub_TPS, coarseImg, ratio, win
  ;;; ratio表示Landsat相对于MODIS的比例, win表示分块儿的窗口大小
  Mod_size = size(coarseImg)
  ns = Mod_size[1]
  nl = Mod_size[2]
  winsub = win-2
  nl_patch = ceil(1.0*nl/winsub)                                                               ; 子块儿数目
  ns_patch = ceil(1.0*ns/winsub)                                                               ; 子块儿数目
  coarse0_TPS = fltarr(ratio*ns, ratio*nl)
  for i=0, ns_patch-1 do begin                                                              ; 对子块儿进行循环
    for j=0, nl_patch-1 do begin
      ai = max([0,i*winsub-1])
      aj = min([ns-1, (i+1)*winsub])
      bi = max([0,j*winsub-1])
      bj = min([nl-1, (j+1)*winsub])
      aad = 0
      bbd = 0
      if aj-ai lt 2 then begin                                                 ; 一定是最靠右的时候像元数目不足
        ai = ai-1
        aad = 1
      endif
      if bj-bi lt 2 then begin                                                 ; 一定是最靠左的时候像元数目不足
        bi = bi-1
        bbd = 1
      endif
      small_dt = coarseImg[ai:aj, bi:bj]

      ;small_dt = coarseImg[i*win:min([ns, (i+1)*win])-1, j*win:min([nl, (j+1)*win])-1]
      xgrid = [ratio/2.0, ratio]
      ygrid = [ratio/2.0, ratio]
      xout = indgen(ratio*(aj-ai+1)) + 0.5
      yout = indgen(ratio*(bj-bi+1)) + 0.5
      tmpTPS = MIN_CURVE_SURF(small_dt, xgrid=xgrid, ygrid=xgrid, /TPS, xout=xout, yout=yout)

      if (aj-ai lt win-1) or (bj-bi lt win-1) then begin                       ; 窗口大小不足win*win
        aii = ratio
        bii = ratio
        if aj-ai lt win-1 then begin                                           ; 列数不足win
          if ai eq 0 then begin                                                ; 靠左
            aii = 0
          endif else begin                                                     ; 靠右
            aii = ratio
          endelse
          if aad eq 1 then begin                                               ; 人为增加了像元数目
            aii = aii + ratio*1
          endif
        endif

        if bj-bi lt win-1 then begin                                           ; 行数不足win
          if bi eq 0 then begin                                                ; 靠上
            bii = 0
          endif else begin                                                     ; 靠下
            bii = ratio
          endelse
          if bbd eq 1 then begin                                               ; 人为增加了像元数目
            bii = bii + ratio*1
          endif
        endif
        tmpTPS = tmpTPS[aii:aii+ratio*min([ns, (i+1)*winsub])-i*winsub*ratio-1, bii:bii+ratio*min([nl, (j+1)*winsub])-j*winsub*ratio-1]
      endif else begin
        tmpTPS = tmpTPS[ratio:ratio*(win-1)-1, ratio:ratio*(win-1)-1]
      endelse

      coarse0_TPS[i*winsub*ratio:ratio*min([ns, (i+1)*winsub])-1, j*winsub*ratio:ratio*min([nl, (j+1)*winsub])-1] = tmpTPS
    endfor
  endfor

  return, coarse0_TPS
end
; -------------------------------------------------------------------------------------------------------------------------




;;; make preparations for unmixing: there are MODIS pixels with the same class but different increment
; -----------------------------------------------------------------------------------------------------
pro pre_unmix, wind_K, wind_frac, n_N, a=A, b=b, num_loc=num_loc, num_cls = num_cls

  win_size = size(wind_K)
  ns = win_size[1]
  nl = win_size[2]
  A=fltarr(ns*nl, n_N)
  b=fltarr(ns*nl)
  num_y =  1+fltarr(n_elements(b))
  tmp_a = fltarr(1, n_N)

  for je=0,nl-1 do begin                           ; 表示的行位置
    for ie=0,ns-1 do begin                         ; 表示的列位置

      ;;; do one more action: decide wheather the new fraction vector is existed in matrix A
      tmp_a[0,*] = wind_frac[ie,je,*]              ; 实际上是对丰度进行遍历
      if total(abs(tmp_a)) eq 0 then continue      ; do not condsider the background
      diff_A = A - tmp_a##(1+fltarr(n_elements(b)))
      sum_diff = total(abs(diff_A), 2)             ; absolute value
      idx_zero = where(sum_diff eq 0, cnt_zero)
      if cnt_zero eq 0 then begin                  ; the new fraction does not exist
        ; 对矩阵进行这样的赋值也是可行的
        A[ie+je*ns,*] = wind_frac[ie,je,*]
        b[ie+je*ns] = wind_K[ie,je]
      endif else begin
        num_y[idx_zero] = num_y[idx_zero] + 1      ; num_y的初始值是1
        b[idx_zero] = b[idx_zero] + wind_K[ie,je]
      endelse
    endfor
  endfor
  b=b/num_y                                   ; 同一个丰度下，有多个增量值时，将增量取平均
  loc_idx = where(b ne 0, num_loc)            ; num_loc 表示有效方程的个数
  cls_idx = where(total(A, 1) ne 0, num_cls)
  
  ;;; 将A和b精简后输出
  A = A[loc_idx, *]
  b = b[loc_idx]
end
;-----------------------------------------------------------------------------------



;;; pixel unmixing: using coarse increments and fraction to unmix coarse pixels
; 在确定背景的时候基于分类结果来确定，不是基于NDVI数据确定
; 实现混合像元分解的功能，避免函数主体过分复杂
; ------------------------------------------------------------------------------------------------------
pro GetTemp_Increment, frac_coarse=frac_coarse, K_fine=K_fine, K_coarse, dlt_loc, $
  w_neighbor, classFine1, fine1, het_index, mask_shrink_c = mask_shrink_c, ratio
  ; 解算出丰度以及增长率,通过传参输出

  ;;; set common two variables
  common V_PUB1, A
  common V_PUB2, b

  tw_neighbor = w_neighbor
  C_size=size(K_coarse)                                             ; 像元在MODIS尺度的增量
  ns_coarse=C_size[1]                                               ; 粗像元尺寸
  nl_coarse=C_size[2]
  f_size=size(classFine1)
  ns=f_size[1]                                                      ; 细像元尺寸
  nl=f_size[2]
  n_N=max(classFine1)
  K_fine=fltarr(ns_coarse,nl_coarse,n_N)                            ; 对每一个Modis像元内的不同地物类型的增量进行确定
  K_fine = K_fine + 255                                             ; 将填充值设置成255
  frac_coarse=fltarr(ns_coarse,nl_coarse,n_N)                       ; n_N表示分类的类别数目
  mask_ex_c = bytarr(ns_coarse,nl_coarse)                           ; MODIS尺度的mask文件，略有扩大
  mask_shrink_c = mask_ex_c                                         ; MODIS尺度的mask文件，略有缩小
  use = 1.0

  ;;; 根据偏移量，生成虚拟的mask文件
  virtual_mask = 255 + bytarr([ns_coarse, nl_coarse]*ratio)
  virtual_mask[dlt_loc[0]:ns+dlt_loc[0]-1, dlt_loc[1]:nl+dlt_loc[1]-1] = classFine1 eq 0   ;;; 不等于0就是有用值

  ; get the abundance of each MODIS pixel 在分类结果上统计丰度，mask文件中背景值是1，有效值是0
  for j=0, nl_coarse-1 do begin
    for i=0, ns_coarse-1 do begin
      tmp_mask_MOD = virtual_mask[ratio*i:ratio*(i+1)-1, ratio*j:ratio*(j+1)-1] 
      idx255 = where(tmp_mask_MOD eq 255, cnt255)                                   ; 取值为背景值255的像元个数                  
      if cnt255 gt 0 then begin                                                     ; 表明这个MOD中有fine像元不在研究区内
        mask_ex_c[i, j] = 1                                                         ; 这个MOD像元先不使用
        mask_shrink_c[i, j] = 1
      endif else begin                                                              ; 像元全都在研究区内
        idx0 = where(tmp_mask_MOD eq 0, cnt0)                                       ; 统计有用像元数目
        if cnt0 gt 0 then begin                                                     ; 是有用的MOD像元
          if (cnt0 lt use*ratio^2) then begin                                       ; 有用fine像元的数目小于0.90
            mask_shrink_c[i, j] = 1                                                 ; 可用fine像元太少
          endif else begin                                                          ; this coarse pixel can be used in unmixing
            ;; 对于有用的MOD像元，统计类型丰度
            lcMODIS = classFine1[ratio*i-dlt_loc[0]:ratio*(i+1)-1-dlt_loc[0], ratio*j-dlt_loc[1]:ratio*(j+1)-1-dlt_loc[1]]
            num=TOTAL(lcMODIS ne 0.0)
            for t=1, n_N do begin
              index=where(lcMODIS eq t, numI)
              frac_coarse[i, j, t-1]=(1.0*numI)/num                                    ; 计算出有用MOD中类型的丰度
            endfor   
                     
          endelse                    
        endif else begin                                                             ; 表示取值全是1，是背景
          mask_ex_c[i, j] = 1
          mask_shrink_c[i, j] = 1
        endelse        
      endelse
    endfor
  endfor

  stdK_coarse=stddev(K_coarse[where(mask_ex_c eq 0)])    ; 粗尺度像元的增量是不规则输入的
  minK_coarse=min(K_coarse[where(mask_ex_c eq 0)])       ; 0 表示有效值
  maxK_coarse=max(K_coarse[where(mask_ex_c eq 0)])  

  ;calculate the growing rate of TM by spectral unmixing
  for j=0,nl_coarse-1,1 do begin
    for i=0,ns_coarse-1,1 do begin                       ; retrieve each target pixel in coarse image
      ; mask的值为0就是数据，MODIS数据一般覆盖比较大，开窗运算应该是可以的
      if (mask_shrink_c[i, j] eq 0) then begin           ; use MODIS pixels which cover full Landsat pixels
        ai=max([0,i-w_neighbor])                         ; the MODIS window location
        bi=min([ns_coarse-1,i+w_neighbor])               ; 预设的窗口大小为3*3，w取值为1
        aj=max([0,j-w_neighbor])
        bj=min([nl_coarse-1,j+w_neighbor])

        wind_K=K_coarse[ai:bi,aj:bj]                     ; K value in the window, increment of MODIS pixels        
        wind_frac=frac_coarse[ai:bi,aj:bj,*]             ; fraction value in the window
        A=fltarr((bi-ai+1)*(bj-aj+1),n_N)
        b=fltarr((bi-ai+1)*(bj-aj+1))        
        pre_unmix, wind_K, wind_frac, n_N, a=A, b=b, num_loc=num_loc, num_cls = num_cls    ; remove background and simplify the linear equations

        while num_loc lt 2*num_cls do begin                    ; the number of equations is less than the number of classes
          tw_neighbor = tw_neighbor + 1
          ai=max([0,i-tw_neighbor])                         ; the MODIS window location
          bi=min([ns_coarse-1,i+tw_neighbor])               ; 预设的窗口大小为3*3，w取值为1
          aj=max([0,j-tw_neighbor])
          bj=min([nl_coarse-1,j+tw_neighbor])
          wind_K=K_coarse[ai:bi,aj:bj]                     ; K value in the window
          wind_frac=frac_coarse[ai:bi,aj:bj,*]        
          pre_unmix, wind_K, wind_frac, n_N, a=A, b=b, num_loc=num_loc, num_cls = num_cls
        endwhile
        tw_neighbor = w_neighbor
                
        xlb=fltarr(n_N,1)+ minK_coarse-stdK_coarse
        xub=fltarr(n_N,1)+ maxK_coarse+stdK_coarse
        
        ;;; 2018.01.01 test constrained_min function
        xbnd = [[xlb], [xub]]
        x = fltarr(n_N) + 1.0/n_N
        nobj = 0
        gcomp = 'SQUARE_ERROR'
        gbnd = [0, n_elements(b)]        
        CONSTRAINED_MIN, x, xbnd, gbnd, nobj, gcomp, inform, EPSTOP = 1.0e-8
        result = x        
        
;        result=IMSL_LINLSQ(b, A, c, bc, bc, contype, Xlb = xlb, Xub = xub, /double)       ; 带约束最小二乘法
        num_nan = finite(result)
        result[where(num_nan eq 0, /null)] = 255         ; 将解混得到的无效值都设置为255,不予考虑
        K_fine[i,j,*]= result[*]                         ; 每个MODIS像元都解算出了一组增量，因此每个窗口中的增量都保留下来了
      endif
    endfor
  endfor
end
;-------------------------------------------------------------------------------------------------------


;;; a sub-function for constrained_min in unmixing
;---------------------------------------------------------------
function SQUARE_ERROR, x
  common V_PUB1
  common V_PUB2
  g = total((A#x - transpose(b))^2)
  return, g
end


;;; get available fine images for current prediction date: only a small subset will be selected based on time range
;;; 这是一个搜索函数，主要是从辅助数据中确定哪些辅助数据可用，将其记录下来：这个函数可以继续优化
;说明：以前的搜索过程是按照闭合的原则进行的，现在不能用这种处理；设定搜索的边界，将单侧最大搜索影像数目设置为3
;-----------------------------------------------------------------------------------------------------
pro GetAncillaryNDVI, band_tm, day_tm, predict, dt_mask, an_time, an_day, t_deta, an_mask, anci_need=anci_need, $
  out_mask = out_mask, cur_anci = cur_anci, loc = loc, time_range
  ; band_tm 表示参照数据的时间，即第几幅，是1-23之间的值
  ; predict 表示预测的波段编号，即第几幅
  ; an_time 表示辅助数据的编号
  ; an_mask 表示辅助数据的mask文件
  ; anci_need 是最终挑选出来的可以使用的辅助数据，记录的是1-23之间的编号
  ; out_mask 是对an_mask进行调整过的mask文件，去掉了交集区域
  ; cur_anci 记录是否当前预测的时间正好处于辅助数据中
  ; loc 用来指示是否左右的mask都有，-1表示只有左边，1表示只有右边，0表示左右都有

  loc = 0
  cur_anci = 0
  an_size = size(an_mask)
  ns = an_size[1]
  nl = an_size[2]
  tot_DOY = [an_day, day_tm]
  tot_time = [an_time, band_tm]                                     ; 辅助数据与参考数据时间的集合
  c_tot_DOY = tot_DOY
  tidx = where(tot_time eq predict, cnt)
  if cnt ne 0 then begin
    c_tot_DOY[tidx] = 10000                                         ; 不可能有10000，排除掉这个值
    cur_anci = predict
  endif  
  
  dlt_DOY = c_tot_DOY - (t_deta/2+t_deta*(predict-1))
  idx_DOY = where(abs(dlt_DOY) le time_range, /null, c_D)           ; get the index of DOY within two months from predictoin date
  t_range = time_range
  while c_D eq 0 do begin
    t_range = t_range + 30                                          ; 将时间增加一个月，继续寻找辅助数据
    idx_DOY = where(abs(dlt_DOY) le t_range, /null, c_D)
  endwhile
    
  tot_time = tot_time[idx_DOY]                                      ; only retain the date within two months from prediction date
;  tidx = where(tot_time eq predict, cnt)
;  if cnt ne 0 then cur_anci = predict                               ; 表明当前景数据正好位于辅助数据中，因为肯定不是参考数据

  ; get the minimum ancillary image 挑选左边的数据 before
  left = tot_time[where(tot_time lt predict, cnt, /null)]           ; 得到的都是负数，不计算0
  left_mask = bytarr(ns, nl) + 255                                  ; 左边的最终mask文件，标记的是波段编号，填充值为255
  if cnt ne 0 then begin                                            ; 只有不为0时才需要判断mask文件
    left_time = left[sort(left)]                                    ; 得到升序排列的数据 [-5，-4，-3] 等等
    for i=0, cnt-1 do begin                                         ; 实际上从远处往近处搜索，这样可以保证时间近的数据能够覆盖时间远的数据
      if left_time[i] eq band_tm then begin
        tmp_mask = dt_mask
      endif else begin
        idx = where(an_time eq left_time[i])
        tmp_mask = an_mask[*, *, idx]                               ; tmp_mask 是一个0-1标记的文件
      endelse

      left_mask = left_mask*tmp_mask                                ; 255作为填充值
      idx_flaged = where(left_mask eq 0)
      left_mask[idx_flaged] = left_time[i]                          ; 最后标记矩阵中剩下的255都是缺值的像元
    endfor
    loc = loc -1                                                    ; 有左边的mask文件    
  endif

  ; 挑选右边的数据 after
  right = tot_time[where(tot_time gt predict, cnt, /null)]           ; 得到的都是正数，不计算0 
  right_mask = bytarr(ns, nl) + 255                                  ; 右边的最终mask文件，标记的是波段编号，初始化为255
  if cnt ne 0 then begin                                             ; 只有不为0时才需要判断mask文件
    right_time = right[reverse(sort(right))]                         ; 最后将波段的编号逆序排列[5,4,3,2]
    for i=0, cnt-1 do begin                                          ; 实际上从远处往近处搜索
      if right_time[i] eq band_tm then begin
        tmp_mask = dt_mask
      endif else begin
        idx = where(an_time eq right_time[i])
        tmp_mask = an_mask[*, *, idx]                               ; tmp_mask 是一个0-1标记的文件
      endelse

      right_mask = right_mask*tmp_mask
      idx_flaged = where(right_mask eq 0)
      right_mask[idx_flaged] = right_time[i]                        ; 最后标记矩阵中剩下的255都是缺值的像元
    endfor
    loc = loc + 1                                                   ; 有右边的mask文件
  endif

  out_mask = bytarr(ns, nl, 2)
  out_mask[*, *, 0] = left_mask                                     ; 将最后找到的mask记录输出
  out_mask[*, *, 1] = right_mask                                    ; 如果left和right有相同的标记，表明当前预测的是受云污染的半景数据
  out_mask[where(out_mask eq cur_anci, /null)] = 255                ; 如果预测数据已经存在于辅助数据中，去掉这景辅助数据

  left_need = left_mask[UNIQ(left_mask, SORT(left_mask))]           ; 挑选矩阵中的惟一值  
  if max(left_need) eq 255 then begin
    left_need = left_need[where(left_need ne 255, /null)]           ; 去掉最大值255，得到左边需要的辅助数据编号
  endif

  right_need = right_mask[UNIQ(right_mask, SORT(right_mask))]       ; 挑选矩阵中的惟一值
  if max(right_need) eq 255 then begin
    right_need = right_need[where(right_need ne 255, /null)]        ; 去掉最大值255，得到右边需要的辅助数据编号
  endif
  anci_need = [left_need, right_need]                               ; 左右两边需要的辅助数据的编号
  base_in = where(anci_need eq band_tm, base_in_num, /null) 
  if (base_in_num eq 0) then anci_need = [anci_need, band_tm]       ; 如果参考数据不在挑选的结果中，将参考数据添加上
end
;--------------------------------------------------------------------------------------------------------------------------


;;; weight combination: combine multi-predictions based on weights derived from coarse pixels
; 使用加权的方式来将不同辅助数据的预测结果进行组合，最终得到一幅完整的预测结果
;------------------------------------------------------------------------------------------------
pro GetCombined_mean, partial_pre, partial_mask, anci_need, predict, modis_name, loc, fine0=fine0, dlt_loc, dt_mask, ratio

  m_size = size(partial_mask)
  ns = m_size[1]                                             ; 细尺度像元的列数
  nl = m_size[2]                                             ; 细尺度像元的行数
  tns = dlt_loc[0]+ns
  tnl = dlt_loc[1]+nl
  fine0 = partial_pre[*, *, 0]
  pre_num = n_elements(anci_need)

  ; 读取预测时刻的MODIS数据
  FileName2 = modis_name+strtrim(predict,1)                   ; 打开第一幅MODIS的NDVI数据
  GetData,ImgData=coarse0, FileName = FileName2, Fid = Fid2, ns = ns_coarse, nl = nl_coarse, data_type = data_type
  if data_type eq 2 then coarse0 = coarse0/10000.0
  coarse0=FLOAT(coarse0)
  envi_file_mng, id=fid2, /remove
  ; 利用shift函数计算9*9的窗口，在窗口内计算权重
  coarse0_win = fltarr(ns_coarse, nl_coarse, 9)               ; 存储3*3窗口的像元
  coarse0_win[*, *, 0] = coarse0
  coarse0_win[*, *, 1] = shift(coarse0, 1, 1)
  coarse0_win[*, *, 2] = shift(coarse0, 1, 0)
  coarse0_win[*, *, 3] = shift(coarse0, 1, -1)
  coarse0_win[*, *, 4] = shift(coarse0, 0, 1)
  coarse0_win[*, *, 5] = shift(coarse0, 0, -1)
  coarse0_win[*, *, 6] = shift(coarse0, -1, 1)
  coarse0_win[*, *, 7] = shift(coarse0, -1, 0)
  coarse0_win[*, *, 8] = shift(coarse0, -1, -1)
  dist_mask = fltarr(ns, nl, pre_num)                         ; 用来存储权重，权重与mask范围保持一致
  
  for pi=0, pre_num-1 do begin

    pre_anci = anci_need[pi]
    tmp_mask = partial_mask[*, *, pi] eq 0                     ; 将0的值转换成1, 1的值转换成0        
    FileName2 = modis_name+strtrim(pre_anci,1)                 ; 打开第一幅Modis 的NDVI数据
    GetData,ImgData=coarse1, FileName = FileName2,Fid = Fid2, data_type = data_type
    if data_type eq 2 then coarse1 = coarse1/10000.0
    coarse1=FLOAT(coarse1)
    envi_file_mng, id=fid2, /remove
    coarse1_win = fltarr(ns_coarse, nl_coarse, 9)
    coarse1_win[*, *, 0] = coarse1
    coarse1_win[*, *, 1] = shift(coarse1, 1, 1)
    coarse1_win[*, *, 2] = shift(coarse1, 1, 0)
    coarse1_win[*, *, 3] = shift(coarse1, 1, -1)
    coarse1_win[*, *, 4] = shift(coarse1, 0, 1)
    coarse1_win[*, *, 5] = shift(coarse1, 0, -1)
    coarse1_win[*, *, 6] = shift(coarse1, -1, 1)
    coarse1_win[*, *, 7] = shift(coarse1, -1, 0)
    coarse1_win[*, *, 8] = shift(coarse1, -1, -1)
    dist_MOD = 1/total(abs(coarse0_win - coarse1_win), 3)       ; 计算出了距离矩阵，实际上是计算的AD值
;    dist_MOD = 1/mean(abs(coarse0_win - coarse1_win)^2, DIMENSION=3)       ; 计算出了距离矩阵，实际上是计算的RMSE值
    dist_fine = fltarr(ns, nl)                                  ; 细尺度下的距离矩阵
    for cj=0, nl_coarse-1 do begin
      for ci=0, ns_coarse-1 do begin

        pixel_c = partial_mask[max([0,ratio*ci-dlt_loc[0]]):min([ratio*(ci+1)-1,tns-1]-dlt_loc[0]), max([0,ratio*cj-dlt_loc[1]]):min([ratio*(cj+1)-1,tnl-1]-dlt_loc[1]), pi]
        if total(pixel_c) lt ratio^2 then begin                 ; there are fine pixels within this coarse pixel
          dist_fine[max([0,ratio*ci-dlt_loc[0]]):min([ratio*(ci+1)-1,tns-1]-dlt_loc[0]), max([0,ratio*cj-dlt_loc[1]]):min([ratio*(cj+1)-1,tnl-1]-dlt_loc[1])]=dist_MOD[ci, cj]
        endif

      endfor
    endfor    
    dist_fine = dist_fine*tmp_mask                              ; 与mask的覆盖范围保持一致的距离矩阵
    dist_mask[*, *, pi] = dist_fine                             ; 将权重矩阵存储到对应的位置

  endfor

  ; 处理方式1：最终生成的加权结果：代码能够保证不会出现死像元，因此可以这样相除
  fine0 = total(partial_pre*dist_mask, 3)/total(dist_mask, 3)
  fine0[where(dt_mask eq 1)] = 0                                ; 控制输出范围最大不能超过研究区域边界
end
;---------------------------------------------------------------------------------


;;; use this function to smooth the prediction results, removing block effects
;-------------------------------------------------------------------------------------------
pro img_smooth, fine1, fine0, classFine1, mask, dtmask, w_neighbor, ratio, num_similar_pixel, DN_min, DN_max

  tmp_dlt = fine0 - fine1
  f_size = size(fine1)
  ns = f_size[1]
  nl = f_size[2]
  for j=0,nl-1,1 do begin
    for i=0,ns-1,1 do begin
      if (mask[i,j] eq 0 ) then begin                                ; 根据mask文件来确定背景，mask取值为0就是数值，不是背景

        tw_neighbor = w_neighbor
        ai=max([0,i-tw_neighbor*ratio])
        bi=min([ns-1,i+tw_neighbor*ratio])                           ; 列
        aj=max([0,j-tw_neighbor*ratio])                              ; 行
        bj=min([nl-1,j+tw_neighbor*ratio])

        ci=i-ai                                                      ; 当前像元在滑动窗口中的位置
        cj=j-aj
        class_t=classFine1[i,j]                                      ; 当前类型值
        numRow=bj-aj+1                                               ; 滑动窗口的大小
        numCol=bi-ai+1                                               ; 滑动窗口的大小
        class_Row=[indgen(numRow,1)+1] ## [intarr(numRow)+1]
        class_Col=[intarr(numRow,1)+1] ## [indgen(numRow)+1]
        D_D_cand=((class_Row-1-cj)^2+(class_Col-1-ci)^2)^0.5 + 0.00000001

        ;select same-class pixels
        class_sel_Win = classFine1[ai:bi,aj:bj]
        tmp_dlt_win = tmp_dlt[ai:bi,aj:bj]
        mask_sel_Win = dtmask[ai:bi,aj:bj] eq 0                        ; 将局部mask进行反转，0表示无值，1表示有值
        class_sel_Win = class_sel_Win*mask_sel_Win
        ind_same_class=where(class_sel_Win eq class_t, num_sameclass)          ; 挑选同类型像元，函数返回的是下标

        ;searching similar pixels from same-class in known fine image
        similarity_sp=fltarr(num_sameclass)
        similarity_sp=abs((fine1[ai:bi,aj:bj])[ind_same_class]-fine1[i,j])               ; 计算相似性
        order_dis=sort(similarity_sp)                                                    ; 对相似性进行排序，从小到大;函数输出的是下标
        number_cand=min([num_similar_pixel,num_sameclass])                               ; 在这里控制相似像元的数目
        indcand=order_dis[0:number_cand-1]                                               ; select the N most similar samples
        D_D_cand=D_D_cand[ind_same_class[indcand]]                                       ; compute weight for each simialr pixel

        ;normalize these distances
        D_D_cand=1.0 + D_D_cand/(tw_neighbor*ratio)                                      ; 这里的用不用理论上来说是没有影响的
        C_D=1.0/D_D_cand
        weight=C_D/total(C_D)

        ;predict the value
        change_cand=tmp_dlt_win[ind_same_class[indcand]]                      ; 挑选出来的相似像元的增量值
        fine0[i,j]=fine1[i,j]+total(weight*change_cand)                                  ; 对增量进行加权平滑

        ; 这里还需要再添加一步异常值检测的处理
        if (fine0[i,j] lt DN_min) then begin
          fine0[i,j] = max([fine1[i,j] + mean(tmp_dlt_win[ind_same_class]), -0.99])
        endif
        if (fine0[i,j] gt DN_max) then begin
          fine0[i,j] = min([fine1[i,j] + mean(tmp_dlt_win[ind_same_class]), 0.99])
        endif

      endif
    endfor
  endfor
end
;------------------------------------------------------------------------------------------------



;;; single date prediction based on temporal increment and spatial increment
; 为了能够处理不规则的影像，这里增加了一个偏移量作为新的参数来进行处理dlt_loc
;-------------------------------------------------------------------------------------------------------------------
pro HysNDVI, fine1, coarse1, c1_tps, coarse0, c0_tps, classFine1, het_index, w_neighbor, mask, fine0=fine0, DN_min, DN_max, dlt_loc, ratio
  ; mask文件用来对插值过程进行控制，只需要处理mask的像元值为0的位置，1表示受污染区域

  ;open the fine image of the first pair
  fine1=float(fine1)
  f_size = size(fine1)
  ns = f_size[1]
  nl = f_size[2]

  ;open the coarse image of the first pair
  c_size = size(coarse1)
  ns_coarse = c_size[1]
  nl_coarse = c_size[2]

  ; 根据TPS插值结果计算空间增量,tps的结果必须先裁剪好
  dlt_c = c0_tps - c1_tps                                                      ; 空间增量计算一次即可
  tns = ns + dlt_loc[0]                                                        ; dlt_loc = [loc_y, loc_x]
  tnl = nl + dlt_loc[1]
  
  ;;; 将t0时间的插值换成F0
  dlt_c = dlt_c[dlt_loc[0]:tns-1, dlt_loc[1]:tnl-1]                            ; 截取有效的空间插值部分
  dlt_c = dlt_c*(mask eq 0)                                                    ; dlt_c中并非所有的值都需要保留下来，只留下mask需要的值即可
  K_coarse=(coarse0-coarse1)                                                   ; MODIS尺度的增量

  ;;; 1、调用函数计算混合像元分解得到的局部窗口内的增量,frac_coarse和K_fine是待求结果。求解时按照最大覆盖范围的MODIS数据计算
  ;---------------------------------------------------------------------------
  GetTemp_Increment, frac_coarse=frac_coarse, K_fine=K_fine, K_coarse, dlt_loc, $
    w_neighbor, classFine1, fine1, het_index, mask_shrink_c = mask_shrink_c, ratio

  ; 在局部区域内对混合像元分解的增量进行修正，避免由于MODIS的问题导致局部解混出错
  s_K_fine = size(K_fine)                                                      ; the filling value of K_fine is 255
  for i=0, s_K_fine[3]-1 do begin                                              ; 分类型进行检查
    K_fine_cls = K_fine[*, *, i]
    effect_incre = K_fine_cls[where(K_fine_cls ne 255.0, /null)]                 ; 这一步挑选已经去除了背景
    mean_incre = mean(effect_incre)                                            ; 整景影像总体而言，增量的均值,这里出现0会有问题
    std_incre = stddev(effect_incre)
    idx_abnormal = where(abs(K_fine_cls - mean_incre) gt 2*std_incre and K_fine_cls ne 255.0, /null, cnt_abnormal)        ; 使用两倍的标准差进行异常值剔除
    if cnt_abnormal gt 0 then begin
      K_fine_cls[idx_abnormal] = mean_incre
      K_fine[*, *, i] = K_fine_cls[*, *]
    endif
  endfor
  K_fine_cls = 0  
  
  ; 计算每个细像元尺度上混合像元分解的增量
  dlt_t = fltarr(ns, nl)                                                        ; 混合像元分解的增量
  for j=0,nl-1,1 do begin
    tj = j+dlt_loc[1]                                                           ; Landsat在MODIS上的位置
    
    for i=0,ns-1,1 do begin                                                     ; 根据像元在TM上的位置，确定MODIS的窗口      
      if mask[i, j] eq 0 then begin                                             ; 如果是数据才需要进行后续处理
        
        ti = i+dlt_loc[0]
        classID=classFine1[i,j]-1                                               ; 分类值
        if (mask_shrink_c[ti/ratio, tj/ratio] eq 1) then begin
          tw_neighbor = w_neighbor + 2                                          ; 说明i，j位置是边缘，扩大窗口
        endif else begin
          tw_neighbor = w_neighbor
        endelse

        K_fine_win=$
          K_fine[max([0,ti/ratio-tw_neighbor]):min([(tns-1)/ratio,ti/ratio+tw_neighbor]), max([0,tj/ratio-tw_neighbor]):min([(tnl-1)/ratio,tj/ratio+tw_neighbor]), classID]
        K_idx=where(K_fine_win ne 255.0, numK_idx)                              ; 只将不为255的变化量那个提取出来，在混合像元分解时填充值设置成了255
        
        if numK_idx gt 3 then begin
          dlt_t[i,j] = mean(K_fine_win[K_idx])                                  ; 使用的是平均的增量值,得到的是Landsat真实mask的对应的像元值
        endif else begin                                                        ; 下面这个是为了，将窗口再扩大
          K_fine_win=$
            K_fine[max([0,ti/ratio-tw_neighbor-1]):min([(tns-1)/ratio,ti/ratio+tw_neighbor+1]), max([0,tj/ratio-tw_neighbor-1]):min([(tnl-1)/ratio,tj/ratio+tw_neighbor+1]), classID]
          K_idx=where(K_fine_win ne 255.0, numK_idx)
          dlt_t[i,j] = mean(K_fine_win[K_idx])          
        endelse

      endif
    endfor
  endfor  
  dlt_t = dlt_t*(mask eq 0)                                                     ; dlt_t中可能出现局部区域解混出现错误的情况


  ;;; 2、根据前面计算得到的空间增量与时间增量，调用CLS方法计算的权重进行组合
  ;-----------------------------------------------------------------
  tK_coarse = K_coarse*(mask_shrink_c eq 0)                                     ; 有效的值才需要参与评估
  dlt_comb = CLS_main(tK_coarse, dlt_c, dlt_t, dlt_loc, ratio, mask_shrink_c)   ; use constrained least square
  Temproal_Ensemble = fine1 + dlt_comb


  ;;; 3、进行误差补偿，有多种方式，平均值，Zhu，纯净度
  ;------------------------------------------------------------------
  fine0 = Temproal_Ensemble                                                         ; 给fine0赋一个初值，避免漏值
  iter = 1
  if iter eq 0 then begin
    for j=0,nl_coarse-1,1 do begin
      for i=0,ns_coarse-1,1 do begin
        ; 只能使用缩小范围的MODIS尺度的mask文件
        if mask_shrink_c[i, j] eq 0 then begin                                                                  ; MODIS尺度的mask，非0表示是数值
          ; 将t1和预测得到的t0时刻的细尺度像元取出
          tmp_fine1 = fine1[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]    ; T1时刻的真实值
          tmp_fine0 = Temproal_Ensemble[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]
          dlt = mean(tmp_fine1) - mean(tmp_fine0) - (coarse1[i,j]-coarse0[i,j])                                 ; 这个增量对于边界区域可能是有问题的
          fine0[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]=$
            tmp_fine0+dlt                                                                                       ; average residuals
        endif
      endfor
    endfor
  endif else begin
    for j=0,nl_coarse-1,1 do begin
      for i=0,ns_coarse-1,1 do begin        
        ; 只能使用缩小范围的MODIS尺度的mask文件
        if mask_shrink_c[i, j] eq 0 then begin                                                                  ; MODIS尺度的mask，非0表示是数值
          ; 将t1和预测得到的t0时刻的细尺度像元取出
          tmp_fine1 = fine1[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]
          tmp_fine0 = Temproal_Ensemble[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]
          tmp_hidx = het_index[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]  ; 像元纯净度
          w_change_tps=(c0_tps[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]])-tmp_fine0
          num_pix = n_elements(tmp_fine1)
          dlt = (coarse0[i,j]-coarse1[i,j]) - (mean(tmp_fine0)-mean(tmp_fine1))                           ; 这个增量对于边界区域可能是有问题的
            
          if (dlt le 0) then begin                                                ; tmp_fine0的估计总体偏大
            ind_noc=where(w_change_tps gt 0, num_noc)                             ; fine0的估计已经偏大了，TPS的结果如果也偏大，就肯定不能用
            if (num_noc gt 0) then begin
              w_change_tps[ind_noc]=0                                             ; 如果待分配的残差是负数，就将所有的TM影像上的正残差强行赋值为0
            endif
          endif else begin
            ind_noc=where(w_change_tps lt 0, num_noc)
            if (num_noc gt 0) then begin
              w_change_tps[ind_noc]=0                                             ; 如果待分配的残值取值为正，将TM上所有的取值为负的残差强行赋值为0
            endif
          endelse

          w_change_tps=abs(w_change_tps)                                          ; 将负的全都转成正的
          w_unform=w_change_tps
          w_unform[*]=abs(dlt)                                                    ; 也转换成为了正数
          w_change=w_change_tps*tmp_hidx+w_unform*(1.0-tmp_hidx)+0.00000001
          w_change=w_change/(mean(w_change)) ;nomalize weight                     ; 这样计算得到最终权重的和等于m，而不等于1
          fine0[max([0,ratio*i-dlt_loc[0]]):min([ratio*(i+1)-1,tns-1])-dlt_loc[0], max([0,ratio*j-dlt_loc[1]]):min([ratio*(j+1)-1,tnl-1])-dlt_loc[1]]=tmp_fine0+w_change*dlt
        endif
      endfor
    endfor
  endelse
  Temproal_Ensemble = 0    
  

  ;;; 对融合结果进行再一次噪声去除，发现融合结果中有一些异常值。
  ;修改：2017年3月27日，使用局部窗口中的2倍标准差的方式来进行异常值检测
  ;-----------------------------------------------------------------------------------------------------------
  ; 2018.12.31 use interger to avoid being lack of memory
  fine0_win = intarr(ns, nl, 9)           
  fine0_win[*, *, 0] = fix(fine0*10000)
  fine0_win[*, *, 1] = fix(shift(fine0, 1, 1)*10000)
  fine0_win[*, *, 2] = fix(shift(fine0, 1, 0)*10000)
  fine0_win[*, *, 3] = fix(shift(fine0, 1, -1)*10000)
  fine0_win[*, *, 4] = fix(shift(fine0, 0, 1)*10000)
  fine0_win[*, *, 5] = fix(shift(fine0, 0, -1)*10000)
  fine0_win[*, *, 6] = fix(shift(fine0, -1, 1)*10000)
  fine0_win[*, *, 7] = fix(shift(fine0, -1, 0)*10000)
  fine0_win[*, *, 8] = fix(shift(fine0, -1, -1)*10000)
  fine0_std = stddev(fine0_win, DIMENSION = 3)/10000.0
  fine0_win = 0                                                    ; 释放内存
  fine0_std = fine0_std*(mask eq 0)                                ; 对标准差进行掩膜处理  
  fine0_filter = mean_filter(fine0, 3)  
  fine0_filter = fine0_filter*(mask eq 0)
  fine0_dlt = abs(fine0 - fine0_filter)
  ind_f = where(fine0_dlt gt 2*fine0_std, cunt)                    ; 在这里设置阈值，对异常值进行检测
  fine0_dlt = 0                                                    ; 释放内存
  if (cunt gt 0) then begin                                           ; 存在潜在的异常值
    ai = ind_f/ns                                                     ; 行数
    aj = ind_f mod ns                                                 ; 列数
    for i=0, cunt-1 do begin
      patch_class = classFine1[max([aj[i]-2, 0]):min([aj[i]+2, ns-1]), max([ai[i]-2, 0]):min([ai[i]+2, nl-1])]
      patch_mask = mask[max([aj[i]-2, 0]):min([aj[i]+2, ns-1]), max([ai[i]-2, 0]):min([ai[i]+2, nl-1])]
      patch_class = patch_class*(patch_mask eq 0)
      patch_class[2, 2] = 0                                           ; 将中心像元置空，避免中心像元的值影响最后的结果
      patch_val = fine0[max([aj[i]-2, 0]):min([aj[i]+2, ns-1]), max([ai[i]-2, 0]):min([ai[i]+2, nl-1])]
      class_id = classFine1[aj[i], ai[i]]                              ; class id of current fine pixel
      idx = where(patch_class eq class_id, nsame, /null)
      if nsame lt 2 then begin                                        ; 同类像元太少，表明可能是孤立的异常值
        fine0[ind_f[i]] = fine0_filter[ind_f[i]]
      endif else begin                                                ; 有足够数量的同类像元，可以进行取平均处理
        fine0[ind_f[i]] = mean(patch_val[idx])
      endelse       
    endfor
  endif
  fine0_filter = 0                                                  ; 释放内存

  ; 最后一步，将图像按照mask的范围输出，超出范围的像元都去掉
  fine0 = fine0*(mask eq 0)
end
;-----------------------------------------------------------------------------------------------------------------


;;; prediction by each fine image within the subset produced by "GetAncillaryNDVI" function
;-----------------------------------------------------------------------------------------------------------------
pro loop_HysNDVIbyAncillary, base, predict, NEWPATH, modis_name, classFine1, het_index, $
  w_neighbor, mask, fine1, fine0=fine0, DN_min, DN_max, dlt_loc, ratio
  ;fine1与base是相对应的，不需要重新读取

  tmp_fine1 = fine1                                            ; 给一个初始值
  fine0 = fine1
  if base lt predict then begin                                ; 预测时间在后，需要进行正向预测
    iter = 1
;    iter = predict - base
  endif else begin                                             ; 预测时间在前，需要进行逆向预测
    iter = -1
;    iter = predict - base
  endelse
  next_pre = base                                              ; 表示下一景预测时间
  
  while next_pre ne predict do begin                           ; 直到预测时间到达才停止    
    current = next_pre                                         ; T1时刻的波段编号
    next_pre = next_pre + iter                                 ; T0时刻的波段编号
    ; 在开始进行预测时，需要先判断前一景预测结果是否已经存在了，需要对已有的临时预测文件进行搜索
    next_file = NEWPATH + '\base_TM_' + strtrim(base,1) + '_' + strtrim(next_pre,1)
    if FILE_TEST(next_file) eq 1 then begin                    ; 如果文件已经存在就需要直接先读取出来
      GetData, ImgData=r_tmp_fine1, FileName = next_file, Fid = Fid1
      r_tmp_fine1 = float(r_tmp_fine1)
      envi_file_mng, id=Fid1, /remove
      tmp_mask = r_tmp_fine1 eq 255.0                          ; 255是填充的背景值
      if total(abs(tmp_mask - mask)) eq 0 then begin           ; 表示两个mask文件是完全一样的,可以不用预测了
        reverse_mask = mask eq 0                               ; mask等于0的标记为1，否则标记为0
        tmp_fine1 = r_tmp_fine1*reverse_mask                   ; 得到消除了255背景以及扩展MODISmask的结果，直接读取出已有的结果
        tmp_fine0 = tmp_fine1                                  ; 如果后续的预测一直不进行，就直接将读取的结果输出
        continue                                               ; 结束本次迭代
      endif
    endif
    
    ; 读取T1时刻的MODIS数据
    FileName2 = modis_name+strtrim(current,1)                ; 打开第一幅Modis 的NDVI数据
    GetData,ImgData=coarse1, FileName = FileName2,Fid = Fid2, Data_Type = Data_Type
    if Data_Type eq 2 then coarse1 = coarse1/10000.0
    coarse1=FLOAT(coarse1)
    envi_file_mng, id=fid2, /remove

    ; 读取T1时刻的TPS插值结果
    GetData, ImgData=c1_tps, FileName = FileName2 + '_TPS', Fid = Fid2, Data_Type = Data_Type
    if Data_Type eq 2 then c1_tps = c1_tps/10000.0
    c1_tps = FLOAT(c1_tps)
    envi_file_mng, id=fid2, /remove

    ; 读取预测时刻的MODIS数据
    FileName3 = modis_name+strtrim(next_pre,1)                ; 打开第二幅MODIS的NDVI数据
    GetData,ImgData=coarse0, FileName = FileName3, Fid = Fid3, Data_Type = Data_Type
    if Data_Type eq 2 then coarse0 = coarse0/10000.0
    coarse0=FLOAT(coarse0)
    envi_file_mng, id=fid3, /remove

    ; 读取对应的TPS插值结果
    GetData, ImgData=c0_tps, FileName = FileName3 + '_TPS', Fid = Fid3, Data_Type = Data_Type
    if Data_Type eq 2 then c0_tps = c0_tps/10000.0
    c0_tps = FLOAT(c0_tps)
    envi_file_mng, id=fid3, /remove

    ; 调用函数进行融合
    HysNDVI, tmp_fine1, coarse1, c1_tps, coarse0, c0_tps, classFine1, het_index, w_neighbor, mask, fine0=tmp_fine0, DN_min, DN_max, dlt_loc, ratio
    tmp_fine1 = tmp_fine0                                    ; 预测的基准值要开始发生改变了
    out_tmp_fine0 = tmp_fine0 + mask*255.0                   ; 在输出之前给数据添加上255的背景值
    Envi_Write_Envi_File, out_tmp_fine0, Out_Name = next_file, r_fid = fid_temp                    ; 输出的是浮点型的数值
    envi_file_mng, id=fid_temp, /remove
  endwhile

  fine0 = tmp_fine0                                            ; 将需要的那景预测结果输出
end
;------------------------------------------------------------------------------------



;;; remove noise in the time series data 
;-------------------------------------------------------------------------------------
pro smooth_NDVIts, filename

  GetData,ImgData=ndvi_ts, ns = ns,nl = nl, FileName = filename, Fid = Fid, Data_Type = Data_Type, map_info = map_info
  envi_file_mng, id = Fid, /remove
  for i=0, nl-1 do begin
    for j=0, ns-1 do begin
      vector_in = ndvi_ts[j, i, *]
      vector_in = reform(vector_in, 23)/10000.0
      savgolFilter = SAVGOL(4, 4, 0, 2)
      rst = CONVOL(vector_in, savgolFilter, /EDGE_TRUNCATE)
      dlt = abs(vector_in - rst)      
      idx = where(dlt lt 0.015, /null, cnt)
;      std = stddev(dlt)
      ;      idx = where(dlt lt mean(dlt)+2*std, /null, cnt)
      rst[idx] = vector_in[idx]
      ndvi_ts[j, i, *] = fix(rst*10000)
    endfor
  endfor  
  
  Envi_Write_Envi_File, ndvi_ts, Out_Name = filename, r_fid = fid_temp, $
    map_info = map_info, ns = ns, nl = nl, out_dt = Data_Type
    
end
;-----------------------------------------------------------------------------------


;-------------------------------------------------------------------------------------
;
;                                main program
; remember: all input images must have correct and consistent geographic coordinates
;-------------------------------------------------------------------------------------

pro  downscaling_IFSDAF_main
  envi, /restore_base_save_files
  envi_batch_init

  t0=systime(1)             ; the initial time of program running

  ; please set the following parameters: all image files should be ENVI standard format - *.hdr
  ; -----------------------------------------------------------------------
  ;; set the time gap of MODIS NDVI, 8(8day composite), 16(16day composite) or 1(daily)
  t_deta = 16 
  ;; current filepath for all image data, please save all data under this path: never use underscore in this path
  Upper_path =  'D:\subset_Shennongjia'
  ;; stacking of MODIS NDVI by the sequense of time
  FileName_MODIS = Upper_path + '\MOD_16Day_ndvi' 
  ;; Landsat NDVI on T0: for example LC81260382015287LGN00_ndvi
  FileName_TM = Upper_path + '\LC81260382015287LGN00_ndvi'  
  ;; Classified map based on clear fine images
  ClassImg = Upper_path + '\landcover'  
  ;; path for partly cloud contaminated fine NDVI along with corresponding Fmasks
  ancillaryFolder = Upper_path + '\Cloud'   ; leave this folder empty when not using partly contaminated images
  ;; name of the final fusion results
  result_name=file_basename(FileName_TM)
  loc = strpos(result_name, '_')
  result_name=Upper_path + '\'+ strmid(result_name, 0, loc)+'_IFSDAF' 
  
  
  ;; you can just leave the following parameters unchanged
  ;----------------------------------------------------------------------
  w_neighbor = 1            ; half window size for unmixing, the whole window is 3*3  
  DN_min = -1               ; the minimum value for NDVI
  DN_max = 1                ; the maximum value for NDVI
  win = 5                   ; window size for TPS interpolation, never change this value
  folder = 'band COUS'      ; folder to store tempory results, this folder will be created within your current data path
  time_range = 60           ; time range for selecting ancillary fine images when prediction current date
  ;------------------------------------------------------------------------

  if FILE_TEST(Upper_path+'\'+folder, /DIRECTORY) eq 0 then begin ; 文件夹不存在时，进行创建
    FILE_MKDIR, Upper_path+'\'+folder, /NOEXPAND_PATH             ; 在当前文件夹中创建一个新的子文件夹
  endif
  NEWPATH = Upper_path+'\'+folder                                 ; filepath for results. 子文件夹的路径
  tm_name = NEWPATH+'\band_tm'
  modis_name = NEWPATH+'\band_mod'                                ; 就是Modis数据在文件夹中的存储名称
  out_day = GetDay(FileName_TM, t_deta)        ; get the exact DOY of clear NDVI
  ; out_day = GetDay2(FileName_TM, t_deta, year, month, day)        ; get the exact DOY of clear NDVI
  band_tm = out_day[0]                                            ; the sequence of clear NDVI
  day_tm = out_day[1]

  ;; open MODIS NDVI stacking file
  ;---------------------------------------------------------------------------------------------------------------
  envi_open_file,FileName_MODIS,r_fid=fid
  envi_file_query,fid,ns=ns_coarse,nl=nl_coarse,nb=nb_coarse,dims=dims_coarse,Data_Type = Data_Type
  map_info_modis = envi_get_map_info(fid=fid)
  proj = map_info_modis.proj
  pix_MOD = map_info_modis.ps[0]                                   ; pixel size of MODIS pixel
  
  ;; open Landsat NDVI file
  ;---------------------------------------------------------------------------------------------------------------
  GetData,ImgData=base_fine1,ns = ns,nl = nl,Data_Type = tData_Type,FileName = FileName_TM, dims=dims, Fid = fid_TM  
  map_info = envi_get_map_info(fid=fid_TM)
  proj=envi_get_projection(fid=fid_TM,pixel_size=pixel_size_result)                    ; 实际上，proj在map_info里面是已经包含的
  ratio = pix_MOD/pixel_size_result[0]                             ; the scale difference between MODIS data and Landsat data
  envi_file_mng, id=fid_TM, /remove                                                    ; remove the base fine image
  Envi_Write_Envi_File, base_fine1, Out_Name = tm_name+strtrim(band_tm,1), map_info = map_info, r_fid = f_fid_base
  n_bands=nb_coarse                                                          ; 用一个新变量将波段数目存储下来
  sfids=LON64ARR(n_bands)                                                    ; 需要将文件id设置长一些，使用长整型
  sfids[band_tm-1] = f_fid_base 
  orig_ns=ns
  orig_nl=nl
  base_fine1 = 0                                       ; clear the storage
  
  ;; 根据Landsat与MODIS的投影信息计算出Landsat在MODIS数据中所处的位置: MC[0]x pixel location, MC[1]y pixel location, MC[2]x map locatoin, MC[3]y map location
  if proj.NAME ne 'Arbitrary' then begin
    loc_x = abs((map_info_modis.mc[2] - map_info_modis.ps[0]*map_info_modis.mc[0])-(map_info.mc[2] - map_info.ps[0]*map_info.mc[0]))/map_info.ps[0]  ; nl
    loc_y = abs((map_info_modis.mc[3] + map_info_modis.ps[1]*map_info_modis.mc[1])-(map_info.mc[3] + map_info.ps[1]*map_info.mc[1]))/map_info.ps[1]  ; ns
    dlt_loc = round([loc_x, loc_y])                                                         ; 以四舍五入的方式定义偏移量    
  endif else begin
    dlt_loc = [0, 0]                                                                        ; 没有投影就直接定义为0偏移
    if total(pixel_size_result) eq 0 then pixel_size_result = pixel_size_result+1
  endelse    
    
  pos=indgen(nb_coarse)
  for iband=0,nb_coarse-1,1 do begin                                         ; 将stack在一起影像文件单独抽取成独立的影像文件
    coarse1 = Envi_Get_Data(Fid = fid, dims = dims_coarse, pos = iband)      ; get MODIS bandi iband
    Envi_Write_Envi_File, coarse1, Out_Name = modis_name+strtrim(iband+1,1), r_fid=r_fid0, map_info = map_info_modis, out_dt = Data_Type
    if FILE_TEST(modis_name+strtrim(iband+1,1)+'_TPS') eq 0 then begin       ; 表明TPS插值文件不存在
      c1_tps = moveWinsub_TPS(coarse1, ratio, win)
      if Data_Type eq 2 then c1_tps = fix(c1_tps)                            ; integer number
      Envi_Write_Envi_File, c1_tps, Out_Name = modis_name+strtrim(iband+1,1)+'_TPS', r_fid=Fid3, out_dt = Data_Type
      envi_file_mng, id=fid3, /remove
    endif
    envi_file_mng, id=r_fid0, /remove                                        ; Modis叠置的NDVI数据分解开后并没有删除掉，后续还会使用到的
  endfor
  coarse1 = 0                                                                ; 去除内存占用
  c1_tps = 0                                                                 ; 去除内存占用
  envi_file_mng, id=fid, /remove                                             ; 将叠置在一起的Modis数据释放掉，避免内存占用的问题


  ;; genereate the prediction order
  ;-------------------------------------------------------------------------------------
  predict=intarr(n_bands-1)                            ; 存储的都是编号
  base=intarr(n_bands-1)                               ; 存储的都是编号
  for i=0,n_bands-band_tm-1,1 do begin
    if i lt 0 then break
    predict[i]=i+band_tm+1                             ; 这个意思实际上是以当前TM影像位置为基础，向后进行预测(正向预测)
    base[i]=predict[i]-1
  endfor
  for j=n_bands-band_tm,n_bands-2,1 do begin           ; 这里的配比非常地具有艺术性
    if j lt 0 then break
    base[j]=n_bands-j
    predict[j]=base[j]-1                               ; 以当前位置为基础进行前向预测（逆向预测）
  endfor
  
  ;; open classified image
  ;--------------------------------------------------------------------------------------------------------------------------------------
  envi_open_file,classImg,r_fid=fid_cls
  classFine1 = envi_get_data(fid=fid_cls,dims=dims,pos=0)      
  cls_id = where(classFine1[uniq(classFine1, sort(classFine1))] ne 0, n_N, /null)   ; n_N is the number of classes  
  dt_mask = classFine1 eq 0                                                         ; 类型等于0的就是背景，设置取值为1，其他取值为0
  max_class=max(classFine1)                                                         ; the maximum number of classes
  if max_class gt n_N then begin                                                    ; there are null classes
    i_new_c=1
    L1_class=intarr(ns,nl)
    for iclass=1, max_class,1 do begin
      ind_ic=where(classFine1 eq iclass, num_ic, /null)
      if (num_ic gt 0) then begin
        L1_class[ind_ic]=i_new_c
        i_new_c=i_new_c+1
      endif
    endfor
    classFine1 = L1_class                                                             ; change the classfied map  
  endif  
  
  ; calculate pixel purity
  if file_test(newpath+'\het_fine') eq 0 then begin
    GetHindex, Class=ClassFine1, Ratio=1.0*ratio, Hindex=het_index                   ; 计算像元纯净度
    Envi_Write_Envi_File, het_index, Out_Name = newpath+'\het_fine', r_fid=fid_het, map_info = map_info
  endif else begin
    GetData, ImgData=het_index, FileName = newpath+'\het_fine', Fid = fid_het
    envi_file_mng, id=fid_het, /remove
  endelse  
  het_index = het_index * (dt_mask eq 0)
  envi_file_mng, id=fid_cls, /remove
  
  ; read partly contaminated data, including NDVI and Fmasks
  if FILE_TEST(ancillaryFolder, /DIRECTORY) eq 0 then begin                        ; 如果连文件夹都没有，表明没有辅助数据
    cnt = 0
  endif else begin
    anciFile = file_search(ancillaryFolder + '\*mask*.hdr', count=cnt)            ; 找到所有带后缀名的文件
  endelse

  if cnt eq 0 then begin
    anciFlag = 0                                                                    ; 表示没有辅助数据
  endif else begin
    anciFlag = 1                                                                    ; 表示有辅助数据
    an_num = n_elements(anciFile)                                                   ; 计算辅助数据的数目，一组算一个，mask及对应的NDVI算一个
    an_time = intarr(an_num)                                                        ; 存储辅助数据在一年中的编号：1,2,3,4,5....
    an_day = intarr(an_num)                                                         ; 存储辅助数据在一年中的天数：1/17/33...  
    an_NDVI = fltarr(ns, nl, an_num)                                                ; 存储辅助数据的NDVI值
    an_mask = bytarr(ns, nl, an_num)                                                ; 存储收缩到MODIS尺度的云mask文件，将边缘都去掉了，只保留了MODIS的范围
    an_mask_orig = an_mask                                                          ; 存储直接从云数据中读取的mask文件，其中1及其他表示云，0表示可用数据    
    for ani=0, an_num-1 do begin      
      out_day = GetDay(anciFile[ani], t_deta)                    ; 计算辅助数据是位于哪一景下面
      ;out_day = GetDay2(anciFile[ani], t_deta, year, month, day)                    ; 计算辅助数据是位于哪一景下面
      an_time[ani] = out_day[0]
      an_day[ani] = out_day[1]            
      GetData, ImgData = tmp_mask, FileName = anciFile[ani], Fid = Fid1             ; 读取mask文件
      envi_file_mng, id=Fid1, /remove
      ; Fmask中0是清晰像元，1是水体背景值是255
      tmp_mask = tmp_mask*(dt_mask eq 0)                                            ; dt_mask 为0表示研究区域，为1表示研究区外的背景
      tmp_mask = tmp_mask gt 1                                                      ; 将mask文件转化成为0-1文件，0表示无云值，1表示有云值
      tmp_mask = tmp_mask + dt_mask                                                 ; 最后填充值其实也直接定成1，即当做云处理
      an_mask_orig[*, *, ani] = byte(tmp_mask)                                      ; 这里的mask是原始的云mask文件
      location1 = strpos(anciFile[ani],'_', /REVERSE_SEARCH)                                         ; 文件命名时需要以下划线为间隔
      Landsat_name = strmid(anciFile[ani], 0, location1)
      ndvi_of_mask = file_search(ancillaryFolder, file_basename(Landsat_name) + '*ndvi*.hdr', count=cnt)   ; cnt = 0 or 1
      ; ndvi_of_mask = file_search(ancillaryFolder, '*'+year+month+day+'*ndvi*.hdr', count=cnt)  ;;;!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ;; NDVI 的获取要重新定义
      NDVIFile = ndvi_of_mask[0]                                                    ; NDVI文件名称，要求数据的文件名比mask的文件名多一个‘mask’
      GetData, ImgData=tmp_ndvi, FileName = NDVIFile, Fid = Fid1, Data_Type = Data_Type          ; 读取NDVI文件
      if Data_Type eq 2 then tmp_ndvi = tmp_ndvi/10000.0      
      ; 2017.04.01对输入的NDVI进行异常值检测
      ab_idx = where(abs(tmp_ndvi) gt 1, /null, ab_cnt)
      mid_tmp_ndvi=median(tmp_ndvi, 5)                                               ; 返回的是一个矩阵
      if ab_cnt gt 0 then tmp_ndvi[ab_idx] = mid_tmp_ndvi[ab_idx]                    ; 有可能中值滤波的结果仍然是异常值      
      an_NDVI[*, *, ani] = tmp_ndvi
      envi_file_mng, id=Fid1, /remove                  
    endfor
    an_mask = an_mask_orig
    print, an_time                                                                  ; 将辅助数据所在的波段输出
    tmp_ndvi = 0
  endelse
  
  for iband=0, n_bands-2 do begin                                                   ; 总共23幅Modis影像，但只进行22次融合
    print, 'it is predicting', predict[iband]
;    if predict[iband] ne 8 then continue
;    if file_test(NEWPATH+'\band_tm'+strtrim(predict[iband], 1)) then continue
    if anciFlag eq 1 then begin                                                     ; 有辅助数据才需要进行搜索
      ; get partly contaminated NDVI within time range
      GetAncillaryNDVI, band_tm, day_tm, predict[iband], dt_mask, an_time, an_day, t_deta,$
           an_mask_orig, anci_need=anci_need, out_mask = out_mask, cur_anci = cur_anci, loc = loc, time_range   
      ; 按照左边和右边的区分，对未知时间点进行预测
      pre_num = n_elements(anci_need)
      partial_pre = fltarr(orig_ns, orig_nl, pre_num)                               ; 存储预测得到的半景结果
      partial_mask = bytarr(orig_ns, orig_nl, pre_num)                              ; 将每幅辅助数据对应的mask文件记录下来
      for pi=0, pre_num-1 do begin                                                  ; 每一景辅助数据单独对数据进行预测        
        pre_anci = anci_need[pi]                                                    ; 本次预测可以用作辅助数据的那一景        
        base_idx = where(an_time eq pre_anci, cnt)
        if cnt eq 0 then begin                                                      ; 表明不在辅助数据中，应该使用作为参考的那景Landsat
          base_band = band_tm    
          FileName1 = tm_name+strtrim(band_tm, 1)
          GetData,ImgData=fine1,ns = ns,nl = nl,Data_Type = Data_Type,FileName = FileName1,Fid = Fid1, dims = dims
          if Data_Type eq 2 then fine1 = fine1/10000.0                              ; 如果输入是整型则需要调整
          fine1 = float(fine1)      
          tmp_mask = dt_mask                                                        ; 设定成研究区的大小
          partial_mask[*, *, pi] = tmp_mask                                         ; 存储mask文件
        endif else begin                                                            ; 表示需要的数据位于辅助数据集合中
          base_band = pre_anci
          fine1 = an_NDVI[*, *, base_idx]         
          tmp_mask = an_mask_orig[*, *, base_idx]                                   ; 使用的是更新过的云mask文件
          partial_mask[*, *, pi] = an_mask_orig[*, *, base_idx]                     ; 最后加权用的是原始的云mask文件
        endelse
        ; 调用函数，计算根据对应的辅助数据得到的预测结果
        loop_HysNDVIbyAncillary, base_band, predict[iband], NEWPATH, modis_name, classFine1,$
            het_index, w_neighbor, tmp_mask, fine1, fine0=fine0, DN_min, DN_max, dlt_loc, ratio
        ; 将预测结果保存下来
        partial_pre[*, *, pi] = fine0
      endfor      
      
      ; 加权处理
      if pre_num eq 1 then begin                                                    ; 表明预测只有一景数据，不用进行加权
        fine0 = partial_pre                                                         ; 如果只有一景不需要进行求和     
      endif else begin                                                              ; 表明有多个预测结果，需要加权
        ; 进行加权处理
        GetCombined_mean, partial_pre, partial_mask, anci_need, predict[iband], modis_name, loc, fine0=fine0, dlt_loc, dt_mask, ratio
      endelse
      partial_pre = 0
      partial_mask = 0
      mask = dt_mask
 
      ; 进行拼接，将预测值与无云区域进行拼接
      if cur_anci ne 0 then begin                                                   ; 预测的正好是辅助波段，需要进行校正
        ; 在最后输出之前,将真实值与校正过的估计值拼接起来
        base_idx = where(an_time eq cur_anci)
        haze_fine = an_NDVI[*, *, base_idx]
        fine0[where(an_mask_orig[*, *, base_idx] eq 0)] = haze_fine[where(an_mask_orig[*, *, base_idx] eq 0)]
        mask = an_mask_orig[*, *, base_idx] eq 0
      endif           
         
    endif else begin                                                                ; 没有辅助数据，按照阶段1处理即可
      ;open the fine image of the first pair
      mask = dt_mask                                                        ; 0表示可以使用的像元
      FileName1 = tm_name+strtrim(base[iband],1)                                    ; 这里打开的实际上是作为参考的那幅TM的NDVI数据
      GetData,ImgData=fine1,ns = ns,nl = nl,Data_Type = Data_Type,FileName = FileName1,Fid = Fid1, dims = dims
      if Data_Type eq 2 then fine1 = fine1/10000.0
      fine1=float(fine1)
      ;    envi_file_mng, id=fid1, /remove                                         ; 这一步清除暂时不能进行，后面叠置需要用到这个数据，最后会单独对其进行释放

      ;open the coarse image of the first pair
      FileName2 = modis_name+strtrim(base[iband],1)                                 ; 打开第一幅Modis 的NDVI数据
      GetData,ImgData=coarse1,ns = ns_coarse,nl = nl_coarse, FileName = FileName2,Fid = Fid2, dims = dims_coarse, Data_Type = Data_Type
      if Data_Type eq 2 then coarse1 = coarse1/10000.0
      coarse1=FLOAT(coarse1)
      envi_file_mng, id=fid2, /remove

      ; 打开TPS插值的结果
      GetData, ImgData=c1_tps, FileName = FileName2 + '_TPS', Fid = Fid2, Data_Type = Data_Type
      if Data_Type eq 2 then c1_tps = c1_tps/10000.0
      c1_tps = FLOAT(c1_tps)
      envi_file_mng, id=fid2, /remove

      ;open the coarse image of the prediction time
      FileName3 = modis_name+strtrim(predict[iband],1)                              ; 打开第二幅MODIS的NDVI数据
      GetData,ImgData=coarse0, FileName = FileName3, Fid = Fid3, Data_Type = Data_Type
      if Data_Type eq 2 then coarse0 = coarse0/10000.0
      coarse0=FLOAT(coarse0)
      envi_file_mng, id=fid3, /remove

      ; 打开TPS插值的结果
      GetData, ImgData=c0_tps, FileName = FileName3 + '_TPS', Fid = Fid3, Data_Type = Data_Type
      if Data_Type eq 2 then c0_tps = c0_tps/10000.0
      c0_tps = FLOAT(c0_tps)
      envi_file_mng, id=fid3, /remove

      ; 调用函数进行融合预测处理：增加了一个偏移量dlt_loc作为新的参数
      HysNDVI, fine1, coarse1, c1_tps, coarse0, c0_tps, classFine1, het_index, w_neighbor, mask, fine0=fine0, DN_min, DN_max, dlt_loc, ratio      
    endelse
    
    ;; pixel smoothing by similar pixels
    FileName1 = tm_name+strtrim(base[iband],1)                         ; 这里打开的实际上是作为参考的那幅TM的NDVI数据
    GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName1,Fid = Fid1
    if Data_Type eq 2 then fine1 = fine1/10000.0
    fine1=float(fine1)
    num_similar_pixel = 20                                             ; 仅仅使用前20个相似像元
    tmp_dlt = fine0 - fine1                                            ; 计算出前面处理得到的增量
    ;;; use patches to increase the processing speed
    patch_size = 1000
    ns_patch = ceil(1.0*ns/patch_size)
    nl_patch = ceil(1.0*nl/patch_size)
    for ip = 0, nl_patch-1 do begin
      for jp = 0, ns_patch-1 do begin
        x1 = ip*patch_size                                ; real location of subset on the whole img
        x2 = min([ns, (ip+1)*patch_size])-1
        y1 = jp*patch_size
        y2 = min([nl, (jp+1)*patch_size])-1
        xx1 = max([0, x1-w_neighbor*ratio])               ; a subset a little bigger
        xx2 = min([ns-1, x2+w_neighbor*ratio])
        yy1 = max([0, y1-w_neighbor*ratio])
        yy2 = min([nl-1, y2+w_neighbor*ratio])
        Lx1 = x1-xx1                                      ; relative location on the bigger subset
        Lx2 = Lx1+x2-x1
        Ly1 = y1-yy1
        Ly2 = Ly1+y2-y1

        t_fine1 = fine1[xx1:xx2, yy1:yy2]
        t_fine0 = fine0[xx1:xx2, yy1:yy2]
        t_class = classFine1[xx1:xx2, yy1:yy2]
        t_mask = mask[xx1:xx2, yy1:yy2]
        if min(t_mask) eq 1 then continue                 ; this patch only contains background value, no smooth
        bgm_mask = dt_mask[xx1:xx2, yy1:yy2]
        img_smooth, t_fine1, t_fine0, t_class, t_mask, bgm_mask, w_neighbor, ratio, num_similar_pixel, DN_min, DN_max
        fine0[x1:x2, y1:y2] = t_fine0[Lx1:Lx2, Ly1:Ly2]
      endfor
    endfor
    t_fine1 = 0
    t_fine0 = 0
    t_class = 0
    t_mask = 0
    tmp_dlt = 0
    
    ;; output the file to tempory folder              
    if Data_Type eq 2 then fine0 = fix(10000*fine0)
    Envi_Write_Envi_File, fine0, Out_Name = tm_name + strtrim(predict[iband],1), r_fid = fid_temp, $
      map_info = map_info, ns = ns, nl = nl, out_dt = Data_Type
    sfids[predict[iband]-1] = fid_temp                                    ; 将新生成的文件ID号存储在数组中
  endfor                                                                  ; 结束了最外层的大循环处理  
  
  bands_name=strarr(n_bands)
  for i = 0,n_bands-1 do begin
    ; 以3位字符的形式将数字输出，并且不足三位时，在前面补0
    bands_name[i] = string(1+i*16,format='(I03)') + '(' + string(i+1, format='(I02)') + ')'    
    envi_open_file, modis_name+strtrim(i+1, 2), r_fid= r_fid              ; 删除MODIS数据
    envi_file_mng, id = r_fid, /remove, /delete    
    envi_open_file, modis_name+strtrim(i+1, 2)+'_TPS', r_fid= r_fid       ; 删除TPS插值数据    
    envi_file_mng, id = r_fid, /remove  
  endfor
  
  if anciFlag eq 1 then begin                                             ; 表示有辅助数据
    baseFile = file_search(NEWPATH + '\base_TM_*.hdr', count=cnt)
    for bi=0, cnt-1 do begin
      envi_open_file, baseFile[bi], r_fid= r_fid                          ; 删除作为基准的辅助数据
      envi_file_mng, id = r_fid, /remove
    endfor
  endif
  
  ;;; 将数据输出
  out_name = result_name
  pos_result=lonarr(n_bands)
  dims_result=lonarr(5,n_bands)
  dims_result[0,*]=-1
  dims_result[1,*]=0
  dims_result[2,*]=orig_ns-1
  dims_result[3,*]=0
  dims_result[4,*]=orig_nl-1
  envi_doit, 'envi_layer_stacking_doit', fid=sfids, pos=pos_result, $
    dims=dims_result,out_dt=Data_Type,out_name=out_name,out_bname=bands_name,out_ps=pixel_size_result,out_proj=proj
  
  for i = 0,n_bands-1 do begin
    r_fid = sfids[i]
    envi_file_mng, id = r_fid, /remove                              ; 释放掉Landsat数据占用的内存
  endfor
  
  smooth_NDVIts, out_name                                           ; smooth the result
  
  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
end
