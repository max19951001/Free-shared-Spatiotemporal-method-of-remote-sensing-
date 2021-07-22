pro image_merge_STARFM_v2
  p = 0 ;规定计算的波段数：0-1、1-2、2-3、3-4
  window_size = 11;窗口大小
  class_num = 1;假定分类数
  A = 4;距离因子
;  compile_opt idl2
;  envi, /restore_base_save_files
;  envi_batch_init, log_file='d:\test\batch.txt'
  ;file_name_F1 = dialog_pickfile(title = 'fine resolusion map of the 1st pair',path = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\LC81190402014203_b6_subset_15m.tif')
  ;file_name_C1 = dialog_pickfile(title = 'coarse resolusion map of the 1st pair',path = 'D:\遥感影像\丽水2014洪水\预处理后数据\npp_2014203_i3_prjed_15m.tif')
  ;file_name_F2 = dialog_pickfile(title = 'fine resolusion map of the 2nd pair',path = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\LC81190402014299_b6_subset_15m.tif')
  ;file_name_C2 = dialog_pickfile(title = 'coarse resolusion map of the 2nd pair',path = 'D:\遥感影像\丽水2014洪水\预处理后数据\npp_2014299_i3_prjed_15m.tif')
  ;file_name_F1 = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\LC81190402014203_b6_subset_15m.tif'
  ;file_name_C1 = 'D:\遥感影像\丽水2014洪水\预处理后数据\npp_2014203_i3_prjed_15m.tif'
  ;file_name_F2 = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\LC81190402014299_b6_subset_15m.tif'
  ;file_name_C2 = 'D:\遥感影像\丽水2014洪水\预处理后数据\npp_2014299_i3_prjed_15m.tif'
  
  file_name_F1 = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\2014123_b6_ref.tif'
  file_name_C1 = 'D:\遥感影像\丽水2014洪水\NPPVIIRS初步剪裁\npp_d20140503_i3_ref_utm_30m.tif'
  file_name_F2 = 'D:\遥感影像\丽水2014洪水\Landsat裁剪\2014299_b6_ref.tif'
  file_name_C2 = 'D:\遥感影像\丽水2014洪水\NPPVIIRS初步剪裁\npp_d20141026_i3_ref_utm_30m.tif'
  
  
  envi_open_file,file_name_F1,r_fid = fid_0203_30
  envi_open_file,file_name_C1,r_fid = fid_0203_150
  envi_open_file,file_name_F2,r_fid = fid_0204_30
  envi_open_file,file_name_C2,r_fid = fid_0204_150
  
  envi_file_query, fid_0203_150, ns=ns, nl=nl, nb=nb , dims = dims
  data_0203_150 = float(ENVI_GET_DATA(fid=fid_0203_150, dims=dims, pos=p))
  
  envi_file_query, fid_0203_30, ns=ns, nl=nl, nb=nb , dims = dims
  data_0203_30 = float(ENVI_GET_DATA(fid=fid_0203_30, dims=dims, pos=p))/10000
  
  envi_file_query, fid_0204_150, ns=ns, nl=nl, nb=nb , dims = dims
  data_0204_150 = float(ENVI_GET_DATA(fid=fid_0204_150, dims=dims, pos=p))
  
  envi_file_query, fid_0204_30, ns=ns, nl=nl, nb=nb , dims = dims
  data_0204_30 = float(ENVI_GET_DATA(fid=fid_0204_30, dims=dims, pos=p))/10000
  
  map_info=envi_get_map_info(fid=fid_0203_30)
  map_info_out=map_info
  
  predicted_data = fltarr(ns,nl)
  predicted_judge = intarr(ns, nl)
  
    ;移动窗口半径
  r = (window_size-1)/2
  
  for  sample = r , ns- r - 1 do begin
    for  line = r , nl- r - 1 do begin
      ;取像元
      ;--------------------------------------------------------------------------------------------------------
      data1 = data_0203_150[(sample-r):(sample+r),(line-r):(line+r)]
      data2 = data_0203_30[(sample-r):(sample+r),(line-r):(line+r)]
      data3 = data_0204_150[(sample-r):(sample+r),(line-r):(line+r)]
      result = 0
      judge = 0
      ;筛选像元
      ;--------------------------------------------------------------------------------------------------------
      ;方差
      similar_pixel_measure = stddev(data2)/class_num
      similar_pixel = intarr(window_size,window_size)
      if similar_pixel_measure ne 0 then begin
        similar_pixel[where(abs(data2-data2[r,r]) lt similar_pixel_measure)] = 1
      endif else begin
        similar_pixel = 1
      endelse
      
      similar_pixel_index = where(similar_pixel eq 1)
      if similar_pixel_index[0] ne -1 then begin
      ;   图像之差
          S_lm = abs(data1[similar_pixel_index] - data2[similar_pixel_index]) > 0.00001
          T_mm = abs(data1[similar_pixel_index] - data3[similar_pixel_index]) > 0.00001
      ;   像素到中心像元的距离 ，并转成相对距离
          sub_index = ARRAY_INDICES([window_size,window_size],similar_pixel_index,/dimensions)
          ;distance = transpose(1+((sub_index[0,*] - (window_size-1)/2)^2 + (sub_index[1,*] - (window_size-1)/2)^2)^0.5/A)
          distance = transpose(1+((sub_index[0,*] - r)^2 + (sub_index[1,*] - r)^2)^0.5*(window_size-1)/2)
          weight = (1/S_lm*T_mm*distance)/total(1/S_lm*T_mm*distance)
          result = total(weight*(data3[similar_pixel_index]+data2[similar_pixel_index]-data1[similar_pixel_index]))
          judge = 1
      endif else begin
          result = data3[r,r]+data2[r,r]-data1[r,r]
          judge = 2 
      endelse
      
      ;先不处理云
;            ;如果反射率result不在0-1之间，进行修正
;      if result ge 0 and result le 1 then begin
;        predicted_data[sample,line] = result
;      endif else begin
;        ;这里使用上一时相的反射率替换，也许并不是最好的方式
;        predicted_data[sample,line] = data2[r,r]
;      endelse

      predicted_data[sample,line] = result
      predicted_judge[sample,line] = judge
    endfor
  endfor
  
  print,'done'
  
  output_filename_judge = 'D:\遥感影像\丽水2014洪水\results\judge_starfm_v4.tif'
  openw,lun,output_filename_judge,/get_lun
  writeu,lun,predicted_judge
  free_lun,lun
  envi_setup_head,fname = output_filename_judge , ns = ns , nl = nl , nb = 1, data_type = 1, offset = 0, interleave = 0, map_info=map_info_out,/write,/open
  
  ;output_filename = file_name_F2 + '_p'
  output_filename = 'D:\遥感影像\丽水2014洪水\results\2014299_b6_starfm_v4.tif'
  openw,lun,output_filename,/get_lun
  writeu,lun,predicted_data
  free_lun,lun
  envi_setup_head,fname = output_filename , ns = ns , nl = nl , nb = 1, data_type = 4, offset = 0, interleave = 0, map_info=map_info_out,/write,/open
  plot, data_0204_30 , predicted_data , TITLE='starfm vs. truth' , XTITLE='real Landsat reflectance' , YTITLE='STARFM-fused reflectance' , psym = 3 , xrange = [0,1] , yrange = [0,1]
  plot,[0,1],[0,1],xrange = [0,1] , yrange = [0,1],/noerase
  
  end