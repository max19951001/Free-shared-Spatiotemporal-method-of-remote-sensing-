
; 说明：分块，istrum最小二乘拟合,残差分配根据阈值，大于阈值的直接认为是土地类型发生变化的，否则，按地物类型未发生变化处理

function cost_fun, g
  common v_pub1,ind_v
  common v_pub2,dep_v
  l=total((ind_v ## g-dep_v)^2)
  return, l
end

;打开数据
pro getdata,imgdata = imgdata,ns = ns,nl = nl,nb = nb,data_type = data_type,$
  filename = filename,map_info = map_info, fid = fid
  filter = ['all file;*.*']
  envi_open_file,filename,r_fid = fid
  envi_file_query,fid,ns = ns,nl = nl,nb = nb,data_type = data_type
  map_info = envi_get_map_info(fid=fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case data_type of
    1:imgdata = bytarr(ns,nl,nb)    ;  byte  byte
    2:imgdata = intarr(ns,nl,nb)    ;  int  integer
    3:imgdata = lonarr(ns,nl,nb)    ;  long  longword integer
    4:imgdata = fltarr(ns,nl,nb)    ;  float  floating point
    5:imgdata = dblarr(ns,nl,nb)    ;  double  double-precision floating
    6:imgdata = complexarr(ns,nl,nb); complex, single-precision, floating-point
    9:imgdata = dcomplexarr(ns,nl,nb);complex, double-precision, floating-point
    12:imgdata = uintarr(ns,nl,nb)   ; unsigned integer vector or array
    13:imgdata = ulonarr(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:imgdata = lon64arr(ns,nl,nb)   ;a 64-bit integer vector or array
    15:imgdata = ulon64arr(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  endcase
  for i = 0,nb-1 do begin
    dt = envi_get_data(fid = fid,dims = dims,pos=i)
    imgdata[*,*,i] = dt[*,*]
  endfor
end
;-----------------------------------------------------------------------------------------------------------------------------------
;主过程
pro ifsdafbymax2
  e=envi()
  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init
  t0=systime(1)

  ;输入参数
  ;-------------------------------------------------------------------------------------------------------------------------------------------
  ;输入数据为反射率值0~1之间
  ;fbase_fp="‪G:\gwydir\landsat20040705.dat"
  fbase_fp1=dialog_pickfile(title='打开基准landsat数据')                    ;基准landsat
  ;cbase_fp="‪G:\gwydir\modis20040705.dat"
  cbase_fp1=dialog_pickfile(title='打开基准mdois数据')                      ;基准modis
  ;cpre_fp="G:\gwydir\modis20040806.dat"
  cpre_fp1=dialog_pickfile(title='打开预测的modis数据')                     ;预测modis
  ;abundance=dialog_pickfile(title='打开丰度图')
  em=dialog_pickfile(title='打开csv文件')
  ;em="‪‪‪E:\OneDrive\IFSDAFBYMAX\svm4.csv"
  scale_factor=10                                                    ;modis和landsat像素比例
  w=15                                                               ;设置窗口的半径为25，则窗口的大小为25*2+1
  coarse_pixel=3                                                      ;粗像素的半径
  dn_min=0.0                                                           ;最小值
  dn_max=1.0                                                           ;最大值
  ;num_similar_pixel=100                                                ;设置相似像素
  ;num_pure=ceil((((coarse_pixel*2)+1)*((coarse_pixel*2)+1))/3)-3                                                           ;每一类最相似的像素
  block_size=50                                                     ;设置分块大小，也就是有多少个粗分辨率像素
  temp_file='e:\temp'                                                  ;临时路径
  ;-------------------------------------------------------------------------------------------------------------------------------------------


  routine_dir=file_dirname(routine_filepath('ifsdafbymax2'))+'\'  ;获取程序目录
  print,routine_dir
  ;1 打开端元光谱文件
  ;*****************************************************************************

  em_spec=read_csv(em,header=em_name)
  em_samples=n_elements(em_name) ;返回列数，也就是端元数
  em_lines  =n_elements(em_spec.(0));em_lines等于波段数
  temp=fltarr(em_samples,em_lines)
  for i=0,em_samples-1 do temp[i,*]=float(em_spec.(i))
  em_spec=temporary(temp);em_spec存储各个波段光谱值

  ;fcls
  cd,routine_dir
  fbaseabd_fp=file_dirname(cpre_fp1)+'\'+file_basename(fbase_fp1)+'_abundance.tif'
  print,fbaseabd_fp
  cmdstr='abundancecaculatemodule.exe '+fbase_fp1+' '+em+' '+fbaseabd_fp
  spawn,cmdstr,/hide
  envi_open_file,fbaseabd_fp,r_fid=fabd_fid
  if fabd_fid eq -1 then begin
    envi_batch_exit
    print,'光谱解混失败'
    return
  endif


   ;打开混合数据

  fbase_fp=dialog_pickfile(title='打开基准landsat混合数据')                    ;基准landsat
  ;cbase_fp="‪G:\gwydir\modis20040705.dat"
  cbase_fp=dialog_pickfile(title='打开基准mdois混合数据')                      ;基准modis
  ;cpre_fp1="G:\gwydir\modis20040806.dat"
  cpre_fp=dialog_pickfile(title='打开预测的modis混合数据')                     ;预测modis
  ;abundance=dialog_pickfile(title='打开丰度图')

  ;打开基准日期的landsat

  envi_open_file,fbase_fp,r_fid=fid0
  envi_file_query,fid0,ns=ns,nl=nl,nb=nb,dims=dims
  map_info = envi_get_map_info(fid = fid0)
  patch_long=block_size*scale_factor
  orig_ns=ns
  orig_nl=nl
  n_ns=ceil(float(orig_ns)/patch_long)
  n_nl=ceil(float(orig_nl)/patch_long) ;获取分块数



  ind_patch1=intarr(4,n_ns*n_nl)
  ind_patch=intarr(4,n_ns*n_nl)
  location=intarr(4,n_ns*n_nl)
  ;分块
  for i_ns=0,n_ns-1,1 do begin
    for i_nl=0,n_nl-1,1 do begin
      ind_patch1[0,n_ns*i_nl+i_ns]=i_ns*patch_long
      ind_patch[0,n_ns*i_nl+i_ns]=max([0,ind_patch1[0,n_ns*i_nl+i_ns]-scale_factor])
      location[0,n_ns*i_nl+i_ns]=ind_patch1[0,n_ns*i_nl+i_ns]-ind_patch[0,n_ns*i_nl+i_ns]

      ind_patch1[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
      ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,ind_patch1[1,n_ns*i_nl+i_ns]+scale_factor])
      location[1,n_ns*i_nl+i_ns]=ind_patch1[1,n_ns*i_nl+i_ns]-ind_patch1[0,n_ns*i_nl+i_ns]+location[0,n_ns*i_nl+i_ns]

      ind_patch1[2,n_ns*i_nl+i_ns]=i_nl*patch_long
      ind_patch[2,n_ns*i_nl+i_ns]=max([0,ind_patch1[2,n_ns*i_nl+i_ns]-scale_factor])
      location[2,n_ns*i_nl+i_ns]=ind_patch1[2,n_ns*i_nl+i_ns]-ind_patch[2,n_ns*i_nl+i_ns]

      ind_patch1[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
      ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,ind_patch1[3,n_ns*i_nl+i_ns]+scale_factor])
      location[3,n_ns*i_nl+i_ns]=ind_patch1[3,n_ns*i_nl+i_ns]-ind_patch1[2,n_ns*i_nl+i_ns]+location[2,n_ns*i_nl+i_ns]
    endfor
  endfor

  tempoutname=temp_file+'\temp_f1'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
    dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
    envi_doit, 'resize_doit', fid=fid0, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
    envi_file_mng, id=r_fid1, /remove
  endfor

  ;对基准的modis操作
  envi_open_file,cbase_fp,r_fid=fid
  tempoutname=temp_file+'\temp_c1'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
    dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
    envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
    envi_file_mng, id=r_fid1, /remove
  endfor
  envi_file_mng, id=fid, /remove

  ;对预测的modis操作
  envi_open_file,cpre_fp,r_fid=fid
  tempoutname=temp_file+'\temp_c0'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
    dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
    envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
    envi_file_mng, id=r_fid1, /remove
  endfor
  envi_file_mng, id=fid, /remove



  ;打开输入的landsat的丰度
  envi_open_file,fbaseabd_fp,r_fid=fid
  envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims
  tempoutname=temp_file+'\class'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
    dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
    envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
    envi_file_mng, id=r_fid1, /remove
  endfor

  envi_file_mng, id=fid, /remove
  envi_file_mng, id=fid0, /remove



  ;------------------------------------------------------------------
  ;处理每一块
  ;-------------------------------------------------------------------
  t0=systime(1)                  ;the initial time of program running

  print,'there are total',n_ns*n_nl,' blocks'

  for isub=0,n_ns*n_nl-1,1 do begin


    ;打开每一块文件

    filename = temp_file+'\temp_f1'
    getdata,imgdata=fine1,ns = ns,nl = nl,nb = nb,data_type = data_type,filename = filename+strtrim(isub+1,1),fid = fid11
    fine1=float(fine1)
    ;  print,ns
    ;  print,nl

    filename = temp_file+'\temp_c1'
    getdata,imgdata=coarse1,ns = ns,nl = nl,nb = nb,data_type = data_type,filename = filename+strtrim(isub+1,1),fid = fid12
    coarse1=float(coarse1)

    filename = temp_file+'\temp_c0'
    getdata,imgdata=coarse2,ns = ns,nl = nl,nb = nb,data_type = data_type,filename = filename+strtrim(isub+1,1),fid = fid13
    coarse2=float(coarse2)

    filename = temp_file+'\class'
    getdata,imgdata=l1_class0,ns = abd_ns,nl = abd_nl,nb=abd_nb,filename = filename+strtrim(isub+1,1),fid = fid14
    ;  print,abd_ns
    ;  print,abd_nl
    ;  print,abd_nb
    ;help,'l1_class0',l1_class0

    ;对输入的丰度图操作
    em_samples=abd_nb ;端元数
    ns_c=floor(ns/scale_factor)
    nl_c=floor(nl/scale_factor) ;获取每一块中modis像素的行列号
    ;  print,ns_c
    ;  print,nl_c
    ;进行丰度聚集，获取modis数据中的丰度图
    cabd_img=fltarr(ns_c,nl_c,em_samples)   ;cabd_img存放 modis丰度值
    for ci=0,nl_c-1 do begin
      for cj=0,ns_c-1 do begin
        ;逐MODIS像素遍历
        win_data=l1_class0[(cj*scale_factor):((cj+1)*scale_factor-1),(ci*scale_factor):((ci+1)*scale_factor-1),*] ;win_data存储的一个MODIS像素内所有Landsat像素的丰度值
        ;      print,win_data
        cabd=fltarr(em_samples);cabd存储的临时粗像素的丰度，也就是从landsat聚集到modis中
        for ib=0,em_samples-1 do begin
          cabd[ib]=mean(win_data[*,*,ib])
        endfor
        ;      index=where(cabd le 0.05,count);将小于0.05丰度的变为0
        ;      if count gt 0 then begin
        ;        cabd[index]=0
        ;      endif
        if total(cabd) lt 0.99 then begin ;将丰度小于1的残差赋值给端元丰度最大的那个
          redi=1.0-total(cabd);残差值
          sortCabd=sort(cabd)
          maxAbundance=cabd[sortCabd[abd_nb-1]]
          maxAbundance+=redi
          cabd[sortCabd[abd_nb-1]]=maxAbundance
        endif
        cabd_img[cj,ci,*]=cabd
      endfor
    endfor
    print,'丰度聚集结束'
    ;help,cabd_img

    ;(1)校正输入landsat过小的值
    for ib=0,nb-1, 1 do begin
      sortindex = sort(fine1[*,*,ib])
      sortindices = (findgen(float(ns)*nl+1))/(float(ns)*nl)
      percentiles=[0.0001, 0.9999]
      dataindices = value_locate(sortindices, percentiles)
      data_1_4= (fine1[*,*,ib])[sortindex[dataindices]]
      ;校正过于小的值
      ind_small=where(fine1[*,*,ib] le data_1_4[0] or fine1[*,*,ib] lt dn_min)
      temp=fine1[*,*,ib]
      temp[ind_small]=min((fine1[*,*,ib])[where(fine1[*,*,ib] gt data_1_4[0] and fine1[*,*,ib] ge dn_min)])
      fine1[*,*,ib]=temp
      ;校正过于大的值
      ind_large=where(fine1[*,*,ib] ge data_1_4[1] or fine1[*,*,ib] gt dn_max)
      temp=fine1[*,*,ib]
      temp[ind_large]=max((fine1[*,*,ib])[where(fine1[*,*,ib] lt data_1_4[1] and fine1[*,*,ib] le dn_max)])
      fine1[*,*,ib]=temp
    endfor
    ;(2)判断base modis数据
    for ib=0,nb-1, 1 do begin
      ;获取大于0.0001小于0.9999的索引位置
      percentiles=[0.0001, 0.9999]
      sortindices = (findgen(float(ns)*nl+1))/(float(ns)*nl);创建一个百分比数组
      dataindices = value_locate(sortindices, percentiles);dataindices返回0.0001和0.9999值在sortindices的位置
      cb_sortindex = sort(coarse1[*,*,ib]);返回从小到大排序的索引
      cb_data= (coarse1[*,*,ib])[cb_sortindex[dataindices]];返回一个两个元素数组，存储允许的最小和最大值
      ;校正过于小的值
      ind_small=where(coarse1[*,*,ib] le cb_data[0] or coarse1[*,*,ib] lt dn_min)
      temp=coarse1[*,*,ib]
      temp[ind_small]=min((coarse1[*,*,ib])[where(coarse1[*,*,ib] gt cb_data[0] and coarse1[*,*,ib] ge dn_min)])
      coarse1[*,*,ib]=temp
      ;校正过于大的值
      ind_large=where(coarse1[*,*,ib] ge cb_data[1] or coarse1[*,*,ib] gt dn_max)
      temp=coarse1[*,*,ib]
      temp[ind_large]=max((coarse1[*,*,ib])[where(coarse1[*,*,ib] lt cb_data[1] and coarse1[*,*,ib] le dn_max)])
      coarse1[*,*,ib]=temp
    endfor
    ;(3)判断pre modis数据
    for ib=0,nb-1, 1 do begin
      ;获取大于0.0001小于0.9999的索引位置
      percentiles=[0.0001, 0.9999]
      sortindices = (findgen(float(ns)*nl+1))/(float(ns)*nl);创建一个百分比数组
      dataindices = value_locate(sortindices, percentiles);dataindices返回0.0001和0.9999值在sortindices的位置
      cp_sortindex = sort(coarse2[*,*,ib]);返回从小到大排序的索引
      cp_data= (coarse2[*,*,ib])[cp_sortindex[dataindices]];返回一个两个元素数组，存储允许的最小和最大值
      ;校正过于小的值
      ind_small=where(coarse2[*,*,ib] le cp_data[0] or coarse2[*,*,ib] lt dn_min)
      temp=coarse2[*,*,ib]
      temp[ind_small]=min((coarse2[*,*,ib])[where(coarse2[*,*,ib] gt cp_data[0] and coarse2[*,*,ib] ge dn_min)])
      coarse2[*,*,ib]=temp
      ;校正过于大的值
      ind_large=where(coarse2[*,*,ib] ge cp_data[1] or coarse2[*,*,ib] gt dn_max)
      temp=coarse2[*,*,ib]
      temp[ind_large]=max((coarse2[*,*,ib])[where(coarse2[*,*,ib] lt cp_data[1] and coarse2[*,*,ib] le dn_max)])
      coarse2[*,*,ib]=temp
    endfor

    ;初始化index_f和index_c索引数组,目的是后面进行选择
    ii=0
    index_f=intarr(ns,nl)
    index_c=intarr(ns_c,nl_c)
    for i=0, ns_c-1, 1 do begin
      for j=0,nl_c-1,1 do begin
        index_f[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]=ii
        index_c[i,j]=ii
        ii=ii+1.0
      endfor
    endfor


    ;landsat行列索引值
    row_ind=intarr(ns,nl)
    col_ind=intarr(ns,nl)
    for i=0,ns-1,1 do begin
      col_ind[i,*]=i
    endfor
    for i=0,nl-1,1 do begin
      row_ind[*,i]=i
    endfor

    ;采样modis像素到modis分辨率，此处是480米
    fine_c1=fltarr(ns_c,nl_c,nb)
    coarse_c1=fltarr(ns_c,nl_c,nb)
    coarse_c2=fltarr(ns_c,nl_c,nb)
    row_c=fltarr(ns_c,nl_c)
    col_c=fltarr(ns_c,nl_c)
    for ic=0,ns_c-1, 1 do begin
      for jc=0,nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic,jc])
        row_c[ic,jc]= mean(row_ind[ind_c])
        col_c[ic,jc]= mean(col_ind[ind_c])
        for ib=0,nb-1,1 do begin
          fine_c1[ic,jc,ib]=mean((fine1[*,*,ib])[ind_c])
          coarse_c1[ic,jc,ib]=mean((coarse1[*,*,ib])[ind_c])
          coarse_c2[ic,jc,ib]=mean((coarse2[*,*,ib])[ind_c])
        endfor
      endfor
    endfor


    ;最小二乘解混
    common v_pub1
    common v_pub2
    gbnd    =[0,100]
    nobj    = 0
    lcomp   = 'cost_fun'

    ;存储的是每一块内每一波段端元的变化值
    c_rate=fltarr(em_samples,nb)
    change_f=make_array(ns,nl,nb,type=data_type);每一landsat像素端元的变化值
    change_c=coarse_c2-coarse_c1

    for ci=0,nl_c-1 do begin
      for cj=0,ns_c-1 do begin
        ai=max([0,ci-coarse_pixel])       ; 移动窗口直径
        bi=min([nl_c-1,ci+coarse_pixel])
        aj=max([0,cj-coarse_pixel])
        bj=min([ns_c-1,cj+coarse_pixel])
        ;        print,"ai",ai
        ;        print,"bi",bi
        ;        print,"aj",aj
        ;        print,"bj",bj
        ;搜索窗口大小
        search_win=(bi-ai+1)*(bj-aj+1)
        ;landsat窗口位置
        fai=ci*scale_factor
        fbi=(ci+1)*scale_factor-1
        faj=cj*scale_factor
        fbj=(cj+1)*scale_factor-1

        min_allow=fltarr(em_samples,nb);每一波段，端元的变化值
        max_allow=fltarr(em_samples,nb)
        ;空间解混开始
        if (search_win gt em_samples) then begin
          fabd_temp=l1_class0[faj:fbj,fai:fbi,*]        ;一个modis像素内landsat像素丰度存在fabd_temp
          ind_v=transpose(reform(cabd_img[aj:bj,ai:bi,*],search_win,em_samples));modis像素每一类的丰度值
          for ib=0,nb-1 do begin
            min_allow[*,ib]=min(coarse_c2[aj:bj,ai:bi,ib]-coarse_c1[aj:bj,ai:bi,ib])-stddev(coarse_c2[aj:bj,ai:bi,ib]-coarse_c1[aj:bj,ai:bi,ib]);对于每一波段，允许变化的最小值
            max_allow[*,ib]=max(coarse_c2[aj:bj,ai:bi,ib]-coarse_c1[aj:bj,ai:bi,ib])+stddev(coarse_c2[aj:bj,ai:bi,ib]-coarse_c1[aj:bj,ai:bi,ib]);对于每一波段，允许变化的最大值
            dep_v = double(reform(change_c[aj:bj,ai:bi,ib],search_win));为移动窗口内所有的modis像素变化值
            x     = fltarr(1,em_samples);最终的端元变化结果存在这个里面
            xbnd  = [[min_allow[*,ib]],[max_allow[*,ib]]];modis像素允许的变化范围
            constrained_min, x, xbnd, gbnd, nobj, lcomp,inform, nstop = 5 ;此处x返回的是每个modis像素内的端元变化范围
            ds_change=fabd_temp*rebin(reform(x,1,1,em_samples),scale_factor,scale_factor,em_samples,/sample);最邻近采样为了保证值不变的情况下进行矩阵相乘
            change_f[faj:fbj,fai:fbi,ib]=total(ds_change,3);change_f为每一landsat的变化值
            ;            x =reform(fine_c1[aj:bj,ai:bi,ib],search_win)
            ;            y =reform(coarse_c1[aj:bj,ai:bi,ib],search_win)
            ;            coefs=linfit(x,y)
            ;            fp_img_temp[faj:fbj,fai:fbi,ib]=fb_img[faj:fbj,fai:fbi,ib]+coefs[1]*change_f[faj:fbj,fai:fbi,ib]
          endfor
        endif
      endfor
    endfor
    ;print,'c_rate',c_rate ;求取的是每一块modis像素的变化值
    print,"光谱解混结束"




    ;调整传感器差异
    fp_img_temp=make_array(ns,nl,nb,type=data_type);时间预测结果
    for ib=0,nb-1  do begin
      x =reform(fine_c1[*,*,ib],ns_c*nl_c)
      y =reform(coarse_c1[*,*,ib],ns_c*nl_c);对窗口内modis分辨率的landsat数据和modis数据进行求解相关系数
      ; print,"x",x
      ;  print,"y",y
      coefs=linfit(x,y)
      fp_img_temp[*,*,ib]=fine1[*,*,ib]+(coefs[1]*change_f[*,*,ib])
      ;print,fp_img_temp[*,*,ib]
    endfor

    ; tempoutname11=temp_file+'\l2_1'
    ;  envi_write_envi_file,fp_img_temp,out_name = tempoutname11

    ;采样时间预测数据到modis分辨率
    coarse_c2_p=fltarr(ns_c,nl_c,nb)
    for ic=0,ns_c-1, 1 do begin
      for jc=0,nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic,jc])
        for ib=0,nb-1,1 do begin
          coarse_c2_p[ic,jc,ib]=mean((fp_img_temp[*,*,ib])[ind_c])
        endfor
      endfor
    endfor

    ;allowed minmum value for each band at tp
    min_allow=fltarr(nb)
    max_allow=fltarr(nb)
    for ib=0,nb-1,1 do begin
      min_allow0=min([min(coarse2[*,*,ib]),min(fp_img_temp[*,*,ib])])
      min_allow[ib]=max([min_allow0,dn_min])
      max_allow0=max([max(coarse2[*,*,ib]),max(fp_img_temp[*,*,ib])])
      max_allow[ib]=min([max_allow0,dn_max])
    endfor



    ;tps插值

    tps1=fltarr(ns,nl,nb)
    for ib=0,nb-1,1 do begin
      tps1[*,*,ib] = min_curve_surf(coarse_c1[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
    endfor

    ;预测日期modis插值
    tps2=fltarr(ns,nl,nb)
    for ib=0,nb-1,1 do begin
      tps2[*,*,ib] = min_curve_surf(coarse_c2[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
    endfor

    similar_tps=fltarr(nb)
    ;进行筛选地物类型发生变化的信息，假设地物类型变化发生在一个较小的区域。那么对于一个移动窗口内，若每一像素的tps差值大于similar_tps,则我们认为是地物类型发生变化的区域
    for ib=0,nb-1 do begin
      similar_tps[ib]=mean(abs(tps2[*,*,ib]-tps1[*,*,ib]))
      print,'simiar_tps====',similar_tps[ib]
    endfor
    ;tps插值
    ;    l2_tps=fltarr(ns,nl,nb)
    ;    for ib=0,nb-1,1 do begin
    ;      l2_tps[*,*,ib] = min_curve_surf(coarse_c2[*,*,ib], col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
    ;    endfor
    ;
    ;    print,'finish tps prediction'

    ;计算fb_img中的相似像素的阈值,后面进行计算同质系数以及加权求和要用
    similar_th=fltarr(nb)
    for iband=0,nb-1,1 do begin
      similar_th[iband]=stddev(fine1[*,*,iband])*2.5/em_samples
    endfor
    ;同质系数计算
    het_index=fltarr(ns,nl,nb)
    for fi=0,nl-1 do begin
      for fj=0,ns-1 do begin
        ai=max([0,fi-w])
        bi=min([nl-1,fi+w]) ;搜索窗口，[fj,fi]为中心像素
        aj=max([0,fj-w])
        bj=min([ns-1,fj+w])
        ;    print,"lajissssssssssssssssssssssssssss",ai,bi,aj,bj
        temp=fltarr(bj-aj+1,bi-ai+1,nb)
        for ib=0,nb-1 do begin
          temp[*,*,ib]=abs(tps2[aj:bj,ai:bi,ib]-tps1[aj:bj,ai:bi,ib])
          index=where(temp[*,*,ib] gt (similar_tps[ib]),nums)
          ;print,'nums======',nums
          if nums gt 0 then begin
            het_index[fj,fi,ib]=nums/((bj-aj+1)*(bi-ai+1)) ;此处he_index表示的是异质性区域的系数，越大的值表示赋值更多给异质性区域
          endif else begin
            het_index[fj,fi,ib]=0.0
          endelse
        endfor
      endfor
    endfor
    ;    het_index=fltarr(ns,nl,nb)
    ;    for fi=0,nl-1 do begin
    ;      for fj=0,ns-1 do begin
    ;        ai=max([0,fi-w])
    ;        bi=min([nl-1,fi+w]) ;搜索窗口，[fj,fi]为中心像素
    ;        aj=max([0,fj-w])
    ;        bj=min([ns-1,fj+w])
    ;        ;    print,"lajissssssssssssssssssssssssssss",ai,bi,aj,bj
    ;        temp=fltarr(bj-aj+1,bi-ai+1,nb)
    ;        for ib=0,nb-1 do begin
    ;          temp[*,*,ib]=fine1[aj:bj,ai:bi,ib]-fine1[fj,fi,ib]
    ;          index=where(temp[*,*,ib] lt similar_th[ib],nums)
    ;          if nums gt 0 then begin
    ;            het_index[fj,fi,ib]=nums/((bj-aj+1)*(bi-ai+1))
    ;          endif else begin
    ;            het_index[fj,fi,ib]=0.0
    ;          endelse
    ;        endfor
    ;      endfor
    ;    endfor

    ;残差分配
    predict_change_c=coarse_c2_p-fine_c1     ;predict change
    real_change_c=coarse_c2-coarse_c1        ;real change
    ;change_r=coarse_c2-coarse_c2_p
    change_r=real_change_c-predict_change_c

    ;redistribute residual
    change_21_c=fltarr(ns_c,nl_c,nb)
    change_21=fltarr(ns,nl,nb)

    for ic=0,ns_c-1, 1 do begin
      ; print,'finish coarse pixel:',ic+1,'column'
      for jc=0,nl_c-1, 1 do begin
        ind_c=where(index_f eq index_c[ic,jc],num_ii)
        for ib=0,nb-1 do begin

          diff_change=change_r[ic,jc,ib];为地物类型发生变化的残差
          w_change_tps=(tps2[*,*,ib])[ind_c]-(fp_img_temp[*,*,ib])[ind_c];为地位类型未发生变化的残差
          ;             print,'diff_change',diff_change
          if (diff_change le 0) then begin
            ind_noc=where(w_change_tps gt 0, num_noc)
            if (num_noc gt 0) then begin
              w_change_tps[ind_noc]=0
            endif
          endif else begin
            ind_noc=where(w_change_tps lt 0, num_noc)
            if (num_noc gt 0) then begin
              w_change_tps[ind_noc]=0
            endif
          endelse
          w_change_tps=abs(w_change_tps)

          w_unform=fltarr(num_ii)     ;evenly distributing residuals to sub-pixels；num_ii为landat像素的个数
          w_unform[*]=abs(diff_change)

;          ;w_change=w_change_tps*het_index[ind_c]+w_unform*(1.0-het_index[ind_c])+0.000001  ;combine these two weights
          w_change=w_unform*het_index[ind_c]+w_change_tps*(1.0-het_index[ind_c])+0.000001
          w_change=w_change/(mean(w_change)) ;nomalize weight

          ;avoid extreme weights
          ind_extrem=where(w_change gt 1.0, num_extrem)
          if (num_extrem gt 0) then begin
            w_change[ind_extrem]=mean(w_change)
          endif
          w_change=w_change/(mean(w_change))
         
;        
          ;distribute residuals according to weight
          temp=change_21[*,*,ib]
          temp[ind_c]=w_change*diff_change
          change_21[*,*,ib]=temp
        endfor
      endfor
    endfor

    ;second prediction: l1+change
    fine2_2=fp_img_temp+change_21
    ;correct abnormal detected change
    for ib=0,nb-1, 1 do begin
      temp=fine2_2[*,*,ib]
      ind_min=where(temp lt min_allow[ib], num_min)
      if (num_min gt 0) then begin
        temp[ind_min]=min_allow[ib]
      endif
      ind_max=where(temp gt max_allow[ib], num_max)
      if (num_max gt 0) then begin
        temp[ind_max]=max_allow[ib]
      endif
      fine2_2[*,*,ib]=temp
    endfor
    change_21=fine2_2-fine1

    print,'finish change prediction step ',isub+1,' block'
    tempoutname1=temp_file+'\temp_change'
    envi_write_envi_file,change_21,out_name = tempoutname1+strtrim(isub+1,1)
    envi_file_mng, id=fid11, /remove
    envi_file_mng, id=fid12, /remove
    envi_file_mng, id=fid13, /remove
    envi_file_mng, id=fid14, /remove

  endfor

  fids = envi_get_file_ids()
  help,fids
  length=n_elements(fids)
  print,length
  for i = 0, length-1 do begin
    envi_file_mng,id = fids[i],/remove
  endfor

  mfid=intarr(n_ns*n_nl)
  mdims=intarr(5,n_ns*n_nl)
  mpos=intarr(nb,n_ns*n_nl)
  pos=indgen(nb)
  x0=intarr(n_ns*n_nl)
  y0=intarr(n_ns*n_nl)

  for isub=0,n_ns*n_nl-1,1 do begin
    envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
    if (sub_fid eq -1) then begin
      envi_batch_exit
      return
    endif
    envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
    mfid[isub] = sub_fid
    mpos[*,isub] = indgen(nb)
    mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
    x0[isub]=ind_patch1[0,isub]
    y0[isub]=ind_patch1[2,isub]
  endfor

  xsize = orig_ns
  ysize = orig_nl
  pixel_size = [1.,1.]

  use_see_through = replicate(1l,n_ns*n_nl)
  see_through_val = replicate(0l,n_ns*n_nl)

  out_name=temp_file+'_change'
  envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,map_info=map_info, $
    out_dt=4, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

  for i=0,n_ns*n_nl-1,1 do begin
    envi_file_mng, id=mfid[i], /remove, /delete
  endfor

  ;权重预测最终结果
  ;;##############step 5: final prediction

  filename6 = out_name
  envi_open_file,filename6,r_fid=fid
  tempoutname=temp_file+'\temp_change'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
    dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
    envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
    envi_file_mng, id=r_fid1, /remove
  endfor
  envi_file_mng, id=fid, /remove,/delete


  for isub=0,n_ns*n_nl-1,1 do begin

    ;open each block image

    filename = temp_file+'\temp_f1'
    getdata,imgdata=fine1,ns = ns,nl = nl,nb = nb,data_type = data_type,filename = filename+strtrim(isub+1,1),fid = fid11
    fine1=float(fine1)

    filename = temp_file+'\temp_c1'
    getdata,imgdata=coarse1,filename = filename+strtrim(isub+1,1),fid = fid12
    coarse1=float(coarse1)

    filename = temp_file+'\temp_c0'
    getdata,imgdata=coarse2,filename = filename+strtrim(isub+1,1),fid = fid13
    coarse2=float(coarse2)

    filename = temp_file+'\class'
    getdata,imgdata=l1_class,filename = filename+strtrim(isub+1,1),fid = fid14

    filename = temp_file+'\temp_change'
    getdata,imgdata=change_21,filename = filename+strtrim(isub+1,1),fid = fid15
    change_21=float(change_21)
    ;place the blended result
    fine2=fine1[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]



    ;compute the distance of each pixel in the window with the target pixel (integrate window)
    d_d_all=((w-indgen(w*2+1)#(intarr(1,w*2+1)+1))^2+(w-(intarr(w*2+1)+1)#indgen(1,w*2+1))^2)^0.5
    d_d_all=reform(d_d_all,(w*2+1)*(w*2+1))

    similar_th=fltarr(nb)
    for iband=0,nb-1,1 do begin
      similar_th[iband,0]=stddev(fine1[*,*,iband])*2.5/em_samples
    endfor

    for i=location[0,isub],location[1,isub],1 do begin           ;retieve each target pixel
      for j=location[2,isub],location[3,isub],1 do begin
        ; if (fine1[i,j,background_band-1] ne background  ) then begin    ;do not process the background

        ai=max([0,i-w])       ; the window location
        bi=min([ns-1,i+w])
        aj=max([0,j-w])
        bj=min([nl-1,j+w])

        ci=i-ai      ;location of target pixel
        cj=j-aj
        col_wind=indgen(bi-ai+1)#(intarr(1,bj-aj+1)+1)*1.0
        row_wind=(intarr(bi-ai+1)+1)#indgen(1,bj-aj+1)*1.0

        ;search similar pixels within window
        similar_cand=fltarr((bi-ai+1)*(bj-aj+1)) ;place the similarity measure between each pixel and the target pixel
        position_cand=intarr((bi-ai+1)*(bj-aj+1))+1  ;place the location of each similar pixel
        for iband=0,nb-1,1 do begin
          cand_band=intarr((bi-ai+1)*(bj-aj+1))
          wind_fine=fine1[ai:bi,aj:bj,iband]
          s_s=abs(wind_fine-wind_fine[ci,cj])
          similar_cand=similar_cand+s_s/(wind_fine[ci,cj]+0.00000001)
          ind_cand=where(s_s lt similar_th[iband])
          cand_band[ind_cand]=1
          position_cand=position_cand*cand_band
        endfor

        indcand=where(position_cand ne 0,number_cand0)  ;select similar pixel initially
        order_dis=sort(similar_cand[indcand])
        ;print,'number-cand0',number_cand0
        ;number_cand=min([number_cand0,num_similar_pixel])
        number_cand=number_cand0
        ind_same_class=indcand[order_dis[0:number_cand-1]]           ; select the n most similar samples

        ;compute weight for each simialr pixel
        ;spatial distance
        if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not an integrate window
          d_d_cand=((ci-col_wind[ind_same_class])^2+(cj-row_wind[ind_same_class])^2)^0.5+0.00001
        endif else begin
          d_d_cand=d_d_all[ind_same_class]      ;integrate window
        endelse

        ;normalize these distances
        d_d_cand=(1.0+d_d_cand/w)*(similar_cand[ind_same_class]+1.0)
        c_d=1.0/d_d_cand
        weight=c_d/total(c_d)

        for iband=0,nb-1,1 do begin
          ;predict the value
          change_cand=(change_21[ai:bi,aj:bj,iband])[ind_same_class]
          fine2[i-location[0,isub],j-location[2,isub],iband]=fine1[i,j,iband]+total(weight*change_cand)

          ;revise the abnormal prediction
          if (fine2[i-location[0,isub],j-location[2,isub],iband] lt dn_min ) then begin ;correct abnomal prediction
            another_predict=max([dn_min, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
            fine2[i-location[0,isub],j-location[2,isub],iband]=min([dn_max,another_predict])
          endif
          if (fine2[i-location[0,isub],j-location[2,isub],iband] gt dn_max ) then begin ;correct abnomal prediction
            another_predict=min([dn_max, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
            fine2[i-location[0,isub],j-location[2,isub],iband]=max([dn_min,another_predict])
          endif
        endfor


        ;endif
      endfor

    endfor

    ; change the type of prediction into the type same as the input image
    case data_type of
      1:fine2 = byte(fine2)    ;  byte  byte
      2:fine2 = fix(fine2)     ;  int  integer
      3:fine2 = long(fine2)    ;  long  longword integer
      4:fine2 = float(fine2)   ;  float  floating point
      5:fine2 = double(fine2)  ;  double  double-precision floating
      6:fine2 = complex(fine2); complex, single-precision, floating-point
      9:fine2 = dcomplex(fine2);complex, double-precision, floating-point
      12:fine2 = uint(fine2)   ; unsigned integer vector or array
      13:fine2 = ulong(fine2)   ;  unsigned longword integer vector or array
      14:fine2 = long64(fine2)   ;a 64-bit integer vector or array
      15:fine2 = ulong64(fine2)   ;an unsigned 64-bit integer vector or array
    endcase

    print,'finish final prediction ',isub+1,' block'
    tempoutname1=temp_file+'\temp_blended'
    envi_write_envi_file,fine2,out_name = tempoutname1+strtrim(isub+1,1)
    envi_file_mng, id=fid11, /remove, /delete
    envi_file_mng, id=fid12, /remove, /delete
    envi_file_mng, id=fid13, /remove, /delete
    envi_file_mng, id=fid14, /remove, /delete
    envi_file_mng, id=fid15, /remove, /delete

  endfor

  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
  ;#########################################
  ;mosiac all the blended patch

  mfid=intarr(n_ns*n_nl)
  mdims=intarr(5,n_ns*n_nl)
  mpos=intarr(nb,n_ns*n_nl)
  pos=indgen(nb)
  x0=intarr(n_ns*n_nl)
  y0=intarr(n_ns*n_nl)

  for isub=0,n_ns*n_nl-1,1 do begin
    envi_open_file, tempoutname1+strtrim(isub+1,1), r_fid= sub_fid
    if (sub_fid eq -1) then begin
      envi_batch_exit
      return
    endif
    envi_file_query,  sub_fid, ns=sub_ns, nl=sub_nl
    mfid[isub] = sub_fid
    mpos[*,isub] = indgen(nb)
    mdims[*,isub] = [-1,0, sub_ns-1,0, sub_nl-1]
    x0[isub]=ind_patch1[0,isub]
    y0[isub]=ind_patch1[2,isub]
  endfor

  xsize = orig_ns
  ysize = orig_nl
  pixel_size = [1.,1.]

  use_see_through = replicate(1l,n_ns*n_nl)
  see_through_val = replicate(0l,n_ns*n_nl)

  out_name=cpre_fp+'_fsdaf20_1'
  envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,map_info=map_info, $
    out_dt=data_type, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

  for i=0,n_ns*n_nl-1,1 do begin
    envi_file_mng, id=mfid[i], /remove, /delete
  endfor




end