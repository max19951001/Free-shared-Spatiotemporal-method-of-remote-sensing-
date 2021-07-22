
;---------------------------------------------------------------------------------
;                          A new spatiotemporal data fusion model        
;                          Using one pairs of fine and coarse images
;                          The program can be used for whole TM scene
;      Note: this version requires users to input pre-classified fine image at t1. 
;      This version is appropriate for areas with complex land cover types and training samples
;      so users can first classify the fine image by supervised classifers like SVM. 
;        
;                   Developed by Zhu Xiaolin,email: zhuxiaolin55@gmail.com
;                        
;             Update history
;             03/09/2016, debug for images with background pixels (e.g. pixels with 0 values)
;                         add parameters "background" and "background_band" to identify background pixels
;             07/01/2017, improve the similar pixel selection process to reduce effect of errors in classification
;             03/11/2018, improve efficiency
;             01/23/2019, debug for when no purest non-change coarse pixels are seleted for computing temporal change
;             02/15/2019, correct the problem of selecting similar pixels when fusing simulation image with many pixels of exactly same values
;                          
;                           Copyright belongs to Xiaolin Zhu
;                           
; please cite: Zhu, X., Helmer E., Liu, D., Chen, J., Gao, F., and Lefsky M. A flexible spatiotemporal
; method for fusing satellite images with different resolutions. Remote Sensing of Environment,doi:10.1016/j.rse.2015.11.016
;---------------------------------------------------------------------------------


;function for open the file

Pro GetData,ImgData = ImgData,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,$
    FileName = FileName,Map_info = map_Info, Fid = Fid
    Filter = ['all file;*.*']
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

;function of multiple linear regression without interception
Function P_OLS, ind_v,dep_v,k,low,high    ;format: ind_v:k*n,dep_v:1*n

  common V_PUB, matrix
  common V_PUB2, y

  nb=n_elements(dep_v)
  x=reform(dep_v)
  varend=fltarr(k)
  yend=fltarr(k,nb)
  for ii=0,k-1 do begin
    varend[ii]=sqrt(total((ind_v[ii,*]-mean(ind_v[ii,*]))^2))
    yend[ii,*]=ind_v[ii,*]
  endfor
  var=sqrt(total((dep_v-mean(dep_v))^2))

  y=dep_v
  matrix=yend
  y=transpose(double(y))

  glow=fltarr(k)+low
  ghigh=fltarr(k)+high
  gintial=fltarr(k)+1.0/float(k)
  gbnd   = [glow, ghigh]
  Lbnd =[0,100]
  nobj   = 1
  g      = gintial
  Lcomp = 'HMBL11'
  nobj=0
  CONSTRAINED_MIN, g, gbnd, Lbnd, nobj, Lcomp, inform
  L =total((matrix ## g- y)^2)
  return,g
END

FUNCTION HMBL11, g
  common V_PUB
  common V_PUB2
  L=total((matrix ## g-y)^2)
  RETURN, L
END



;-------------------------------------------------------------------
;                       main program
;-------------------------------------------------------------------

pro  FSDAF_preclassification

 
;please set the following parameters
;----------------------------------------------------------------------
 w=20                 ;set the half window size, if 25, the window size is 25*2+1=51
 num_similar_pixel=20         ;set number of similar pixels
 num_pure=100                 ;number of most purest coarse pixels in each class selected fro change calculation 
 DN_min=0.0                   ;set the range of DN value of the image,If byte, 0 and 255
 DN_max=10000.0
 scale_factor=16              ;set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
 block_size=30                ;set the size of block, e.g., 20 means 20*20 coarse pixels, if process whole ETM scene, set 30~50
 background=0                 ;set the value of background pixels. 0 means that pixels will be considered as background if one of its bands= 0
 background_band=3            ;which band with value = background indicating background pixels. Sometimes, background pixels have different values in different bands 
 temp_file='D:\temp'          ;set path of a folder to store temporary files
;------------------------------------------------------------------------


;open the fine image of the first pair
FileName1 = Dialog_PickFile(title = 'open the fine image of the first pair:')
envi_open_file,FileName1,r_fid=fid0
envi_file_query,fid0,ns=ns,nl=nl,nb=nb,dims=dims
map_info = envi_get_map_info(fid = fid0)
patch_long=block_size*scale_factor
orig_ns=ns
orig_nl=nl
n_ns=ceil(float(orig_ns)/patch_long)
n_nl=ceil(float(orig_nl)/patch_long)

  ind_patch1=intarr(4,n_ns*n_nl)           ;divide the whole scene into 1000*1000 block
  ind_patch=intarr(4,n_ns*n_nl)
  location=intarr(4,n_ns*n_nl)

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

tempoutname=temp_file+'\temp_F1'

pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid0, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor




;open the coarse image of the first pair
;-----------------------------------------------------------
FileName2 = Dialog_PickFile(title = 'open the coarse image of the first pair:')
envi_open_file,FileName2,r_fid=fid
tempoutname=temp_file+'\temp_C1'
pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor
envi_file_mng, id=fid, /remove


;open the coarse image of the prediction time
;-----------------------------------------------------------
FileName5 = Dialog_PickFile(title = 'open the coarse image of the prediction time:')
envi_open_file,FileName5,r_fid=fid
tempoutname=temp_file+'\temp_C0'
pos=indgen(nb)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor
envi_file_mng, id=fid, /remove


;open the class image
;-----------------------------------------------------------
out_name = Dialog_PickFile(title = 'open the class image of fine image in the 1st pair:')
envi_open_file,out_name,r_fid=fid
tempoutname=temp_file+'\class'
pos=indgen(1)
for isub=0,n_ns*n_nl-1,1 do begin
  dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
  envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
    out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
  envi_file_mng, id=r_fid1, /remove
endfor

envi_file_mng, id=fid, /remove
envi_file_mng, id=fid0, /remove

;

;------------------------------------------------------------------
;process  each block
;-------------------------------------------------------------------
t0=systime(1)                  ;the initial time of program running

print,'there are total',n_ns*n_nl,' blocks'

for isub=0,n_ns*n_nl-1,1 do begin

  ;open each block image

  FileName = temp_file+'\temp_F1'
  GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid11
  fine1=float(fine1)

  FileName = temp_file+'\temp_C1'
  GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid12
  coarse1=FLOAT(coarse1)
  
  FileName = temp_file+'\temp_C0'
  GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid13
  coarse2=FLOAT(coarse2)

  FileName = temp_file+'\class'
  GetData,ImgData=L1_class0,FileName = FileName+strtrim(isub+1,1),Fid = Fid14
  
  num_class=max(L1_class0)
  ;recode the classification map if the subset does not have all classes
  i_new_c=0
  L1_class=intarr(ns,nl)
  for iclass=0, num_class-1,1 do begin
    ind_ic=where(L1_class0 eq iclass+1 and fine1[*,*, background_band-1] ne background, num_ic)
    if (num_ic gt 0) then begin
      L1_class[ind_ic]=i_new_c+1
      i_new_c=i_new_c+1
    endif
  endfor

  num_class=max(L1_class)
  
  
if (num_class gt 0) then begin   ;do not process if the whole subset is background
 

;correct extreme noise in fine1 becase extreme values will affect the allowed data range
for ib=0,nb-1, 1 do begin
  sortIndex = Sort(fine1[*,*,ib])
  sortIndices = (Findgen(float(ns)*nl+1))/(float(ns)*nl)
  Percentiles=[0.0001, 0.9999]
  dataIndices = Value_Locate(sortIndices, Percentiles)
  data_1_4= (fine1[*,*,ib])[sortIndex[dataIndices]]
  ;correct too small values
  ind_small=where(fine1[*,*,ib] le data_1_4[0] or fine1[*,*,ib] lt DN_min)
  temp=fine1[*,*,ib]
  temp[ind_small]=min((fine1[*,*,ib])[where(fine1[*,*,ib] gt data_1_4[0] and fine1[*,*,ib] ge DN_min)])
  fine1[*,*,ib]=temp
  ;correct too large values
  ind_large=where(fine1[*,*,ib] ge data_1_4[1] or fine1[*,*,ib] gt DN_max)
  temp=fine1[*,*,ib]
  temp[ind_large]=max((fine1[*,*,ib])[where(fine1[*,*,ib] lt data_1_4[1] and fine1[*,*,ib] le DN_max)])
  fine1[*,*,ib]=temp
endfor


;get index image between coarse and fine resolutions
ii=0
ns_c=floor(ns/scale_factor)
nl_c=floor(nl/scale_factor)
index_f=intarr(ns,nl)
index_c=intarr(ns_c,nl_c)
for i=0, ns_c-1, 1 do begin
  for j=0,nl_c-1,1 do begin
    index_f[i*scale_factor:(i+1)*scale_factor-1, j*scale_factor:(j+1)*scale_factor-1]=ii
    index_c[i,j]=ii
    ii=ii+1.0
  endfor
endfor

        
  ;col and row index
  row_ind=intarr(ns,nl)
  col_ind=intarr(ns,nl)
  for i=0,ns-1,1 do begin
    col_ind[i,*]=i
  endfor
  for i=0,nl-1,1 do begin
    row_ind[*,i]=i
  endfor
   
  ;resample coarse image to coarse resolution
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
  

;print,'number of class:',num_class

;step 2: get fracture of each class within each coarse pixel at t1
   Fraction1=fltarr(ns_c,nl_c,num_class)
   for ic=0,ns_c-1, 1 do begin
     for jc=0,nl_c-1, 1 do begin
       ind_c=where(index_f eq index_c[ic,jc], num_c)
       L1_class_c=L1_class[ind_c]
       for iclass=0, num_class-1,1 do begin
         ind_ic=where(L1_class_c eq iclass+1, num_ic)
         Fraction1[ic,jc,iclass]=float(num_ic)/float(num_c)         
       endfor
       if (total(Fraction1[ic,jc,*]) le 0.999) then begin   ;avoild pixels have background fine pixels
           Fraction1[ic,jc,*]=0
       endif
     endfor
   endfor

;get the heterogenity of each fine pixel
het_index=fltarr(ns,nl)
scale_d=w

for i=0,ns-1, 1 do begin
  for j=0,nl-1, 1 do begin  
    ai=max([0,i-scale_d])       ; the window location
    bi=min([ns-1,i+scale_d])
    aj=max([0,j-scale_d])
    bj=min([nl-1,j+scale_d])
    class_t=L1_class[i,j]
    ;select same-class pixels
    ind_same_class=where(L1_class[ai:bi,aj:bj] eq class_t, num_sameclass)
    het_index[i,j]=float(num_sameclass)/((bi-ai+1.0)*(bj-aj+1.0))
  endfor
endfor

;tempoutname11=temp_file+'\het_index'
;Envi_Write_Envi_File,het_index,Out_Name = tempoutname11



;step 3: METHOD2:estimate average spectral change of each class using pixels without land cover change
c_rate=fltarr(num_class,nb)

;allowed change value for each band
min_allow=fltarr(nb)
max_allow=fltarr(nb)
for ib=0,nb-1,1 do begin
  min_allow[ib]=min(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])-stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
  max_allow[ib]=max(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])+stddev(coarse_c2[*,*,ib]-coarse_c1[*,*,ib])
endfor

for ib=0,nb-1, 1 do begin
  X_matrix=fltarr(num_class,num_pure*num_class)
  y_matrix=fltarr(1,num_pure*num_class)
  ii=0
  for ic=1, num_class, 1 do begin
    order=sort(Fraction1[*,*,ic-1])
    order=reverse(order)
    ind_f=where(Fraction1[*,*,ic-1] gt 0.01, num_f) ;make sure all selected modis pixel contain class i
    num_pure1=min([num_f,num_pure])
    change_c=(coarse_c2[*,*,ib])[order[0:(num_pure1-1)]]-(coarse_c1[*,*,ib])[order[0:(num_pure1-1)]]        
   
    ;only use 0.1-0.9 samples to exclude the land cover change pixels
    sortIndex = Sort(change_c)
    sortIndices = (Findgen(num_pure1+1))/float(num_pure1)
    ;Percentiles=[0.25, 0.75] ; use this one if large land cover change existed
    Percentiles=[0.1, 0.9]
    dataIndices = Value_Locate(sortIndices, Percentiles)
    data_1_4= change_c[sortIndex[dataIndices]]   
    ind_nonchange=where(change_c ge data_1_4[0] and change_c le data_1_4[1],num_nonc) 
    if (num_nonc gt 0) then begin          
      y_matrix[0,ii:ii+num_nonc-1]=change_c[ind_nonchange]  
      for icc=0, num_class-1,1 do begin          
          f_c=(Fraction1[*,*,icc])[order[0:(num_pure1-1)]]      
          X_matrix[icc,ii:ii+num_nonc-1]=f_c[ind_nonchange]
      endfor
      ii=ii+num_nonc
   endif
  endfor
  X_matrix=X_matrix[*,0:ii-1]
  y_matrix=y_matrix[0,0:ii-1]
  result=P_OLS(X_matrix,y_matrix,num_class,min_allow[ib],max_allow[ib])  ;format: ind_v:k*n,dep_v:1*n
  c_rate[*,ib]=result[*]
endfor


;;step4: predict L2 assuming no land cover change
L2_1=fine1
for ic=1, num_class, 1 do begin
  ind_L1_class=where(L1_class eq ic)
  for ib=0,nb-1, 1 do begin
    temp=L2_1[*,*,ib]
    temp[ind_L1_class]=(fine1[*,*,ib])[ind_L1_class]+c_rate[ic-1,ib]
    L2_1[*,*,ib]=temp
  endfor
endfor

;tempoutname11=temp_file+'\L2_1'
;Envi_Write_Envi_File,L2_1,Out_Name = tempoutname11


;resample L2_1 image to coarse resolution

coarse_c2_p=fltarr(ns_c,nl_c,nb)
for ic=0,ns_c-1, 1 do begin
  for jc=0,nl_c-1, 1 do begin
    ind_c=where(index_f eq index_c[ic,jc])
    for ib=0,nb-1,1 do begin
      coarse_c2_p[ic,jc,ib]=mean((L2_1[*,*,ib])[ind_c])      
    endfor
  endfor
endfor



;allowed minmum value for each band at t2
min_allow=fltarr(nb)
max_allow=fltarr(nb)
for ib=0,nb-1,1 do begin
  min_allow0=min([min(coarse2[*,*,ib]),min(L2_1[*,*,ib])])
  min_allow[ib]=max([min_allow0,DN_min])
  max_allow0=max([max(coarse2[*,*,ib]),max(L2_1[*,*,ib])])
  max_allow[ib]=min([max_allow0,DN_max])
endfor



;step5: predict L2 using TPS 
L2_tps=fltarr(ns,nl,nb)
for ib=0,nb-1,1 do begin
  L2_tps[*,*,ib] = MIN_CURVE_SURF(coarse_c2[*,*,ib], col_c, row_c,/TPS, XPOUT=col_ind,  YPOUT=row_ind)
endfor

print,'finish TPS prediction'

;tempoutname11=temp_file+'\L2_tps'
;Envi_Write_Envi_File,L2_tps,Out_Name = tempoutname11



;step 6: redistribute residual
;change residual
predict_change_c=coarse_c2_p-fine_c1     ;predict change
real_change_c=coarse_c2-coarse_c1        ;real change
;change_R=coarse_c2-coarse_c2_p
change_R=real_change_c-predict_change_c

     ;redistribute residual
     change_21_c=fltarr(ns_c,nl_c,nb)
     change_21=fltarr(ns,nl,nb)
  
    for ic=0,ns_c-1, 1 do begin
        ; print,'finish coarse pixel:',ic+1,'column'
       for jc=0,nl_c-1, 1 do begin
  
         ind_c=where(index_f eq index_c[ic,jc],num_ii)
        
           for ib=0,nb-1,1 do begin
             
              diff_change=change_R[ic,jc,ib]  
              w_change_tps=(L2_tps[*,*,ib])[ind_c]-(L2_1[*,*,ib])[ind_c]
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
  
             w_unform=fltarr(num_ii)     ;evenly distributing residuals to sub-pixels 
             w_unform[*]=abs(diff_change)
              
             w_change=w_change_tps*het_index[ind_c]+w_unform*(1.0-het_index[ind_c])+0.000001  ;combine these two weights
             
             w_change=w_change/(mean(w_change)) ;nomalize weight
             
             ;avoid extreme weights
             ind_extrem=where(w_change gt 10, num_extrem)
             if (num_extrem gt 0) then begin
               w_change[ind_extrem]=mean(w_change)
             endif
             w_change=w_change/(mean(w_change))
          
                      
             ;distribute residuals according to WEIGHT           
             temp=change_21[*,*,ib]
             temp[ind_c]=w_change*diff_change            
             change_21[*,*,ib]=temp

           endfor

       endfor
     endfor
     
;     tempoutname11=temp_file+'\change_21'
;     Envi_Write_Envi_File,change_21,Out_Name = tempoutname11
;        
     ;second prediction: L1+change
     fine2_2=L2_1+change_21
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

 endif else begin
    change_21=fine1-fine1
 endelse
;
 change_21=change_21[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]
 change_21=float(change_21)

print,'finish change prediction step ',isub+1,' block'
tempoutname1=temp_file+'\temp_change'
Envi_Write_Envi_File,change_21,Out_Name = tempoutname1+strtrim(isub+1,1)
envi_file_mng, id=Fid11, /remove
envi_file_mng, id=Fid12, /remove
envi_file_mng, id=Fid13, /remove
envi_file_mng, id=Fid14, /remove

endfor

;
;;--------------------------------------------------------------------------------------
;mosiac all the change patch

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

  use_see_through = replicate(1L,n_ns*n_nl)
  see_through_val = replicate(0L,n_ns*n_nl)

    out_name=temp_file+'_change'
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
    out_dt=4, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

    for i=0,n_ns*n_nl-1,1 do begin
      envi_file_mng, id=mfid[i], /remove, /delete
    endfor



;;##############step 5: final prediction

FileName6 = out_name
envi_open_file,FileName6,r_fid=fid
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

  FileName = temp_file+'\temp_F1'
  GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid11
  fine1=float(fine1)

  FileName = temp_file+'\temp_C1'
  GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid12
  coarse1=FLOAT(coarse1)

  FileName = temp_file+'\temp_C0'
  GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid13
  coarse2=FLOAT(coarse2)

  FileName = temp_file+'\class'
  GetData,ImgData=L1_class,FileName = FileName+strtrim(isub+1,1),Fid = Fid14

  FileName = temp_file+'\temp_change'
  GetData,ImgData=change_21,FileName = FileName+strtrim(isub+1,1),Fid = Fid15
  change_21=FLOAT(change_21)
  ;place the blended result
  fine2=fine1[location[0,isub]:location[1,isub],location[2,isub]:location[3,isub],*]



;compute the distance of each pixel in the window with the target pixel (integrate window)
D_D_all=((w-indgen(w*2+1)#(intarr(1,w*2+1)+1))^2+(w-(intarr(w*2+1)+1)#indgen(1,w*2+1))^2)^0.5
D_D_all=reform(D_D_all,(w*2+1)*(w*2+1))

similar_th=fltarr(nb)
for iband=0,nb-1,1 do begin
  similar_th[iband,0]=stddev(fine1[*,*,iband])*2.0/num_class   
endfor

for i=location[0,isub],location[1,isub],1 do begin           ;retieve each target pixel
  for j=location[2,isub],location[3,isub],1 do begin
    if (fine1[i,j,background_band-1] ne background  ) then begin    ;do not process the background

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
           S_S=abs(wind_fine-wind_fine[ci,cj])
           similar_cand=similar_cand+S_S/(wind_fine[ci,cj]+0.00000001)
           ind_cand=where(S_S lt similar_th[iband])
           cand_band[ind_cand]=1
           position_cand=position_cand*cand_band
        endfor
        
        indcand=where(position_cand ne 0,number_cand0)  ;select similar pixel initially
        ; place the spatial distance measure between each calidated pixesls and the target pixel 
        if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not an integrate window
          distance_cand=((ci-col_wind)^2+(cj-row_wind)^2)^0.5+0.00001
        endif else begin
          distance_cand=D_D_all      ;integrate window
        endelse
        ;add a small weight from spatial distance to spectral distance to avoid all pixels in the window with same similarity. This happens in perfect simulated images
        combine_similar_cand=(similar_cand+0.00001)*(10.0+distance_cand/w);spatial distance has very small effect
        order_dis=sort(combine_similar_cand[indcand])
        number_cand=min([number_cand0,num_similar_pixel])
        ind_same_class=indcand[order_dis[0:number_cand-1]]           ; select the N most similar samples
        
                    
      ;normalize these distances
       D_D_cand=distance_cand[ind_same_class]
       C_D=(1.0+D_D_cand/w)*(similar_cand[ind_same_class]+1.0)
       C_D=1.0/C_D
       weight=C_D/total(C_D)

      for iband=0,nb-1,1 do begin
        ;predict the value
        change_cand=(change_21[ai:bi,aj:bj,iband])[ind_same_class]
        fine2[i-location[0,isub],j-location[2,isub],iband]=fine1[i,j,iband]+total(weight*change_cand)

        ;revise the abnormal prediction
        if (fine2[i-location[0,isub],j-location[2,isub],iband] lt DN_min ) then begin ;correct abnomal prediction
          another_predict=max([DN_min, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
          fine2[i-location[0,isub],j-location[2,isub],iband]=min([DN_max,another_predict])
        endif
        if (fine2[i-location[0,isub],j-location[2,isub],iband] gt DN_max ) then begin ;correct abnomal prediction
          another_predict=min([DN_max, fine1[i,j,iband]+(coarse2[i,j,iband]-coarse1[i,j,iband])])
          fine2[i-location[0,isub],j-location[2,isub],iband]=max([DN_min,another_predict])
        endif
      endfor


    endif
  endfor

endfor

; change the type of prediction into the type same as the input image
case Data_Type Of
  1:fine2 = Byte(fine2)    ;  BYTE  Byte
  2:fine2 = FIX(fine2)     ;  INT  Integer
  3:fine2 = LONG(fine2)    ;  LONG  Longword integer
  4:fine2 = FLOAT(fine2)   ;  FLOAT  Floating point
  5:fine2 = DOUBLE(fine2)  ;  DOUBLE  Double-precision floating
  6:fine2 = COMPLEX(fine2); complex, single-precision, floating-point
  9:fine2 = DCOMPLEX(fine2);complex, double-precision, floating-point
  12:fine2 = UINT(fine2)   ; unsigned integer vector or array
  13:fine2 = ULONG(fine2)   ;  unsigned longword integer vector or array
  14:fine2 = LONG64(fine2)   ;a 64-bit integer vector or array
  15:fine2 = ULONG64(fine2)   ;an unsigned 64-bit integer vector or array
EndCase

print,'finish final prediction ',isub+1,' block'
tempoutname1=temp_file+'\temp_blended'
Envi_Write_Envi_File,fine2,Out_Name = tempoutname1+strtrim(isub+1,1)
envi_file_mng, id=Fid11, /remove, /delete
envi_file_mng, id=Fid12, /remove, /delete
envi_file_mng, id=Fid13, /remove, /delete
envi_file_mng, id=Fid14, /remove, /delete
envi_file_mng, id=Fid15, /remove, /delete

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

use_see_through = replicate(1L,n_ns*n_nl)
see_through_val = replicate(0L,n_ns*n_nl)

out_name=FileName5+'_FSDAF'
envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
  dims=mdims, out_name=out_name, xsize=xsize, $
  ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
  out_dt=Data_Type, pixel_size=pixel_size, $
  background=0, see_through_val=see_through_val, $
  use_see_through=use_see_through

for i=0,n_ns*n_nl-1,1 do begin
  envi_file_mng, id=mfid[i], /remove, /delete
endfor


end