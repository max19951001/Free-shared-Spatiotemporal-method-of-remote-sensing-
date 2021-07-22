
;---------------------------------------------------------------------------------
;                            ESTARFM PROGRAM
;               Using two pairs of fine and coarse images
;         the program can be used for whole TM scene and VI index product
;   A fast version of ESTARFM code: improve the effeciency of original ESTARFM 
;   through controling the number of similar pixels
;        Developed by (1) Zhu Xiaolin,email: zhuxiaolin55@gmail.com
;             Department of Land Surveying and Geo-Informatics
;             The Hong Kong Polytechnic University
;             
;            Debugging history: 
;            1)5/12/2012 correct the abnormal prediction
;            2)9/29/2013 correct the spatial distance calculation for integrate window 
;            3)7/10/2014 correct the abnormal value of spectral distance and use all bands to indentify background
;            4)2/13/2017 add one parameter to specify the value of pixels of background or missing
;            5)1/1/2018  improve efficiency and modified the weight for fusing VI index
;            6)3/11/2008 correct a bug in spatial distance caculation
;            7)7/27/2018 correct a bug of abnormal coversion coefficent estimation when the two input pairs are too similar
;            8)7/27/2018 improve the prediction when no enough simular pixels are selected
;                        
;Please cite the reference: Xiaolin Zhu, Jin Chen, Feng Gao, & Jeffrey G Masek.
;An enhanced spatial and temporal adaptive reflectance fusion model for complex
;heterogeneous regions. Remote Sensing of Environment,2010,114,2610-2623
;
;                     Copyright belongs to Xiaolin Zhu
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

;-------------------------------------------------------------------
;                       main program
;-------------------------------------------------------------------

pro  ESTARFM_FAST

 ;please set the following parameters
;----------------------------------------------------------------------
 w=25.0                 ;set the haif window size, if 25, the window size is 25*2+1=51 fine pixels
 num_class=4.0          ;set the estimated number of classes, please set a larger value if blending images with very few bands
 num_similar_pixel=20   ;set number of similar pixels, a smaller value is faster but accuracy may be lower
 DN_min=0               ;set the range of DN value of the image,If byte, 0 and 255
 DN_max=10000.0
 background=-9999       ;the value of background and missng pixels in both MODIS and Landsat images
 patch_long=500         ;set the size of each block,if process whole ETM scene, set 500-1000
 temp_file='D:\temp'    ;set the temporary file location, temporary files will be deleted after the work
;------------------------------------------------------------------------


 ;open the fine image of the first pair
  FileName1 = Dialog_PickFile(title = 'open the fine image of the first pair:')
  envi_open_file,FileName1,r_fid=fid
  envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims
  map_info = envi_get_map_info(fid=fid)
  orig_ns=ns
  orig_nl=nl
  n_ns=ceil(float(ns)/patch_long)
  n_nl=ceil(float(nl)/patch_long)

  ind_patch=intarr(4,n_ns*n_nl)
  for i_ns=0,n_ns-1,1 do begin
    for i_nl=0,n_nl-1,1 do begin
        ind_patch[0,n_ns*i_nl+i_ns]=i_ns*patch_long
        ind_patch[1,n_ns*i_nl+i_ns]=min([ns-1,(i_ns+1)*patch_long-1])
        ind_patch[2,n_ns*i_nl+i_ns]=i_nl*patch_long
        ind_patch[3,n_ns*i_nl+i_ns]=min([nl-1,(i_nl+1)*patch_long-1])
    endfor
  endfor

  tempoutname=temp_file+'\temp_F1'

  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
  endfor

  envi_file_mng, id=fid, /remove


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


    ;open the fine image of the second pair
  ;-----------------------------------------------------------
 FileName3 = Dialog_PickFile(title = 'open the fine image of the second pair:')
  envi_open_file,FileName3,r_fid=fid
  tempoutname=temp_file+'\temp_F2'
  pos=indgen(nb)
  for isub=0,n_ns*n_nl-1,1 do begin
      dims=[-1,ind_patch[0,isub],ind_patch[1,isub],ind_patch[2,isub],ind_patch[3,isub]]
      envi_doit, 'resize_doit', fid=fid, pos=pos, dims=dims, interp=0, rfact=[1,1], $
      out_name=tempoutname+strtrim(isub+1,1), r_fid=r_fid1
      envi_file_mng, id=r_fid1, /remove
  endfor
   envi_file_mng, id=fid, /remove

  ;open the coarse image of the second pair
  ;-----------------------------------------------------------

  FileName4 = Dialog_PickFile(title = 'open the coarse image of the second pair:')
  envi_open_file,FileName4,r_fid=fid
  tempoutname=temp_file+'\temp_C2'
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


;------------------------------------------------------------------
        ;process  each block
;-------------------------------------------------------------------
 t0=systime(1)                  ;the initial time of program running

print,'there are total',n_ns*n_nl,' blocks'

for isub=0,n_ns*n_nl-1,1 do begin

;open each block image

    FileName = temp_file+'\temp_F1'
     GetData,ImgData=fine1,ns = ns,nl = nl,nb = nb,Data_Type = Data_Type,FileName = FileName+strtrim(isub+1,1),Fid = Fid1
    fine1=float(fine1)

    FileName = temp_file+'\temp_C1'
    GetData,ImgData=coarse1,FileName = FileName+strtrim(isub+1,1),Fid = Fid2
    coarse1=FLOAT(coarse1)

    FileName = temp_file+'\temp_F2'
    GetData,ImgData=fine2,FileName = FileName+strtrim(isub+1,1),Fid = Fid3
    fine2=FLOAT(fine2)

    FileName = temp_file+'\temp_C2'
    GetData,ImgData=coarse2,FileName = FileName+strtrim(isub+1,1),Fid = Fid4
    coarse2=FLOAT(coarse2)

    FileName = temp_file+'\temp_C0'
    GetData,ImgData=coarse0,FileName = FileName+strtrim(isub+1,1),Fid = Fid5
    coarse0=FLOAT(coarse0)

    fine0=fltarr(ns,nl,nb)     ;place the blended result


    ;row index of images
    row_index=intarr(ns,nl)
    for i=0,nl-1,1 do begin
      row_index[*,i]=i
    endfor
    ;column index of images
    col_index=intarr(ns,nl)
    for i=0,ns-1,1 do begin
      col_index[i,*]=i
    endfor
    
   ;compute the uncertainty,0.2% of each band is uncertain
    uncertain=(DN_max*0.002)*(2^0.5)

     similar_th=fltarr(nb,2)          ;compute the threshold of similar pixel seeking

     for iband=0,nb-1,1 do begin
       similar_th[iband,0]=stddev(fine1[*,*,iband])*2.0/num_class   ;pair 1
       similar_th[iband,1]=stddev(fine2[*,*,iband])*2.0/num_class   ;pair 2
     endfor

     ;compute the distance of each pixel in the window with the target pixel (integrate window)
      D_D_all=1.0+((w-indgen(w*2+1)#(intarr(1,w*2+1)+1))^2+(w-(intarr(w*2+1)+1)#indgen(1,w*2+1))^2)^0.5/float(w)
    
     ;find interaction of valid pixels of all input images: exclude missing pixels and background
     valid_index=bytarr(ns,nl)
     ind_valid=where(fine1[*,*,0] ne background and fine2[*,*,0] ne background and coarse1[*,*,0] ne background $
      and coarse2[*,*,0] ne background and coarse0[*,*,0] ne background,num_valid)
     if (num_valid gt 0) then valid_index[ind_valid]=1   ;mark good pixels in all images
     
     for j=0,nl-1,1 do begin              ;retieve each target pixel
       for i=0,ns-1,1 do begin
       
        if (valid_index[i,j] eq 1) then begin    ;do not process the background

          ai=max([0,i-w])       ; the window location
          bi=min([ns-1,i+w])
          aj=max([0,j-w])
          bj=min([nl-1,j+w])

          ind_wind_valid=where(valid_index[ai:bi,aj:bj] eq 1)
          position_cand=intarr((bi-ai+1)*(bj-aj+1))+1  ;place the location of each similar pixel
          similar_cand=fltarr((bi-ai+1)*(bj-aj+1)) ;pleace the similarity measure between each pixel and the target pixel
          row_wind=row_index[ai:bi,aj:bj]
          col_wind=col_index[ai:bi,aj:bj]
          
          ;searching for similar pixels
          for ipair=0,1,1 do begin
             for iband=0,nb-1,1 do begin
                 cand_band=intarr((bi-ai+1)*(bj-aj+1))
                 case ipair of
                  0:S_S=abs(fine1[ai:bi,aj:bj,iband]-fine1[i,j,iband])
                  1:S_S=abs(fine2[ai:bi,aj:bj,iband]-fine2[i,j,iband])
                 endcase
                 similar_cand=similar_cand+S_S/(similar_th[iband,ipair]+0.00000001)
                 ind_cand=where(S_S lt similar_th[iband,ipair])
                 cand_band[ind_cand]=1
                 position_cand=position_cand*cand_band
             endfor
          endfor
          cand_band=0
        
          indcand0=where(position_cand ne 0 and valid_index[ai:bi,aj:bj] eq 1,number_cand0)  ;select similar pixel initially
          order_dis=sort(similar_cand[indcand0])
          number_cand=min([number_cand0,num_similar_pixel])
          indcand=indcand0[order_dis[0:number_cand-1]]           ; select the N most similar samples
        
          if (number_cand gt 5) then begin

             S_D_cand=fltarr(number_cand)                ;compute the correlation
             x_cand=col_wind[indcand]
             y_cand=row_wind[indcand]
             finecand=fltarr(number_cand,nb*2)
             coasecand=fltarr(number_cand,nb*2)
             for ib=0,nb-1, 1 do begin
               finecand[*,ib]=(fine1[ai:bi,aj:bj,ib])[indcand]
               finecand[*,ib+nb]=(fine2[ai:bi,aj:bj,ib])[indcand]
               coasecand[*,ib]=(coarse1[ai:bi,aj:bj,ib])[indcand]
               coasecand[*,ib+nb]=(coarse2[ai:bi,aj:bj,ib])[indcand]
             endfor
             
             if (nb eq 1) then begin  ; for images with one band, like NDVI
               S_D_cand=1.0-0.5*(abs((finecand[*,0]-coasecand[*,0])/(finecand[*,0]+coasecand[*,0]))+abs((finecand[*,1]-coasecand[*,1])/(finecand[*,1]+coasecand[*,1])))
             endif else begin   
              ; for images with multiple bands             
               sdx=stddev(finecand,DIMENSION=2)
               sdy=stddev(coasecand,DIMENSION=2)         
               meanx=mean(finecand,DIMENSION=2)
               meany=mean(coasecand,DIMENSION=2)
               x_meanx=fltarr(number_cand,nb*2)
               y_meany=fltarr(number_cand,nb*2)
               for ib=0,nb*2-1, 1 do begin
                 x_meanx[*,ib]=finecand[*,ib]-meanx
                 y_meany[*,ib]=coasecand[*,ib]-meany
               endfor     
               S_D_cand=nb*2.0*mean(x_meanx*y_meany,DIMENSION=2)/(sdx*sdy)/(nb*2.0-1) 
             endelse
               ind_nan=where(S_D_cand ne S_D_cand,num_nan)
               if (num_nan gt 0) then S_D_cand[ind_nan]=0.5 ;correct the NaN value of correlation

              D_D_cand=fltarr(number_cand)        ;spatial distance
              if ((bi-ai+1)*(bj-aj+1) lt (w*2.0+1)*(w*2.0+1)) then begin   ;not an integrate window
                 D_D_cand=1.0+((i-x_cand)^2+(j-y_cand)^2)^0.5/float(w)              
              endif else begin
                 D_D_cand[0:number_cand-1]=D_D_all[indcand]      ;integrate window
              endelse
              C_D=(1.0-S_D_cand)*D_D_cand+0.0000001            ;combined distance
              weight=(1.0/C_D)/total(1.0/C_D)

              for iband=0,nb-1,1 do begin     ;compute V
                  fine_cand=[(fine1[ai:bi,aj:bj,iband])[indcand],(fine2[ai:bi,aj:bj,iband])[indcand]]
                  corse_cand=[(coarse1[ai:bi,aj:bj,iband])[indcand],(coarse2[ai:bi,aj:bj,iband])[indcand]]
                  coarse_change=abs(mean((coarse1[ai:bi,aj:bj,iband])[indcand])-mean((coarse2[ai:bi,aj:bj,iband])[indcand]))         
                  if ( coarse_change ge DN_max*0.02) then begin ;to ensure changes in coarse image large enough to obtain the conversion coefficient
                        regress_result=regress(corse_cand,fine_cand,FTEST=fvalue)
                        sig=1.0-f_pdf(fvalue,1,number_cand*2-2)
                       ;correct the result with no significancy or inconsistent change or too large value  
                        if (sig le 0.05 and regress_result[0] gt 0 and regress_result[0] le 5) then begin
                             V_cand=regress_result[0]
                        endif else begin
                             V_cand=1.0
                        endelse
                  endif else begin
                        V_cand=1.0
                  endelse

                    ; compute the temporal weight
                     difc_pair1=abs(mean((coarse0[ai:bi,aj:bj,iband])[ind_wind_valid])-mean((coarse1[ai:bi,aj:bj,iband])[ind_wind_valid]))+0.01^5
                     difc_pair2=abs(mean((coarse0[ai:bi,aj:bj,iband])[ind_wind_valid])-mean((coarse2[ai:bi,aj:bj,iband])[ind_wind_valid]))+0.01^5
                     T_weight1=(1.0/difc_pair1)/(1.0/difc_pair1+1.0/difc_pair2)
                     T_weight2=(1.0/difc_pair2)/(1.0/difc_pair1+1.0/difc_pair2)

                    ;predict from pair1
                     coase0_cand=(coarse0[ai:bi,aj:bj,iband])[indcand]
                     coase1_cand=(coarse1[ai:bi,aj:bj,iband])[indcand]
                     fine01=fine1[i,j,iband]+total(weight*V_cand*(coase0_cand-coase1_cand))
                     ;predict from pair2
                     coase2_cand=(coarse2[ai:bi,aj:bj,iband])[indcand]
                     fine02=fine2[i,j,iband]+total(weight*V_cand*(coase0_cand-coase2_cand))
                     ;the final prediction
                     fine0[i,j,iband]=T_weight1*fine01+T_weight2*fine02
                     ;revise the abnormal prediction
                     if (fine0[i,j,iband] le DN_min or fine0[i,j,iband] ge DN_max) then begin
                        fine01=total(weight*(fine1[ai:bi,aj:bj,iband])[indcand])
                        fine02=total(weight*(fine2[ai:bi,aj:bj,iband])[indcand])  
                        fine0[i,j,iband]=T_weight1*fine01+T_weight2*fine02
                     endif
                  endfor
               endif else begin   ;for the case of no enough similar pixel selected
                     for iband=0,nb-1,1 do begin  
                     ; compute the temporal weight
                        difc_pair1=mean((coarse0[ai:bi,aj:bj,iband])[ind_wind_valid])-mean((coarse1[ai:bi,aj:bj,iband])[ind_wind_valid])+0.01^5
                        difc_pair1_a=abs(difc_pair1)
                        difc_pair2=mean((coarse0[ai:bi,aj:bj,iband])[ind_wind_valid])-mean((coarse2[ai:bi,aj:bj,iband])[ind_wind_valid])+0.01^5
                        difc_pair2_a=abs(difc_pair2)
                        T_weight1=(1.0/difc_pair1_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        T_weight2=(1.0/difc_pair2_a)/(1.0/difc_pair1_a+1.0/difc_pair2_a)
                        fine0[i,j,iband]=T_weight1*(fine1[i,j,iband]+difc_pair1)+T_weight2*(fine2[i,j,iband]+difc_pair2)
                     endfor
               endelse
              endif             
             endfor
            endfor

     ; change the type of prediction into the type same as the input image
    case Data_Type Of
        1:fine0 = Byte(fine0)    ;  BYTE  Byte
        2:fine0 = FIX(fine0)     ;  INT  Integer
        3:fine0 = LONG(fine0)    ;  LONG  Longword integer
        4:fine0 = FLOAT(fine0)   ;  FLOAT  Floating point
        5:fine0 = DOUBLE(fine0)  ;  DOUBLE  Double-precision floating
        6:fine0 = COMPLEX(fine0); complex, single-precision, floating-point
        9:fine0 = DCOMPLEX(fine0);complex, double-precision, floating-point
        12:fine0 = UINT(fine0)   ; unsigned integer vector or array
        13:fine0 = ULONG(fine0)   ;  unsigned longword integer vector or array
        14:fine0 = LONG64(fine0)   ;a 64-bit integer vector or array
        15:fine0a = ULONG64(fine0)   ;an unsigned 64-bit integer vector or array
    EndCase

       print,'finished ',isub+1,' block'
         tempoutname1=temp_file+'\temp_blended'
         Envi_Write_Envi_File,fine0,Out_Name = tempoutname1+strtrim(isub+1,1)
         envi_file_mng, id=Fid1, /remove, /delete
         envi_file_mng, id=Fid2, /remove, /delete
         envi_file_mng, id=Fid3, /remove, /delete
         envi_file_mng, id=Fid4, /remove, /delete
         envi_file_mng, id=Fid5, /remove, /delete
endfor

;;--------------------------------------------------------------------------------------
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
      x0[isub]=ind_patch[0,isub]
      y0[isub]=ind_patch[2,isub]
  endfor

    xsize = orig_ns
    ysize = orig_nl
    pixel_size = [1.,1.]

    use_see_through = replicate(1L,n_ns*n_nl)
    see_through_val = replicate(0L,n_ns*n_nl)

    out_name=FileName5+'_ESTARFM_FAST'
    envi_doit, 'mosaic_doit', fid=mfid, pos=mpos, $
    dims=mdims, out_name=out_name, xsize=xsize, $
    ysize=ysize, x0=x0, y0=y0, georef=0,MAP_INFO=map_info, $
    out_dt=Data_Type, pixel_size=pixel_size, $
    background=0, see_through_val=see_through_val, $
    use_see_through=use_see_through

    for i=0,n_ns*n_nl-1,1 do begin
      envi_file_mng, id=mfid[i], /remove, /delete
    endfor

print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'


end