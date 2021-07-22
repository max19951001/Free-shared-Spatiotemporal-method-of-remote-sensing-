

clear all
clc
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%

DN_min = 0;  %% minimal DN value
DN_max = 10000;   %% maximal DN value; if using reflectance, use 1 as DN_max

class=4;  %% number of classes

 
%%%%%%%%%

[filename, pathname] = uigetfile( '*.tif', 'input the fine image of the first pair' );
FR_T1 = imread( [pathname, filename] );

%%% for output %%%
infilename = ( [pathname, filename] );

[filename, pathname] = uigetfile( '*.tif', 'input the coarse image of the first pair' );
CR_T1 = imread( [pathname, filename] );


[filename, pathname] = uigetfile( '*.tif', 'input the coarse image at the prediction date' );
CR_T2 = imread( [pathname, filename] );
output_name = filename;

% [filename, pathname] = uigetfile( '*.tif', 'input the fine image at the prediction date for validation' );
% FR_T2 = imread( [pathname, filename] );


CR_T1=double(CR_T1);
CR_T2=double(CR_T2);

FR_T1=double(FR_T1);
% FR_T2=double(FR_T2);

[xH,yH,bands]=size(FR_T1);
[xL,yL,bands]=size(CR_T1);
Factor=xH/xL;
 

%%%%%%%%%%%%%%%%%   Clustering  %%%%%%%%%%%%%%%%

disp('Cluster T1 FR image')

iteration=10; %%  number of iteration in K-means clustering 

minimum_n=xH*yH/500;  %% the minimal numer of pixels in a cluster
FR_T1_V = zeros(xH*yH,bands);
for i=1:bands
    FR_T1_V(:,i) = reshape(FR_T1(:,:,i),xH*yH,1);
end  
[centroid, result] = Clustering(FR_T1_V, 'kmeans', class , iteration); 
map_T1=reshape(result,xH,yH);
clear iteration minimum_n FR_T1_V centroid result

%%%%%%%%%%%%%%%%%%%%%% Create T1 FR fraction images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Create T1 FR fraction images')


n_HRpixels = zeros(class,1);
for c=1:class
    n_HRpixels(c) = sum(sum(map_T1==c));
end

meanV = zeros(bands,class);
Cov_M = zeros(bands,bands,class);
for c=1:class
    n=0;
    data_p = zeros(n_HRpixels(c),bands);
    for i=1:xH
        for j=1:yH
            if map_T1(i,j)==c
                n=n+1;
                data_p(n,:) = squeeze(FR_T1(i,j,:))';
            end
        end
    end
    Cov_M(:,:,c) = cov(data_p);
    meanV(:,c) = mean(data_p)';
end

FR_fraction_T1=zeros(xH,yH,class);
v_inv=zeros(bands,bands,class);
for c=1:class
    v_inv(:,:,c) = inv(Cov_M(:,:,c));
end

a=zeros(class,1);
for i=1:xH
    for j=1:yH
        s=squeeze(FR_T1(i,j,:));
        for c=1:class
           u = meanV(:,c);
           a(c) = ( (s-u)'*v_inv(:,:,c)*(s-u) ).^(-1);
        end
         FR_fraction_T1(i,j,:)=a/sum(a);
    end
end

clear  n_HRpixels meanV u a s v_inv Cov_M meanV n data_p
%%%%%%%%%%%%%%%% Degrade T1 FR fraction images to CR scale %%%%%%%%%%%%%%%%%%%%

 CR_fraction_T1=zeros(xL,yL,class);
 for c=1:class
     for i=1:xL
         for j=1:yL
             tmp = sum(sum(  FR_fraction_T1( (i-1)*Factor+1:i*Factor,(j-1)*Factor+1:j*Factor,c) ))/Factor^2;
             CR_fraction_T1(i,j,c)=tmp;
         end
     end
 end
 clear tmp  
%%%%%%%%%%%%%%%%%%%%  Estimate T2  CR endmembers   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimate T2 CR endmembers')

CR_endmember_T2 =  Estimate_endmember_T2(CR_fraction_T1 ,CR_T2 ,CR_T1 ,xL,yL,class,bands);

%%%%%%%%%%%%%%%%%%%%%%%  Estimate T1 to T2 FR endmember change  %%%%%%%%%%%%%%%%%%%%%%%

disp('Estimate T1 to T2 FR endmember change')

FR_endmember_change = Estimate_endmember_change(CR_fraction_T1 ,CR_T2 ,CR_T1 ,xL,yL,class,bands); 

%%%%%%%%%%%%%%%%%%%%%%  unmix T2 CR image to fraction images   %%%%%%%%%%%%%

disp('Unmix T2 CR image to fraction images')

n_iteration=100;  %% iteration number used in LSE
CR_fraction_T2=CLMM_Unmix(CR_T2,CR_endmember_T2,class,bands,xL,yL,n_iteration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% predict T2 time FR fraction images  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Predict T2 time FR fraction images')


fractionchange_Interpolation  =  Downscale_FractionChangeImages(CR_fraction_T1,CR_fraction_T2,FR_fraction_T1,class,xH,yH,Factor);
FR_fraction_T2 = fractionchange_Interpolation + FR_fraction_T1;

FR_fraction_T2(FR_fraction_T2>1)=1;
FR_fraction_T2(FR_fraction_T2<0)=0;

temp = sum(FR_fraction_T2,3);
for c=1:class
    FR_fraction_T2(:,:,c) = FR_fraction_T2(:,:,c) ./temp;
end
clear temp
%%%%%%%%%%%%%%%%%%%%%%% Estimate T1 and T2 endmembers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimate T1 and T2 endmembers')


FR_endmember_T1= Estimate_endmember_T1(FR_fraction_T1,FR_T1,xH,yH,class,bands);
FR_endmember_T2 = FR_endmember_T1 + FR_endmember_change;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create T2 FR temporal prediciton image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Create T2 FR temporal prediciton image ')


temp=zeros(xH,yH,bands);
for b=1:bands
    for c=1:class
        temp(:,:,b) = temp(:,:,b) + FR_endmember_T2(b,c) * FR_fraction_T2(:,:,c) - FR_endmember_T1(b,c) * FR_fraction_T1(:,:,c) ;
    end
end
F2TP=FR_T1+temp; 
clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%   Spatial interpolation of the T2 CR image  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Spatial interpolation of the T2 CR image ')


F2SP=zeros(xH,yH,bands);
  for b=1:bands
        F2SP(:,:,b) = SpatialInterpolation(CR_T2(:,:,b),Factor,'cubic');
  end
  
%%%%%%%%%%%%  Distribute the residuals %%%%%%%%%%%
%%%%%%%%%%%%  The same method used in FSDAF is used here  %%%%%%%%%%% 

disp('Distribute the residuals ')


m = Factor^2;
HI = zeros(xH,yH);
Result=zeros(xH+2*Factor,yH+2*Factor);       
Result(Factor+1:xH+Factor,Factor+1:yH+Factor)=map_T1; 
for i=Factor+1:Factor+xH
    for j=Factor+1:Factor+yH
        SD = Result(i-Factor/2:i+Factor/2,j-Factor/2:j+Factor/2 );
        HI(i-Factor,j-Factor) = sum(sum(SD==Result(i,j)))/(Factor+1)^2;
    end
end 

Eho = zeros(xH,yH,bands);
Eho = F2SP - F2TP;

Ehe = zeros(xH,yH,bands);
R = zeros(xL,yL,bands);
for i=1:xL
    for j=1:yL
        tmp1 = F2TP( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: );
        tmp2 = FR_T1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: );
        R(i,j,:) = squeeze(CR_T2(i,j,:))-squeeze(CR_T1(i,j,:)) - (  squeeze(sum(sum(  tmp1 ))) - squeeze(sum(sum( tmp2 ))) ) / Factor^2;
    end
end
 
 for i=1:xL
     for j=1:yL
         for b=1:bands
             Ehe((i-1)*Factor+1:i*Factor,(j-1)*Factor+1:j*Factor,b) = R(i,j,b);
         end          
     end
 end
 
 CW = zeros(xH,yH,bands);
for b=1:bands
    CW(:,:,b) = abs(Eho(:,:,b)) .* HI + abs(Ehe(:,:,b)) .* (ones(xH,yH)-HI) + 0.0000001;
end

W = zeros(xH,yH,bands);
for b = 1:bands
    for i=1:xH
        for j=1:yH
            W(i,j,b) = CW(i,j,b) / sum(sum( abs( CW( ( ceil(i/Factor)-1)*Factor+1: ceil(i/Factor)*Factor, (ceil(j/Factor)-1)*Factor+1:ceil(j/Factor)*Factor,b  ) )));
        end
    end
end

r = zeros(xH,yH,bands);
for i=1:xH
    for j=1:yH
        for b=1:bands
              r(i,j,b) =  m * R( ceil(i/Factor),ceil(j/Factor), b ) * W(i,j,b) ;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DeltaF=zeros(xH,yH,bands);
for b=1:bands
    for c=1:class
        DeltaF(:,:,b) = DeltaF(:,:,b) + FR_endmember_T2(b,c) * FR_fraction_T2(:,:,c) - FR_endmember_T1(b,c) * FR_fraction_T1(:,:,c) ;
    end
end
DeltaF=DeltaF+r;

clear Result r m  SD R Ehe Eho W tmp1 tmp2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Using similar pixels in the final prediction   %%%%%%%%%%%%%%%%%%%%%%%%

disp('Using similar pixels in the final prediction')


FR_T1_ext = zeros(xH+4*Factor,yH+4*Factor,bands);       
FR_T1_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor,:)=FR_T1; 

DeltaF_ext = zeros(xH+4*Factor,yH+4*Factor,bands);       
DeltaF_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor,:)=DeltaF; 

%%%%%%%%%%%%  map_T1 is used to distribute the error;%%%%%%%%%%%%  
%%%%%%  improvements to use other map with change information could be developed %%%%

map_T1_ext = zeros(xH+4*Factor,yH+4*Factor);      
map_T1_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor)=map_T1;
HI_ext = zeros(xH+4*Factor,yH+4*Factor);      
HI_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor)=HI; 

Sk = zeros(4*Factor+1,4*Factor+1);
n=20; %% number of similar pixels
w = Factor;

Matrix = zeros(2*Factor+1,2*Factor+1);
for i=1:2*Factor+1
    for j=1:2*Factor+1
         Matrix(i,j) =  sqrt(  (i-(Factor+1))^2 + (j-(Factor+1))^2  ) ;
    end
end

A=2;
Dk = zeros(2*Factor+1,2*Factor+1);
for i=1:2*Factor+1
    for j=1:2*Factor+1
         Dk(i,j) =  1 + Matrix(i,j)/ (w/A);
    end
end

wk_T1 = zeros(2*Factor+1,2*Factor+1);
wk_T2 = zeros(2*Factor+1,2*Factor+1);
wk = zeros(2*Factor+1,2*Factor+1);

SFSDAF_prediction = zeros(xH,yH,bands);
for i=2*Factor+1:2*Factor+xH
    for j=2*Factor+1:2*Factor+yH
        SD = FR_T1_ext( i-Factor:i+Factor,j-Factor:j+Factor,:  );
        for b=1:bands
            SD(:,:,b)  = SD(:,:,b) - FR_T1_ext(i,j,b) +0.0000000001;
        end
        SD=abs(SD);
                
        Sk_FR = sum(SD,3);     
        
            Sk = Sk_FR;
            Sk_row = reshape(Sk,(2*Factor+1)*(2*Factor+1),1);
            tmp = sort(Sk_row);
            threhold_n = tmp(n+1);
            
           SD_map = Sk<=threhold_n;
            wk = (  (1./Dk.* SD_map) / sum(sum(1./Dk.* SD_map))  );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for b=1:bands
            SFSDAF_prediction(i-2*Factor,j-2*Factor,b) = FR_T1(i-2*Factor,j-2*Factor,b) + sum(sum(DeltaF_ext(i-Factor:i+Factor,j-Factor:j+Factor,b) .* wk)); 
            
            if SFSDAF_prediction(i-2*Factor,j-2*Factor,b)< DN_min
                another_prediction = max( DN_min, FR_T1(i-2*Factor,j-2*Factor,b) + ( CR_T2( ceil(i/Factor)-2, ceil(j/Factor)-2 , b)  - CR_T1( ceil(i/Factor)-2, ceil(j/Factor)-2 , b) ) );
                SFSDAF_prediction(i-2*Factor,j-2*Factor,b) = min( DN_max, another_prediction );
            end
            
            if SFSDAF_prediction(i-2*Factor,j-2*Factor,b)> DN_max
                another_prediction = min( DN_max, FR_T1(i-2*Factor,j-2*Factor,b) + ( CR_T2( ceil(i/Factor)-2, ceil(j/Factor)-2 , b)  - CR_T1( ceil(i/Factor)-2, ceil(j/Factor)-2 , b) ) );
                SFSDAF_prediction(i-2*Factor,j-2*Factor,b) = max( DN_min, another_prediction );
            end
            %%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
end



%%%%%%%%%%%%%% output %%%%%%%%%%%%%

disp('Outputting')

info = imfinfo(infilename);

if ~isfield(info,'GeoAsciiParamsTag')   
    
    for b=1:bands
         name = strcat('SFSDAF_prediction_b',num2str(b),'.txt');
         WriteGridASCII( name,yH,xH,SFSDAF_prediction(:,:,b) )
    end  
    
else    
   
    [A,R] = geotiffread(infilename);
    info = geotiffinfo(infilename);
    geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
    tiffTags = struct('TileLength',1024,'TileWidth',1024);
    outfilename = strcat('SFSDAF_', output_name);
    geotiffwrite(outfilename,SFSDAF_prediction,R,'GeoKeyDirectoryTag',geoTags, ...
    'TiffType','bigtiff','TiffTags',tiffTags) 

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except SFSDAF_prediction

% save('SFSDAF_prediction.mat', 'SFSDAF_prediction')
