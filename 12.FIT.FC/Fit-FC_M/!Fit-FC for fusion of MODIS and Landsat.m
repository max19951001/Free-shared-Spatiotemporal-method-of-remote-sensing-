%%%%This is the code for Fit-FC produced by Dr Qunming Wang; Email: wqm11111@126.com
%%%%Copyright belongs to Qunming Wang
%%%%When using the code, please cite the fowllowing papers
%%%%Q. Wang, P. M. Atkinson. Spatio-temporal fusion for daily Sentinel-2 images. Remote Sensing of Environment, 2018, 204: 31–42.
clear all;
load MODIS_t1;
load MODIS_t2;
load Landsat_t1;
s=16;%%%s is the zoom ratio between MODIS and Landsat
%%%%%%MODIS need to be interpolated to ensure the zoom ratio is 16
%%%%%%E.g., if Landsat is 1600*1600 pixels, MODIS needs to be interpolated to 100*100 pixels

J1=Landsat_t1;
I_MS1=MODIS_t1;
I_MS2=MODIS_t2;
[Number_row,Number_col,Number_bands]=size(Landsat_t1);

W=input('Enter the radius of the window size in local regression: ');%%%E.g., if W=1, the window size if 3*3
for i=1:Number_bands
    [a,b,q]=guidedfilter_MS_low(I_MS1(:,:,i), I_MS2(:,:,i), W);
    I_MS0(:,:,i)=a.*I_MS1(:,:,i)+b;
    Z1(:,:,i)=imresize(a,s,'nearest').*J1(:,:,i)+imresize(b,s,'nearest');
    RB(:,:,i)=I_MS2(:,:,i)-I_MS0(:,:,i);%%%RB is the residual of the regression model
end

Z2_interpolated=Z1+imresize(RB,s,'bicubic');

w0=20;
N_S=20;
A=(2*w0+1)/2;
Z0=zeros(size(Z1));
B1=input('Enter the number of the red band in the Landsta cube: ');
B2=input('Enter the number of the NIR band in the Landsta cube: ');

tic
for i=1:Number_bands
    Z(:,:,i)=STARFM_fast_2016_v2(Z0(:,:,i),Z0(:,:,i),Z2_interpolated(:,:,i),J1(:,:,B1),J1(:,:,B2),w0,N_S,A);
end
alltime=toc
FalseColorf=Z(:,:,[4 3 2]);xf=imadjust(FalseColorf,stretchlim(FalseColorf),[]);figure,imshow(xf);
